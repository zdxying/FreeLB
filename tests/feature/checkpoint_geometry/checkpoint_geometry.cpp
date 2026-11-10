/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// checkpoint_geometry.cpp
// A snapshot now carries the geometry identity: grid, origin/extent, per-cell
// entity, boundary patches and a record digest.  The test writes such a
// snapshot from a 3D record, checks the header round-trip, then checks that the
// digest rejects a different record and that a matching restart continues
// bit-identically.
//
// The legacy path is covered too: a snapshot written without a record cannot be
// validated, so checkCompatible() (which has no record) must refuse it instead
// of resuming into a silently different geometry.

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>

#include "freelb.h"
#include "freelb.hh"
#include "../../test_util.h"
#include "utils/field_checksum.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t BulkFlag = std::uint8_t(2);
constexpr int Nx = 12;
constexpr int Ny = 12;
constexpr int Nz = 8;
constexpr int kPre = 100;
constexpr int kPost = 40;

using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
using CELL = Cell<T, LatSet, FIELDS>;
using BulkTask =
  tmp::Key_TypePair<BulkFlag,
                    collision::BGK<moment::rhoU<CELL>,
                                   equilibrium::SecondOrder<CELL>>>;
using TaskCollection = tmp::TupleWrapper<BulkTask>;
using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;

// One embedded solid box, so the record has a band, an entity and a patch.
flbgeo::VoxelGeometry<T> makeRecord(T solidShift) {
  const AABB<T, 3> domain(Vector<T, 3>{T(0), T(0), T(0)},
                          Vector<T, 3>{T(Nx), T(Ny), T(Nz)});
  flbgeo::VoxelGeometryBuilder<T> builder(Nx, Ny, Nz, domain, T(1), T(1.5));
  builder.addBox(domain, 1, 1, flbgeo::VoxelRole::Domain);
  const AABB<T, 3> solid(Vector<T, 3>{T(4) + solidShift, T(4), T(2)},
                         Vector<T, 3>{T(8) + solidShift, T(8), T(6)});
  builder.addBox(solid, 2, 2, flbgeo::VoxelRole::Solid);
  // The band lives on both sides of the surface, so the selector must reach
  // past the box to pick up the fluid cells that carry the entity.
  builder.addPatch(2, flbgeo::PatchKind::Wall,
                   solid.getExtended(Vector<T, 3>{T(1.5)}));
  return builder.build();
}

struct Sim {
  BlockGeometry3D<T> Geo;
  BlockFieldManager<FLAG, T, 3> FlagFM;
  BaseConverter<T> Conv;
  std::unique_ptr<BlockLatticeManager<T, LatSet, FIELDS>> Lattice;

  explicit Sim(const flbgeo::VoxelGeometry<T>& record)
      : Geo(record.getNx(), record.getNy(), record.getNz(), 1,
            record.getDomain(), record.getVoxelSize(), 1),
        FlagFM(Geo, VoidFlag), Conv(LatSet::cs2) {
    flbgeo::writeVoxelFlags(record, FlagFM, BulkFlag, VoidFlag);
    Conv.SimplifiedConverterFromRT(record.getNx(), T(0.1), T(0.56));
    ValuePack initValues(Conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
    Lattice.reset(
      new BlockLatticeManager<T, LatSet, FIELDS>(Geo, initValues, Conv));
  }

  void step() {
    Lattice->Stream();
    Lattice->template ApplyCellDynamics<NSTask>(FlagFM);
    Lattice->NormalFullCommunicate();
  }

  std::array<double, 3> popChecksum() {
    auto& popf =
      Lattice->getBlockLat(0).template getField<POP<T, LatSet::q>>();
    return frelb_diag::PopChecksumHost(popf, Lattice->getBlockLat(0).getN(),
                                       LatSet::q);
  }
};

}  // namespace

int main() {
  const flbgeo::VoxelGeometry<T> record = makeRecord(T(0));
  CHECK(record.countRole(flbgeo::VoxelRole::Domain) > 0);
  CHECK(record.countRole(flbgeo::VoxelRole::Solid) > 0);
  CHECK_EQ(record.getPatchCount(), std::size_t{1});
  CHECK(!record.getPatchCells(0).empty());

  const std::string dir = "ckpt";
  Sim A(record);
  for (int step = 0; step < kPre; ++step) A.step();
  flbckpt::writeWithGeometry(dir, kPre, record, *A.Lattice, A.FlagFM);

  // ---- header: the identity is stored and readable ----
  flbckpt::CheckpointHeader hdr;
  CHECK(flbckpt::readHeader(dir, -1, hdr));
  CHECK_EQ(hdr.version, flbckpt::kFormatVersion);
  CHECK(hdr.geometry.present);
  CHECK_EQ(hdr.geometry.nx, std::uint32_t(Nx));
  CHECK_EQ(hdr.geometry.ny, std::uint32_t(Ny));
  CHECK_EQ(hdr.geometry.nz, std::uint32_t(Nz));
  const auto& domain = record.getDomain();
  for (int i = 0; i < 3; ++i) {
    CHECK_NEAR(hdr.geometry.origin[i], static_cast<double>(domain.getMin()[i]),
               1e-12);
    CHECK_NEAR(hdr.geometry.max[i], static_cast<double>(domain.getMax()[i]),
               1e-12);
  }
  CHECK_NEAR(hdr.geometry.voxelSize, static_cast<double>(record.getVoxelSize()),
             1e-12);
  CHECK_EQ(hdr.geometry.digest, flbgeo::geometryDigest(record));
  CHECK_EQ(hdr.geometry.entity.size(), std::size_t(Nx * Ny * Nz));
  for (std::size_t id = 0; id < record.size(); ++id) {
    CHECK_EQ(hdr.geometry.entity[id], record.getEntity(id));
  }
  CHECK_EQ(hdr.geometry.patches.size(), std::size_t{1});
  CHECK_EQ(hdr.geometry.patches[0].entity, std::uint16_t{2});
  CHECK_EQ(static_cast<int>(hdr.geometry.patches[0].kind),
           static_cast<int>(flbgeo::PatchKind::Wall));
  const auto& patchCells = record.getPatchCells(0);
  CHECK_EQ(hdr.geometry.patches[0].cells.size(), patchCells.size());
  for (std::size_t i = 0; i < patchCells.size(); ++i) {
    CHECK_EQ(hdr.geometry.patches[0].cells[i],
             static_cast<std::uint64_t>(patchCells[i]));
  }

  // ---- validation: matching record passes, different record fails ----
  CHECK(flbckpt::checkCompatibleWithGeometry(dir, -1, record, *A.Lattice,
                                             A.FlagFM));
  // No record means no validation, so the legacy check must refuse the file.
  CHECK(!flbckpt::checkCompatible(dir, -1, *A.Lattice, A.FlagFM));
  const flbgeo::VoxelGeometry<T> moved = makeRecord(T(1));
  CHECK(flbgeo::geometryDigest(moved) != flbgeo::geometryDigest(record));
  CHECK(!flbckpt::checkCompatibleWithGeometry(dir, -1, moved, *A.Lattice,
                                              A.FlagFM));

  // ---- restart with the matching record continues bit-identically ----
  Sim B(record);
  const std::int64_t restored =
    flbckpt::readWithGeometry(dir, -1, record, *B.Lattice, B.FlagFM);
  CHECK_EQ(restored, std::int64_t(kPre));

  for (int step = 0; step < kPost; ++step) A.step();
  for (int step = 0; step < kPost; ++step) B.step();
  const std::array<double, 3> cA = A.popChecksum();
  const std::array<double, 3> cB = B.popChecksum();
  CHECK_EQ(cA[0], cB[0]);
  CHECK_EQ(cA[1], cB[1]);
  CHECK_EQ(cA[2], cB[2]);

  std::cout << "[checkpoint_geometry] digest=" << hdr.geometry.digest
            << " cells=" << hdr.geometry.entity.size()
            << " patches=" << hdr.geometry.patches.size()
            << " patch_cells=" << patchCells.size() << " sum=" << cA[0]
            << " sumsq=" << cA[1] << std::endl;

  return flbtest::summary("feature/checkpoint_geometry");
}
