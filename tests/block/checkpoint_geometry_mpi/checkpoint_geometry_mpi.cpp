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

// checkpoint_geometry_mpi.cpp
// Two-rank geometry identity: every rank writes the same VoxelGeometry record
// into its own checkpoint file, and a restart reads all rank files and checks
// that they agree before validating the caller's record.  The two-rank layout
// also keeps the block table rank-independent, so the restarted run can be
// compared bit-for-bit on each rank.
//
// The rank-count change itself is covered by block/checkpoint_decomp; this test
// covers the part that needs more than one rank: the per-rank identity blocks
// and their cross-file agreement.

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
constexpr int Nx = 16;
constexpr int Ny = 16;
constexpr int Nz = 10;
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

flbgeo::VoxelGeometry<T> makeRecord() {
  const AABB<T, 3> domain(Vector<T, 3>{T(0), T(0), T(0)},
                          Vector<T, 3>{T(Nx), T(Ny), T(Nz)});
  flbgeo::VoxelGeometryBuilder<T> builder(Nx, Ny, Nz, domain, T(1), T(1.5));
  builder.addBox(domain, 1, 1, flbgeo::VoxelRole::Domain);
  const AABB<T, 3> solid(Vector<T, 3>{T(5), T(5), T(2)},
                         Vector<T, 3>{T(10), T(10), T(8)});
  builder.addBox(solid, 2, 2, flbgeo::VoxelRole::Solid);
  builder.addPatch(2, flbgeo::PatchKind::Wall,
                   solid.getExtended(Vector<T, 3>{T(1.5)}));
  return builder.build();
}

struct Sim {
  BlockGeometryHelper3D<T> Helper;
  std::unique_ptr<BlockGeometry3D<T>> Geo;
  std::unique_ptr<BlockFieldManager<FLAG, T, 3>> FlagFM;
  BaseConverter<T> Conv;
  BlockLatticeManager<T, LatSet, FIELDS>* Lattice;

  explicit Sim(const flbgeo::VoxelGeometry<T>& record)
      : Helper(record.getNx(), record.getNy(), record.getNz(),
               record.getDomain(), record.getVoxelSize(), 8, 1, 0, false),
        Conv(LatSet::cs2), Lattice(nullptr) {
    Helper.CreateBlocks(2, 2, 2);
    Helper.LoadBalancing(mpi().getSize());
    Geo.reset(new BlockGeometry3D<T>(Helper));
    FlagFM.reset(new BlockFieldManager<FLAG, T, 3>(*Geo, VoidFlag));
    flbgeo::writeVoxelFlags(record, *FlagFM, BulkFlag, VoidFlag);
    Geo->InitVoxelMap(*FlagFM);
    Conv.SimplifiedConverterFromRT(record.getNx(), T(0.1), T(0.56));
    ValuePack initValues(Conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
    Lattice = new BlockLatticeManager<T, LatSet, FIELDS>(*Geo, initValues, Conv);
  }

  ~Sim() { delete Lattice; }

  void step() {
    Lattice->template ApplyCellDynamics<NSTask>(*FlagFM);
    Lattice->Stream();
    Lattice->NormalFullCommunicate();
  }

  // local checksum over every owned block; the same rank restarts with the
  // same block set, so the comparison is bit-exact
  std::array<double, 3> popChecksum() {
    std::array<double, 3> total{0, 0, 0};
    for (std::size_t ib = 0; ib < Geo->getBlocks().size(); ++ib) {
      auto& popf =
        Lattice->getBlockLat(static_cast<int>(ib))
          .template getField<POP<T, LatSet::q>>();
      const std::array<double, 3> part = frelb_diag::PopChecksumHost(
        popf, Lattice->getBlockLat(static_cast<int>(ib)).getN(), LatSet::q);
      for (int i = 0; i < 3; ++i) total[i] += part[i];
    }
    return total;
  }
};

}  // namespace

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);
  CHECK_EQ(mpi().getSize(), 2);

  const flbgeo::VoxelGeometry<T> record = makeRecord();
  const std::string dir = "ckpt";

  Sim A(record);
  CHECK(A.Geo->getBlockNum() > 0);
  for (int step = 0; step < kPre; ++step) A.step();
  flbckpt::writeWithGeometry(dir, kPre, record, *A.Lattice, *A.FlagFM);

  CHECK(flbckpt::checkCompatibleWithGeometry(dir, -1, record, *A.Lattice,
                                             *A.FlagFM));
  CHECK(!flbckpt::checkCompatible(dir, -1, *A.Lattice, *A.FlagFM));

  Sim B(record);
  const std::int64_t restored =
    flbckpt::readWithGeometry(dir, -1, record, *B.Lattice, *B.FlagFM);
  CHECK_EQ(restored, std::int64_t(kPre));

  for (int step = 0; step < kPost; ++step) A.step();
  for (int step = 0; step < kPost; ++step) B.step();
  const std::array<double, 3> cA = A.popChecksum();
  const std::array<double, 3> cB = B.popChecksum();
  CHECK_EQ(cA[0], cB[0]);
  CHECK_EQ(cA[1], cB[1]);
  CHECK_EQ(cA[2], cB[2]);

  if (mpi().getRank() == 0) {
    std::cout << "[checkpoint_geometry_mpi] ranks=" << mpi().getSize()
              << " digest=" << flbgeo::geometryDigest(record)
              << " sum=" << cA[0] << " sumsq=" << cA[1] << std::endl;
  }

  return flbtest::summary("block/checkpoint_geometry_mpi");
}
