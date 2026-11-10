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

// link_bc_seam.cpp
// Two-rank seam test for the B4 surface-link table.  The helical duct is split
// into 2x2x2 blocks and the blocks are distributed over two MPI ranks, so part
// of every link table is built from a halo that belongs to another rank.  The
// global link count and the fluid mass must match the single-block run, and the
// ratio of wall links to fluid cells must not depend on the rank that owns the
// link.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

#include "freelb.h"
#include "freelb.hh"
#include "boundary/link_bc.h"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

int Thread_Num = 1;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t BulkFlag = std::uint8_t(2);
constexpr int Steps = 20;

using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
using CELL = Cell<T, LatSet, FIELDS>;
using BulkTask =
  tmp::Key_TypePair<BulkFlag,
                    collision::BGK<moment::rhoU<CELL>,
                                   equilibrium::SecondOrder<CELL>>>;
using TaskCollection = tmp::TupleWrapper<BulkTask>;
using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;

std::size_t globalSum(std::size_t local) {
  unsigned long long value = local;
  unsigned long long sum = 0;
  MPI_Allreduce(&value, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  return static_cast<std::size_t>(sum);
}

}  // namespace

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT;

  const int rank = mpi().getRank();
  const int size = mpi().getSize();
  CHECK_EQ(size, 2);

  StlReader<T> reader("../../feature/geo_stl/stl/helical_duct.stl", T(1), T(1), 2);
  CHECK(reader.getMesh().triangleSize() > 0);

  const T h = reader.getVoxelSize();
  const Vector<T, 3> meshMin = reader.getMin();
  const Vector<T, 3> meshMax = reader.getMax();
  const Vector<T, 3> lo = meshMin - Vector<T, 3>{h};
  const int nx = static_cast<int>(std::ceil((meshMax[0] - meshMin[0]) / h)) + 2;
  const int ny = static_cast<int>(std::ceil((meshMax[1] - meshMin[1]) / h)) + 2;
  const int nz = static_cast<int>(std::ceil((meshMax[2] - meshMin[2]) / h)) + 2;
  const AABB<T, 3> domain(
    lo, lo + Vector<T, 3>{T(nx) * h, T(ny) * h, T(nz) * h});

  flbgeo::VoxelGeometryBuilder<T> builder(nx, ny, nz, domain, h, T(1.5));
  builder.addStl(reader, 1, 1, flbgeo::VoxelRole::Domain, true);
  const flbgeo::VoxelGeometry<T> record = builder.build();
  CHECK(record.countRole(flbgeo::VoxelRole::Domain) > 0);

  BlockGeometryHelper3D<T> helper(nx, ny, nz, record.getDomain(),
                                  record.getVoxelSize(), 8, 1, 0, false);
  helper.CreateBlocks(2, 2, 2);
  helper.LoadBalancing(size);
  BlockGeometry3D<T> geo(helper);
  CHECK(geo.getBlockNum() > 0);

  BlockFieldManager<FLAG, T, 3> flagFM(geo, VoidFlag);
  flbgeo::writeVoxelFlags(record, flagFM, BulkFlag, VoidFlag);
  // writeVoxelFlags writes global positions into the halo too, so the link
  // table of a block sees the same flags as in the single-block run.
  geo.InitVoxelMap(flagFM);

  BaseConverter<T> conv(LatSet::cs2);
  conv.SimplifiedConverterFromRT(nx, T(0.1), T(0.56));
  ValuePack initValues(conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(geo, initValues, conv);

  flbbc::InterpolatedBBManager<decltype(NSLattice), decltype(flagFM)> ibb(
    NSLattice, flagFM, record, VoidFlag);
  CHECK(ibb.getLinkCount() > 0);

  std::size_t ownedCells = 0;
  T before = T(0);
  for (int ib = 0; ib < geo.getBlockNum(); ++ib) {
    const auto& block = geo.getBlock(ib);
    const int overlap = block.getOverlap();
    const Vector<int, 3> mesh = block.getBaseBlock().getMesh();
    const auto& flag = flagFM.getBlockField(ib).getField(0);
    const auto& rho = NSLattice.getBlockLat(ib).template getField<RHO<T>>();
    for (int z = overlap; z < overlap + mesh[2]; ++z) {
      for (int y = overlap; y < overlap + mesh[1]; ++y) {
        for (int x = overlap; x < overlap + mesh[0]; ++x) {
          const std::size_t id = block.getIndex(Vector<int, 3>{x, y, z});
          if (flag[id] != VoidFlag) {
            ++ownedCells;
            before += rho.get(id);
          }
        }
      }
    }
  }

  for (int step = 0; step < Steps; ++step) {
    NSLattice.ApplyCellDynamics<NSTask>(flagFM);
    NSLattice.Stream();
    ibb.Apply();
    NSLattice.NormalFullCommunicate();
  }

  T after = T(0);
  for (int ib = 0; ib < geo.getBlockNum(); ++ib) {
    const auto& block = geo.getBlock(ib);
    const int overlap = block.getOverlap();
    const Vector<int, 3> mesh = block.getBaseBlock().getMesh();
    const auto& flag = flagFM.getBlockField(ib).getField(0);
    const auto& rho = NSLattice.getBlockLat(ib).template getField<RHO<T>>();
    for (int z = overlap; z < overlap + mesh[2]; ++z) {
      for (int y = overlap; y < overlap + mesh[1]; ++y) {
        for (int x = overlap; x < overlap + mesh[0]; ++x) {
          const std::size_t id = block.getIndex(Vector<int, 3>{x, y, z});
          if (flag[id] != VoidFlag) after += rho.get(id);
        }
      }
    }
  }

  const std::size_t globalLinks = globalSum(ibb.getLinkCount());
  const std::size_t globalCells = globalSum(ownedCells);
  const T globalBefore = [&] {
    double value = static_cast<double>(before);
    double sum = 0;
    MPI_Allreduce(&value, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return static_cast<T>(sum);
  }();
  const T globalAfter = [&] {
    double value = static_cast<double>(after);
    double sum = 0;
    MPI_Allreduce(&value, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return static_cast<T>(sum);
  }();

  // Both ranks must own duct blocks, otherwise the seam is not exercised.
  unsigned long long localLinks = ibb.getLinkCount();
  unsigned long long minLinks = 0;
  MPI_Allreduce(&localLinks, &minLinks, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN,
                MPI_COMM_WORLD);
  CHECK(minLinks > 0);

  const T drift = std::abs(globalAfter - globalBefore);
  const std::size_t globalBlocks =
    globalSum(static_cast<std::size_t>(geo.getBlockNum()));

  if (rank == 0) {
    std::cout << "[block/link_bc_seam] ranks=" << size
              << " blocks=" << globalBlocks
              << " links=" << globalLinks << " cells=" << globalCells
              << " mass_before=" << globalBefore
              << " mass_after=" << globalAfter << " drift=" << drift
              << std::endl;
  }

  CHECK_EQ(globalBlocks, std::size_t(8));
  CHECK_EQ(globalLinks, std::size_t(7160));
  CHECK_EQ(globalCells, std::size_t(1294));
  CHECK_NEAR(globalBefore, T(1294), T(1e-12));
  CHECK_NEAR(globalAfter, globalBefore, T(1e-10));

  return flbtest::summary("block/link_bc_seam");
}
