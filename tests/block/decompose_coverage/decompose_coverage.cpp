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

// decompose_coverage.cpp
// B3 coverage contract: every active cell of the geometry record is owned by
// exactly one block, the block table is reproducible from the geometry alone,
// and the same table is seen on every MPI rank.  The split pass runs on the
// helper before LoadBalancing(), so the table itself must not depend on the
// rank count.

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "freelb.h"
#include "freelb.hh"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

int Thread_Num = 1;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t ActiveFlag = std::uint8_t(2);
constexpr int N = 16;
constexpr int MaxBlocks = 24;

std::uint64_t hashLayout(const std::vector<BlockLayoutEntry>& layout) {
  std::uint64_t h = 1469598103934665603ull;
  auto mix = [&](std::uint64_t v) {
    h ^= v;
    h *= 1099511628211ull;
  };
  for (const BlockLayoutEntry& entry : layout) {
    mix(static_cast<std::uint64_t>(entry.BlockId));
    mix(static_cast<std::uint64_t>(entry.Overlap));
    mix(static_cast<std::uint64_t>(entry.Level));
    for (int d = 0; d < 3; ++d) {
      mix(static_cast<std::uint64_t>(entry.IndexMin[d]) +
          0x9e3779b97f4a7c15ull);
      mix(static_cast<std::uint64_t>(entry.Mesh[d]));
    }
  }
  return h;
}

}  // namespace

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);

  const AABB<T, 3> domain(Vector<T, 3>{}, Vector<T, 3>{T(N), T(N), T(N)});
  flbgeo::VoxelGeometryBuilder<T> builder(N, N, N, domain, T(1), T(1.5));
  builder.addBox(domain, 1, 1, flbgeo::VoxelRole::Domain);
  builder.addBox(AABB<T, 3>(Vector<T, 3>{T(4), T(4), T(4)},
                            Vector<T, 3>{T(12), T(12), T(12)}),
                 2, 2, flbgeo::VoxelRole::Solid);
  const flbgeo::VoxelGeometry<T> record = builder.build();
  const std::size_t fluid = record.countRole(flbgeo::VoxelRole::Domain);
  CHECK(fluid > 0);
  CHECK(fluid < record.size());

  auto prepare = [&](BlockGeometryHelper3D<T>& helper) {
    helper.getFlagField().Init(record, ActiveFlag, VoidFlag);
    helper.CreateBlocks(2, 2, 2);
    helper.FillAwareDecomposition(MaxBlocks, T(0), false);
    helper.ValidateActiveCoverage(VoidFlag, false);
  };

  BlockGeometryHelper3D<T> helper(N, N, N, domain, T(1), 8, 1, 0, false);
  prepare(helper);
  const std::vector<BlockLayoutEntry> layout = helper.GetBlockLayout();
  CHECK(layout.size() >= 8);
  CHECK(layout.size() <= static_cast<std::size_t>(MaxBlocks));
  for (std::size_t i = 0; i < layout.size(); ++i) {
    CHECK_EQ(layout[i].BlockId, static_cast<int>(i));
    CHECK_EQ(layout[i].Overlap, 1);
  }

  // The table is a function of the geometry and the limits, not of the order
  // the split search happened to find its candidates.
  BlockGeometryHelper3D<T> again(N, N, N, domain, T(1), 8, 1, 0, false);
  prepare(again);
  CHECK(layout == again.GetBlockLayout());
  CHECK_EQ(hashLayout(layout), hashLayout(again.GetBlockLayout()));

  helper.LoadBalancing(mpi().getSize());
  BlockGeometry3D<T> geo(helper);
  CHECK(geo.getBlockNum() >= 1);

  BlockFieldManager<FLAG, T, 3> FlagFM(geo, VoidFlag);
  flbgeo::writeVoxelFlags(record, FlagFM, ActiveFlag, VoidFlag);
  geo.InitVoxelMap(FlagFM, VoidFlag);

  const std::size_t owned = geo.getTotalOwnedCellNum(FlagFM, VoidFlag);
  const std::size_t active = geo.getTotalActiveCellNum(FlagFM, VoidFlag);
  CHECK_EQ(owned, fluid);
  CHECK(active >= owned);

  // Every rank builds the same block table; only the rank assignment differs.
  const std::uint64_t localHash = hashLayout(layout);
  unsigned long long hashMin = 0;
  unsigned long long hashMax = 0;
#ifdef MPI_ENABLED
  MPI_Allreduce(&localHash, &hashMin, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(&localHash, &hashMax, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                MPI_COMM_WORLD);
#else
  hashMin = localHash;
  hashMax = localHash;
#endif
  CHECK_EQ(hashMin, hashMax);

  if (mpi().getRank() == 0) {
    std::cout << "[block/decompose_coverage] rank=" << mpi().getSize()
              << " blocks=" << layout.size() << " fluid=" << fluid
              << " owned=" << owned << " active=" << active
              << " layout_hash=" << localHash << std::endl;
  }

  return flbtest::summary("block/decompose_coverage");
}
