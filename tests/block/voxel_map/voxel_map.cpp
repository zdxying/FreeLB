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
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// voxel_map.cpp
// P0: compile and exercise the compact iteration path, fill-ratio diagnostics,
// and the empty-block semantics of VoxelMap.

#include <chrono>
#include <cstdint>
#include <thread>
#include <utility>

#include "freelb.h"
#include "freelb.hh"
#include "data_struct/voxel_map.h"
#include "utils/memory_usage.h"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

int Thread_Num = 1;

struct CountDynamics {
  static inline std::size_t count = 0;

  template <typename CELL>
  static void Execute(std::uint8_t, CELL&) {
    ++count;
  }
};

int main() {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t ActiveFlag = std::uint8_t(2);
  constexpr int Nx = 12;
  constexpr int Ny = 12;
  constexpr int Nz = 12;

  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;

  BaseConverter<T> BaseConv(LatSet::cs2);
  AABB<T, 3> domain(Vector<T, 3>{}, Vector<T, 3>{T(Nx), T(Ny), T(Nz)});
  BlockGeometry3D<T> Geo(Nx, Ny, Nz, 1, domain, T(1), 1);
  BlockFieldManager<FLAG, T, 3> FlagFM(Geo, VoidFlag);

  Block3D<T>& block = Geo.getBlock(0);
  auto& flags = FlagFM.getBlockField(0).getField(0);
  const int overlap = block.getOverlap();
  const std::size_t allocated = block.getN();

  auto resetFlags = [&]() {
    for (std::size_t id = 0; id < flags.size(); ++id) flags[id] = VoidFlag;
  };

  ValuePack initValues(BaseConv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(Geo, initValues, BaseConv);

  // --- thin plate: one plane in the extended block -------------------------
  resetFlags();
  std::size_t expected = 0;
  for (std::size_t id = 0; id < flags.size(); ++id) {
    const Vector<int, 3> loc = block.getLoc(id);
    if (loc[0] == overlap + Nx / 2) {
      flags[id] = ActiveFlag;
      ++expected;
    }
  }
  const std::uint64_t generationBefore = Geo.getGeometryGeneration();
  Geo.InitVoxelMap(FlagFM);
  CHECK_EQ(Geo.getGeometryGeneration(), generationBefore + 1);
  CHECK_EQ(block.getVoxelMap().getVoxNum(), expected);
  CHECK_EQ(block.getVoxNum(), expected);
  CHECK_EQ(Geo.getTotalActiveCellNum(FlagFM, VoidFlag), expected);
  CHECK_EQ(Geo.getTotalOwnedCellNum(FlagFM, VoidFlag), std::size_t(Ny * Nz));
  CHECK_EQ(block.getVoxelMap().getOwnedNum(), std::size_t(Ny * Nz));
  CHECK_NEAR(Geo.getFillRatio(FlagFM, VoidFlag),
             static_cast<double>(expected) / static_cast<double>(allocated), 1e-12);
  CHECK_NEAR(Geo.getOwnedFillRatio(FlagFM, VoidFlag),
             static_cast<double>(Ny * Nz) / static_cast<double>(allocated), 1e-12);
  bool sawActive = false;
  bool sawVoid = false;
  for (std::size_t id = 0; id < allocated; ++id) {
    if (block.getVoxelMap().isActive(id)) {
      CHECK_EQ(block.getVoxelMap()[block.getVoxelMap().aabbToSlot(id)], id);
      sawActive = true;
    } else {
      sawVoid = true;
    }
  }
  CHECK(sawActive);
  CHECK(sawVoid);
  for (std::size_t i = 0; i < expected; ++i) {
    CHECK_EQ(block.getLoc(block.getVoxelMap()[i])[0], overlap + Nx / 2);
  }
  CountDynamics::count = 0;
  NSLattice.ApplyCellDynamics<CountDynamics>(FlagFM);
  CHECK_EQ(CountDynamics::count, expected);
  Geo.PrintFillRatio(FlagFM, VoidFlag);

  // --- diagonal sheet: exercises non-contiguous active ids -----------------
  resetFlags();
  expected = 0;
  for (std::size_t id = 0; id < flags.size(); ++id) {
    const Vector<int, 3> loc = block.getLoc(id);
    if (loc[0] + loc[1] == (Nx + 2 * overlap - 1)) {
      flags[id] = ActiveFlag;
      ++expected;
    }
  }
  Geo.InitVoxelMap(FlagFM);
  CHECK_EQ(block.getVoxNum(), expected);
  CHECK_EQ(Geo.getTotalActiveCellNum(FlagFM, VoidFlag), expected);
  for (std::size_t i = 0; i < expected; ++i) {
    const Vector<int, 3> loc = block.getLoc(block.getVoxelMap()[i]);
    CHECK_EQ(loc[0] + loc[1], Nx + 2 * overlap - 1);
  }
  CountDynamics::count = 0;
  NSLattice.ApplyCellDynamics<CountDynamics>(FlagFM);
  CHECK_EQ(CountDynamics::count, expected);

  // --- fully void block: built but empty, not a fallback to getN() ---------
  resetFlags();
  Geo.InitVoxelMap(FlagFM);
  CHECK(block.getVoxelMap().isInitialized());
  CHECK(block.getVoxelMap().empty());
  CHECK_EQ(block.getVoxNum(), 0);
  CHECK_EQ(Geo.getTotalActiveCellNum(FlagFM, VoidFlag), 0);
  CHECK_EQ(Geo.getTotalOwnedCellNum(FlagFM, VoidFlag), 0);
  CHECK_EQ(block.getVoxelMap().getOwnedNum(), 0);
  CHECK_NEAR(Geo.getFillRatio(FlagFM, VoidFlag), 0., 1e-12);
  CHECK_NEAR(Geo.getOwnedFillRatio(FlagFM, VoidFlag), 0., 1e-12);
  CountDynamics::count = 0;
  NSLattice.ApplyCellDynamics<CountDynamics>(FlagFM);
  CHECK_EQ(CountDynamics::count, 0);

  // --- dense block: every cell is visited ---------------------------------
  resetFlags();
  for (std::size_t id = 0; id < flags.size(); ++id) flags[id] = ActiveFlag;
  Geo.InitVoxelMap(FlagFM);
  CHECK_EQ(block.getVoxNum(), allocated);
  CHECK_EQ(block.getVoxelMap().getOwnedNum(), std::size_t(Nx * Ny * Nz));
  CHECK_EQ(Geo.getTotalOwnedCellNum(FlagFM, VoidFlag), std::size_t(Nx * Ny * Nz));
  CHECK_NEAR(Geo.getFillRatio(FlagFM, VoidFlag), 1., 1e-12);
  CHECK_NEAR(Geo.getOwnedFillRatio(FlagFM, VoidFlag),
             static_cast<double>(Nx * Ny * Nz) / static_cast<double>(allocated), 1e-12);
  CountDynamics::count = 0;
  NSLattice.ApplyCellDynamics<CountDynamics>(FlagFM);
  CHECK_EQ(CountDynamics::count, allocated);

  // --- 2D twin: fill metrics and the same compact iteration path -----------
  AABB<T, 2> domain2(Vector<T, 2>{}, Vector<T, 2>{T(8), T(8)});
  BlockGeometry2D<T> Geo2(8, 8, 1, domain2, T(1), 1);
  BlockFieldManager<FLAG, T, 2> FlagFM2(Geo2, VoidFlag);
  Block2D<T>& block2 = Geo2.getBlock(0);
  auto& flags2 = FlagFM2.getBlockField(0).getField(0);
  std::size_t expected2 = 0;
  for (std::size_t id = 0; id < flags2.size(); ++id) {
    const Vector<int, 2> loc = block2.getLoc(id);
    if (loc[0] == block2.getOverlap() + 4) {
      flags2[id] = ActiveFlag;
      ++expected2;
    }
  }
  Geo2.InitVoxelMap(FlagFM2);
  CHECK_EQ(block2.getVoxNum(), expected2);
  CHECK_EQ(Geo2.getTotalActiveCellNum(FlagFM2, VoidFlag), expected2);
  CHECK_EQ(block2.getVoxelMap().getOwnedNum(), std::size_t(8));
  CHECK_EQ(Geo2.getTotalOwnedCellNum(FlagFM2, VoidFlag), std::size_t(8));
  CHECK_NEAR(Geo2.getFillRatio(FlagFM2, VoidFlag),
             static_cast<double>(expected2) / static_cast<double>(block2.getN()), 1e-12);
  CHECK_NEAR(Geo2.getOwnedFillRatio(FlagFM2, VoidFlag),
             static_cast<double>(8) / static_cast<double>(block2.getN()), 1e-12);
  for (std::size_t i = 0; i < expected2; ++i) {
    CHECK_EQ(block2.getLoc(block2.getVoxelMap()[i])[0], block2.getOverlap() + 4);
  }

  // --- VoxelMap value semantics and rebuild --------------------------------
  VoxelMap a(3);
  a.getVox_to_AABB()[0] = 7;
  a.getVox_to_AABB()[1] = 9;
  a.getVox_to_AABB()[2] = 11;
  VoxelMap b(a);
  CHECK_EQ(b.getVoxNum(), std::size_t(3));
  CHECK_EQ(b[1], std::size_t(9));
  VoxelMap c(std::move(b));
  CHECK_EQ(c.getVoxNum(), std::size_t(3));
  CHECK_EQ(c[2], std::size_t(11));
  CHECK(!b.isInitialized());
  CHECK_EQ(b.getVoxNum(), std::size_t(0));
  a.Init(1);
  CHECK_EQ(a.getVoxNum(), std::size_t(1));
  CHECK_EQ(a[0], std::size_t(0));
  a.Init(0);
  CHECK(a.isInitialized());
  CHECK(a.empty());

  // --- Init(BlockGeometryHelper&) rebuild bumps the generation -------------
  BlockGeometryHelper3D<T> helper(Nx, Ny, Nz, domain, T(1), 4, 1);
  helper.CreateBlocks(2, 2, 1);
  helper.LoadBalancing(1);
  BlockGeometry3D<T> helperGeo(helper);
  CHECK_EQ(helperGeo.getBlockNum(), 4);
  const std::uint64_t helperGeneration = helperGeo.getGeometryGeneration();
  helperGeo.Init(helper);
  CHECK(helperGeo.getGeometryGeneration() > helperGeneration);
  CHECK_EQ(helperGeo.getBlockNum(), 4);

  // --- diagnostics presence ------------------------------------------------
#if defined(__linux__)
  CHECK(flbmem::peakRssBytes() > 0);
#endif
  flbmem::printMemoryUsage("block/voxel_map");

  Timer timer;
  ++timer;
  std::this_thread::sleep_for(std::chrono::milliseconds(2));
  timer.Print_MainLoopPerformance(Geo.getTotalCellNum(),
                                  Geo.getTotalActiveCellNum(FlagFM, VoidFlag));

  return flbtest::summary("block/voxel_map");
}
