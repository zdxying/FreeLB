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

// geo_stl.cpp
// P0 fixture: a closed, warped, square helical duct (see stl/helical_duct.stl).
// The test does not use any private core asset.  It checks inside/outside
// polarity, then measures the active-cell fill ratio at two voxel sizes.

#include <cmath>
#include <cstdint>

#include "freelb.h"
#include "freelb.hh"
#include "utils/memory_usage.h"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

int Thread_Num = 1;

namespace {

struct FillResult {
  std::size_t active{};
  std::size_t owned{};
  std::size_t allocated{};
  double fill{};
};

// Build a block geometry from the STL at one voxel size.  The helper creates
// block cells, tags those that touch the solid surface, merges them into
// blocks, and trims fully void layers.  The resulting FLAG field is then
// copied into a BlockGeometry3D field so the P0 metrics can be exercised.
FillResult measure(const char* path, T voxelSize, int maxBlocks = 0,
                   bool print = true, int rcbBlocks = 0) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t InsideFlag = std::uint8_t(2);

  StlReader<T> reader(path, voxelSize, T(1), 2);
  BlockGeometryHelper3D<T> helper(reader, 8, 1, 1, true);
  helper.getFlagField().Init(reader, 1, InsideFlag, VoidFlag);
  helper.TagBlockCells(VoidFlag);
  if (rcbBlocks > 0) {
    helper.RCBOptimization(rcbBlocks);
  } else {
    helper.CreateBlocks(true, false);
    helper.RemoveUnusedCells(VoidFlag, false);
    if (maxBlocks > 0) helper.FillAwareDecomposition(maxBlocks, T(0), false);
  }
  helper.LoadBalancing(1);

  BlockGeometry3D<T> geo(helper);
  BlockFieldManager<FLAG, T, 3> FlagFM(geo, VoidFlag);
  FlagFM.forEach([&](auto& blockfield, std::size_t id) {
    blockfield.SetField(id, helper.getFlagField()[blockfield.getBlock().getLoc_t(id)]);
  });

  geo.InitVoxelMap(FlagFM);
  if (print) geo.PrintFillRatio(FlagFM, VoidFlag);

  FillResult result;
  result.active = geo.getTotalActiveCellNum(FlagFM, VoidFlag);
  result.owned = geo.getTotalOwnedCellNum(FlagFM, VoidFlag);
  result.allocated = geo.getTotalCellNum();
  result.fill = geo.getOwnedFillRatio(FlagFM, VoidFlag);
  return result;
}

}  // namespace

int main() {
  constexpr T pi = T(3.14159265358979323846);
  constexpr T R = T(10);
  constexpr T pitch = T(30);
  constexpr T turns = T(1.25);

  const char* path = "./stl/helical_duct.stl";
  StlReader<T> reader(path, T(1), T(1), 2);
  CHECK(reader.getMesh().triangleSize() > 0);
  CHECK(reader.getMesh().getMax_Min()[0] > T(20));
  CHECK(reader.getMesh().getMax_Min()[2] > T(30));

  auto centerline = [&](T t) {
    const T a = T(2) * pi * turns * t;
    return Vector<T, 3>{R * std::cos(a), R * std::sin(a), pitch * t};
  };
  auto isInside = [&](const Vector<T, 3>& p) {
    Octree<T>* node = reader.getTree()->find(p);
    return node != nullptr && node->isLeaf() && node->getInside();
  };

  for (T t : {T(0.05), T(0.25), T(0.50), T(0.75), T(0.95)}) {
    CHECK(isInside(centerline(t)));
  }
  const Vector<T, 3> outside =
    reader.getMesh().getMin() - Vector<T, 3>{T(2), T(2), T(2)};
  CHECK(!isInside(outside));

  const FillResult coarse = measure(path, T(1));
  CHECK(coarse.active > 0);
  CHECK(coarse.allocated > coarse.active);
  CHECK(coarse.active < coarse.allocated / 2);
  CHECK(coarse.fill > 0.01);
  CHECK(coarse.fill < 0.30);

  const FillResult fine = measure(path, T(0.5));
  CHECK(fine.active > coarse.active * 4);
  CHECK(fine.active < coarse.active * 12);
  CHECK(fine.allocated > fine.active);
  CHECK(fine.fill > 0.01);
  CHECK(fine.fill < 0.40);
  CHECK(std::abs(fine.fill - coarse.fill) < 0.20);

  // B3 diagnostics: the helical duct is thin and curved, so this is the hard
  // end of the fill budget.  The split pass only adds blocks; it cannot merge
  // the existing 11/22-block table, so the cap decides how far it can go.
  for (int maxBlocks : {16, 24, 32}) {
    const FillResult r = measure(path, T(0.5), maxBlocks, false);
    std::cout << "[feature/geo_stl] h=0.5 max_blocks=" << maxBlocks
              << " active=" << r.active << " owned=" << r.owned
              << " allocated=" << r.allocated << " owned_fill=" << r.fill
              << std::endl;
    CHECK_EQ(r.owned, fine.owned);
    CHECK(r.fill >= fine.fill - 1e-12);
  }
  for (int maxBlocks : {16, 24}) {
    const FillResult r = measure(path, T(1), maxBlocks, false);
    std::cout << "[feature/geo_stl] h=1.0 max_blocks=" << maxBlocks
              << " active=" << r.active << " owned=" << r.owned
              << " allocated=" << r.allocated << " owned_fill=" << r.fill
              << std::endl;
    CHECK_EQ(r.owned, coarse.owned);
    CHECK(r.fill >= coarse.fill - 1e-12);
  }

  // Alternative table shape at the same budget: RCB on the same flag field.
  // This tells us whether the remaining gap is the split topology or the
  // dense AABB storage itself.
  for (int n : {16, 32}) {
    const FillResult r = measure(path, T(0.5), 0, false, n);
    std::cout << "[feature/geo_stl] h=0.5 rcb=" << n
              << " active=" << r.active << " owned=" << r.owned
              << " allocated=" << r.allocated << " owned_fill=" << r.fill
              << std::endl;
    CHECK_EQ(r.owned, fine.owned);
  }
  for (int n : {8, 16}) {
    const FillResult r = measure(path, T(1), 0, false, n);
    std::cout << "[feature/geo_stl] h=1.0 rcb=" << n
              << " active=" << r.active << " owned=" << r.owned
              << " allocated=" << r.allocated << " owned_fill=" << r.fill
              << std::endl;
    CHECK_EQ(r.owned, coarse.owned);
  }

  flbmem::printMemoryUsage("feature/geo_stl");
#if defined(__linux__)
  CHECK(flbmem::peakRssBytes() > 0);
#endif

  return flbtest::summary("feature/geo_stl");
}
