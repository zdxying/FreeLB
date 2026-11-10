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

// geo_bulk.cpp
// B3/B5 bulk fixture: a 32^3 domain with a 24^3 solid block.  The raw fluid
// fraction is 0.578; the test also runs the record through writeVoxelFlags and
// InitVoxelMap on a one-block baseline and checks the owned fill budget.

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "freelb.h"
#include "freelb.hh"

#include "../../test_util.h"

using T = FLOAT;

int Thread_Num = 1;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t BulkFlag = std::uint8_t(2);

struct DecompResult {
  std::size_t blocks{};
  std::size_t allocated{};
  std::size_t owned{};
  double fill{};
  std::vector<BlockLayoutEntry> layout{};
};

// Decompose the record into a uniform base and let FillAwareDecomposition
// tighten it.  The helper is rebuilt from the record on every call so the
// layout is a function of the geometry and the limits only.
DecompResult measureDecomp(const flbgeo::VoxelGeometry<T>& geometry, int bx,
                           int by, int bz, int maxBlocks, T target,
                           bool verbose = false) {
  BlockGeometryHelper3D<T> helper(
    geometry.getNx(), geometry.getNy(), geometry.getNz(), geometry.getDomain(),
    geometry.getVoxelSize(), 8, 1, 0, false);
  helper.getFlagField().Init(geometry, BulkFlag, VoidFlag);
  helper.CreateBlocks(bx, by, bz);
  helper.FillAwareDecomposition(maxBlocks, target, verbose);
  helper.LoadBalancing(1);

  BlockGeometry3D<T> blockgeo(helper);
  BlockFieldManager<FLAG, T, 3> flagFM(blockgeo, VoidFlag);
  flbgeo::writeVoxelFlags(geometry, flagFM, BulkFlag, VoidFlag);
  blockgeo.InitVoxelMap(flagFM, VoidFlag);

  DecompResult result;
  result.blocks = blockgeo.getBlockNum();
  result.allocated = blockgeo.getTotalCellNum();
  result.owned = blockgeo.getTotalOwnedCellNum(flagFM, VoidFlag);
  result.fill = result.allocated ? static_cast<double>(result.owned) / result.allocated : 0.;
  result.layout = helper.GetBlockLayout();
  return result;
}

}  // namespace

int main() {
  const flbgeo::GeometryFile<T> data = flbgeo::readGeometryFile<T>("geo_bulk.flbgeo");
  const flbgeo::VoxelGeometry<T>& geometry = data.Geometry;

  CHECK_EQ(geometry.getNx(), 32);
  CHECK_EQ(geometry.getNy(), 32);
  CHECK_EQ(geometry.getNz(), 32);

  const std::size_t fluid = geometry.countRole(flbgeo::VoxelRole::Domain);
  const std::size_t solid = geometry.countRole(flbgeo::VoxelRole::Solid);
  const std::size_t center = geometry.findCellId(Vector<T, 3>{T(16), T(16), T(16)});
  CHECK_EQ(fluid + solid, geometry.size());
  CHECK_EQ(solid, geometry.countEntity(2));
  CHECK_EQ(geometry.getEntity(center), std::uint16_t{2});
  CHECK_EQ(geometry.getMaterial(center), std::uint16_t{2});
  CHECK(geometry.isSolid(center));
  CHECK(geometry.isDomain(geometry.findCellId(Vector<T, 3>{T(1.5), T(1.5), T(1.5)})));

  const double rawFill = static_cast<double>(fluid) / geometry.size();
  CHECK(rawFill > 0.55);
  CHECK(rawFill < 0.60);

  // One-block baseline: this tests the record -> mask -> voxel-map path; B3
  // owns the multi-block decomposition and its budget.
  BlockGeometry3D<T> blockgeo(32, 32, 32, 1, geometry.getDomain(),
                              geometry.getVoxelSize(), 1);
  CHECK_EQ(blockgeo.getBlockNum(), 1);
  BlockFieldManager<FLAG, T, 3> flagFM(blockgeo, VoidFlag);
  flbgeo::writeVoxelFlags(geometry, flagFM, BulkFlag, VoidFlag);
  blockgeo.InitVoxelMap(flagFM, VoidFlag);

  const std::size_t active = blockgeo.getTotalActiveCellNum(flagFM, VoidFlag);
  const std::size_t owned = blockgeo.getTotalOwnedCellNum(flagFM, VoidFlag);
  const std::size_t allocated = blockgeo.getTotalCellNum();
  const double ownedFill = static_cast<double>(owned) / allocated;

  std::cout << "[feature/geo_bulk] domain=" << fluid << " solid=" << solid
            << " allocated=" << allocated << " owned=" << owned
            << " active=" << active << " raw_fill=" << rawFill
            << " owned_fill=" << ownedFill
            << " blocks=" << blockgeo.getBlockNum() << std::endl;

  CHECK_EQ(active, fluid);
  CHECK(ownedFill > 0.45);
  CHECK(ownedFill < 0.60);

  // The record -> helper flag field bridge is the decomposition input; it must
  // see every domain cell exactly once.
  {
    BlockGeometryHelper3D<T> bridge(32, 32, 32, geometry.getDomain(),
                                    geometry.getVoxelSize(), 8, 1, 0, false);
    bridge.getFlagField().Init(geometry, BulkFlag, VoidFlag);
    CHECK_EQ(CountActiveCells(bridge.getBaseBlock(), bridge.getFlagField(),
                              VoidFlag),
             fluid);
  }

  // B3: uniform base plus fill-aware split.  Starting at 8 blocks gives the
  // split a table to work from; max_blocks bounds the result.
  for (int maxBlocks : {8, 16, 24, 32}) {
    const DecompResult r =
      measureDecomp(geometry, 2, 2, 2, maxBlocks, T(0), false);
    std::cout << "[feature/geo_bulk] max_blocks=" << maxBlocks
              << " blocks=" << r.blocks << " allocated=" << r.allocated
              << " owned=" << r.owned << " owned_fill=" << r.fill
              << std::endl;
    CHECK_EQ(r.owned, fluid);
    CHECK(r.blocks <= static_cast<std::size_t>(maxBlocks));
    // B3 bulk budget: the 16-block result must clear the target floor, and
    // using the cap for the slab pass must not lose footprint.
    if (maxBlocks == 16) {
      CHECK(r.fill > 0.45);
      CHECK(r.fill < 0.60);
    }
    if (maxBlocks == 24) {
      CHECK(r.fill > 0.48);
    }
  }

  // The layout is a function of the geometry and the limits, not of the call
  // order: two independent helpers must produce the same table.
  {
    const DecompResult a = measureDecomp(geometry, 2, 2, 2, 16, T(0));
    const DecompResult b = measureDecomp(geometry, 2, 2, 2, 16, T(0));
    CHECK(a.layout == b.layout);
  }

  return flbtest::summary("feature/geo_bulk");
}
