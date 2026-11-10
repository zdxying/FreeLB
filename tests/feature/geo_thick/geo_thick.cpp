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

// geo_thick.cpp
// B3/B4 thick-wall, multi-entity fixture: two 4-cell-thick plates in a 32^3
// fluid domain.  The test checks the wall thickness, both entity/material
// pairs, the surface band, and the record -> mask -> voxel-map path.

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

std::size_t cell(const flbgeo::VoxelGeometry<T>& geo, int x, int y, int z) {
  return geo.getIndex(Vector<int, 3>{x, y, z});
}

struct DecompResult {
  std::size_t blocks{};
  std::size_t allocated{};
  std::size_t owned{};
  double fill{};
};

// Start from a uniform table and let the fill-aware split tighten it.  Rebuilt
// from the record on every call so the table does not depend on call order.
DecompResult measureDecomp(const flbgeo::VoxelGeometry<T>& geometry, int bx,
                           int by, int bz, int maxBlocks) {
  BlockGeometryHelper3D<T> helper(
    geometry.getNx(), geometry.getNy(), geometry.getNz(), geometry.getDomain(),
    geometry.getVoxelSize(), 8, 1, 0, false);
  helper.getFlagField().Init(geometry, BulkFlag, VoidFlag);
  helper.CreateBlocks(bx, by, bz);
  helper.FillAwareDecomposition(maxBlocks, T(0), false);
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
  return result;
}

}  // namespace

int main() {
  const flbgeo::GeometryFile<T> data = flbgeo::readGeometryFile<T>("geo_thick.flbgeo");
  const flbgeo::VoxelGeometry<T>& geometry = data.Geometry;

  CHECK_EQ(geometry.getNx(), 32);
  CHECK_EQ(geometry.getNy(), 32);
  CHECK_EQ(geometry.getNz(), 32);

  const std::size_t fluid = geometry.countRole(flbgeo::VoxelRole::Domain);
  const std::size_t solid = geometry.countRole(flbgeo::VoxelRole::Solid);
  CHECK_EQ(fluid + solid, geometry.size());
  CHECK_EQ(solid, geometry.countEntity(2) + geometry.countEntity(3));
  CHECK_EQ(geometry.countEntity(2), std::size_t{4 * 16 * 16});
  CHECK_EQ(geometry.countEntity(3), std::size_t{4 * 16 * 16});

  // Each plate is four cells thick along x; the 12-cell gap is fluid.
  for (int x = 8; x < 12; ++x) {
    const std::size_t id = cell(geometry, x, 16, 16);
    CHECK(geometry.isSolid(id));
    CHECK_EQ(geometry.getEntity(id), std::uint16_t{2});
    CHECK_EQ(geometry.getMaterial(id), std::uint16_t{2});
  }
  for (int x = 12; x < 20; ++x) {
    CHECK(geometry.isDomain(cell(geometry, x, 16, 16)));
  }
  for (int x = 20; x < 24; ++x) {
    const std::size_t id = cell(geometry, x, 16, 16);
    CHECK(geometry.isSolid(id));
    CHECK_EQ(geometry.getEntity(id), std::uint16_t{3});
    CHECK_EQ(geometry.getMaterial(id), std::uint16_t{3});
  }
  CHECK(geometry.isDomain(cell(geometry, 7, 16, 16)));
  CHECK(geometry.isDomain(cell(geometry, 24, 16, 16)));

  // The band covers the fluid cells adjacent to both plates.
  CHECK(!geometry.getBandCells().empty());
  const std::size_t nearWall = cell(geometry, 7, 16, 16);
  CHECK(geometry.hasSurfaceSample(nearWall));
  CHECK_EQ(geometry.getSurfaceSample(nearWall).Entity, std::uint16_t{2});

  BlockGeometry3D<T> blockgeo(32, 32, 32, 1, geometry.getDomain(),
                              geometry.getVoxelSize(), 1);
  BlockFieldManager<FLAG, T, 3> flagFM(blockgeo, VoidFlag);
  flbgeo::writeVoxelFlags(geometry, flagFM, BulkFlag, VoidFlag);
  blockgeo.InitVoxelMap(flagFM, VoidFlag);

  const std::size_t active = blockgeo.getTotalActiveCellNum(flagFM, VoidFlag);
  const std::size_t owned = blockgeo.getTotalOwnedCellNum(flagFM, VoidFlag);
  const std::size_t allocated = blockgeo.getTotalCellNum();
  const double ownedFill = static_cast<double>(owned) / allocated;

  std::cout << "[feature/geo_thick] domain=" << fluid << " solid=" << solid
            << " entity2=" << geometry.countEntity(2)
            << " entity3=" << geometry.countEntity(3)
            << " allocated=" << allocated << " owned=" << owned
            << " active=" << active << " owned_fill=" << ownedFill
            << " blocks=" << blockgeo.getBlockNum() << std::endl;

  CHECK_EQ(active, fluid);
  CHECK(owned > 0);
  CHECK(ownedFill > 0.70);

  // Thick walls leave more fluid between blocks than the bulk shell, so the
  // split pass should keep the fill at or above the 8-block uniform baseline
  // (0.658 in the plan's baseline table) while the cap grows.
  for (int maxBlocks : {8, 16, 24}) {
    const DecompResult r = measureDecomp(geometry, 2, 2, 2, maxBlocks);
    std::cout << "[feature/geo_thick] max_blocks=" << maxBlocks
              << " blocks=" << r.blocks << " allocated=" << r.allocated
              << " owned=" << r.owned << " owned_fill=" << r.fill
              << std::endl;
    CHECK_EQ(r.owned, fluid);
    CHECK(r.blocks <= static_cast<std::size_t>(maxBlocks));
    CHECK(r.fill >= 0.65);
  }

  return flbtest::summary("feature/geo_thick");
}
