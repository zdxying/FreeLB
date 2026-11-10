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

// voxel_mask.cpp
// B1a bridge: flbgeo::writeVoxelFlags turns a voxel record into the active mask
// of a 3D block geometry.  The test rebuilds the shared helical-duct geometry
// through the existing helper flag path (feature/geo_stl), writes the mask from
// a VoxelGeometry on the same lattice, and requires both masks to agree cell by
// cell, overlap halo included.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

#include "freelb.h"
#include "freelb.hh"

#include "../../test_util.h"

using T = FLOAT;

int Thread_Num = 1;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t InsideFlag = std::uint8_t(2);

}  // namespace

int main() {
  StlReader<T> reader("../geo_stl/stl/helical_duct.stl", T(1), T(1), 2);
  CHECK(reader.getMesh().triangleSize() > 0);

  // Reference path: the helper flag field used by feature/geo_stl.
  BlockGeometryHelper3D<T> helper(reader, 8, 1, 1, true);
  helper.getFlagField().Init(reader, 1, InsideFlag, VoidFlag);
  helper.TagBlockCells(VoidFlag);
  helper.CreateBlocks(true, false);
  helper.RemoveUnusedCells(VoidFlag, false);
  helper.LoadBalancing(1);

  BlockGeometry3D<T> blockgeo(helper);
  BlockFieldManager<FLAG, T, 3> helperFM(blockgeo, VoidFlag);
  helperFM.forEach([&](auto& blockField, std::size_t id) {
    blockField.SetField(
      id, helper.getFlagField()[blockField.getBlock().getLoc_t(id)]);
  });

  // Voxel record on the same lattice.  The helper stores its voxel centers at
  // mesh.min - h/2 + loc * h, so the record's lower corner is one voxel below
  // mesh.min and its cell centers land on the helper's centers.  The grid size
  // mirrors GeometryFlagField::Init with ext = 1.
  const T h = reader.getVoxelSize();
  const Vector<T, 3> meshMin = reader.getMin();
  const Vector<T, 3> meshMax = reader.getMax();
  const Vector<T, 3> lo = meshMin - Vector<T, 3>{h};
  const int nx = static_cast<int>(std::ceil((meshMax[0] - meshMin[0]) / h)) + 2;
  const int ny = static_cast<int>(std::ceil((meshMax[1] - meshMin[1]) / h)) + 2;
  const int nz = static_cast<int>(std::ceil((meshMax[2] - meshMin[2]) / h)) + 2;
  const AABB<T, 3> domain(lo,
                          lo + Vector<T, 3>{T(nx) * h, T(ny) * h, T(nz) * h});

  flbgeo::VoxelGeometryBuilder<T> builder(nx, ny, nz, domain, h, T(1.5));
  builder.addStl(reader, 1, 1, flbgeo::VoxelRole::Domain);
  const flbgeo::VoxelGeometry<T> vgeo = builder.build();

  BlockFieldManager<FLAG, T, 3> geoFM(blockgeo, VoidFlag);
  flbgeo::writeVoxelFlags(vgeo, geoFM, InsideFlag, VoidFlag);

  std::size_t cells = 0;
  std::size_t mismatch = 0;
  std::size_t activeHelper = 0;
  std::size_t activeGeo = 0;
  for (std::size_t ib = 0; ib < helperFM.size(); ++ib) {
    const auto& hfield = helperFM.getBlockField(static_cast<int>(ib)).getField(0);
    const auto& gfield = geoFM.getBlockField(static_cast<int>(ib)).getField(0);
    cells += hfield.size();
    for (std::size_t id = 0; id < hfield.size(); ++id) {
      if (hfield[id] == InsideFlag) ++activeHelper;
      if (gfield[id] == InsideFlag) ++activeGeo;
      if (hfield[id] != gfield[id]) ++mismatch;
    }
  }

  std::cout << "[feature/voxel_mask] cells=" << cells
            << " mismatch=" << mismatch << " active_helper=" << activeHelper
            << " active_geo=" << activeGeo
            << " voxel_domain=" << vgeo.countRole(flbgeo::VoxelRole::Domain)
            << " voxel_solid=" << vgeo.countRole(flbgeo::VoxelRole::Solid)
            << std::endl;

  CHECK_EQ(mismatch, std::size_t{0});
  CHECK_EQ(activeGeo, activeHelper);
  CHECK(activeHelper > 0);

  // End to end: the voxel mask is a valid source for the active-cell map.
  blockgeo.InitVoxelMap(geoFM, VoidFlag);
  CHECK_EQ(blockgeo.getTotalActiveCellNum(geoFM, VoidFlag), activeHelper);

  return flbtest::summary("feature/voxel_mask");
}
