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

// voxel_entity.cpp
// B1a: VoxelGeometryBuilder combines a domain box and the closed helical-duct STL
// into one voxel record, keeps the solid's entity/material under a repeated
// domain box, and overlays a wall patch without changing occupancy.  The STL is
// the shared fixture of feature/geo_stl; no private asset is used.

#include <algorithm>
#include <cmath>
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

struct VoxelCounts {
  std::size_t domain{};
  std::size_t solid{};
  std::size_t outside{};
  std::size_t entity1{};
  std::size_t entity2{};
  std::size_t material1{};
  std::size_t material2{};
};

VoxelCounts tally(const flbgeo::VoxelGeometry<T>& geo) {
  VoxelCounts counts;
  for (std::size_t id = 0; id < geo.size(); ++id) {
    if (geo.isDomain(id)) ++counts.domain;
    if (geo.isSolid(id)) ++counts.solid;
    if (geo.isOutside(id)) ++counts.outside;
    if (geo.getEntity(id) == std::uint16_t{1}) ++counts.entity1;
    if (geo.getEntity(id) == std::uint16_t{2}) ++counts.entity2;
    if (geo.getMaterial(id) == std::uint16_t{1}) ++counts.material1;
    if (geo.getMaterial(id) == std::uint16_t{2}) ++counts.material2;
  }
  return counts;
}

void checkSameCounts(const VoxelCounts& got, const VoxelCounts& expected) {
  CHECK_EQ(got.domain, expected.domain);
  CHECK_EQ(got.solid, expected.solid);
  CHECK_EQ(got.outside, expected.outside);
  CHECK_EQ(got.entity1, expected.entity1);
  CHECK_EQ(got.entity2, expected.entity2);
  CHECK_EQ(got.material1, expected.material1);
  CHECK_EQ(got.material2, expected.material2);
}

void printCounts(const char* label, const VoxelCounts& counts) {
  std::cout << "[feature/voxel_entity] " << label << " domain=" << counts.domain
            << " solid=" << counts.solid << " outside=" << counts.outside
            << " entity1=" << counts.entity1 << " entity2=" << counts.entity2
            << " material1=" << counts.material1 << " material2=" << counts.material2
            << std::endl;
}

}  // namespace

int main() {
  constexpr T pi = T(3.14159265358979323846);
  constexpr T R = T(10);
  constexpr T pitch = T(30);
  constexpr T turns = T(1.25);
  constexpr T h = T(1);
  constexpr T bandWidth = T(1.5);
  constexpr int nx = 27;
  constexpr int ny = 27;
  constexpr int nz = 36;

  const AABB<T, 3> domain(Vector<T, 3>{T(-13), T(-13), T(-3)},
                          Vector<T, 3>{T(14), T(14), T(33)});
  CHECK_EQ(domain.getExtension()[0], T(nx) * h);
  CHECK_EQ(domain.getExtension()[1], T(ny) * h);
  CHECK_EQ(domain.getExtension()[2], T(nz) * h);

  StlReader<T> reader("../geo_stl/stl/helical_duct.stl", T(1), T(1), 2);
  CHECK(reader.getMesh().triangleSize() > 0);

  flbgeo::VoxelGeometryBuilder<T> builder(nx, ny, nz, domain, h, bandWidth);

  // The builder has no read accessor, so build() doubles as an occupancy
  // snapshot: it sorts the band and applies any pending patches, then returns
  // a copy.  Repeating it rebuilds the same state.
  builder.addBox(domain, 1, 1, flbgeo::VoxelRole::Domain);
  builder.addStl(reader, 2, 2, flbgeo::VoxelRole::Solid);
  const VoxelCounts afterStl = tally(builder.build());
  printCounts("after box+stl", afterStl);

  CHECK(afterStl.domain > 0);
  CHECK(afterStl.solid > 0);
  CHECK_EQ(afterStl.outside, std::size_t{0});
  CHECK_EQ(afterStl.domain, afterStl.entity1);
  CHECK_EQ(afterStl.domain, afterStl.material1);
  CHECK_EQ(afterStl.solid, afterStl.entity2);
  CHECK_EQ(afterStl.solid, afterStl.material2);

  // A domain box added after the STL must not overwrite solid cells.
  builder.addBox(domain, 1, 1, flbgeo::VoxelRole::Domain);
  const VoxelCounts afterSecondBox = tally(builder.build());
  printCounts("after second box", afterSecondBox);
  checkSameCounts(afterSecondBox, afterStl);

  builder.addPatch(2, flbgeo::PatchKind::Wall, domain);
  const VoxelCounts beforeBuild = tally(builder.build());
  auto geo = builder.build();
  const VoxelCounts afterBuild = tally(geo);
  printCounts("after build", afterBuild);
  checkSameCounts(afterBuild, beforeBuild);
  CHECK_EQ(afterBuild.domain, geo.countRole(flbgeo::VoxelRole::Domain));
  CHECK_EQ(afterBuild.solid, geo.countRole(flbgeo::VoxelRole::Solid));
  CHECK_EQ(afterBuild.outside, geo.countRole(flbgeo::VoxelRole::Outside));
  CHECK_EQ(afterBuild.domain, geo.countEntity(1));
  CHECK_EQ(afterBuild.solid, geo.countEntity(2));

  // Surface band: |distance| <= bandWidth and unit normals.
  const std::vector<std::size_t>& band = geo.getBandCells();
  CHECK(!band.empty());
  for (std::size_t id : band) {
    CHECK(geo.hasSurfaceSample(id));
    const flbgeo::SurfaceSample<T>& sample = geo.getSurfaceSample(id);
    CHECK(std::abs(sample.Distance) <= bandWidth + T(1e-9));
    CHECK(std::abs(sample.Normal.getnorm() - T(1)) < T(1e-6));
  }

  // The wall patch is an overlay on domain band cells of entity 2.
  CHECK_EQ(geo.getPatchCount(), std::size_t{1});
  if (geo.getPatchCount() == 1) {
    const flbgeo::BoundaryPatch& patch = geo.getPatch(0);
    const std::vector<std::size_t>& cells = geo.getPatchCells(0);
    CHECK_EQ(patch.Entity, std::uint16_t{2});
    CHECK(patch.Kind == flbgeo::PatchKind::Wall);
    CHECK(!cells.empty());
    for (std::size_t id : cells) {
      CHECK(geo.isDomain(id));
      CHECK_EQ(geo.getSurfaceSample(id).Entity, std::uint16_t{2});
      CHECK(std::find(band.begin(), band.end(), id) != band.end());
    }
  }

  // Centerline point at t = 0.5 lies inside the STL, hence in the solid.
  const T a = T(2) * pi * turns * T(0.5);
  const Vector<T, 3> point{R * std::cos(a), R * std::sin(a), pitch * T(0.5)};
  const std::size_t centerId = geo.findCellId(point);
  CHECK(centerId != flbgeo::kInvalidCell);
  if (centerId != flbgeo::kInvalidCell) CHECK(geo.isSolid(centerId));

  // The first cell center is the domain's lower corner offset by half a voxel.
  const Vector<T, 3> corner = domain.getMin() + Vector<T, 3>{h / T(2)};
  const std::size_t cornerId = geo.findCellId(corner);
  CHECK(cornerId != flbgeo::kInvalidCell);
  if (cornerId != flbgeo::kInvalidCell) {
    CHECK(geo.isDomain(cornerId));
    CHECK_EQ(geo.getEntity(cornerId), std::uint16_t{1});
    CHECK_EQ(geo.getMaterial(cornerId), std::uint16_t{1});
  }

  // Index <-> location round trips.
  const std::size_t probes[] = {std::size_t{0},
                                std::size_t{1},
                                std::size_t{nx} - 1,
                                std::size_t{nx},
                                std::size_t{nx} * std::size_t{ny} - 1,
                                geo.size() / 2,
                                geo.size() - 1};
  for (std::size_t id : probes) {
    CHECK_EQ(geo.getIndex(geo.getLoc(id)), id);
    CHECK_EQ(geo.findCellId(geo.getPosition(id)), id);
  }

  return flbtest::summary("feature/voxel_entity");
}
