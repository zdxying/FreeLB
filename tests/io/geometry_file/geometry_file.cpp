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

// geometry_file.cpp
// B1b: the FLBGEO file must round-trip a full VoxelGeometry and its recipe.
// The fixture is the shared helical duct; occupancy, entities, materials,
// patches, the surface band, and the source list are compared cell by cell.

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>

#include "freelb.h"
#include "freelb.hh"
#include "io/geometry_file.hh"

#include "../../test_util.h"

using T = FLOAT;

int Thread_Num = 1;

namespace {

template <typename T>
void checkVector(const Vector<T, 3>& got, const Vector<T, 3>& expected) {
  for (int i = 0; i < 3; ++i) CHECK_EQ(got[i], expected[i]);
}

}  // namespace

int main() {
  constexpr T h = T(1);
  constexpr T bandWidth = T(1.5);
  constexpr int nx = 27;
  constexpr int ny = 27;
  constexpr int nz = 36;
  const char* stlPath = "../../feature/geo_stl/stl/helical_duct.stl";
  const char* filePath = "./geometry_roundtrip.flbgeo";

  const AABB<T, 3> domain(Vector<T, 3>{T(-13), T(-13), T(-3)},
                          Vector<T, 3>{T(14), T(14), T(33)});

  StlReader<T> reader(stlPath, T(1), T(1), 2);
  CHECK(reader.getMesh().triangleSize() > 0);

  flbgeo::VoxelGeometryBuilder<T> builder(nx, ny, nz, domain, h, bandWidth);
  builder.addBox(domain, 1, 1, flbgeo::VoxelRole::Domain);
  builder.addStl(reader, 2, 2, flbgeo::VoxelRole::Solid);
  builder.addPatch(2, flbgeo::PatchKind::Wall, domain);
  const flbgeo::VoxelGeometry<T> geometry = builder.build();

  CHECK_EQ(geometry.size(), std::size_t(nx) * std::size_t(ny) * std::size_t(nz));
  CHECK(geometry.countRole(flbgeo::VoxelRole::Domain) > 0);
  CHECK(geometry.countRole(flbgeo::VoxelRole::Solid) > 0);
  CHECK(!geometry.getBandCells().empty());
  CHECK_EQ(geometry.getPatchCount(), std::size_t{1});
  CHECK(!geometry.getPatchCells(0).empty());

  flbgeo::GeometryRecipe<T> recipe;
  recipe.BandWidth = bandWidth;

  flbgeo::GeometrySource<T> boxSource;
  boxSource.Kind = flbgeo::GeometrySourceKind::Box;
  boxSource.Box = domain;
  boxSource.Entity = 1;
  boxSource.Material = 1;
  boxSource.Role = flbgeo::VoxelRole::Domain;
  boxSource.RequireClosed = false;

  flbgeo::GeometrySource<T> stlSource;
  stlSource.Kind = flbgeo::GeometrySourceKind::Stl;
  stlSource.Path = "stl/helical_duct.stl";
  stlSource.Entity = 2;
  stlSource.Material = 2;
  stlSource.Role = flbgeo::VoxelRole::Solid;
  stlSource.RequireClosed = true;

  recipe.Sources = {boxSource, stlSource};

  flbgeo::writeGeometryFile(filePath, geometry, recipe);

  {
    std::ifstream file(filePath, std::ios::binary | std::ios::ate);
    CHECK(file.good());
    CHECK(file.tellg() > std::streamoff{0});
  }

  const flbgeo::GeometryFile<T> loaded =
    flbgeo::readGeometryFile<T>(filePath);
  const flbgeo::VoxelGeometry<T>& got = loaded.Geometry;

  CHECK_EQ(got.getNx(), geometry.getNx());
  CHECK_EQ(got.getNy(), geometry.getNy());
  CHECK_EQ(got.getNz(), geometry.getNz());
  CHECK_EQ(got.size(), geometry.size());
  checkVector(got.getDomain().getMin(), geometry.getDomain().getMin());
  checkVector(got.getDomain().getMax(), geometry.getDomain().getMax());
  CHECK_EQ(got.getVoxelSize(), geometry.getVoxelSize());

  for (std::size_t id = 0; id < geometry.size(); ++id) {
    CHECK_EQ(static_cast<int>(got.getRole(id)),
             static_cast<int>(geometry.getRole(id)));
    CHECK_EQ(got.getEntity(id), geometry.getEntity(id));
    CHECK_EQ(got.getMaterial(id), geometry.getMaterial(id));
  }

  CHECK_EQ(got.getPatchCount(), geometry.getPatchCount());
  for (std::size_t p = 0; p < geometry.getPatchCount(); ++p) {
    CHECK_EQ(got.getPatch(p).Entity, geometry.getPatch(p).Entity);
    CHECK_EQ(static_cast<int>(got.getPatch(p).Kind),
             static_cast<int>(geometry.getPatch(p).Kind));
    CHECK_EQ(got.getPatchCells(p).size(), geometry.getPatchCells(p).size());
    for (std::size_t i = 0; i < geometry.getPatchCells(p).size(); ++i) {
      CHECK_EQ(got.getPatchCells(p)[i], geometry.getPatchCells(p)[i]);
    }
  }

  const std::vector<std::size_t>& srcBand = geometry.getBandCells();
  const std::vector<std::size_t>& dstBand = got.getBandCells();
  CHECK_EQ(got.getBandSize(), geometry.getBandSize());
  CHECK_EQ(dstBand.size(), srcBand.size());
  for (std::size_t i = 0; i < srcBand.size(); ++i) {
    CHECK_EQ(dstBand[i], srcBand[i]);
    const flbgeo::SurfaceSample<T>& expected = geometry.getSurfaceSample(srcBand[i]);
    const flbgeo::SurfaceSample<T>& actual = got.getSurfaceSample(dstBand[i]);
    CHECK_EQ(actual.Distance, expected.Distance);
    checkVector(actual.Normal, expected.Normal);
    CHECK_EQ(actual.Entity, expected.Entity);
  }

  CHECK_EQ(loaded.Recipe.BandWidth, recipe.BandWidth);
  CHECK_EQ(loaded.Recipe.Sources.size(), recipe.Sources.size());
  for (std::size_t i = 0; i < recipe.Sources.size(); ++i) {
    const flbgeo::GeometrySource<T>& expected = recipe.Sources[i];
    const flbgeo::GeometrySource<T>& actual = loaded.Recipe.Sources[i];
    CHECK_EQ(static_cast<int>(actual.Kind), static_cast<int>(expected.Kind));
    CHECK_EQ(actual.Entity, expected.Entity);
    CHECK_EQ(actual.Material, expected.Material);
    CHECK_EQ(static_cast<int>(actual.Role), static_cast<int>(expected.Role));
    CHECK_EQ(actual.RequireClosed, expected.RequireClosed);

    switch (expected.Kind) {
      case flbgeo::GeometrySourceKind::Stl:
        CHECK_EQ(actual.Path, expected.Path);
        break;
      case flbgeo::GeometrySourceKind::Box:
        checkVector(actual.Box.getMin(), expected.Box.getMin());
        checkVector(actual.Box.getMax(), expected.Box.getMax());
        break;
      case flbgeo::GeometrySourceKind::Sphere:
        checkVector(actual.Center, expected.Center);
        CHECK_EQ(actual.Radius, expected.Radius);
        break;
    }
  }

  std::remove(filePath);
  return flbtest::summary("io/geometry_file");
}
