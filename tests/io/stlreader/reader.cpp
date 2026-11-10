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

// reader.cpp test

#include <algorithm>
#include <cmath>
#include <limits>

#include "freelb.h"
#include "freelb.hh"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D3Q27<T>;

int main() {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t AABBFlag = std::uint8_t(2);
  constexpr std::uint8_t BdFlag = std::uint8_t(4);

  // StlReader<T> reader("./stl/radiator.stl", 0.5);
  StlReader<T> reader("./stl/cylinder100_20_20finer.stl", 2., 2.);

  // Surface queries: signedDistance must match the nearest triangle and take
  // its sign from the octree; evalSurfaceNormal must be unit length and point
  // from the solid into the fluid.
  {
    const std::vector<Triangle<T>>& triangles = reader.getMesh().getTriangles();
    CHECK(!triangles.empty());

    auto insideAt = [&](const Vector<T, 3>& p) {
      Octree<T>* node = reader.getTree()->find(p);
      return node != nullptr && node->isLeaf() && node->getInside();
    };
    auto nearestDistance = [&](const Vector<T, 3>& p) {
      T best = std::numeric_limits<T>::max();
      for (const Triangle<T>& tri : triangles) {
        const Vector<T, 3> q = tri.closestPtPointTriangle(p);
        best = std::min(best, (q - p).getnorm());
      }
      return best;
    };

    const Vector<T, 3> mid = (reader.getMin() + reader.getMax()) * T(0.5);
    const Vector<T, 3> outside =
      reader.getMin() - Vector<T, 3>(reader.getVoxelSize() * T(2));
    const T sdMid = reader.signedDistance(mid);
    const T sdOutside = reader.signedDistance(outside);

    CHECK_EQ(sdMid < T(0), insideAt(mid));
    CHECK(std::abs(std::abs(sdMid) - nearestDistance(mid)) < T(1e-9));
    CHECK(sdOutside > T(0));
    CHECK(!insideAt(outside));
    CHECK(std::abs(sdOutside - nearestDistance(outside)) < T(1e-9));

    Triangle<T>& tri = reader.getMesh().getTri(0);
    const Vector<T, 3> center = tri.getCenter();
    const Vector<T, 3> probe =
      center + tri.getNormal() * (reader.getVoxelSize() * T(0.5));
    const Vector<T, 3> outward =
      insideAt(probe) ? center - probe : probe - center;
    const Vector<T, 3> normal = reader.evalSurfaceNormal(probe);
    CHECK(std::abs(normal.getnorm() - T(1)) < T(1e-6));
    CHECK(normal * outward > T(0.9) * outward.getnorm());
  }

  // reader.getTree()->write("octree");
  // geometry.Setup<LatSet>();

  // Geometry3D
  // Geometry3D<T> geometry(reader, 1);
  // geometry.SetupBoundary<LatSet>(1,2);

  // vtkWriter::FieldFlagWriter<std::uint8_t> flagwriter(
  //   "flag", geometry.getGeoFlagField().getField().getdataPtr(),
  //   geometry.getGeoFlagField().getField().size());
  // vtkStruPointsWriter<T, LatSet::d> writer("Geometry", geometry);
  // writer.addtoWriteList(&flagwriter);
  // writer.Write();


  // BlockGeometry3D
  BlockGeometryHelper3D<T> GeoHelper(reader, 10);
  GeoHelper.CreateBlocks(true);
  GeoHelper.AdaptiveOptimization(4);
  GeoHelper.LoadBalancing();
  BlockGeometry3D<T> blockgeometry(GeoHelper);
  CHECK_EQ(blockgeometry.getBlockNum(), 4);
  CHECK_EQ(blockgeometry.getTotalCellNum(), std::size_t(52272));
  // BlockGeometry3D<T> blockgeometry(reader, 4);
  BlockFieldManager<FlagField, T, LatSet::d> FlagFM(blockgeometry, VoidFlag);
  FlagFM.ReadOctree(reader.getTree(), AABBFlag);
  FlagFM.SetupBoundary<LatSet>(AABBFlag, VoidFlag, BdFlag);

  vtmwriter::ScalarWriter GeoFlagWriter("flag", FlagFM);
  vtmwriter::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", blockgeometry);
  GeoWriter.addWriterSet(GeoFlagWriter);
  GeoWriter.WriteBinary();

  return flbtest::summary("io/stlreader");

  
  // radius height centre
  // Cylinder<T> cylinder(T(15),
  //                      Vector<T, 3>{T(0), T(0), T(30)},
  //                      Vector<T, 3>{T(20), T(20), T(5)});
  // Cylinder<T> cylinder(T(5),
  //                      Vector<T, 3>{T(0), T(0), T(20)},
  //                      Vector<T, 3>{T(0), T(0), T(10)});
  // geometry.setFlag(cylinder, 2);
  // AABB<T,3> cavity(Vector<T, 3>(T(-2), T(-2), T(-2)),
  //                Vector<T, 3>(T(80), T(40), T(80)));
}