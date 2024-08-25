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

#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;
using LatSet = D3Q7<T>;

int main() {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t AABBFlag = std::uint8_t(2);
  constexpr std::uint8_t BdFlag = std::uint8_t(4);

  // StlReader<T> reader("./stl/radiator.stl", 0.5);
  StlReader<T> reader("./stl/cylinder100_20_20finer.stl", 2., 2.);
  // reader.getTree()->write("octree");
  // geometry.Setup<LatSet>();

  // Geometry3D
  Geometry3D<T> geometry(reader, 1);
  geometry.SetupBoundary<LatSet>(1,2);

  vtkWriter::FieldFlagWriter<std::uint8_t> flagwriter(
    "flag", geometry.getGeoFlagField().getField().getdataPtr(),
    geometry.getGeoFlagField().getField().size());
  vtkStruPointsWriter<T, LatSet::d> writer("Geometry", geometry);
  writer.addtoWriteList(&flagwriter);
  writer.Write();


  // BlockGeometry3D
  BlockGeometry3D<T> blockgeometry(reader, 2);
  BlockFieldManager<FlagField, T, LatSet::d> FlagFM(blockgeometry, VoidFlag);
  FlagFM.ReadOctree(reader.getTree(), AABBFlag);
  FlagFM.SetupBoundary<LatSet>(AABBFlag, VoidFlag, BdFlag);

  vtmwriter::ScalarWriter GeoFlagWriter("flag", FlagFM);
  vtmwriter::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", blockgeometry);
  GeoWriter.addWriterSet(GeoFlagWriter);
  GeoWriter.WriteBinary();

  
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


  return 0;
}