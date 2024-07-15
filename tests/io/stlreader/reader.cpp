// reader.cpp test

#include "head.h"
#include "freelb.h"
#include "freelb.hh"
#include "io/stlreader.h"
#include "io/stlreader.hh"
#include "geometry/geometry3d.h"

int Thread_Num = 1;

using T = FLOAT;
using LatSet = D3Q7<T>;
int main() {
  // StlReader<T> reader("./stl/radiator.stl", 0.5);
  StlReader<T> reader("./stl/cylinder100_20_20finer.stl", 2., 2.);
  // reader.getTree()->write("octree");
  // geometry.Setup<LatSet>();
  Geometry3D<T> geometry(reader, 1);
  geometry.SetupBoundary<LatSet>(1,2);

  vtkWriter::FieldFlagWriter<std::uint8_t> flagwriter(
      "flag", geometry.getGeoFlagField().getField().getdataPtr(),
      geometry.getGeoFlagField().getField().size());
  vtkStruPointsWriter<T, LatSet::d> writer("Geometry", geometry);
  writer.addtoWriteList(&flagwriter);
  writer.Write();

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