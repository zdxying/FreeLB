#include "freelb.h"
#include "freelb.hh"

int Thread_Num = 1;

using T = FLOAT;
using LatSet = D2Q9<T>;

int Ni = 100;
int Nj = 100;
T Cell_Len = 1.;

int main() {
  std::uint8_t cavflag = std::uint8_t(1);
  std::uint8_t domainflag = std::uint8_t(2);
  std::uint8_t circleflag = std::uint8_t(4);
  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  Circle<T> circle(
      T(Ni / 10), Vector<T, 2>(T(Ni * Cell_Len * 0.5), T(Nj * Cell_Len * 0.5)));
  Circle<T> Largercircle(T(Ni / 10) + 1, Vector<T, 2>(T(Ni * Cell_Len * 0.5),
                                                      T(Nj * Cell_Len * 0.5)));
  Geometry2D<T> Geo(Ni, Nj, cavity, Cell_Len);
  Geo.setFlag(circle, circleflag);

  RefinedGeometry2D<T> RefGeo(Geo, 2);
  RefGeo.Setup(Geo.getGeoFlagField().getField(), circleflag, Largercircle,
               2);
  RefGeo.UpdateAbstractTree();
  RefGeo.UpdateRefineLevelsVtk();
  RefGeo.UpdateFlag(Largercircle, circleflag);
  RefGeo.UpdateFlagsVtk();
  RefGeo.UpdateBlockAABB(cavflag, 0.6);
  RefGeo.UpdateBlockAABBVtk();

  // -------------------writers-------------------
  
  vtkWriter::FlagWriter GeoFlagwriter("flag",
                                      Geo.getGeoFlagField().getField());
  vtkWriter::FlagWriter GeoLevelwriter("Level",
                                       RefGeo.getRefineLevels().getField());
  vtkWriter::ScalerWriter GeoBlockwriter("Block",
                                       RefGeo.getBlockIndex().getField());
  vtkStruPointsWriter<T, LatSet::d> GeoWriter("CavGeo", Geo);
  GeoWriter.addtoWriteList(&GeoFlagwriter, &GeoLevelwriter, &GeoBlockwriter);
  GeoWriter.Write();

  vtkWriter::UnStruFlagWriter levelwriter("level", RefGeo.getRefineLevelsVtk());
  vtkWriter::UnStruFlagWriter flagwriter("flag", RefGeo.getFlagsVtk());
  vtkPolyWriter<T, LatSet::d> writer("Geometry");
  writer.addtoWriteList(&levelwriter, &flagwriter);
  writer.Write(RefGeo.getUpdatedAbstractTree());

  vtkUnStruGridWriter<T, LatSet::d> gridwriter("Grid");
  gridwriter.addtoWriteList(&levelwriter, &flagwriter);
  gridwriter.Write(RefGeo.getUpdatedAbstractTree());
}