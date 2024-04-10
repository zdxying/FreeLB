#include <mpi.h>

#include "freelb.h"
#include "freelb.hh"

int Thread_Num = 1;

using T = float;
using LatSet = D2Q9<T>;

int Ni = 100;
int Nj = 100;
T Cell_Len = 1.;

T U_Max = 0.1;
T Kine_Visc = 0.01;

int main() {
  std::uint8_t cavflag = std::uint8_t(1);
  std::uint8_t domainflag = std::uint8_t(2);
  std::uint8_t circleflag = std::uint8_t(4);

  // converter
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConvertFromViscosity(Ni, U_Max, Kine_Visc);

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  Circle<T> circle(
      T(Ni / 3), Vector<T, 2>(T(Ni * Cell_Len * 0.5), T(Nj * Cell_Len * 0.5)));

  BlockGeometry2D<T> Geo(Ni, Nj, 4, cavity, Cell_Len);
  Geo.ReadAABBs(cavity, cavflag);
  Geo.setFlag(circle, circleflag);

  BlockVectFieldAOS<T, 2> BlockVelo(Geo.getBlockSizes(), T(0));

  BlockLatticeManager<T, LatSet> LatMan(Geo, BaseConv, BlockVelo);
  // test
  BlockScalerField<T> BlockField(Geo.getBlockSizes(), T(1));
  Geo.forEachVoxel(circle, [&](int id, int blockid) {
    BlockField.getBlockField(blockid).SetField(id, T(2));
  });

  // -------------------writers-------------------

  vtmwriter::ScalerWriter GeoFlagWriter("flag", Geo.getGeoFlags());
  vtmwriter::ScalerWriter Swriter("S", BlockField);

  vtmwriter::vtmWriter<T, LatSet::d> VtmWriter("GeoFlag", Geo);
  VtmWriter.addWriterSet(&GeoFlagWriter);
  VtmWriter.addWriterSet(&Swriter);
  VtmWriter.WriteBinary();

  // vtkWriter::FlagWriter GeoFlagwriter("flag",
  //                                     Geo.getGeoFlagField().getField());

  // vtkStruPointsWriter<T, LatSet::d> GeoWriter("CavGeo", Geo);
  // GeoWriter.addtoWriteList(&GeoFlagwriter);
  // GeoWriter.Write();

  return 0;
}