// cavblock2d.cpp

// Lid-driven cavity flow 2d
// this is a benchmark for the freeLB library

// the top wall is set with a constant velocity,
// while the other walls are set with a no-slip boundary condition
// Bounce-Back-like method is used:
// Bounce-Back-Moving-Wall method for the top wall
// Bounce-Back method for the other walls

// block data structure is used

#include "freelb.h"
#include "freelb.hh"

// int Total_Macro_Step = 0;
using T = FLOAT;
using LatSet = D2Q9<T>;

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/

int Thread_Num;
int Ni;
int Nj;
T Cell_Len;
int BlockNum;

void readParam() {
  /*reader*/
  iniReader param_reader("divrefblock2d.ini");
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockNum = param_reader.getValue<int>("Mesh", "BlockNum");
}

int main() {
  std::uint8_t VoidFlag = std::uint8_t(1);
  std::uint8_t AABBFlag = std::uint8_t(2);
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  readParam();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)), Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  // use 0.5 for refined block
  AABB<T, 2> toplid(Vector<T, 2>(Cell_Len / 2, T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T((Ni - 0.5) * Cell_Len), T(Nj * Cell_Len)));
  // use for refined block
  AABB<T, 2> outercavity(Vector<T, 2>(T(Ni * Cell_Len) / 8 + 1, T(Nj * Cell_Len) / 8 + 1),
                         Vector<T, 2>(T(Ni * Cell_Len) * 7 / 8 - 1, T(Nj * Cell_Len) * 7 / 8 - 1));
  AABB<T, 2> innercavity(Vector<T, 2>(T(Ni * Cell_Len) * 3 / 8 + 1, T(Nj * Cell_Len) * 3 / 8 + 1),
                         Vector<T, 2>(T(Ni * Cell_Len) * 5 / 8 - 1, T(Nj * Cell_Len) * 5 / 8 - 1));

  // geometry helper
  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, Ni / 8, cavity, Cell_Len);
  GeoHelper.forEachBlockCell([&](BasicBlock<T, 2>& block) {
    if (!isOverlapped(block, outercavity)) {
      block.refine();
    }
  });
  GeoHelper.forEachBlockCell([&](BasicBlock<T, 2>& block) {
    if (isOverlapped(block, innercavity)) {
      block.refine();
    }
  });
  GeoHelper.CreateBlocks();
  // GeoHelper.Optimize(12);
  GeoHelper.AdaptiveOptimization(16);

  // std::cout << "StdDev: " << GeoHelper.ComputeStdDev() << std::endl;

  BlockGeometry2D<T> Geo(GeoHelper, cavity, AABBFlag, VoidFlag);
  Geo.SetupBoundary<LatSet>(AABBFlag, BouncebackFlag);
  Geo.setFlag(toplid, BouncebackFlag, BBMovingWallFlag);

  vtmno::ScalerWriter GeoFlagWriterNO("flagno", Geo.getGeoFlags(), Geo.getGeoMeshes());
  vtmno::vtmWriter<T, LatSet::d> GeoWriterNO("GeoFlagNO", Geo);
  GeoWriterNO.addWriterSet(&GeoFlagWriterNO);
  GeoWriterNO.WriteBinary();

  // vtmwriter::ScalerWriter GeoFlagWriter("flag", Geo.getGeoFlags());
  // vtmwriter::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", Geo);
  // GeoWriter.addWriterSet(&GeoFlagWriter);
  // GeoWriter.WriteBinary();

  return 0;
}