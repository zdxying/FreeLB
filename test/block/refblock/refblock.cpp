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
int Ni;
int Nj;
T Cell_Len;
int BlockNum;

void readParam() {
  /*reader*/
  iniReader param_reader("refblock.ini");
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

  // Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)), Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> toplid(Vector<T, 2>(Cell_Len, T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T((Ni - 1) * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> Coarsecavity(Vector<T, 2>(T(Ni * Cell_Len) / 3, T(Nj * Cell_Len) / 3),
                          Vector<T, 2>(T(Ni * Cell_Len) * 2 / 3, T(Nj * Cell_Len) * 2 / 3));
  // test geocomm
  Circle<T> circle(T(Ni * Cell_Len / 8), Vector<T, 2>(T(Ni * Cell_Len / 2), T(Nj * Cell_Len / 2)));
  // refblock
  AABB<T, 2> top13cavity(Vector<T, 2>(T(0), T(Nj * Cell_Len / 3 + 1)),
                         Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> top23cavity(Vector<T, 2>(T(0), T(Nj * Cell_Len * 2 / 3 + 1)),
                         Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));

  BlockGeometry2D<T> Geo(Ni, Nj, BlockNum, cavity, Cell_Len, AABBFlag, VoidFlag, 1, true);
  // Geo.RefineBlock(Coarsecavity);
  // use block.refine() instead
  Geo.forEachBlock([&](Block2D<T>& block) {
    if (isOverlapped(block.getBaseBlock(), top13cavity)) {
      block.Refine();
    }
  });
  Geo.forEachBlock([&](Block2D<T>& block) {
    if (isOverlapped(block.getBaseBlock(), top23cavity)) {
      block.Refine();
    }
  });
  Geo.InitAllComm();
  Geo.ReadAABBs(cavity, AABBFlag);
  Geo.SetupBoundary<LatSet>(AABBFlag, BouncebackFlag);
  Geo.setFlag(toplid, BouncebackFlag, BBMovingWallFlag);
  Geo.setFlag(circle, AABBFlag, BBMovingWallFlag);


  vtmwriter::ScalerWriter GeoFlagWriter("flag", Geo.getGeoFlags());
  vtmwriter::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", Geo);
  GeoWriter.addWriterSet(&GeoFlagWriter);
  GeoWriter.WriteBinary();

  vtmno::ScalerWriter GeoFlagWriterNO("flagno", Geo.getGeoFlags(), Geo.getGeoMeshes());
  vtmno::vtmWriter<T, LatSet::d> GeoWriterNO("GeoFlagNO", Geo);
  GeoWriterNO.addWriterSet(&GeoFlagWriterNO);
  GeoWriterNO.WriteBinary();
}