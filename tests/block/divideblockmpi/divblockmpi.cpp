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

void readParam() {
  /*reader*/
  iniReader param_reader("divblockm.ini");
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
}

int main(int argc, char* argv[]) {
  std::uint8_t VoidFlag = std::uint8_t(1);
  std::uint8_t AABBFlag = std::uint8_t(2);
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  // Printer::Print_BigBanner(std::string("Initializing..."));

  Mpi().init(&argc, &argv);

  MPI_DEBUG_WAIT

  int world_size = Mpi().getSize();

  readParam();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)), Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> toplid(Vector<T, 2>(Cell_Len, T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T((Ni - 1) * Cell_Len), T(Nj * Cell_Len)));
  // test geocomm
  Circle<T> circle(T(Ni * Cell_Len / 8), Vector<T, 2>(T(Ni * Cell_Len / 2), T(Nj * Cell_Len / 2)));

  BlockGeometryMPIHelper2D<T> GeoHelper(Ni, Nj, world_size, cavity, Cell_Len);
  const BasicBlock<T, 2>& Block = GeoHelper.getBlock(Mpi().getRank());
  const AABB<int, 2>& idxblock = Block.getIdxBlock();
  const AABB<T, 2>& aabb = Block.getAABB();

  Block2D<T> Geo(aabb, idxblock, Mpi().getRank(), Cell_Len, AABBFlag, VoidFlag);
  GeoHelper.InitMPIBlockCommStru(Geo.getMPIBlockComm());
  mpi::MPIBlockBufferInit(Geo.getMPIBlockComm(), Geo.getMPIBlockBuffer());

  
}