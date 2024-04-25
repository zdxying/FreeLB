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
T RT;
int Thread_Num;

/*physical property*/
T rho_ref;    // g/mm^3
T Dyna_Visc;  // Pa·s Dynamic viscosity of the liquid
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
T Ra;         // Rayleigh number
/*init conditions*/
Vector<T, 2> U_Ini;  // mm/s
T U_Max;
T P_char;

/*bcs*/
Vector<T, 2> U_Wall;  // mm/s

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  /*reader*/
  iniReader param_reader("cavityparam2d.ini");
  // Thread_Num = param_reader.getValue<int>("OMP", "Thread_Num");
  /*mesh*/
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  /*CA mesh*/
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  /*physical property*/
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Dyna_Visc = param_reader.getValue<T>("Physical_Property", "Dyna_Visc");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");
  /*init conditions*/
  U_Ini[0] = param_reader.getValue<T>("Init_Conditions", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("Init_Conditions", "U_Ini1");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  P_char = param_reader.getValue<T>("Init_Conditions", "P_char");
  /*bcs*/
  U_Wall[0] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall1");
  // LB
  RT = param_reader.getValue<T>("LB", "RT");
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");
  tol = param_reader.getValue<T>("tolerance", "tol");

  /*output to console*/

  MPI_RANK(0)

  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:"
            << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
            << "Tolerance:         " << tol << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
            << "----------------------------------------------" << std::endl;
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

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConvertFromViscosity(Ni, U_Max, Kine_Visc);
  // BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max,
  // Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)), Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> toplid(Vector<T, 2>(Cell_Len, T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T((Ni - 1) * Cell_Len), T(Nj * Cell_Len)));

  BlockGeometryMPIHelper2D<T> GeoHelper(Ni, Nj, world_size, cavity, Cell_Len);
  // TODO: create BlockGeometry instead of Block from GeoHelper
  const BasicBlock<T, 2>& Block = GeoHelper.getBlock(Mpi().getRank());
  const AABB<int, 2>& idxblock = Block.getIdxBlock();
  const AABB<T, 2>& aabb = Block.getAABB();

  Block2D<T> Geo(aabb, idxblock, Mpi().getRank(), Cell_Len, AABBFlag, VoidFlag);
  Geo.ReadAABBs(cavity, AABBFlag);
  Geo.SetupBoundary<LatSet>(cavity, AABBFlag, BouncebackFlag);
  GeoHelper.InitMPIBlockCommStru(Geo.getMPIBlockComm());
  mpi::MPIBlockBufferInit(Geo.getMPIBlockComm(), Geo.getMPIBlockBuffer());
  Geo.GeoFlagMPIComm();
  Geo.setFlag(toplid, BouncebackFlag, BBMovingWallFlag);
  Geo.GeoFlagMPIComm();

  vtmwriter::ScalerWriter GeoFlagWriter("flag", Geo.getGeoFlagField());
  vtmwriter::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", Geo.getSelfBlock());
  GeoWriter.addWriterSet(&GeoFlagWriter);
  GeoWriter.MPIWriteBinary();

  // ------------------ define lattice ------------------
  // velocity field
  VectorFieldAOS<T, LatSet::d> Velocity(Geo.getN());
  // set initial value of field
  Geo.forEach(toplid, Geo.getGeoFlagField().getField(), BBMovingWallFlag,
              [&Velocity](int id) { Velocity.SetField(id, U_Wall); });
  // lattice
  BlockLattice<T, LatSet> BlockLat(Geo, BaseConv, Velocity);
  mpi::MPIBlockBufferInit(Geo.getMPIBlockComm(), BlockLat.getMPIBlockBuffer(), LatSet::q);
  BlockLat.EnableToleranceU();
  // bcs
  BBLikeFixedBlockBdManager<T, LatSet, BounceBackLikeMethod<T, LatSet>::normal_bounceback> NS_BB(
    "NS_BB", {&BlockLat}, BouncebackFlag, VoidFlag);
  BBLikeFixedBlockBdManager<T, LatSet, BounceBackLikeMethod<T, LatSet>::movingwall_bounceback> NS_BBMW(
    "NS_BBMW", {&BlockLat}, BBMovingWallFlag, VoidFlag);
  BlockBoundaryManager BM(&NS_BB, &NS_BBMW);

  T res = 1;  // res of each process

  // writers
  vtmwriter::ScalerWriter RhoWriter("Rho", BlockLat.getRhoField());
  vtmwriter::VectorWriter VecWriter("Velocity", Velocity);
  vtmwriter::vtmWriter<T, LatSet::d> NSWriter("cavityblock2d", Geo.getSelfBlock());
  NSWriter.addWriterSet(&RhoWriter, &VecWriter);

  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;
  NSWriter.MPIWriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && res > tol) {
    ++MainLoopTimer;
    ++OutputTimer;

    BlockLat.UpdateRho(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag));
    BlockLat.UpdateU(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag));
    BlockLat.template BGK<Equilibrium<T, LatSet>::SecondOrder>(
      Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag | BouncebackFlag | BBMovingWallFlag));
    BlockLat.Stream();
    BlockLat.MPIcommunicate(Geo.getMPIBlockComm());

    BM.Apply(MainLoopTimer());
    Mpi().barrier();

    if (MainLoopTimer() % OutputStep == 0) {
      res = BlockLat.getToleranceU();
      Mpi().barrier();
      Mpi().reduceAndBcast(res, MPI_MAX);
      OutputTimer.Print_InnerLoopPerformance(GeoHelper.getN(), OutputStep);
      Printer::Print_Res<T>(res);
      NSWriter.MPIWriteBinary(MainLoopTimer());
    }
  }
  NSWriter.MPIWriteBinary(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(GeoHelper.getN());

  // Mpi().~MpiManager();
}