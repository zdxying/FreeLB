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
  iniReader param_reader("cavityref2d.ini");
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

int main() {
  std::uint8_t VoidFlag = std::uint8_t(1);
  std::uint8_t AABBFlag = std::uint8_t(2);
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  Printer::Print_BigBanner(std::string("Initializing..."));

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
  // use 0.5 for refined block
  AABB<T, 2> toplid(Vector<T, 2>(Cell_Len / 2, T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T((Ni - 0.5) * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> Coarsecavity(Vector<T, 2>(T(Ni * Cell_Len) / 3 + 1, T(Nj * Cell_Len) / 3 + 1),
                          Vector<T, 2>(T(Ni * Cell_Len) * 2 / 3 - 1, T(Nj * Cell_Len) * 2 / 3 - 1));

  BlockGeometry2D<T> Geo(Ni, Nj, Thread_Num, cavity, Cell_Len, AABBFlag, VoidFlag, 1, true);
  Geo.forEachBlock([&](Block2D<T>& block) {
    if (!isOverlapped(block.getBaseBlock(), Coarsecavity)) {
      // BasicBlock<T, 2> RefBlock = block.getBaseBlock().getRefinedBlock();
      // block = Block2D<T>(RefBlock, AABBFlag, VoidFlag, 2);
      block.Refine();
    }
  });
  Geo.InitAllComm();
  Geo.ReadAABBs(cavity, AABBFlag);
  Geo.SetupBoundary<LatSet>(AABBFlag, BouncebackFlag);
  Geo.setFlag(toplid, BouncebackFlag, BBMovingWallFlag);

  // vtmwriter::ScalerWriter GeoFlagWriter("flag", Geo.getGeoFlags());
  // vtmwriter::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", Geo);
  // GeoWriter.addWriterSet(&GeoFlagWriter);
  // GeoWriter.WriteBinary();

  vtmno::ScalerWriter GeoFlagWriterNO("flagno", Geo.getGeoFlags(), Geo.getGeoMeshes());
  vtmno::vtmWriter<T, LatSet::d> GeoWriterNO("GeoFlagNO", Geo);
  GeoWriterNO.addWriterSet(&GeoFlagWriterNO);
  GeoWriterNO.WriteBinary();

  // ------------------ define lattice ------------------
  // velocity field
  BlockVectFieldAOS<T, 2> Velocity(Geo.getBlockSizes());
  // set initial value of field
  Geo.forEachVoxel(toplid, BBMovingWallFlag, [&Velocity](int id, int blockid) {
    Velocity.getBlockField(blockid).SetField(id, U_Wall);
  });
  // lattice
  BlockLatticeManager<T, LatSet> LatMan(Geo, BaseConv, Velocity);
  LatMan.EnableToleranceU();
  T res = 1;

  // bcs
  BBLikeBlockFixedBdManager<T, LatSet, BounceBackLikeMethod<T, LatSet>::normal_bounceback> NS_BB(
    "NS_BB", LatMan, BouncebackFlag, VoidFlag);
  BBLikeBlockFixedBdManager<T, LatSet, BounceBackLikeMethod<T, LatSet>::movingwall_bounceback> NS_BBMW(
    "NS_BBMW", LatMan, BBMovingWallFlag, VoidFlag);
  BlockBoundaryManager BM(&NS_BB, &NS_BBMW);

  // writers
  vtmno::ScalerWriter RhoWriterNO("Rho", LatMan.getRhoField(), Geo.getGeoMeshes());
  vtmno::VectorWriter VecWriterNO("Velocity", Velocity, Geo.getGeoMeshes());
  vtmno::vtmWriter<T, LatSet::d> NSWriterNO("cavref2dNO", Geo, 1);
  NSWriterNO.addWriterSet(&RhoWriterNO, &VecWriterNO);

  vtmwriter::ScalerWriter RhoWriter("Rho", LatMan.getRhoField());
  vtmwriter::VectorWriter VecWriter("Velocity", Velocity);
  vtmwriter::vtmWriter<T, LatSet::d> NSWriter("cavref2d", Geo);
  NSWriter.addWriterSet(&RhoWriter, &VecWriter);


  // /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;
  NSWriterNO.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && res > tol) {

    LatMan.UpdateRho(MainLoopTimer(), std::uint8_t(AABBFlag));
    LatMan.UpdateU(MainLoopTimer(), std::uint8_t(AABBFlag));
    LatMan.template BGK<Equilibrium<T, LatSet>::SecondOrder>(MainLoopTimer(),
      std::uint8_t(AABBFlag | BouncebackFlag | BBMovingWallFlag));

    // to communicate before stream, non-overlapped vtm writer should be used
    LatMan.Stream(MainLoopTimer());

    LatMan.Communicate(MainLoopTimer());

    BM.Apply(MainLoopTimer());

    ++MainLoopTimer;
    ++OutputTimer;
    if (MainLoopTimer() % OutputStep == 0) {
      res = LatMan.getToleranceU(1);
      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Print_Res<T>(res);
      NSWriterNO.WriteBinary(MainLoopTimer());
    }
  }
  NSWriterNO.WriteBinary(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());

  return 0;
}