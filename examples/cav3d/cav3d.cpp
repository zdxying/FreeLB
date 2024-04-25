// Cavity.cpp
// Lid-driven cavity flow 3d
// this is a benchmark for the freeLB library

// the top wall is set with a constant velocity,
// while the other walls are set with a no-slip boundary condition
// Bounce-Back-like method is used:
// Bounce-Back-Moving-Wall method for the top wall
// Bounce-Back method for the other walls

#include "freelb.h"
#include "freelb.hh"

// int Total_Macro_Step = 0;
using T = FLOAT;
using LatSet = D3Q19<T>;

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
int Nk;
T Cell_Len;
T RT;
int Thread_Num;

/*physical property*/
T rho_ref;    // g/mm^3
T Dyna_Visc;  // Pa·s Dynamic viscosity of the liquid
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
T Ra;         // Rayleigh number
/*init conditions*/
Vector<T, 3> U_Ini;  // mm/s
T U_Max;
T P_char;

/*bcs*/
Vector<T, 3> U_Wall;  // mm/s

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  /*reader*/
  iniReader param_reader("Cavityparam3D.ini");
  // Thread_Num = param_reader.getValue<int>("OMP", "Thread_Num");
  /*mesh*/
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  /*CA mesh*/
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Nk = param_reader.getValue<int>("Mesh", "Nk");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  /*physical property*/
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Dyna_Visc = param_reader.getValue<T>("Physical_Property", "Dyna_Visc");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");
  /*init conditions*/
  U_Ini[0] = param_reader.getValue<T>("Init_Conditions", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("Init_Conditions", "U_Ini1");
  U_Ini[2] = param_reader.getValue<T>("Init_Conditions", "U_Ini2");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  P_char = param_reader.getValue<T>("Init_Conditions", "P_char");
  /*bcs*/
  U_Wall[0] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall1");
  U_Wall[2] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall2");
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
            << "----------------------------------------------\n"
            << std::endl;
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

  // define geometry
  AABB<T, 3> cavity(
      Vector<T, 3>(T(0), T(0), T(0)),
      Vector<T, 3>(T(Ni * Cell_Len), T(Nj * Cell_Len), T(Nk * Cell_Len)));
  AABB<T, 3> toplid(Vector<T, 3>(Cell_Len, Cell_Len, T((Nk - 1) * Cell_Len)),
                    Vector<T, 3>(T((Ni - 1) * Cell_Len), T((Nj - 1) * Cell_Len),
                                 T(Nk * Cell_Len)));
  Geometry3D<T> Geo(Ni, Nj, Nk, cavity, Cell_Len);
  // boundary flag: 2 for no-slip, 3 for moving wall
  Geo.SetupBoundary<LatSet>();
  Geo.setFlag(toplid, BouncebackFlag, BBMovingWallFlag);

  vtkWriter::FieldFlagWriter<std::uint8_t> flagwriter(
      "flag", Geo.getGeoFlagField().getField().getdata(),
      Geo.getGeoFlagField().getField().size());
  vtkStruPointsWriter<T, LatSet::d> GeoWriter("CavGeo", Geo);
  GeoWriter.addtoWriteList(&flagwriter);
  GeoWriter.Write();

  // ------------------ define lattice ------------------
  // velocity field
  VectorFieldAOS<T, LatSet::d> Velocity(Geo.getVoxelsNum());
  // set initial value of field
  Geo.forEachVoxel(toplid, BBMovingWallFlag,
                   [&Velocity](int id) { Velocity.SetField(id, U_Wall); });
  // lattice
  BasicLattice<T, LatSet> NSLattice(Geo, BaseConv, Velocity);
  NSLattice.EnableToleranceU();
  // bcs
  BBLikeFixedBoundary<
      T, LatSet, BounceBackLikeMethod<T, LatSet>::normal_bounceback>
      NS_BB("NS_BB", NSLattice, BouncebackFlag);
  BBLikeFixedBoundary<
      T, LatSet, BounceBackLikeMethod<T, LatSet>::movingwall_bounceback>
      NS_BBMW("NS_BBMW", NSLattice, BBMovingWallFlag);
  BoundaryManager BM(&NS_BB, &NS_BBMW);


  T res = 1;
  // BasicLBM<T, LatSet> NS(NSLattice, BM);

  vtkWriter::FieldScalerWriter<T> RhoWriter(
      "rho", NSLattice.getRhoField().getField().getdata(),
      NSLattice.getRhoField().getField().size());
  vtkWriter::FieldVectorWriter_AOS<T, LatSet::d> VelocityWriter(
      "velocity", NSLattice.getVelocityField().getField().getdata(),
      NSLattice.getVelocityField().getField().size());
  vtkStruPointsWriter<T, LatSet::d> NSWriter("NS", Geo);
  NSWriter.addtoWriteList(&RhoWriter, &VelocityWriter);

  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  while (MainLoopTimer() < MaxStep && res > tol)
  {
    ++MainLoopTimer;
    ++OutputTimer;

    // NSLattice.UpdateRho(NSLattice.getIndex());
    // NSLattice.UpdateU(NSLattice.getInnerIndex());
    // NS.Apply<Equilibrium<T, LatSet>::Feq_secondOrder>();
    NSLattice.UpdateRho(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag|BouncebackFlag|BBMovingWallFlag));
    NSLattice.UpdateU(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag));
    // NSLattice.BGK<Equilibrium<T, LatSet>::SecondOrder>();
    NSLattice.BGK<Equilibrium<T, LatSet>::SecondOrder>(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag|BouncebackFlag|BBMovingWallFlag));

    NSLattice.Stream();
    BM.Apply();

    if (MainLoopTimer() % OutputStep == 0) {
      res = NSLattice.getToleranceU();
      OutputTimer.Print_InnerLoopPerformance(NSLattice.getN(), OutputStep);
      Printer::Print_Res<T>(res);
    }
  }
  
  NSWriter.Write(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(NSLattice.getN());

  return 0;
}