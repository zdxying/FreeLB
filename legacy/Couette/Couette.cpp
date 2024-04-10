// Couette.cpp
// Couette flow 2d
// this is a benchmark for the freeLB library

// the top wall is set with a constant velocity with Bounce-Back-Moving-Wall
// the inlet and outlet are set with a Periodic boundary condition
// the bottom wall is set with a no-slip boundary condition with Bounce-Back

// note that the four corners are set with Bounce-Back-Like method
// do not set with periodic boundary condition

#include "freelb.h"
#include "freelb.hh"

// int Total_Macro_Step = 0;
using T = FLOAT;
using LatSet = D2Q9<T>;
using LatStru = lat::D2Q9;
/*----------------------------------------------
                        Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
T Cell_Len;
T RT;
int Thread_Num;
/*---------------------
        physical param
---------------------*/

/*physical property*/
T rho_ref;    // g/mm^3
T Dyna_Visc;  // Pa·s Dynamic viscosity of the liquid
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
T Ra;         // Rayleigh number
/*init conditions*/
T U_Ini[2];  // m/s
T U_Max;
T P_char;

/*bcs*/
T U_Wall[2];
/*---------------------
        LB param
---------------------*/

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

// geometry
T radius;
T positionx;
T positiony;

void readParam() {
  /*reader*/
  iniReader param_reader("Couetteparam.ini");
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
  // geometry
  radius = param_reader.getValue<T>("Geometry", "radius");
  positionx = param_reader.getValue<T>("Geometry", "positionx");
  positiony = param_reader.getValue<T>("Geometry", "positiony");

  /*output to console*/
  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:"
            << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
            << "tolerance:         " << tol << "\n"
            << "----------------------------------------------\n"
            << std::endl;
}

int main() {
  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni - 2, U_Max, RT);
  // Conv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // lattice structure D2Q9
  LatStru Lattice(Ni, Nj);
  Lattice.SetPeriod(0);

  // geometry
  Geometry2DLegacy<T> Geo(Ni, Nj);
  // Geo.circle(radius, positionx, positiony);
  // Geo.quickset_boundary(1);
  Geo.set_left(3);   // inlety
  Geo.set_right(4);  // outlet
  Geo.set_top(2);    // moving plate
  Geo.set_bottom(1);

  Velocity2D<T> Field(U_Ini, BaseConv, Geo);
  Field.Setby_PhysU(2, U_Wall);

  // bcs
  BounceBack<T, LatSet, LatStru> NS_BB(Lattice);
  NS_BB.SetFromGeo(1, Geo);
  Periodic<T, LatSet, LatStru> NS_P(Lattice);
  NS_P.SetFromGeo(3, Geo);
  NS_P.SetFromGeo(4, Geo);
  // BounceBackMovingWall<T, LatSet, LatStru> NS_BBMW(Lattice, Field);
  BounceBackLike<T, LatSet, LatStru,
                 bouncebackmethod<T, LatSet>::movingwall_bounceback>
      NS_BBMW(Lattice, Field, "Bounce-Back-Moving-Wall");
  NS_BBMW.SetFromGeo(2, Geo);
  // Field.SetUwallFromVector(NS_BBMW.GetBd_Ids());

  std::vector<BasicBoundary<T, LatSet, LatStru> *> Boundary{&NS_BB, &NS_BBMW,
                                                            &NS_P};
  BoundaryManager2D<T, LatSet, LatStru> NS_BoundaryMan(Boundary);

  // Field.SetUwallFromVector(NS_GP.GetBoundaryCellIds());

  // lbm method
  lbm2D<T, LatSet, LatStru> NS(Field, BaseConv, NS_BoundaryMan, "NS", "rho");

  // Cell index for Rho update
  IndexManager2D CellIdxRho(Ni, Nj, NS.get_InIdx());
  Geo.get_Idx(CellIdxRho.Get(), 3);
  Geo.get_Idx(CellIdxRho.Get(), 4);
  // Cell index for U update
  IndexManager2D CellIdxU(Ni, Nj, NS.get_InIdx());
  Geo.get_Idx(CellIdxU.Get(), 1);
  Geo.get_Idx(CellIdxU.Get(), 3);
  Geo.get_Idx(CellIdxU.Get(), 4);

  // lbm manager
  LBManager2D<T> LBM(Field, &NS);

  // writer
  DataWriter2D<T> datWriter(work_dir, Field, ConvManager, &LBM);
  // datWriter.Write_Geometry();

  // tolerance
  LBM.SetupToleranceU(tol);
  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && LBM.NOTConvergedU) {
    ++MainLoopTimer;
    ++OutputTimer;

    // NS.apply<Equilibrium<T, LatSet>::Feq_secondOrder,
    //          &lbm2D<T, LatSet, LatStru>::poststream_innerrho,
    //          &lbm2D<T, LatSet, LatStru>::poststream_inneru>();

    NS.apply<Equilibrium<T, LatSet>::Feq_secondOrder>(CellIdxRho.Get(),
                                                      CellIdxU.Get());
    if (MainLoopTimer() % OutputStep == 0) {
      LBM.calcToleranceU();
      OutputTimer.Print_InnerLoopPerformance(Ni, Nj, OutputStep);
      // datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());
      Printer::Print_Res<T>(LBM.Current_Res);
      // getchar();
    }
  }
  datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni, Nj);

  return 0;
}