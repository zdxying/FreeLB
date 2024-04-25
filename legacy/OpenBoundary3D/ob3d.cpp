// ob3d.cpp
// open boundary 3d
// this is a benchmark for the freeLB library

// the inlet is set with a constant velocity
// Bounce-Back-Moving-Wall method is used
// the outlet is set with a constant pressure
// Anti-Bounce-Back with pressure method is used

// note that extra IndexManager is used: CellIdxRho and CellIdxU
// CellIdxRho is based on popInIdx but add the inlet boundary cells

// CellIdxU is based on popInIdx but add the outlet boundary cells
// the outlet velocity will be updated but the inlet remains constant
// this way the bounce-back-moving-wall is performed

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
/*---------------------
        physical param
---------------------*/

/*physical property*/
T rho_ref;    // g/mm^3
T Dyna_Visc;  // Pa·s Dynamic viscosity of the liquid
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
T Ra;         // Rayleigh number
/*init conditions*/
Vector<T, 3> U_Ini;  // mm/s
T U_Max;

/*bcs*/
Vector<T, 3> U_Wall;  // mm/s

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  /*reader*/
  iniReader param_reader("ob3dparam.ini");
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
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");
  /*init conditions*/
  U_Ini[0] = param_reader.getValue<T>("Init_Conditions", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("Init_Conditions", "U_Ini1");
  U_Ini[2] = param_reader.getValue<T>("Init_Conditions", "U_Ini2");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
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
            << "----------------------------------------------\n"
            << std::endl;
}

int main() {
  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Nj, U_Max, RT);
  // Conv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // define geometry
  Cylinder<T> cylinder(T(Ni / 2), Vector<T, 3>(T(0), T(Nj), T(0)),
                       Vector<T, 3>(T(Ni / 2), T(0), T(Nk / 2)));
  Cylinder<T> Inletc(T((Ni - 2) * Cell_Len / 2),
                     Vector<T, 3>(T(0), T(Cell_Len), T(0)),
                     Vector<T, 3>(T(Ni * Cell_Len / 2), T(-0.5 * Cell_Len),
                                  T(Nk * Cell_Len / 2)));
  Cylinder<T> Outletc(
      T((Ni - 2) * Cell_Len / 2), Vector<T, 3>(T(0), T(Cell_Len), T(0)),
      Vector<T, 3>(T(Ni * Cell_Len / 2), T(Nj * Cell_Len - 0.5 * Cell_Len),
                   T(Nk * Cell_Len / 2)));

  // AABB<T,3> Inlet(
  //     Vector<T, 3>(T(0), T(-0.5 * Cell_Len), T(0)),
  //     Vector<T, 3>(T(Ni * Cell_Len), T(0.5 * Cell_Len), T(Nk * Cell_Len)));
  // AABB<T,3> Outlet(
  //     Vector<T, 3>(T(0), T(Nj * Cell_Len - 0.5 * Cell_Len), T(0)),
  //     Vector<T, 3>(T(Ni * Cell_Len), T(Nj * Cell_Len + 0.5 * Cell_Len),
  //                  T(Nk * Cell_Len)));
  VoxelGeometry3D<T> Geo(Ni, Nj, Nk, cylinder, Cell_Len);
  Geo.Setup<LatSet>();
  Geo.setFlag(Inletc, 1, 2);
  Geo.setFlag(Outletc, 1, 3);
  Geo.WriteStruPoints();

  VelocityField3D<T> Field(BaseConv, U_Ini, Geo);
  Field.setVelocity(U_Wall, 2);

  // bcs
  BounceBackLike3D<T, LatSet, BBlikemethod<T, LatSet>::normal_bounceback>
      NS_BB(1, "NS_BB");
  BounceBackLike3D<T, LatSet, BBlikemethod<T, LatSet>::movingwall_bounceback>
      NS_BBMW(2, "NS_BBMW");
  BounceBackLike3D<T, LatSet,
                   BBlikemethod<T, LatSet>::anti_bounceback_pressure>
      NS_ABBP(3, "Anti-Bounce-Back");
  BoundaryManager3D<T, LatSet> BM(&NS_BB, &NS_BBMW, &NS_ABBP);

  // lbm method
  lbm3D<T, LatSet> NS(Field, BaseConv, BM, "NS", "rho");
  // NS.addtoPostRhoIdx(2);
  NS.addtoPostUIdx(3);
  NS.EnableToleranceU();
  T res = 1;

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && res > tol) {
    ++MainLoopTimer;
    ++OutputTimer;
    NS.Run<Equilibrium<T, LatSet>::Feq_secondOrder, true, true>();
    // NS.Run<&Equilibrium3D<T, LatSet>::Feq_secondOrder, true, true>();
    if (MainLoopTimer() % OutputStep == 0) {
      res = NS.getToleranceU();
      OutputTimer.Print_InnerLoopPerformance(Ni * Nj * Nk, OutputStep);
      Printer::Print_Res<T>(res);
    }
  }
  NS.WriteStruPoints(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni * Nj * Nk);
  return 0;
}