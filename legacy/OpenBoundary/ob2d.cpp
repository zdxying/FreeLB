// ob2d.cpp
// open boundary 2d
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

// physical properties
T rho_ref;    // g/mm^3
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
T Ra;         // Rayleigh number
// init conditions
T U_Ini[2];  // m/s
T U_Max;
T P_char;

// bcs
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
  
  iniReader param_reader("ob2dparam.ini");
  // Thread_Num = param_reader.getValue<int>("OMP", "Thread_Num");
  // mesh
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  // physical properties
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");
  // init conditions
  U_Ini[0] = param_reader.getValue<T>("Init_Conditions", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("Init_Conditions", "U_Ini1");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  P_char = param_reader.getValue<T>("Init_Conditions", "P_char");
  // bcs
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
  // Lattice.SetPeriod(0);

  // geometry
  Geometry2DLegacy<T> Geo(Ni, Nj);
  // Geo.circle(radius, positionx, positiony, -1);
  Geo.quickset_boundary(1);
  Geo.set_left(2);   // inlet
  Geo.set_right(3);  // outlet

  Velocity2D<T> Field(U_Ini, BaseConv, Geo);
  Field.Setby_PhysU(2, U_Wall);

  // bcs
  BounceBack<T, LatSet, LatStru> NS_BB(Lattice);
  NS_BB.SetFromGeo(1, Geo);
  // GeneralisedPeriodic<T, LatSet, LatStru> NS_GP(Lattice, Field);
  // NS_GP.SetFromGeo(2, Geo);
  // NS_GP.SetFromGeo(3, Geo);
  // inlet: BBMW
  // BounceBackMovingWall<T, LatSet, LatStru> NS_BBMW(Lattice, Field);
  BounceBackLike<T, LatSet, LatStru,
                 bouncebackmethod<T, LatSet>::movingwall_bounceback>
      NS_BBMW(Lattice, Field, "Bounce-Back-Moving-Wall");
  NS_BBMW.SetFromGeo(2, Geo);
  // outlet: ABB
  // AntiBounceBack<T, LatSet, LatStru> NS_ABB(Lattice);
  BounceBackLike<T, LatSet, LatStru,
                 bouncebackmethod<T, LatSet>::anti_bounceback_pressure>
      NS_ABBP(Lattice, Field, "Anti-Bounce-Back");
  NS_ABBP.SetFromGeo(3, Geo);

  std::vector<BasicBoundary<T, LatSet, LatStru> *> Boundary{&NS_BB, &NS_BBMW,
                                                            &NS_ABBP};
  BoundaryManager2D<T, LatSet, LatStru> NS_BoundaryMan(Boundary);

  // Field.SetUwallFromVector(NS_GP.GetBoundaryCellIds());

  // lbm method
  lbm2D<T, LatSet, LatStru> NS(Field, BaseConv, NS_BoundaryMan, "NS", "rho");

  // set inlet rho may not necessary ----2023.11.7
  // NS.setPopRho_From_VoxelFlag(2, getRhoInlet<T, LatSet>(BaseConv, Ni, Nj));

  // Cell index for Rho update
  IndexManager2D CellIdxRho(Ni, Nj, NS.get_InIdx());
  Geo.get_Idx(CellIdxRho.Get(), 2);
  // Cell index for U update
  IndexManager2D CellIdxU(Ni, Nj, NS.get_InIdx());
  Geo.get_Idx(CellIdxU.Get(), 3);

  // lbm manager
  LBManager2D<T> LBM(Field, &NS);

  // writer
  DataWriter2D<T> datWriter(work_dir, Field, ConvManager, &LBM);
  // datWriter.Write_Geometry();

  // tolerance
  LBM.SetupToleranceU(tol);
  // count and timer
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