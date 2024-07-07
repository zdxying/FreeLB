// Poiseuille.cpp
// poiseuille flow 2d


// pressure driven flow in a channel
// pressure gradient is applied in the x direction, i.e. along the channel

// the left and right boundary(inlet and outlet) are set to be
// generalised-periodic boundary: periodic boundary with pressure variations.
// the top and bottom boundary are set to be bounce-back boundary.

// p = cs^2 * rho
// In incompressible flows, the absolute pressure value is defined
// up to an arbitrary constant. Thereby, we set p_out to 1 (simulation units)
// and: p_in = p_out + delta_p,
// where delta_p can be calculated by:
// delta_p = 8 * Conv.Lattice_VisKine * Conv.Lattice_charU * Ni / Nj / Nj
// the rho in the inlet and outlet are set to be constant:
// rho_out = 1, rho_in = 1 + delta_p * LatSet::InvCs2;

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
T Dyna_Visc;  // Pa·s Dynamic viscosity of the liquid
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
  
  iniReader param_reader("Poiseuilleparam.ini");
  // mesh
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  // physical properties
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Dyna_Visc = param_reader.getValue<T>("Physical_Property", "Dyna_Visc");
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
  Lattice.SetPeriod(0);

  // geometry
  Geometry2DLegacy<T> Geo(Ni, Nj);
  // Geo.circle(radius, positionx, positiony, -1);
  Geo.quickset_boundary(1);
  Geo.set_left(2);   // inlet
  Geo.set_right(3);  // outlet

  Velocity2D<T> Field(U_Ini, BaseConv, Geo);

  // bcs
  BounceBack<T, LatSet, LatStru> NS_BB(Lattice);
  NS_BB.SetFromGeo(1, Geo);
  GeneralisedPeriodic<T, LatSet, LatStru> NS_GP(Lattice, Field);
  NS_GP.SetFromGeo(2, Geo);
  NS_GP.SetFromGeo(3, Geo);
  std::vector<BasicBoundary<T, LatSet, LatStru> *> NS_Boundaries{&NS_BB,
                                                                 &NS_GP};
  BoundaryManager2D<T, LatSet, LatStru> NS_BoundaryMan(NS_Boundaries);
  // BoundaryManager2D<T, LatSet, LatStru> NS_BoundaryMan(&NS_BB, nullptr,
  // nullptr,
  //                                                      nullptr, &NS_GP);

  // Field.SetUwallFromVector(NS_GP.GetBoundaryCellIds());

  // lbm method
  lbm2D<T, LatSet, LatStru> NS(Field, BaseConv, NS_BoundaryMan,
                               std::string("NS"), std::string("rho"));
  NS.setPopRho_From_VoxelFlag(2, getRhoInlet<T, LatSet>(BaseConv, Ni, Nj));

  IndexManager2D CellIdxU(Ni, Nj, NS.get_InIdx());
  // updated inlet and outlet velocity
  // top and bottom velocity are set to be zero and will not be updated
  Geo.get_Idx(CellIdxU.Get(), 2);
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
    // NS.apply<Equilibrium<T,LatSet>::Feq_secondOrder>();// use second order eq
    // to avoid divergence

    // NS.apply<Equilibrium<T, LatSet>::Feq_secondOrder,
    //          &lbm2D<T, LatSet, LatStru>::poststream_innerrho,
    //          &lbm2D<T, LatSet, LatStru>::poststream_u>();

    NS.apply<Equilibrium<T, LatSet>::Feq_secondOrder>(NS.get_InIdx(),
                                                      CellIdxU.Get());

    if (MainLoopTimer() % OutputStep == 0) {
      LBM.calcToleranceU();
      OutputTimer.Print_InnerLoopPerformance(Ni, Nj, OutputStep);
      Printer::Print_Res<T>(LBM.Current_Res);
      // datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());
      // getchar();
    }
  }
  datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni, Nj);

  return 0;
}