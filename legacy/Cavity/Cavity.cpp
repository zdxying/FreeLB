// Cavity.cpp
// Lid-driven cavity flow 2d
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
using LatSet = D2Q9<T>;
using LatStru = lat::D2Q9;

// namespace lat
// {
// 	using LatSet = D2Q9;
// }
// namespace lat
// using LatInfo =

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
T Cell_Len;
T RT;
int Thread_Num;

/*physical property*/
T rho_ref;   // g/mm^3
T Dyna_Visc; // Pa·s Dynamic viscosity of the liquid
T Kine_Visc; // mm^2/s kinematic viscosity of the liquid
T Ra;        // Rayleigh number
/*init conditions*/
T U_Ini[2]; // m/s
T U_Max;
T P_char;

/*bcs*/
T U_Wall[2];

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam()
{
  /*reader*/
  iniReader param_reader("Cavityparam.ini");
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
  std::cout << "------------Simulation Parameters:-------------\n"
            << std::endl;
  std::cout << "[Simulation_Settings]:"
            << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
            << "Tolerance:         " << tol << "\n"
            #ifdef _OPENMP
            << "Running on " << Thread_Num <<  " threads\n"
            #endif
            << "----------------------------------------------\n"
            << std::endl;
}

int main()
{
  Printer::Print_BigBanner(std::string("Initializing..."));
  
  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConvertFromViscosity(Ni - 2, U_Max, Kine_Visc);
  // ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // lattice structure D2Q9
  LatStru Lattice(Ni, Nj);

  Geometry2DLegacy<T> Geo(Ni, Nj);
  Geo.quickset_boundary(1);
  Geo.set_top(2);

  Velocity2D<T> Field(U_Ini, BaseConv, Geo);
  // bcs
  BounceBack<T, LatSet, LatStru> NS_BB(Lattice);
  NS_BB.SetFromGeo(1, Geo);
  BounceBackMovingWall<T, LatSet, LatStru> NS_BBMW(Lattice, Field);
  NS_BBMW.SetFromGeo(2, Geo);
  BoundaryManager2D<T, LatSet, LatStru> BM(&NS_BB, nullptr, &NS_BBMW);

  // set initial value of field
  Field.Setby_PhysU(NS_BBMW.GetBd_Ids(), U_Wall);

  // lbm method
  lbm2D<T, LatSet, LatStru> NS(Field, BaseConv, BM,
                               std::string("NS"), std::string("rho"));

  // lbm manager
  LBManager2D<T> LBM(Field, &NS);
  LBM.SetupToleranceU(tol);

  // writer
  DataWriter2D<T> datWriter(work_dir, Field, ConvManager, &LBM);
  // datWriter.Write_Geometry();

  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  // datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  /*write dist func*/
  // std::ofstream distfunc;
  // distfunc.open(datWriter.LBCA_dir + "distfunc" + ".dat");
  while (MainLoopTimer() < MaxStep && LBM.NOTConvergedU)
  // while (Total_Step < MaxStep)
  {
    ++MainLoopTimer;
    ++OutputTimer;
    NS.apply<Equilibrium<T, LatSet>::Feq_secondOrder>();
    // datWriter.Write_Step(distfunc, Total_Step);
    // datWriter.Write_distribution_function(distfunc, (Ni/2 + Nj/2*Ni), 0,
    // LB.ADT);
    if (MainLoopTimer() % OutputStep == 0)
    {
      LBM.calcToleranceU();
      OutputTimer.Print_InnerLoopPerformance(Ni, Nj, OutputStep);
      Printer::Print_Res<T>(LBM.Current_Res);
      // datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());
      // getchar();
    }
  }
  // datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());

  // distfunc.close();

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni, Nj);

  return 0;
}