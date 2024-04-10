
#include "ca/zhu_stefanescu2d.h"
#include "ca/zhu_stefanescu2d.hh"
#include "freelb.h"
#include "freelb.hh"

// int Total_Macro_Step = 0;
using T = FLOAT;
using LatSet = D2Q5<T>;
using LatStru = lat::D2Q5;
/*----------------------------------------------
              Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
T Cell_Len;
T RT;

T TimeStep;
int Thread_Num;
/*---------------------
        physical param
---------------------*/
/*physical property*/
T rho_ref;    // g/mm^3
T Diff_Liq;   // mm^2/s Diffusion coefficient in liquid
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
T Ra;         // Rayleigh number
/*init conditions*/
T Conc_Ini;  // wt.%
T U_Ini[2];  // mm/s
T U_Max;
T P_char;

/*bcs*/
T Conc_Wall;  // wt.%
T U_Wall[2];
/*---------------------
        LB param
---------------------*/
T Ch;  // char high conc
T Cl;  // char low conc

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  /*reader*/
  iniReader param_reader("SDT.ini");
  /*mesh*/
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  /*CA mesh*/
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  /*physical property*/
  rho_ref = param_reader.getValue<T>("Phys_Prop", "rho_ref");
  Diff_Liq = param_reader.getValue<T>("Phys_Prop", "Diff_Liq");
  Kine_Visc = param_reader.getValue<T>("Phys_Prop", "Kine_Visc");
  Ra = param_reader.getValue<T>("Phys_Prop", "Ra");
  // Kine_Visc = Dyna_Visc / rho_ref;
  /*init conditions*/
  Conc_Ini = param_reader.getValue<T>("ICs", "Conc_Ini");
  U_Ini[0] = param_reader.getValue<T>("ICs", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("ICs", "U_Ini1");
  U_Max = param_reader.getValue<T>("ICs", "U_Max");
  P_char = param_reader.getValue<T>("ICs", "P_char");
  Ch = param_reader.getValue<T>("ICs", "Ch");
  Cl = param_reader.getValue<T>("ICs", "Cl");
  /*bcs*/
  Conc_Wall = param_reader.getValue<T>("BCs", "Conc_Wall");
  U_Wall[0] = param_reader.getValue<T>("BCs", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("BCs", "Velo_Wall1");
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
            << "tolerance:         " << tol << "\n"
            << "----------------------------------------------\n"
            << std::endl;
}

int main() {
  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  BaseConverter<T> BaseConv(LatSet::cs2);
  // Conv.SimplifiedConvertFromViscosity(Ni - 2, U_Max, Kine_Visc);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max,
                         Kine_Visc);

  ConcConverter<T> ConcConv(LatSet::cs2, BaseConv, Conc_Ini);
  ConcConv.ConvertConc(Cl, Ch, Diff_Liq);
  // ConcConv.Enable_Non_Uniform_TimeStep(100);

  UnitConvManager<T> ConvManager(&BaseConv, nullptr, &ConcConv);
  ConvManager.Check_and_Print();

  // lattice structure
  LatStru Lattice_D2Q9(Ni, Nj);

  Geometry2DLegacy<T> Geo(Ni, Nj);
  Geo.quickset_boundary(1);

  Velocity2D<T> Field(BaseConv, Geo, U_Ini);

  // --------------------- BCs ---------------------
  // bcs for NS
  BounceBack<T, LatSet, LatStru> NSSO_BB(Lattice_D2Q9);
  NSSO_BB.SetFromGeo(1, Geo);
  BoundaryManager2D<T, LatSet, LatStru> NS_BM(&NSSO_BB);

  // bcs for Conc
  // to use the same BB BCs with NS, the latset and latstru must be the same
  // BoundaryManager2D<T, LatSet, LatStru> SO_BM(&NSSO_BB);

  // --------------------- lbm ---------------------
  lbm2D<T, LatSet, LatStru> SO(Field, ConcConv, NS_BM, std::string("SO"),
                               std::string("C"));
  LBManager2D<T> LBM(Field, nullptr, nullptr, &SO);

  // give a fluctuation to SO field
  SO.setPopRho(Ni / 2, Nj / 2, ConcConv.getLatticeRho(Conc_Ini * 2));

  // writer
  DataWriter2D<T> datWriter(work_dir, Field, ConvManager, &LBM);
  datWriter.Write_Geometry();
  FLBplot<T> plotPop(work_dir, "/FLBplotPop");
  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  datWriter.Write_Rho_Phys_Only(MainLoopTimer());
  plotPop.write_plot(MainLoopTimer(), SO.getAverRho());
  while (MainLoopTimer() < MaxStep) {

    SO.applyRho<Equilibrium<T, LatSet>::Feq_firstOrder,
                &lbm2D<T, LatSet, LatStru>::poststream_rho>();

    // SO.Collide_BGK<Equilibrium<T, LatSet>::Feq_firstOrder>(SO.get_Idx());
    // SO.Stream(SO.get_InIdx());
    // SO.poststream_rho();

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      OutputTimer.Print_InnerLoopPerformance(Ni, Nj, OutputStep);
      datWriter.Write_Rho_Phys_Only(MainLoopTimer());
      plotPop.write_plot(MainLoopTimer(), SO.getAverRho());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni, Nj);

  return 0;
}