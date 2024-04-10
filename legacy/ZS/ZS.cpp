
#include "ca/zhu_stefanescu2d.h"
#include "ca/zhu_stefanescu2d.hh"
#include "freelb.h"
#include "freelb.hh"
// Known bugs: Segmentation fault may occur, but running rhe executable again
// may resolve this without re-compile
//  this may be caused by parallel

// int Total_Macro_Step = 0;
using T = FLOAT;
using LatSet0 = D2Q9<T>;
using LatSet1 = D2Q9<T>;

using LatStru0 = lat::D2Q9;
using LatStru1 = lat::D2Q9;
using LatStruca = calat::D2Q8;
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

/*nucleation and growth*/
T GT_Coeff;  // mm*K Gibbs-Thomson coefficient
T Delta;     // anisotropy coefficient
T pref_Orine;

/*Phase diagram*/
T T_Melt;      // K
T T_Eute;      // K
T m_Liquidus;  // abs slope of The liquidus;
T m_Solidus;   // abs slope of The solidus;

/*physical property*/
T rho_ref;              // g/mm^3
T Solutal_Expan_Coeff;  // wt.%^-1 Solutal expansion coefficient
T Thermal_Expan_Coeff;  // K^-1 Thermal expansion coefficient
T SHeatCap_Liq;         // J·g^−1·K^−1 specific heat capacity of liquid
T SHeatCap_Soli;        // J·g^−1·K^−1 specific heat capacity of solid
T LatHeat;              // J·g^−1 Enthalpy of fusion
T T_Cond_Liq;           // W·mm^−1·K^−1 Thermal conductivity of liquid
T T_Cond_Soli;          // W·mm^−1·K^−1 Thermal conductivity of solid
T T_Cond_Amb;           // W·mm^−1·K^−1 Thermal conductivity of ambient
T Diff_Soli;            // mm^2/s Diffusion coefficient of solute
T Diff_Liq;             // mm^2/s Diffusion coefficient in liquid
T Dyna_Visc;            // Pa·s Dynamic viscosity of the liquid
T Kine_Visc;            // mm^2/s kinematic viscosity of the liquid
T TDiff;                // mm^2/s Thermal diffusivity of the liquid
T Ra;                   // Rayleigh number
/*init conditions*/
T Temp_Ini;  // K
T Conc_Ini;  // wt.%
T U_Ini[2];  // mm/s
T U_Max;
T P_char;

/*bcs*/
T Temp_Wall;  // K
T Conc_Wall;  // wt.%
T U_Wall[2];
/*---------------------
        LB param
---------------------*/
T Th;  // char high temp
T Tl;  // char low temp
T Ch;  // char high conc
T Cl;  // char low conc

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  /*reader*/
  iniReader param_reader("ZSparam.ini");
  /*mesh*/
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  /*CA mesh*/
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  /*nucleation and growth*/
  GT_Coeff = param_reader.getValue<T>("Nuc_and_Growth", "GT_Coeff");
  Delta = param_reader.getValue<T>("Nuc_and_Growth", "Delta");
  pref_Orine = param_reader.getValue<T>("Nuc_and_Growth", "pref_Orine");
  /*Phase diagram*/
  T_Melt = param_reader.getValue<T>("Phase_Diagram", "T_Melt");
  T_Eute = param_reader.getValue<T>("Phase_Diagram", "T_Eute");
  m_Liquidus = param_reader.getValue<T>("Phase_Diagram", "m_Liquidus");
  m_Solidus = param_reader.getValue<T>("Phase_Diagram", "m_Solidus");
  /*physical property*/
  rho_ref = param_reader.getValue<T>("Phys_Prop", "rho_ref");
  Solutal_Expan_Coeff =
      param_reader.getValue<T>("Phys_Prop", "Solutal_Expan_Coeff");
  Thermal_Expan_Coeff =
      param_reader.getValue<T>("Phys_Prop", "Thermal_Expan_Coeff");
  SHeatCap_Liq = param_reader.getValue<T>("Phys_Prop", "SHeatCap_Liq");
  SHeatCap_Soli = param_reader.getValue<T>("Phys_Prop", "SHeatCap_Soli");
  LatHeat = param_reader.getValue<T>("Phys_Prop", "LatHeat");
  T_Cond_Liq = param_reader.getValue<T>("Phys_Prop", "T_Cond_Liq");
  T_Cond_Soli = param_reader.getValue<T>("Phys_Prop", "T_Cond_Soli");
  T_Cond_Amb = param_reader.getValue<T>("Phys_Prop", "T_Cond_Amb");
  Diff_Soli = param_reader.getValue<T>("Phys_Prop", "Diff_Soli");
  Diff_Liq = param_reader.getValue<T>("Phys_Prop", "Diff_Liq");
  Dyna_Visc = param_reader.getValue<T>("Phys_Prop", "Dyna_Visc");
  Kine_Visc = param_reader.getValue<T>("Phys_Prop", "Kine_Visc");
  Ra = param_reader.getValue<T>("Phys_Prop", "Ra");
  // Kine_Visc = Dyna_Visc / rho_ref;
  TDiff = param_reader.getValue<T>("Phys_Prop", "TDiff");
  /*init conditions*/
  Temp_Ini = param_reader.getValue<T>("ICs", "Temp_Ini");
  Th = param_reader.getValue<T>("ICs", "Th");
  Tl = param_reader.getValue<T>("ICs", "Tl");
  Conc_Ini = param_reader.getValue<T>("ICs", "Conc_Ini");
  U_Ini[0] = param_reader.getValue<T>("ICs", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("ICs", "U_Ini1");
  U_Max = param_reader.getValue<T>("ICs", "U_Max");
  P_char = param_reader.getValue<T>("ICs", "P_char");
  /*bcs*/
  Conc_Wall = param_reader.getValue<T>("BCs", "Conc_Wall");
  Temp_Wall = param_reader.getValue<T>("BCs", "Temp_Wall");
  U_Wall[0] = param_reader.getValue<T>("BCs", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("BCs", "Velo_Wall1");
  // LB
  RT = param_reader.getValue<T>("LB", "RT");
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");
  tol = param_reader.getValue<T>("tolerance", "tol");

  Cl = 0;
  Ch = (T_Melt - T_Eute) / m_Liquidus;

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

  BaseConverter<T> BaseConv(LatSet0::cs2);
  // Conv.SimplifiedConvertFromViscosity(Ni - 2, U_Max, Kine_Visc);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max,
                         Kine_Visc);

  TempConverter<T> TempConv(LatSet1::cs2, BaseConv, Temp_Ini);
  // Conv.ConvertTempFromTDiff_with_Ra(Tl, Th, TDiff, Ra);
  // Conv.ConvertTempFromSHeatCap_and_TCond_with_Ra(Tl, Th, T_Cond_Liq,
  // SHeatCap_Liq, Ra);
  TempConv.ConvertTempFromSHeatCap_and_TCond_with_Texpan(
      Tl, Th, T_Cond_Liq, SHeatCap_Liq, Thermal_Expan_Coeff);

  ConcConverter<T> ConcConv(LatSet0::cs2, BaseConv, Conc_Ini);
  ConcConv.ConvertConc_withCExpan(Cl, Ch, Diff_Liq, Solutal_Expan_Coeff);
  // ConcConv.Enable_Non_Uniform_TimeStep(100);

  ZSConverter<T> CAConv(BaseConv, TempConv, ConcConv, T_Melt, T_Eute, m_Solidus,
                        m_Liquidus, GT_Coeff);
  // Conv.ConvertLatentHeat(LatHeat);

  // PhaseDiagramConverter<T> PDConv(TempConv, ConcConv, T_Melt, T_Eute,
  // m_Solidus, m_Liquidus);

  UnitConvManager<T> ConvManager(&BaseConv, &TempConv, &ConcConv, &CAConv);
  ConvManager.Check_and_Print();

  // lattice structure
  LatStru0 Lattice_D2Q9(Ni, Nj);
  LatStru1 Lattice_D2Q5(Ni, Nj);

  Geometry2DLegacy<T> Geo(Ni, Nj);
  Geo.quickset_boundary(1);
  //// Grow with forced convection
  // BounceBack<T, LatSet0, LatStru0> SO_BB(Lattice_D2Q9);
  AntiBounceBack<T, LatSet0, LatStru0> SO_BB(Lattice_D2Q9);
  SO_BB.SetFromGeo(1, Geo);
  std::vector<BasicBoundary<T, LatSet0, LatStru0> *> SOBoundary {&SO_BB};
  BoundaryManager2D<T, LatSet0, LatStru0> SO_BM(SOBoundary);

  Geo.set_left(2);   // inlet
  Geo.set_right(3);  // outlet
  // index for poststream update
  std::vector<int> inlet;
  std::vector<int> outlet;
  Geo.get_Idx(inlet, 2);
  Geo.get_Idx(outlet, 3);

  ////

  Velocity2D<T> Field(U_Ini, BaseConv, Geo);
  //// Grow with forced convection
  Field.Setby_PhysU(2, U_Wall);
  ////
  CAZSField2D<T> CAField(Geo);

  // --------------------- BCs ---------------------
  // bcs for NS
  BounceBack<T, LatSet0, LatStru0> NS_BB(Lattice_D2Q9);
  NS_BB.SetFromGeo(1, Geo);

  //// Grow with forced convection
  BounceBackLike<T, LatSet0, LatStru0,
                 bouncebackmethod<T, LatSet0>::movingwall_bounceback>
      NS_BBMW(Lattice_D2Q9, Field, "Bounce-Back-Moving-Wall");
  NS_BBMW.SetFromGeo(2, Geo);
  // outlet: ABB
  // AntiBounceBack<T, LatSet0, LatStru0> NS_ABB(Lattice);
  BounceBackLike<T, LatSet0, LatStru0,
                 bouncebackmethod<T, LatSet0>::anti_bounceback_pressure>
      NS_ABBP(Lattice_D2Q9, Field, "Anti-Bounce-Back");
  NS_ABBP.SetFromGeo(3, Geo);

  std::vector<BasicBoundary<T, LatSet0, LatStru0> *> NSBoundary{
      &NS_BB, &NS_BBMW, &NS_ABBP};

  BoundaryManager2D<T, LatSet0, LatStru0> NS_BM(NSBoundary);
  ////

  // bcs for Conc
  // to use the same BB BCs with NS, the latset and latstru must be the same
  // BoundaryManager2D<T, LatSet0, LatStru0> SO_BM(&NSSO_BB);

  // bcs for thermal
  BounceBack<T, LatSet1, LatStru1> TH_BB(Lattice_D2Q5);
  Geo.set_left(2);
  Geo.set_right(3);
  TH_BB.SetFromGeo(1, Geo);
  AntiBounceBack<T, LatSet1, LatStru1> TH_ABB(Lattice_D2Q5);
  TH_ABB.SetFromGeo(2, Geo);
  TH_ABB.SetFromGeo(3, Geo);
  BoundaryManager2D<T, LatSet1, LatStru1> TH_BM(&TH_BB, &TH_ABB);

  // cell communicator, get_method: erase(small mesh size)
  CellIndexManager2D<LatSet0, LatStru0> CellComm(Ni, Nj, CAField.State,
                                                 SO_BB.GetBdCell_Dirs());

  // --------------------- lbm ---------------------
  // lbm method
  lbm2D<T, LatSet0, LatStru0> NS(Field, BaseConv, NS_BM, std::string("NS"),
                                 std::string("rho"), CellComm.Get_Cells(),
                                 CellComm.Get_InCells());
  NS.EnableForce();

  lbm2D<T, LatSet0, LatStru0> SO(Field, ConcConv, SO_BM, std::string("SO"),
                                 std::string("C"), CellComm.Get_Cells(),
                                 CellComm.Get_InCells());

  lbm2D<T, LatSet1, LatStru1> TH(Field, TempConv, TH_BM, std::string("TH"),
                                 std::string("T"));
  // TH.setPopRho_From_VoxelFlag(2, T(1));
  // TH.setPopRho_From_VoxelFlag(3, T(0));
  TH.AddtoInIdx(1);

  LBManager2D<T> LBM(Field, &NS, &TH, &SO);

  // --------------------- CA ---------------------
  ZhuStefanescu2D<T, LatStruca> CA(CAField, CAConv, LBM, Delta, pref_Orine,
                                   Index2D::GetId(Ni / 2, Nj / 2, Ni), 8);

  // writer
  DataWriter2D<T> datWriter(work_dir, Field, ConvManager, &LBM, &CAField);
  datWriter.Write_Geometry();
  FLBplot<T> plotZS(work_dir, "/FLBplotZS");
  // FLBplot<T> plotPop(work_dir, "/FLBplotPop");
  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  // datWriter.Write_Fluid_Rho_Phys_And_ZSCA(MainLoopTimer());
  // datWriter.Write_MovingBCs(MainLoopTimer(), NSSO_BB.GetBd_Ids());
  // datWriter.Wrtie_FromFlag(MainLoopTimer(), CA.getInterface(),
  //                          std::string("/Interface_"));
  // datWriter.Write_ZSCA_Only(MainLoopTimer());
  // datWriter.Write_ZSCA(MainLoopTimer(), CAField);
  // datWriter.Write_Rho_Phys_Only(MainLoopTimer());
  datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());
  plotZS.write_plot(MainLoopTimer(), SO.getAverRho());
  // plotPop.write_plot(MainLoopTimer(), SO.getAverRho());

  CellComm.Get(CA.Trans_GetInterface());

  while (MainLoopTimer() < MaxStep) {
    // lbm
    NS.applyF<Equilibrium<T, LatSet0>::Feq_secondOrder,
              ForceMethod2D<T, LatSet0>::VelocityO2_SpatialO2,
              &lbm2D<T, LatSet0, LatStru0>::poststream_innerrho,
              &lbm2D<T, LatSet0, LatStru0>::poststream_inneru>();

    //// Grow with forced convection
    // NS.apply<Equilibrium<T, LatSet0>::Feq_secondOrder,
    //          &lbm2D<T, LatSet0, LatStru0>::poststream_innerrho,
    //          &lbm2D<T, LatSet0, LatStru0>::poststream_inneru>();
    // //// additional post stream
    // NS.compute_rho(inlet);
    // NS.compute_u(outlet);
    //// end of Grow with forced convection

    SO.applyRho<Equilibrium<T, LatSet0>::Feq_firstOrder,
                &lbm2D<T, LatSet0, LatStru0>::poststream_innerrho>();
    // SO.applyRho_Partial<Equilibrium<T, LatSet0>::Feq_firstOrder,
    //                     &lbm2D<T, LatSet0,
    //                     LatStru0>::poststream_rho>(CAField.f);
    // TH.applyRho<Equilibrium<T, LatSet1>::Feq_secondOrder>();
    // datWriter.Write_Rho_Phys_Only(MainLoopTimer()+1);
    // Force term
    LBM.set_SolutalBuoyancy();

    // source term

    // must enable non-uniform time step
    // if (MainLoopTimer() % ConcConv.TimeStepCoeff == 0) {
    //   SO.applyRho<Equilibrium<T, LatSet0>::Feq_secondOrder,
    //               &lbm2D<T, LatSet1, LatStru0>::poststream_innerrho>();
    // }

    // ca
    CA.apply_SLICapture();
    // CA.apply_SimpleCapture();
    // CA.apply_MixedCapture();
    // lbm bcs
    // get cells
    CellComm.Get(CA.Trans_GetInterface());

    ++MainLoopTimer;
    ++OutputTimer;
    // plotZS.write_plot(MainLoopTimer(), CA.getStatisticalPopRho());
    if (MainLoopTimer() % OutputStep == 0) {
      // datWriter.Wrtie_FromFlag(MainLoopTimer(), CA.getInterface(),
      //                          std::string("/Interface_"));
      // datWriter.Write_MovingBCs(MainLoopTimer(), NSSO_BB.GetBd_Ids());
      // datWriter.Write_ZSCA_Only(MainLoopTimer());
      datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());
      // datWriter.Write_Rho_Phys_Only(MainLoopTimer());
      // datWriter.Write_ZSCA(MainLoopTimer(), CAField);
      plotZS.write_plot(MainLoopTimer(), CA.getStatisticalPopRho());
      // plotPop.write_plot(MainLoopTimer(), SO.getAverRho());

      OutputTimer.Print_InnerLoopPerformance(Ni, Nj, OutputStep);
      Printer::Print_SolidFraction<T>(CellComm.getSolidFraction<T>());
      Printer::Print<int>("Interface Cells", CA.Trans_GetInterface().size());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni, Nj);
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  return 0;
}

// attention! DO NOT call std::srand(std::time(0)) in a while loop
// If the while loop doesn't take more than a second to complete,
// then time(0) will return the same value every time we call srand(time(0)),
// which means we are setting the same seed for rand() repeatedly.