
#include "ca/zhu_stefanescu3d.h"
#include "ca/zhu_stefanescu3d.hh"
#include "freelb.h"
#include "freelb.hh"
// Known bugs: Segmentation fault may occur, but running rhe executable again
// may resolve this without re-compile
//  this may be caused by parallel

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

T TimeStep;
int Thread_Num;
/*---------------------
        physical param
---------------------*/

/*nucleation and growth*/
T GT_Coeff;  // mm*K Gibbs-Thomson coefficient
T Epsilon;   // anisotropy coefficient
// preferred growth angle around z-axis 1st rotation
T Psi;
// preferred growth angle around x-axis 2nd rotation
T Theta;
// preferred growth angle around z-axis 3rd rotation
T Phi;

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
T Temp_Ini;          // K
T Conc_Ini;          // wt.%
Vector<T, 3> U_Ini;  // mm/s
T U_Max;
T P_char;

/*bcs*/
T Temp_Wall;          // K
T Conc_Wall;          // wt.%
Vector<T, 3> U_Wall;  // mm/s
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
std::string work_dir;

void readParam() {
  /*reader*/
  iniReader param_reader("ZS3Dparam.ini");
  /*mesh*/
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  /*CA mesh*/
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Nk = param_reader.getValue<int>("Mesh", "Nk");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  /*nucleation and growth*/
  GT_Coeff = param_reader.getValue<T>("Nuc_and_Growth", "GT_Coeff");
  Epsilon = param_reader.getValue<T>("Nuc_and_Growth", "Epsilon");
  Psi = param_reader.getValue<T>("Nuc_and_Growth", "Psi");
  Theta = param_reader.getValue<T>("Nuc_and_Growth", "Theta");
  Phi = param_reader.getValue<T>("Nuc_and_Growth", "Phi");
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
  U_Ini[2] = param_reader.getValue<T>("ICs", "U_Ini2");
  U_Max = param_reader.getValue<T>("ICs", "U_Max");
  P_char = param_reader.getValue<T>("ICs", "P_char");
  /*bcs*/
  Conc_Wall = param_reader.getValue<T>("BCs", "Conc_Wall");
  Temp_Wall = param_reader.getValue<T>("BCs", "Temp_Wall");
  U_Wall[0] = param_reader.getValue<T>("BCs", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("BCs", "Velo_Wall1");
  U_Wall[2] = param_reader.getValue<T>("BCs", "Velo_Wall2");
  // LB
  RT = param_reader.getValue<T>("LB", "RT");
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  Cl = 0;
  Ch = (T_Melt - T_Eute) / m_Liquidus;

  /*output to console*/
  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:"
            << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
            << "----------------------------------------------\n"
            << std::endl;
}

int main() {
  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // ---Unit Converters---
  BaseConverter<T> BaseConv(LatSet::cs2);
  // Conv.SimplifiedConvertFromViscosity(Ni - 2, U_Max, Kine_Visc);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max,
                         Kine_Visc);

  TempConverter<T> TempConv(LatSet::cs2, BaseConv, Temp_Ini);
  // Conv.ConvertTempFromTDiff_with_Ra(Tl, Th, TDiff, Ra);
  // Conv.ConvertTempFromSHeatCap_and_TCond_with_Ra(Tl, Th, T_Cond_Liq,
  // SHeatCap_Liq, Ra);
  TempConv.ConvertTempFromSHeatCap_and_TCond_with_Texpan(
      Tl, Th, T_Cond_Liq, SHeatCap_Liq, Thermal_Expan_Coeff);

  ConcConverter<T> ConcConv(LatSet::cs2, BaseConv, Conc_Ini);
  ConcConv.ConvertConc_withCExpan(Cl, Ch, Diff_Liq, Solutal_Expan_Coeff);
  // ConcConv.Enable_Non_Uniform_TimeStep(100);

  ZSConverter<T> CAConv(BaseConv, TempConv, ConcConv, T_Melt, T_Eute, m_Solidus,
                        m_Liquidus, GT_Coeff);
  // Conv.ConvertLatentHeat(LatHeat);

  // PhaseDiagramConverter<T> PDConv(TempConv, ConcConv, T_Melt, T_Eute,
  // m_Solidus, m_Liquidus);

  UnitConvManager<T> ConvManager(&BaseConv, &TempConv, &ConcConv, &CAConv);
  ConvManager.Check_and_Print();

  // ---Geometry---
  AABB<T, 3> cavity(
      Vector<T, 3>(T(0), T(0), T(0)),
      Vector<T, 3>(T(Ni * Cell_Len), T(Nj * Cell_Len), T(Nk * Cell_Len)));
  VoxelGeometry3D<T> Geo(Ni, Nj, Nk, cavity, Cell_Len);
  Geo.Setup<LatSet>();

  // velocity field
  VelocityField3D<T> Field(BaseConv, U_Ini, Geo);
  // CA field
  // CAZSField3D<T> CAField(Field);
  // CAField.SetState(1, 1);

  // --------------------- BCs ---------------------
  // bcs for NS
  BounceBackLike3D<T, LatSet, BBlikemethod<T, LatSet>::normal_bounceback> NS_BB(
      1, "NS_BB");
  BoundaryManager3D<T, LatSet> NS_BM(&NS_BB);
  // moving boundary
  BounceBackLike3D<T, LatSet, BBlikemethod<T, LatSet>::normal_bounceback>
      NS_BBMoving(1, "NS_BBMoving");
  MovingBoundaryManager3D<T, LatSet> MovingNS_BM(NS_BBMoving, Geo,
                                                 CAField.getStates());

  // bcs for Conc
  BounceBackLike3D<T, LatSet, BBlikemethod<T, LatSet>::normal_bounceback> SO_BB(
      1, "SO_BB");
  BoundaryManager3D<T, LatSet> SO_BM(&SO_BB);
  // BoundaryManager2D<T, LatSet> SO_BM(&NSSO_BB);

  // bcs for thermal
  BounceBackLike3D<T, LatSet, BBlikemethod<T, LatSet>::normal_bounceback> TH_BB(
      1, "TH_BB");
  BoundaryManager3D<T, LatSet> TH_BM(&TH_BB);

  // cell communicator, get_method: erase(small mesh size)
  // CellIndexManager2D<LatSet> CellComm(Ni, Nj, CAField.State,
  //                                     SO_BB.GetBdCell_Dirs());

  // --------------------- lbm ---------------------
  // lbm method
  lbm3D<T, LatSet> NS(Field, BaseConv, NS_BM, "NS", "rho",
                      &(CAField.getStates()), &NS_BBMoving);

  lbm3D<T, LatSet> SO(Field, ConcConv, SO_BM, "SO", "C", &(CAField.getStates()),
                      &NS_BBMoving);

  lbm3D<T, LatSet> TH(Field, TempConv, TH_BM, "TH", "T");

  // --------------------- CA ---------------------
  CA::ZhuStefanescu3D<T, LatSet> CA(
      Field, CAConv, TH, SO, Psi, Theta, Phi, Epsilon,
      Index3D::GetId(Ni / 2, Nj / 2, Nk / 2, Ni, Ni * Nj));

  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  MovingNS_BM.Communicate(CA.getInterface());
  // CellComm.Get(CA.Trans_GetInterface());

  // writer
  vtkWriter::FieldScalerWriter<T> WMCWriter("WMC", Geo.getGlobalIdx(),
                                            CA.getWMC());
  vtkWriter::FieldScalerWriter<T> DeltaFsWriter("DeltaFs", Geo.getGlobalIdx(),
                                                CA.getDeltaFs());
  vtkWriter::FieldScalerWriter<T> FsWriter("Fs", Geo.getGlobalIdx(),
                                           CA.getFs());
  vtkWriter::FieldFlagWriter<CA::CAState> StateWriter(
      "State", Geo.getGlobalIdx(), CA.getState());
  vtkStruPointsWriter<T, LatSet::d> vtkwriter("ZS3D", Geo.getVoxelSize(),
                                              Geo.getMin(), Ni, Nj, Nk);
  vtkwriter.addtoWriteList(&WMCWriter, &DeltaFsWriter, &FsWriter, &StateWriter);

  SO.WriteRho(MainLoopTimer());
  vtkwriter.Write(MainLoopTimer());

  while (MainLoopTimer() < MaxStep) {
    // lbm
    // NS.Run_Commun<Equilibrium<T, LatSet>::Feq_secondOrder, true, true>();

    //// Grow with forced convection
    // NS.apply<Equilibrium<T, LatSet>::Feq_secondOrder,
    //          &lbm2D<T, LatSet>::poststream_innerrho,
    //          &lbm2D<T, LatSet>::poststream_inneru>();
    // //// additional post stream
    // NS.compute_rho(inlet);
    // NS.compute_u(outlet);
    //// end of Grow with forced convection

    SO.Run_Commun<Equilibrium<T, LatSet>::Feq_secondOrder, true, true>();
    // SO.applyRho_Partial<Equilibrium<T, LatSet>::Feq_firstOrder,
    //                     &lbm2D<T, LatSet,
    //                     LatStru0>::poststream_rho>(CAField.f);
    // TH.applyRho<Equilibrium<T, LatSet>::Feq_secondOrder>();
    // datWriter.Write_Rho_Phys_Only(MainLoopTimer()+1);
    // Force term
    // LBM.set_SolutalBuoyancy();

    CA.apply_SimpleCapture();

    MovingNS_BM.Communicate(CA.getInterface());

    ++MainLoopTimer;
    ++OutputTimer;
    // plotZS.write_plot(MainLoopTimer(), CA.getStatisticalPopRho());
    if (MainLoopTimer() % OutputStep == 0) {
      OutputTimer.Print_InnerLoopPerformance(Ni * Nj * Nk, OutputStep);
      Printer::Print_SolidFraction<T>(CA.getSolidFraction());
      Printer::Print<int>("Interface Cells", CA.getInterface().size());
      SO.WriteRho(MainLoopTimer());
      vtkwriter.Write(MainLoopTimer());
      // CAField.WriteStates(MainLoopTimer());
      //   CAField.WriteFs(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni * Nj * Nk);
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  return 0;
}

// attention! DO NOT call std::srand(std::time(0)) in a while loop
// If the while loop doesn't take more than a second to complete,
// then time(0) will return the same value every time we call srand(time(0)),
// which means we are setting the same seed for rand() repeatedly.