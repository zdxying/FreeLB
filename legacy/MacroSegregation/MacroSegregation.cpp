
#include "freelb.h"
#include "freelb.hh"

// int Total_Macro_Step = 0;
using T = FLOAT;
using LatSet0 = D2Q9<T>;
using LatSet1 = D2Q5<T>;

using LatStru0 = lat::D2Q9;
using LatStru1 = lat::D2Q5;
using LatStruca = calat::D2Q8;
/*----------------------------------------------
                        Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
T Cell_Len;
T RT;

T TimeStep;

/*---------------------
        physical param
---------------------*/

/*nucleation and growth*/
T Nuc_Dens_Surf;  // m^-2  = 1 mm^-2   -------10^3 m^-1
T Nuc_Dens_Bulk;  // m^-3  = 0.1 mm^-3  -----1.2*10^4 m^-2
T DT_Mean_Surf;   // K Average undercooling
T DT_Std_Surf;    // K Standard deviation of undercooling
T DT_Mean_Bulk;   // K Average undercooling
T DT_Std_Bulk;    // K Standard deviation of undercooling
T Growth_Para;
T GT_Coeff;  // m*K Gibbs-Thomson coefficient

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
T U_Ini[2];  // m/s
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
  iniReader param_reader("MSparam.ini");
  /*mesh*/
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  /*CA mesh*/
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  /*nucleation and growth*/
  Nuc_Dens_Surf =
      param_reader.getValue<T>("Nuc_and_Growth", "Nuc_Dens_Surf");
  Nuc_Dens_Bulk =
      param_reader.getValue<T>("Nuc_and_Growth", "Nuc_Dens_Bulk");
  DT_Mean_Surf =
      param_reader.getValue<T>("Nuc_and_Growth", "DT_Mean_Surf");
  DT_Std_Surf = param_reader.getValue<T>("Nuc_and_Growth", "DT_Std_Surf");
  DT_Mean_Bulk =
      param_reader.getValue<T>("Nuc_and_Growth", "DT_Mean_Bulk");
  DT_Std_Bulk = param_reader.getValue<T>("Nuc_and_Growth", "DT_Std_Bulk");
  Growth_Para = param_reader.getValue<T>("Nuc_and_Growth", "Growth_Para");
  GT_Coeff = param_reader.getValue<T>("Nuc_and_Growth", "GT_Coeff");
  /*Phase diagram*/
  T_Melt = param_reader.getValue<T>("Phase_Diagram", "T_Melt");
  T_Eute = param_reader.getValue<T>("Phase_Diagram", "T_Eute");
  m_Liquidus = param_reader.getValue<T>("Phase_Diagram", "m_Liquidus");
  m_Solidus = param_reader.getValue<T>("Phase_Diagram", "m_Solidus");
  /*physical property*/
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Solutal_Expan_Coeff =
      param_reader.getValue<T>("Physical_Property", "Solutal_Expan_Coeff");
  Thermal_Expan_Coeff =
      param_reader.getValue<T>("Physical_Property", "Thermal_Expan_Coeff");
  SHeatCap_Liq = param_reader.getValue<T>("Physical_Property", "SHeatCap_Liq");
  SHeatCap_Soli =
      param_reader.getValue<T>("Physical_Property", "SHeatCap_Soli");
  LatHeat = param_reader.getValue<T>("Physical_Property", "LatHeat");
  T_Cond_Liq = param_reader.getValue<T>("Physical_Property", "T_Cond_Liq");
  T_Cond_Soli = param_reader.getValue<T>("Physical_Property", "T_Cond_Soli");
  T_Cond_Amb = param_reader.getValue<T>("Physical_Property", "T_Cond_Amb");
  Diff_Soli = param_reader.getValue<T>("Physical_Property", "Diff_Soli");
  Diff_Liq = param_reader.getValue<T>("Physical_Property", "Diff_Liq");
  Dyna_Visc = param_reader.getValue<T>("Physical_Property", "Dyna_Visc");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");
  Ra = param_reader.getValue<T>("Physical_Property", "Ra");
  // Kine_Visc = Dyna_Visc / rho_ref;
  TDiff = param_reader.getValue<T>("Physical_Property", "TDiff");
  /*init conditions*/
  Temp_Ini = param_reader.getValue<T>("Init_Conditions", "Temp_Ini");
  Th = param_reader.getValue<T>("Init_Conditions", "Th");
  Tl = param_reader.getValue<T>("Init_Conditions", "Tl");
  Conc_Ini = param_reader.getValue<T>("Init_Conditions", "Conc_Ini");
  U_Ini[0] = param_reader.getValue<T>("Init_Conditions", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("Init_Conditions", "U_Ini1");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  P_char = param_reader.getValue<T>("Init_Conditions", "P_char");
  /*bcs*/
  Conc_Wall = param_reader.getValue<T>("Boundary_Conditions", "Conc_Wall");
  Temp_Wall = param_reader.getValue<T>("Boundary_Conditions", "Temp_Wall");
  U_Wall[0] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall1");
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

  TempConverter<T> TempConv(LatSet1::cs2, BaseConv);
  // Conv.ConvertTempFromTDiff_with_Ra(Tl, Th, TDiff, Ra);
  // Conv.ConvertTempFromSHeatCap_and_TCond_with_Ra(Tl, Th, T_Cond_Liq,
  // SHeatCap_Liq, Ra);
  TempConv.ConvertTempFromSHeatCap_and_TCond_with_Texpan(
      Tl, Th, T_Cond_Liq, SHeatCap_Liq, Thermal_Expan_Coeff);

  ConcConverter<T> ConcConv(LatSet1::cs2, BaseConv);
  ConcConv.ConvertConc_withCExpan(Cl, Ch, Diff_Liq, Solutal_Expan_Coeff);
  ConcConv.Enable_Non_Uniform_TimeStep(100);

  GandinConverter<T> CAConv(BaseConv, TempConv, ConcConv, T_Melt, T_Eute,
  m_Solidus, m_Liquidus);
  CAConv.ConvertCA(DT_Mean_Bulk, DT_Std_Bulk, DT_Mean_Surf, DT_Std_Surf,
                   Nuc_Dens_Bulk, Nuc_Dens_Surf, Growth_Para);
  // Conv.ConvertLatentHeat(LatHeat);

  // PhaseDiagramConverter<T> PDConv(TempConv, ConcConv, T_Melt, T_Eute, m_Solidus,
                                  // m_Liquidus);

  UnitConvManager<T> ConvManager(&BaseConv, &TempConv, &ConcConv, &CAConv);
  ConvManager.Check_and_Print();

  // lattice structure
  LatStru0 Lattice_D2Q9(Ni, Nj);
  LatStru1 Lattice_D2Q5(Ni, Nj);

  Geometry2DLegacy<T> Geo(Ni, Nj);
  Geo.quickset_boundary(1);

  Velocity2D<T> Field(Geo, U_Ini, U_Wall);

  CAGField2D<T> CAField(Geo);

  // cell communicator, get_method: erase(small mesh size)
  GetGCells2D CellComm(Ni, Nj, CAField.State, 1);

  // --------------------- BCs ---------------------
  // bcs for NS
  BounceBack<T, LatSet0, LatStru0> NS_BB(Lattice_D2Q9);
  NS_BB.SetFromGeo(1, Geo);
  BoundaryManager2D<T, LatSet0, LatStru0> NS_BM(&NS_BB);

  // bcs for Conc
  BounceBack<T, LatSet1, LatStru1> SO_BB(Lattice_D2Q5);
  SO_BB.SetFromGeo(1, Geo);
  BoundaryManager2D<T, LatSet1, LatStru1> SO_BM(&SO_BB);

  // bcs for thermal
  BounceBack<T, LatSet1, LatStru1> TH_BB(Lattice_D2Q5);
  Geo.set_left(2);
  Geo.set_right(3);
  TH_BB.SetFromGeo(1, Geo);
  AntiBounceBack<T, LatSet1, LatStru1> TH_ABB(Lattice_D2Q5);
  TH_ABB.SetFromGeo(2, Geo);
  TH_ABB.SetFromGeo(3, Geo);
  BoundaryManager2D<T, LatSet1, LatStru1> TH_BM(&TH_BB, &TH_ABB);

  CellComm.enable_boundflag(NS_BB.GetBdCell_Ids());

  // --------------------- lbm ---------------------
  // lbm method
  lbm2D<T, LatSet0, LatStru0> NS(T(1), Field, BaseConv, NS_BM,
                                 CellComm.Get_Cells(), std::string("NS"),
                                 std::string("rho"));
  NS.EnableForce();

  lbm2D<T, LatSet1, LatStru1> TH(T(1), Field, TempConv, TH_BM,
                                 std::string("TH"), std::string("T"));
  TH.setPopRho_From_VoxelFlag(2, T(1));
  TH.setPopRho_From_VoxelFlag(3, T(0));
  TH.AddtoInIdx(1);

  lbm2D<T, LatSet1, LatStru1> SO(T(1), Field, ConcConv, SO_BM,
                                 CellComm.Get_Cells(), std::string("SO"),
                                 std::string("C"));

  LBManager2D<T> LBM(Field, &NS, &TH, &SO);

  // --------------------- CA ---------------------
  GandinCA2D<T, LatStruca> CAM(CAField, CAConv, LBM, CellComm);

  // writer
  DataWriter2D<T> datWriter(work_dir, Field, ConvManager, &LBM);

  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  // datWriter.Write_Fluid_Phys_And_Rho_Phys(MainLoopTimer());

  CAM.SingleNuc(Ni/2, Nj/2, T(0));

  while (MainLoopTimer() < MaxStep) {
    // get cells
    CellComm.Get();
    // lbm
    NS.applyF<Equilibrium<T, LatSet0>::Feq_secondOrder,
              ForceMethod2D<T, LatSet0>::VelocityO2_SpatialO2>();
    TH.applyRho<Equilibrium<T, LatSet1>::Feq_secondOrder,
                &lbm2D<T, LatSet1, LatStru1>::poststream_innerrho>();
    // Force term
    LBM.set_ThermalandSolutalBuoyancy();

    // source term

    // must enable non-uniform time step
    if (MainLoopTimer() % ConcConv.TimeStepCoeff == 0) {
      SO.applyRho<Equilibrium<T, LatSet1>::Feq_secondOrder,
                  &lbm2D<T, LatSet1, LatStru1>::poststream_innerrho>();
      SO.UpdateBCs(CAM.getNewGrowings(), CellComm.Get_isbound());
    }

    //ca
    CAM.apply_SingleNuc();

    // lbm bcs
    NS.UpdateBCs(CAM.getNewGrowings(), CellComm.Get_isbound());


    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      OutputTimer.Print_InnerLoopPerformance(Ni, Nj, OutputStep);
      Printer::Print_Res<T>(LBM.Current_Res);
      datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());
    }
  }

  datWriter.Write_Fluid_Rho_Phys(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni, Nj);

  return 0;
}

// attention! DO NOT call std::srand(std::time(0)) in a while loop
// If the while loop doesn't take more than a second to complete,
// then time(0) will return the same value every time we call srand(time(0)),
// which means we are setting the same seed for rand() repeatedly.