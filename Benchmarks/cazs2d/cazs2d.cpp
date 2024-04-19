
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
using LatSet1 = D2Q5<T>;

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
T Temp_Ini;          // K
T Conc_Ini;          // wt.%
Vector<T, 2> U_Ini;  // mm/s
T U_Max;
T P_char;

/*bcs*/
T Temp_Wall;          // K
T Conc_Wall;          // wt.%
Vector<T, 2> U_Wall;  // mm/s
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
  iniReader param_reader("ZS2Dparam.ini");
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

  Cl = 0;
  Ch = (T_Melt - T_Eute) / m_Liquidus;

  /*output to console*/
  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:"
            << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
            << "----------------------------------------------" << std::endl;
}

int main() {
  std::uint8_t voidflag = std::uint8_t(1);
  std::uint8_t AABBflag = std::uint8_t(2);
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t FI_Flag =
      static_cast<std::uint8_t>(CA::CAType::Fluid | CA::CAType::Interface);

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

  ConcConverter<T> ConcConv(LatSet1::cs2, BaseConv, Conc_Ini);
  ConcConv.ConvertConc_withCExpan(Cl, Ch, Diff_Liq, Solutal_Expan_Coeff);
  // ConcConv.Enable_Non_Uniform_TimeStep(100);

  ZSConverter<T> CAConv(BaseConv, TempConv, ConcConv, T_Melt, T_Eute, m_Solidus,
                        m_Liquidus, GT_Coeff);
  // Conv.ConvertLatentHeat(LatHeat);

  // PhaseDiagramConverter<T> PDConv(TempConv, ConcConv, T_Melt, T_Eute,
  // m_Solidus, m_Liquidus);

  UnitConvManager<T> ConvManager(&BaseConv, &TempConv, &ConcConv, &CAConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  // NS geometry
  Geometry2D<T> Geo0(Ni, Nj, cavity, Cell_Len, Vector<T, 2>{}, AABBflag,
                     voidflag);
  Geo0.SetupBoundary<LatSet0>(AABBflag, BouncebackFlag);
  // thermal geometry
  Geometry2D<T> Geo1(Ni, Nj, cavity, Cell_Len, Vector<T, 2>{}, AABBflag,
                     voidflag);
  Geo1.SetupBoundary<LatSet1>(AABBflag, BouncebackFlag);

  // ------------------ define lattice ------------------
  // velocity field
  VectorFieldAOS<T, LatSet0::d> Velocity(Geo0.getVoxelsNum());
	// lbm
  BasicLattice<T, LatSet0> NSLattice(Geo0, BaseConv, Velocity);

  BasicLattice<T, LatSet1> SOLattice(Geo1, ConcConv, Velocity);

  BasicLattice<T, LatSet1> THLattice(Geo1, TempConv, Velocity);

	  // --------------------- CA ---------------------
  CA::ZhuStefanescu2D<T, LatSet0> CA(
      CAConv, SOLattice, THLattice, NSLattice, Delta, pref_Orine,
      Geo0.findCellId(
          Vector<T, 2>{T(Ni * Cell_Len / 2), T(Nj * Cell_Len / 2)}));

  // --------------------- BCs ---------------------
  // NS
  BBLikeFixedBoundary<
      T, LatSet0, BounceBackLikeMethod<T, LatSet0>::normal_bounceback>
      NS_BB("NS_BB", NSLattice, BouncebackFlag);
  BBLikeMovingBoundary<
      T, LatSet0, BounceBackLikeMethod<T, LatSet0>::normal_bounceback>
      NS_MBB("NS_MBB", NSLattice, CA.getInterface(), CA::CAType::Solid);

  // Conc
  BBLikeFixedBoundary<
      T, LatSet1, BounceBackLikeMethod<T, LatSet1>::normal_bounceback>
      SO_BB("SO_BB", SOLattice, BouncebackFlag);
  BBLikeMovingBoundary<
      T, LatSet1, BounceBackLikeMethod<T, LatSet1>::normal_bounceback>
      SO_MBB("SO_MBB", SOLattice, CA.getInterface(), CA::CAType::Solid);


  Buoyancy<T, LatSet0> Force(NSLattice, Velocity);
  Force.AddSource(&SOLattice);
  Force.AddSource(&THLattice);



  // writer

//   vtkWriter::PhysScalerWriter RhoWriter(
//       "rho", NSLattice.getRhoField().getField(), BaseConv);
  vtkWriter::PhysScalerWriter CWriter(
      "C", SOLattice.getRhoField().getField(), ConcConv);
//   vtkWriter::FlagWriter CellTypwWriter(
//       "CellType", CA.getState().getField());
//   vtkWriter::PhysVelocityWriter_AOS<T, LatSet0::d> VelocityWriter(
//       "velocity", NSLattice.getVelocityField().getField(), BaseConv);
  //
//   vtkWriter::ScalerWriter CurvWriter(
//       "Curv", CA.getCurvature().getField());
//   vtkWriter::ScalerWriter DFsWriter(
//       "DFs", CA.getDeltaFs().getField());
    //   vtkWriter::ScalerWriter FsWriter(
    //   "Fs", CA.getFs().getField());
  vtkStruPointsWriter<T, LatSet0::d> NCWriter("CAZS2D", Geo0);
  NCWriter.addtoWriteList(&CWriter);

  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  NCWriter.Write(MainLoopTimer());
  while (MainLoopTimer() < MaxStep) {
    ++MainLoopTimer;
    ++OutputTimer;

    NSLattice.UpdateRho(CA.getState().getField(), FI_Flag);
    SOLattice.UpdateRho_Source(CA.getState().getField(), FI_Flag, CA.getExcessC_().getField());

    CA.apply_SimpleCapture();

    Force.GetBuoyancy(CA.getState().getField(), FI_Flag);

    Force.BGK_U<Equilibrium<T, LatSet0>::SecondOrder>(
        CA.getState().getField(), CA::CAType::Fluid);
    Force.BGK<Equilibrium<T, LatSet0>::SecondOrder>(
        CA.getState().getField(), CA::CAType::Interface);

    SOLattice.BGK_Source<Equilibrium<T, LatSet1>::SecondOrder>(
        CA.getState().getField(), FI_Flag, CA.getExcessC_().getField());

    NSLattice.Stream();
    SOLattice.Stream();

    NS_BB.Apply();
    NS_MBB.Apply();
    SO_BB.Apply();
    SO_MBB.Apply();

    if (MainLoopTimer() % OutputStep == 0) {
      OutputTimer.Print_InnerLoopPerformance(NSLattice.getN(), OutputStep);
      // Printer::Print_SolidFraction<T>(CellComm.getSolidFraction<T>());
      Printer::Print<int>("Interface", CA.getInterface().size());
      Printer::Print<T>("Solid%", CA.getSolidCountFracton()*100);
      NCWriter.Write(MainLoopTimer());
    }
  }
  NCWriter.Write(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni, Nj);
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();
  return 0;
}

// attention! DO NOT call std::srand(std::time(0)) in a while loop
// If the while loop doesn't take more than a second to complete,
// then time(0) will return the same value every time we call srand(time(0)),
// which means we are setting the same seed for rand() repeatedly.