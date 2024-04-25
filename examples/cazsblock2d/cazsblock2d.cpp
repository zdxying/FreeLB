/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// cazsblock2d.cpp


#include "ca/zhu_stefanescu2d.h"
#include "ca/zhu_stefanescu2d.hh"
#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;
using LatSet0 = D2Q9<T>;
using LatSet1 = D2Q5<T>;
using LatSetCA = CA::D2Q8<T>;

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
  iniReader param_reader("cazsblock2dparam.ini");
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
  Solutal_Expan_Coeff = param_reader.getValue<T>("Phys_Prop", "Solutal_Expan_Coeff");
  Thermal_Expan_Coeff = param_reader.getValue<T>("Phys_Prop", "Thermal_Expan_Coeff");
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
  std::uint8_t VoidFlag = std::uint8_t(1);
  std::uint8_t AABBFlag = std::uint8_t(2);
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t FI_Flag =
    static_cast<std::uint8_t>(CA::CAType::Fluid | CA::CAType::Interface);

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // ------------------ define converters ------------------
  BaseConverter<T> BaseConv(LatSet0::cs2);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);

  TempConverter<T> TempConv(LatSet1::cs2, BaseConv, Temp_Ini);
  // Conv.ConvertTempFromTDiff_with_Ra(Tl, Th, TDiff, Ra);
  // Conv.ConvertTempFromSHeatCap_and_TCond_with_Ra(Tl, Th, T_Cond_Liq,
  // SHeatCap_Liq, Ra);
  TempConv.ConvertTempFromSHeatCap_and_TCond_with_Texpan(Tl, Th, T_Cond_Liq, SHeatCap_Liq,
                                                         Thermal_Expan_Coeff);

  ConcConverter<T> ConcConv(LatSet1::cs2, BaseConv, Conc_Ini);
  ConcConv.ConvertConc_withCExpan(Cl, Ch, Diff_Liq, Solutal_Expan_Coeff);

  ZSConverter<T> CAConv(BaseConv, TempConv, ConcConv, T_Melt, T_Eute, m_Solidus,
                        m_Liquidus, GT_Coeff);

  UnitConvManager<T> ConvManager(&BaseConv, &TempConv, &ConcConv, &CAConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  // NS geometry
  BlockGeometry2D<T> Geo(Ni, Nj, Thread_Num, cavity, Cell_Len);

  // ------------------ define flag field ------------------
  BlockFieldManager<FlagField, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(cavity,
                 [&](FlagField& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet0>(cavity, BouncebackFlag);

  vtmo::ScalerWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag", Geo, 1);
  GeoWriter.addWriterSet(&FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define lattice ------------------
  // velocity field
  BlockFieldManager<VectorFieldAOS<T, 2>, T, 2> VelocityFM(Geo);
  // lbm
  BlockLatticeManager<T, LatSet0> NSLattice(Geo, BaseConv, VelocityFM);

  BlockLatticeManager<T, LatSet1> SOLattice(Geo, ConcConv, VelocityFM);

  BlockLatticeManager<T, LatSet1> THLattice(Geo, TempConv, VelocityFM);

  // --------------------- CA ---------------------
  CA::BlockZhuStefanescu2DManager<T, LatSetCA> CA(VelocityFM, CAConv, SOLattice,
                                                  THLattice, Delta, pref_Orine);
  // set CA State field
  CA.getStateFM().forEach(FlagFM, AABBFlag, [&](auto& field, std::size_t id) {
    field.SetField(id, CA::CAType::Fluid);
  });
  CA.Setup(Geo.getIndex(Vector<int, 2>{Ni / 2, Nj / 2}));
  // --------------------- BCs ---------------------
  // NS
  BBLikeFixedBlockBdManager<T, LatSet0,
                            BounceBackLikeMethod<T, LatSet0>::normal_bounceback>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  BBLikeMovingBlockBdManager<T, LatSet0,
                             BounceBackLikeMethod<T, LatSet0>::normal_bounceback>
    NS_MBB("NS_MBB", NSLattice, CA.getInterfaces(), FlagFM, CA::CAType::Solid);

  // Conc
  BBLikeFixedBlockBdManager<T, LatSet1,
                            BounceBackLikeMethod<T, LatSet1>::normal_bounceback>
    SO_BB("SO_BB", SOLattice, FlagFM, BouncebackFlag, VoidFlag);

  BBLikeMovingBlockBdManager<T, LatSet1,
                             BounceBackLikeMethod<T, LatSet1>::normal_bounceback>
    SO_MBB("SO_MBB", SOLattice, CA.getInterfaces(), FlagFM, CA::CAType::Solid);


  BlockBuoyancyManager<T, LatSet0> Force(NSLattice, VelocityFM);
  Force.AddSource(SOLattice);
  Force.AddSource(THLattice);

  // writer
  vtmo::ScalerWriter CWriter("Conc", SOLattice.getRhoFM());
  // vtmo::ScalerWriter FsWriter("Fs", CA.getFsFM());
  // vtmo::ScalerWriter EXCWriter("ExC", CA.getExcessCFM());
  vtmo::ScalerWriter StateWriter("State", CA.getStateFM());
  vtmo::VectorWriter VecWriter("Velocity", VelocityFM);
  vtmo::vtmWriter<T, 2> MainWriter("cazsblock2d", Geo, 1);
  MainWriter.addWriterSet(&CWriter, &StateWriter, &VecWriter);

  vtmo::ScalerWriter RhoWriter("Rho", NSLattice.getRhoFM());
  MainWriter.addWriterSet(&RhoWriter);

  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  MainWriter.WriteBinary(MainLoopTimer());

  while (MainLoopTimer() < MaxStep) {
    NSLattice.UpdateRho(MainLoopTimer(), FI_Flag, CA.getStateFM());
    SOLattice.UpdateRho_Source(MainLoopTimer(), FI_Flag, CA.getStateFM(),
                               CA.getExcessCFM());

    CA.Apply_SimpleCapture();

    Force.GetBuoyancy(MainLoopTimer(), FI_Flag, CA.getStateFM());

    Force.BGK_U<Equilibrium<T, LatSet0>::SecondOrder>(MainLoopTimer(), CA::CAType::Fluid,
                                                      CA.getStateFM());
    Force.BGK<Equilibrium<T, LatSet0>::SecondOrder>(
      MainLoopTimer(), CA::CAType::Interface, CA.getStateFM());

    SOLattice.BGK_Source<Equilibrium<T, LatSet1>::SecondOrder>(
      MainLoopTimer(), FI_Flag, CA.getStateFM(), CA.getExcessCFM());

    NSLattice.Stream(MainLoopTimer());
    SOLattice.Stream(MainLoopTimer());

    NS_BB.Apply(MainLoopTimer());
    NS_MBB.Apply(MainLoopTimer());
    SO_BB.Apply(MainLoopTimer());
    SO_MBB.Apply(MainLoopTimer());

    NSLattice.Communicate(MainLoopTimer());
    SOLattice.Communicate(MainLoopTimer());

    // SOLattice.getRhoFM().CommunicateAll(MainLoopTimer());

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      VelocityFM.CommunicateAll();
      SOLattice.getRhoFM().CommunicateAll();
      NSLattice.getRhoFM().CommunicateAll();

      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Print<std::size_t>("Interface", CA.getInterfaceNum());
      Printer::Print<T>("Solid%", T(CA.getSolidCount()) * 100 / Geo.getTotalCellNum());
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());
    }
  }

  MainWriter.WriteBinary(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();
  return 0;
}