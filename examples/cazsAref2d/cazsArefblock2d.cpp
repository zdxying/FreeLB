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

// cazsArefblock2d.cpp
// lattice Boltzmann method coupled with cellular automata(CA) for 2D solidification with adaptive mesh refinement

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
int BlockCellNx;
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

// physical properties
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
// init conditions
T Temp_Ini;          // K
T Conc_Ini;          // wt.%
Vector<T, 2> U_Ini;  // mm/s
T U_Max;

// bcs
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
int RefineCheckStep;

std::string work_dir;

void readParam(std::vector<T>& refThold, std::vector<T>& coaThold) {
  iniReader param_reader("cazsArefblock2dparam.ini");
  // mesh
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellNx = param_reader.getValue<int>("Mesh", "BlockCellNx");
  /*nucleation and growth*/
  GT_Coeff = param_reader.getValue<T>("Nuc_and_Growth", "GT_Coeff");
  Delta = param_reader.getValue<T>("Nuc_and_Growth", "Delta");
  pref_Orine = param_reader.getValue<T>("Nuc_and_Growth", "pref_Orine");
  /*Phase diagram*/
  T_Melt = param_reader.getValue<T>("Phase_Diagram", "T_Melt");
  T_Eute = param_reader.getValue<T>("Phase_Diagram", "T_Eute");
  m_Liquidus = param_reader.getValue<T>("Phase_Diagram", "m_Liquidus");
  m_Solidus = param_reader.getValue<T>("Phase_Diagram", "m_Solidus");
  // physical properties
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
  // init conditions
  Temp_Ini = param_reader.getValue<T>("ICs", "Temp_Ini");
  Th = param_reader.getValue<T>("ICs", "Th");
  Tl = param_reader.getValue<T>("ICs", "Tl");
  Conc_Ini = param_reader.getValue<T>("ICs", "Conc_Ini");
  U_Ini[0] = param_reader.getValue<T>("ICs", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("ICs", "U_Ini1");
  U_Max = param_reader.getValue<T>("ICs", "U_Max");
  // bcs
  Conc_Wall = param_reader.getValue<T>("BCs", "Conc_Wall");
  Temp_Wall = param_reader.getValue<T>("BCs", "Temp_Wall");
  U_Wall[0] = param_reader.getValue<T>("BCs", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("BCs", "Velo_Wall1");
  // LB
  RT = param_reader.getValue<T>("LB", "RT");
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");
  RefineCheckStep = param_reader.getValue<int>("Simulation_Settings", "RefineCheckStep");
  // threshold
  param_reader.getVector<T>("RefTHOLD", refThold);
  param_reader.getVector<T>("CoaTHOLD", coaThold);

  Cl = 0;
  Ch = (T_Melt - T_Eute) / m_Liquidus;


  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:" << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
            << "----------------------------------------------" << std::endl;
}

int main() {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t AABBFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
  constexpr std::uint8_t FI_Flag = CA::CAType::Fluid | CA::CAType::Interface;

  Printer::Print_BigBanner(std::string("Initializing..."));

  std::vector<T> RefThold;
  std::vector<T> CoaThold;
  readParam(RefThold, CoaThold);

  // ------------------ define converters ------------------
  BaseConverter<T> BaseConv(LatSet0::cs2);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);

  TempConverter<T> TempConv(LatSet1::cs2, BaseConv, Temp_Ini);

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
  AABB<T, 2> seedcavity(Vector<T, 2>(T((Ni / 2 - Ni / BlockCellNx) * Cell_Len),
                                     T((Nj / 2 - Ni / BlockCellNx) * Cell_Len)),
                        Vector<T, 2>(T((Ni / 2 + Ni / BlockCellNx) * Cell_Len),
                                     T((Nj / 2 + Ni / BlockCellNx) * Cell_Len)));

  // geometry helper
  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, cavity, Cell_Len, Ni / BlockCellNx);
  GeoHelper.forEachBlockCell([&](BasicBlock<T, 2>& block) {
    if (isOverlapped(block, seedcavity)) {
      block.refine();
    }
  });
  GeoHelper.forEachBlockCell([&](BasicBlock<T, 2>& block) {
    if (isOverlapped(block, seedcavity)) {
      block.refine();
    }
  });
  GeoHelper.CheckRefine();
  GeoHelper.CreateBlocks();
  GeoHelper.AdaptiveOptimization(Thread_Num);
  GeoHelper.LoadBalancing();

  // geometry
  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(cavity,
                 [&](auto& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet0>(cavity, BouncebackFlag);

  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define lattice ------------------
  using NSFIELDS = TypePack<RHO<T>, VELOCITY<T, 2>, POP<T, LatSet0::q>, SCALARFORCE<T>>;
  ValuePack NSInitValues(BaseConv.getLatRhoInit(), Vector<T, 2>{}, T{}, T{});
  using NSCELL = Cell<T, LatSet0, NSFIELDS>;
  BlockLatticeManager<T, LatSet0, NSFIELDS> NSLattice(Geo, NSInitValues, BaseConv);

  using CONCFIELDS = TypePack<CONC<T>, POP<T, LatSet1::q>, RHOINIT<T>, GBETA<T>>;
  using CONCFIELDREFS = TypePack<VELOCITY<T, 2>, CA::EXCESSC<T>>;
  BlockFieldManager<CA::EXCESSC<T>, T, 2>* tempExcessC = nullptr;
  using CONCFIELDPACK = TypePack<CONCFIELDS, CONCFIELDREFS>;
  ValuePack CONCInitValues(ConcConv.getLatRhoInit(), T{}, ConcConv.getLatRhoInit(),
                           ConcConv.getLattice_gbeta());
  using CONCCELL = Cell<T, LatSet1, ExtractFieldPack<CONCFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet1, CONCFIELDPACK> SOLattice(
    Geo, CONCInitValues, ConcConv, &NSLattice.getField<VELOCITY<T, 2>>(), tempExcessC);

  using TEMPFIELDS = TypePack<TEMP<T>, POP<T, LatSet1::q>, RHOINIT<T>, GBETA<T>>;
  using TEMPFIELDREFS = TypePack<VELOCITY<T, 2>>;
  using TEMPFIELDPACK = TypePack<TEMPFIELDS, TEMPFIELDREFS>;
  ValuePack TEMPInitValues(TempConv.getLatRhoInit(), T{}, TempConv.getLatRhoInit(),
                           TempConv.getLattice_gbeta());
  using TEMPCELL = Cell<T, LatSet1, ExtractFieldPack<TEMPFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet1, TEMPFIELDPACK> THLattice(
    Geo, TEMPInitValues, TempConv, &NSLattice.getField<VELOCITY<T, 2>>());

  // --------------------- dynamic lattice ---------------------
  DynamicBlockLatticeHelper2D<T, LatSet0, NSFIELDS> NSDynLatHelper(NSLattice, GeoHelper,
                                                                   RefThold, CoaThold, 2);
  DynamicBlockLatticeHelper2D<T, LatSet1, CONCFIELDPACK> SODynLatHelper(
    SOLattice, GeoHelper, RefThold, CoaThold, 2);

  // --------------------- CA ---------------------
  ValuePack CAInitValues(CA::CAType::Boundary, T{}, T{}, T{}, T{}, T{}, T{},
                         TempConv.getLatRhoInit(), ConcConv.getLatRhoInit());
  CA::BlockZhuStefanescu2DManager<T, LatSetCA> CA(
    Geo, CAConv, Delta, pref_Orine, CAInitValues, &NSLattice.getField<VELOCITY<T, 2>>(),
    &SOLattice.getField<CONC<T>>(), &THLattice.getField<TEMP<T>>());

  SOLattice.template addField<CA::EXCESSC<T>>(CA.template getField<CA::EXCESSC<T>>());
  // set CA State field
  CA.getField<CA::STATE>().forEach(
    FlagFM, AABBFlag | BouncebackFlag,
    [&](auto& field, std::size_t id) { field.SetField(id, CA::CAType::Fluid); });
  CA.Setup(Geo.getIndex(Vector<int, 2>{Ni / 2, Nj / 2}));

  // --------------------- BCs ---------------------
  // NS
  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,
                            BlockLatticeManager<T, LatSet0, NSFIELDS>,
                            BlockFieldManager<FLAG, T, 2>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  BBLikeMovingBlockBdManager<bounceback::normal<NSCELL>,
                             BlockLatticeManager<T, LatSet0, NSFIELDS>,
                             BlockFieldManager<FLAG, T, 2>>
    NS_MBB("NS_MBB", NSLattice, CA.getInterfaces(), FlagFM, CA::CAType::Solid);

  // Conc
  BBLikeFixedBlockBdManager<bounceback::normal<CONCCELL>,
                            BlockLatticeManager<T, LatSet1, CONCFIELDPACK>,
                            BlockFieldManager<FLAG, T, 2>>
    SO_BB("SO_BB", SOLattice, FlagFM, BouncebackFlag, VoidFlag);

  BBLikeMovingBlockBdManager<bounceback::normal<CONCCELL>,
                             BlockLatticeManager<T, LatSet1, CONCFIELDPACK>,
                             BlockFieldManager<FLAG, T, 2>>
    SO_MBB("SO_MBB", SOLattice, CA.getInterfaces(), FlagFM, CA::CAType::Solid);


  // define task/ dynamics:
  // NS task
  // bulk task
  using NSBulkTask =
    tmp::Key_TypePair<CA::CAType::Fluid,
                      collision::BGKForce_Feq_RhoU<equilibrium::SecondOrder<NSCELL>,
                                                   force::ScalarForce<NSCELL>, true>>;
  // wall task
  using NSWallTask =
    tmp::Key_TypePair<CA::CAType::Interface,
                      collision::BGKForce_Feq<equilibrium::SecondOrder<NSCELL>,
                                              force::ScalarForce<NSCELL>>>;

  using NSTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSBulkTask, NSWallTask>;

  // SO task
  using SOTask =
    tmp::Key_TypePair<FI_Flag,
                      collision::BGKSource_Feq_Rho<equilibrium::SecondOrder<CONCCELL>,
                                                   CA::EXCESSC<T>, true>>;

  using SOTaskSelector = TaskSelector<std::uint8_t, CONCCELL, SOTask>;

  // buoyancy
  using SObuoyancyTask = tmp::Key_TypePair<FI_Flag, force::Buoyancy<NSCELL, CONCCELL>>;
  using SObuoyancyTaskSelector =
    CoupledTaskSelector<std::uint8_t, NSCELL, CONCCELL, SObuoyancyTask>;
  BlockLatManagerCoupling SObuoyancy(NSLattice, SOLattice);
  using THbuoyancyTask = tmp::Key_TypePair<FI_Flag, force::Buoyancy<NSCELL, TEMPCELL>>;
  using THbuoyancyTaskSelector =
    CoupledTaskSelector<std::uint8_t, NSCELL, TEMPCELL, THbuoyancyTask>;
  BlockLatManagerCoupling THbuoyancy(NSLattice, THLattice);

  // writer
  vtmo::ScalarWriter CWriter("Conc", SOLattice.getField<CONC<T>>());
  vtmo::ScalarWriter StateWriter("State", CA.getField<CA::STATE>());
  vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::vtmWriter<T, 2> MainWriter("cazsblock2d", Geo, 1);
  MainWriter.addWriterSet(CWriter, StateWriter, VecWriter);

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  MainWriter.WriteBinary(MainLoopTimer());

  while (MainLoopTimer() < MaxStep) {
    // clear force
    NSLattice.getField<SCALARFORCE<T>>().InitValue(T{});
    // get buoyancy
    SObuoyancy.ApplyCellDynamics<SObuoyancyTaskSelector>(MainLoopTimer(),
                                                         CA.getField<CA::STATE>());
    THbuoyancy.ApplyCellDynamics<THbuoyancyTaskSelector>(MainLoopTimer(),
                                                         CA.getField<CA::STATE>());

    // NS task
    NSLattice.ApplyCellDynamics<NSTaskSelector>(MainLoopTimer(),
                                                CA.getField<CA::STATE>());
    NSLattice.Stream(MainLoopTimer());
    NS_BB.Apply(MainLoopTimer());
    NS_MBB.Apply(MainLoopTimer());
    NSLattice.Communicate(MainLoopTimer());

    // SO task
    SOLattice.ApplyCellDynamics<SOTaskSelector>(MainLoopTimer(),
                                                CA.getField<CA::STATE>());
    SOLattice.Stream(MainLoopTimer());
    SO_BB.Apply(MainLoopTimer());
    SO_MBB.Apply(MainLoopTimer());
    SOLattice.Communicate(MainLoopTimer());

    CA.Apply_SimpleCapture();


    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      // Velocity and Conc Field Communication
      NSLattice.getField<VELOCITY<T, 2>>().CommunicateAll();
      SOLattice.getField<CONC<T>>().CommunicateAll();
      CA.Communicate();

      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Print<std::size_t>("Interface", CA.getInterfaceNum());
      Printer::Print<T>("Solid%", T(CA.getSolidCount()) * 100 / Geo.getTotalCellNum());
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());
    }

    // ----- adaptive mesh refinement -----
    if (MainLoopTimer() % RefineCheckStep == 0) {
      if (CA.WillRefineBlockCells(GeoHelper)) {
        SODynLatHelper.GeoRefine(Thread_Num);

        // Geo Init
        Geo.Init(GeoHelper);

        // field reconstruction
        FlagFM.Init(VoidFlag);
        FlagFM.forEach(cavity, [&](FlagField& field, std::size_t id) {
          field.SetField(id, AABBFlag);
        });
        FlagFM.template SetupBoundary<LatSet0>(cavity, BouncebackFlag);

        CA.getField<CA::STATE>().InitCopy(
          GeoHelper, CA::CAType::Boundary, FlagFM, AABBFlag,
          [&](auto& field, std::size_t id) { field.SetField(id, CA::CAType::Fluid); });
        CA.getField<CA::STATE>().NormalCommunicate();

        CA.CAFieldDataInit(GeoHelper);

        NSLattice.Init(GeoHelper, NSInitValues.values);
        SOLattice.Init(GeoHelper, CONCInitValues.values);
        THLattice.Init(GeoHelper, TEMPInitValues.values);

        NSDynLatHelper.PopFieldInit();
        SODynLatHelper.PopFieldInit();

        CA.Init();

        // Bcs init
        NS_BB.Init();
        NS_MBB.Init();
        SO_BB.Init();
        SO_MBB.Init();

        CWriter.Init(SOLattice.getField<CONC<T>>());
        StateWriter.Init(CA.getField<CA::STATE>());
        VecWriter.Init(NSLattice.getField<VELOCITY<T, 2>>());
        MainWriter.Init();
        MainWriter.addWriterSet(CWriter, StateWriter, VecWriter);
      }
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();
  return 0;
}