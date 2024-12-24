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

// dambreak2d.cpp


#include "freelb.h"
#include "freelb.hh"
#include "lbm/freeSurface.h"


using T = FLOAT;
using LatSet = D2Q9<T>;

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
T Cell_Len;
T RT;
int Thread_Num;
int BlockCellLen;

// physical properties
T rho_ref;    // g/mm^3
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
/*free surface*/
// surface tension: N/m = kg/s^2
// fluid: 0.0728 N/m at 20 C = 0.0728 kg/s^2 = 72.8 g/s^2
T surface_tension_coefficient;
// Anti jitter value
T VOF_Trans_Threshold;
// When to remove lonely cells
T LonelyThreshold;
// init conditions
Vector<T, 2> U_Ini;  // mm/s
T U_Max;

// bcs
Vector<T, 2> U_Wall;  // mm/s

// power-law
T BehaviorIndex;
T MInViscCoef;
T MaxViscCoef;

// LES
T Smagorinsky;

// Simulation settings
int MaxStep;
int OutputStep;


void readParam() {
  iniReader param_reader("dambreak2d.ini");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");
  // physical properties
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");

  /*free surface*/
  surface_tension_coefficient =
    param_reader.getValue<T>("Free_Surface", "surface_tension_coefficient");
  VOF_Trans_Threshold = param_reader.getValue<T>("Free_Surface", "VOF_Trans_Threshold");
  LonelyThreshold = param_reader.getValue<T>("Free_Surface", "LonelyThreshold");

  // power law
  BehaviorIndex = param_reader.getValue<T>("PowerLaw", "BehaviorIndex");
  MInViscCoef = param_reader.getValue<T>("PowerLaw", "MInViscCoef");
  MaxViscCoef = param_reader.getValue<T>("PowerLaw", "MaxViscCoef");

  // init conditions
  U_Ini[0] = param_reader.getValue<T>("Init_Conditions", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("Init_Conditions", "U_Ini1");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  // bcs
  U_Wall[0] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall1");
  // LB
  RT = param_reader.getValue<T>("LB", "RT");
  // LES
  Smagorinsky = param_reader.getValue<T>("LES", "Smagorinsky");
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  MPI_RANK(0)
  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:" << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
            << "----------------------------------------------" << std::endl;
}

int getOutputStep(const Timer& timer) {
  if (timer() > 1990) return 1;
  else return 1000;
}

int main(int argc, char* argv[]) {
  // constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  // constexpr std::uint8_t AABBFlag = std::uint8_t(2);
  // constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);

  mpi().init(&argc, &argv);

  MPI_DEBUG_WAIT

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // define geometry
  AABB<T, 2> cavity(Vector<T, 2>{},
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> fluid(Vector<T, 2>{},
                   Vector<T, 2>(T(int(Ni / 2) * Cell_Len), T(int(Nj / 2) * Cell_Len)));
  
  // method 1: dirctly define geometry, [serial][openmp]
  // BlockGeometry2D<T> Geo(Ni, Nj, Thread_Num, cavity, Cell_Len, 2);
  // end method 1

  // method 2: use geohelper for complex geometry, [mpi]
  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, cavity, Cell_Len, BlockCellLen, 2);
  GeoHelper.CreateBlocks();
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);
  // end method 2

  // ------------------ define flag field ------------------
  // BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  // FlagFM.forEach(cavity,
  //                [&](auto& field, std::size_t id) { field.SetField(id, AABBFlag); });
  // FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);

  // vtmo::ScalarWriter Flagvtm("Flag", FlagFM);
  // vtmo::vtmWriter<T, LatSet::d> FlagWriter("Flag", Geo, 1);
  // FlagWriter.addWriterSet(Flagvtm);
  // FlagWriter.WriteBinary();

  // ------------------ define lattice ------------------
  using NSFIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>, CONSTFORCE<T, LatSet::d>, OMEGA<T>, SMAGORINSKY<T>>;

  using ALLFIELDS = MergeFieldPack<NSFIELDS, olbfs::FSFIELDS<T, LatSet>, olbfs::FSPARAMS<T>>::mergedpack;
  // using ALLNSFS_FIELDS = MergeFieldPack<NSFIELDS, olbfs::FSFIELDS<T, LatSet>, olbfs::FSPARAMS<T>>::mergedpack;
  // using ALLFIELDS = MergeFieldPack<ALLNSFS_FIELDS, PowerLawPARAMS<T>>::mergedpack;

  // a conversion factor of unit s^2 / g
  // [surface_tension_coefficient_factor * surface_tension_coefficient] = [1]
  // (LatRT_ - T(0.5)) * cs2 * deltaX_ * deltaX_ / VisKine_

  T surface_tension_coefficient_factor =
    BaseConv.Conv_Time * BaseConv.Conv_Time / (rho_ref * std::pow(BaseConv.Conv_L, 3));

  ValuePack NSInitValues(BaseConv.getLatRhoInit(), Vector<T, 2>{}, T{}, Vector<T, 2>{T{}, -BaseConv.Lattice_g}, BaseConv.getOMEGA(), Smagorinsky);
  ValuePack FSInitValues(olbfs::FSType::Void, olbfs::FSFlag::None, T{}, T{}, Vector<T, LatSet::q>{}, Vector<T, 2>{});
  ValuePack FSParamsInitValues(LonelyThreshold, VOF_Trans_Threshold, true, surface_tension_coefficient_factor* surface_tension_coefficient);
  IF_MPI_RANK(0)
  {std::cout << "surface: " << surface_tension_coefficient_factor* surface_tension_coefficient << std::endl;}
  // power-law dynamics for non-Newtonian fluid
  // ValuePack PowerLawInitValues(BaseConv.Lattice_VisKine, BehaviorIndex - 1, BaseConv.Lattice_VisKine*MInViscCoef, BaseConv.Lattice_VisKine*MaxViscCoef);

  auto ALLValues = mergeValuePack(NSInitValues, FSInitValues, FSParamsInitValues);
  // auto ALLNSFSValues = mergeValuePack(NSInitValues, FSInitValues, FSParamsInitValues);
  // auto ALLValues = mergeValuePack(ALLNSFSValues, PowerLawInitValues);

  using NSCELL = Cell<T, LatSet, ALLFIELDS>;
  using NSLAT = BlockLatticeManager<T, LatSet, ALLFIELDS>;
  using NSBlockLat = BlockLattice<T, LatSet, ALLFIELDS>;
  using NSBlockLatMan = BlockLatticeManager<T, LatSet, ALLFIELDS>;
  BlockLatticeManager<T, LatSet, ALLFIELDS> NSLattice(Geo, ALLValues, BaseConv);

  //// free surface

  // set cell state
  NSLattice.getField<olbfs::STATE>().forEach(
    cavity, [&](auto& field, std::size_t id) { field.SetField(id, olbfs::FSType::Gas); });
  // set fluid
  NSLattice.getField<olbfs::STATE>().forEach(
    fluid, [&](auto& field, std::size_t id) { field.SetField(id, olbfs::FSType::Fluid); });

  NSLattice.getField<olbfs::STATE>().template SetupBoundary<LatSet>(cavity, olbfs::FSType::Wall);

  olbfs::FreeSurfaceHelper<NSLAT>::Init(NSLattice);

  //// end free surface

  // ------------------ define task/ dynamics ------------------
  // NS task  PowerLaw_BGKForce
  // openlb used BGK dynamics for Gas cells
  using NSBulkTask =
    tmp::Key_TypePair<olbfs::FSType::Fluid | olbfs::FSType::Interface | olbfs::FSType::Gas,
                      collision::SmagorinskyForceBGK<moment::forcerhoU<NSCELL, force::ConstForce<NSCELL>, true>, 
                      equilibrium::SecondOrder<NSCELL>, force::ConstForce<NSCELL>>>;
  using NSWallTask = tmp::Key_TypePair<olbfs::FSType::Wall, collision::BounceBack<NSCELL>>;

  using NSTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSBulkTask, NSWallTask>;

  using ForceRhoUTask = tmp::Key_TypePair<olbfs::FSType::Fluid | olbfs::FSType::Interface,
                                          moment::rhoU<NSCELL, true>>;
  using RhoUTaskSelector = TaskSelector<std::uint8_t, NSCELL, ForceRhoUTask>;

  // bcs
  // BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,
  //                         BlockLatticeManager<T, LatSet, ALLFIELDS>,
  //                         BlockFieldManager<FLAG, T, LatSet::d>>
  // NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  // ------------------ define writers ------------------
  vtmo::ScalarWriter rhovtm("rho", NSLattice.getField<RHO<T>>());
  vtmo::ScalarWriter MassWriter("Mass", NSLattice.getField<olbfs::MASS<T>>());
  vtmo::VectorWriter VeloWriter("Velo", NSLattice.getField<VELOCITY<T, LatSet::d>>());
  vtmo::ScalarWriter VOFWriter("VOF", NSLattice.getField<olbfs::VOLUMEFRAC<T>>());
  vtmo::ScalarWriter StateWriter("State", NSLattice.getField<olbfs::STATE>());
  vtmo::vtmWriter<T, LatSet::d> Writer("dambreak2d", Geo, 1);
  Writer.addWriterSet(rhovtm, MassWriter, VOFWriter, VeloWriter, StateWriter);

  FieldStatistics RhoStat(NSLattice.getField<RHO<T>>());
  FieldStatistics MassStat(NSLattice.getField<olbfs::MASS<T>>());

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  Writer.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  Printer::PrintTitle("Step: 0");
  Printer::Print("Average Rho", RhoStat.getAverage());
  Printer::Print("Average Mass", MassStat.getAverage());
  Printer::Print("Max Mass", MassStat.getMax());
  Printer::Print("Min Mass", MassStat.getMin());
  Printer::Endl();

  while (MainLoopTimer() < MaxStep) {
    ++MainLoopTimer;
    ++OutputTimer;

    NSLattice.ApplyCellDynamics<NSTaskSelector>(NSLattice.getField<olbfs::STATE>());
    NSLattice.Stream();
    // NS_BB.Apply(MainLoopTimer());
    NSLattice.NormalAllCommunicate();

    olbfs::FreeSurfaceApply<NSBlockLatMan>::Apply(NSLattice, MainLoopTimer());

    if (MainLoopTimer() % OutputStep == 0) {
      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Print("Average Rho", RhoStat.getAverage());
      Printer::Print("Average Mass", MassStat.getAverage());
      Printer::Print("Max Mass", MassStat.getMax());
      Printer::Print("Min Mass", MassStat.getMin());
      Printer::Endl();
      Writer.WriteBinary(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());

  return 0;
}