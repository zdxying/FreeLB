// dambreak2d.cpp

// Lid-driven cavity flow 2d


// the top wall is set with a constant velocity,
// while the other walls are set with a no-slip boundary condition
// Bounce-Back-like method is used:
// Bounce-Back-Moving-Wall method for the top wall
// Bounce-Back method for the other walls

#include "freelb.h"
#include "freelb.hh"
#include "lbm/freeSurface.h"
#include "lbm/freeSurface.hh"

// int Total_Macro_Step = 0;
using T = FLOAT;
using LatSet = D3Q27<T>;

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
int Nk;
T Cell_Len;
T RT;
int Thread_Num;

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
Vector<T, 3> U_Ini;  // mm/s
T U_Max;

// bcs
Vector<T, 3> U_Wall;  // mm/s

// power-law
T BehaviorIndex;
T MInViscCoef;
T MaxViscCoef;

// Simulation settings
int MaxStep;
int OutputStep;

std::string work_dir;

void readParam() {
  iniReader param_reader("dambreak3d.ini");
  // Thread_Num = param_reader.getValue<int>("OMP", "Thread_Num");
  // mesh
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Nk = param_reader.getValue<int>("Mesh", "Nk");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
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
  U_Ini[2] = param_reader.getValue<T>("Init_Conditions", "U_Ini2");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  // bcs
  U_Wall[0] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall1");
  U_Wall[2] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall2");
  // LB
  RT = param_reader.getValue<T>("LB", "RT");
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");


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

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // define geometry
  AABB<T, 3> cavity(Vector<T, 3>{},
                    Vector<T, 3>(T(Ni * Cell_Len), T(Nj * Cell_Len), T(Nk * Cell_Len)));
  AABB<T, 3> fluid(Vector<T, 3>{},
                   Vector<T, 3>(T(int(Ni / 2) * Cell_Len), T(int(Nj) * Cell_Len), T(int(Nk/2) * Cell_Len)));
  BlockGeometry3D<T> Geo(Ni, Nj, Nk, Thread_Num, cavity, Cell_Len);

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
  using NSFIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>, SCALARCONSTFORCE<T>>;

  // using ALLFIELDS = MergeFieldPack<NSFIELDS, FS::FSFIELDS<T, LatSet>, FS::FSPARAMS<T>>::mergedpack;
  using ALLNSFS_FIELDS = MergeFieldPack<NSFIELDS, FS::FSFIELDS<T, LatSet>, FS::FSPARAMS<T>>::mergedpack;
  using ALLFIELDS = MergeFieldPack<ALLNSFS_FIELDS, PowerLawPARAMS<T>>::mergedpack;

  // a conversion factor of unit s^2 / g
  // [surface_tension_coefficient_factor * surface_tension_coefficient] = [1]
  // (LatRT_ - T(0.5)) * cs2 * deltaX_ * deltaX_ / VisKine_

  T surface_tension_coefficient_factor =
    BaseConv.Conv_Time * BaseConv.Conv_Time / (rho_ref * std::pow(BaseConv.Conv_L, 3));

  ValuePack NSInitValues(BaseConv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{}, -BaseConv.Lattice_g);
  ValuePack FSInitValues(FS::FSType::Solid, T{}, T{}, T{});
  ValuePack FSParamsInitValues(LonelyThreshold, VOF_Trans_Threshold, true, surface_tension_coefficient_factor* surface_tension_coefficient);
  // power-law dynamics for non-Newtonian fluid
  ValuePack PowerLawInitValues(BaseConv.Lattice_VisKine, BehaviorIndex - 1, BaseConv.Lattice_VisKine*MInViscCoef, BaseConv.Lattice_VisKine*MaxViscCoef);

  // auto ALLValues = mergeValuePack(NSInitValues, FSInitValues, FSParamsInitValues);
  auto ALLNSFSValues = mergeValuePack(NSInitValues, FSInitValues, FSParamsInitValues);
  auto ALLValues = mergeValuePack(ALLNSFSValues, PowerLawInitValues);

  using NSCELL = Cell<T, LatSet, ALLFIELDS>;
  using NSLAT = BlockLatticeManager<T, LatSet, ALLFIELDS>;
  using NSBlockLat = BlockLattice<T, LatSet, ALLFIELDS>;
  BlockLatticeManager<T, LatSet, ALLFIELDS> NSLattice(Geo, ALLValues, BaseConv);

  //// free surface

  // set cell state
  NSLattice.getField<FS::STATE>().forEach(
    cavity, [&](auto& field, std::size_t id) { field.SetField(id, FS::FSType::Gas); });
  // set fluid
  NSLattice.getField<FS::STATE>().forEach(
    fluid, [&](auto& field, std::size_t id) { field.SetField(id, FS::FSType::Fluid); });

  NSLattice.getField<FS::STATE>().template SetupBoundary<LatSet>(cavity, FS::FSType::Wall);

  FS::FreeSurfaceHelper<NSLAT>::Init(NSLattice);

  //// end free surface

  // define task/ dynamics:
  // NS task
  using NSBulkTask =
    tmp::Key_TypePair<FS::FSType::Fluid | FS::FSType::Interface,
                      collision::BGKForce_Feq_RhoU<equilibrium::SecondOrder<NSCELL>,
                                                   force::ScalarConstForce<NSCELL>, true>>;
  using NSWallTask = tmp::Key_TypePair<FS::FSType::Wall, collision::BounceBack<NSCELL>>;

  using NSTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSBulkTask, NSWallTask>;

  // bcs
  // BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,
  //                         BlockLatticeManager<T, LatSet, ALLFIELDS>,
  //                         BlockFieldManager<FLAG, T, LatSet::d>>
  // NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  vtmo::ScalarWriter rhovtm("rho", NSLattice.getField<RHO<T>>());
  vtmo::ScalarWriter MassWriter("Mass", NSLattice.getField<FS::MASS<T>>());
  vtmo::VectorWriter VeloWriter("Velo", NSLattice.getField<VELOCITY<T, LatSet::d>>());
  vtmo::ScalarWriter VOFWriter("VOF", NSLattice.getField<FS::VOLUMEFRAC<T>>());
  vtmo::ScalarWriter StateWriter("State", NSLattice.getField<FS::STATE>());
  vtmo::vtmWriter<T, LatSet::d> Writer("dambreak3d", Geo, 1);
  Writer.addWriterSet(StateWriter, VeloWriter);

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  Writer.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  while (MainLoopTimer() < MaxStep) {
    ++MainLoopTimer;
    ++OutputTimer;

    NSLattice.ApplyCellDynamics<NSTaskSelector>(MainLoopTimer(), NSLattice.getField<FS::STATE>());
    NSLattice.Stream(MainLoopTimer());
    // NS_BB.Apply(MainLoopTimer());
    NSLattice.Communicate(MainLoopTimer());

    FS::FreeSurfaceApply<BlockLatticeManager<T, LatSet, ALLFIELDS>>::Apply(NSLattice, MainLoopTimer());


    if (MainLoopTimer() % OutputStep == 0) {
      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      Writer.WriteBinary(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());

  return 0;
}