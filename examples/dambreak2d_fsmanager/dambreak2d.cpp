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
using LatSet = D2Q9<T>;

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
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
Vector<T, 2> U_Ini;  // mm/s
T U_Max;

// bcs
Vector<T, 2> U_Wall;  // mm/s

// Simulation settings
int MaxStep;
int OutputStep;

std::string work_dir;

void readParam() {
  iniReader param_reader("dambreak2d.ini");
  // Thread_Num = param_reader.getValue<int>("OMP", "Thread_Num");
  // mesh
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  // physical properties
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");

  /*free surface*/
  surface_tension_coefficient =
    param_reader.getValue<T>("Free_Surface", "surface_tension_coefficient");
  VOF_Trans_Threshold = param_reader.getValue<T>("Free_Surface", "VOF_Trans_Threshold");
  LonelyThreshold = param_reader.getValue<T>("Free_Surface", "LonelyThreshold");
  // init conditions
  U_Ini[0] = param_reader.getValue<T>("Init_Conditions", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("Init_Conditions", "U_Ini1");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  // bcs
  U_Wall[0] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall1");
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
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> fluid(Vector<T, 2>(T(0), T(0)),
                   Vector<T, 2>(T(int(Ni / 2) * Cell_Len), T(int(Nj / 2) * Cell_Len)));
  BlockGeometry2D<T> Geo(Ni, Nj, 1, cavity, Cell_Len);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(cavity,
                 [&](auto& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);

  vtmo::ScalarWriter Flagvtm("Flag", FlagFM);
  vtmo::vtmWriter<T, LatSet::d> FlagWriter("Flag", Geo, 1);
  FlagWriter.addWriterSet(Flagvtm);
  FlagWriter.WriteBinary();

  // ------------------ define lattice ------------------
  using NSFIELDS = TypePack<RHO<T>, VELOCITY<T, 2>, POP<T, LatSet::q>, SCALARCONSTFORCE<T>>;
  ValuePack NSInitValues(BaseConv.getLatRhoInit(), Vector<T, 2>{}, T{}, -BaseConv.Lattice_g);
  using NSCELL = Cell<T, LatSet, NSFIELDS>;
  BlockLatticeManager<T, LatSet, NSFIELDS> NSLattice(Geo, NSInitValues, BaseConv);

  //// free surface
  // a conversion factor of unit s^2 / g
  // [surface_tension_coefficient_factor * surface_tension_coefficient] = [1]
  // (LatRT_ - T(0.5)) * cs2 * deltaX_ * deltaX_ / VisKine_
  T surface_tension_coefficient_factor =
    BaseConv.Conv_Time * BaseConv.Conv_Time / (rho_ref * std::pow(BaseConv.Conv_L, 3));

  ValuePack FSInitValues(FS::FSType::Solid, T{}, T{}, T{});

  FS::FreeSurface2DManager<T, LatSet, NSFIELDS> FreeSurface(NSLattice, FSInitValues, LonelyThreshold, VOF_Trans_Threshold);
  // set cell state
  FreeSurface.getField<FS::STATE>().forEach(
    cavity, [&](auto& field, std::size_t id) { field.SetField(id, FS::FSType::Gas); });
  // set fluid
  FreeSurface.getField<FS::STATE>().forEach(
    fluid, [&](auto& field, std::size_t id) { field.SetField(id, FS::FSType::Fluid); });
  // set interface
  FreeSurface.Init();
  //// end free surface

  // define task/ dynamics:
  // NS task
  using NSBulkTask =
    tmp::Key_TypePair<FS::FSType::Fluid | FS::FSType::Interface,
                      collision::BGKForce_Feq_RhoU<equilibrium::SecondOrder<NSCELL>,
                                                   force::ScalarConstForce<NSCELL>, true>>;
  
  using NSTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSBulkTask>;

  // bcs
    BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,
                            BlockLatticeManager<T, LatSet, NSFIELDS>,
                            BlockFieldManager<FLAG, T, 2>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  vtmo::ScalarWriter rhovtm("rho", NSLattice.getField<RHO<T>>());
  vtmo::ScalarWriter MassWriter("Mass", FreeSurface.getField<FS::MASS<T>>());
  vtmo::VectorWriter VeloWriter("Velo", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter VOFWriter("VOF", FreeSurface.getField<FS::VOLUMEFRAC<T>>());
  vtmo::ScalarWriter StateWriter("State", FreeSurface.getField<FS::STATE>());
  vtmo::vtmWriter<T, LatSet::d> Writer("dambreak2d", Geo, 1);
  Writer.addWriterSet(rhovtm, MassWriter, VOFWriter, StateWriter, VeloWriter);

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  Writer.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  while (MainLoopTimer() < MaxStep) {
    ++MainLoopTimer;
    ++OutputTimer;

    NSLattice.ApplyCellDynamics<NSTaskSelector>(MainLoopTimer(), FreeSurface.getField<FS::STATE>());

    NSLattice.Stream(MainLoopTimer());
    NS_BB.Apply(MainLoopTimer());

    FreeSurface.Apply();

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