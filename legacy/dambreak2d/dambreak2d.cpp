// dambreak2d.cpp

// Lid-driven cavity flow 2d
// this is a benchmark for the freeLB library

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
T Ra;         // Rayleigh number
/*free surface*/
// surface tension: N/m = kg/s^2
// water: 0.0728 N/m at 20 C = 0.0728 kg/s^2 = 72.8 g/s^2
T surface_tension_coefficient;
// Anti jitter value
T transitionThreshold;
// When to remove lonely cells
T lonelyThreshold;
// init conditions
Vector<T, 2> U_Ini;  // mm/s
T U_Max;
T P_char;

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
  transitionThreshold =
      param_reader.getValue<T>("Free_Surface", "transitionThreshold");
  lonelyThreshold = param_reader.getValue<T>("Free_Surface", "lonelyThreshold");
  // init conditions
  U_Ini[0] = param_reader.getValue<T>("Init_Conditions", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("Init_Conditions", "U_Ini1");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  P_char = param_reader.getValue<T>("Init_Conditions", "P_char");
  // bcs
  U_Wall[0] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall1");
  // LB
  RT = param_reader.getValue<T>("LB", "RT");
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  
  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:"
            << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
            << "----------------------------------------------"
            << std::endl;
}

int main() {
  std::uint8_t BouncebackFlag = std::uint8_t(3);

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  // BaseConv.SimplifiedConvertFromViscosity(Ni, U_Max, Kine_Visc);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max,
                         Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // define geometry
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> water(
      Vector<T, 2>(T(0), T(0)),
      Vector<T, 2>(T(int(Ni / 4) * Cell_Len), T(int(1 * Nj / 3) * Cell_Len)));

  Geometry2D<T> Geo(Ni, Nj, cavity, Cell_Len, Vector<T, 2>{}, Type::Gas,
                    Type::Solid);
  Geo.SetupBoundary<LatSet>(Type::Gas, BouncebackFlag);
  Geo.setFlag(water, Type::Gas, Type::Fluid);

  vtkWriter::FieldFlagWriter<std::uint8_t> flagwriter(
      "flag", Geo.getGeoFlagField().getField().getdata(),
      Geo.getGeoFlagField().getField().size());
  vtkStruPointsWriter<T, LatSet::d> GeoWriter("dambreakGeo", Geo);
  GeoWriter.addtoWriteList(&flagwriter);
  GeoWriter.Write();

  // ------------------ define lattice ------------------
  // velocity field
  VectorFieldAOS<T, LatSet::d> Velocity(Geo.getVoxelsNum());

  // bcs
  BBLikeFixedBoundary<
      T, LatSet, BounceBackLikeMethod<T, LatSet>::normal_bounceback>
      NS_BB("NS_BB", Geo, BouncebackFlag);
  BoundaryManager<T, LatSet> NS_BM(&NS_BB);

  // lattice
  BasicLattice<T, LatSet> NSLattice(Geo, BaseConv, Velocity);
  // force 
  ConstBulkForce<T, LatSet> Gravity(NSLattice, Vector<T, LatSet::d>{T(0), -BaseConv.Lattice_g}, Velocity);

  //// free surface
  // a conversion factor of unit s^2 / g
  // [surface_tension_coefficient_factor * surface_tension_coefficient] = [1]
  // (LatRT_ - T(0.5)) * cs2 * deltaX_ * deltaX_ / VisKine_
  T surface_tension_coefficient_factor =
      BaseConv.Conv_Time * BaseConv.Conv_Time /
      (rho_ref * std::pow(BaseConv.Conv_L, 3));
  FreeSurface2D<T, LatSet> FreeSurface(NSLattice);
  FreeSurface.setSrufaceTension(surface_tension_coefficient_factor *
                                surface_tension_coefficient);
  // set cell type for free surface
  // set solid
  Geo.forEachVoxel(Type::Solid, [&FreeSurface](int id) {
    FreeSurface.getType().SetField(id, Type::Solid);
  });
  // set gas
  Geo.forEachVoxel(Type::Gas, [&FreeSurface](int id) {
    FreeSurface.getType().SetField(id, Type::Gas);
  });
  // set fluid
  Geo.forEachVoxel(water, [&FreeSurface](int id) {
    FreeSurface.getType().SetField(id, Type::Fluid);
  });
  // set interface
  FreeSurface.getType().getField().for_isflag(
      Type::Fluid, [&FreeSurface, &Geo](std::size_t id) {
        if (Geo.template hasNeighborFlag<LatSet>(id, Type::Gas)) {
          FreeSurface.getType().SetField(id, Type::Interface);
        }
      });
  // set pops
  FreeSurface.getType().getField().for_isflag(
      Type::Fluid, [&NSLattice](std::size_t id) {
        NSLattice.InitPop(id, NSLattice.getLatRhoInit());
      });
  FreeSurface.getType().getField().for_isflag(
      Type::Interface, [&NSLattice](std::size_t id) {
        NSLattice.InitPop(id, NSLattice.getLatRhoInit() / 2);
      });
  FreeSurface.getType().getField().for_isflag(
      Type::Gas, [&NSLattice](std::size_t id) { NSLattice.InitPop(id, T(0)); });

  FreeSurface.Initialize();
  //// end free surface

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  // vtk writer
  vtkWriter::FieldScalerWriter<T> RhoWriter(
      "rho", NSLattice.getRhoField().getField().getdata(),
      NSLattice.getRhoField().getField().size());
  vtkWriter::FieldVectorWriter_AOS<T, LatSet::d> VelocityWriter(
      "velocity", NSLattice.getVelocityField().getField().getdata(),
      NSLattice.getVelocityField().getField().size());
  vtkWriter::FieldFlagWriter<Type> Typewriter(
      "type", FreeSurface.getType().getField().getdata(),
      FreeSurface.getType().getField().size());
  vtkWriter::FieldFlagWriter<Flag> Flagwriter(
      "flag", FreeSurface.getFlag().getField().getdata(),
      FreeSurface.getFlag().getField().size());
  vtkWriter::FieldScalerWriter<T> Masswriter(
      "mass", FreeSurface.getMass().getField().getdata(),
      FreeSurface.getMass().getField().size());
  vtkWriter::FieldScalerWriter<T> Fractionwriter(
      "mass", FreeSurface.getF().getField().getdata(),
      FreeSurface.getF().getField().size());
  vtkStruPointsWriter<T, LatSet::d> FSWriter("FS", Geo);
  FSWriter.addtoWriteList(&Typewriter, &Flagwriter, &RhoWriter,
                          &Masswriter, &Fractionwriter, &VelocityWriter);

  FSWriter.Write(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  while (MainLoopTimer() < MaxStep) {
    ++MainLoopTimer;
    ++OutputTimer;

    FreeSurface.Prepare();
    NSLattice.UpdateRho(NSLattice.getIndex());
    Gravity.BGK<Equilibrium<T, LatSet>::SecondOrder>(FreeSurface.getFluidIdx());
    Gravity.BGK<Equilibrium<T, LatSet>::SecondOrder>(FreeSurface.getInterfaceIdx());
    NSLattice.Stream();
    NS_BM.Apply(NSLattice);
    FreeSurface.Apply();

    if (MainLoopTimer() % OutputStep == 0) {
      OutputTimer.Print_InnerLoopPerformance(Ni * Nj, OutputStep);
      FSWriter.Write(MainLoopTimer());
    }
  }
  FSWriter.Write(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni * Nj);

  return 0;
}