// cavity2d.cpp

// Lid-driven cavity flow 2d
// this is a benchmark for the freeLB library

// the top wall is set with a constant velocity,
// while the other walls are set with a no-slip boundary condition
// Bounce-Back-like method is used:
// Bounce-Back-Moving-Wall method for the top wall
// Bounce-Back method for the other walls

#include "freelb.h"
#include "freelb.hh"

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

/*physical property*/
T rho_ref;    // g/mm^3
T Dyna_Visc;  // Pa·s Dynamic viscosity of the liquid
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
T Ra;         // Rayleigh number
/*init conditions*/
Vector<T, 2> U_Ini;  // mm/s
T U_Max;
T P_char;

/*bcs*/
Vector<T, 2> U_Wall;  // mm/s

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  /*reader*/
  iniReader param_reader("cavityparam2d.ini");
  // Thread_Num = param_reader.getValue<int>("OMP", "Thread_Num");
  /*mesh*/
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  /*CA mesh*/
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  /*physical property*/
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Dyna_Visc = param_reader.getValue<T>("Physical_Property", "Dyna_Visc");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");
  /*init conditions*/
  U_Ini[0] = param_reader.getValue<T>("Init_Conditions", "U_Ini0");
  U_Ini[1] = param_reader.getValue<T>("Init_Conditions", "U_Ini1");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  P_char = param_reader.getValue<T>("Init_Conditions", "P_char");
  /*bcs*/
  U_Wall[0] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall0");
  U_Wall[1] = param_reader.getValue<T>("Boundary_Conditions", "Velo_Wall1");
  // LB
  RT = param_reader.getValue<T>("LB", "RT");
  // Simulation settings
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");
  tol = param_reader.getValue<T>("tolerance", "tol");

  /*output to console*/
  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:"
            << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
            << "Tolerance:         " << tol << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
            << "----------------------------------------------"
            << std::endl;
}

void getDiffVelo(VectorFieldAOS<T, LatSet::d>& diff, const VectorFieldAOS<T, LatSet::d>& v1, const VectorFieldAOS<T, LatSet::d>& v2) {
  for (int i = 0; i < v1.getField().size(); ++i) {
    diff.SetField(i, v1.get(i) - v2.get(i));
  }
}

void getDiffRho(ScalerField<T>& diff, const ScalerField<T>& r1, const ScalerField<T>& r2) {
  for (int i = 0; i < r1.getField().size(); ++i) {
    diff.SetField(i, r1.get(i) - r2.get(i));
  }
}

int main() {
  std::uint8_t VoidFlag = std::uint8_t(1);
  std::uint8_t AABBFlag = std::uint8_t(2);
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConvertFromViscosity(Ni, U_Max, Kine_Visc);
  // BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max,
  // Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> toplid(Vector<T, 2>(Cell_Len, T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T((Ni-1) * Cell_Len),
                                 T(Nj * Cell_Len)));
  Geometry2D<T> Geo(Ni, Nj, cavity, Cell_Len);
  // boundary flag: 2 for no-slip, 3 for moving wall
  Geo.SetupBoundary<LatSet>();
  Geo.setFlag(toplid, BouncebackFlag, BBMovingWallFlag);

  vtkWriter::FieldFlagWriter<std::uint8_t> flagwriter(
      "flag", Geo.getGeoFlagField().getField().getdata(),
      Geo.getGeoFlagField().getField().size());
  vtkStruPointsWriter<T, LatSet::d> GeoWriter("CavGeo", Geo);
  GeoWriter.addtoWriteList(&flagwriter);
  GeoWriter.Write();

  // ------------------ define lattice ------------------
  // velocity field
  VectorFieldAOS<T, LatSet::d> Velocity(Geo.getVoxelsNum());
  // set initial value of field
  Geo.forEachVoxel(toplid, BBMovingWallFlag,
                     [&Velocity](int id) { Velocity.SetField(id, U_Wall); });
  // lattice
  BasicLattice<T, LatSet> NSLattice(Geo, BaseConv, Velocity);
  NSLattice.EnableToleranceU();
  // bcs
  BBLikeFixedBoundary<
      T, LatSet, BounceBackLikeMethod<T, LatSet>::normal_bounceback>
      NS_BB("NS_BB", NSLattice, BouncebackFlag);
  BBLikeFixedBoundary<
      T, LatSet, BounceBackLikeMethod<T, LatSet>::movingwall_bounceback>
      NS_BBMW("NS_BBMW", NSLattice, BBMovingWallFlag);
  BoundaryManager BM(&NS_BB, &NS_BBMW);

  T res = 1;
  // BasicLBM<T, LatSet> NS(NSLattice, BM);

  // --------------------lattice 2---------------------
  // velocity field
  // VectorFieldAOS<T, LatSet::d> Velocity2(Geo.getVoxelsNum());
  // Velocity2 = Velocity;
  // BasicLattice<T, LatSet> NSLattice2(Geo, BaseConv, Velocity2);
  // NSLattice2.EnableToleranceU();
  // // bcs
  // BBLikeFixedBoundary<
  //     T, LatSet, BounceBackLikeMethod<T, LatSet>::normal_bounceback>
  //     NS_BB2("NS_BB", NSLattice2, BouncebackFlag);
  // BBLikeFixedBoundary<
  //     T, LatSet, BounceBackLikeMethod<T, LatSet>::movingwall_bounceback>
  //     NS_BBMW2("NS_BBMW", NSLattice2, BBMovingWallFlag);
  // BoundaryManager BM2(&NS_BB2, &NS_BBMW2);

  // T res2 = 1;

  // // comparation
  // VectorFieldAOS<T, LatSet::d> DiffVelocity(Geo.getVoxelsNum(), Vector<T,LatSet::d>{});
  // ScalerField<T> DiffRho(Geo.getVoxelsNum(), 0);
  // --------------------lattice 2---------------------


  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  // vtkWriter::FieldScalerWriter<T> RhoWriter(
  //     "rho", NSLattice.getRhoField().getField().getdata(),
  //     NSLattice.getRhoField().getField().size());
  // vtkWriter::FieldVectorWriter_AOS<T, LatSet::d> VelocityWriter(
  //     "velocity", NSLattice.getVelocityField().getField().getdata(),
  //     NSLattice.getVelocityField().getField().size());
  // vtkStruPointsWriter<T, LatSet::d> NSWriter("NS", Geo);
  // NSWriter.addtoWriteList(&RhoWriter, &VelocityWriter);

  // vti writer
  vtiwriter::VectorWriter VeloWriter("Velocity", Velocity.getField());
  vtiwriter::ScalerWriter RhoWriter("Rho", NSLattice.getRhoField().getField());
  vtiwriter::vtiManager NSWriter("NS", Geo.getVoxelSize(), Geo.getMin(), Vector<int,2>{Ni+1,Nj+1});
  NSWriter.addWriter(&RhoWriter, &VeloWriter);
  NSWriter.WriteBinary(MainLoopTimer());
  
  // vtiwriter::VectorWriter DiffVeloWriter("DiffVelocity", DiffVelocity.getField());
  // vtiwriter::ScalerWriter DiffRhoWriter("DiffRho", DiffRho.getField());
  // vtiwriter::vtiManager VTIWriter("DiffV", Geo.getVoxelSize(), Geo.getMin(), Vector<int,2>{Ni+1,Nj+1});
  // VTIWriter.addWriter(&DiffRhoWriter,&DiffVeloWriter);
  // VTIWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && res > tol)
  {
    ++MainLoopTimer;
    ++OutputTimer;

    // NSLattice.UpdateRho(NSLattice.getIndex());
    // NSLattice.UpdateU(NSLattice.getInnerIndex());
    // NSLattice.BGK<Equilibrium<T, LatSet>::SecondOrder>();
    // NS.Apply<Equilibrium<T, LatSet>::Feq_secondOrder>();

    NSLattice.UpdateRho(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag));
    NSLattice.UpdateU(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag));
    // NSLattice.BGK<Equilibrium<T, LatSet>::SecondOrder>();
    NSLattice.BGK<Equilibrium<T, LatSet>::SecondOrder>(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag|BouncebackFlag|BBMovingWallFlag));
    NSLattice.Stream();
    BM.Apply();

    // test code
    // NSLattice2.UpdateRho(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag));
    // NSLattice2.UpdateU(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag));
    // // NSLattice2.BGK<Equilibrium<T, LatSet>::SecondOrder>();
    // NSLattice2.BGK<Equilibrium<T, LatSet>::SecondOrder>(Geo.getGeoFlagField().getField(), std::uint8_t(AABBFlag|BouncebackFlag|BBMovingWallFlag));
    // NSLattice2.Stream();
    // BM2.Apply();
    

    if (MainLoopTimer() % OutputStep == 0) {
      res = NSLattice.getTolU(1);
      // res2 = NSLattice2.getTolU(1);
      OutputTimer.Print_InnerLoopPerformance(NSLattice.getN(), OutputStep);
      Printer::Print_Res<T>(res);
      // Printer::Print_Res<T>(res2);
      // getDiffVelo(DiffVelocity, Velocity, Velocity2);
      // getDiffRho(DiffRho, NSLattice.getRhoField(), NSLattice2.getRhoField());
      NSWriter.WriteBinary(MainLoopTimer());
      // VTIWriter.WriteBinary(MainLoopTimer());
    }
  }
  NSWriter.WriteBinary(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(NSLattice.getN());

  return 0;
}

/*
UpdateRho(AABBFlag);
UpdateU(AABBFlag);
BGK() and BGK(AABBFlag|BouncebackFlag|BBMovingWallFlag) have the same result 

this is because the existance of the corner cells
in the present bcs method, the bcs of the direction where 2 neighbor void cells are in the same line is not considered
the corner cells will be "contaminated" by the void cells
*/
/*
update rho at boundary will lead to derease of rho and average rho

*/