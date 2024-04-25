
#include "freelb.h"
#include "freelb.hh"

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

int Thread_Num;
/*---------------------
        physical param
---------------------*/

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
Vector<T, 2> U_Wall;  // mm/s
/*---------------------
        LB param
---------------------*/
T Th;  // char high temp
T Tl;  // char low temp

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  /*reader*/
  iniReader param_reader("ncparam2d.ini");
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
  Solutal_Expan_Coeff = param_reader.getValue<T>("Physical_Property", "Solutal_Expan_Coeff");
  Thermal_Expan_Coeff = param_reader.getValue<T>("Physical_Property", "Thermal_Expan_Coeff");
  SHeatCap_Liq = param_reader.getValue<T>("Physical_Property", "SHeatCap_Liq");
  SHeatCap_Soli = param_reader.getValue<T>("Physical_Property", "SHeatCap_Soli");
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
            << "tolerance:         " << tol << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
            << "----------------------------------------------" << std::endl;
}

int main() {
  std::uint8_t voidflag = std::uint8_t(1);
  std::uint8_t AABBflag = std::uint8_t(2);
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t leftFlag = std::uint8_t(8);
  std::uint8_t rightFlag = std::uint8_t(16);
  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet0::cs2);
  BaseConv.SimplifiedConvertFromViscosity(Ni - 2, U_Max, Kine_Visc);
  // Conv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni*Cell_Len, U_Max, Kine_Visc);
  TempConverter<T> TempConv(LatSet1::cs2, BaseConv, Temp_Ini);
  TempConv.ConvertTempFromTDiff_with_Ra(Tl, Th, TDiff, Ra);
  // Conv.ConvertTempFromSHeatCap_and_TCond_with_Ra(Tl, Th, T_Cond_Liq,
  // SHeatCap_Liq, Ra); Conv.ConvertTempFromSHeatCap_and_TCond_with_Texpan(Tl,
  // Th, T_Cond_Liq, SHeatCap_Liq, Thermal_Expan_Coeff);
  UnitConvManager<T> ConvManager(&BaseConv, &TempConv);
  ConvManager.Check_and_Print();

  // define geometry
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)), Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> leftside(Vector<T, 2>(T(0), T(0)), Vector<T, 2>(T(Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> rightside(Vector<T, 2>(T((Ni - 1) * Cell_Len), T(0)),
                       Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  // NS geometry
  Geometry2D<T> Geo0(Ni, Nj, cavity, Cell_Len, Vector<T, 2>{}, AABBflag, voidflag);
  Geo0.SetupBoundary<LatSet0>(AABBflag, BouncebackFlag);
  // thermal geometry
  Geometry2D<T> Geo1(Ni, Nj, cavity, Cell_Len, Vector<T, 2>{}, AABBflag, voidflag);
  Geo1.SetupBoundary<LatSet1>(AABBflag, BouncebackFlag);
  Geo1.setFlag(leftside, BouncebackFlag, leftFlag);
  Geo1.setFlag(rightside, BouncebackFlag, rightFlag);

  vtkWriter::FieldFlagWriter<std::uint8_t> geo0flagwriter(
    "flag0", Geo0.getGeoFlagField().getField().getdata(), Geo0.getGeoFlagField().getField().size());
  vtkWriter::FieldFlagWriter<std::uint8_t> geo1flagwriter(
    "flag1", Geo1.getGeoFlagField().getField().getdata(), Geo1.getGeoFlagField().getField().size());
  vtkStruPointsWriter<T, LatSet0::d> GeoWriter("CavGeo", Geo0);
  GeoWriter.addtoWriteList(&geo0flagwriter, &geo1flagwriter);
  GeoWriter.Write();

  // ------------------ define lattice ------------------
  // velocity field
  VectorFieldAOS<T, LatSet0::d> Velocity(Geo0.getVoxelsNum());
  // LATTICE
  BasicLattice<T, LatSet0> NSLattice(Geo0, BaseConv, Velocity);
  //   NSLattice.EnableToleranceU();
  BasicLattice<T, LatSet1> THLattice(Geo1, TempConv, Velocity);
  // bcs
  // NS BCS
  std::vector<int> NS_BCS_Index;
  Geo0.forEachVoxel(BouncebackFlag, [&NS_BCS_Index](int id) { NS_BCS_Index.push_back(id); });
  BBLikeFixedBoundary<T, LatSet0, BounceBackLikeMethod<T, LatSet0>::normal_bounceback> NS_BB(
    "NS_BB", NSLattice, BouncebackFlag);
  BoundaryManager NS_BM(&NS_BB);
  // thermal BCS
  BBLikeFixedBoundary<T, LatSet1, BounceBackLikeMethod<T, LatSet1>::normal_bounceback> TH_BB(
    "TH_BB", THLattice, BouncebackFlag);
  BBLikeFixedBoundary<T, LatSet1, BounceBackLikeMethod<T, LatSet1>::anti_bounceback_O2> TH_ABBL(
    "TH_ABBL", THLattice, leftFlag);
  BBLikeFixedBoundary<T, LatSet1, BounceBackLikeMethod<T, LatSet1>::anti_bounceback_O2> TH_ABBR(
    "TH_ABBR", THLattice, rightFlag);
  BoundaryManager TH_BM(&TH_BB, &TH_ABBL, &TH_ABBR);

  THLattice.EnableToleranceRho();
  Geo1.forEachVoxel(leftFlag, [&THLattice](int id) { THLattice.getRhoField().SetField(id, T(1)); });
  Geo1.forEachVoxel(rightFlag,
                    [&THLattice](int id) { THLattice.getRhoField().SetField(id, T(0)); });

  T res = 1;

  Buoyancy<T, LatSet0> Force(NSLattice, Velocity);
  Force.AddSource(&THLattice);

  vtkWriter::FieldScalerWriter<T> RhoWriter("rho", NSLattice.getRhoField().getField().getdata(),
                                            NSLattice.getRhoField().getField().size());
  vtkWriter::FieldScalerWriter<T> TWriter("T", THLattice.getRhoField().getField().getdata(),
                                          THLattice.getRhoField().getField().size());
  vtkWriter::FieldVectorWriter_AOS<T, LatSet0::d> VelocityWriter(
    "velocity", NSLattice.getVelocityField().getField().getdata(),
    NSLattice.getVelocityField().getField().size());
  vtkStruPointsWriter<T, LatSet0::d> NCWriter("FS", Geo0);
  NCWriter.addtoWriteList(&RhoWriter, &TWriter, &VelocityWriter);
  /*count and timer*/
  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  //   NCWriter.Write(MainLoopTimer());

  while (MainLoopTimer() < MaxStep && res > tol) {
    ++MainLoopTimer;
    ++OutputTimer;

    NSLattice.UpdateRho(NSLattice.getIndex());
    THLattice.UpdateRho(THLattice.getInnerIndex());
    TH_BB.UpdateRho();

    Force.GetBuoyancy();

    Force.BGK_U<Equilibrium<T, LatSet0>::SecondOrder>(NSLattice.getInnerIndex());
    Force.BGK<Equilibrium<T, LatSet0>::SecondOrder>(NS_BCS_Index);
    NSLattice.Stream();

    THLattice.BGK<Equilibrium<T, LatSet1>::SecondOrder>(THLattice.getIndex());
    THLattice.Stream();

    NS_BM.Apply();
    TH_BM.Apply();

    if (MainLoopTimer() % OutputStep == 0) {
      res = THLattice.getToleranceRho();
      OutputTimer.Print_InnerLoopPerformance(NSLattice.getN(), OutputStep);
      Printer::Print_Res<T>(res);
      //   NCWriter.Write(MainLoopTimer());
    }
  }

  NCWriter.Write(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Ni, Nj);

  return 0;
}
