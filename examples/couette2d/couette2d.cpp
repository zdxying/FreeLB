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

// couette2d.cpp
// Couette flow 2d
// this is a benchmark for the freeLB library

// the top wall is set with a constant velocity with Bounce-Back-Moving-Wall
// the inlet and outlet are set with a Periodic boundary condition
// the bottom wall is set with a no-slip boundary condition with Bounce-Back

// note that the four corners are set with Bounce-Back-Like method
// do not set with periodic boundary condition

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
/*---------------------
        physical param
---------------------*/

// physical properties
T rho_ref;    // g/mm^3
T Dyna_Visc;  // Pa·s Dynamic viscosity of the liquid
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
T Ra;         // Rayleigh number
// init conditions
Vector<T, 2> U_Ini;  // mm/s
T U_Max;

// bcs
Vector<T, 2> U_Wall;  // mm/s

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  iniReader param_reader("Couetteparam.ini");
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
  Dyna_Visc = param_reader.getValue<T>("Physical_Property", "Dyna_Visc");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");
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
  tol = param_reader.getValue<T>("tolerance", "tol");


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
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t PeriodicFlag1 = std::uint8_t(8);
  std::uint8_t PeriodicFlag2 = std::uint8_t(16);
  std::uint8_t BBMovingWallFlag = std::uint8_t(32);

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni, U_Max, RT);
  // Conv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> left(Vector<T, 2>(T(0), Cell_Len),
                  Vector<T, 2>(Cell_Len, T((Nj - 1) * Cell_Len)));
  AABB<T, 2> right(Vector<T, 2>(T((Ni - 1) * Cell_Len), Cell_Len),
                   Vector<T, 2>(T(Ni * Cell_Len), T((Nj - 1) * Cell_Len)));
  AABB<T, 2> toplid(Vector<T, 2>(T(0), T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  Geometry2D<T> Geo(Ni, Nj, cavity, Cell_Len);
  Geo.SetupBoundary<LatSet>();
  Geo.setFlag(left, BouncebackFlag, PeriodicFlag1);
  Geo.setFlag(right, BouncebackFlag, PeriodicFlag2);
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
  BBLikeFixedBoundary<T, LatSet, BounceBackLikeMethod<T, LatSet>::normal_bounceback>
    NS_BB("NS_BB", NSLattice, BouncebackFlag);
  BBLikeFixedBoundary<T, LatSet, BounceBackLikeMethod<T, LatSet>::movingwall_bounceback>
    NS_BBMW("NS_BBMW", NSLattice, BBMovingWallFlag);
  FixedPeriodicBoundary<T, LatSet> NS_Pleft(NSLattice, left, right, PeriodicFlag1);
  FixedPeriodicBoundary<T, LatSet> NS_Pright(NSLattice, right, left, PeriodicFlag2);
  BoundaryManager BM1(&NS_BB, &NS_BBMW);
  BoundaryManager BM2(&NS_Pleft, &NS_Pright);

  T res = 1;

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  vtkWriter::FieldScalarWriter<T> RhoWriter("rho",
                                            NSLattice.getRhoField().getField().getdata(),
                                            NSLattice.getRhoField().getField().size());
  vtkWriter::FieldVectorWriter_AOS<T, LatSet::d> VelocityWriter(
    "velocity", NSLattice.getVelocityField().getField().getdata(),
    NSLattice.getVelocityField().getField().size());
  vtkStruPointsWriter<T, LatSet::d> NSWriter("NS", Geo);
  NSWriter.addtoWriteList(&RhoWriter, &VelocityWriter);

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  // NSWriter.Write(MainLoopTimer());

  while (MainLoopTimer() < MaxStep && res > tol) {
    ++MainLoopTimer;
    ++OutputTimer;

    NSLattice.UpdateRho(NSLattice.getIndex());
    NSLattice.UpdateU(NSLattice.getInnerIndex());
    BM2.UpdateU();
    NSLattice.BGK<Equilibrium<T, LatSet>::SecondOrder>();
    BM2.Apply();
    NSLattice.Stream();
    BM1.Apply();

    if (MainLoopTimer() % OutputStep == 0) {
      res = NSLattice.getToleranceU();
      OutputTimer.Print_InnerLoopPerformance(NSLattice.getN(), OutputStep);
      Printer::Print_Res<T>(res);
      Printer::Endl();
      // NSWriter.Write(MainLoopTimer());
    }
  }
  NSWriter.Write(MainLoopTimer());
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(NSLattice.getN());

  return 0;
}