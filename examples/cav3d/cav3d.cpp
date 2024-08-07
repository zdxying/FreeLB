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

// cav3d.cpp

// Lid-driven cavity flow 3d

// the top wall is set with a constant velocity,
// while the other walls are set with a no-slip boundary condition
// Bounce-Back-like method is used:
// Bounce-Back-Moving-Wall method for the top wall
// Bounce-Back method for the other walls

#include "freelb.h"
#include "freelb.hh"

// int Total_Macro_Step = 0;
using T = FLOAT;
using LatSet = D3Q19<T>;

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
// init conditions
Vector<T, 3> U_Ini;  // mm/s
T U_Max;

// bcs
Vector<T, 3> U_Wall;  // mm/s

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  // reader
  iniReader param_reader("cavity3d.ini");
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // mesh
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Nk = param_reader.getValue<int>("Mesh", "Nk");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  // physical properties
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");
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
  tol = param_reader.getValue<T>("tolerance", "tol");
  // parallel
#ifdef _OPENMP
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
#endif

  
  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:"
            << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
            << "Tolerance:         " << tol << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
            << "----------------------------------------------" << std::endl;
}

int main() {
  // std::uint8_t VoidFlag = std::uint8_t(1);
  std::uint8_t AABBFlag = std::uint8_t(2);
  std::uint8_t BouncebackFlag = std::uint8_t(4);
  std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 3> cavity(Vector<T, 3>{},
                    Vector<T, 3>(T(Ni * Cell_Len), T(Nj * Cell_Len), T(Nk * Cell_Len)));
  AABB<T, 3> toplid(
    Vector<T, 3>(Cell_Len, Cell_Len, T((Nk - 1) * Cell_Len)),
    Vector<T, 3>(T((Ni - 1) * Cell_Len), T((Nj - 1) * Cell_Len), T(Nk * Cell_Len)));

  Geometry3D<T> Geo(Ni, Nj, Nk, cavity, Cell_Len);

  Geo.SetupBoundary<LatSet>();
  Geo.setFlag(toplid, BouncebackFlag, BBMovingWallFlag);

  vtkWriter::FieldFlagWriter<std::uint8_t> flagwriter(
    "flag", Geo.getGeoFlagField().getField().getdataPtr(),
    Geo.getGeoFlagField().getField().size());
  vtkStruPointsWriter<T, LatSet::d> GeoWriter("CavGeo", Geo);
  GeoWriter.addtoWriteList(&flagwriter);
  GeoWriter.Write();

  // ------------------ define lattice ------------------
  // velocity field
  VectorFieldAOS<T, LatSet::d> Velocity(Geo.getVoxelsNum());
  // set initial value of field
  Vector<T, 3> LatU_Wall = BaseConv.getLatticeU(U_Wall);
  Geo.forEachVoxel(toplid, BBMovingWallFlag,
                   [&](int id) { Velocity.SetField(id, LatU_Wall); });
  // lattice
  PopLattice<T, LatSet> NSLattice(Geo, BaseConv, Velocity);
  NSLattice.EnableToleranceU();
  // bcs
  BBLikeFixedBoundary<T, LatSet, BounceBackLikeMethod<T, LatSet>::normal_bounceback>
    NS_BB("NS_BB", NSLattice, BouncebackFlag);
  BBLikeFixedBoundary<T, LatSet, BounceBackLikeMethod<T, LatSet>::movingwall_bounceback>
    NS_BBMW("NS_BBMW", NSLattice, BBMovingWallFlag);
  BoundaryManager BM(&NS_BB, &NS_BBMW);

  T res = 1;

  // writers

  // vtkWriter::FieldScalarWriter<T> RhoWriter("rho",
  //                                           NSLattice.getRhoField().getField().getdataPtr(),
  //                                           NSLattice.getRhoField().getField().size());
  // vtkWriter::FieldVectorWriter_AOS<T, LatSet::d> VelocityWriter(
  //   "velocity", NSLattice.getVelocityField().getField().getdataPtr(),
  //   NSLattice.getVelocityField().getField().size());
  // vtkStruPointsWriter<T, LatSet::d> NSWriter("NS", Geo);
  // NSWriter.addtoWriteList(&RhoWriter, &VelocityWriter);

  // vti writer
  vtiwriter::VectorWriter Velovti("Velocity", Velocity.getField());
  vtiwriter::ScalarWriter Rhovti("Rho", NSLattice.getRhoField().getField());
  vtiwriter::vtiManager NSvtim("NSvti", Geo.getVoxelSize(), Geo.getMin(),
                                 Vector<int, 3>{Ni + 1, Nj + 1, Nk + 1});
  NSvtim.addWriter(Rhovti, Velovti);

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  NSvtim.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  while (MainLoopTimer() < MaxStep && res > tol) {
    NSLattice.UpdateRho(Geo.getGeoFlagField().getField(), AABBFlag);
    NSLattice.UpdateU(Geo.getGeoFlagField().getField(), AABBFlag);
    NSLattice.BGK<Equilibrium<T, LatSet>::SecondOrder>(
      Geo.getGeoFlagField().getField(),
      std::uint8_t(AABBFlag | BouncebackFlag | BBMovingWallFlag));
    NSLattice.Stream();
    BM.Apply();

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      res = NSLattice.getTolU(1);
      OutputTimer.Print_InnerLoopPerformance(NSLattice.getN(), OutputStep);
      Printer::Print_Res<T>(res);
      Printer::Endl();
      NSvtim.WriteBinary(MainLoopTimer());
    }
  }

  // NSWriter.Write(MainLoopTimer());
  
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(NSLattice.getN());

  return 0;
}