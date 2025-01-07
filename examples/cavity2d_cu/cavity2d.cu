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

// cavblock3d.cpp

// Lid-driven cavity flow 3d
// this is a benchmark for the freeLB

// the top wall is set with a constant velocity,
// while the other walls are set with a no-slip boundary condition
// Bounce-Back-like method is used:
// Bounce-Back-Moving-Wall method for the top wall
// Bounce-Back method for the other walls

// block data structure is used

#include "freelb.h"
#include "freelb.hh"

using T = float;
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
T Dyna_Visc;  // Pa·s Dynamic viscosity of the liquid
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
// init conditions
Vector<T, LatSet::d> U_Ini;  // mm/s
T U_Max;

// bcs
Vector<T, LatSet::d> U_Wall;  // mm/s

// Simulation settings
int MaxStep;
int OutputStep;
T tol;
std::string work_dir;

void readParam() {
  iniReader param_reader("cavity2d.ini");
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
  std::cout << "[Simulation_Settings]:" << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
            << "Tolerance:         " << tol << "\n"
            << "----------------------------------------------" << std::endl;
}

int main() {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t AABBFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
  constexpr std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> toplid(Vector<T, 2>(Cell_Len, T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T((Ni - 1) * Cell_Len), T(Nj * Cell_Len)));
  BlockGeometry2D<T> Geo(Ni, Nj, 1, cavity, Cell_Len);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(cavity,
                 [&](FLAG& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
  FlagFM.forEach(toplid, [&](FLAG& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, BBMovingWallFlag);
  });
  // do not forget to copy to device
  FlagFM.copyToDevice();

  // vtmwriter::ScalarWriter FlagWriter("flag", FlagFM);
  // vtmwriter::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", Geo);
  // GeoWriter.addWriterSet(FlagWriter);
  // GeoWriter.WriteBinary();

  // ------------------ define lattice ------------------
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  using cudevFIELDS = typename ExtractCudevFieldPack<FIELDS>::cudev_pack;

  using HostCell = Cell<T, LatSet, FIELDS>;
  using CELL = cudev::Cell<T, LatSet, cudevFIELDS>;

  ValuePack InitValues(BaseConv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  // lattice
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(Geo, InitValues, BaseConv);
  NSLattice.EnableToleranceU();
  T res = 1;
  // set initial value of field
  Vector<T, 2> LatU_Wall = BaseConv.getLatticeU(U_Wall);
  NSLattice.getField<VELOCITY<T, LatSet::d>>().forEach(
    toplid, FlagFM, BBMovingWallFlag,
    [&](auto& field, std::size_t id) { field.SetField(id, LatU_Wall); }); 

  // bcs
  // BBLikeFixedBlockBdManager<bounceback::normal<CELL>, BlockLatticeManager<T, LatSet, FIELDS>, BlockFieldManager<FLAG, T, LatSet::d>>
  //   NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);
  // BBLikeFixedBlockBdManager<bounceback::movingwall<CELL>, BlockLatticeManager<T, LatSet, FIELDS>, BlockFieldManager<FLAG, T, LatSet::d>>
  //   NS_BBMW("NS_BBMW", NSLattice, FlagFM, BBMovingWallFlag, VoidFlag);
  // BlockBoundaryManager BM(&NS_BB, &NS_BBMW);

  // define task/ dynamics:
  // bulk task
  using BulkTask = tmp::Key_TypePair<AABBFlag, collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  // wall task
  using WallTask = tmp::Key_TypePair<BouncebackFlag | BBMovingWallFlag, collision::BGK<moment::useFieldrhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  // BCs task as a collision process, if used, bcs will be handled in the collision process
  using BBTask = tmp::Key_TypePair<BouncebackFlag, collision::BounceBack<CELL>>;
  using BBMVTask = tmp::Key_TypePair<BBMovingWallFlag, collision::BounceBackMovingWall<CELL>>;
  // task collection
  // using TaskCollection = tmp::TupleWrapper<BulkTask, WallTask>;
  using TaskCollection = tmp::TupleWrapper<BulkTask, BBTask, BBMVTask>;
  // task executor
  using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;

  // task: update rho and u
  using RhoUTask = tmp::Key_TypePair<AABBFlag, moment::rhoU<CELL, true>>;
  using TaskCollectionRhoU = tmp::TupleWrapper<RhoUTask>;
  using TaskSelectorRhoU = tmp::TaskSelector<TaskCollectionRhoU, std::uint8_t, CELL>;

  // writers
  vtmwriter::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());
  vtmwriter::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, LatSet::d>>());
  vtmwriter::vtmWriter<T, LatSet::d> NSWriter("cavblock2dcu", Geo);
  NSWriter.addWriterSet(RhoWriter, VecWriter);

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  std::cout << "Total Cells: " << Geo.getTotalCellNum() << std::endl;

  NSLattice.getField<POP<T, LatSet::q>>().copyToDevice();
  NSLattice.getField<RHO<T>>().copyToDevice();
  NSLattice.getField<VELOCITY<T, LatSet::d>>().copyToDevice();

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;
  
  NSWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && res > tol) {

    // NSLattice.ApplyCellDynamics<NSTask>(FlagFM);
    NSLattice.CuDevApplyCellDynamics<NSTask>(FlagFM);
    // NSLattice.CuDevApplyCellDynamics<collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>();
    
    // NSLattice.Stream();
    NSLattice.CuDevStream();
    // BM.Apply(MainLoopTimer());

    // NSLattice.Communicate(MainLoopTimer());

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      NSLattice.CuDevApplyCellDynamics<TaskSelectorRhoU>(FlagFM);
      cudaDeviceSynchronize();
      NSLattice.getBlockLat(0).getField<RHO<T>>().copyToHost();
      NSLattice.getBlockLat(0).getField<VELOCITY<T, LatSet::d>>().copyToHost();
      
      res = NSLattice.getToleranceU(-1);
      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Print_Res<T>(res);
      Printer::Endl();
      NSWriter.WriteBinary(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();

  return 0;
}