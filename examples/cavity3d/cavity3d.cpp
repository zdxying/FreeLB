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

// cavity3d.cpp

// Lid-driven cavity flow 3d

// the top wall is set with a constant velocity,
// while the other walls are set with a no-slip boundary condition
// Bounce-Back-Moving-Wall method for the top wall
// Bounce-Back method for the other walls

#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;
using LatSet = D3Q19<T>;

/*----------------------------------------------
                Simulation Parameters
-----------------------------------------------*/
T RT;
// geometry
int Ni;
int Nj;
int Nk;
T Cell_Len;
int BlockCellLen;
int Thread_Num;
// physical properties
T rho_ref;    // g/mm^3
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

void readParam() {
  iniReader param_reader("cavity3d.ini");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  // mesh
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Nk = param_reader.getValue<int>("Mesh", "Nk");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");
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

  MPI_RANK(0)
  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:" << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
            << "Tolerance:         " << tol << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
            << "----------------------------------------------" << std::endl;
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t AABBFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
  constexpr std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  mpi().init(&argc, &argv);

  MPI_DEBUG_WAIT

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 3> cavity(
    Vector<T, 3>{}, Vector<T, 3>{T(Ni * Cell_Len), T(Nj * Cell_Len), T(Nk * Cell_Len)});
  AABB<T, 3> toplid(Vector<T, 3>{Cell_Len, Cell_Len, T((Nk - 1) * Cell_Len)},
    Vector<T, 3>{T((Ni - 1) * Cell_Len), T((Nj - 1) * Cell_Len), T(Nk * Cell_Len)});

  // method 1: dirctly define geometry, [serial][openmp]
  BlockGeometry3D<T> Geo(Ni, Nj, Nk, Thread_Num, cavity, Cell_Len);
  // end method 1

  // method 2: use geohelper for complex geometry, [mpi]
  // BlockGeometryHelper3D<T> GeoHelper(Ni, Nj, Nk, cavity, Cell_Len, BlockCellLen);
  // GeoHelper.CreateBlocks();
  // GeoHelper.AdaptiveOptimization(mpi().getSize());
  // GeoHelper.LoadBalancing(mpi().getSize());
  // BlockGeometry3D<T> Geo(GeoHelper);
  // end method 2

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(
    cavity, [&](FLAG& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
  FlagFM.forEach(toplid, [&](FLAG& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, BBMovingWallFlag);
  });
  // write flag field
  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define lattice ------------------
  // alias for collection of all fields
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  // alias for cell interface
  using CELL = Cell<T, LatSet, FIELDS>;
  // initial values for all fields
  ValuePack InitValues(BaseConv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  // lattice
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(Geo, InitValues, BaseConv);
  NSLattice.EnableToleranceU();
  T res = 1;

  // set initial value of field
  Vector<T, LatSet::d> LatU_Wall = BaseConv.getLatticeU(U_Wall);
  NSLattice.getField<VELOCITY<T, LatSet::d>>().forEach(toplid, FlagFM, BBMovingWallFlag,
    [&](auto& field, std::size_t id) { field.SetField(id, LatU_Wall); });

  // bcs
  BBLikeFixedBlockBdManager<bounceback::normal<CELL>,
    BlockLatticeManager<T, LatSet, FIELDS>, BlockFieldManager<FLAG, T, LatSet::d>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::movingwall<CELL>,
    BlockLatticeManager<T, LatSet, FIELDS>, BlockFieldManager<FLAG, T, LatSet::d>>
    NS_BBMW("NS_BBMW", NSLattice, FlagFM, BBMovingWallFlag, VoidFlag);
  BlockBoundaryManager BM(&NS_BB, &NS_BBMW);

  // ------------------ define task/ dynamics ------------------
  // bulk task
  using BulkTask = tmp::Key_TypePair<AABBFlag,
    collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  // BCs task as a collision process
  using BBTask = tmp::Key_TypePair<BouncebackFlag, collision::BounceBack<CELL>>;
  using BBMVTask =
    tmp::Key_TypePair<BBMovingWallFlag, collision::BounceBackMovingWall<CELL>>;
  // task collection
  using TaskCollection = tmp::TupleWrapper<BulkTask, BBTask, BBMVTask>;
  // task executor
  using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;

  // task: update rho and u
  using RhoUTask = tmp::Key_TypePair<AABBFlag, moment::rhoU<CELL, true>>;
  using TaskCollectionRhoU = tmp::TupleWrapper<RhoUTask>;
  using TaskSelectorRhoU = tmp::TaskSelector<TaskCollectionRhoU, std::uint8_t, CELL>;

  // ------------------ define writers ------------------
  // vtmo::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());
  vtmo::PhysScalarWriter physRhoWriter("physRho", NSLattice.getField<RHO<T>>(),
    std::bind(&BaseConverter<T>::getPhysRho, &BaseConv, std::placeholders::_1));
  // vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T,
  // LatSet::d>>());
  vtmo::PhysVectorWriter physVecWriter("physVelocity",
    NSLattice.getField<VELOCITY<T, LatSet::d>>(),
    std::bind(&BaseConverter<T>::getPhysU<LatSet::d>, &BaseConv, std::placeholders::_1));
  vtmo::vtmWriter<T, LatSet::d> NSWriter("cavity3d", Geo);
  NSWriter.addWriterSet(physRhoWriter, physVecWriter);

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  NSWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && res > tol) {
    NSLattice.ApplyCellDynamics<NSTask>(FlagFM);
    NSLattice.Stream();
    // BM.Apply(MainLoopTimer());
    NSLattice.NormalCommunicate();

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      NSLattice.ApplyCellDynamics<TaskSelectorRhoU>(FlagFM);

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