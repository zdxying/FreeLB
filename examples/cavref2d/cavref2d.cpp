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

// cavref2d.cpp

// Lid-driven cavity flow 2d with refined block structure

// the top wall is set with a constant velocity,
// while the other walls are set with a no-slip boundary condition
// Bounce-Back-like method is used:
// Bounce-Back-Moving-Wall method for the top wall
// Bounce-Back method for the other walls


#include "freelb.h"
#include "freelb.hh"


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
int BlockCellNx;

// physical properties
T rho_ref;    // g/mm^3
T Dyna_Visc;  // Pa·s Dynamic viscosity of the liquid
T Kine_Visc;  // mm^2/s kinematic viscosity of the liquid
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
  iniReader param_reader("cavityref2d.ini");
  // mesh
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  // parallel
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  // mesh
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellNx = param_reader.getValue<int>("Mesh", "BlockCellNx");
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
            << "Tolerance:         " << tol << "\n"
#ifdef _OPENMP
            << "Running on " << Thread_Num << " threads\n"
#endif
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
  // use 0.5 for refined block
  AABB<T, 2> toplid(Vector<T, 2>(Cell_Len / 2, T((Nj - 1) * Cell_Len)),
                    Vector<T, 2>(T((Ni - 0.5) * Cell_Len), T(Nj * Cell_Len)));
  // for refined block
  AABB<T, 2> innercavity(
    Vector<T, 2>(T(Ni * Cell_Len) / BlockCellNx, T(Nj * Cell_Len) / BlockCellNx),
    Vector<T, 2>(T(Ni * Cell_Len) * (BlockCellNx - 1) / BlockCellNx,
                 T(Nj * Cell_Len) * (BlockCellNx - 1) / BlockCellNx));

  // geometry helper
  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, cavity, Cell_Len, Ni / BlockCellNx);
  GeoHelper.forEachBlockCell([&](BasicBlock<T, 2>& block) {
    if (!isOverlapped(block, innercavity)) {
      block.refine();
    }
  });
  GeoHelper.CreateBlocks();
  GeoHelper.AdaptiveOptimization(Thread_Num);
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(cavity,
                 [&](auto& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
  FlagFM.forEach(toplid, [&](auto& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, BBMovingWallFlag);
  });

  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define lattice ------------------
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  using FIELDREFS = TypePack<FLAG>;
  using FIELDSPACK = TypePack<FIELDS, FIELDREFS>;
  using ALLFIELDS = ExtractFieldPack<FIELDSPACK>::mergedpack;
  using CELL = Cell<T, LatSet, ALLFIELDS>;
  ValuePack InitValues(BaseConv.getLatRhoInit(), Vector<T, 2>{}, T{});
  // lattice
  BlockLatticeManager<T, LatSet, FIELDSPACK> NSLattice(Geo, InitValues, BaseConv, &FlagFM);
  NSLattice.EnableToleranceU();
  T res = 1;
  // set initial value of field
  Vector<T, 2> LatU_Wall = BaseConv.getLatticeU(U_Wall);
  NSLattice.getField<VELOCITY<T, LatSet::d>>().forEach(FlagFM, BBMovingWallFlag,
                     [&](auto& field, std::size_t id) { field.SetField(id, LatU_Wall); });

  // bcs
  BBLikeFixedBlockBdManager<bounceback::normal<CELL>, BlockLatticeManager<T, LatSet, FIELDSPACK>, BlockFieldManager<FLAG, T, 2>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::movingwall<CELL>, BlockLatticeManager<T, LatSet, FIELDSPACK>, BlockFieldManager<FLAG, T, 2>>
    NS_BBMW("NS_BBMW", NSLattice, FlagFM, BBMovingWallFlag, VoidFlag);
  BlockBoundaryManager BM(&NS_BB, &NS_BBMW);

  // define task/ dynamics:
  // to use refined/multi-level block structure, macroscopic fields should be updated each time step for pop conversion
  // bulk task
  using BulkTask = tmp::Key_TypePair<AABBFlag, collision::BGK_Feq_RhoU<equilibrium::SecondOrder<CELL>, true>>;
  // wall task
  using WallTask = tmp::Key_TypePair<BouncebackFlag | BBMovingWallFlag, collision::BGK_Feq<equilibrium::SecondOrder<CELL>, true>>;
  // task collection
  using TaskCollection = tmp::TupleWrapper<BulkTask, WallTask>;
  // task executor
  using TaskSelector = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;
  // task: update rho and u
  // using RhoUTask = tmp::Key_TypePair<AABBFlag, moment::rhou<CELL>>;
  // using TaskCollectionRhoU = tmp::TupleWrapper<RhoUTask>;
  // using TaskSelectorRhoU = tmp::TaskSelector<TaskCollectionRhoU, std::uint8_t, CELL>;

  // writers
  vtmo::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());
  vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::vtmWriter<T, LatSet::d> NSWriter("cavref2d", Geo, 1);
  NSWriter.addWriterSet(RhoWriter, VecWriter);

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;
  NSWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && res > tol) {
    NSLattice.ApplyCellDynamics<TaskSelector>(MainLoopTimer(), FlagFM);
    // streaming
    NSLattice.Stream(MainLoopTimer());
    // boundary conditions
    BM.Apply(MainLoopTimer());
    // block communication
    NSLattice.Communicate(MainLoopTimer());

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      res = NSLattice.getToleranceU(-1);
      OutputTimer.Print_InnerLoopPerformance(GeoHelper.getTotalBaseCellNum(), OutputStep);
      Printer::Print_Res<T>(res);
      Printer::Endl();
      NSWriter.WriteBinary(MainLoopTimer());
    }
  }
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(GeoHelper.getTotalBaseCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();

  return 0;
}