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

// the top wall is set with a constant velocity with Bounce-Back-Moving-Wall
// the inlet and outlet are set with a Periodic boundary condition
// the bottom wall is set with a no-slip boundary condition with Bounce-Back

// note that the four corners are set with Bounce-Back-Like method
// do not set with periodic boundary condition

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
/*---------------------
        physical param
---------------------*/

// physical properties
T rho_ref;    // g/mm^3
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
  iniReader param_reader("couetteparam.ini");
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
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t AABBFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
  constexpr std::uint8_t BBMovingWallFlag = std::uint8_t(8);
  constexpr std::uint8_t PeriodicFlag = std::uint8_t(16);

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
  BlockGeometry2D<T> Geo(Ni, Nj, Thread_Num, cavity, Cell_Len);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(cavity,
                 [&](FLAG& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
  FlagFM.forEach(toplid, [&](FLAG& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, BBMovingWallFlag);
  });
  FlagFM.forEach(left, [&](FLAG& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, PeriodicFlag);
  });
  FlagFM.forEach(right, [&](FLAG& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, PeriodicFlag);
  });

  vtmwriter::ScalarWriter FlagWriter("flag", FlagFM);
  vtmwriter::vtmWriter<T, 2> GeoWriter("GeoFlag", Geo);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define lattice ------------------
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>, StrainRateMag<T>>;
  using CELL = Cell<T, LatSet, FIELDS>;
  ValuePack InitValues(BaseConv.getLatRhoInit(), Vector<T, 2>{}, T{}, T{});
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
  BBLikeFixedBlockBdManager<bounceback::normal<CELL>, BlockLatticeManager<T, LatSet, FIELDS>, BlockFieldManager<FLAG, T, 2>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::movingwall<CELL>, BlockLatticeManager<T, LatSet, FIELDS>, BlockFieldManager<FLAG, T, 2>>
    NS_BBMW("NS_BBMW", NSLattice, FlagFM, BBMovingWallFlag, VoidFlag);
  FixedPeriodicBoundaryManager<BlockLatticeManager<T, LatSet, FIELDS>, BlockFieldManager<FLAG, T, 2>>
    NS_Peri("NS_Peri", NSLattice, FlagFM, PeriodicFlag, VoidFlag);
  // direction should be handled carefully
  NS_Peri.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  BlockBoundaryManager BM(&NS_BB, &NS_BBMW);
  

  // define task/ dynamics:
  // bulk task
  using BulkTask = tmp::Key_TypePair<AABBFlag | PeriodicFlag, collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  // wall task
  using WallTask = tmp::Key_TypePair<BouncebackFlag | BBMovingWallFlag, collision::BGK<moment::useFieldrhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  // task collection
  using TaskCollection = tmp::TupleWrapper<BulkTask, WallTask>;
  // task executor
  using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;

  using Momenta = moment::MomentaTuple<moment::rhoU<CELL, true>, moment::shearRateMag<CELL, true>>;
  using RhoUTask = tmp::Key_TypePair<AABBFlag|PeriodicFlag, Momenta>;
  using TaskCollectionRhoU = tmp::TupleWrapper<RhoUTask>;
  using TaskSelectorRhoU = tmp::TaskSelector<TaskCollectionRhoU, std::uint8_t, CELL>;

  // writers
  vtmo::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());
  vtmo::PhysVectorWriter PhysVecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>(),
    std::bind(&BaseConverter<T>::getPhysU<2>, &BaseConv, std::placeholders::_1));
  vtmo::PhysScalarWriter PhysShearRateWriter("ShearRate", NSLattice.getField<StrainRateMag<T>>(),
    std::bind(&BaseConverter<T>::getPhysStrainRate, &BaseConv, std::placeholders::_1));
  vtmo::vtmWriter<T, LatSet::d> NSWriter("couette2d", Geo);
  NSWriter.addWriterSet(RhoWriter, PhysVecWriter, PhysShearRateWriter);
  
  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;

  NSWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && res > tol) {

    NSLattice.ApplyCellDynamics<NSTask>(FlagFM);
    NS_Peri.Apply();
    NSLattice.Stream();
    NSLattice.NormalFullCommunicate();
    BM.Apply();

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