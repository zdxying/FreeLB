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
std::string work_dir;

void readParam() {
  iniReader param_reader("openbd2dparam.ini");
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
  std::cout << "[Simulation_Settings]:" << "TotalStep:         " << MaxStep << "\n"
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
  constexpr std::uint8_t InletFlag = std::uint8_t(8);
  constexpr std::uint8_t OutletFlag = std::uint8_t(16);

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, LatSet::d> cavity(Vector<T, LatSet::d>{},
                            Vector<T, LatSet::d>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, LatSet::d> left(Vector<T, LatSet::d>{},
                          Vector<T, LatSet::d>(T(1), T(Nj - 1) * Cell_Len));
  AABB<T, LatSet::d> right(Vector<T, LatSet::d>(T(Ni - 1) * Cell_Len, T(1)),
                           Vector<T, LatSet::d>(T(Ni * Cell_Len), T(Nj - 1) * Cell_Len));
  BlockGeometry2D<T> Geo(Ni, Nj, Thread_Num, cavity, Cell_Len);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(cavity,
                 [&](FLAG& field, std::size_t id) { field.SetField(id, AABBFlag); });
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
  FlagFM.forEach(left, [&](FLAG& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, InletFlag);
  });
  FlagFM.forEach(right, [&](FLAG& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, OutletFlag);
  });

  vtmwriter::ScalarWriter FlagWriter("flag", FlagFM);
  vtmwriter::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag", Geo);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define lattice ------------------
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  using CELL = Cell<T, LatSet, FIELDS>;
  ValuePack InitValues(BaseConv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  // lattice
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(Geo, InitValues, BaseConv);
  NSLattice.EnableToleranceU();
  T res = 1;
  // set initial value of field
  Vector<T, LatSet::d> LatU_Wall = BaseConv.getLatticeU(U_Wall);
  NSLattice.getField<VELOCITY<T, LatSet::d>>().forEach(
    FlagFM, InletFlag,
    [&](auto& field, std::size_t id) { field.SetField(id, LatU_Wall); });

  // bcs
  BBLikeFixedBlockBdManager<bounceback::normal<CELL>,
                            BlockLatticeManager<T, LatSet, FIELDS>,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::movingwall<CELL>,
                            BlockLatticeManager<T, LatSet, FIELDS>,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    NS_Inlet("NS_BBMW", NSLattice, FlagFM, InletFlag, VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::anti_pressure<CELL>,
                            BlockLatticeManager<T, LatSet, FIELDS>,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    NS_Outlet("NS_BBMW", NSLattice, FlagFM, OutletFlag, VoidFlag);
  BlockBoundaryManager BM(&NS_BB, &NS_Inlet, &NS_Outlet);

  // define task/ dynamics:
  // rho u
  using RhoUTask = tmp::Key_TypePair<AABBFlag|OutletFlag, moment::rhou<CELL, true>>;
  using rhoTask = tmp::Key_TypePair<InletFlag, moment::rho<CELL>>;
  using RhoUTaskSelector = TaskSelector<std::uint8_t, CELL, RhoUTask>;
  // collision
  using collisionTask = collision::BGK<moment::UseFieldRhoU<CELL>, equilibrium::SecondOrder<CELL>>;

  // writers
  vtmwriter::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());
  vtmwriter::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmwriter::vtmWriter<T, LatSet::d> NSWriter("cavblock2d", Geo);
  NSWriter.addWriterSet(RhoWriter, VecWriter);

  // count and timer
  Timer MainLoopTimer;
  Timer OutputTimer;
  NSWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && res > tol) {

    NSLattice.ApplyCellDynamics<RhoUTaskSelector>(MainLoopTimer(), FlagFM);
    NSLattice.ApplyCellDynamics<collisionTask>(MainLoopTimer());

    NSLattice.Stream(MainLoopTimer());

    BM.Apply(MainLoopTimer());

    NSLattice.Communicate(MainLoopTimer());
    // use CommunicateAll() (less efficient) instead, 
    // because bounceback::normal used only inner cells to apply Bounceback
    // however, Bounceback is not applied to the outer cells which are to be communicated 
    // and partial(efficient) communication only communicate pops in certain directions
    // the bounce-back direction might not be communicated, becomming invalid pop
    // NSLattice.getField<POP<T, LatSet::q>>().CommunicateAll(MainLoopTimer());

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
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