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
// checkpoint / restart
int Restart;
int CheckpointStep;
std::string RestartDir;

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
  // checkpoint / restart: off by default, so old .ini files keep working
  Restart = param_reader.getValueOr<int>("Simulation_Settings", "Restart", 0);
  CheckpointStep =
    param_reader.getValueOr<int>("Simulation_Settings", "CheckpointStep", 0);
  RestartDir =
    param_reader.getValueOr<std::string>("Simulation_Settings", "RestartDir", "./ckpt");

  MPI_RANK(0)
  std::cout << "------------Simulation Parameters:-------------\n" << std::endl;
  std::cout << "[Simulation_Settings]:" << "TotalStep:         " << MaxStep << "\n"
            << "OutputStep:        " << OutputStep << "\n"
            << "Tolerance:         " << tol << "\n"
            << "Restart:           " << (Restart ? "on" : "off") << "\n"
            << "CheckpointStep:    " << CheckpointStep << "\n"
            << "RestartDir:        " << RestartDir << "\n"
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

  MPI_DEBUG_WAIT;

  Printer::Print_BigBanner(std::string("Initializing..."));

  readParam();

  // ------------------ geometry file ------------------
  // The grid/domain/voxel size come from the repository recipe; [Mesh] must
  // agree, otherwise the converter and the block geometry would describe a
  // different lattice than the mask.
  flbgeo::GeometryFile<T> geometryFile =
    flbgeo::readGeometryFile<T>("cavity3d.flbgeo");
  const flbgeo::VoxelGeometry<T>& record = geometryFile.Geometry;
  if (Ni != record.getNx() || Nj != record.getNy() || Nk != record.getNz() ||
      Cell_Len != record.getVoxelSize()) {
    IF_MPI_RANK(0) {
      std::cerr << "[cavity3d] geometry file " << record.getNx() << "x"
                << record.getNy() << "x" << record.getNz() << " (cell "
                << record.getVoxelSize() << ") does not match [Mesh] " << Ni
                << "x" << Nj << "x" << Nk << " (cell " << Cell_Len << ")"
                << std::endl;
    }
    return 1;
  }

  // converters
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, 3> cavity = record.getDomain();
  AABB<T, 3> toplid(Vector<T, 3>{Cell_Len, Cell_Len, T((Nk - 1) * Cell_Len)},
    Vector<T, 3>{T((Ni - 1) * Cell_Len), T((Nj - 1) * Cell_Len), T(Nk * Cell_Len)});

  // method 1: dirctly define geometry, [serial][openmp]
  // BlockGeometry3D<T> Geo(Ni, Nj, Nk, Thread_Num, cavity, Cell_Len);
  // end method 1

  // method 2: use geohelper for complex geometry, [mpi]
  BlockGeometryHelper3D<T> GeoHelper(Ni, Nj, Nk, cavity, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks();
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry3D<T> Geo(GeoHelper);
  // end method 2

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  // occupancy from the record: Domain cells are bulk, solid/outside are void
  flbgeo::writeVoxelFlags(record, FlagFM, AABBFlag, VoidFlag);
  // domain-edge walls have no surface band, so the cavity box still selects
  // them; solid surfaces are selected by their patches below
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
  FlagFM.forEach(toplid, [&](FLAG& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, BBMovingWallFlag);
  });
  // the obstacle's Wall patch drives its bounce-back cells
  flbbc::PatchFlagMap patchFlags;
  patchFlags.wall = BouncebackFlag;
  const flbbc::PatchBindStats patchStats =
    flbbc::bindPatches(record, FlagFM, patchFlags, AABBFlag);
  IF_MPI_RANK(0) {
    std::cout << "[cavity3d] obstacle patch: " << patchStats.wall
              << " cells, unmapped " << patchStats.unmappedCells << std::endl;
  }
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
  // a patch cell keeps its bulk bit and gains the wall bit
  using BBPatchTask = tmp::Key_TypePair<std::uint8_t(AABBFlag | BouncebackFlag),
                                        collision::BounceBack<CELL>>;
  // task collection
  using TaskCollection =
    tmp::TupleWrapper<BulkTask, BBTask, BBMVTask, BBPatchTask>;
  // task executor
  using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;

  // task: update rho and u
  using RhoUTask = tmp::Key_TypePair<AABBFlag, moment::rhoU<CELL, true>>;
  using RhoUPatchTask =
    tmp::Key_TypePair<std::uint8_t(AABBFlag | BouncebackFlag),
                      moment::rhoU<CELL, true>>;
  using TaskCollectionRhoU = tmp::TupleWrapper<RhoUTask, RhoUPatchTask>;
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

  // ------------------ restart: restore the snapshot ------------------
  // The snapshot carries the geometry identity (B7a) and is validated against
  // the record above; flags and BCs are replayed from the same file.
  // Must run before the main loop; the BC managers only hold references, so
  // nothing has to be rebuilt after a restore.
  Timer MainLoopTimer;
  if (Restart) {
    if (!flbckpt::checkCompatibleWithGeometry(RestartDir, -1, record, NSLattice,
                                              FlagFM)) {
      Printer::Print_BigBanner(std::string("Restart aborted: incompatible checkpoint"));
      return 1;
    }
    MainLoopTimer.set(
      flbckpt::readWithGeometry(RestartDir, -1, record, NSLattice, FlagFM));
    // Re-seed the tolerance history: UOld/RhoOld are not part of the snapshot
    // (they are overwritten after every evaluation and do not affect the physics)
    NSLattice.EnableToleranceU();
    IF_MPI_RANK(0) {
      std::cout << "[Restart] resumed at step " << MainLoopTimer() << std::endl;
    }
  }

  // count and timer
  Timer OutputTimer;

  NSWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep && res > tol) {
    // normal comm
    // NSLattice.ApplyCellDynamics<NSTask>(FlagFM);
    // NSLattice.Stream();
    // NSLattice.NormalCommunicate();
    //

    // hidden comm
    NSLattice.Stream();
    NSLattice.ApplyCellDynamicsWithHiddenComm<NSTask>(FlagFM);
    //

    // BM.Apply(MainLoopTimer());


    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      // normal comm
      //
      // hidden comm
      NSLattice.NormalFullCommunicate();
      //
      NSLattice.ApplyCellDynamics<TaskSelectorRhoU>(FlagFM);

      res = NSLattice.getToleranceU(-1);
      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Print_Res<T>(res);
      Printer::Endl();
      // NSWriter.WriteBinary(MainLoopTimer());

      // Snapshot point: communication is done, before the next Stream (quiescent)
      if (CheckpointStep > 0 && MainLoopTimer() % CheckpointStep == 0)
        flbckpt::writeWithGeometry(RestartDir, MainLoopTimer(), record, NSLattice,
                                   FlagFM);
    }
  }

  // Write one more at the end, so a run stopped early by the tolerance can resume
  if (CheckpointStep > 0) {
    NSLattice.NormalFullCommunicate();
    flbckpt::writeWithGeometry(RestartDir, MainLoopTimer(), record, NSLattice, FlagFM);
  }

#ifdef XCORE_ENABLED
  xcore::MemoryPool::getInstance().print_status();
#endif

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();

  return 0;
}