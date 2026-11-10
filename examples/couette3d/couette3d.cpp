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

// couette3d.cpp
// Plane Couette flow in a channel periodic along x and y.  Both walls are solid
// slabs; the bottom is static and the top moves along x.  The two Wall patches
// have the same kind, so the per-entity flag override in bindPatches is what
// gives them different boundary tasks: entity 2 -> bounce-back, entity 3 ->
// moving-wall bounce-back.  The .flbgeo supplies the geometry and the patches.

#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;
using LatSet = D3Q19<T>;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t BulkFlag = std::uint8_t(2);
constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
constexpr std::uint8_t MovingWallFlag = std::uint8_t(8);
constexpr std::uint8_t PeriodicFlag = std::uint8_t(16);

int Ni;
int Nj;
int Nk;
T Cell_Len;
int BlockCellLen;
int Thread_Num;
T RT;
T rho_ref;
T Kine_Visc;
T U_Max;
T U_Wall;
int MaxStep;
int OutputStep;
T tol;

void readParam() {
  iniReader param_reader("couette3d.ini");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");
  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Nk = param_reader.getValue<int>("Mesh", "Nk");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");
  RT = param_reader.getValue<T>("LB", "RT");
  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");
  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");
  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");
  U_Wall = param_reader.getValue<T>("Boundary_Conditions", "U_Wall");
  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");
  tol = param_reader.getValue<T>("tolerance", "tol");
}

}  // namespace

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT;
  readParam();

  flbgeo::GeometryFile<T> geometryFile =
    flbgeo::readGeometryFile<T>("couette3d.flbgeo");
  const flbgeo::VoxelGeometry<T>& record = geometryFile.Geometry;
  if (Ni != record.getNx() || Nj != record.getNy() || Nk != record.getNz() ||
      Cell_Len != record.getVoxelSize()) {
    IF_MPI_RANK(0) {
      std::cerr << "[couette3d] geometry file does not match [Mesh]" << std::endl;
    }
    return 1;
  }

  BaseConverter<T> conv(LatSet::cs2);
  conv.ConvertFromRT(Cell_Len, RT, rho_ref, Ni * Cell_Len, U_Max, Kine_Visc);

  const AABB<T, 3> domain = record.getDomain();
  BlockGeometryHelper3D<T> helper(Ni, Nj, Nk, domain, Cell_Len, BlockCellLen);
  helper.CreateBlocks();
  helper.LoadBalancing(mpi().getSize());
  BlockGeometry3D<T> geo(helper);

  BlockFieldManager<FLAG, T, 3> FlagFM(geo, VoidFlag);
  flbgeo::writeVoxelFlags(record, FlagFM, BulkFlag, VoidFlag);

  // periodic x/y faces: the solid slabs span the whole cross-section, so only
  // fluid cells are painted
  const AABB<T, 3> left(Vector<T, 3>{T(0), T(0), T(0)},
                        Vector<T, 3>{T(1), T(Nj), T(Nk)});
  const AABB<T, 3> right(Vector<T, 3>{T(Ni - 1), T(0), T(0)},
                         Vector<T, 3>{T(Ni), T(Nj), T(Nk)});
  const AABB<T, 3> front(Vector<T, 3>{T(0), T(0), T(0)},
                         Vector<T, 3>{T(Ni), T(1), T(Nk)});
  const AABB<T, 3> back(Vector<T, 3>{T(0), T(Nj - 1), T(0)},
                        Vector<T, 3>{T(Ni), T(Nj), T(Nk)});
  auto makePeriodic = [&](const AABB<T, 3>& box) {
    FlagFM.forEach(box, [&](FLAG& field, std::size_t id) {
      if (util::isFlag(field.get(id), BulkFlag)) {
        field.SetField(id, std::uint8_t(BulkFlag | PeriodicFlag));
      }
    });
  };
  makePeriodic(left);
  makePeriodic(right);
  makePeriodic(front);
  makePeriodic(back);

  // both patches are Wall; the entity override gives the top slab its own task
  flbbc::PatchFlagMap patchFlags;
  patchFlags.wall = BouncebackFlag;
  std::unordered_map<std::uint16_t, std::uint8_t> byEntity;
  byEntity[3] = MovingWallFlag;
  const flbbc::PatchBindStats patchStats = flbbc::bindPatches(
    record, FlagFM, patchFlags, BulkFlag, &byEntity, PeriodicFlag);
  IF_MPI_RANK(0) {
    std::cout << "[couette3d] patches: " << patchStats.wall
              << " cells (unmapped " << patchStats.unmappedCells << ")"
              << std::endl;
  }

  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  using CELL = Cell<T, LatSet, FIELDS>;
  ValuePack initValues(conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> lattice(geo, initValues, conv);

  // the moving-wall task reads the cell velocity
  Vector<T, LatSet::d> wallVelocity =
    conv.getLatticeU(Vector<T, LatSet::d>{U_Wall, T(0), T(0)});
  auto& velFM = lattice.getField<VELOCITY<T, LatSet::d>>();
  for (std::size_t ib = 0; ib < FlagFM.getBlockFields().size(); ++ib) {
    auto& flagBlock = FlagFM.getBlockField(static_cast<int>(ib)).getField(0);
    auto& velBlock = velFM.getBlockField(static_cast<int>(ib)).getField(0);
    for (std::size_t id = 0; id < flagBlock.size(); ++id) {
      if (util::isFlag(flagBlock[id], MovingWallFlag)) {
        velBlock.set(id, wallVelocity);
      }
    }
  }

  using BulkTask = tmp::Key_TypePair<BulkFlag,
    collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  // The boundary flags replaced the bulk bit, so the keys are the boundary
  // bits themselves; isFlag() is a bit test, so the most specific task goes
  // first and Bulk (bit 2) last.
  using PeriodicBulkTask = tmp::Key_TypePair<PeriodicFlag,
    collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  using BBTask = tmp::Key_TypePair<BouncebackFlag, collision::BounceBack<CELL>>;
  using MovingTask = tmp::Key_TypePair<MovingWallFlag,
    collision::BounceBackMovingWall<CELL>>;
  using TaskCollection =
    tmp::TupleWrapper<MovingTask, BBTask, PeriodicBulkTask, BulkTask>;
  using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;
  using RhoUTask = tmp::Key_TypePair<BulkFlag, moment::rhoU<CELL, true>>;
  using PeriodicRhoUTask = tmp::Key_TypePair<PeriodicFlag,
                                             moment::rhoU<CELL, true>>;
  using TaskCollectionRhoU =
    tmp::TupleWrapper<RhoUTask, PeriodicRhoUTask>;
  using TaskSelectorRhoU =
    tmp::TaskSelector<TaskCollectionRhoU, std::uint8_t, CELL>;

  FixedPeriodicBoundaryManager<BlockLatticeManager<T, LatSet, FIELDS>,
                               BlockFieldManager<FLAG, T, 3>>
    NS_Peri("NS_Peri", lattice, FlagFM, PeriodicFlag, VoidFlag);
  NS_Peri.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  NS_Peri.Setup(front, NbrDirection::YN, back, NbrDirection::YP);

  vtmo::PhysVectorWriter velocityWriter(
    "physVelocity", velFM,
    std::bind(&BaseConverter<T>::getPhysU<LatSet::d>, &conv,
              std::placeholders::_1));
  vtmo::vtmWriter<T, LatSet::d> writer("couette3d", geo);
  writer.addWriterSet(velocityWriter);

  for (int step = 0; step < MaxStep; ++step) {
    lattice.ApplyCellDynamics<NSTask>(FlagFM);
    NS_Peri.Apply();
    lattice.Stream();
    lattice.NormalFullCommunicate();
    if (OutputStep > 0 && (step + 1) % OutputStep == 0) {
      writer.WriteBinary(step + 1);
    }
  }

  // refresh rho/u from the populations before reading the profile
  lattice.ApplyCellDynamics<TaskSelectorRhoU>(FlagFM);

  // Centerline profile at x = y = mid.  The x/y periodic faces carry no wall
  // task (they are refreshed by the periodic copy), so a plane average includes
  // that face strip; the centreline is the meaningful profile.
  const int midX = Ni / 2;
  const int midY = Nj / 2;
  std::vector<double> ux(Nk, 0.0);
  std::vector<char> have(Nk, 0);
  for (std::size_t ib = 0; ib < geo.getBlocks().size(); ++ib) {
    const auto& block = geo.getBlock(static_cast<int>(ib));
    const int overlap = block.getOverlap();
    const auto& vel = velFM.getBlockField(static_cast<int>(ib)).getField(0);
    for (int k = 0; k < Nk; ++k) {
      if (have[k]) continue;
      const auto loc = Vector<int, 3>{midX + overlap, midY + overlap, k + overlap};
      const Vector<T, 3> pos = block.getVoxel(loc);
      if (std::abs(pos[0] - (T(midX) + T(0.5))) > T(1e-9) ||
          std::abs(pos[1] - (T(midY) + T(0.5))) > T(1e-9)) {
        continue;
      }
      ux[k] = static_cast<double>(vel[block.getIndex(loc)][0]);
      have[k] = 1;
    }
  }
  IF_MPI_RANK(0) {
    std::cout << "[couette3d] lattice wall u = " << wallVelocity[0]
              << "\n[couette3d] z  ux  couette (lattice units)" << std::endl;
    const double wallSpeed = static_cast<double>(wallVelocity[0]);
    for (int k = 1; k + 1 < Nk; ++k) {
      if (!have[k]) continue;
      const double zc = static_cast<double>(k) + 0.5;
      const double exact = wallSpeed * (zc - 1.0) / static_cast<double>(Nk - 2);
      std::cout << zc << ' ' << ux[k] << ' ' << exact << std::endl;
    }
  }

  return 0;
}
