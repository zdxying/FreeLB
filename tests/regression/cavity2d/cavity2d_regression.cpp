/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 */

// Numerical regression: a small, serial lid-driven cavity (D2Q9, BGK).
//
// The point is to pin the *solver* itself, not a container: the whole CPU LBM
// pipeline (BlockGeometry2D -> BlockLatticeManager -> Stream + ApplyCellDynamics
// with BGK / bounce-back / bounce-back-moving-wall tasks) runs for a fixed
// number of steps on a fixed grid and its distribution-function checksum is
// compared against a pinned baseline.  A change in the scheme, the equilibrium,
// the boundary treatment or the streaming order moves the checksum and fails
// here instead of silently changing every simulation.
//
// The grid is tiny and single-block so the test runs in well under a second and
// is independent of MPI / OpenMP scheduling.

#include "freelb.h"
#include "freelb.hh"
#include "utils/field_checksum.h"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

// Pinned baseline (FLOAT_TYPE=double).  Relative tolerance below absorbs the
// last-bits differences between compilers / optimization levels; a real change
// to the solver moves sumsq far more than that.
static constexpr double kBaselineSumSq = 7728.702031;  // 32x32 lid cavity, 1000 steps
static constexpr double kRelTol = 1e-6;

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t AABBFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
  constexpr std::uint8_t BBMovingWallFlag = std::uint8_t(8);

  mpi().init(&argc, &argv);

  const int Ni = 32;
  const int Nj = 32;
  const T Cell_Len = T(1);
  const int NSTEPS = 1000;

  // converters: rho = 1, characteristic u = 0.1, RT = 0.56
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni, T(0.1), T(0.56));

  // geometry (single block)
  AABB<T, 2> cavity(Vector<T, 2>{}, Vector<T, 2>{T(Ni * Cell_Len), T(Nj * Cell_Len)});
  AABB<T, 2> toplid(Vector<T, 2>{Cell_Len, T((Nj - 1) * Cell_Len)},
                    Vector<T, 2>{T((Ni - 1) * Cell_Len), T(Nj * Cell_Len)});
  BlockGeometry2D<T> Geo(Ni, Nj, 1, cavity, Cell_Len);

  // flag field
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(cavity, [&](FLAG& field, std::size_t id) {
    field.SetField(id, AABBFlag);
  });
  FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
  FlagFM.forEach(toplid, [&](FLAG& field, std::size_t id) {
    if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, BBMovingWallFlag);
  });

  // lattice
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  using CELL = Cell<T, LatSet, FIELDS>;
  ValuePack InitValues(BaseConv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(Geo, InitValues, BaseConv);

  // moving-lid velocity
  Vector<T, LatSet::d> U_Wall(T{0.1}, T{0});
  Vector<T, LatSet::d> LatU_Wall = BaseConv.getLatticeU(U_Wall);
  NSLattice.getField<VELOCITY<T, LatSet::d>>().forEach(
      toplid, FlagFM, BBMovingWallFlag,
      [&](auto& field, std::size_t id) { field.SetField(id, LatU_Wall); });

  // tasks: bulk BGK + bounce-back + bounce-back moving wall
  using BulkTask = tmp::Key_TypePair<
      AABBFlag, collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  using BBTask = tmp::Key_TypePair<BouncebackFlag, collision::BounceBack<CELL>>;
  using BBMVTask =
      tmp::Key_TypePair<BBMovingWallFlag, collision::BounceBackMovingWall<CELL>>;
  using TaskCollection = tmp::TupleWrapper<BulkTask, BBTask, BBMVTask>;
  using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;
  using RhoUTask = tmp::Key_TypePair<AABBFlag, moment::rhoU<CELL, true>>;
  using TaskSelectorRhoU =
      tmp::TaskSelector<tmp::TupleWrapper<RhoUTask>, std::uint8_t, CELL>;

  for (int step = 0; step < NSTEPS; ++step) {
    NSLattice.Stream();
    NSLattice.ApplyCellDynamics<NSTask>(FlagFM);
  }
  NSLattice.ApplyCellDynamics<TaskSelectorRhoU>(FlagFM);

  // distribution-function checksum
  auto& popf = NSLattice.getBlockLat(0).getField<POP<T, LatSet::q>>();
  const std::size_t n = NSLattice.getBlockLat(0).getN();
  const std::array<double, 3> s = frelb_diag::PopChecksumHost(popf, n, LatSet::q);
  frelb_diag::PrintChecksum("cavity2d", s);

  // no blow-up / NaN
  CHECK(std::isfinite(s[0]) && std::isfinite(s[1]) && std::isfinite(s[2]));
  // mass conservation: sum of all pops == N_cells * rho_init (rho_init = 1)
  CHECK_NEAR(s[0], static_cast<double>(Geo.getTotalCellNum()), 1e-9);
  // pinned fingerprint of the field's fine structure
  CHECK_NEAR(s[1], kBaselineSumSq, kBaselineSumSq * kRelTol);

  return flbtest::summary("regression/cavity2d");
}
