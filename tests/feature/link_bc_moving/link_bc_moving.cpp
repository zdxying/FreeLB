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
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// link_bc_moving.cpp
// B4 moving-wall baseline: Couette flow in the same periodic channel as
// link_bc_accuracy.  The bottom wall (entity 2) is static, the top wall
// (entity 3) moves along x with Uw; the analytic profile is linear in z.
// The 6000-step comparison is a short-horizon baseline: the moving-wall
// scheme still drifts at very long times, so this is not a converged value.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <unordered_map>

#include "freelb.h"
#include "freelb.hh"
#include "boundary/link_bc.h"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

int Thread_Num = 1;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t BulkFlag = std::uint8_t(2);
constexpr std::uint8_t PeriodicFlag = std::uint8_t(4);
constexpr int Nx = 16;
constexpr int Ny = 16;
constexpr int Nz = 8;
constexpr T ZDomain0 = T(0.25);
constexpr T ZDomain1 = T(8.25);
constexpr T Wall0 = T(0.5);
constexpr T Wall1 = T(8.5);
constexpr T Uw = T(0.01);
constexpr int Steps = 6000;

using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
using CELL = Cell<T, LatSet, FIELDS>;
using BulkTask =
  tmp::Key_TypePair<BulkFlag,
                    collision::BGK<moment::rhoU<CELL>,
                                   equilibrium::SecondOrder<CELL>>>;
using PeriodicBulkTask =
  tmp::Key_TypePair<BulkFlag | PeriodicFlag,
                    collision::BGK<moment::rhoU<CELL>,
                                   equilibrium::SecondOrder<CELL>>>;
using TaskCollection = tmp::TupleWrapper<BulkTask, PeriodicBulkTask>;
using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;
using RhoUTask = tmp::Key_TypePair<BulkFlag, moment::rhoU<CELL, true>>;
using PeriodicRhoUTask =
  tmp::Key_TypePair<BulkFlag | PeriodicFlag, moment::rhoU<CELL, true>>;
using RhoUTaskCollection = tmp::TupleWrapper<RhoUTask, PeriodicRhoUTask>;
using TaskSelectorRhoU =
  tmp::TaskSelector<RhoUTaskCollection, std::uint8_t, CELL>;

}  // namespace

int main() {
  const AABB<T, 3> domain(Vector<T, 3>{T(0), T(0), ZDomain0},
                          Vector<T, 3>{T(Nx), T(Ny), ZDomain1});
  flbgeo::VoxelGeometryBuilder<T> builder(Nx, Ny, Nz, domain, T(1), T(1.5));
  builder.addBox(domain, 1, 1, flbgeo::VoxelRole::Domain);
  StlReader<T> low("../link_bc/stl/slab_low.stl", T(1), T(1), 2);
  StlReader<T> high("../link_bc/stl/slab_high.stl", T(1), T(1), 2);
  builder.addStl(low, 2, 2, flbgeo::VoxelRole::Solid, true);
  builder.addStl(high, 3, 3, flbgeo::VoxelRole::Solid, true);
  const flbgeo::VoxelGeometry<T> record = builder.build();
  CHECK(record.getBandSize() > 0);

  BlockGeometry3D<T> geo(Nx, Ny, Nz, 1, record.getDomain(),
                         record.getVoxelSize(), 1);
  BlockFieldManager<FLAG, T, 3> flagFM(geo, VoidFlag);
  flbgeo::writeVoxelFlags(record, flagFM, BulkFlag, VoidFlag);

  const AABB<T, 3> left(Vector<T, 3>{T(0), T(0), ZDomain0},
                        Vector<T, 3>{T(1), T(Ny), ZDomain1});
  const AABB<T, 3> right(Vector<T, 3>{T(Nx - 1), T(0), ZDomain0},
                         Vector<T, 3>{T(Nx), T(Ny), ZDomain1});
  const AABB<T, 3> bottom(Vector<T, 3>{T(0), T(0), ZDomain0},
                          Vector<T, 3>{T(Nx), T(1), ZDomain1});
  const AABB<T, 3> top(Vector<T, 3>{T(0), T(Ny - 1), ZDomain0},
                       Vector<T, 3>{T(Nx), T(Ny), ZDomain1});
  auto makePeriodic = [&](const AABB<T, 3>& box) {
    flagFM.forEach(box, [&](FLAG& field, std::size_t id) {
      if (util::isFlag(field.get(id), BulkFlag)) {
        field.SetField(id, std::uint8_t(BulkFlag | PeriodicFlag));
      }
    });
  };
  makePeriodic(left);
  makePeriodic(right);
  makePeriodic(bottom);
  makePeriodic(top);

  BaseConverter<T> conv(LatSet::cs2);
  conv.SimplifiedConverterFromRT(Nx, T(0.1), T(0.56));
  ValuePack initValues(conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(geo, initValues, conv);

  FixedPeriodicBoundaryManager<BlockLatticeManager<T, LatSet, FIELDS>,
                               BlockFieldManager<FLAG, T, 3>>
    NS_Peri("NS_Peri", NSLattice, flagFM, PeriodicFlag, VoidFlag);
  NS_Peri.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  NS_Peri.Setup(bottom, NbrDirection::YN, top, NbrDirection::YP);
  // The x and y faces share the corner column, whose diagonal populations
  // come from the diagonal partner; a face pair alone cannot supply them.
  const AABB<T, 3> cornerLL(Vector<T, 3>{T(0), T(0), ZDomain0},
                            Vector<T, 3>{T(1), T(1), ZDomain1});
  const AABB<T, 3> cornerLR(Vector<T, 3>{T(Nx - 1), T(0), ZDomain0},
                            Vector<T, 3>{T(Nx), T(1), ZDomain1});
  const AABB<T, 3> cornerUL(Vector<T, 3>{T(0), T(Ny - 1), ZDomain0},
                            Vector<T, 3>{T(1), T(Ny), ZDomain1});
  const AABB<T, 3> cornerUR(Vector<T, 3>{T(Nx - 1), T(Ny - 1), ZDomain0},
                            Vector<T, 3>{T(Nx), T(Ny), ZDomain1});
  NS_Peri.Setup(cornerLL, NbrDirection::XN | NbrDirection::YN, cornerUR,
                NbrDirection::XP | NbrDirection::YP);
  NS_Peri.Setup(cornerLR, NbrDirection::XP | NbrDirection::YN, cornerUL,
                NbrDirection::XN | NbrDirection::YP);

  std::unordered_map<std::uint16_t, Vector<T, 3>> wallVelocities;
  wallVelocities.emplace(std::uint16_t{3}, Vector<T, 3>{Uw, T(0), T(0)});
  flbbc::InterpolatedBBManager<decltype(NSLattice), decltype(flagFM)> ibb(
    NSLattice, flagFM, record, VoidFlag, PeriodicFlag, &wallVelocities, nullptr,
    std::uint8_t{3});
  CHECK(ibb.getLinkCount() > 0);

  const auto& block = geo.getBlock(0);
  const int overlap = block.getOverlap();
  auto cellId = [&](int x, int y, int z) {
    return block.getIndex(
      Vector<int, 3>{x + overlap, y + overlap, z + overlap});
  };
  auto findWallLink = [&](std::size_t id, int cz) -> const flbbc::SurfaceLink<T>* {
    for (const flbbc::SurfaceLink<T>& link : ibb.getLinks(0)) {
      if (link.CellId != id) continue;
      const auto c = latset::c<LatSet>(link.Direction);
      if (c[2] == cz) return &link;
    }
    return nullptr;
  };
  const flbbc::SurfaceLink<T>* bottomLink = findWallLink(cellId(8, 8, 0), -1);
  const flbbc::SurfaceLink<T>* topLink =
    findWallLink(cellId(8, 8, Nz - 1), 1);
  CHECK(bottomLink != nullptr);
  CHECK(topLink != nullptr);
  CHECK_NEAR(bottomLink->Fraction, T(0.25), T(1e-9));
  CHECK_NEAR(topLink->Fraction, T(0.75), T(1e-9));
  CHECK_EQ(bottomLink->Entity, std::uint16_t{2});
  CHECK_EQ(topLink->Entity, std::uint16_t{3});
  CHECK_NEAR(bottomLink->WallVelocity[0], T(0), T(1e-12));
  CHECK_NEAR(topLink->WallVelocity[0], Uw, T(1e-12));

  for (int step = 0; step < Steps; ++step) {
    NSLattice.ApplyCellDynamics<NSTask>(flagFM);
    NS_Peri.Apply();
    NSLattice.Stream();
    ibb.Apply();
  }
  NSLattice.ApplyCellDynamics<TaskSelectorRhoU>(flagFM);

  const auto& vel =
    NSLattice.getBlockLat(0).template getField<VELOCITY<T, LatSet::d>>();
  T err{};
  // least-squares slope over the interior layers, relative to the analytic
  // Couette slope: separates a uniform deficit from a near-wall kink
  T n{}, sz{}, szz{}, su{}, szu{};
  for (int z = 0; z < Nz; ++z) {
    const T zc = ZDomain0 + T(z) + T(0.5);
    const T expected = Uw * (zc - Wall0) / (Wall1 - Wall0);
    const T d = (vel.get(cellId(8, 8, z))[0] - expected) / Uw;
    err += d * d;
    const T u = vel.get(cellId(8, 8, z))[0];
    n += T(1);
    sz += zc;
    szz += zc * zc;
    su += u;
    szu += zc * u;
  }
  err = std::sqrt(err / T(Nz));
  const T slope = (n * szu - sz * su) / (n * szz - sz * sz);
  const T slopeRatio = slope / (Uw / (Wall1 - Wall0));

  std::cout << "[feature/link_bc_moving] links=" << ibb.getLinkCount()
            << " bottom_q=" << bottomLink->Fraction
            << " top_q=" << topLink->Fraction << " profile_l2=" << err
            << " slope_ratio=" << slopeRatio << std::endl;


  CHECK(err < T(0.10));

  return flbtest::summary("feature/link_bc_moving");
}
