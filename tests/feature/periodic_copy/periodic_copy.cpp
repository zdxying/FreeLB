/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress on FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 *
 */

// periodic_copy.cpp
// The periodic exchange must place every population that leaves a face on the
// opposite halo cell.  Stream() pulls each direction from (cell - c), so a copy
// that writes the face cell itself lands one cell inside the domain instead.
// A fully periodic box has no wall and therefore no other mass or momentum
// source, which makes the exchange testable: after one Apply() every halo cell
// must equal its periodic image, and collide + copy + stream must conserve
// mass and momentum to roundoff.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "freelb.h"
#include "freelb.hh"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

int Thread_Num = 1;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t BulkFlag = std::uint8_t(2);
constexpr std::uint8_t PeriodicFlag = std::uint8_t(4);
constexpr int Nx = 8;
constexpr int Ny = 8;
constexpr int Steps = 200;
constexpr T Amplitude = T(0.02);
// a uniform mean flow gives the momentum check a non-zero value to defend
constexpr T MeanU = T(0.01);

using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
using CELL = Cell<T, LatSet, FIELDS>;
using BulkTask = tmp::Key_TypePair<
  BulkFlag, collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
using PeriodicBulkTask = tmp::Key_TypePair<
  std::uint8_t(BulkFlag | PeriodicFlag),
  collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
using TaskCollection = tmp::TupleWrapper<BulkTask, PeriodicBulkTask>;
using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;

}  // namespace

int main() {
  const AABB<T, 2> domain(Vector<T, 2>{}, Vector<T, 2>{T(Nx), T(Ny)});
  BlockGeometry2D<T> geo(Nx, Ny, 1, domain, T(1));
  BlockFieldManager<FLAG, T, 2> flagFM(geo, VoidFlag);
  flagFM.forEach(domain, [&](FLAG& field, std::size_t id) {
    field.SetField(id, BulkFlag);
  });

  const AABB<T, 2> left(Vector<T, 2>{}, Vector<T, 2>{T(1), T(Ny)});
  const AABB<T, 2> right(Vector<T, 2>{T(Nx - 1), T(0)},
                         Vector<T, 2>{T(Nx), T(Ny)});
  const AABB<T, 2> bottom(Vector<T, 2>{}, Vector<T, 2>{T(Nx), T(1)});
  const AABB<T, 2> top(Vector<T, 2>{T(0), T(Ny - 1)},
                       Vector<T, 2>{T(Nx), T(Ny)});
  const AABB<T, 2> cornerLL(Vector<T, 2>{}, Vector<T, 2>{T(1), T(1)});
  const AABB<T, 2> cornerLR(Vector<T, 2>{T(Nx - 1), T(0)},
                            Vector<T, 2>{T(Nx), T(1)});
  const AABB<T, 2> cornerUL(Vector<T, 2>{T(0), T(Ny - 1)},
                            Vector<T, 2>{T(1), T(Ny)});
  const AABB<T, 2> cornerUR(Vector<T, 2>{T(Nx - 1), T(Ny - 1)},
                            Vector<T, 2>{T(Nx), T(Ny)});
  auto makePeriodic = [&](const AABB<T, 2>& box) {
    flagFM.forEach(box, [&](FLAG& field, std::size_t id) {
      field.SetField(id, std::uint8_t(BulkFlag | PeriodicFlag));
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
                               BlockFieldManager<FLAG, T, 2>>
    NS_Peri("NS_Peri", NSLattice, flagFM, PeriodicFlag, VoidFlag);
  NS_Peri.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  NS_Peri.Setup(bottom, NbrDirection::YN, top, NbrDirection::YP);
  // A corner column is the periodic image of the diagonal corner, so the two
  // face pairs cannot supply its diagonal populations.
  NS_Peri.Setup(cornerLL, NbrDirection::XN | NbrDirection::YN, cornerUR,
                NbrDirection::XP | NbrDirection::YP);
  NS_Peri.Setup(cornerLR, NbrDirection::XP | NbrDirection::YN, cornerUL,
                NbrDirection::XN | NbrDirection::YP);

  const auto& block = geo.getBlock(0);
  const int overlap = block.getOverlap();
  auto cellId = [&](int x, int y) {
    return block.getIndex(Vector<int, 2>{x + overlap, y + overlap});
  };
  auto& blat = NSLattice.getBlockLat(0);

  // A divergence-free mean flow plus a density wave, written as populations:
  // neither is a periodic image of the other, so a copy that reads the wrong
  // cell or the wrong direction shows up immediately.
  constexpr T TwoPi = T(6.283185307179586);
  std::array<T, LatSet::q> feq;
  for (int y = 0; y < Ny; ++y) {
    for (int x = 0; x < Nx; ++x) {
      const T rho = T(1) + Amplitude * std::sin(TwoPi * T(x) / T(Nx)) *
                               std::sin(TwoPi * T(y) / T(Ny));
      const Vector<T, 2> u{MeanU + Amplitude * std::sin(TwoPi * T(y) / T(Ny)),
                           Amplitude * std::cos(TwoPi * T(x) / T(Nx))};
      equilibrium::SecondOrder<CELL>::apply(feq, rho, u);
      CELL cell(cellId(x, y), blat);
      for (unsigned int k = 0; k < LatSet::q; ++k) cell[k] = feq[k];
    }
  }

  // ---- the copy itself: every halo cell must hold its periodic image ----
  NS_Peri.Apply();
  auto checkHalo = [&](int hx, int hy, int sx, int sy, NbrDirection dir) -> int {
    std::vector<unsigned int> dirs;
    getPeriodicPopDir<LatSet>(dir, dirs);
    CELL halo(cellId(hx, hy), blat);
    CELL image(cellId(sx, sy), blat);
    int checked = 0;
    for (unsigned int k : dirs) {
      CHECK_NEAR(halo[k], image[k], T(0));
      ++checked;
    }
    return checked;
  };
  int haloChecks = 0;
  for (int y = 0; y < Ny; ++y) {
    haloChecks += checkHalo(-1, y, Nx - 1, y, NbrDirection::XN);
    haloChecks += checkHalo(Nx, y, 0, y, NbrDirection::XP);
  }
  for (int x = 0; x < Nx; ++x) {
    haloChecks += checkHalo(x, -1, x, Ny - 1, NbrDirection::YN);
    haloChecks += checkHalo(x, Ny, x, 0, NbrDirection::YP);
  }
  haloChecks += checkHalo(-1, -1, Nx - 1, Ny - 1,
                          NbrDirection::XN | NbrDirection::YN);
  haloChecks += checkHalo(Nx, -1, 0, Ny - 1,
                          NbrDirection::XP | NbrDirection::YN);
  haloChecks += checkHalo(-1, Ny, Nx - 1, 0,
                          NbrDirection::XN | NbrDirection::YP);
  haloChecks += checkHalo(Nx, Ny, 0, 0,
                          NbrDirection::XP | NbrDirection::YP);

  // ---- collide + copy + stream conserves mass and momentum ----
  auto totals = [&](T& mass, T& px, T& py) {
    mass = px = py = T(0);
    for (int y = 0; y < Ny; ++y) {
      for (int x = 0; x < Nx; ++x) {
        CELL cell(cellId(x, y), blat);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          const T f = cell[k];
          mass += f;
          px += f * T(latset::c<LatSet>(k)[0]);
          py += f * T(latset::c<LatSet>(k)[1]);
        }
      }
    }
  };
  T mass0{}, px0{}, py0{};
  totals(mass0, px0, py0);
  for (int step = 0; step < Steps; ++step) {
    NSLattice.ApplyCellDynamics<NSTask>(flagFM);
    NS_Peri.Apply();
    NSLattice.Stream();
  }
  T mass1{}, px1{}, py1{};
  totals(mass1, px1, py1);

  std::cout << "[feature/periodic_copy] halo_checks=" << haloChecks
            << " mass=" << mass0 << "->" << mass1
            << " px=" << px0 << "->" << px1 << " py=" << py0 << "->" << py1
            << std::endl;

  CHECK_NEAR(mass1, mass0, T(1e-11) * mass0);
  CHECK_NEAR(px1, px0, T(1e-11) * std::max(std::abs(px0), T(1)));
  CHECK_NEAR(py1, py0, T(1e-11) * std::max(std::abs(py0), T(1)));
  // the flow must survive: a copy that drops the image decays it.
  CHECK(std::abs(px1) > T(0.5) * std::abs(px0));

  return flbtest::summary("feature/periodic_copy");
}
