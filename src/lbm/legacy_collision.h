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

// legacy_collision.h
//
// Old version of the BGK collision operator, only used by the legacy
// PopLattice path (data_struct/lattice.hh, lbm/buoyancy.h).
// Kept out of lbm/collision.h so active umbrella builds do not compile it.

#pragma once

#include "lbm/collision.h"  // PopCell, moment, latset


namespace legacy {

template <typename T, typename LatSet>
struct BGK {
  // BGK collision operator
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  __any__ static void apply(PopCell<T, LatSet>& cell) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }

  // BGK collision operator with force
  // template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  // __any__ static void applyForce(PopCell<T, LatSet>& cell, const Vector<T, LatSet::d>& force) {
  //   std::array<T, LatSet::q> feq{};
  //   GetFeq(feq, cell.getVelocity(), cell.getRho());

  //   std::array<T, LatSet::q> fi{};
  //   force::ForcePop<T, LatSet>::compute(fi, cell.getVelocity(), force);

  //   const T omega = cell.getOmega();
  //   const T _omega = cell.get_Omega();
  //   const T fomega = cell.getfOmega();

  //   for (int i = 0; i < LatSet::q; ++i) {
  //     cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
  //   }
  // }

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  __any__ static void applySource(PopCell<T, LatSet>& cell, const std::array<T, LatSet::q>& fi) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  __any__ static void applySource(PopCell<T, LatSet>& cell, const T S) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * S * latset::w<LatSet>(i);
    }
  }
};

}  // namespace legacy
