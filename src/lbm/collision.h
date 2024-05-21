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

// collision.h

#pragma once

#include "data_struct/cell.h"
#include "lbm/force.h"

namespace collision {

// a typical BGK collision process with:
// macroscopic variables updated
// equilibrium distribution function calculated
template <typename EquilibriumScheme, bool WriteToField = false>
struct BGK_Feq_RhoU {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;

  static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    moment::template rhou<CELL, WriteToField>::apply(cell, rho, u);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(cell, feq, rho, u);
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }
  static void apply(CELL& cell, const Vector<T, LatSet::d>& force) {
    
  }
};

// a typical BGK collision process with:
// equilibrium distribution function calculated
template <typename EquilibriumScheme, bool WriteToField = false>
struct BGK_Feq {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;

  static void apply(CELL& cell) {
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(cell, feq, cell.getRho(), cell.getVelocity());
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }
};

// a typical BGK collision process
template <typename CELL, bool WriteToField = false>
struct bgk {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell, const std::array<T, LatSet::q>& feq) {
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }
};

template <typename T, typename LatSet>
struct BGK {
  // BGK collision operator
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  static void apply(Cell<T, LatSet>& cell) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();

    for (int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  static void apply(BCell<T, LatSet>& cell) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();

    for (int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }

  struct Source {
    static void apply(BCell<T, LatSet>& cell, const std::array<T, LatSet::q>& fi) {
      const T fomega = cell.getfOmega();
      for (int i = 0; i < LatSet::q; ++i) {
        cell[i] += fomega * fi[i];
      }
    }
  };


  // BGK collision operator with force
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  static void applyForce(Cell<T, LatSet>& cell, const Vector<T, LatSet::d>& force) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    std::array<T, LatSet::q> fi{};
    force::Force<T, LatSet>::ComputeForcePop(fi, cell.getVelocity(), force);

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  static void applyForce(BCell<T, LatSet>& cell, const Vector<T, LatSet::d>& force) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    std::array<T, LatSet::q> fi{};
    force::Force<T, LatSet>::ComputeForcePop(fi, cell.getVelocity(), force);

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  static void applySource(Cell<T, LatSet>& cell, const std::array<T, LatSet::q>& fi) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  static void applySource(BCell<T, LatSet>& cell, const std::array<T, LatSet::q>& fi) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  static void applySource(Cell<T, LatSet>& cell, const T S) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * S * LatSet::w[i];
    }
  }
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  static void applySource(BCell<T, LatSet>& cell, const T S) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * S * LatSet::w[i];
    }
  }
};

// full way bounce back, could be regarded as a mpdified collision process
// swap the populations in the opposite direction
// LatSet must have rest population(D2Q4 is not supported)
template <typename T, typename LatSet>
struct BounceBack {
  static constexpr int halfq = LatSet::q / 2 + 1;

  static void apply(BCell<T, LatSet>& cell) {
    for (int i = 1; i < halfq; ++i) {
      T temp = cell[i];
      cell[i] = cell[LatSet::opp[i]];
      cell[LatSet::opp[i]] = temp;
    }
  };
};

}  // namespace collision