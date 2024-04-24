/* This file is part of FreeLB
 * 
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 * 
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 * 
 */

// collision.h

#pragma once

#include "data_struct/cell.h"
#include "lbm/force.h"

namespace collision {

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

}  // namespace collision