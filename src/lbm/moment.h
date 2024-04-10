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

// moment.h

#pragma once

#include "data_struct/cell.h"

namespace moment {

template <typename T, typename LatSet>
struct Rho {
  // return rho
  static T get(const BasicCell<T, LatSet>& cell) {
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho += cell[i];
    return rho;
  }

  // return rho with source
  static T get(const BasicCell<T, LatSet>& cell, const T source) {
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho += cell[i];
    rho += source * T(0.5);
    return rho;
  }

  // compute rho
  static void apply(const BasicCell<T, LatSet>& cell, T& rho) {
    rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho += cell[i];
  }

  // compute rho with source: C = sum(f) + q/2
  static void apply(const BasicCell<T, LatSet>& cell, T& rho, const T source) {
    rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho += cell[i];
    rho += source * T(0.5);
  }
};


template <typename T, typename LatSet>
struct Velocity {
  // get velocity
  static Vector<T, LatSet::d> get(const BasicCell<T, LatSet>& cell) {
    Vector<T, LatSet::d> u;
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + LatSet::c[i] * cell[i];
    }
    return u / rho;
  }

  // get velocity with force
  static Vector<T, LatSet::d> get(const BasicCell<T, LatSet>& cell,
                                  const std::array<T, LatSet::q>& fi) {
    Vector<T, LatSet::d> u;
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + (LatSet::c[i] * (cell[i] + T(0.5) * fi[i]));
    }
    return u / rho;
  }

  // compute velocity
  static void apply(const BasicCell<T, LatSet>& cell, Vector<T, LatSet::d>& u) {
    u.clear();
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + LatSet::c[i] * cell[i];
    }
    u = u / rho;
  }

  // compute velocity with force
  static void apply(const BasicCell<T, LatSet>& cell, Vector<T, LatSet::d>& u,
                    const std::array<T, LatSet::q>& fi) {
    u.clear();
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + (LatSet::c[i] * (cell[i] + T(0.5) * fi[i]));
    }
    u = u / rho;
  }
};

template <typename T, typename LatSet>
struct RhoVelocity{
  // compute rho and velocity
  static void apply(const BasicCell<T, LatSet>& cell, T& rho, Vector<T, LatSet::d>& u) {
    rho = T(0);
    u.clear();
    for (int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + LatSet::c[i] * cell[i];
    }
    u = u / rho;
  }
};

} // namespace moment