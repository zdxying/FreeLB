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

// moment.h

#pragma once

#include "data_struct/cell.h"

namespace moment {

template <typename CELL, bool Write_To_Field = false>
struct rho {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static T get(CELL& cell) {
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
    if constexpr (Write_To_Field) cell.getRho() = rho_value;
    return rho_value;
  }
  static void apply(CELL& cell) {
    T& rho_value = cell.getRho();
    rho_value = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
  }
  static void apply(CELL& cell, T& rho_value) {
    rho_value = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
    if constexpr (Write_To_Field) cell.getRho() = rho_value;
  }

  // apply with source term
};

template <typename CELL, bool Write_To_Field = false>
struct u {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static Vector<T, LatSet::d> get(CELL& cell) {
    Vector<T, LatSet::d> u_value;
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value /= rho_value;
    if constexpr (Write_To_Field) cell.getVelocity() = u_value;
    return u_value;
  }
  static void apply(CELL& cell) {
    Vector<T, LatSet::d>& u_value = cell.getVelocity();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value /= rho_value;
  }
  static void apply(CELL& cell, Vector<T, LatSet::d>& u_value) {
    u_value.clear();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value /= rho_value;
    if constexpr (Write_To_Field) cell.getVelocity() = u_value;
  }
};

template <typename CELL, bool Write_To_Field = false>
struct rhou {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell, T& rho_value, Vector<T, LatSet::d>& u_value) {
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value /= rho_value;
    if constexpr (Write_To_Field) {
      cell.getRho() = rho_value;
      cell.getVelocity() = u_value;
    }
  }
  static void apply(CELL& cell) {
    T& rho_value = cell.getRho();
    Vector<T, LatSet::d>& u_value = cell.getVelocity();
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value /= rho_value;
  }
};

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
struct RhoVelocity {
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

}  // namespace moment