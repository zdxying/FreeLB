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

// equlibrium.h

#pragma once
// lattice boltzmann method implementations

#include "data_struct/Vector.h"

template <typename T, typename LatSet>
class BasicCell;

template <typename T, typename LatSet>
class Cell;

// calc equilibrium distribution function
// sum(feq_i) = rho, for both first and second order
template <typename T, typename LatSet>
struct Equilibrium {
  static inline T Order1(int k, const Vector<T, LatSet::d> &u, T rho) {
    return LatSet::w[k] * rho * (T(1) + LatSet::InvCs2 * (u * LatSet::c[k]));
  }

  static inline T Order1_Incompresible(int k, const Vector<T, LatSet::d> &u, T rho) {
    return LatSet::w[k] * (rho + LatSet::InvCs2 * (u * LatSet::c[k]));
  }

  static inline T Order2(int k, const Vector<T, LatSet::d> &u, T rho, T u2) {
    T uc = u * LatSet::c[k];
    return LatSet::w[k] * rho *
           (T(1) + LatSet::InvCs2 * uc + uc * uc * T(0.5) * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T(0.5));
  }

  // static inline T Order2_Incompresible(int k, const Vector<T, LatSet::d> &u, T rho, T
  // u2) {
  //   T uc = u * LatSet::c[k];
  //   return LatSet::w[k] * (rho + LatSet::InvCs2 * uc + uc * uc * T(0.5) *
  //   LatSet::InvCs4 -
  //                          LatSet::InvCs2 * u2 * T(0.5));
  // }

  static void Feq_firstOrder(T *feq, const Vector<T, LatSet::d> &u, T rho) {
    for (int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order1(k, u, rho);
    }
  }
  static void FirstOrder(std::array<T, LatSet::q> &feq, const Vector<T, LatSet::d> &u,
                         T rho) {
    for (int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order1(k, u, rho);
    }
  }

  static void Feq_firstOrder_Incompresible(T *feq, const Vector<T, LatSet::d> &u, T rho) {
    for (int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order1_Incompresible(k, u, rho);
    }
  }
  static void FirstOrder_Incompresible(std::array<T, LatSet::q> &feq,
                                       const Vector<T, LatSet::d> &u, T rho) {
    for (int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order1_Incompresible(k, u, rho);
    }
  }

  static void Feq_secondOrder(T *feq, const Vector<T, LatSet::d> &u, T rho) {
    T u2 = u.getnorm2();
    for (int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order2(k, u, rho, u2);
    }
  }
  static void SecondOrder(std::array<T, LatSet::q> &feq, const Vector<T, LatSet::d> &u,
                          T rho) {
    T u2 = u.getnorm2();
    for (int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order2(k, u, rho, u2);
    }
  }
  // static void SecondOrder_Incompresible(std::array<T, LatSet::q> &feq,
  //                                       const Vector<T, LatSet::d> &u, T rho) {
  //   T u2 = u.getnorm2();
  //   for (int k = 0; k < LatSet::q; ++k) {
  //     feq[k] = Order2_Incompresible(k, u, rho, u2);
  //   }
  // }

  // init cell to equilibrium
  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T)>
  static void InitEquilibrium(BasicCell<T, LatSet> &cell, const T rho,
                              const Vector<T, LatSet::d> &u) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, u, rho);
    for (int i = 0; i < LatSet::q; ++i) cell[i] = feq[i];
  }
  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T)>
  static void InitEquilibrium(Cell<T, LatSet> &cell) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());
    for (int i = 0; i < LatSet::q; ++i) cell[i] = feq[i];
  }
};
