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

// #include "data_struct/Vector.h"
#include "lbm/lattice_set.h"

template <typename T, typename LatSet, typename TypePack>
class Cell;

// calc equilibrium distribution function
// sum(feq_i) = rho, for both first and second order
template <typename T, typename LatSet>
struct Equilibrium {
  __any__ static inline T Order1(int k, const Vector<T, LatSet::d> &u, T rho) {
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }

  __any__ static inline T Order1_Incompresible(int k, const Vector<T, LatSet::d> &u, T rho) {
    return latset::w<LatSet>(k) * (rho + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }

  __any__ static inline T Order2(int k, const Vector<T, LatSet::d> &u, T rho, T u2) {
    T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho *
           (T(1) + LatSet::InvCs2 * uc + uc * uc * T(0.5) * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T(0.5));
  }

  // __any__ static inline T Order2_Incompresible(int k, const Vector<T, LatSet::d> &u, T rho, T
  // u2) {
  //   T uc = u * latset::c<LatSet>(k);
  //   return latset::w<LatSet>(k) * (rho + LatSet::InvCs2 * uc + uc * uc * T(0.5) *
  //   LatSet::InvCs4 -
  //                          LatSet::InvCs2 * u2 * T(0.5));
  // }

  __any__ static void FirstOrder_Incompresible(std::array<T, LatSet::q> &feq, const Vector<T, LatSet::d> &u, T rho) {
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order1_Incompresible(k, u, rho);
    }
  }
  __any__ static void SecondOrder(std::array<T, LatSet::q> &feq, const Vector<T, LatSet::d> &u, T rho) {
    T u2 = u.getnorm2();
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order2(k, u, rho, u2);
    }
  }
};

namespace equilibrium {

template <typename CELL>
struct SecondOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;

  using GenericRho = typename CELL::GenericRho;

  __any__ static inline T get(int k, const Vector<T, LatSet::d> &u, T rho, T u2) {
    const T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho *
           (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T{0.5});
  }

  __any__ static void apply(std::array<T, LatSet::q> &feq, T rho, const Vector<T, LatSet::d> &u) {
    const T u2 = u.getnorm2();
#ifdef UNROLLFOR
    apply_impl(feq, rho, u, u2, std::make_index_sequence<LatSet::q>());
#else
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = get(k, u, rho, u2);
    }
#endif
  }

  __any__ static void apply(std::array<T, LatSet::q> &feq, const CELL &cell) {
    const T rho = cell.template get<GenericRho>();
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
    const T u2 = u.getnorm2();
#ifdef UNROLLFOR
    apply_impl(feq, rho, u, u2, std::make_index_sequence<LatSet::q>());
#else
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = get(k, u, rho, u2);
    }
#endif
  }

  // unroll the for loop
  template <std::size_t... Is>
  __any__ static inline void apply_impl(std::array<T, LatSet::q> &feq, T rho, const Vector<T, LatSet::d> &u, T u2, std::index_sequence<Is...>) {
    // Unpack the sequence and call `get` for each index at compile time
    ((feq[Is] = get<Is>(rho, u, u2)), ...);
  }
  // use template to enable compile-time evaluation of latset::c and latset::w
  template <unsigned int k>
  __any__ static inline T get(T rho, const Vector<T, LatSet::d> &u, T u2) {
    const T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho *
           (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T{0.5});
  }
};

template <typename CELL>
struct FirstOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;

  using GenericRho = typename CELL::GenericRho;

  __any__ static inline T get(int k, const Vector<T, LatSet::d> &u, T rho) {
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }

  __any__ static void apply(std::array<T, LatSet::q> &feq, T rho, const Vector<T, LatSet::d> &u) {
#ifdef UNROLLFOR
    apply_impl(feq, rho, u, std::make_index_sequence<LatSet::q>());
#else
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = get(k, u, rho);
    }
#endif
  }

  __any__ static void apply(std::array<T, LatSet::q> &feq, const CELL &cell) {
    const T rho = cell.template get<GenericRho>();
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
#ifdef UNROLLFOR
    apply_impl(feq, rho, u, std::make_index_sequence<LatSet::q>());
#else
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = get(k, u, rho);
    }
#endif
  }

  // unroll the for loop
  template <std::size_t... Is>
  __any__ static inline void apply_impl(std::array<T, LatSet::q> &feq, T rho, const Vector<T, LatSet::d> &u, std::index_sequence<Is...>) {
    // Unpack the sequence and call `get` for each index at compile time
    ((feq[Is] = get<Is>(rho, u)), ...);
  }
  // use template to enable compile-time evaluation of latset::c and latset::w
  template <unsigned int k>
  __any__ static inline T get(T rho, const Vector<T, LatSet::d> &u) {
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }
};


}  // namespace equilibrium
