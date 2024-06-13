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

// powerlaw.h

#pragma once

#include "lbm/moment.h"

// m
struct ConsistencyIndexBase : public FieldBase<1> {};
// n-1
struct BehaviorIndexMinus1Base : public FieldBase<1> {};
// max nu
struct MaxNuBase : public FieldBase<1> {};
// min nu
struct MinNuBase : public FieldBase<1> {};

template <typename T>
using ConsistencyIndex = Array<T, ConsistencyIndexBase>;
template <typename T>
using BehaviorIndexMinus1 = Array<T, BehaviorIndexMinus1Base>;
template <typename T>
using MaxNu = Array<T, MaxNuBase>;
template <typename T>
using MinNu = Array<T, MinNuBase>;

template <typename T>
using PowerLawPARAMS = TypePack<ConsistencyIndex<T>, BehaviorIndexMinus1<T>, MinNu<T>, MaxNu<T>>;
// power-law dynamics for non-Newtonian fluid

// in power-law dynamics: compute omega from magnitude of shear rate
template <typename CELLTYPE>
struct PowerLaw_Omega {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static inline T get(T gamma, CELL& cell) {
    T m = cell.template get<ConsistencyIndex<T>>();
    T n_minus1 = cell.template get<BehaviorIndexMinus1<T>>();
    // Ostwald-de Waele/ power-law
    // nu = m * gamma^(n-1)
    T nu = m * std::pow(gamma, n_minus1);
    
    nu = std::min(nu, cell.template get<MaxNu<T>>());
    nu = std::max(nu, cell.template get<MinNu<T>>());
    
    // nu = cs^2 * (tau - 0.5) = (1/omega - 0.5)/3
    T omega = T{1} / (nu * LatSet::InvCs2 + T{0.5});
    return omega;
  }
};

namespace collision {

template <typename EquilibriumScheme, bool WriteToField = false>
struct PowerLaw_BGK_Feq_RhoU {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;

  static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    moment::template rhou<CELL, WriteToField>::apply(cell, rho, u);
    // strain rate
    std::array<T, util::SymmetricMatrixSize<LatSet::d>()> strain_rate{};
    moment::template strainRate<CELL>::get(cell, strain_rate, rho, u);
    // magnitude of shear rate
    T gamma = moment::template shearRateMag<CELL>::get(strain_rate);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // BGK collision
    const T omega = PowerLaw_Omega<CELL>::get(gamma, cell);
    const T _omega = T{1} - omega;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }
};

template <typename EquilibriumScheme, typename ForceScheme, bool WriteToField = false,
          unsigned int dir = 2>
struct PowerLaw_BGKForce_Feq_RhoU {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    moment::template forceRhou<CELL, WriteToField, dir>::apply(
      cell, ForceScheme::getForce(cell), rho, u);
    // strain rate
    std::array<T, util::SymmetricMatrixSize<LatSet::d>()> strain_rate{};
    moment::template strainRate<CELL>::get(cell, strain_rate, rho, u);
    // magnitude of shear rate
    T gamma = moment::template shearRateMag<CELL>::get(strain_rate);
    // compute force term
    std::array<T, LatSet::q> fi{};
    ForceScheme::apply(u, ForceScheme::getForce(cell), fi);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // PowerLaw-BGK collision
    const T omega = PowerLaw_Omega<CELL>::get(gamma, cell);
    const T _omega = T{1} - omega;
    const T fomega = T{1} - T{0.5} * omega;

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }
};

}  // namespace collision