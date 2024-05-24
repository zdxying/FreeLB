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

// power-law dynamics for non-Newtonian fluid

// in power-law dynamics: compute omega from magnitude of shear rate
template <typename CELLTYPE>
struct omegaFromShearRate {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static inline T computeOmega(const T gamma, const T m, const T n_minus1) {
    // Ostwald-de Waele/ power-law
    // nu = m * gamma^(n-1)
    T nu = m * std::pow(gamma, n_minus1);
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

  static void apply(CELL& cell, const T m, const T n_minus1) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    moment::template rhou<CELL, WriteToField>::apply(cell, rho, u);
    // strain rate
    std::array<T, util::SymmetricMatrixSize<LatSet::d>> strain_rate{};
    moment::template strainRate<CELL>::get(cell, strain_rate, rho, u);
    // magnitude of shear rate
    T gamma = moment::template shearRateMag<CELL>::get(strain_rate);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(cell, feq, rho, u);
    // BGK collision
    const T omega = omegaFromShearRate<CELL>::computeOmega(gamma, m, n_minus1);
    const T _omega = T{1} - omega;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }
};

}  // namespace collision