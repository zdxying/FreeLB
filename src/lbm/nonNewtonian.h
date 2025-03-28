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

// nonNewtonian.h

#pragma once

#include "lbm/moment.ur.h"


template <typename T>
using PowerLawPARAMS = TypePack<PL_k<T>, PL_m<T>, MinShearRate<T>, MaxShearRate<T>, OMEGA<T>>;

// power-law dynamics: compute omega from magnitude of shear rate
template <typename CELLTYPE>
struct PowerLaw_Omega {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // in lattice unit
  static inline T get(T gamma, CELL& cell) {
    const T gamma_min = cell.template get<MinShearRate<T>>();
    const T gamma_max = cell.template get<MaxShearRate<T>>();
    const T k = cell.template get<PL_k<T>>();
    const T m = cell.template get<PL_m<T>>();

    if (gamma < gamma_min) gamma = gamma_min;
    else if (gamma > gamma_max) gamma = gamma_max;
    // Ostwald-de Waele/ power-law: nu = k * gamma^(m)
    const T nu = k * std::pow(gamma, m);
    // lat_nu = cs^2 * (tau - 0.5) = (1/omega - 0.5)/3
    const T omega = T{1} / (nu * LatSet::InvCs2 + T{0.5});
    cell.template get<OMEGA<T>>() = omega;
    return omega;
  }
};


template <typename T>
using CrossPARAMS = TypePack<Cross_eta0<T>, Cross_t<T>, Cross_m<T>, MinShearRate<T>, MaxShearRate<T>, OMEGA<T>>;

template <typename CELLTYPE>
struct Cross_Omega {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // in lattice unit
  static inline T get(T gamma, CELL& cell) {
    const T gamma_min = cell.template get<MinShearRate<T>>();
    const T gamma_max = cell.template get<MaxShearRate<T>>();
    const T cross_eta0 = cell.template get<Cross_eta0<T>>();
    const T cross_t = cell.template get<Cross_t<T>>();
    const T cross_m = cell.template get<Cross_m<T>>();

    if (gamma < gamma_min) gamma = gamma_min;
    else if (gamma > gamma_max) gamma = gamma_max;
    // cross
    const T nu = cross_eta0/(1+std::pow(cross_t*gamma, cross_m));
    // lat_nu = cs^2 * (tau - 0.5) = (1/omega - 0.5)/3
    const T omega = T{1} / (nu * LatSet::InvCs2 + T{0.5});
    cell.template get<OMEGA<T>>() = omega;
    return omega;
  }
};


namespace collision {

template <typename MomentaScheme, typename EquilibriumScheme, typename RheologyOemga>
struct Rheology_BGK {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;

  static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    MomentaScheme::apply(cell, rho, u);
    // strain rate
    std::array<T, util::SymmetricMatrixSize<LatSet::d>()> strain_rate{};
    moment::template strainRate<CELL>::apply(cell, rho, u, strain_rate);
    // magnitude of shear rate
    T gamma = moment::template shearRateMag<CELL>::get(strain_rate);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // Rheology omega
    const T omega = RheologyOemga::get(gamma, cell);
    const T _omega = T{1} - omega;

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }
};

template <typename MomentaScheme, typename EquilibriumScheme, typename ForceScheme, typename RheologyOemga>
struct Rheology_BGKForce {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    const auto force = ForceScheme::getForce(cell);
    MomentaScheme::apply(cell, force, rho, u);
    // strain rate
    std::array<T, util::SymmetricMatrixSize<LatSet::d>()> strain_rate{};
    moment::template strainRate<CELL>::apply(cell, rho, u, strain_rate);
    // magnitude of shear rate
    T gamma = moment::template shearRateMag<CELL>::get(strain_rate);
    // compute force term
    std::array<T, LatSet::q> fi{};
    ForceScheme::apply(u, force, fi);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // Rheology omega
    const T omega = RheologyOemga::get(gamma, cell);
    const T _omega = T{1} - omega;
    const T fomega = T{1} - T{0.5} * omega;

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }
};

}  // namespace collision