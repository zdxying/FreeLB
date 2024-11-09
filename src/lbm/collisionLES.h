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

// collisionLES.h

#pragma once

#include "lbm/moment.ur.h"

namespace collision {

template <typename MomentaScheme, typename EquilibriumScheme, bool WriteToField = false>
struct SmagorinskyBGK {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;

  static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    MomentaScheme::apply(cell, rho, u);

    // get Smagorinsky Effective Omega
    // second moment of non-equilibrium part of the distribution function
    std::array<T, util::SymmetricMatrixSize<LatSet::d>()> Pi_ab{};
    moment::template Pi_ab_neq<CELL>::get(cell, Pi_ab, rho, u);
    // nrom of Pi_ab
    T Pi_ab_norm = std::sqrt(util::NormSquare<T, LatSet::d>(Pi_ab));
    // coefficient
    const T smagorinsky = cell.template get<SMAGORINSKY<T>>();
    T coeff = smagorinsky * smagorinsky * LatSet::InvCs2 * LatSet::InvCs2 * 2 * std::sqrt(2) / cell.template get<RHO<T>>();
    // Molecular realaxation time
    const T RT_mol = T{1} / cell.getOmega();
    // Turbulent effective time
    const T RT_tur = T{0.5} * (std::sqrt(RT_mol * RT_mol + coeff * Pi_ab_norm) - RT_mol);
    // effective omega
    const T Eomega = T{1} / (RT_mol + RT_tur);

    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // BGK collision
    const T _Eomega = T{1} - Eomega;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = Eomega * feq[i] + _Eomega * cell[i];
    }
  }
};

template <typename MomentaScheme, typename EquilibriumScheme, typename ForceScheme, bool WriteToField = false>
struct SmagorinskyForceBGK {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;

  static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    const auto force = ForceScheme::getForce(cell);
    MomentaScheme::apply(cell, force, rho, u);
    // compute force term
    std::array<T, LatSet::q> fi{};
    ForceScheme::apply(u, force, fi);

    // get Smagorinsky Effective Omega
    // second moment of non-equilibrium part of the distribution function
    std::array<T, util::SymmetricMatrixSize<LatSet::d>()> Pi_ab{};
    moment::template forcePi_ab_neq<CELL>::get(cell, Pi_ab, rho, u, force);
    // nrom of Pi_ab
    T Pi_ab_norm = std::sqrt(util::NormSquare<T, LatSet::d>(Pi_ab));
    // coefficient
    const T smagorinsky = cell.template get<SMAGORINSKY<T>>();
    const T coeff = smagorinsky * smagorinsky * LatSet::InvCs2 * LatSet::InvCs2 * 2 * std::sqrt(2) / rho;
    // const T coeff = smagorinsky * smagorinsky * LatSet::InvCs2 * std::sqrt(2);
    // Molecular realaxation time
    const T RT_mol = T{1} / cell.getOmega(); // cell.template get<OMEGA<T>>();
    // Turbulent effective time
    const T RT_tur = T{0.5} * (std::sqrt(RT_mol * RT_mol + coeff * Pi_ab_norm) - RT_mol);
    // const T RT_tur = coeff * Pi_ab_norm;
    // effective omega
    const T Eomega = T{1} / (RT_mol + RT_tur);
    // cell.template get<OMEGA<T>>() = Eomega;

    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // BGK collision
    const T _Eomega = T{1} - Eomega;
    const T fEomega = T{1} - T{0.5} * Eomega;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = Eomega * feq[i] + _Eomega * cell[i] + fEomega * fi[i];
    }
  }
};

template <typename EquilibriumScheme, typename ForceScheme, bool WriteToField = false>
struct PalabosSmagorinskyForceBGK_Feq_RhoU {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;

  static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    moment::template forceRhou<CELL, ForceScheme, WriteToField>::apply(cell, ForceScheme::getForce(cell), rho, u);
    // compute force term
    std::array<T, LatSet::q> fi{};
    ForceScheme::apply(u, ForceScheme::getForce(cell), fi);

    // get Smagorinsky Effective Omega
    // second moment of non-equilibrium part of the distribution function
    std::array<T, util::SymmetricMatrixSize<LatSet::d>()> Pi_ab{};
    moment::template Pi_ab_neq<CELL>::get(cell, Pi_ab, rho, u);
    // nrom of Pi_ab
    const T Pi_ab_sqr = util::NormSquare<T, LatSet::d>(Pi_ab);
    const T Pi_ab_norm = std::sqrt(Pi_ab_sqr);
    // coefficient
    const T omega = cell.getOmega();
    const T smagorinsky = cell.template get<SMAGORINSKY<T>>();
    const T coeff = (smagorinsky * LatSet::InvCs2 * omega) * (smagorinsky * LatSet::InvCs2 * omega) * T{0.5} / rho;
    const T linearTerm = coeff * Pi_ab_norm;
    const T squareTerm = T{2} * coeff * coeff * Pi_ab_sqr;

    const T Eomega = omega * (1 - linearTerm + squareTerm);
    // cell.template get<OMEGA<T>>() = Eomega;

    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // BGK collision
    const T _Eomega = T{1} - Eomega;
    const T fEomega = T{1} - T{0.5} * Eomega;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = Eomega * feq[i] + _Eomega * cell[i] + fEomega * fi[i];
    }
  }
};

}