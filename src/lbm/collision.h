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


#include "lbm/moment.h"


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
  using GenericRho = typename CELL::GenericRho;

  static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    moment::template rhou<CELL, WriteToField>::apply(cell, rho, u);
    // moment::template rhou<CELL, WriteToField>::apply(cell);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // EquilibriumScheme::apply(feq, cell);
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
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
  using GenericRho = typename CELL::GenericRho;

  static void apply(CELL& cell) {
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, cell);
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }
};

// a typical BGK collision process with:
// macroscopic variables updated
// equilibrium distribution function calculated
// force term
template <typename EquilibriumScheme, typename ForceScheme, bool WriteToField = false,
          unsigned int dir = 2>
struct BGKForce_Feq_RhoU {
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
    // compute force term
    std::array<T, LatSet::q> fi{};
    ForceScheme::apply(u, ForceScheme::getForce(cell), fi);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }
};

// a typical BGK collision process with:
// equilibrium distribution function calculated
// force term
template <typename EquilibriumScheme, typename ForceScheme>
struct BGKForce_Feq {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  static void apply(CELL& cell) {
    // compute force term
    std::array<T, LatSet::q> fi{};
    ForceScheme::apply(cell, fi);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, cell);
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }
};

// a typical BGK collision process with:
// macroscopic variables updated
// equilibrium distribution function calculated
// force term
template <typename EquilibriumScheme, typename SOURCE, bool WriteToField = false>
struct BGKSource_Feq_Rho {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const auto source = cell.template get<SOURCE>();
    moment::template sourceRho<CELL, WriteToField>::apply(cell, rho, source);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * source * LatSet::w[i];
    }
  }
};

// full way bounce back, could be regarded as a mpdified collision process
// swap the populations in the opposite direction
// LatSet must have rest population(D2Q4 is not supported)
template <typename CELLTYPE>
struct BounceBack {
  using CELL = CELLTYPE;
  using LatSet = typename CELL::LatticeSet;
  using T = typename CELL::FloatType;
  using GenericRho = typename CELL::GenericRho;
  static constexpr int startdir = LatSet::q % 2 == 0 ? 0 : 1;

  static void apply(CELL& cell) {
    for (int i = startdir; i < LatSet::q; i += 2) {
      T temp = cell[i];
      int iopp = i + 1;
      cell[i] = cell[iopp];
      cell[iopp] = temp;
    }
  }
  
};

  // static inline void apply(CELL &cell, unsigned int k) {
  //   cell[k] = cell.getPrevious(LatSet::opp[k]) +
  //             2 * LatSet::InvCs2 * LatSet::w[k] * cell.template get<GenericRho>() *
  //               (cell.template get<VELOCITY<T, LatSet::d>>() * LatSet::c[k]);
  // }

// full way bounce back with moving wall, could be regarded as a modified collision process
// swap the populations in the opposite direction
// LatSet must have rest population(D2Q4 is not supported)
template <typename CELLTYPE>
struct BounceBackMovingWall {
  using CELL = CELLTYPE;
  using LatSet = typename CELL::LatticeSet;
  using T = typename CELL::FloatType;
  using GenericRho = typename CELL::GenericRho;
  static constexpr int startdir = LatSet::q % 2 == 0 ? 0 : 1;

  static void apply(CELL& cell) {
    T rhox = 2 * LatSet::InvCs2 * cell.template get<GenericRho>();
    for (int i = startdir; i < LatSet::q; i += 2) {
      T temp = cell[i];
      int iopp = i + 1;
      T uc = cell.template get<VELOCITY<T, LatSet::d>>() * LatSet::c[i] * LatSet::w[i] * rhox;
      cell[i] = cell[iopp] + uc;
      cell[iopp] = temp - uc;
    }
  }
  
};


// old version of BGK collision

template <typename T, typename LatSet>
struct BGK {
  // BGK collision operator
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  static void apply(PopCell<T, LatSet>& cell) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();

    for (int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }

  // BGK collision operator with force
  // template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  // static void applyForce(PopCell<T, LatSet>& cell, const Vector<T, LatSet::d>& force) {
  //   std::array<T, LatSet::q> feq{};
  //   GetFeq(feq, cell.getVelocity(), cell.getRho());

  //   std::array<T, LatSet::q> fi{};
  //   force::ForcePop<T, LatSet>::compute(fi, cell.getVelocity(), force);

  //   const T omega = cell.getOmega();
  //   const T _omega = cell.get_Omega();
  //   const T fomega = cell.getfOmega();

  //   for (int i = 0; i < LatSet::q; ++i) {
  //     cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
  //   }
  // }

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  static void applySource(PopCell<T, LatSet>& cell, const std::array<T, LatSet::q>& fi) {
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
  static void applySource(PopCell<T, LatSet>& cell, const T S) {
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