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


#include "lbm/moment.ur.h"


namespace collision {


template <typename MomentaScheme, typename EquilibriumScheme>
struct BGK {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    MomentaScheme::apply(cell, rho, u);
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

template <typename MomentaScheme, typename EquilibriumScheme, typename ForceScheme>
struct BGKForce {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    const auto force = ForceScheme::getForce(cell);
    MomentaScheme::apply(cell, force, rho, u);
    // compute force term
    std::array<T, LatSet::q> fi{};
    ForceScheme::apply(u, force, fi);
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

  __any__ static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const auto source = cell.template get<SOURCE>();
    moment::template sourceRho<CELL, SOURCE, WriteToField>::apply(cell, rho, source);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * source * latset::w<LatSet>(i);
    }
  }
};

// full way bounce back, could be regarded as a mpdified collision process
// swap the populations in the opposite direction
template <typename CELLTYPE>
struct BounceBack {
  using CELL = CELLTYPE;
  using LatSet = typename CELL::LatticeSet;
  using T = typename CELL::FloatType;
  using GenericRho = typename CELL::GenericRho;
  static constexpr unsigned int startdir = LatSet::q % 2 == 0 ? 0 : 1;

  __any__ static void apply(CELL& cell) {
    for (unsigned int i = startdir; i < LatSet::q; i += 2) {
      const T temp = cell[i];
      const unsigned int iopp = i + 1;
      cell[i] = cell[iopp];
      cell[iopp] = temp;
    }
  }
  
};

  // static inline void apply(CELL &cell, unsigned int k) {
  //   cell[k] = cell.getPrevious(latset::opp<LatSet>(k)) +
  //             2 * LatSet::InvCs2 * latset::w<LatSet>(k) * cell.template get<GenericRho>() *
  //               (cell.template get<VELOCITY<T, LatSet::d>>() * latset::c<LatSet>(k));
  // }

// full way bounce back with moving wall, could be regarded as a modified collision process
// swap the populations in the opposite direction
template <typename CELLTYPE>
struct BounceBackMovingWall {
  using CELL = CELLTYPE;
  using LatSet = typename CELL::LatticeSet;
  using T = typename CELL::FloatType;
  using GenericRho = typename CELL::GenericRho;
  static constexpr unsigned int startdir = LatSet::q % 2 == 0 ? 0 : 1;

  __any__ static void apply(CELL& cell) {
    const T rhox = 2 * LatSet::InvCs2 * cell.template get<GenericRho>();
    for (unsigned int i = startdir; i < LatSet::q; i += 2) {
      const T temp = cell[i];
      const unsigned int iopp = i + 1;
      const T uc = cell.template get<VELOCITY<T, LatSet::d>>() * latset::c<LatSet>(i) * latset::w<LatSet>(i) * rhox;
      cell[i] = cell[iopp] + uc;
      cell[iopp] = temp - uc;
    }
  }
  
};

}  // namespace collision


// old version of BGK collision

namespace legacy {

template <typename T, typename LatSet>
struct BGK {
  // BGK collision operator
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  __any__ static void apply(PopCell<T, LatSet>& cell) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }

  // BGK collision operator with force
  // template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  // __any__ static void applyForce(PopCell<T, LatSet>& cell, const Vector<T, LatSet::d>& force) {
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
  __any__ static void applySource(PopCell<T, LatSet>& cell, const std::array<T, LatSet::q>& fi) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  __any__ static void applySource(PopCell<T, LatSet>& cell, const T S) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * S * latset::w<LatSet>(i);
    }
  }
};

}  // namespace collision