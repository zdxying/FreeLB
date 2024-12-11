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

// collisionMRT.h

#pragma once


#include "lbm/moment.ur.h"
#include "lbm/collisionMRT.h"

#ifdef _UNROLLFOR

namespace collision {

#ifdef __CUDA_ARCH__
template <typename T, typename LatSet, typename TypePack>
using CELL = cudev::Cell<T, LatSet, TypePack>;

#else
template <typename T, typename LatSet, typename TypePack>
using CELL = Cell<T, LatSet, TypePack>;
#endif

// a typical MRT collision process with:
// macroscopic variables updated
// equilibrium distribution function calculated
template <typename T, typename TypePack, bool WriteToField>
struct MRT_Feq_RhoU<CELL<T, D2Q9<T> ,TypePack>, WriteToField> {
  using LatSet = D2Q9<T>;
  using equilibriumscheme = equilibrium::SecondOrder<CELL<T, D2Q9<T> ,TypePack>>;
  using GenericRho = typename CELL<T, D2Q9<T> ,TypePack>::GenericRho;

  __any__ static void apply(CELL<T, D2Q9<T> ,TypePack>& cell) {
    const T omega = cell.getOmega();
    // relaxation time vector
    const T rtvec[LatSet::q] {T{}, T{11./10.}, T{11./10.}, T{}, T{11./10.}, T{}, T{11./10.}, omega, omega};
    // relaxation time matrix
    T InvM_S[LatSet::q][LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        InvM_S[i][j] = mrt::InvM<LatSet>(i, j) * rtvec[j];
      }
    }

    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    moment::template rhoU<CELL<T, D2Q9<T> ,TypePack>, WriteToField>::apply(cell, rho, u);

    // compute Momenta
    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      momenta[i] = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];
      }
    }
    // compute Equilibrium
    // std::array<T, LatSet::q> feq{};
    // equilibrium::SecondOrder<CELL<T, D2Q9<T> ,TypePack>>::apply(feq, rho, u);
    // T momentaEq[LatSet::q] {};
    // for(unsigned int i = 0; i < LatSet::q; ++i) {
    //   momentaEq[i] = T{};
    //   for (unsigned int j = 0; j < LatSet::q; ++j) {
    //     momentaEq[i] += mrt::M<LatSet>(i, j) * feq[j];
    //   }
    // }
    // T momentaEq[LatSet::q] {};
    // const T ux2 = u[0] * u[0];
    // const T uy2 = u[1] * u[1];
    // const T rhoux = rho * u[0];
    // const T rhouy = rho * u[1];
    // const T cse0 = 3 * rho * (ux2 + uy2);
    // momentaEq[0] = rho;
    // momentaEq[1] = -2 * rho + cse0;
    // momentaEq[2] = 9 * rho * ux2 * uy2 - cse0 + rho;
    // momentaEq[3] = rhoux;
    // momentaEq[4] = rhoux * (3 * uy2 - 1);
    // momentaEq[5] = rhouy;
    // momentaEq[6] = rhouy * (3 * ux2 - 1);
    // momentaEq[7] = rho * (ux2 - uy2);
    // momentaEq[8] = rho * u[0] * u[1];

    // compute off-Equilibrium part: momenta - momentaEq
    T delmomenta[LatSet::q] {};
    const T ux2 = u[0] * u[0];
    const T uy2 = u[1] * u[1];
    const T rhoux = rho * u[0];
    const T rhouy = rho * u[1];
    const T cse0 = 3 * rho * (ux2 + uy2);
    delmomenta[0] = momenta[0] - rho;
    delmomenta[1] = momenta[1] + 2 * rho - cse0;
    delmomenta[2] = momenta[2] - 9 * rho * ux2 * uy2 + cse0 - rho;
    delmomenta[3] = momenta[3] - rhoux;
    delmomenta[4] = momenta[4] - rhoux * (3 * uy2 - 1);
    delmomenta[5] = momenta[5] - rhouy;
    delmomenta[6] = momenta[6] - rhouy * (3 * ux2 - 1);
    delmomenta[7] = momenta[7] - rho * (ux2 - uy2);
    delmomenta[8] = momenta[8] - rho * u[0] * u[1];

    // MRT collision
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T momentaEqS{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEqS += InvM_S[i][j] * delmomenta[j];
      }
      cell[i] -= momentaEqS;
    }
  }

};

}

#endif