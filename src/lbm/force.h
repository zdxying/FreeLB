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

// force.h
#pragma once

#include "lbm/moment.ur.h"
#include "utils/util.h"

namespace force {

template <typename T, typename LatSet>
struct ForcePopImpl {
  __any__ static inline void compute(std::array<T, LatSet::q> &Fi, const Vector<T, LatSet::d> &u, const Vector<T, LatSet::d> &F) {
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      Fi[i] = latset::w<LatSet>(i) * F * ((latset::c<LatSet>(i) - u) * LatSet::InvCs2 + (latset::c<LatSet>(i) * u * LatSet::InvCs4) * latset::c<LatSet>(i));
    }
  }
};

template <typename T, typename LatSet, unsigned int d>
struct ScalarForcePopImpl {
  __any__ static inline void compute(std::array<T, LatSet::q> &Fi, const Vector<T, LatSet::d> &u, const T F) {
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      const T v1 = (latset::c<LatSet>(i)[d] - u[d]) * LatSet::InvCs2;
      const T v2 = (latset::c<LatSet>(i) * u * LatSet::InvCs4) * latset::c<LatSet>(i)[d];
      Fi[i] = latset::w<LatSet>(i) * F * (v1 + v2);
    }
  }
};

template <typename T, typename LatSet>
struct ForcePop {

  // calculate discrete force
  __any__ static inline void compute(std::array<T, LatSet::q> &Fi, const Vector<T, LatSet::d> &u, const Vector<T, LatSet::d> &F) {
    ForcePopImpl<T, LatSet>::compute(Fi, u, F);
  }

  template <unsigned int d>
  __any__ static inline void computeScalar(std::array<T, LatSet::q> &Fi, const Vector<T, LatSet::d> &u, const T F) {
    ScalarForcePopImpl<T, LatSet, d>::compute(Fi, u, F);
  }
};

template <typename CELL>
struct Force {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // scalardir will NOT be used in Force<>
  static constexpr unsigned int scalardir = 0;

  __any__ static void apply(CELL &cell, std::array<T, LatSet::q> &Fi) {
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
    const Vector<T, LatSet::d> &F = cell.template get<FORCE<T, LatSet::d>>();
    ForcePop<T, LatSet>::compute(Fi, u, F);
  }
  __any__ static void apply(const Vector<T, LatSet::d> &u, const Vector<T, LatSet::d> &F,
                    std::array<T, LatSet::q> &Fi) {
    ForcePop<T, LatSet>::compute(Fi, u, F);
  }
  __any__ static Vector<T, LatSet::d> &getVectorForce(CELL &cell) {
    return cell.template get<FORCE<T, LatSet::d>>();
  }
  __any__ static auto &getForce(CELL &cell) { return cell.template get<FORCE<T, LatSet::d>>(); }
};

template <typename CELL>
struct ConstForce {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // scalardir will NOT be used in ConstForce<>
  static constexpr unsigned int scalardir = 0;

  __any__ static void apply(CELL &cell, std::array<T, LatSet::q> &Fi) {
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
    const Vector<T, LatSet::d> &F = cell.template get<CONSTFORCE<T, LatSet::d>>();
    ForcePop<T, LatSet>::compute(Fi, u, F);
  }
  __any__ static void apply(const Vector<T, LatSet::d> &u, const Vector<T, LatSet::d> &F, std::array<T, LatSet::q> &Fi) {
    ForcePop<T, LatSet>::compute(Fi, u, F);
  }
  __any__ static Vector<T, LatSet::d> &getVectorForce(CELL &cell) {
    return cell.template get<CONSTFORCE<T, LatSet::d>>();
  }
  __any__ static auto &getForce(CELL &cell) {
    return cell.template get<CONSTFORCE<T, LatSet::d>>();
  }
};

template <typename CELL, unsigned int dir = 2>
struct ScalarForce {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

  __any__ static void apply(CELL &cell, std::array<T, LatSet::q> &Fi) {
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
    const T F = cell.template get<SCALARFORCE<T>>();
    ForcePop<T, LatSet>::template computeScalar<scalardir>(Fi, u, F);
  }
  __any__ static void apply(const Vector<T, LatSet::d> &u, T F, std::array<T, LatSet::q> &Fi) {
    ForcePop<T, LatSet>::template computeScalar<scalardir>(Fi, u, F);
  }
  __any__ static Vector<T, LatSet::d> getVectorForce(const CELL &cell) {
    if constexpr (LatSet::d == 2) {
      return Vector<T, LatSet::d>{T{}, cell.template get<SCALARFORCE<T>>()};
    } else if constexpr (LatSet::d == 3) {
      return Vector<T, LatSet::d>{T{}, T{}, cell.template get<SCALARFORCE<T>>()};
    }
  }
  __any__ static auto getForce(const CELL &cell) {
    return cell.template get<SCALARFORCE<T>>();
  }
};

template <typename CELL, unsigned int dir = 2>
struct ScalarConstForce {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

  __any__ static void apply(CELL &cell, std::array<T, LatSet::q> &Fi) {
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
    const T F = cell.template get<SCALARCONSTFORCE<T>>();
    ForcePop<T, LatSet>::template computeScalar<scalardir>(Fi, u, F);
  }
  __any__ static void apply(const Vector<T, LatSet::d> &u, T F, std::array<T, LatSet::q> &Fi) {
    ForcePop<T, LatSet>::template computeScalar<scalardir>(Fi, u, F);
  }
  __any__ static Vector<T, LatSet::d> getVectorForce(const CELL &cell) {
    if constexpr (LatSet::d == 2) {
      return Vector<T, LatSet::d>{T{}, cell.template get<SCALARCONSTFORCE<T>>()};
    } else if constexpr (LatSet::d == 3) {
      return Vector<T, LatSet::d>{T{}, T{}, cell.template get<SCALARCONSTFORCE<T>>()};
    }
  }
  __any__ static auto getForce(const CELL &cell) {
    return cell.template get<SCALARCONSTFORCE<T>>();
  }
};

template <typename CELLType0, typename CELLType1, bool ScalarForce = true>
struct Buoyancy {
  using CELL = CELLType0;
  using CELL0 = CELLType0;
  using CELL1 = CELLType1;

  using T = typename CELL0::FloatType;
  using LatSet0 = typename CELL0::LatticeSet;
  using LatSet1 = typename CELL1::LatticeSet;
  static constexpr unsigned int d = LatSet0::d - 1;

  using GenericRho = typename CELL1::GenericRho;

  // CELL0: NS cell, CELL1: thermal or solute cell
  __any__ static void apply(CELL0 &cell0, CELL1 &cell1) {
    const T rho = cell1.template get<GenericRho>();
    const T rhoInit = cell1.template get<RHOINIT<T>>();
    const T gbeta = cell1.template get<GBETA<T>>();
    if constexpr (ScalarForce) {
      cell0.template get<SCALARFORCE<T>>() += gbeta * (rho - rhoInit);
    } else {
      cell0.template get<FORCE<T, LatSet0::d>>()[d] += gbeta * (rho - rhoInit);
    }
  }
};


}  // namespace force