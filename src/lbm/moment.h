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

// moment.h

#pragma once

#include "data_struct/cell.h"

namespace moment {

template <typename CELLTYPE, bool WriteToField = false>
struct rho {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  using GenericRho = typename CELL::GenericRho;

  static inline T get(CELL& cell) {
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
    if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
    return rho_value;
  }
  static inline void apply(CELL& cell, T& rho_value) {
    rho_value = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
    if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
  }
  // always write to field
  static inline void apply(CELL& cell) {
    T& rho_value = cell.template get<GenericRho>();
    rho_value = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
  }
};

// update rho(usually temperature or concentration in advection-diffusion problems) with
// source term, no need to preprocess the source term
template <typename CELLTYPE, bool WriteToField = false>
struct sourceRho {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  using GenericRho = typename CELL::GenericRho;

  static inline T get(CELL& cell, T source) {
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
    // fOmega: avoid lattice artifact
    rho_value += source * T{0.5};
    if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
    return rho_value;
  }
  static inline void apply(CELL& cell, T& rho_value, T source) {
    rho_value = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
    rho_value += source * T{0.5};
    if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
  }
  // always write to field
  static inline void apply(CELL& cell, T source) {
    T& rho_value = cell.template get<GenericRho>();
    rho_value = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
    rho_value += source * T{0.5};
  }
};

template <typename CELLTYPE, bool WriteToField = false>
struct u {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static inline Vector<T, LatSet::d> get(CELL& cell) {
    Vector<T, LatSet::d> u_value;
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value /= rho_value;
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
    return u_value;
  }
  static inline void apply(CELL& cell, Vector<T, LatSet::d>& u_value) {
    u_value.clear();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value /= rho_value;
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
  }
  // always write to field
  static inline void apply(CELL& cell) {
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value /= rho_value;
  }
};

template <typename CELLTYPE, bool WriteToField = false, unsigned int dir = 2>
struct forceU {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

  static inline Vector<T, LatSet::d> get(CELL& cell,
                                         const Vector<T, LatSet::d>& f_alpha) {
    Vector<T, LatSet::d> u_value;
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value += f_alpha * T{0.5};
    u_value /= rho_value;
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
    return u_value;
  }
  static inline Vector<T, LatSet::d> get(CELL& cell, T f) {
    Vector<T, LatSet::d> u_value;
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value[scalardir] += f * T{0.5};
    u_value /= rho_value;
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
    return u_value;
  }
  static inline void apply(CELL& cell, Vector<T, LatSet::d>& u_value,
                           const Vector<T, LatSet::d>& f_alpha) {
    u_value.clear();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value += f_alpha * T{0.5};
    u_value /= rho_value;
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
  }
  static inline void apply(CELL& cell, Vector<T, LatSet::d>& u_value, T f) {
    u_value.clear();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value[scalardir] += f * T{0.5};
    u_value /= rho_value;
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
  }
  // always write to field
  static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha) {
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value += f_alpha * T{0.5};
    u_value /= rho_value;
  }
  static inline void apply(CELL& cell, T f) {
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value[scalardir] += f * T{0.5};
    u_value /= rho_value;
  }
};

template <typename CELLTYPE, bool WriteToField = false>
struct rhou {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL& cell, T& rho_value, Vector<T, LatSet::d>& u_value) {
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value /= rho_value;
    if constexpr (WriteToField) {
      cell.template get<GenericRho>() = rho_value;
      cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
    }
  }
  // will write to field regardless of WriteToField value
  static inline void apply(CELL& cell) {
    T& rho_value = cell.template get<GenericRho>();
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value /= rho_value;
  }
};

template <typename CELLTYPE, bool WriteToField = false, unsigned int dir = 2>
struct forceRhou {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  using GenericRho = typename CELL::GenericRho;

  static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

  static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value,
                           Vector<T, LatSet::d>& u_value) {
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value += f_alpha * T{0.5};
    u_value /= rho_value;
    if constexpr (WriteToField) {
      cell.template get<GenericRho>() = rho_value;
      cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
    }
  }
  // for scalar force
  static inline void apply(CELL& cell, T f, T& rho_value, T u_value) {
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value[scalardir] += f * T{0.5};
    u_value /= rho_value;
    if constexpr (WriteToField) {
      cell.template get<GenericRho>() = rho_value;
      cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
    }
  }
  // always write to field
  static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha) {
    T& rho_value = cell.template get<GenericRho>();
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value += f_alpha * T{0.5};
    u_value /= rho_value;
  }
  static inline void apply(CELL& cell, T f) {
    T& rho_value = cell.template get<GenericRho>();
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += LatSet::c[i] * cell[i];
    }
    u_value[scalardir] += f * T{0.5};
    u_value /= rho_value;
  }
};

// strain rate tensor/ rate of deformation matrix
template <typename CELLTYPE>
struct strainRate {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // S_ab = (-3/(2*tau))*SUM_i(f_i^(1)*c_ia*c_ib)
  // f_i^(1) = f_i - f_i^eq
  // PI_ab^eq = SUM_i(f_i^eq*c_ia*c_ib) = rho*U_a*U_b + rho*Cs^2*delta_ab
  static inline void get(CELL& cell,
                         std::array<T, util::SymmetricMatrixSize<LatSet::d>>& tensor,
                         const T rho, const Vector<T, LatSet::d>& u) {
    unsigned int i{};
    const T coeff = T{-1.5} * cell.getOmega();
    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        T value{};
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          value += LatSet::c[k][alpha] * LatSet::c[k][beta] * cell[k];
        }
        // remove the equilibrium part: PI_ab^eq
        value -= rho * u[alpha] * u[beta];
        if (alpha == beta) value -= rho * LatSet::cs2;
        // multiply by the coefficient
        value *= coeff;
        tensor[i] = value;
        ++i;
      }
    }
  }
};

// magnitude of shear rate
template <typename CELLTYPE>
struct shearRateMag {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // \dot{gamma} = sqrt(2*trace(S^2)) = sqrt(2*SUM_{a,b}S_ab^2)
  static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>>& tensor) {
    T value{};
    for (unsigned int i = 0; i < util::SymmetricMatrixSize<LatSet::d>; ++i) {
      value += tensor[i] * tensor[i];
    }
    return std::sqrt(T{2} * value);
  }
};

template <typename T, typename LatSet>
struct Rho {
  // return rho
  static T get(const BasicCell<T, LatSet>& cell) {
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho += cell[i];
    return rho;
  }

  // return rho with source
  static T get(const BasicCell<T, LatSet>& cell, const T source) {
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho += cell[i];
    rho += source * T{0.5};
    return rho;
  }

  // compute rho
  static void apply(const BasicCell<T, LatSet>& cell, T& rho) {
    rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho += cell[i];
  }

  // compute rho with source: C = sum(f) + q/2
  static void apply(const BasicCell<T, LatSet>& cell, T& rho, const T source) {
    rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho += cell[i];
    rho += source * T{0.5};
  }
};

template <typename T, typename LatSet>
struct Velocity {
  // get velocity
  static Vector<T, LatSet::d> get(const BasicCell<T, LatSet>& cell) {
    Vector<T, LatSet::d> u;
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + LatSet::c[i] * cell[i];
    }
    return u / rho;
  }

  // get velocity with force
  static Vector<T, LatSet::d> get(const BasicCell<T, LatSet>& cell,
                                  const std::array<T, LatSet::q>& fi) {
    Vector<T, LatSet::d> u;
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + (LatSet::c[i] * (cell[i] + T{0.5} * fi[i]));
    }
    return u / rho;
  }

  // compute velocity
  static void apply(const BasicCell<T, LatSet>& cell, Vector<T, LatSet::d>& u) {
    u.clear();
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + LatSet::c[i] * cell[i];
    }
    u = u / rho;
  }

  // compute velocity with force
  static void apply(const BasicCell<T, LatSet>& cell, Vector<T, LatSet::d>& u,
                    const std::array<T, LatSet::q>& fi) {
    u.clear();
    T rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + (LatSet::c[i] * (cell[i] + T{0.5} * fi[i]));
    }
    u = u / rho;
  }
};

template <typename T, typename LatSet>
struct RhoVelocity {
  // compute rho and velocity
  static void apply(const BasicCell<T, LatSet>& cell, T& rho, Vector<T, LatSet::d>& u) {
    rho = T(0);
    u.clear();
    for (int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + LatSet::c[i] * cell[i];
    }
    u = u / rho;
  }
};

}  // namespace moment