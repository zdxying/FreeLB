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

#pragma once

#include "data_struct/cell.h"

namespace moment {


// ----------------------------------------------------------------------------
// empty momenta
template <typename CELLTYPE>
struct NoMomenta {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell, T& rho, Vector<T, LatSet::d>& u) {
  }
};


// ----------------------------------------------------------------------------
// moment from field value

template <typename CELLTYPE>
struct useFieldrho {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline void apply(CELL& cell, T& rho) {
    rho = cell.template get<GenericRho>();
  }
};

template <typename CELLTYPE>
struct useFieldU {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline void apply(CELL& cell, Vector<T, LatSet::d>& u) {
    u = cell.template get<VELOCITY<T, LatSet::d>>();
  }
  __any__ static inline void apply(CELL& cell, T& rho, Vector<T, LatSet::d>& u) {
    u = cell.template get<VELOCITY<T, LatSet::d>>();
  }
};

template <typename CELLTYPE>
struct useFieldrhoU {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline void apply(CELL& cell, T& rho, Vector<T, LatSet::d>& u) {
    rho = cell.template get<GenericRho>();
    u = cell.template get<VELOCITY<T, LatSet::d>>();
  }
  // for compatible with force scheme
  __any__ static inline void apply(CELL& cell, const T force, T& rho, Vector<T, LatSet::d>& u) {
    apply(cell, rho, u);
  }
  // for compatible with force scheme
  __any__ static inline void apply(CELL& cell, const Vector<T, LatSet::d>& force, T& rho, Vector<T, LatSet::d>& u) {
    apply(cell, rho, u);
  }
};


// ----------------------------------------------------------------------------
// moment from const value

template <typename CELLTYPE, typename CONSTRHOTYPE = CONSTRHO<typename CELLTYPE::FloatType>, bool WriteToField = false>
struct constrho {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline T get(CELL& cell) {
    const T rho_value = cell.template get<CONSTRHOTYPE>();
    if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
    return rho_value;
  }
  __any__ static inline void apply(CELL& cell, T& rho_value) {
    rho_value = cell.template get<CONSTRHOTYPE>();
    if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
  }
  __any__ static inline void apply(CELL& cell) {
    static_assert(WriteToField, "constrho::apply(CELL& cell) must write to field");
    cell.template get<GenericRho>() = cell.template get<CONSTRHOTYPE>();
  }
};

template <typename CELLTYPE, bool WriteToField = false>
struct constU {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline Vector<T, LatSet::d> get(CELL& cell) {
    const Vector<T, LatSet::d> u_value = cell.template get<CONSTU<T, LatSet::d>>();
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
    return u_value;
  }
  __any__ static inline void apply(CELL& cell, Vector<T, LatSet::d>& u_value) {
    u_value = cell.template get<CONSTU<T, LatSet::d>>();
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
  }
  __any__ static inline void apply(CELL& cell) {
    static_assert(WriteToField, "constU::apply(CELL& cell) must write to field");
    cell.template get<VELOCITY<T, LatSet::d>>() = cell.template get<CONSTU<T, LatSet::d>>();
  }
};


// ----------------------------------------------------------------------------
// moment calculation

template <typename CELLTYPE, bool WriteToField>
struct rhoImpl {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline void apply(CELL& cell, T& rho_value) {
    rho_value = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
    if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
  }
};

template <typename CELLTYPE, bool WriteToField = false>
struct rho {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline T get(CELL& cell) {
    T rho_value{};
    apply(cell, rho_value);
    return rho_value;
  }
  __any__ static inline void apply(CELL& cell) {
    static_assert(WriteToField, "rho::apply(CELL& cell) must write to field");
    T& rho_value = cell.template get<GenericRho>();
    apply(cell, rho_value);
  }

  __any__ static inline void apply(CELL& cell, T& rho_value) {
    rhoImpl<CELL, WriteToField>::apply(cell, rho_value);
  }
  __any__ static inline void apply(CELL& cell, T& rho_value, Vector<T, LatSet::d>& u_value) {
    rhoImpl<CELL, WriteToField>::apply(cell, rho_value);
  }
};


template <typename CELLTYPE, typename SOURCE, bool WriteToField>
struct sourcerhoImpl {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline void apply(CELL& cell, const T source, T& rho_value) {
    rho_value = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
    // fOmega: avoid lattice artifact
    rho_value += source * T{0.5};
    if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
  }
};

// update rho(usually temperature or concentration in advection-diffusion problems) with
// source term, no need to preprocess the source term
template <typename CELLTYPE, typename SOURCE, bool WriteToField = false>
struct sourcerho {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline T get(CELL& cell) {
    const auto source = cell.template get<SOURCE>();
    return get(cell, source);
  }
  __any__ static inline T get(CELL& cell, const T source) {
    T rho_value{};
    apply(cell, source, rho_value);
    return rho_value;
  }
  __any__ static inline void apply(CELL& cell) {
    static_assert(WriteToField, "sourcerho::apply(CELL& cell) must write to field");
    const auto source = cell.template get<SOURCE>();
    apply(cell, source);
  }
  __any__ static inline void apply(CELL& cell, const T source) {
    static_assert(WriteToField, "sourcerho::apply(CELL& cell, T source) must write to field");
    T& rho_value = cell.template get<GenericRho>();
    apply(cell, source, rho_value);
  }

  __any__ static inline void apply(CELL& cell, const T source, T& rho_value) {
    sourcerhoImpl<CELL, SOURCE, WriteToField>::apply(cell, source, rho_value);
  }
};


template <typename CELLTYPE, bool WriteToField>
struct UImpl {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline void apply(CELL& cell, Vector<T, LatSet::d>& u_value) {
    u_value.clear();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += latset::c<LatSet>(i) * cell[i];
    }
    u_value /= rho_value;
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
  }
};

template <typename CELLTYPE, bool WriteToField = false>
struct U {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline Vector<T, LatSet::d> get(CELL& cell) {
    Vector<T, LatSet::d> u_value{};
    apply(cell, u_value);
    return u_value;
  }
  __any__ static inline void apply(CELL& cell) {
    static_assert(WriteToField, "U::apply(CELL& cell) must write to field");
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(cell, u_value);
  }

  __any__ static inline void apply(CELL& cell, Vector<T, LatSet::d>& u_value) {
    UImpl<CELL, WriteToField>::apply(cell, u_value);
  }
};


template <typename CELLTYPE, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forceUImpl {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

  __any__ static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha, Vector<T, LatSet::d>& u_value) {
    u_value.clear();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += latset::c<LatSet>(i) * cell[i];
    }
    u_value += f_alpha * T{0.5};
    u_value /= rho_value;
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
  }
  // for scalar force
  __any__ static inline void apply(CELL& cell, const T f, Vector<T, LatSet::d>& u_value) {
    u_value.clear();
    T rho_value{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += latset::c<LatSet>(i) * cell[i];
    }
    u_value[scalardir] += f * T{0.5};
    u_value /= rho_value;
    if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
  }
};

template <typename CELLTYPE, typename ForceScheme, bool WriteToField = false, unsigned int dir = 2>
struct forceU {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

  __any__ static inline Vector<T, LatSet::d> get(CELL& cell) {
    const auto force = ForceScheme::getForce(cell);
    return get(cell, force);
  }
  __any__ static inline Vector<T, LatSet::d> get(CELL& cell,
                                         const Vector<T, LatSet::d>& f_alpha) {
    Vector<T, LatSet::d> u_value;
    apply(cell, f_alpha, u_value);
    return u_value;
  }
  __any__ static inline Vector<T, LatSet::d> get(CELL& cell, const T f) {
    Vector<T, LatSet::d> u_value;
    apply(cell, f, u_value);
    return u_value;
  }
  __any__ static inline void apply(CELL& cell) {
    static_assert(WriteToField, "forceU::apply(CELL& cell) must write to field");
    const auto force = ForceScheme::getForce(cell);
    apply(cell, force);
  }
  __any__ static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha) {
    static_assert(WriteToField, "forceU::apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha) must write to field");
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(cell, f_alpha, u_value);
  }
  __any__ static inline void apply(CELL& cell, const T f) {
    static_assert(WriteToField, "forceU::apply(CELL& cell, T f) must write to field");
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(cell, f, u_value);
  }

  __any__ static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha, Vector<T, LatSet::d>& u_value) {
    forceUImpl<CELL, ForceScheme, WriteToField, dir>::apply(cell, f_alpha, u_value);
  }
  // for scalar force
  __any__ static inline void apply(CELL& cell, const T f, Vector<T, LatSet::d>& u_value) {
    forceUImpl<CELL, ForceScheme, WriteToField, dir>::apply(cell, f, u_value);
  }
};


template <typename CELLTYPE, bool WriteToField>
struct rhoUImpl {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell, T& rho_value, Vector<T, LatSet::d>& u_value) {
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += latset::c<LatSet>(i) * cell[i];
    }
    u_value /= rho_value;
    if constexpr (WriteToField) {
      cell.template get<GenericRho>() = rho_value;
      cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
    }
  }
};

template <typename CELLTYPE, bool WriteToField = false>
struct rhoU {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    static_assert(WriteToField, "rhoU::apply(CELL& cell) must write to field");
    T& rho_value = cell.template get<GenericRho>();
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(cell, rho_value, u_value);
  }

  __any__ static void apply(CELL& cell, T& rho_value, Vector<T, LatSet::d>& u_value) {
    rhoUImpl<CELL, WriteToField>::apply(cell, rho_value, u_value);
  }
};


template <typename CELLTYPE, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forcerhoUImpl {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

  __any__ static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value,
                           Vector<T, LatSet::d>& u_value) {
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += latset::c<LatSet>(i) * cell[i];
    }
    u_value += f_alpha * T{0.5};
    u_value /= rho_value;
    if constexpr (WriteToField) {
      cell.template get<GenericRho>() = rho_value;
      cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
    }
  }
  // for scalar force
  __any__ static inline void apply(CELL& cell, const T f, T& rho_value, Vector<T, LatSet::d>& u_value) {
    rho_value = T{};
    u_value.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho_value += cell[i];
      u_value += latset::c<LatSet>(i) * cell[i];
    }
    u_value[scalardir] += f * T{0.5};
    u_value /= rho_value;
    if constexpr (WriteToField) {
      cell.template get<GenericRho>() = rho_value;
      cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
    }
  }
};

template <typename CELLTYPE, typename ForceScheme, bool WriteToField = false, unsigned int dir = 2>
struct forcerhoU {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

  __any__ static inline void apply(CELL& cell) {
    static_assert(WriteToField, "forcerhoU::apply(CELL& cell) must write to field");
    const auto force = ForceScheme::getForce(cell);
    apply(cell, force);
  }
  __any__ static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha) {
    static_assert(WriteToField, "forcerhoU::apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha) must write to field");
    T& rho_value = cell.template get<GenericRho>();
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(cell, f_alpha, rho_value, u_value);
  }
  __any__ static inline void apply(CELL& cell, const T f) {
    static_assert(WriteToField, "forcerhoU::apply(CELL& cell, T f) must write to field");
    T& rho_value = cell.template get<GenericRho>();
    Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(cell, f, rho_value, u_value);
  }

  __any__ static inline void apply(CELL& cell, T& rho_value, Vector<T, LatSet::d>& u_value) {
    const auto force = ForceScheme::getForce(cell);
    apply(cell, force, rho_value, u_value);
  }

  __any__ static inline void apply(CELL& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value,
                           Vector<T, LatSet::d>& u_value) {
    forcerhoUImpl<CELL, ForceScheme, WriteToField, dir>::apply(cell, f_alpha, rho_value, u_value);
  }
  // for scalar force
  __any__ static inline void apply(CELL& cell, const T f, T& rho_value, Vector<T, LatSet::d>& u_value) {
    forcerhoUImpl<CELL, ForceScheme, WriteToField, dir>::apply(cell, f, rho_value, u_value);
  }
};


// second moment of non-equilibrium part of the distribution function
template <typename CELLTYPE>
struct Pi_ab_neq {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // f_i^(1) = f_i - f_i^eq
  // PI_ab^neq = SUM_i((f_i - f_i^eq)*c_ia*c_ib)
  // PI_ab^eq = SUM_i(f_i^eq*c_ia*c_ib) = rho*U_a*U_b + rho*Cs^2*delta_ab
  __any__ static inline void apply(CELL& cell, const T rho, const Vector<T, LatSet::d>& u, 
  std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor) {
    unsigned int i{};
    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        T value{};
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          value += latset::c<LatSet>(k)[alpha] * latset::c<LatSet>(k)[beta] * cell[k];
        }
        // remove the equilibrium part: PI_ab^eq
        value -= rho * u[alpha] * u[beta];
        if (alpha == beta) value -= rho * LatSet::cs2;
        tensor[i] = value;
        ++i;
      }
    }
  }
};

// second moment of non-equilibrium part of the distribution function with force
template <typename CELLTYPE>
struct forcePi_ab_neq {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // f_i^(1) = f_i - f_i^eq
  // PI_ab^neq = SUM_i((f_i - f_i^eq)*c_ia*c_ib)
  // PI_ab^eq = SUM_i(f_i^eq*c_ia*c_ib) = rho*U_a*U_b + rho*Cs^2*delta_ab
  __any__ static inline void apply(CELL& cell, const T rho, const Vector<T, LatSet::d>& u, const Vector<T, LatSet::d>& f_alpha,
                         std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor) {
    unsigned int i{};
    // remove force term in u
    Vector<T, LatSet::d> unew = u - f_alpha * T{0.5}; //(T{0.5} / rho);

    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        T value{};
        const T force = T{0.5} * (f_alpha[alpha] * unew[beta] + f_alpha[beta] * unew[alpha]);  //* rho
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          value += latset::c<LatSet>(k)[alpha] * latset::c<LatSet>(k)[beta] * cell[k];
        }
        value += force;
        // remove the equilibrium part: PI_ab^eq
        value -= rho * unew[alpha] * unew[beta];
        if (alpha == beta) value -= rho * LatSet::cs2;
        tensor[i] = value;
        ++i;
      }
    }
  }
};

// stress tensor
template <typename CELLTYPE>
struct stress {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // sigma_ab = -(1 - 1/(2*tau))*SUM_i(f_i^(1)*c_ia*c_ib)
  // f_i^(1) = f_i - f_i^eq
  // PI_ab^eq = SUM_i(f_i^eq*c_ia*c_ib) = rho*U_a*U_b + rho*Cs^2*delta_ab
  __any__ static inline void apply(CELL& cell, const T rho, const Vector<T, LatSet::d>& u, 
  std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& stress_tensor) {
    unsigned int i{};
    const T coeff = T{0.5} * cell.getOmega() - T{1};
    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        T value{};
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          value += latset::c<LatSet>(k)[alpha] * latset::c<LatSet>(k)[beta] * cell[k];
        }
        // remove the equilibrium part: PI_ab^eq
        value -= rho * u[alpha] * u[beta];
        if (alpha == beta) value -= rho * LatSet::cs2;
        // multiply by the coefficient
        value *= coeff;
        stress_tensor[i] = value;
        ++i;
      }
    }
  }
};

// strain rate tensor/ rate of deformation matrix
template <typename CELLTYPE>
struct strainRate {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // sigma_ab = -(1 - delta_T/(2*tau)) * SUM_i(f_i^(1)*c_ia*c_ib)
  // sigma_ab = 2 * mu * S_ab = 2 * rho * nu * S_ab
  // S_ab = simga_ab / (2 * rho * Cs^2 * (tau - 0.5))
  // S_ab = (1/rho)(-Cs^2/(2*tau))*SUM_i(f_i^(1)*c_ia*c_ib)
  // f_i^(1) = f_i - f_i^eq
  // PI_ab^eq = SUM_i(f_i^eq*c_ia*c_ib) = rho*U_a*U_b + rho*Cs^2*delta_ab
  __any__ static inline void apply(CELL& cell, const T rho, const Vector<T, LatSet::d>& u, 
  std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& strain_rate_tensor) {
    unsigned int i{};
    // use cell.template get<OMEGA<T>>() here instead of cell.getOmega()
    const T coeff = T{-1.5} * cell.template get<OMEGA<T>>() / cell.template get<typename CELL::GenericRho>();
    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        T value{};
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          value += latset::c<LatSet>(k)[alpha] * latset::c<LatSet>(k)[beta] * cell[k];
        }
        // remove the equilibrium part: PI_ab^eq
        value -= rho * u[alpha] * u[beta];
        if (alpha == beta) value -= rho * LatSet::cs2;
        // multiply by the coefficient
        value *= coeff;
        strain_rate_tensor[i] = value;
        ++i;
      }
    }
  }
};


template <typename CELLTYPE>
struct shearRateMagImpl {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& strain_rate_tensor) {
    T value{};
    unsigned int i{};
    for (unsigned int alpha = 0; alpha < LatSet::d; ++alpha) {
      for (unsigned int beta = alpha; beta < LatSet::d; ++beta) {
        T sq = strain_rate_tensor[i] * strain_rate_tensor[i];
        if (alpha != beta) sq *= T{2};
        value += sq;
        ++i;
      }
    }
    return std::sqrt(T{2} * value);
  }
};

// magnitude of shear rate
template <typename CELLTYPE, bool WriteToField = false>
struct shearRateMag {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    const T rho_value = cell.template get<GenericRho>();
    const Vector<T, LatSet::d>& u_value = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(cell, rho_value, u_value);
  }
  __any__ static void apply(CELL& cell, const T rho_value, const Vector<T, LatSet::d>& u_value) {
    static_assert(WriteToField, "shearRateMag::apply(CELL& cell, T& rho_value, Vector<T, LatSet::d>& u_value) must write to field");
    std::array<T, util::SymmetricMatrixSize<LatSet::d>()> strain_rate_tensor;
    strainRate<CELL>::apply(cell, rho_value, u_value, strain_rate_tensor);
    cell.template get<StrainRateMag<T>>() = get(strain_rate_tensor);
  }

  __any__ static inline T get(CELL& cell, const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& strain_rate_tensor) {
    const T gamma = get(strain_rate_tensor);
    if constexpr (WriteToField) cell.template get<StrainRateMag<T>>() = gamma;
    return gamma;
  }
  // \dot{gamma} = sqrt(2*trace(S^2)) = sqrt(2*SUM_{a,b}S_ab^2)
  __any__ static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& strain_rate_tensor) {
    return shearRateMagImpl<CELL>::get(strain_rate_tensor);
  }
};


// a tuple containing multiple momenta
template <typename... Momenta>
struct MomentaTuple {
  // get the first template parameter 
  using CELL = typename tmp::FirstParam<Momenta...>::type::CELL;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    (Momenta::apply(cell), ...);
  }
  __any__ static void apply(CELL& cell, T& rho_value, Vector<T, LatSet::d>& u_value) {
    (Momenta::apply(cell, rho_value, u_value), ...);
  }
};

// ----------------- legacy code -----------------

template <typename T, typename LatSet>
struct Rho {
  // return rho
  static T get(const BasicPopCell<T, LatSet>& cell) {
    T rho = T(0);
    for (unsigned int i = 0; i < LatSet::q; ++i) rho += cell[i];
    return rho;
  }

  // return rho with source
  static T get(const BasicPopCell<T, LatSet>& cell, const T source) {
    T rho = T(0);
    for (unsigned int i = 0; i < LatSet::q; ++i) rho += cell[i];
    rho += source * T{0.5};
    return rho;
  }

  // compute rho
  static void apply(const BasicPopCell<T, LatSet>& cell, T& rho) {
    rho = T(0);
    for (unsigned int i = 0; i < LatSet::q; ++i) rho += cell[i];
  }

  // compute rho with source: C = sum(f) + q/2
  static void apply(const BasicPopCell<T, LatSet>& cell, T& rho, const T source) {
    rho = T(0);
    for (unsigned int i = 0; i < LatSet::q; ++i) rho += cell[i];
    rho += source * T{0.5};
  }
};

template <typename T, typename LatSet>
struct Velocity {
  // get velocity
  static Vector<T, LatSet::d> get(const BasicPopCell<T, LatSet>& cell) {
    Vector<T, LatSet::d> u;
    T rho = T(0);
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + latset::c<LatSet>(i) * cell[i];
    }
    return u / rho;
  }

  // get velocity with force
  static Vector<T, LatSet::d> get(const BasicPopCell<T, LatSet>& cell,
                                  const std::array<T, LatSet::q>& fi) {
    Vector<T, LatSet::d> u;
    T rho = T(0);
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + (latset::c<LatSet>(i) * (cell[i] + T{0.5} * fi[i]));
    }
    return u / rho;
  }

  // compute velocity
  static void apply(const BasicPopCell<T, LatSet>& cell, Vector<T, LatSet::d>& u) {
    u.clear();
    T rho = T(0);
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + latset::c<LatSet>(i) * cell[i];
    }
    u = u / rho;
  }

  // compute velocity with force
  static void apply(const BasicPopCell<T, LatSet>& cell, Vector<T, LatSet::d>& u,
                    const std::array<T, LatSet::q>& fi) {
    u.clear();
    T rho = T(0);
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + (latset::c<LatSet>(i) * (cell[i] + T{0.5} * fi[i]));
    }
    u = u / rho;
  }
};

template <typename T, typename LatSet>
struct RhoVelocity {
  // compute rho and velocity
  static void apply(const BasicPopCell<T, LatSet>& cell, T& rho, Vector<T, LatSet::d>& u) {
    rho = T(0);
    u.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u = u + latset::c<LatSet>(i) * cell[i];
    }
    u = u / rho;
  }
};

}  // namespace moment