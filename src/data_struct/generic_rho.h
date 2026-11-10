/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 *
 */

// generic_rho.h
//
// GenericRho selection: find the RHO/TEMP/CONC field alias inside a field
// TypePack.  Lives in its own header because two independent lattice headers
// need it and neither should have to be included before the other:
//   * data_struct/block_lattice_base.h  (host BlockLatticeBase, global set)
//   * data_struct/cuda_block_lattice.h  (cudev::BlockLatticeBase, cudev set)
//
// Kept out of utils/tmp.h: these need the field aliases from utils/alias.h,
// which the generic TMP header should not depend on.  Both sets are needed --
// the CUDA cell lives in namespace cudev and names the cudev version.

#pragma once

#include "utils/alias.h"  // RHO / TEMP / CONC
#include "utils/tmp.h"    // TypePack, isOneOf, isTypeInTuple

template <typename T, typename Pack>
struct GetGenericRhoType {
  using type = std::conditional_t<
    isTypeInTuple<RHO<T>, Pack>::value, RHO<T>,
    std::conditional_t<
      isTypeInTuple<TEMP<T>, Pack>::value, TEMP<T>,
      std::conditional_t<isTypeInTuple<CONC<T>, Pack>::value, CONC<T>, void>>>;
};

template <typename T, typename First>
struct FindGenericRhoType;

template <typename T, typename First, typename... Rest>
struct FindGenericRhoType<T, TypePack<First, Rest...>> {
  using type =
    std::conditional_t<isOneOf<First, RHO<T>, TEMP<T>, CONC<T>>::value, First,
                       typename FindGenericRhoType<T, TypePack<Rest...>>::type>;
};

template <typename T, typename First>
struct FindGenericRhoType<T, TypePack<First>> {
  using type =
    std::conditional_t<isOneOf<First, RHO<T>, TEMP<T>, CONC<T>>::value, First, void>;
};

namespace cudev {

#ifdef __CUDACC__

template <typename T, typename Pack>
struct FindGenericRhoType;

template <typename T, typename First, typename... Rest>
struct FindGenericRhoType<T, TypePack<First, Rest...>> {
  using type =
    std::conditional_t<isOneOf<First, RHO<T>, TEMP<T>, CONC<T>>::value, First,
                       typename FindGenericRhoType<T, TypePack<Rest...>>::type>;
};

template <typename T, typename First>
struct FindGenericRhoType<T, TypePack<First>> {
  using type =
    std::conditional_t<isOneOf<First, RHO<T>, TEMP<T>, CONC<T>>::value, First, void>;
};

#endif

}  // namespace cudev
