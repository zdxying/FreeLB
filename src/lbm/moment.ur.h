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

// moment.ur.h

// unroll for loop version of moment.h with template specialization

#pragma once

#include "lbm/moment.h"

#ifdef _UNROLLFOR

namespace moment {

#ifdef __CUDA_ARCH__
template <typename T, typename LatSet, typename TypePack>
using CELL = cudev::Cell<T, LatSet, TypePack>;

#else
template <typename T, typename LatSet, typename TypePack>
using CELL = Cell<T, LatSet, TypePack>;
#endif

template <typename T, typename TypePack, bool WriteToField>
struct rho<CELL<T, D2Q9<T> ,TypePack>, WriteToField> {
  using GenericRho = typename CELL<T, D2Q9<T> ,TypePack>::GenericRho;

  __any__ static inline T get(CELL<T, D2Q9<T> ,TypePack>& cell) {
    if constexpr (WriteToField) {
			T rho_value{};
			rho_value = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
			cell.template get<GenericRho>() = rho_value;
			return rho_value;
		} else {
			return cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
		}
  }
  __any__ static inline void apply(CELL<T, D2Q9<T> ,TypePack>& cell, T& rho_value) {
		rho_value = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
    if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
  }
  // always write to field
  __any__ static inline void apply(CELL<T, D2Q9<T> ,TypePack>& cell) {
    T& rho_value = cell.template get<GenericRho>();
    rho_value = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
  }
};

template <typename T, typename TypePack, bool WriteToField>
struct rho<CELL<T, D3Q19<T> ,TypePack>, WriteToField> {
  using GenericRho = typename CELL<T, D3Q19<T> ,TypePack>::GenericRho;

  __any__ static inline T get(CELL<T, D3Q19<T> ,TypePack>& cell) {
    if constexpr (WriteToField) {
			T rho_value{};
			rho_value = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] 
			+ cell[9] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18];
			cell.template get<GenericRho>() = rho_value;
			return rho_value;
		} else {
			return cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8]
			+ cell[9] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18];
		}
  }
  __any__ static inline void apply(CELL<T, D3Q19<T> ,TypePack>& cell, T& rho_value) {
		rho_value = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8]
		+ cell[9] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18];
    if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
  }
  // always write to field
  __any__ static inline void apply(CELL<T, D3Q19<T> ,TypePack>& cell) {
    T& rho_value = cell.template get<GenericRho>();
    rho_value = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8]
		+ cell[9] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18];
  }
};


template <typename T, typename TypePack, bool WriteToField>
struct rhou<CELL<T, D2Q9<T> ,TypePack>, WriteToField> {
  using GenericRho = typename CELL<T, D2Q9<T> ,TypePack>::GenericRho;

  __any__ static void apply(CELL<T, D2Q9<T> ,TypePack>& cell, T& rho_value, Vector<T, 2>& u_value) {
    rho_value = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
    u_value[0] = (cell[1] - cell[2] + cell[5] - cell[6] + cell[7] - cell[8]) / rho_value;
		u_value[1] = (cell[3] - cell[4] + cell[5] - cell[6] - cell[7] + cell[8]) / rho_value;
    if constexpr (WriteToField) {
      cell.template get<GenericRho>() = rho_value;
      cell.template get<VELOCITY<T, 2>>() = u_value;
    }
  }
  // will write to field regardless of WriteToField value
  __any__ static void apply(CELL<T, D2Q9<T> ,TypePack>& cell) {
    T& rho_value = cell.template get<GenericRho>();
    Vector<T, 2>& u_value = cell.template get<VELOCITY<T, 2>>();
		rho_value = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
    u_value[0] = (cell[1] - cell[2] + cell[5] - cell[6] + cell[7] - cell[8]) / rho_value;
		u_value[1] = (cell[3] - cell[4] + cell[5] - cell[6] - cell[7] + cell[8]) / rho_value;
  }
};

template <typename T, typename TypePack, bool WriteToField>
struct rhou<CELL<T, D3Q19<T> ,TypePack>, WriteToField> {
  using GenericRho = typename CELL<T, D3Q19<T> ,TypePack>::GenericRho;

  __any__ static void apply(CELL<T, D3Q19<T> ,TypePack>& cell, T& rho_value, Vector<T, 3>& u_value) {
    rho_value = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18];
    u_value[0] = (cell[1] - cell[2] + cell[7] - cell[8] + cell[9] - cell[10] + cell[13] - cell[14] + cell[15] - cell[16]) / rho_value;
		u_value[1] = (cell[3] - cell[4] + cell[7] - cell[8] + cell[11] - cell[12] - cell[13] + cell[14] + cell[17] - cell[18]) / rho_value;
		u_value[1] = (cell[5] - cell[6] + cell[9] - cell[10] + cell[11] - cell[12] - cell[15] + cell[16] - cell[17] + cell[18]) / rho_value;
    if constexpr (WriteToField) {
      cell.template get<GenericRho>() = rho_value;
      cell.template get<VELOCITY<T, 3>>() = u_value;
    }
  }
  // will write to field regardless of WriteToField value
  __any__ static void apply(CELL<T, D3Q19<T> ,TypePack>& cell) {
    T& rho_value = cell.template get<GenericRho>();
    Vector<T, 3>& u_value = cell.template get<VELOCITY<T, 3>>();
		rho_value = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18];
    u_value[0] = (cell[1] - cell[2] + cell[7] - cell[8] + cell[9] - cell[10] + cell[13] - cell[14] + cell[15] - cell[16]) / rho_value;
		u_value[1] = (cell[3] - cell[4] + cell[7] - cell[8] + cell[11] - cell[12] - cell[13] + cell[14] + cell[17] - cell[18]) / rho_value;
		u_value[1] = (cell[5] - cell[6] + cell[9] - cell[10] + cell[11] - cell[12] - cell[15] + cell[16] - cell[17] + cell[18]) / rho_value;
  }
};

}  // namespace moment

#endif