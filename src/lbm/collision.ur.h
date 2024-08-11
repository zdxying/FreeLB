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

// collision.ur.h

// unroll for loop version of collision.h with template specialization

#pragma once

#include "lbm/collision.h"

#ifdef _UNROLLFOR

namespace collision {

#ifdef __CUDA_ARCH__
template <typename T, typename LatSet, typename TypePack>
using CELL = cudev::Cell<T, LatSet, TypePack>;

#else
template <typename T, typename LatSet, typename TypePack>
using CELL = Cell<T, LatSet, TypePack>;
#endif

// full way bounce back, could be regarded as a mpdified collision process
// swap the populations in the opposite direction
template <typename T, typename TypePack>
struct BounceBack<CELL<T, D2Q9<T> ,TypePack>> {
  __any__ static void apply(CELL<T, D2Q9<T> ,TypePack>& cell) {
		const T celltemp1 = cell[1];
		cell[1] = cell[2];
		cell[2] = celltemp1;
		const T celltemp3 = cell[3];
		cell[3] = cell[4];
		cell[4] = celltemp3;
		const T celltemp5 = cell[5];
		cell[5] = cell[6];
		cell[6] = celltemp5;
		const T celltemp7 = cell[7];
		cell[7] = cell[8];
		cell[8] = celltemp7;
  }
};

template <typename T, typename TypePack>
struct BounceBack<CELL<T, D3Q19<T> ,TypePack>> {
  __any__ static void apply(CELL<T, D3Q19<T> ,TypePack>& cell) {
		const T celltemp1 = cell[1];
		cell[1] = cell[2];
		cell[2] = celltemp1;
		const T celltemp3 = cell[3];
		cell[3] = cell[4];
		cell[4] = celltemp3;
		const T celltemp5 = cell[5];
		cell[5] = cell[6];
		cell[6] = celltemp5;
		const T celltemp7 = cell[7];
		cell[7] = cell[8];
		cell[8] = celltemp7;
		const T celltemp9 = cell[9];
		cell[9] = cell[10];
		cell[10] = celltemp9;
		const T celltemp11 = cell[11];
		cell[11] = cell[12];
		cell[12] = celltemp11;
		const T celltemp13 = cell[13];
		cell[13] = cell[14];
		cell[14] = celltemp13;
		const T celltemp15 = cell[15];
		cell[15] = cell[16];
		cell[16] = celltemp15;
		const T celltemp17 = cell[17];
		cell[17] = cell[18];
		cell[18] = celltemp17;
  }
};

// full way bounce back with moving wall, could be regarded as a modified collision
// process swap the populations in the opposite direction
template <typename T, typename TypePack>
struct BounceBackMovingWall<CELL<T, D2Q9<T> ,TypePack>> {
  using LatSet = typename CELL<T, D2Q9<T> ,TypePack>::LatticeSet;
  using GenericRho = typename CELL<T, D2Q9<T> ,TypePack>::GenericRho;

  __any__ static void apply(CELL<T, D2Q9<T> ,TypePack>& cell) {
		constexpr T InvCs2_2 = 2 * LatSet::InvCs2;
    const T rhox = InvCs2_2 * cell.template get<GenericRho>();
		const T u0 = cell.template get<VELOCITY<T, 2>>()[0];
		const T u1 = cell.template get<VELOCITY<T, 2>>()[1];

		const T celltemp1 = cell[1];
		const T uc1 = u0 * latset::w<LatSet>(1) * rhox;
		cell[1] = cell[2] + uc1;
		cell[2] = celltemp1 - uc1;

		const T celltemp3 = cell[3];
		const T uc3 = u1 * latset::w<LatSet>(3) * rhox;
		cell[3] = cell[4] + uc3;
		cell[4] = celltemp3 - uc3;

		const T celltemp5 = cell[5];
		const T uc5 = (u0 + u1) * latset::w<LatSet>(5) * rhox;
		cell[5] = cell[6] + uc5;
		cell[6] = celltemp5 - uc5;

		const T celltemp7 = cell[7];
		const T uc7 = (u0 - u1) * latset::w<LatSet>(7) * rhox;
		cell[7] = cell[8] + uc7;
		cell[8] = celltemp7 - uc7;
  }
};

template <typename T, typename TypePack>
struct BounceBackMovingWall<CELL<T, D3Q19<T> ,TypePack>> {
  using LatSet = typename CELL<T, D3Q19<T> ,TypePack>::LatticeSet;
  using GenericRho = typename CELL<T, D3Q19<T> ,TypePack>::GenericRho;

  __any__ static void apply(CELL<T, D3Q19<T> ,TypePack>& cell) {
		constexpr T InvCs2_2 = 2 * LatSet::InvCs2;
    const T rhox = InvCs2_2 * cell.template get<GenericRho>();
		const T u0 = cell.template get<VELOCITY<T, 3>>()[0];
		const T u1 = cell.template get<VELOCITY<T, 3>>()[1];
		const T u2 = cell.template get<VELOCITY<T, 3>>()[2];

		const T celltemp1 = cell[1];
		const T uc1 = u0 * latset::w<LatSet>(1) * rhox;
		cell[1] = cell[2] + uc1;
		cell[2] = celltemp1 - uc1;

		const T celltemp3 = cell[3];
		const T uc3 = u1 * latset::w<LatSet>(3) * rhox;
		cell[3] = cell[4] + uc3;
		cell[4] = celltemp3 - uc3;

		const T celltemp5 = cell[5];
		const T uc5 = u2 * latset::w<LatSet>(5) * rhox;
		cell[5] = cell[6] + uc5;
		cell[6] = celltemp5 - uc5;

		const T celltemp7 = cell[7];
		const T uc7 = (u0 + u1) * latset::w<LatSet>(7) * rhox;
		cell[7] = cell[8] + uc7;
		cell[8] = celltemp7 - uc7;

		const T celltemp9 = cell[9];
		const T uc9 = (u0 + u2) * latset::w<LatSet>(9) * rhox;
		cell[9] = cell[10] + uc9;
		cell[10] = celltemp9 - uc9;

		const T celltemp11 = cell[11];
		const T uc11 = (u1 + u2) * latset::w<LatSet>(11) * rhox;
		cell[11] = cell[12] + uc11;
		cell[12] = celltemp11 - uc11;

		const T celltemp13 = cell[13];
		const T uc13 = (u0 - u1) * latset::w<LatSet>(13) * rhox;
		cell[13] = cell[14] + uc13;
		cell[14] = celltemp13 - uc13;

		const T celltemp15 = cell[15];
		const T uc15 = (u0 - u2) * latset::w<LatSet>(15) * rhox;
		cell[15] = cell[16] + uc15;
		cell[16] = celltemp15 - uc15;

		const T celltemp17 = cell[17];
		const T uc17 = (u1 - u2) * latset::w<LatSet>(17) * rhox;
		cell[17] = cell[18] + uc17;
		cell[18] = celltemp17 - uc17;
  }
};


}  // namespace collision

#endif