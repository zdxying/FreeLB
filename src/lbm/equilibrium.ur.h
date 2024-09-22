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

// equlibrium.h

#pragma once

#include "lbm/equilibrium.h"

#ifdef _UNROLLFOR

namespace equilibrium {

#ifdef __CUDA_ARCH__
template <typename T, typename LatSet, typename TypePack>
using CELL = cudev::Cell<T, LatSet, TypePack>;

#else
template <typename T, typename LatSet, typename TypePack>
using CELL = Cell<T, LatSet, TypePack>;
#endif

template <typename T, typename TypePack>
struct SecondOrder<CELL<T, D2Q9<T> ,TypePack>> {
  using LatSet = typename CELL<T, D2Q9<T> ,TypePack>::LatticeSet;
	using CELLTYPE = CELL<T, D2Q9<T> ,TypePack>;


  __any__ static inline T get(int k, const Vector<T, LatSet::d> &u, T rho, T u2) {
    const T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 - LatSet::InvCs2 * u2 * T{0.5});
  }

  __any__ static void apply(std::array<T, LatSet::q> &feq, T rho, const Vector<T, LatSet::d> &u) {
		constexpr T InvCs4_ = T{0.5} * LatSet::InvCs4;
    const T u0 = u[0];
		const T u1 = u[1];
		
		const T u0p1 = u0 + u1;
		const T u0m1 = u0 - u1;
		const T InvCs4_u0p1_2  = InvCs4_ * u0p1 * u0p1;
		const T InvCs4_u0m1_2 = InvCs4_ * u0m1 * u0m1;

		const T u0_2 = u0 * u0;
		const T u1_2 = u1 * u1;

		const T _InvCs2_u2_ = T{1} - LatSet::InvCs2 * (u0_2 + u1_2) * T{0.5};

		const T InvCs2u0 = LatSet::InvCs2 * u0;
		const T InvCs2u1 = LatSet::InvCs2 * u1;
		const T InvCs2u0p1 = LatSet::InvCs2 * u0p1;
		const T InvCs2u0m1 = LatSet::InvCs2 * u0m1;
		const T InvCs4_u0_2 = InvCs4_ * u0_2;
		const T InvCs4_u1_2 = InvCs4_ * u1_2;

		feq[0] = latset::w<LatSet>(0) * rho * _InvCs2_u2_;
		feq[1] = latset::w<LatSet>(1) * rho * (_InvCs2_u2_ + InvCs2u0   + InvCs4_u0_2);
		feq[2] = latset::w<LatSet>(2) * rho * (_InvCs2_u2_ - InvCs2u0   + InvCs4_u0_2);
		feq[3] = latset::w<LatSet>(3) * rho * (_InvCs2_u2_ + InvCs2u1   + InvCs4_u1_2);
		feq[4] = latset::w<LatSet>(4) * rho * (_InvCs2_u2_ - InvCs2u1   + InvCs4_u1_2);
		feq[5] = latset::w<LatSet>(5) * rho * (_InvCs2_u2_ + InvCs2u0p1 + InvCs4_u0p1_2);
		feq[6] = latset::w<LatSet>(6) * rho * (_InvCs2_u2_ - InvCs2u0p1 + InvCs4_u0p1_2);
		feq[7] = latset::w<LatSet>(7) * rho * (_InvCs2_u2_ + InvCs2u0m1 + InvCs4_u0m1_2);
		feq[8] = latset::w<LatSet>(8) * rho * (_InvCs2_u2_ - InvCs2u0m1 + InvCs4_u0m1_2);
  }

	using GenericRho = typename CELLTYPE::GenericRho;
	__any__ static void apply(std::array<T, LatSet::q> &feq, const CELLTYPE &cell) {
    const T rho = cell.template get<GenericRho>();
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
		apply(feq, rho, u);
	}
};

template <typename T, typename TypePack>
struct SecondOrder<CELL<T, D3Q19<T> ,TypePack>> {
  using LatSet = typename CELL<T, D3Q19<T> ,TypePack>::LatticeSet;
	using CELLTYPE = CELL<T, D3Q19<T> ,TypePack>;


  __any__ static inline T get(int k, const Vector<T, LatSet::d> &u, T rho, T u2) {
    const T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 - LatSet::InvCs2 * u2 * T{0.5});
  }

  __any__ static void apply(std::array<T, LatSet::q> &feq, T rho, const Vector<T, LatSet::d> &u) {
		constexpr T InvCs4_ = T{0.5} * LatSet::InvCs4;
    const T u0 = u[0];
		const T u1 = u[1];
		const T u2 = u[2];
		
		const T u0p1 = u0 + u1;
		const T u0p2 = u0 + u2;
		const T u1p2 = u1 + u2;
		const T u0m1 = u0 - u1;
		const T u0m2 = u0 - u2;
		const T u1m2 = u1 - u2;

		const T u0_2 = u0 * u0;
		const T u1_2 = u1 * u1;
		const T u2_2 = u2 * u2;

		const T _InvCs2_u2_ = T{1} - LatSet::InvCs2 * (u0_2 + u1_2 + u2_2) * T{0.5};

		const T InvCs2u0 = LatSet::InvCs2 * u0;
		const T InvCs2u1 = LatSet::InvCs2 * u1;
		const T InvCs2u2 = LatSet::InvCs2 * u2;

		const T InvCs2u0p1 = LatSet::InvCs2 * u0p1;
		const T InvCs2u0p2 = LatSet::InvCs2 * u0p2;
		const T InvCs2u1p2 = LatSet::InvCs2 * u1p2;
		const T InvCs2u0m1 = LatSet::InvCs2 * u0m1;
		const T InvCs2u0m2 = LatSet::InvCs2 * u0m2;
		const T InvCs2u1m2 = LatSet::InvCs2 * u1m2;

		const T InvCs4_u0_2   = InvCs4_ * u0_2;
		const T InvCs4_u1_2   = InvCs4_ * u1_2;
		const T InvCs4_u2_2   = InvCs4_ * u2_2;
		const T InvCs4_u0p1_2 = InvCs4_ * u0p1 * u0p1;
		const T InvCs4_u0p2_2 = InvCs4_ * u0p2 * u0p2;
		const T InvCs4_u1p2_2 = InvCs4_ * u1p2 * u1p2;
		const T InvCs4_u0m1_2 = InvCs4_ * u0m1 * u0m1;
		const T InvCs4_u0_2_2 = InvCs4_ * u0m2 * u0m2;
		const T InvCs4_u1_2_2 = InvCs4_ * u1m2 * u1m2;

		feq[0] = latset::w<LatSet>(0) * rho * _InvCs2_u2_;

		feq[1] = latset::w<LatSet>(1) * rho * (_InvCs2_u2_ + InvCs2u0   + InvCs4_u0_2);
		feq[2] = latset::w<LatSet>(2) * rho * (_InvCs2_u2_ - InvCs2u0   + InvCs4_u0_2);
		feq[3] = latset::w<LatSet>(3) * rho * (_InvCs2_u2_ + InvCs2u1   + InvCs4_u1_2);
		feq[4] = latset::w<LatSet>(4) * rho * (_InvCs2_u2_ - InvCs2u1   + InvCs4_u1_2);
		feq[5] = latset::w<LatSet>(5) * rho * (_InvCs2_u2_ + InvCs2u2   + InvCs4_u2_2);
		feq[6] = latset::w<LatSet>(6) * rho * (_InvCs2_u2_ - InvCs2u2   + InvCs4_u2_2);

		feq[7]  = latset::w<LatSet>(7)  * rho * (_InvCs2_u2_ + InvCs2u0p1  + InvCs4_u0p1_2);
		feq[8]  = latset::w<LatSet>(8)  * rho * (_InvCs2_u2_ - InvCs2u0p1  + InvCs4_u0p1_2);
		feq[9]  = latset::w<LatSet>(9)  * rho * (_InvCs2_u2_ + InvCs2u0p2  + InvCs4_u0p2_2);
		feq[10] = latset::w<LatSet>(10) * rho * (_InvCs2_u2_ - InvCs2u0p2  + InvCs4_u0p2_2);
		feq[11] = latset::w<LatSet>(11) * rho * (_InvCs2_u2_ + InvCs2u1p2  + InvCs4_u1p2_2);
		feq[12] = latset::w<LatSet>(12) * rho * (_InvCs2_u2_ - InvCs2u1p2  + InvCs4_u1p2_2);

		feq[13] = latset::w<LatSet>(13) * rho * (_InvCs2_u2_ + InvCs2u0m1 + InvCs4_u0m1_2);
		feq[14] = latset::w<LatSet>(14) * rho * (_InvCs2_u2_ - InvCs2u0m1 + InvCs4_u0m1_2);
		feq[15] = latset::w<LatSet>(15) * rho * (_InvCs2_u2_ + InvCs2u0m2 + InvCs4_u0_2_2);
		feq[16] = latset::w<LatSet>(16) * rho * (_InvCs2_u2_ - InvCs2u0m2 + InvCs4_u0_2_2);
		feq[17] = latset::w<LatSet>(17) * rho * (_InvCs2_u2_ + InvCs2u1m2 + InvCs4_u1_2_2);
		feq[18] = latset::w<LatSet>(18) * rho * (_InvCs2_u2_ - InvCs2u1m2 + InvCs4_u1_2_2);
  }

	using GenericRho = typename CELLTYPE::GenericRho;
	__any__ static void apply(std::array<T, LatSet::q> &feq, const CELLTYPE &cell) {
    const T rho = cell.template get<GenericRho>();
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
		apply(feq, rho, u);
	}
};

}  // namespace equilibrium

#endif