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


//------------------------------------

template <typename T, typename TypePack>
struct SecondOrderImpl<CELL<T, D2Q5<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;

__any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u){
constexpr T InvCs2x = T{0.5} * LatSet::InvCs2;
const T var0 = T{1} - (u[0]*u[0] + u[1]*u[1]) * InvCs2x;
const T rhowk0 = latset::w<LatSet>(0) * rho;
const T rhowk1 = latset::w<LatSet>(1) * rho;
const T InvCs2uck1 = LatSet::InvCs2 * (u[0]);
const T var0InvCs4uc2k1 = var0 + InvCs2uck1*InvCs2uck1*T{0.5};
const T InvCs2uck3 = LatSet::InvCs2 * (u[1]);
const T var0InvCs4uc2k3 = var0 + InvCs2uck3*InvCs2uck3*T{0.5};
feq[0] = rhowk0 * var0;
feq[1] = rhowk1 * (var0InvCs4uc2k1 + InvCs2uck1);
feq[2] = rhowk1 * (var0InvCs4uc2k1 - InvCs2uck1);
feq[3] = rhowk1 * (var0InvCs4uc2k3 + InvCs2uck3);
feq[4] = rhowk1 * (var0InvCs4uc2k3 - InvCs2uck3);
}
};

template <typename T, typename TypePack>
struct SecondOrderImpl<CELL<T, D2Q9<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;

__any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u){
constexpr T InvCs2x = T{0.5} * LatSet::InvCs2;
const T var0 = T{1} - (u[0]*u[0] + u[1]*u[1]) * InvCs2x;
const T rhowk0 = latset::w<LatSet>(0) * rho;
const T rhowk1 = latset::w<LatSet>(1) * rho;
const T rhowk5 = latset::w<LatSet>(5) * rho;
const T InvCs2uck1 = LatSet::InvCs2 * (u[0]);
const T var0InvCs4uc2k1 = var0 + InvCs2uck1*InvCs2uck1*T{0.5};
const T InvCs2uck3 = LatSet::InvCs2 * (u[1]);
const T var0InvCs4uc2k3 = var0 + InvCs2uck3*InvCs2uck3*T{0.5};
const T InvCs2uck5 = LatSet::InvCs2 * (u[0]+u[1]);
const T var0InvCs4uc2k5 = var0 + InvCs2uck5*InvCs2uck5*T{0.5};
const T InvCs2uck7 = LatSet::InvCs2 * (u[0]-u[1]);
const T var0InvCs4uc2k7 = var0 + InvCs2uck7*InvCs2uck7*T{0.5};
feq[0] = rhowk0 * var0;
feq[1] = rhowk1 * (var0InvCs4uc2k1 + InvCs2uck1);
feq[2] = rhowk1 * (var0InvCs4uc2k1 - InvCs2uck1);
feq[3] = rhowk1 * (var0InvCs4uc2k3 + InvCs2uck3);
feq[4] = rhowk1 * (var0InvCs4uc2k3 - InvCs2uck3);
feq[5] = rhowk5 * (var0InvCs4uc2k5 + InvCs2uck5);
feq[6] = rhowk5 * (var0InvCs4uc2k5 - InvCs2uck5);
feq[7] = rhowk5 * (var0InvCs4uc2k7 + InvCs2uck7);
feq[8] = rhowk5 * (var0InvCs4uc2k7 - InvCs2uck7);
}
};

template <typename T, typename TypePack>
struct SecondOrderImpl<CELL<T, D3Q7<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;

__any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u){
constexpr T InvCs2x = T{0.5} * LatSet::InvCs2;
const T var0 = T{1} - (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) * InvCs2x;
const T rhowk0 = latset::w<LatSet>(0) * rho;
const T rhowk1 = latset::w<LatSet>(1) * rho;
const T InvCs2uck1 = LatSet::InvCs2 * (u[0]);
const T var0InvCs4uc2k1 = var0 + InvCs2uck1*InvCs2uck1*T{0.5};
const T InvCs2uck3 = LatSet::InvCs2 * (u[1]);
const T var0InvCs4uc2k3 = var0 + InvCs2uck3*InvCs2uck3*T{0.5};
const T InvCs2uck5 = LatSet::InvCs2 * (u[2]);
const T var0InvCs4uc2k5 = var0 + InvCs2uck5*InvCs2uck5*T{0.5};
feq[0] = rhowk0 * var0;
feq[1] = rhowk1 * (var0InvCs4uc2k1 + InvCs2uck1);
feq[2] = rhowk1 * (var0InvCs4uc2k1 - InvCs2uck1);
feq[3] = rhowk1 * (var0InvCs4uc2k3 + InvCs2uck3);
feq[4] = rhowk1 * (var0InvCs4uc2k3 - InvCs2uck3);
feq[5] = rhowk1 * (var0InvCs4uc2k5 + InvCs2uck5);
feq[6] = rhowk1 * (var0InvCs4uc2k5 - InvCs2uck5);
}
};

template <typename T, typename TypePack>
struct SecondOrderImpl<CELL<T, D3Q15<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;

__any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u){
constexpr T InvCs2x = T{0.5} * LatSet::InvCs2;
const T var0 = T{1} - (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) * InvCs2x;
const T rhowk0 = latset::w<LatSet>(0) * rho;
const T rhowk1 = latset::w<LatSet>(1) * rho;
const T rhowk7 = latset::w<LatSet>(7) * rho;
const T InvCs2uck1 = LatSet::InvCs2 * (u[0]);
const T var0InvCs4uc2k1 = var0 + InvCs2uck1*InvCs2uck1*T{0.5};
const T InvCs2uck3 = LatSet::InvCs2 * (u[1]);
const T var0InvCs4uc2k3 = var0 + InvCs2uck3*InvCs2uck3*T{0.5};
const T InvCs2uck5 = LatSet::InvCs2 * (u[2]);
const T var0InvCs4uc2k5 = var0 + InvCs2uck5*InvCs2uck5*T{0.5};
const T InvCs2uck7 = LatSet::InvCs2 * (u[0]+u[1]+u[2]);
const T var0InvCs4uc2k7 = var0 + InvCs2uck7*InvCs2uck7*T{0.5};
const T InvCs2uck9 = LatSet::InvCs2 * (u[0]+u[1]-u[2]);
const T var0InvCs4uc2k9 = var0 + InvCs2uck9*InvCs2uck9*T{0.5};
const T InvCs2uck11 = LatSet::InvCs2 * (u[0]-u[1]+u[2]);
const T var0InvCs4uc2k11 = var0 + InvCs2uck11*InvCs2uck11*T{0.5};
const T InvCs2uck13 = LatSet::InvCs2 * (-u[0]+u[1]+u[2]);
const T var0InvCs4uc2k13 = var0 + InvCs2uck13*InvCs2uck13*T{0.5};
feq[0] = rhowk0 * var0;
feq[1] = rhowk1 * (var0InvCs4uc2k1 + InvCs2uck1);
feq[2] = rhowk1 * (var0InvCs4uc2k1 - InvCs2uck1);
feq[3] = rhowk1 * (var0InvCs4uc2k3 + InvCs2uck3);
feq[4] = rhowk1 * (var0InvCs4uc2k3 - InvCs2uck3);
feq[5] = rhowk1 * (var0InvCs4uc2k5 + InvCs2uck5);
feq[6] = rhowk1 * (var0InvCs4uc2k5 - InvCs2uck5);
feq[7] = rhowk7 * (var0InvCs4uc2k7 + InvCs2uck7);
feq[8] = rhowk7 * (var0InvCs4uc2k7 - InvCs2uck7);
feq[9] = rhowk7 * (var0InvCs4uc2k9 + InvCs2uck9);
feq[10] = rhowk7 * (var0InvCs4uc2k9 - InvCs2uck9);
feq[11] = rhowk7 * (var0InvCs4uc2k11 + InvCs2uck11);
feq[12] = rhowk7 * (var0InvCs4uc2k11 - InvCs2uck11);
feq[13] = rhowk7 * (var0InvCs4uc2k13 + InvCs2uck13);
feq[14] = rhowk7 * (var0InvCs4uc2k13 - InvCs2uck13);
}
};

template <typename T, typename TypePack>
struct SecondOrderImpl<CELL<T, D3Q19<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;

__any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u){
constexpr T InvCs2x = T{0.5} * LatSet::InvCs2;
const T var0 = T{1} - (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) * InvCs2x;
const T rhowk0 = latset::w<LatSet>(0) * rho;
const T rhowk1 = latset::w<LatSet>(1) * rho;
const T rhowk7 = latset::w<LatSet>(7) * rho;
const T InvCs2uck1 = LatSet::InvCs2 * (u[0]);
const T var0InvCs4uc2k1 = var0 + InvCs2uck1*InvCs2uck1*T{0.5};
const T InvCs2uck3 = LatSet::InvCs2 * (u[1]);
const T var0InvCs4uc2k3 = var0 + InvCs2uck3*InvCs2uck3*T{0.5};
const T InvCs2uck5 = LatSet::InvCs2 * (u[2]);
const T var0InvCs4uc2k5 = var0 + InvCs2uck5*InvCs2uck5*T{0.5};
const T InvCs2uck7 = LatSet::InvCs2 * (u[0]+u[1]);
const T var0InvCs4uc2k7 = var0 + InvCs2uck7*InvCs2uck7*T{0.5};
const T InvCs2uck9 = LatSet::InvCs2 * (u[0]+u[2]);
const T var0InvCs4uc2k9 = var0 + InvCs2uck9*InvCs2uck9*T{0.5};
const T InvCs2uck11 = LatSet::InvCs2 * (u[1]+u[2]);
const T var0InvCs4uc2k11 = var0 + InvCs2uck11*InvCs2uck11*T{0.5};
const T InvCs2uck13 = LatSet::InvCs2 * (u[0]-u[1]);
const T var0InvCs4uc2k13 = var0 + InvCs2uck13*InvCs2uck13*T{0.5};
const T InvCs2uck15 = LatSet::InvCs2 * (u[0]-u[2]);
const T var0InvCs4uc2k15 = var0 + InvCs2uck15*InvCs2uck15*T{0.5};
const T InvCs2uck17 = LatSet::InvCs2 * (u[1]-u[2]);
const T var0InvCs4uc2k17 = var0 + InvCs2uck17*InvCs2uck17*T{0.5};
feq[0] = rhowk0 * var0;
feq[1] = rhowk1 * (var0InvCs4uc2k1 + InvCs2uck1);
feq[2] = rhowk1 * (var0InvCs4uc2k1 - InvCs2uck1);
feq[3] = rhowk1 * (var0InvCs4uc2k3 + InvCs2uck3);
feq[4] = rhowk1 * (var0InvCs4uc2k3 - InvCs2uck3);
feq[5] = rhowk1 * (var0InvCs4uc2k5 + InvCs2uck5);
feq[6] = rhowk1 * (var0InvCs4uc2k5 - InvCs2uck5);
feq[7] = rhowk7 * (var0InvCs4uc2k7 + InvCs2uck7);
feq[8] = rhowk7 * (var0InvCs4uc2k7 - InvCs2uck7);
feq[9] = rhowk7 * (var0InvCs4uc2k9 + InvCs2uck9);
feq[10] = rhowk7 * (var0InvCs4uc2k9 - InvCs2uck9);
feq[11] = rhowk7 * (var0InvCs4uc2k11 + InvCs2uck11);
feq[12] = rhowk7 * (var0InvCs4uc2k11 - InvCs2uck11);
feq[13] = rhowk7 * (var0InvCs4uc2k13 + InvCs2uck13);
feq[14] = rhowk7 * (var0InvCs4uc2k13 - InvCs2uck13);
feq[15] = rhowk7 * (var0InvCs4uc2k15 + InvCs2uck15);
feq[16] = rhowk7 * (var0InvCs4uc2k15 - InvCs2uck15);
feq[17] = rhowk7 * (var0InvCs4uc2k17 + InvCs2uck17);
feq[18] = rhowk7 * (var0InvCs4uc2k17 - InvCs2uck17);
}
};

template <typename T, typename TypePack>
struct SecondOrderImpl<CELL<T, D3Q27<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;

__any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u){
constexpr T InvCs2x = T{0.5} * LatSet::InvCs2;
const T var0 = T{1} - (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) * InvCs2x;
const T rhowk0 = latset::w<LatSet>(0) * rho;
const T rhowk1 = latset::w<LatSet>(1) * rho;
const T rhowk7 = latset::w<LatSet>(7) * rho;
const T rhowk19 = latset::w<LatSet>(19) * rho;
const T InvCs2uck1 = LatSet::InvCs2 * (u[0]);
const T var0InvCs4uc2k1 = var0 + InvCs2uck1*InvCs2uck1*T{0.5};
const T InvCs2uck3 = LatSet::InvCs2 * (u[1]);
const T var0InvCs4uc2k3 = var0 + InvCs2uck3*InvCs2uck3*T{0.5};
const T InvCs2uck5 = LatSet::InvCs2 * (u[2]);
const T var0InvCs4uc2k5 = var0 + InvCs2uck5*InvCs2uck5*T{0.5};
const T InvCs2uck7 = LatSet::InvCs2 * (u[0]+u[1]);
const T var0InvCs4uc2k7 = var0 + InvCs2uck7*InvCs2uck7*T{0.5};
const T InvCs2uck9 = LatSet::InvCs2 * (u[0]+u[2]);
const T var0InvCs4uc2k9 = var0 + InvCs2uck9*InvCs2uck9*T{0.5};
const T InvCs2uck11 = LatSet::InvCs2 * (u[1]+u[2]);
const T var0InvCs4uc2k11 = var0 + InvCs2uck11*InvCs2uck11*T{0.5};
const T InvCs2uck13 = LatSet::InvCs2 * (u[0]-u[1]);
const T var0InvCs4uc2k13 = var0 + InvCs2uck13*InvCs2uck13*T{0.5};
const T InvCs2uck15 = LatSet::InvCs2 * (u[0]-u[2]);
const T var0InvCs4uc2k15 = var0 + InvCs2uck15*InvCs2uck15*T{0.5};
const T InvCs2uck17 = LatSet::InvCs2 * (u[1]-u[2]);
const T var0InvCs4uc2k17 = var0 + InvCs2uck17*InvCs2uck17*T{0.5};
const T InvCs2uck19 = LatSet::InvCs2 * (u[0]+u[1]+u[2]);
const T var0InvCs4uc2k19 = var0 + InvCs2uck19*InvCs2uck19*T{0.5};
const T InvCs2uck21 = LatSet::InvCs2 * (u[0]+u[1]-u[2]);
const T var0InvCs4uc2k21 = var0 + InvCs2uck21*InvCs2uck21*T{0.5};
const T InvCs2uck23 = LatSet::InvCs2 * (u[0]-u[1]+u[2]);
const T var0InvCs4uc2k23 = var0 + InvCs2uck23*InvCs2uck23*T{0.5};
const T InvCs2uck25 = LatSet::InvCs2 * (-u[0]+u[1]+u[2]);
const T var0InvCs4uc2k25 = var0 + InvCs2uck25*InvCs2uck25*T{0.5};
feq[0] = rhowk0 * var0;
feq[1] = rhowk1 * (var0InvCs4uc2k1 + InvCs2uck1);
feq[2] = rhowk1 * (var0InvCs4uc2k1 - InvCs2uck1);
feq[3] = rhowk1 * (var0InvCs4uc2k3 + InvCs2uck3);
feq[4] = rhowk1 * (var0InvCs4uc2k3 - InvCs2uck3);
feq[5] = rhowk1 * (var0InvCs4uc2k5 + InvCs2uck5);
feq[6] = rhowk1 * (var0InvCs4uc2k5 - InvCs2uck5);
feq[7] = rhowk7 * (var0InvCs4uc2k7 + InvCs2uck7);
feq[8] = rhowk7 * (var0InvCs4uc2k7 - InvCs2uck7);
feq[9] = rhowk7 * (var0InvCs4uc2k9 + InvCs2uck9);
feq[10] = rhowk7 * (var0InvCs4uc2k9 - InvCs2uck9);
feq[11] = rhowk7 * (var0InvCs4uc2k11 + InvCs2uck11);
feq[12] = rhowk7 * (var0InvCs4uc2k11 - InvCs2uck11);
feq[13] = rhowk7 * (var0InvCs4uc2k13 + InvCs2uck13);
feq[14] = rhowk7 * (var0InvCs4uc2k13 - InvCs2uck13);
feq[15] = rhowk7 * (var0InvCs4uc2k15 + InvCs2uck15);
feq[16] = rhowk7 * (var0InvCs4uc2k15 - InvCs2uck15);
feq[17] = rhowk7 * (var0InvCs4uc2k17 + InvCs2uck17);
feq[18] = rhowk7 * (var0InvCs4uc2k17 - InvCs2uck17);
feq[19] = rhowk19 * (var0InvCs4uc2k19 + InvCs2uck19);
feq[20] = rhowk19 * (var0InvCs4uc2k19 - InvCs2uck19);
feq[21] = rhowk19 * (var0InvCs4uc2k21 + InvCs2uck21);
feq[22] = rhowk19 * (var0InvCs4uc2k21 - InvCs2uck21);
feq[23] = rhowk19 * (var0InvCs4uc2k23 + InvCs2uck23);
feq[24] = rhowk19 * (var0InvCs4uc2k23 - InvCs2uck23);
feq[25] = rhowk19 * (var0InvCs4uc2k25 + InvCs2uck25);
feq[26] = rhowk19 * (var0InvCs4uc2k25 - InvCs2uck25);
}
};


//------------------------------------


}  // namespace equilibrium

#endif