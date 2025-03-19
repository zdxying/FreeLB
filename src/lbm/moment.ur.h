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



//------------------------------------

template <typename T, typename TypePack, bool WriteToField>
struct rhoImpl<CELL<T, D2Q5<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4];
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};

template <typename T, typename TypePack, bool WriteToField>
struct rhoImpl<CELL<T, D2Q9<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8];
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};

template <typename T, typename TypePack, bool WriteToField>
struct rhoImpl<CELL<T, D3Q7<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6];
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};

template <typename T, typename TypePack, bool WriteToField>
struct rhoImpl<CELL<T, D3Q15<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14];
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};

template <typename T, typename TypePack, bool WriteToField>
struct rhoImpl<CELL<T, D3Q19<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18];
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};

template <typename T, typename TypePack, bool WriteToField>
struct rhoImpl<CELL<T, D3Q27<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26];
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};


//------------------------------------


//------------------------------------

template <typename T, typename TypePack, typename SOURCE, bool WriteToField>
struct sourcerhoImpl<CELL<T, D2Q5<T>, TypePack>, SOURCE, WriteToField>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, const T source, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+source * T{0.5};
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};

template <typename T, typename TypePack, typename SOURCE, bool WriteToField>
struct sourcerhoImpl<CELL<T, D2Q9<T>, TypePack>, SOURCE, WriteToField>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, const T source, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+source * T{0.5};
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};

template <typename T, typename TypePack, typename SOURCE, bool WriteToField>
struct sourcerhoImpl<CELL<T, D3Q7<T>, TypePack>, SOURCE, WriteToField>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, const T source, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+source * T{0.5};
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};

template <typename T, typename TypePack, typename SOURCE, bool WriteToField>
struct sourcerhoImpl<CELL<T, D3Q15<T>, TypePack>, SOURCE, WriteToField>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, const T source, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+source * T{0.5};
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};

template <typename T, typename TypePack, typename SOURCE, bool WriteToField>
struct sourcerhoImpl<CELL<T, D3Q19<T>, TypePack>, SOURCE, WriteToField>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, const T source, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18]+source * T{0.5};
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};

template <typename T, typename TypePack, typename SOURCE, bool WriteToField>
struct sourcerhoImpl<CELL<T, D3Q27<T>, TypePack>, SOURCE, WriteToField>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, const T source, T& rho_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]+source * T{0.5};
if constexpr (WriteToField) cell.template get<GenericRho>() = rho_value;
}
};


//------------------------------------


//------------------------------------

template <typename T, typename TypePack, bool WriteToField>
struct UImpl<CELL<T, D2Q5<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4];
u_value[0] = (cell[1]-cell[2]) / rho_value;
u_value[1] = (cell[3]-cell[4]) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};

template <typename T, typename TypePack, bool WriteToField>
struct UImpl<CELL<T, D2Q9<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8];
u_value[0] = (cell[1]-cell[2]+cell[5]-cell[6]+cell[7]-cell[8]) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[5]-cell[6]-cell[7]+cell[8]) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};

template <typename T, typename TypePack, bool WriteToField>
struct UImpl<CELL<T, D3Q7<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6];
u_value[0] = (cell[1]-cell[2]) / rho_value;
u_value[1] = (cell[3]-cell[4]) / rho_value;
u_value[2] = (cell[5]-cell[6]) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};

template <typename T, typename TypePack, bool WriteToField>
struct UImpl<CELL<T, D3Q15<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[11]-cell[12]-cell[13]+cell[14]) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[9]-cell[10]-cell[11]+cell[12]+cell[13]-cell[14]) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[7]-cell[8]-cell[9]+cell[10]+cell[11]-cell[12]+cell[13]-cell[14]) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};

template <typename T, typename TypePack, bool WriteToField>
struct UImpl<CELL<T, D3Q19<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16]) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18]) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18]) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};

template <typename T, typename TypePack, bool WriteToField>
struct UImpl<CELL<T, D3Q27<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static inline void apply(CELLTYPE& cell, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16]+cell[19]-cell[20]+cell[21]-cell[22]+cell[23]-cell[24]-cell[25]+cell[26]) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18]+cell[19]-cell[20]+cell[21]-cell[22]-cell[23]+cell[24]+cell[25]-cell[26]) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18]+cell[19]-cell[20]-cell[21]+cell[22]+cell[23]-cell[24]+cell[25]-cell[26]) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};


//------------------------------------


//------------------------------------

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forceUImpl<CELL<T, D2Q5<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4];
u_value[0] = (cell[1]-cell[2] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4] + f_alpha[1]*T{0.5}) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
__any__ static inline void apply(CELLTYPE& cell, const T f, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4];
u_value[0] = cell[1]-cell[2];
u_value[1] = cell[3]-cell[4];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forceUImpl<CELL<T, D2Q9<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8];
u_value[0] = (cell[1]-cell[2]+cell[5]-cell[6]+cell[7]-cell[8] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[5]-cell[6]-cell[7]+cell[8] + f_alpha[1]*T{0.5}) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
__any__ static inline void apply(CELLTYPE& cell, const T f, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8];
u_value[0] = cell[1]-cell[2]+cell[5]-cell[6]+cell[7]-cell[8];
u_value[1] = cell[3]-cell[4]+cell[5]-cell[6]-cell[7]+cell[8];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forceUImpl<CELL<T, D3Q7<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6];
u_value[0] = (cell[1]-cell[2] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4] + f_alpha[1]*T{0.5}) / rho_value;
u_value[2] = (cell[5]-cell[6] + f_alpha[2]*T{0.5}) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
__any__ static inline void apply(CELLTYPE& cell, const T f, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6];
u_value[0] = cell[1]-cell[2];
u_value[1] = cell[3]-cell[4];
u_value[2] = cell[5]-cell[6];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
u_value[2] /= rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forceUImpl<CELL<T, D3Q15<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[11]-cell[12]-cell[13]+cell[14] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[9]-cell[10]-cell[11]+cell[12]+cell[13]-cell[14] + f_alpha[1]*T{0.5}) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[7]-cell[8]-cell[9]+cell[10]+cell[11]-cell[12]+cell[13]-cell[14] + f_alpha[2]*T{0.5}) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
__any__ static inline void apply(CELLTYPE& cell, const T f, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14];
u_value[0] = cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[11]-cell[12]-cell[13]+cell[14];
u_value[1] = cell[3]-cell[4]+cell[7]-cell[8]+cell[9]-cell[10]-cell[11]+cell[12]+cell[13]-cell[14];
u_value[2] = cell[5]-cell[6]+cell[7]-cell[8]-cell[9]+cell[10]+cell[11]-cell[12]+cell[13]-cell[14];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
u_value[2] /= rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forceUImpl<CELL<T, D3Q19<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18] + f_alpha[1]*T{0.5}) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18] + f_alpha[2]*T{0.5}) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
__any__ static inline void apply(CELLTYPE& cell, const T f, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18];
u_value[0] = cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16];
u_value[1] = cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18];
u_value[2] = cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
u_value[2] /= rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forceUImpl<CELL<T, D3Q27<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16]+cell[19]-cell[20]+cell[21]-cell[22]+cell[23]-cell[24]-cell[25]+cell[26] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18]+cell[19]-cell[20]+cell[21]-cell[22]-cell[23]+cell[24]+cell[25]-cell[26] + f_alpha[1]*T{0.5}) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18]+cell[19]-cell[20]-cell[21]+cell[22]+cell[23]-cell[24]+cell[25]-cell[26] + f_alpha[2]*T{0.5}) / rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
__any__ static inline void apply(CELLTYPE& cell, const T f, Vector<T, LatSet::d>& u_value){
const T rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26];
u_value[0] = cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16]+cell[19]-cell[20]+cell[21]-cell[22]+cell[23]-cell[24]-cell[25]+cell[26];
u_value[1] = cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18]+cell[19]-cell[20]+cell[21]-cell[22]-cell[23]+cell[24]+cell[25]-cell[26];
u_value[2] = cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18]+cell[19]-cell[20]-cell[21]+cell[22]+cell[23]-cell[24]+cell[25]-cell[26];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
u_value[2] /= rho_value;
if constexpr (WriteToField) cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}
};


//------------------------------------


//------------------------------------

template <typename T, typename TypePack, bool WriteToField>
struct rhoUImpl<CELL<T, D2Q5<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static void apply(CELLTYPE& cell, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4];
u_value[0] = (cell[1]-cell[2]) / rho_value;
u_value[1] = (cell[3]-cell[4]) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};

template <typename T, typename TypePack, bool WriteToField>
struct rhoUImpl<CELL<T, D2Q9<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static void apply(CELLTYPE& cell, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8];
u_value[0] = (cell[1]-cell[2]+cell[5]-cell[6]+cell[7]-cell[8]) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[5]-cell[6]-cell[7]+cell[8]) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};

template <typename T, typename TypePack, bool WriteToField>
struct rhoUImpl<CELL<T, D3Q7<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static void apply(CELLTYPE& cell, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6];
u_value[0] = (cell[1]-cell[2]) / rho_value;
u_value[1] = (cell[3]-cell[4]) / rho_value;
u_value[2] = (cell[5]-cell[6]) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};

template <typename T, typename TypePack, bool WriteToField>
struct rhoUImpl<CELL<T, D3Q15<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static void apply(CELLTYPE& cell, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[11]-cell[12]-cell[13]+cell[14]) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[9]-cell[10]-cell[11]+cell[12]+cell[13]-cell[14]) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[7]-cell[8]-cell[9]+cell[10]+cell[11]-cell[12]+cell[13]-cell[14]) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};

template <typename T, typename TypePack, bool WriteToField>
struct rhoUImpl<CELL<T, D3Q19<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static void apply(CELLTYPE& cell, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16]) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18]) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18]) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};

template <typename T, typename TypePack, bool WriteToField>
struct rhoUImpl<CELL<T, D3Q27<T>, TypePack>, WriteToField>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;
using GenericRho = typename CELLTYPE::GenericRho;

__any__ static void apply(CELLTYPE& cell, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16]+cell[19]-cell[20]+cell[21]-cell[22]+cell[23]-cell[24]-cell[25]+cell[26]) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18]+cell[19]-cell[20]+cell[21]-cell[22]-cell[23]+cell[24]+cell[25]-cell[26]) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18]+cell[19]-cell[20]-cell[21]+cell[22]+cell[23]-cell[24]+cell[25]-cell[26]) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};


//------------------------------------


//------------------------------------

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forcerhoUImpl<CELL<T, D2Q5<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4];
u_value[0] = (cell[1]-cell[2] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4] + f_alpha[1]*T{0.5}) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
__any__ static inline void apply(CELLTYPE& cell, const T f, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4];
u_value[0] = cell[1]-cell[2];
u_value[1] = cell[3]-cell[4];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forcerhoUImpl<CELL<T, D2Q9<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8];
u_value[0] = (cell[1]-cell[2]+cell[5]-cell[6]+cell[7]-cell[8] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[5]-cell[6]-cell[7]+cell[8] + f_alpha[1]*T{0.5}) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
__any__ static inline void apply(CELLTYPE& cell, const T f, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8];
u_value[0] = cell[1]-cell[2]+cell[5]-cell[6]+cell[7]-cell[8];
u_value[1] = cell[3]-cell[4]+cell[5]-cell[6]-cell[7]+cell[8];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forcerhoUImpl<CELL<T, D3Q7<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6];
u_value[0] = (cell[1]-cell[2] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4] + f_alpha[1]*T{0.5}) / rho_value;
u_value[2] = (cell[5]-cell[6] + f_alpha[2]*T{0.5}) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
__any__ static inline void apply(CELLTYPE& cell, const T f, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6];
u_value[0] = cell[1]-cell[2];
u_value[1] = cell[3]-cell[4];
u_value[2] = cell[5]-cell[6];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
u_value[2] /= rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forcerhoUImpl<CELL<T, D3Q15<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[11]-cell[12]-cell[13]+cell[14] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[9]-cell[10]-cell[11]+cell[12]+cell[13]-cell[14] + f_alpha[1]*T{0.5}) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[7]-cell[8]-cell[9]+cell[10]+cell[11]-cell[12]+cell[13]-cell[14] + f_alpha[2]*T{0.5}) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
__any__ static inline void apply(CELLTYPE& cell, const T f, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14];
u_value[0] = cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[11]-cell[12]-cell[13]+cell[14];
u_value[1] = cell[3]-cell[4]+cell[7]-cell[8]+cell[9]-cell[10]-cell[11]+cell[12]+cell[13]-cell[14];
u_value[2] = cell[5]-cell[6]+cell[7]-cell[8]-cell[9]+cell[10]+cell[11]-cell[12]+cell[13]-cell[14];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
u_value[2] /= rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forcerhoUImpl<CELL<T, D3Q19<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18] + f_alpha[1]*T{0.5}) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18] + f_alpha[2]*T{0.5}) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
__any__ static inline void apply(CELLTYPE& cell, const T f, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18];
u_value[0] = cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16];
u_value[1] = cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18];
u_value[2] = cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
u_value[2] /= rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};

template <typename T, typename TypePack, typename ForceScheme, bool WriteToField, unsigned int dir>
struct forcerhoUImpl<CELL<T, D3Q27<T>, TypePack>, ForceScheme, WriteToField, dir>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;
using GenericRho = typename CELLTYPE::GenericRho;

static constexpr unsigned int scalardir = dir >= 2 ? LatSet::d - 1 : dir;

__any__ static inline void apply(CELLTYPE& cell, const Vector<T, LatSet::d>& f_alpha, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26];
u_value[0] = (cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16]+cell[19]-cell[20]+cell[21]-cell[22]+cell[23]-cell[24]-cell[25]+cell[26] + f_alpha[0]*T{0.5}) / rho_value;
u_value[1] = (cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18]+cell[19]-cell[20]+cell[21]-cell[22]-cell[23]+cell[24]+cell[25]-cell[26] + f_alpha[1]*T{0.5}) / rho_value;
u_value[2] = (cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18]+cell[19]-cell[20]-cell[21]+cell[22]+cell[23]-cell[24]+cell[25]-cell[26] + f_alpha[2]*T{0.5}) / rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
__any__ static inline void apply(CELLTYPE& cell, const T f, T& rho_value, Vector<T, LatSet::d>& u_value){
rho_value = cell[0]+cell[1]+cell[2]+cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26];
u_value[0] = cell[1]-cell[2]+cell[7]-cell[8]+cell[9]-cell[10]+cell[13]-cell[14]+cell[15]-cell[16]+cell[19]-cell[20]+cell[21]-cell[22]+cell[23]-cell[24]-cell[25]+cell[26];
u_value[1] = cell[3]-cell[4]+cell[7]-cell[8]+cell[11]-cell[12]-cell[13]+cell[14]+cell[17]-cell[18]+cell[19]-cell[20]+cell[21]-cell[22]-cell[23]+cell[24]+cell[25]-cell[26];
u_value[2] = cell[5]-cell[6]+cell[9]-cell[10]+cell[11]-cell[12]-cell[15]+cell[16]-cell[17]+cell[18]+cell[19]-cell[20]-cell[21]+cell[22]+cell[23]-cell[24]+cell[25]-cell[26];
u_value[scalardir] += f * T{0.5};
u_value[0] /= rho_value;
u_value[1] /= rho_value;
u_value[2] /= rho_value;
if constexpr (WriteToField) {cell.template get<GenericRho>() = rho_value; cell.template get<VELOCITY<T, LatSet::d>>() = u_value;}
}
};


//------------------------------------


//------------------------------------

template <typename T, typename TypePack>
struct Pi_ab_neq<CELL<T, D2Q5<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
tensor[0] = cell[1]+cell[2]- rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = - rho*(u[0]*u[1]);
tensor[2] = cell[3]+cell[4]- rho*(u[1]*u[1]+LatSet::cs2);
}
};

template <typename T, typename TypePack>
struct Pi_ab_neq<CELL<T, D2Q9<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
tensor[0] = cell[1]+cell[2]+cell[5]+cell[6]+cell[7]+cell[8]- rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = cell[5]+cell[6]-cell[7]-cell[8]- rho*(u[0]*u[1]);
tensor[2] = cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]- rho*(u[1]*u[1]+LatSet::cs2);
}
};

template <typename T, typename TypePack>
struct Pi_ab_neq<CELL<T, D3Q7<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
tensor[0] = cell[1]+cell[2]- rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = - rho*(u[0]*u[1]);
tensor[2] = - rho*(u[0]*u[2]);
tensor[3] = cell[3]+cell[4]- rho*(u[1]*u[1]+LatSet::cs2);
tensor[4] = - rho*(u[1]*u[2]);
tensor[5] = cell[5]+cell[6]- rho*(u[2]*u[2]+LatSet::cs2);
}
};

template <typename T, typename TypePack>
struct Pi_ab_neq<CELL<T, D3Q15<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
tensor[0] = cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]- rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = cell[7]+cell[8]+cell[9]+cell[10]-cell[11]-cell[12]-cell[13]-cell[14]- rho*(u[0]*u[1]);
tensor[2] = cell[7]+cell[8]-cell[9]-cell[10]+cell[11]+cell[12]-cell[13]-cell[14]- rho*(u[0]*u[2]);
tensor[3] = cell[3]+cell[4]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]- rho*(u[1]*u[1]+LatSet::cs2);
tensor[4] = cell[7]+cell[8]-cell[9]-cell[10]-cell[11]-cell[12]+cell[13]+cell[14]- rho*(u[1]*u[2]);
tensor[5] = cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]- rho*(u[2]*u[2]+LatSet::cs2);
}
};

template <typename T, typename TypePack>
struct Pi_ab_neq<CELL<T, D3Q19<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
tensor[0] = cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[13]+cell[14]+cell[15]+cell[16]- rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = cell[7]+cell[8]-cell[13]-cell[14]- rho*(u[0]*u[1]);
tensor[2] = cell[9]+cell[10]-cell[15]-cell[16]- rho*(u[0]*u[2]);
tensor[3] = cell[3]+cell[4]+cell[7]+cell[8]+cell[11]+cell[12]+cell[13]+cell[14]+cell[17]+cell[18]- rho*(u[1]*u[1]+LatSet::cs2);
tensor[4] = cell[11]+cell[12]-cell[17]-cell[18]- rho*(u[1]*u[2]);
tensor[5] = cell[5]+cell[6]+cell[9]+cell[10]+cell[11]+cell[12]+cell[15]+cell[16]+cell[17]+cell[18]- rho*(u[2]*u[2]+LatSet::cs2);
}
};

template <typename T, typename TypePack>
struct Pi_ab_neq<CELL<T, D3Q27<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
tensor[0] = cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[13]+cell[14]+cell[15]+cell[16]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]- rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = cell[7]+cell[8]-cell[13]-cell[14]+cell[19]+cell[20]+cell[21]+cell[22]-cell[23]-cell[24]-cell[25]-cell[26]- rho*(u[0]*u[1]);
tensor[2] = cell[9]+cell[10]-cell[15]-cell[16]+cell[19]+cell[20]-cell[21]-cell[22]+cell[23]+cell[24]-cell[25]-cell[26]- rho*(u[0]*u[2]);
tensor[3] = cell[3]+cell[4]+cell[7]+cell[8]+cell[11]+cell[12]+cell[13]+cell[14]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]- rho*(u[1]*u[1]+LatSet::cs2);
tensor[4] = cell[11]+cell[12]-cell[17]-cell[18]+cell[19]+cell[20]-cell[21]-cell[22]-cell[23]-cell[24]+cell[25]+cell[26]- rho*(u[1]*u[2]);
tensor[5] = cell[5]+cell[6]+cell[9]+cell[10]+cell[11]+cell[12]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]- rho*(u[2]*u[2]+LatSet::cs2);
}
};


//------------------------------------


//------------------------------------

template <typename T, typename TypePack>
struct forcePi_ab_neq<CELL<T, D2Q5<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u_f, const Vector<T, LatSet::d>& f_alpha, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const Vector<T, LatSet::d> u = u_f - f_alpha * T{0.5};
tensor[0] = cell[1]+cell[2]+ T{0.5}*(f_alpha[0]*u[0]+f_alpha[0]*u[0]) - rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = T{0.5}*(f_alpha[0]*u[1]+f_alpha[1]*u[0]) - rho*(u[0]*u[1]);
tensor[2] = cell[3]+cell[4]+ T{0.5}*(f_alpha[1]*u[1]+f_alpha[1]*u[1]) - rho*(u[1]*u[1]+LatSet::cs2);
}
};

template <typename T, typename TypePack>
struct forcePi_ab_neq<CELL<T, D2Q9<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u_f, const Vector<T, LatSet::d>& f_alpha, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const Vector<T, LatSet::d> u = u_f - f_alpha * T{0.5};
tensor[0] = cell[1]+cell[2]+cell[5]+cell[6]+cell[7]+cell[8]+ T{0.5}*(f_alpha[0]*u[0]+f_alpha[0]*u[0]) - rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = cell[5]+cell[6]-cell[7]-cell[8]+ T{0.5}*(f_alpha[0]*u[1]+f_alpha[1]*u[0]) - rho*(u[0]*u[1]);
tensor[2] = cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]+ T{0.5}*(f_alpha[1]*u[1]+f_alpha[1]*u[1]) - rho*(u[1]*u[1]+LatSet::cs2);
}
};

template <typename T, typename TypePack>
struct forcePi_ab_neq<CELL<T, D3Q7<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u_f, const Vector<T, LatSet::d>& f_alpha, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const Vector<T, LatSet::d> u = u_f - f_alpha * T{0.5};
tensor[0] = cell[1]+cell[2]+ T{0.5}*(f_alpha[0]*u[0]+f_alpha[0]*u[0]) - rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = T{0.5}*(f_alpha[0]*u[1]+f_alpha[1]*u[0]) - rho*(u[0]*u[1]);
tensor[2] = T{0.5}*(f_alpha[0]*u[2]+f_alpha[2]*u[0]) - rho*(u[0]*u[2]);
tensor[3] = cell[3]+cell[4]+ T{0.5}*(f_alpha[1]*u[1]+f_alpha[1]*u[1]) - rho*(u[1]*u[1]+LatSet::cs2);
tensor[4] = T{0.5}*(f_alpha[1]*u[2]+f_alpha[2]*u[1]) - rho*(u[1]*u[2]);
tensor[5] = cell[5]+cell[6]+ T{0.5}*(f_alpha[2]*u[2]+f_alpha[2]*u[2]) - rho*(u[2]*u[2]+LatSet::cs2);
}
};

template <typename T, typename TypePack>
struct forcePi_ab_neq<CELL<T, D3Q15<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u_f, const Vector<T, LatSet::d>& f_alpha, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const Vector<T, LatSet::d> u = u_f - f_alpha * T{0.5};
tensor[0] = cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+ T{0.5}*(f_alpha[0]*u[0]+f_alpha[0]*u[0]) - rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = cell[7]+cell[8]+cell[9]+cell[10]-cell[11]-cell[12]-cell[13]-cell[14]+ T{0.5}*(f_alpha[0]*u[1]+f_alpha[1]*u[0]) - rho*(u[0]*u[1]);
tensor[2] = cell[7]+cell[8]-cell[9]-cell[10]+cell[11]+cell[12]-cell[13]-cell[14]+ T{0.5}*(f_alpha[0]*u[2]+f_alpha[2]*u[0]) - rho*(u[0]*u[2]);
tensor[3] = cell[3]+cell[4]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+ T{0.5}*(f_alpha[1]*u[1]+f_alpha[1]*u[1]) - rho*(u[1]*u[1]+LatSet::cs2);
tensor[4] = cell[7]+cell[8]-cell[9]-cell[10]-cell[11]-cell[12]+cell[13]+cell[14]+ T{0.5}*(f_alpha[1]*u[2]+f_alpha[2]*u[1]) - rho*(u[1]*u[2]);
tensor[5] = cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]+ T{0.5}*(f_alpha[2]*u[2]+f_alpha[2]*u[2]) - rho*(u[2]*u[2]+LatSet::cs2);
}
};

template <typename T, typename TypePack>
struct forcePi_ab_neq<CELL<T, D3Q19<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u_f, const Vector<T, LatSet::d>& f_alpha, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const Vector<T, LatSet::d> u = u_f - f_alpha * T{0.5};
tensor[0] = cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[13]+cell[14]+cell[15]+cell[16]+ T{0.5}*(f_alpha[0]*u[0]+f_alpha[0]*u[0]) - rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = cell[7]+cell[8]-cell[13]-cell[14]+ T{0.5}*(f_alpha[0]*u[1]+f_alpha[1]*u[0]) - rho*(u[0]*u[1]);
tensor[2] = cell[9]+cell[10]-cell[15]-cell[16]+ T{0.5}*(f_alpha[0]*u[2]+f_alpha[2]*u[0]) - rho*(u[0]*u[2]);
tensor[3] = cell[3]+cell[4]+cell[7]+cell[8]+cell[11]+cell[12]+cell[13]+cell[14]+cell[17]+cell[18]+ T{0.5}*(f_alpha[1]*u[1]+f_alpha[1]*u[1]) - rho*(u[1]*u[1]+LatSet::cs2);
tensor[4] = cell[11]+cell[12]-cell[17]-cell[18]+ T{0.5}*(f_alpha[1]*u[2]+f_alpha[2]*u[1]) - rho*(u[1]*u[2]);
tensor[5] = cell[5]+cell[6]+cell[9]+cell[10]+cell[11]+cell[12]+cell[15]+cell[16]+cell[17]+cell[18]+ T{0.5}*(f_alpha[2]*u[2]+f_alpha[2]*u[2]) - rho*(u[2]*u[2]+LatSet::cs2);
}
};

template <typename T, typename TypePack>
struct forcePi_ab_neq<CELL<T, D3Q27<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u_f, const Vector<T, LatSet::d>& f_alpha, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const Vector<T, LatSet::d> u = u_f - f_alpha * T{0.5};
tensor[0] = cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[13]+cell[14]+cell[15]+cell[16]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]+ T{0.5}*(f_alpha[0]*u[0]+f_alpha[0]*u[0]) - rho*(u[0]*u[0]+LatSet::cs2);
tensor[1] = cell[7]+cell[8]-cell[13]-cell[14]+cell[19]+cell[20]+cell[21]+cell[22]-cell[23]-cell[24]-cell[25]-cell[26]+ T{0.5}*(f_alpha[0]*u[1]+f_alpha[1]*u[0]) - rho*(u[0]*u[1]);
tensor[2] = cell[9]+cell[10]-cell[15]-cell[16]+cell[19]+cell[20]-cell[21]-cell[22]+cell[23]+cell[24]-cell[25]-cell[26]+ T{0.5}*(f_alpha[0]*u[2]+f_alpha[2]*u[0]) - rho*(u[0]*u[2]);
tensor[3] = cell[3]+cell[4]+cell[7]+cell[8]+cell[11]+cell[12]+cell[13]+cell[14]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]+ T{0.5}*(f_alpha[1]*u[1]+f_alpha[1]*u[1]) - rho*(u[1]*u[1]+LatSet::cs2);
tensor[4] = cell[11]+cell[12]-cell[17]-cell[18]+cell[19]+cell[20]-cell[21]-cell[22]-cell[23]-cell[24]+cell[25]+cell[26]+ T{0.5}*(f_alpha[1]*u[2]+f_alpha[2]*u[1]) - rho*(u[1]*u[2]);
tensor[5] = cell[5]+cell[6]+cell[9]+cell[10]+cell[11]+cell[12]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]+ T{0.5}*(f_alpha[2]*u[2]+f_alpha[2]*u[2]) - rho*(u[2]*u[2]+LatSet::cs2);
}
};


//------------------------------------


//------------------------------------

template <typename T, typename TypePack>
struct stress<CELL<T, D2Q5<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const T coeff = T{0.5} * cell.getOmega() - T{1};
tensor[0] = (cell[1]+cell[2]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (-rho*(u[0]*u[1])) * coeff;
tensor[2] = (cell[3]+cell[4]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
}
};

template <typename T, typename TypePack>
struct stress<CELL<T, D2Q9<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const T coeff = T{0.5} * cell.getOmega() - T{1};
tensor[0] = (cell[1]+cell[2]+cell[5]+cell[6]+cell[7]+cell[8]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (cell[5]+cell[6]-cell[7]-cell[8]-rho*(u[0]*u[1])) * coeff;
tensor[2] = (cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
}
};

template <typename T, typename TypePack>
struct stress<CELL<T, D3Q7<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const T coeff = T{0.5} * cell.getOmega() - T{1};
tensor[0] = (cell[1]+cell[2]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (-rho*(u[0]*u[1])) * coeff;
tensor[2] = (-rho*(u[0]*u[2])) * coeff;
tensor[3] = (cell[3]+cell[4]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
tensor[4] = (-rho*(u[1]*u[2])) * coeff;
tensor[5] = (cell[5]+cell[6]-rho*(u[2]*u[2]+LatSet::cs2)) * coeff;
}
};

template <typename T, typename TypePack>
struct stress<CELL<T, D3Q15<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const T coeff = T{0.5} * cell.getOmega() - T{1};
tensor[0] = (cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (cell[7]+cell[8]+cell[9]+cell[10]-cell[11]-cell[12]-cell[13]-cell[14]-rho*(u[0]*u[1])) * coeff;
tensor[2] = (cell[7]+cell[8]-cell[9]-cell[10]+cell[11]+cell[12]-cell[13]-cell[14]-rho*(u[0]*u[2])) * coeff;
tensor[3] = (cell[3]+cell[4]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
tensor[4] = (cell[7]+cell[8]-cell[9]-cell[10]-cell[11]-cell[12]+cell[13]+cell[14]-rho*(u[1]*u[2])) * coeff;
tensor[5] = (cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]-rho*(u[2]*u[2]+LatSet::cs2)) * coeff;
}
};

template <typename T, typename TypePack>
struct stress<CELL<T, D3Q19<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const T coeff = T{0.5} * cell.getOmega() - T{1};
tensor[0] = (cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[13]+cell[14]+cell[15]+cell[16]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (cell[7]+cell[8]-cell[13]-cell[14]-rho*(u[0]*u[1])) * coeff;
tensor[2] = (cell[9]+cell[10]-cell[15]-cell[16]-rho*(u[0]*u[2])) * coeff;
tensor[3] = (cell[3]+cell[4]+cell[7]+cell[8]+cell[11]+cell[12]+cell[13]+cell[14]+cell[17]+cell[18]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
tensor[4] = (cell[11]+cell[12]-cell[17]-cell[18]-rho*(u[1]*u[2])) * coeff;
tensor[5] = (cell[5]+cell[6]+cell[9]+cell[10]+cell[11]+cell[12]+cell[15]+cell[16]+cell[17]+cell[18]-rho*(u[2]*u[2]+LatSet::cs2)) * coeff;
}
};

template <typename T, typename TypePack>
struct stress<CELL<T, D3Q27<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
const T coeff = T{0.5} * cell.getOmega() - T{1};
tensor[0] = (cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[13]+cell[14]+cell[15]+cell[16]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (cell[7]+cell[8]-cell[13]-cell[14]+cell[19]+cell[20]+cell[21]+cell[22]-cell[23]-cell[24]-cell[25]-cell[26]-rho*(u[0]*u[1])) * coeff;
tensor[2] = (cell[9]+cell[10]-cell[15]-cell[16]+cell[19]+cell[20]-cell[21]-cell[22]+cell[23]+cell[24]-cell[25]-cell[26]-rho*(u[0]*u[2])) * coeff;
tensor[3] = (cell[3]+cell[4]+cell[7]+cell[8]+cell[11]+cell[12]+cell[13]+cell[14]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
tensor[4] = (cell[11]+cell[12]-cell[17]-cell[18]+cell[19]+cell[20]-cell[21]-cell[22]-cell[23]-cell[24]+cell[25]+cell[26]-rho*(u[1]*u[2])) * coeff;
tensor[5] = (cell[5]+cell[6]+cell[9]+cell[10]+cell[11]+cell[12]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]-rho*(u[2]*u[2]+LatSet::cs2)) * coeff;
}
};


//------------------------------------


//------------------------------------

template <typename T, typename TypePack>
struct strainRate<CELL<T, D2Q5<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
T omega{};
if constexpr(cell.template hasField<OMEGA<T>>()) omega = cell.template get<OMEGA<T>>();
else omega = cell.getOmega();
const T coeff = T{-1.5} * omega / cell.template get<typename CELLTYPE::GenericRho>();
tensor[0] = (cell[1]+cell[2]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (-rho*(u[0]*u[1])) * coeff;
tensor[2] = (cell[3]+cell[4]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
}
};

template <typename T, typename TypePack>
struct strainRate<CELL<T, D2Q9<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
T omega{};
if constexpr(cell.template hasField<OMEGA<T>>()) omega = cell.template get<OMEGA<T>>();
else omega = cell.getOmega();
const T coeff = T{-1.5} * omega / cell.template get<typename CELLTYPE::GenericRho>();
tensor[0] = (cell[1]+cell[2]+cell[5]+cell[6]+cell[7]+cell[8]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (cell[5]+cell[6]-cell[7]-cell[8]-rho*(u[0]*u[1])) * coeff;
tensor[2] = (cell[3]+cell[4]+cell[5]+cell[6]+cell[7]+cell[8]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
}
};

template <typename T, typename TypePack>
struct strainRate<CELL<T, D3Q7<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
T omega{};
if constexpr(cell.template hasField<OMEGA<T>>()) omega = cell.template get<OMEGA<T>>();
else omega = cell.getOmega();
const T coeff = T{-1.5} * omega / cell.template get<typename CELLTYPE::GenericRho>();
tensor[0] = (cell[1]+cell[2]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (-rho*(u[0]*u[1])) * coeff;
tensor[2] = (-rho*(u[0]*u[2])) * coeff;
tensor[3] = (cell[3]+cell[4]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
tensor[4] = (-rho*(u[1]*u[2])) * coeff;
tensor[5] = (cell[5]+cell[6]-rho*(u[2]*u[2]+LatSet::cs2)) * coeff;
}
};

template <typename T, typename TypePack>
struct strainRate<CELL<T, D3Q15<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
T omega{};
if constexpr(cell.template hasField<OMEGA<T>>()) omega = cell.template get<OMEGA<T>>();
else omega = cell.getOmega();
const T coeff = T{-1.5} * omega / cell.template get<typename CELLTYPE::GenericRho>();
tensor[0] = (cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (cell[7]+cell[8]+cell[9]+cell[10]-cell[11]-cell[12]-cell[13]-cell[14]-rho*(u[0]*u[1])) * coeff;
tensor[2] = (cell[7]+cell[8]-cell[9]-cell[10]+cell[11]+cell[12]-cell[13]-cell[14]-rho*(u[0]*u[2])) * coeff;
tensor[3] = (cell[3]+cell[4]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
tensor[4] = (cell[7]+cell[8]-cell[9]-cell[10]-cell[11]-cell[12]+cell[13]+cell[14]-rho*(u[1]*u[2])) * coeff;
tensor[5] = (cell[5]+cell[6]+cell[7]+cell[8]+cell[9]+cell[10]+cell[11]+cell[12]+cell[13]+cell[14]-rho*(u[2]*u[2]+LatSet::cs2)) * coeff;
}
};

template <typename T, typename TypePack>
struct strainRate<CELL<T, D3Q19<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
T omega{};
if constexpr(cell.template hasField<OMEGA<T>>()) omega = cell.template get<OMEGA<T>>();
else omega = cell.getOmega();
const T coeff = T{-1.5} * omega / cell.template get<typename CELLTYPE::GenericRho>();
tensor[0] = (cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[13]+cell[14]+cell[15]+cell[16]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (cell[7]+cell[8]-cell[13]-cell[14]-rho*(u[0]*u[1])) * coeff;
tensor[2] = (cell[9]+cell[10]-cell[15]-cell[16]-rho*(u[0]*u[2])) * coeff;
tensor[3] = (cell[3]+cell[4]+cell[7]+cell[8]+cell[11]+cell[12]+cell[13]+cell[14]+cell[17]+cell[18]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
tensor[4] = (cell[11]+cell[12]-cell[17]-cell[18]-rho*(u[1]*u[2])) * coeff;
tensor[5] = (cell[5]+cell[6]+cell[9]+cell[10]+cell[11]+cell[12]+cell[15]+cell[16]+cell[17]+cell[18]-rho*(u[2]*u[2]+LatSet::cs2)) * coeff;
}
};

template <typename T, typename TypePack>
struct strainRate<CELL<T, D3Q27<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;

__any__ static inline void apply(CELLTYPE& cell, const T rho, const Vector<T, LatSet::d>& u, std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& tensor){
T omega{};
if constexpr(cell.template hasField<OMEGA<T>>()) omega = cell.template get<OMEGA<T>>();
else omega = cell.getOmega();
const T coeff = T{-1.5} * omega / cell.template get<typename CELLTYPE::GenericRho>();
tensor[0] = (cell[1]+cell[2]+cell[7]+cell[8]+cell[9]+cell[10]+cell[13]+cell[14]+cell[15]+cell[16]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]-rho*(u[0]*u[0]+LatSet::cs2)) * coeff;
tensor[1] = (cell[7]+cell[8]-cell[13]-cell[14]+cell[19]+cell[20]+cell[21]+cell[22]-cell[23]-cell[24]-cell[25]-cell[26]-rho*(u[0]*u[1])) * coeff;
tensor[2] = (cell[9]+cell[10]-cell[15]-cell[16]+cell[19]+cell[20]-cell[21]-cell[22]+cell[23]+cell[24]-cell[25]-cell[26]-rho*(u[0]*u[2])) * coeff;
tensor[3] = (cell[3]+cell[4]+cell[7]+cell[8]+cell[11]+cell[12]+cell[13]+cell[14]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]-rho*(u[1]*u[1]+LatSet::cs2)) * coeff;
tensor[4] = (cell[11]+cell[12]-cell[17]-cell[18]+cell[19]+cell[20]-cell[21]-cell[22]-cell[23]-cell[24]+cell[25]+cell[26]-rho*(u[1]*u[2])) * coeff;
tensor[5] = (cell[5]+cell[6]+cell[9]+cell[10]+cell[11]+cell[12]+cell[15]+cell[16]+cell[17]+cell[18]+cell[19]+cell[20]+cell[21]+cell[22]+cell[23]+cell[24]+cell[25]+cell[26]-rho*(u[2]*u[2]+LatSet::cs2)) * coeff;
}
};


//------------------------------------


//------------------------------------

template <typename T, typename TypePack>
struct shearRateMagImpl<CELL<T, D2Q5<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q5<T>, TypePack>;
using LatSet = D2Q5<T>;

__any__ static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& srtensor){
T value = srtensor[0]*srtensor[0] + srtensor[1]*srtensor[1]*2 + srtensor[2]*srtensor[2] ;
return std::sqrt(2 * value);
}
};

template <typename T, typename TypePack>
struct shearRateMagImpl<CELL<T, D2Q9<T>, TypePack>>{
using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
using LatSet = D2Q9<T>;

__any__ static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& srtensor){
T value = srtensor[0]*srtensor[0] + srtensor[1]*srtensor[1]*2 + srtensor[2]*srtensor[2] ;
return std::sqrt(2 * value);
}
};

template <typename T, typename TypePack>
struct shearRateMagImpl<CELL<T, D3Q7<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q7<T>, TypePack>;
using LatSet = D3Q7<T>;

__any__ static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& srtensor){
T value = srtensor[0]*srtensor[0] + srtensor[1]*srtensor[1]*2 + srtensor[2]*srtensor[2]*2 + srtensor[3]*srtensor[3] + srtensor[4]*srtensor[4]*2 + srtensor[5]*srtensor[5] ;
return std::sqrt(2 * value);
}
};

template <typename T, typename TypePack>
struct shearRateMagImpl<CELL<T, D3Q15<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q15<T>, TypePack>;
using LatSet = D3Q15<T>;

__any__ static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& srtensor){
T value = srtensor[0]*srtensor[0] + srtensor[1]*srtensor[1]*2 + srtensor[2]*srtensor[2]*2 + srtensor[3]*srtensor[3] + srtensor[4]*srtensor[4]*2 + srtensor[5]*srtensor[5] ;
return std::sqrt(2 * value);
}
};

template <typename T, typename TypePack>
struct shearRateMagImpl<CELL<T, D3Q19<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q19<T>, TypePack>;
using LatSet = D3Q19<T>;

__any__ static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& srtensor){
T value = srtensor[0]*srtensor[0] + srtensor[1]*srtensor[1]*2 + srtensor[2]*srtensor[2]*2 + srtensor[3]*srtensor[3] + srtensor[4]*srtensor[4]*2 + srtensor[5]*srtensor[5] ;
return std::sqrt(2 * value);
}
};

template <typename T, typename TypePack>
struct shearRateMagImpl<CELL<T, D3Q27<T>, TypePack>>{
using CELLTYPE = CELL<T, D3Q27<T>, TypePack>;
using LatSet = D3Q27<T>;

__any__ static inline T get(const std::array<T, util::SymmetricMatrixSize<LatSet::d>()>& srtensor){
T value = srtensor[0]*srtensor[0] + srtensor[1]*srtensor[1]*2 + srtensor[2]*srtensor[2]*2 + srtensor[3]*srtensor[3] + srtensor[4]*srtensor[4]*2 + srtensor[5]*srtensor[5] ;
return std::sqrt(2 * value);
}
};


//------------------------------------



}  // namespace moment

#endif