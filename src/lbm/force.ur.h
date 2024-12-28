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

// force.ur.h
#pragma once

#include "lbm/force.h"

#ifdef _UNROLLFOR

namespace force {


//------------------------------------

template <typename T>
struct ForcePopImpl<T, D2Q5<T>>{
__any__ static inline void compute(std::array<T, 5> &Fi, const Vector<T, 2> &u, const Vector<T, 2> &F){
const T u0 = u[0];
const T u1 = u[1];
const T F0 = F[0];
const T F1 = F[1];
const T T1_u0 = T{1} - u0;
const T _T1_u0 = T{-1} - u0;
const T _u0 = -u0;
const T T1_u1 = T{1} - u1;
const T _T1_u1 = T{-1} - u1;
const T _u1 = -u1;
const T cuInvCs2_1 = D2Q5<T>::InvCs2 * (u0);
const T cuInvCs2_3 = D2Q5<T>::InvCs2 * (u1);
Fi[0] = latset::w<D2Q5<T>>(0)*D2Q5<T>::InvCs2 * (-(F0*u0+F1*u1));
Fi[1] = latset::w<D2Q5<T>>(1)*D2Q5<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_1) + F1*(_u1));
Fi[2] = latset::w<D2Q5<T>>(2)*D2Q5<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_1) + F1*(_u1));
Fi[3] = latset::w<D2Q5<T>>(3)*D2Q5<T>::InvCs2 * (F0*(_u0) + F1*(T1_u1+cuInvCs2_3));
Fi[4] = latset::w<D2Q5<T>>(4)*D2Q5<T>::InvCs2 * (F0*(_u0) + F1*(_T1_u1+cuInvCs2_3));
}
};

template <typename T>
struct ForcePopImpl<T, D2Q9<T>>{
__any__ static inline void compute(std::array<T, 9> &Fi, const Vector<T, 2> &u, const Vector<T, 2> &F){
const T u0 = u[0];
const T u1 = u[1];
const T F0 = F[0];
const T F1 = F[1];
const T T1_u0 = T{1} - u0;
const T _T1_u0 = T{-1} - u0;
const T _u0 = -u0;
const T T1_u1 = T{1} - u1;
const T _T1_u1 = T{-1} - u1;
const T _u1 = -u1;
const T cuInvCs2_1 = D2Q9<T>::InvCs2 * (u0);
const T cuInvCs2_3 = D2Q9<T>::InvCs2 * (u1);
const T cuInvCs2_5 = D2Q9<T>::InvCs2 * (u0+u1);
const T cuInvCs2_7 = D2Q9<T>::InvCs2 * (u0-u1);
Fi[0] = latset::w<D2Q9<T>>(0)*D2Q9<T>::InvCs2 * (-(F0*u0+F1*u1));
Fi[1] = latset::w<D2Q9<T>>(1)*D2Q9<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_1) + F1*(_u1));
Fi[2] = latset::w<D2Q9<T>>(2)*D2Q9<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_1) + F1*(_u1));
Fi[3] = latset::w<D2Q9<T>>(3)*D2Q9<T>::InvCs2 * (F0*(_u0) + F1*(T1_u1+cuInvCs2_3));
Fi[4] = latset::w<D2Q9<T>>(4)*D2Q9<T>::InvCs2 * (F0*(_u0) + F1*(_T1_u1+cuInvCs2_3));
Fi[5] = latset::w<D2Q9<T>>(5)*D2Q9<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_5) + F1*(T1_u1+cuInvCs2_5));
Fi[6] = latset::w<D2Q9<T>>(6)*D2Q9<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_5) + F1*(_T1_u1+cuInvCs2_5));
Fi[7] = latset::w<D2Q9<T>>(7)*D2Q9<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_7) + F1*(_T1_u1-cuInvCs2_7));
Fi[8] = latset::w<D2Q9<T>>(8)*D2Q9<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_7) + F1*(T1_u1-cuInvCs2_7));
}
};

template <typename T>
struct ForcePopImpl<T, D3Q7<T>>{
__any__ static inline void compute(std::array<T, 7> &Fi, const Vector<T, 3> &u, const Vector<T, 3> &F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T F0 = F[0];
const T F1 = F[1];
const T F2 = F[2];
const T T1_u0 = T{1} - u0;
const T _T1_u0 = T{-1} - u0;
const T _u0 = -u0;
const T T1_u1 = T{1} - u1;
const T _T1_u1 = T{-1} - u1;
const T _u1 = -u1;
const T T1_u2 = T{1} - u2;
const T _T1_u2 = T{-1} - u2;
const T _u2 = -u2;
const T cuInvCs2_1 = D3Q7<T>::InvCs2 * (u0);
const T cuInvCs2_3 = D3Q7<T>::InvCs2 * (u1);
const T cuInvCs2_5 = D3Q7<T>::InvCs2 * (u2);
Fi[0] = latset::w<D3Q7<T>>(0)*D3Q7<T>::InvCs2 * (-(F0*u0+F1*u1+F2*u2));
Fi[1] = latset::w<D3Q7<T>>(1)*D3Q7<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_1) + F1*(_u1) + F2*(_u2));
Fi[2] = latset::w<D3Q7<T>>(2)*D3Q7<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_1) + F1*(_u1) + F2*(_u2));
Fi[3] = latset::w<D3Q7<T>>(3)*D3Q7<T>::InvCs2 * (F0*(_u0) + F1*(T1_u1+cuInvCs2_3) + F2*(_u2));
Fi[4] = latset::w<D3Q7<T>>(4)*D3Q7<T>::InvCs2 * (F0*(_u0) + F1*(_T1_u1+cuInvCs2_3) + F2*(_u2));
Fi[5] = latset::w<D3Q7<T>>(5)*D3Q7<T>::InvCs2 * (F0*(_u0) + F1*(_u1) + F2*(T1_u2+cuInvCs2_5));
Fi[6] = latset::w<D3Q7<T>>(6)*D3Q7<T>::InvCs2 * (F0*(_u0) + F1*(_u1) + F2*(_T1_u2+cuInvCs2_5));
}
};

template <typename T>
struct ForcePopImpl<T, D3Q15<T>>{
__any__ static inline void compute(std::array<T, 15> &Fi, const Vector<T, 3> &u, const Vector<T, 3> &F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T F0 = F[0];
const T F1 = F[1];
const T F2 = F[2];
const T T1_u0 = T{1} - u0;
const T _T1_u0 = T{-1} - u0;
const T _u0 = -u0;
const T T1_u1 = T{1} - u1;
const T _T1_u1 = T{-1} - u1;
const T _u1 = -u1;
const T T1_u2 = T{1} - u2;
const T _T1_u2 = T{-1} - u2;
const T _u2 = -u2;
const T cuInvCs2_1 = D3Q15<T>::InvCs2 * (u0);
const T cuInvCs2_3 = D3Q15<T>::InvCs2 * (u1);
const T cuInvCs2_5 = D3Q15<T>::InvCs2 * (u2);
const T cuInvCs2_7 = D3Q15<T>::InvCs2 * (u0+u1+u2);
const T cuInvCs2_9 = D3Q15<T>::InvCs2 * (u0+u1-u2);
const T cuInvCs2_11 = D3Q15<T>::InvCs2 * (u0-u1+u2);
const T cuInvCs2_13 = D3Q15<T>::InvCs2 * (-u0+u1+u2);
Fi[0] = latset::w<D3Q15<T>>(0)*D3Q15<T>::InvCs2 * (-(F0*u0+F1*u1+F2*u2));
Fi[1] = latset::w<D3Q15<T>>(1)*D3Q15<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_1) + F1*(_u1) + F2*(_u2));
Fi[2] = latset::w<D3Q15<T>>(2)*D3Q15<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_1) + F1*(_u1) + F2*(_u2));
Fi[3] = latset::w<D3Q15<T>>(3)*D3Q15<T>::InvCs2 * (F0*(_u0) + F1*(T1_u1+cuInvCs2_3) + F2*(_u2));
Fi[4] = latset::w<D3Q15<T>>(4)*D3Q15<T>::InvCs2 * (F0*(_u0) + F1*(_T1_u1+cuInvCs2_3) + F2*(_u2));
Fi[5] = latset::w<D3Q15<T>>(5)*D3Q15<T>::InvCs2 * (F0*(_u0) + F1*(_u1) + F2*(T1_u2+cuInvCs2_5));
Fi[6] = latset::w<D3Q15<T>>(6)*D3Q15<T>::InvCs2 * (F0*(_u0) + F1*(_u1) + F2*(_T1_u2+cuInvCs2_5));
Fi[7] = latset::w<D3Q15<T>>(7)*D3Q15<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_7) + F1*(T1_u1+cuInvCs2_7) + F2*(T1_u2+cuInvCs2_7));
Fi[8] = latset::w<D3Q15<T>>(8)*D3Q15<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_7) + F1*(_T1_u1+cuInvCs2_7) + F2*(_T1_u2+cuInvCs2_7));
Fi[9] = latset::w<D3Q15<T>>(9)*D3Q15<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_9) + F1*(T1_u1+cuInvCs2_9) + F2*(_T1_u2-cuInvCs2_9));
Fi[10] = latset::w<D3Q15<T>>(10)*D3Q15<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_9) + F1*(_T1_u1+cuInvCs2_9) + F2*(T1_u2-cuInvCs2_9));
Fi[11] = latset::w<D3Q15<T>>(11)*D3Q15<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_11) + F1*(_T1_u1-cuInvCs2_11) + F2*(T1_u2+cuInvCs2_11));
Fi[12] = latset::w<D3Q15<T>>(12)*D3Q15<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_11) + F1*(T1_u1-cuInvCs2_11) + F2*(_T1_u2+cuInvCs2_11));
Fi[13] = latset::w<D3Q15<T>>(13)*D3Q15<T>::InvCs2 * (F0*(_T1_u0-cuInvCs2_13) + F1*(T1_u1+cuInvCs2_13) + F2*(T1_u2+cuInvCs2_13));
Fi[14] = latset::w<D3Q15<T>>(14)*D3Q15<T>::InvCs2 * (F0*(T1_u0-cuInvCs2_13) + F1*(_T1_u1+cuInvCs2_13) + F2*(_T1_u2+cuInvCs2_13));
}
};

template <typename T>
struct ForcePopImpl<T, D3Q19<T>>{
__any__ static inline void compute(std::array<T, 19> &Fi, const Vector<T, 3> &u, const Vector<T, 3> &F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T F0 = F[0];
const T F1 = F[1];
const T F2 = F[2];
const T T1_u0 = T{1} - u0;
const T _T1_u0 = T{-1} - u0;
const T _u0 = -u0;
const T T1_u1 = T{1} - u1;
const T _T1_u1 = T{-1} - u1;
const T _u1 = -u1;
const T T1_u2 = T{1} - u2;
const T _T1_u2 = T{-1} - u2;
const T _u2 = -u2;
const T cuInvCs2_1 = D3Q19<T>::InvCs2 * (u0);
const T cuInvCs2_3 = D3Q19<T>::InvCs2 * (u1);
const T cuInvCs2_5 = D3Q19<T>::InvCs2 * (u2);
const T cuInvCs2_7 = D3Q19<T>::InvCs2 * (u0+u1);
const T cuInvCs2_9 = D3Q19<T>::InvCs2 * (u0+u2);
const T cuInvCs2_11 = D3Q19<T>::InvCs2 * (u1+u2);
const T cuInvCs2_13 = D3Q19<T>::InvCs2 * (u0-u1);
const T cuInvCs2_15 = D3Q19<T>::InvCs2 * (u0-u2);
const T cuInvCs2_17 = D3Q19<T>::InvCs2 * (u1-u2);
Fi[0] = latset::w<D3Q19<T>>(0)*D3Q19<T>::InvCs2 * (-(F0*u0+F1*u1+F2*u2));
Fi[1] = latset::w<D3Q19<T>>(1)*D3Q19<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_1) + F1*(_u1) + F2*(_u2));
Fi[2] = latset::w<D3Q19<T>>(2)*D3Q19<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_1) + F1*(_u1) + F2*(_u2));
Fi[3] = latset::w<D3Q19<T>>(3)*D3Q19<T>::InvCs2 * (F0*(_u0) + F1*(T1_u1+cuInvCs2_3) + F2*(_u2));
Fi[4] = latset::w<D3Q19<T>>(4)*D3Q19<T>::InvCs2 * (F0*(_u0) + F1*(_T1_u1+cuInvCs2_3) + F2*(_u2));
Fi[5] = latset::w<D3Q19<T>>(5)*D3Q19<T>::InvCs2 * (F0*(_u0) + F1*(_u1) + F2*(T1_u2+cuInvCs2_5));
Fi[6] = latset::w<D3Q19<T>>(6)*D3Q19<T>::InvCs2 * (F0*(_u0) + F1*(_u1) + F2*(_T1_u2+cuInvCs2_5));
Fi[7] = latset::w<D3Q19<T>>(7)*D3Q19<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_7) + F1*(T1_u1+cuInvCs2_7) + F2*(_u2));
Fi[8] = latset::w<D3Q19<T>>(8)*D3Q19<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_7) + F1*(_T1_u1+cuInvCs2_7) + F2*(_u2));
Fi[9] = latset::w<D3Q19<T>>(9)*D3Q19<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_9) + F1*(_u1) + F2*(T1_u2+cuInvCs2_9));
Fi[10] = latset::w<D3Q19<T>>(10)*D3Q19<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_9) + F1*(_u1) + F2*(_T1_u2+cuInvCs2_9));
Fi[11] = latset::w<D3Q19<T>>(11)*D3Q19<T>::InvCs2 * (F0*(_u0) + F1*(T1_u1+cuInvCs2_11) + F2*(T1_u2+cuInvCs2_11));
Fi[12] = latset::w<D3Q19<T>>(12)*D3Q19<T>::InvCs2 * (F0*(_u0) + F1*(_T1_u1+cuInvCs2_11) + F2*(_T1_u2+cuInvCs2_11));
Fi[13] = latset::w<D3Q19<T>>(13)*D3Q19<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_13) + F1*(_T1_u1-cuInvCs2_13) + F2*(_u2));
Fi[14] = latset::w<D3Q19<T>>(14)*D3Q19<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_13) + F1*(T1_u1-cuInvCs2_13) + F2*(_u2));
Fi[15] = latset::w<D3Q19<T>>(15)*D3Q19<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_15) + F1*(_u1) + F2*(_T1_u2-cuInvCs2_15));
Fi[16] = latset::w<D3Q19<T>>(16)*D3Q19<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_15) + F1*(_u1) + F2*(T1_u2-cuInvCs2_15));
Fi[17] = latset::w<D3Q19<T>>(17)*D3Q19<T>::InvCs2 * (F0*(_u0) + F1*(T1_u1+cuInvCs2_17) + F2*(_T1_u2-cuInvCs2_17));
Fi[18] = latset::w<D3Q19<T>>(18)*D3Q19<T>::InvCs2 * (F0*(_u0) + F1*(_T1_u1+cuInvCs2_17) + F2*(T1_u2-cuInvCs2_17));
}
};

template <typename T>
struct ForcePopImpl<T, D3Q27<T>>{
__any__ static inline void compute(std::array<T, 27> &Fi, const Vector<T, 3> &u, const Vector<T, 3> &F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T F0 = F[0];
const T F1 = F[1];
const T F2 = F[2];
const T T1_u0 = T{1} - u0;
const T _T1_u0 = T{-1} - u0;
const T _u0 = -u0;
const T T1_u1 = T{1} - u1;
const T _T1_u1 = T{-1} - u1;
const T _u1 = -u1;
const T T1_u2 = T{1} - u2;
const T _T1_u2 = T{-1} - u2;
const T _u2 = -u2;
const T cuInvCs2_1 = D3Q27<T>::InvCs2 * (u0);
const T cuInvCs2_3 = D3Q27<T>::InvCs2 * (u1);
const T cuInvCs2_5 = D3Q27<T>::InvCs2 * (u2);
const T cuInvCs2_7 = D3Q27<T>::InvCs2 * (u0+u1);
const T cuInvCs2_9 = D3Q27<T>::InvCs2 * (u0+u2);
const T cuInvCs2_11 = D3Q27<T>::InvCs2 * (u1+u2);
const T cuInvCs2_13 = D3Q27<T>::InvCs2 * (u0-u1);
const T cuInvCs2_15 = D3Q27<T>::InvCs2 * (u0-u2);
const T cuInvCs2_17 = D3Q27<T>::InvCs2 * (u1-u2);
const T cuInvCs2_19 = D3Q27<T>::InvCs2 * (u0+u1+u2);
const T cuInvCs2_21 = D3Q27<T>::InvCs2 * (u0+u1-u2);
const T cuInvCs2_23 = D3Q27<T>::InvCs2 * (u0-u1+u2);
const T cuInvCs2_25 = D3Q27<T>::InvCs2 * (-u0+u1+u2);
Fi[0] = latset::w<D3Q27<T>>(0)*D3Q27<T>::InvCs2 * (-(F0*u0+F1*u1+F2*u2));
Fi[1] = latset::w<D3Q27<T>>(1)*D3Q27<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_1) + F1*(_u1) + F2*(_u2));
Fi[2] = latset::w<D3Q27<T>>(2)*D3Q27<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_1) + F1*(_u1) + F2*(_u2));
Fi[3] = latset::w<D3Q27<T>>(3)*D3Q27<T>::InvCs2 * (F0*(_u0) + F1*(T1_u1+cuInvCs2_3) + F2*(_u2));
Fi[4] = latset::w<D3Q27<T>>(4)*D3Q27<T>::InvCs2 * (F0*(_u0) + F1*(_T1_u1+cuInvCs2_3) + F2*(_u2));
Fi[5] = latset::w<D3Q27<T>>(5)*D3Q27<T>::InvCs2 * (F0*(_u0) + F1*(_u1) + F2*(T1_u2+cuInvCs2_5));
Fi[6] = latset::w<D3Q27<T>>(6)*D3Q27<T>::InvCs2 * (F0*(_u0) + F1*(_u1) + F2*(_T1_u2+cuInvCs2_5));
Fi[7] = latset::w<D3Q27<T>>(7)*D3Q27<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_7) + F1*(T1_u1+cuInvCs2_7) + F2*(_u2));
Fi[8] = latset::w<D3Q27<T>>(8)*D3Q27<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_7) + F1*(_T1_u1+cuInvCs2_7) + F2*(_u2));
Fi[9] = latset::w<D3Q27<T>>(9)*D3Q27<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_9) + F1*(_u1) + F2*(T1_u2+cuInvCs2_9));
Fi[10] = latset::w<D3Q27<T>>(10)*D3Q27<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_9) + F1*(_u1) + F2*(_T1_u2+cuInvCs2_9));
Fi[11] = latset::w<D3Q27<T>>(11)*D3Q27<T>::InvCs2 * (F0*(_u0) + F1*(T1_u1+cuInvCs2_11) + F2*(T1_u2+cuInvCs2_11));
Fi[12] = latset::w<D3Q27<T>>(12)*D3Q27<T>::InvCs2 * (F0*(_u0) + F1*(_T1_u1+cuInvCs2_11) + F2*(_T1_u2+cuInvCs2_11));
Fi[13] = latset::w<D3Q27<T>>(13)*D3Q27<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_13) + F1*(_T1_u1-cuInvCs2_13) + F2*(_u2));
Fi[14] = latset::w<D3Q27<T>>(14)*D3Q27<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_13) + F1*(T1_u1-cuInvCs2_13) + F2*(_u2));
Fi[15] = latset::w<D3Q27<T>>(15)*D3Q27<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_15) + F1*(_u1) + F2*(_T1_u2-cuInvCs2_15));
Fi[16] = latset::w<D3Q27<T>>(16)*D3Q27<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_15) + F1*(_u1) + F2*(T1_u2-cuInvCs2_15));
Fi[17] = latset::w<D3Q27<T>>(17)*D3Q27<T>::InvCs2 * (F0*(_u0) + F1*(T1_u1+cuInvCs2_17) + F2*(_T1_u2-cuInvCs2_17));
Fi[18] = latset::w<D3Q27<T>>(18)*D3Q27<T>::InvCs2 * (F0*(_u0) + F1*(_T1_u1+cuInvCs2_17) + F2*(T1_u2-cuInvCs2_17));
Fi[19] = latset::w<D3Q27<T>>(19)*D3Q27<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_19) + F1*(T1_u1+cuInvCs2_19) + F2*(T1_u2+cuInvCs2_19));
Fi[20] = latset::w<D3Q27<T>>(20)*D3Q27<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_19) + F1*(_T1_u1+cuInvCs2_19) + F2*(_T1_u2+cuInvCs2_19));
Fi[21] = latset::w<D3Q27<T>>(21)*D3Q27<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_21) + F1*(T1_u1+cuInvCs2_21) + F2*(_T1_u2-cuInvCs2_21));
Fi[22] = latset::w<D3Q27<T>>(22)*D3Q27<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_21) + F1*(_T1_u1+cuInvCs2_21) + F2*(T1_u2-cuInvCs2_21));
Fi[23] = latset::w<D3Q27<T>>(23)*D3Q27<T>::InvCs2 * (F0*(T1_u0+cuInvCs2_23) + F1*(_T1_u1-cuInvCs2_23) + F2*(T1_u2+cuInvCs2_23));
Fi[24] = latset::w<D3Q27<T>>(24)*D3Q27<T>::InvCs2 * (F0*(_T1_u0+cuInvCs2_23) + F1*(T1_u1-cuInvCs2_23) + F2*(_T1_u2+cuInvCs2_23));
Fi[25] = latset::w<D3Q27<T>>(25)*D3Q27<T>::InvCs2 * (F0*(_T1_u0-cuInvCs2_25) + F1*(T1_u1+cuInvCs2_25) + F2*(T1_u2+cuInvCs2_25));
Fi[26] = latset::w<D3Q27<T>>(26)*D3Q27<T>::InvCs2 * (F0*(T1_u0-cuInvCs2_25) + F1*(_T1_u1+cuInvCs2_25) + F2*(_T1_u2+cuInvCs2_25));
}
};


//------------------------------------


//------------------------------------

template <typename T>
struct ScalarForcePopImpl<T, D2Q5<T>, 0>{
__any__ static inline void compute(std::array<T, 5> &Fi, const Vector<T, 2> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T T1_ud = T{1} - u0;
const T _T1_ud = T{-1} - u0;
const T _ud = -u0;
const T cuInvCs2_1 = D2Q5<T>::InvCs2 * (u0);
Fi[0] = latset::w<D2Q5<T>>(0) * D2Q5<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D2Q5<T>>(1) * D2Q5<T>::InvCs2 * F * (T1_ud+cuInvCs2_1);
Fi[2] = latset::w<D2Q5<T>>(2) * D2Q5<T>::InvCs2 * F * (_T1_ud+cuInvCs2_1);
Fi[3] = latset::w<D2Q5<T>>(3) * D2Q5<T>::InvCs2 * F * (_ud);
Fi[4] = latset::w<D2Q5<T>>(4) * D2Q5<T>::InvCs2 * F * (_ud);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D2Q5<T>, 1>{
__any__ static inline void compute(std::array<T, 5> &Fi, const Vector<T, 2> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T T1_ud = T{1} - u1;
const T _T1_ud = T{-1} - u1;
const T _ud = -u1;
const T cuInvCs2_3 = D2Q5<T>::InvCs2 * (u1);
Fi[0] = latset::w<D2Q5<T>>(0) * D2Q5<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D2Q5<T>>(1) * D2Q5<T>::InvCs2 * F * (_ud);
Fi[2] = latset::w<D2Q5<T>>(2) * D2Q5<T>::InvCs2 * F * (_ud);
Fi[3] = latset::w<D2Q5<T>>(3) * D2Q5<T>::InvCs2 * F * (T1_ud+cuInvCs2_3);
Fi[4] = latset::w<D2Q5<T>>(4) * D2Q5<T>::InvCs2 * F * (_T1_ud+cuInvCs2_3);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D2Q9<T>, 0>{
__any__ static inline void compute(std::array<T, 9> &Fi, const Vector<T, 2> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T T1_ud = T{1} - u0;
const T _T1_ud = T{-1} - u0;
const T _ud = -u0;
const T cuInvCs2_1 = D2Q9<T>::InvCs2 * (u0);
const T cuInvCs2_5 = D2Q9<T>::InvCs2 * (u0+u1);
const T cuInvCs2_7 = D2Q9<T>::InvCs2 * (u0-u1);
Fi[0] = latset::w<D2Q9<T>>(0) * D2Q9<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D2Q9<T>>(1) * D2Q9<T>::InvCs2 * F * (T1_ud+cuInvCs2_1);
Fi[2] = latset::w<D2Q9<T>>(2) * D2Q9<T>::InvCs2 * F * (_T1_ud+cuInvCs2_1);
Fi[3] = latset::w<D2Q9<T>>(3) * D2Q9<T>::InvCs2 * F * (_ud);
Fi[4] = latset::w<D2Q9<T>>(4) * D2Q9<T>::InvCs2 * F * (_ud);
Fi[5] = latset::w<D2Q9<T>>(5) * D2Q9<T>::InvCs2 * F * (T1_ud+cuInvCs2_5);
Fi[6] = latset::w<D2Q9<T>>(6) * D2Q9<T>::InvCs2 * F * (_T1_ud+cuInvCs2_5);
Fi[7] = latset::w<D2Q9<T>>(7) * D2Q9<T>::InvCs2 * F * (T1_ud+cuInvCs2_7);
Fi[8] = latset::w<D2Q9<T>>(8) * D2Q9<T>::InvCs2 * F * (_T1_ud+cuInvCs2_7);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D2Q9<T>, 1>{
__any__ static inline void compute(std::array<T, 9> &Fi, const Vector<T, 2> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T T1_ud = T{1} - u1;
const T _T1_ud = T{-1} - u1;
const T _ud = -u1;
const T cuInvCs2_3 = D2Q9<T>::InvCs2 * (u1);
const T cuInvCs2_5 = D2Q9<T>::InvCs2 * (u0+u1);
const T cuInvCs2_7 = D2Q9<T>::InvCs2 * (u0-u1);
Fi[0] = latset::w<D2Q9<T>>(0) * D2Q9<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D2Q9<T>>(1) * D2Q9<T>::InvCs2 * F * (_ud);
Fi[2] = latset::w<D2Q9<T>>(2) * D2Q9<T>::InvCs2 * F * (_ud);
Fi[3] = latset::w<D2Q9<T>>(3) * D2Q9<T>::InvCs2 * F * (T1_ud+cuInvCs2_3);
Fi[4] = latset::w<D2Q9<T>>(4) * D2Q9<T>::InvCs2 * F * (_T1_ud+cuInvCs2_3);
Fi[5] = latset::w<D2Q9<T>>(5) * D2Q9<T>::InvCs2 * F * (T1_ud+cuInvCs2_5);
Fi[6] = latset::w<D2Q9<T>>(6) * D2Q9<T>::InvCs2 * F * (_T1_ud+cuInvCs2_5);
Fi[7] = latset::w<D2Q9<T>>(7) * D2Q9<T>::InvCs2 * F * (_T1_ud-cuInvCs2_7);
Fi[8] = latset::w<D2Q9<T>>(8) * D2Q9<T>::InvCs2 * F * (T1_ud-cuInvCs2_7);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q7<T>, 0>{
__any__ static inline void compute(std::array<T, 7> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u0;
const T _T1_ud = T{-1} - u0;
const T _ud = -u0;
const T cuInvCs2_1 = D3Q7<T>::InvCs2 * (u0);
Fi[0] = latset::w<D3Q7<T>>(0) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q7<T>>(1) * D3Q7<T>::InvCs2 * F * (T1_ud+cuInvCs2_1);
Fi[2] = latset::w<D3Q7<T>>(2) * D3Q7<T>::InvCs2 * F * (_T1_ud+cuInvCs2_1);
Fi[3] = latset::w<D3Q7<T>>(3) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[4] = latset::w<D3Q7<T>>(4) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[5] = latset::w<D3Q7<T>>(5) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[6] = latset::w<D3Q7<T>>(6) * D3Q7<T>::InvCs2 * F * (_ud);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q7<T>, 1>{
__any__ static inline void compute(std::array<T, 7> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u1;
const T _T1_ud = T{-1} - u1;
const T _ud = -u1;
const T cuInvCs2_3 = D3Q7<T>::InvCs2 * (u1);
Fi[0] = latset::w<D3Q7<T>>(0) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q7<T>>(1) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[2] = latset::w<D3Q7<T>>(2) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[3] = latset::w<D3Q7<T>>(3) * D3Q7<T>::InvCs2 * F * (T1_ud+cuInvCs2_3);
Fi[4] = latset::w<D3Q7<T>>(4) * D3Q7<T>::InvCs2 * F * (_T1_ud+cuInvCs2_3);
Fi[5] = latset::w<D3Q7<T>>(5) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[6] = latset::w<D3Q7<T>>(6) * D3Q7<T>::InvCs2 * F * (_ud);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q7<T>, 2>{
__any__ static inline void compute(std::array<T, 7> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u2;
const T _T1_ud = T{-1} - u2;
const T _ud = -u2;
const T cuInvCs2_5 = D3Q7<T>::InvCs2 * (u2);
Fi[0] = latset::w<D3Q7<T>>(0) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q7<T>>(1) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[2] = latset::w<D3Q7<T>>(2) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[3] = latset::w<D3Q7<T>>(3) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[4] = latset::w<D3Q7<T>>(4) * D3Q7<T>::InvCs2 * F * (_ud);
Fi[5] = latset::w<D3Q7<T>>(5) * D3Q7<T>::InvCs2 * F * (T1_ud+cuInvCs2_5);
Fi[6] = latset::w<D3Q7<T>>(6) * D3Q7<T>::InvCs2 * F * (_T1_ud+cuInvCs2_5);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q15<T>, 0>{
__any__ static inline void compute(std::array<T, 15> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u0;
const T _T1_ud = T{-1} - u0;
const T _ud = -u0;
const T cuInvCs2_1 = D3Q15<T>::InvCs2 * (u0);
const T cuInvCs2_7 = D3Q15<T>::InvCs2 * (u0+u1+u2);
const T cuInvCs2_9 = D3Q15<T>::InvCs2 * (u0+u1-u2);
const T cuInvCs2_11 = D3Q15<T>::InvCs2 * (u0-u1+u2);
const T cuInvCs2_13 = D3Q15<T>::InvCs2 * (-u0+u1+u2);
Fi[0] = latset::w<D3Q15<T>>(0) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q15<T>>(1) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_1);
Fi[2] = latset::w<D3Q15<T>>(2) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_1);
Fi[3] = latset::w<D3Q15<T>>(3) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[4] = latset::w<D3Q15<T>>(4) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[5] = latset::w<D3Q15<T>>(5) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[6] = latset::w<D3Q15<T>>(6) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[7] = latset::w<D3Q15<T>>(7) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_7);
Fi[8] = latset::w<D3Q15<T>>(8) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_7);
Fi[9] = latset::w<D3Q15<T>>(9) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_9);
Fi[10] = latset::w<D3Q15<T>>(10) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_9);
Fi[11] = latset::w<D3Q15<T>>(11) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_11);
Fi[12] = latset::w<D3Q15<T>>(12) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_11);
Fi[13] = latset::w<D3Q15<T>>(13) * D3Q15<T>::InvCs2 * F * (_T1_ud-cuInvCs2_13);
Fi[14] = latset::w<D3Q15<T>>(14) * D3Q15<T>::InvCs2 * F * (T1_ud-cuInvCs2_13);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q15<T>, 1>{
__any__ static inline void compute(std::array<T, 15> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u1;
const T _T1_ud = T{-1} - u1;
const T _ud = -u1;
const T cuInvCs2_3 = D3Q15<T>::InvCs2 * (u1);
const T cuInvCs2_7 = D3Q15<T>::InvCs2 * (u0+u1+u2);
const T cuInvCs2_9 = D3Q15<T>::InvCs2 * (u0+u1-u2);
const T cuInvCs2_11 = D3Q15<T>::InvCs2 * (u0-u1+u2);
const T cuInvCs2_13 = D3Q15<T>::InvCs2 * (-u0+u1+u2);
Fi[0] = latset::w<D3Q15<T>>(0) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q15<T>>(1) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[2] = latset::w<D3Q15<T>>(2) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[3] = latset::w<D3Q15<T>>(3) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_3);
Fi[4] = latset::w<D3Q15<T>>(4) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_3);
Fi[5] = latset::w<D3Q15<T>>(5) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[6] = latset::w<D3Q15<T>>(6) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[7] = latset::w<D3Q15<T>>(7) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_7);
Fi[8] = latset::w<D3Q15<T>>(8) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_7);
Fi[9] = latset::w<D3Q15<T>>(9) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_9);
Fi[10] = latset::w<D3Q15<T>>(10) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_9);
Fi[11] = latset::w<D3Q15<T>>(11) * D3Q15<T>::InvCs2 * F * (_T1_ud-cuInvCs2_11);
Fi[12] = latset::w<D3Q15<T>>(12) * D3Q15<T>::InvCs2 * F * (T1_ud-cuInvCs2_11);
Fi[13] = latset::w<D3Q15<T>>(13) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_13);
Fi[14] = latset::w<D3Q15<T>>(14) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_13);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q15<T>, 2>{
__any__ static inline void compute(std::array<T, 15> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u2;
const T _T1_ud = T{-1} - u2;
const T _ud = -u2;
const T cuInvCs2_5 = D3Q15<T>::InvCs2 * (u2);
const T cuInvCs2_7 = D3Q15<T>::InvCs2 * (u0+u1+u2);
const T cuInvCs2_9 = D3Q15<T>::InvCs2 * (u0+u1-u2);
const T cuInvCs2_11 = D3Q15<T>::InvCs2 * (u0-u1+u2);
const T cuInvCs2_13 = D3Q15<T>::InvCs2 * (-u0+u1+u2);
Fi[0] = latset::w<D3Q15<T>>(0) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q15<T>>(1) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[2] = latset::w<D3Q15<T>>(2) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[3] = latset::w<D3Q15<T>>(3) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[4] = latset::w<D3Q15<T>>(4) * D3Q15<T>::InvCs2 * F * (_ud);
Fi[5] = latset::w<D3Q15<T>>(5) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_5);
Fi[6] = latset::w<D3Q15<T>>(6) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_5);
Fi[7] = latset::w<D3Q15<T>>(7) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_7);
Fi[8] = latset::w<D3Q15<T>>(8) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_7);
Fi[9] = latset::w<D3Q15<T>>(9) * D3Q15<T>::InvCs2 * F * (_T1_ud-cuInvCs2_9);
Fi[10] = latset::w<D3Q15<T>>(10) * D3Q15<T>::InvCs2 * F * (T1_ud-cuInvCs2_9);
Fi[11] = latset::w<D3Q15<T>>(11) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_11);
Fi[12] = latset::w<D3Q15<T>>(12) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_11);
Fi[13] = latset::w<D3Q15<T>>(13) * D3Q15<T>::InvCs2 * F * (T1_ud+cuInvCs2_13);
Fi[14] = latset::w<D3Q15<T>>(14) * D3Q15<T>::InvCs2 * F * (_T1_ud+cuInvCs2_13);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q19<T>, 0>{
__any__ static inline void compute(std::array<T, 19> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u0;
const T _T1_ud = T{-1} - u0;
const T _ud = -u0;
const T cuInvCs2_1 = D3Q19<T>::InvCs2 * (u0);
const T cuInvCs2_7 = D3Q19<T>::InvCs2 * (u0+u1);
const T cuInvCs2_9 = D3Q19<T>::InvCs2 * (u0+u2);
const T cuInvCs2_13 = D3Q19<T>::InvCs2 * (u0-u1);
const T cuInvCs2_15 = D3Q19<T>::InvCs2 * (u0-u2);
Fi[0] = latset::w<D3Q19<T>>(0) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q19<T>>(1) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_1);
Fi[2] = latset::w<D3Q19<T>>(2) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_1);
Fi[3] = latset::w<D3Q19<T>>(3) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[4] = latset::w<D3Q19<T>>(4) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[5] = latset::w<D3Q19<T>>(5) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[6] = latset::w<D3Q19<T>>(6) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[7] = latset::w<D3Q19<T>>(7) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_7);
Fi[8] = latset::w<D3Q19<T>>(8) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_7);
Fi[9] = latset::w<D3Q19<T>>(9) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_9);
Fi[10] = latset::w<D3Q19<T>>(10) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_9);
Fi[11] = latset::w<D3Q19<T>>(11) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[12] = latset::w<D3Q19<T>>(12) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[13] = latset::w<D3Q19<T>>(13) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_13);
Fi[14] = latset::w<D3Q19<T>>(14) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_13);
Fi[15] = latset::w<D3Q19<T>>(15) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_15);
Fi[16] = latset::w<D3Q19<T>>(16) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_15);
Fi[17] = latset::w<D3Q19<T>>(17) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[18] = latset::w<D3Q19<T>>(18) * D3Q19<T>::InvCs2 * F * (_ud);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q19<T>, 1>{
__any__ static inline void compute(std::array<T, 19> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u1;
const T _T1_ud = T{-1} - u1;
const T _ud = -u1;
const T cuInvCs2_3 = D3Q19<T>::InvCs2 * (u1);
const T cuInvCs2_7 = D3Q19<T>::InvCs2 * (u0+u1);
const T cuInvCs2_11 = D3Q19<T>::InvCs2 * (u1+u2);
const T cuInvCs2_13 = D3Q19<T>::InvCs2 * (u0-u1);
const T cuInvCs2_17 = D3Q19<T>::InvCs2 * (u1-u2);
Fi[0] = latset::w<D3Q19<T>>(0) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q19<T>>(1) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[2] = latset::w<D3Q19<T>>(2) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[3] = latset::w<D3Q19<T>>(3) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_3);
Fi[4] = latset::w<D3Q19<T>>(4) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_3);
Fi[5] = latset::w<D3Q19<T>>(5) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[6] = latset::w<D3Q19<T>>(6) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[7] = latset::w<D3Q19<T>>(7) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_7);
Fi[8] = latset::w<D3Q19<T>>(8) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_7);
Fi[9] = latset::w<D3Q19<T>>(9) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[10] = latset::w<D3Q19<T>>(10) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[11] = latset::w<D3Q19<T>>(11) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_11);
Fi[12] = latset::w<D3Q19<T>>(12) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_11);
Fi[13] = latset::w<D3Q19<T>>(13) * D3Q19<T>::InvCs2 * F * (_T1_ud-cuInvCs2_13);
Fi[14] = latset::w<D3Q19<T>>(14) * D3Q19<T>::InvCs2 * F * (T1_ud-cuInvCs2_13);
Fi[15] = latset::w<D3Q19<T>>(15) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[16] = latset::w<D3Q19<T>>(16) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[17] = latset::w<D3Q19<T>>(17) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_17);
Fi[18] = latset::w<D3Q19<T>>(18) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_17);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q19<T>, 2>{
__any__ static inline void compute(std::array<T, 19> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u2;
const T _T1_ud = T{-1} - u2;
const T _ud = -u2;
const T cuInvCs2_5 = D3Q19<T>::InvCs2 * (u2);
const T cuInvCs2_9 = D3Q19<T>::InvCs2 * (u0+u2);
const T cuInvCs2_11 = D3Q19<T>::InvCs2 * (u1+u2);
const T cuInvCs2_15 = D3Q19<T>::InvCs2 * (u0-u2);
const T cuInvCs2_17 = D3Q19<T>::InvCs2 * (u1-u2);
Fi[0] = latset::w<D3Q19<T>>(0) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q19<T>>(1) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[2] = latset::w<D3Q19<T>>(2) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[3] = latset::w<D3Q19<T>>(3) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[4] = latset::w<D3Q19<T>>(4) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[5] = latset::w<D3Q19<T>>(5) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_5);
Fi[6] = latset::w<D3Q19<T>>(6) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_5);
Fi[7] = latset::w<D3Q19<T>>(7) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[8] = latset::w<D3Q19<T>>(8) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[9] = latset::w<D3Q19<T>>(9) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_9);
Fi[10] = latset::w<D3Q19<T>>(10) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_9);
Fi[11] = latset::w<D3Q19<T>>(11) * D3Q19<T>::InvCs2 * F * (T1_ud+cuInvCs2_11);
Fi[12] = latset::w<D3Q19<T>>(12) * D3Q19<T>::InvCs2 * F * (_T1_ud+cuInvCs2_11);
Fi[13] = latset::w<D3Q19<T>>(13) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[14] = latset::w<D3Q19<T>>(14) * D3Q19<T>::InvCs2 * F * (_ud);
Fi[15] = latset::w<D3Q19<T>>(15) * D3Q19<T>::InvCs2 * F * (_T1_ud-cuInvCs2_15);
Fi[16] = latset::w<D3Q19<T>>(16) * D3Q19<T>::InvCs2 * F * (T1_ud-cuInvCs2_15);
Fi[17] = latset::w<D3Q19<T>>(17) * D3Q19<T>::InvCs2 * F * (_T1_ud-cuInvCs2_17);
Fi[18] = latset::w<D3Q19<T>>(18) * D3Q19<T>::InvCs2 * F * (T1_ud-cuInvCs2_17);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q27<T>, 0>{
__any__ static inline void compute(std::array<T, 27> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u0;
const T _T1_ud = T{-1} - u0;
const T _ud = -u0;
const T cuInvCs2_1 = D3Q27<T>::InvCs2 * (u0);
const T cuInvCs2_7 = D3Q27<T>::InvCs2 * (u0+u1);
const T cuInvCs2_9 = D3Q27<T>::InvCs2 * (u0+u2);
const T cuInvCs2_13 = D3Q27<T>::InvCs2 * (u0-u1);
const T cuInvCs2_15 = D3Q27<T>::InvCs2 * (u0-u2);
const T cuInvCs2_19 = D3Q27<T>::InvCs2 * (u0+u1+u2);
const T cuInvCs2_21 = D3Q27<T>::InvCs2 * (u0+u1-u2);
const T cuInvCs2_23 = D3Q27<T>::InvCs2 * (u0-u1+u2);
const T cuInvCs2_25 = D3Q27<T>::InvCs2 * (-u0+u1+u2);
Fi[0] = latset::w<D3Q27<T>>(0) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q27<T>>(1) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_1);
Fi[2] = latset::w<D3Q27<T>>(2) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_1);
Fi[3] = latset::w<D3Q27<T>>(3) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[4] = latset::w<D3Q27<T>>(4) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[5] = latset::w<D3Q27<T>>(5) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[6] = latset::w<D3Q27<T>>(6) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[7] = latset::w<D3Q27<T>>(7) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_7);
Fi[8] = latset::w<D3Q27<T>>(8) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_7);
Fi[9] = latset::w<D3Q27<T>>(9) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_9);
Fi[10] = latset::w<D3Q27<T>>(10) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_9);
Fi[11] = latset::w<D3Q27<T>>(11) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[12] = latset::w<D3Q27<T>>(12) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[13] = latset::w<D3Q27<T>>(13) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_13);
Fi[14] = latset::w<D3Q27<T>>(14) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_13);
Fi[15] = latset::w<D3Q27<T>>(15) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_15);
Fi[16] = latset::w<D3Q27<T>>(16) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_15);
Fi[17] = latset::w<D3Q27<T>>(17) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[18] = latset::w<D3Q27<T>>(18) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[19] = latset::w<D3Q27<T>>(19) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_19);
Fi[20] = latset::w<D3Q27<T>>(20) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_19);
Fi[21] = latset::w<D3Q27<T>>(21) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_21);
Fi[22] = latset::w<D3Q27<T>>(22) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_21);
Fi[23] = latset::w<D3Q27<T>>(23) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_23);
Fi[24] = latset::w<D3Q27<T>>(24) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_23);
Fi[25] = latset::w<D3Q27<T>>(25) * D3Q27<T>::InvCs2 * F * (_T1_ud-cuInvCs2_25);
Fi[26] = latset::w<D3Q27<T>>(26) * D3Q27<T>::InvCs2 * F * (T1_ud-cuInvCs2_25);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q27<T>, 1>{
__any__ static inline void compute(std::array<T, 27> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u1;
const T _T1_ud = T{-1} - u1;
const T _ud = -u1;
const T cuInvCs2_3 = D3Q27<T>::InvCs2 * (u1);
const T cuInvCs2_7 = D3Q27<T>::InvCs2 * (u0+u1);
const T cuInvCs2_11 = D3Q27<T>::InvCs2 * (u1+u2);
const T cuInvCs2_13 = D3Q27<T>::InvCs2 * (u0-u1);
const T cuInvCs2_17 = D3Q27<T>::InvCs2 * (u1-u2);
const T cuInvCs2_19 = D3Q27<T>::InvCs2 * (u0+u1+u2);
const T cuInvCs2_21 = D3Q27<T>::InvCs2 * (u0+u1-u2);
const T cuInvCs2_23 = D3Q27<T>::InvCs2 * (u0-u1+u2);
const T cuInvCs2_25 = D3Q27<T>::InvCs2 * (-u0+u1+u2);
Fi[0] = latset::w<D3Q27<T>>(0) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q27<T>>(1) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[2] = latset::w<D3Q27<T>>(2) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[3] = latset::w<D3Q27<T>>(3) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_3);
Fi[4] = latset::w<D3Q27<T>>(4) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_3);
Fi[5] = latset::w<D3Q27<T>>(5) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[6] = latset::w<D3Q27<T>>(6) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[7] = latset::w<D3Q27<T>>(7) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_7);
Fi[8] = latset::w<D3Q27<T>>(8) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_7);
Fi[9] = latset::w<D3Q27<T>>(9) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[10] = latset::w<D3Q27<T>>(10) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[11] = latset::w<D3Q27<T>>(11) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_11);
Fi[12] = latset::w<D3Q27<T>>(12) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_11);
Fi[13] = latset::w<D3Q27<T>>(13) * D3Q27<T>::InvCs2 * F * (_T1_ud-cuInvCs2_13);
Fi[14] = latset::w<D3Q27<T>>(14) * D3Q27<T>::InvCs2 * F * (T1_ud-cuInvCs2_13);
Fi[15] = latset::w<D3Q27<T>>(15) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[16] = latset::w<D3Q27<T>>(16) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[17] = latset::w<D3Q27<T>>(17) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_17);
Fi[18] = latset::w<D3Q27<T>>(18) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_17);
Fi[19] = latset::w<D3Q27<T>>(19) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_19);
Fi[20] = latset::w<D3Q27<T>>(20) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_19);
Fi[21] = latset::w<D3Q27<T>>(21) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_21);
Fi[22] = latset::w<D3Q27<T>>(22) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_21);
Fi[23] = latset::w<D3Q27<T>>(23) * D3Q27<T>::InvCs2 * F * (_T1_ud-cuInvCs2_23);
Fi[24] = latset::w<D3Q27<T>>(24) * D3Q27<T>::InvCs2 * F * (T1_ud-cuInvCs2_23);
Fi[25] = latset::w<D3Q27<T>>(25) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_25);
Fi[26] = latset::w<D3Q27<T>>(26) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_25);
}
};

template <typename T>
struct ScalarForcePopImpl<T, D3Q27<T>, 2>{
__any__ static inline void compute(std::array<T, 27> &Fi, const Vector<T, 3> &u, const T F){
const T u0 = u[0];
const T u1 = u[1];
const T u2 = u[2];
const T T1_ud = T{1} - u2;
const T _T1_ud = T{-1} - u2;
const T _ud = -u2;
const T cuInvCs2_5 = D3Q27<T>::InvCs2 * (u2);
const T cuInvCs2_9 = D3Q27<T>::InvCs2 * (u0+u2);
const T cuInvCs2_11 = D3Q27<T>::InvCs2 * (u1+u2);
const T cuInvCs2_15 = D3Q27<T>::InvCs2 * (u0-u2);
const T cuInvCs2_17 = D3Q27<T>::InvCs2 * (u1-u2);
const T cuInvCs2_19 = D3Q27<T>::InvCs2 * (u0+u1+u2);
const T cuInvCs2_21 = D3Q27<T>::InvCs2 * (u0+u1-u2);
const T cuInvCs2_23 = D3Q27<T>::InvCs2 * (u0-u1+u2);
const T cuInvCs2_25 = D3Q27<T>::InvCs2 * (-u0+u1+u2);
Fi[0] = latset::w<D3Q27<T>>(0) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[1] = latset::w<D3Q27<T>>(1) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[2] = latset::w<D3Q27<T>>(2) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[3] = latset::w<D3Q27<T>>(3) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[4] = latset::w<D3Q27<T>>(4) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[5] = latset::w<D3Q27<T>>(5) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_5);
Fi[6] = latset::w<D3Q27<T>>(6) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_5);
Fi[7] = latset::w<D3Q27<T>>(7) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[8] = latset::w<D3Q27<T>>(8) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[9] = latset::w<D3Q27<T>>(9) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_9);
Fi[10] = latset::w<D3Q27<T>>(10) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_9);
Fi[11] = latset::w<D3Q27<T>>(11) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_11);
Fi[12] = latset::w<D3Q27<T>>(12) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_11);
Fi[13] = latset::w<D3Q27<T>>(13) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[14] = latset::w<D3Q27<T>>(14) * D3Q27<T>::InvCs2 * F * (_ud);
Fi[15] = latset::w<D3Q27<T>>(15) * D3Q27<T>::InvCs2 * F * (_T1_ud-cuInvCs2_15);
Fi[16] = latset::w<D3Q27<T>>(16) * D3Q27<T>::InvCs2 * F * (T1_ud-cuInvCs2_15);
Fi[17] = latset::w<D3Q27<T>>(17) * D3Q27<T>::InvCs2 * F * (_T1_ud-cuInvCs2_17);
Fi[18] = latset::w<D3Q27<T>>(18) * D3Q27<T>::InvCs2 * F * (T1_ud-cuInvCs2_17);
Fi[19] = latset::w<D3Q27<T>>(19) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_19);
Fi[20] = latset::w<D3Q27<T>>(20) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_19);
Fi[21] = latset::w<D3Q27<T>>(21) * D3Q27<T>::InvCs2 * F * (_T1_ud-cuInvCs2_21);
Fi[22] = latset::w<D3Q27<T>>(22) * D3Q27<T>::InvCs2 * F * (T1_ud-cuInvCs2_21);
Fi[23] = latset::w<D3Q27<T>>(23) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_23);
Fi[24] = latset::w<D3Q27<T>>(24) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_23);
Fi[25] = latset::w<D3Q27<T>>(25) * D3Q27<T>::InvCs2 * F * (T1_ud+cuInvCs2_25);
Fi[26] = latset::w<D3Q27<T>>(26) * D3Q27<T>::InvCs2 * F * (_T1_ud+cuInvCs2_25);
}
};


//------------------------------------



}  // namespace force

#endif