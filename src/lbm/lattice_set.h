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

// lattice set for lbm and ca
// lattice_set.h

#pragma once

#include "data_struct/Vector.h"
#include "utils/util.h"


// lattice set variable template s
namespace latsetdata {
// lattice set for discrete velocity
template <unsigned int D, unsigned int Q>
__constexpr__ Vector<int, D> c[Q] = {};

// lattice set for weight
// we have to use Fraction here not template typename T
// cause partial specialization of variable template is not allowed
template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> w[Q] = {};

// lattice set for opposite direction
template <unsigned int D, unsigned int Q>
__constexpr__ int opp[Q] = {};

}  // namespace latsetdata

namespace latset {

// lattice set functions
template <unsigned int D, unsigned int Q>
constexpr const Vector<int, D>& c(unsigned int i) {
  return latsetdata::c<D, Q>[i];
}

template <typename T, unsigned int D, unsigned int Q>
constexpr T w(unsigned int i) {
  return latsetdata::w<D, Q>[i].template operator()<T>();
}

template <unsigned int D, unsigned int Q>
constexpr int opp(unsigned int i) {
  return latsetdata::opp<D, Q>[i];
}

// lattice set functions using LatSet template
template <typename LatSet>
constexpr const Vector<int, LatSet::d>& c(unsigned int i) {
#ifdef __CUDA_ARCH__
  return c<LatSet::d, LatSet::q>(i);
#else
  return LatSet::c[i];
#endif
}

template <typename LatSet>
constexpr typename LatSet::FloatType w(unsigned int i) {
#ifdef __CUDA_ARCH__
  return w<typename LatSet::FloatType, LatSet::d, LatSet::q>(i);
#else
  return LatSet::w[i];
#endif
}

template <typename LatSet>
constexpr int opp(unsigned int i) {
#ifdef __CUDA_ARCH__
  return opp<LatSet::d, LatSet::q>(i);
#else
  return LatSet::opp[i];
#endif
}

}  // namespace latset


// lattice set, including dimension, velocity, weight, and direction
// basic lattice set, DdQq
// lattice set: DdQq(e.g., D2Q9)
// the member of DdQq struct can be accessed without inheritance,
// e.g., struct D1Q3 : public Basic_Lattice_Set<1, 3>
template <unsigned int D, unsigned int Q>
struct Basic_Lattice_Set {
  // DdQq
  static constexpr unsigned int d = D;
  static constexpr unsigned int q = Q;
};

// D1Q3
template <typename T>
struct D1Q3 : public Basic_Lattice_Set<1, 3> {
  using FloatType = T;
  static constexpr int c[q][d] = {{0}, {1}, {-1}};
  static constexpr T w[q] = {T(1) / T(3), T(1) / T(6), T(1) / T(6)};
  static constexpr int opp[q] = {0, 2, 1};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
};

// D2Q4
template <typename T>
struct D2Q4 : public Basic_Lattice_Set<2, 4> {
  using FloatType = T;
  static constexpr Vector<int, 2> c[q] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
  static constexpr T w[q] = {T(1) / T(4), T(1) / T(4), T(1) / T(4), T(1) / T(4)};
  static constexpr int opp[q] = {1, 0, 3, 2};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
};

namespace latsetdata {
template <>
__constexpr__ Vector<int, 2> c<2, 4>[4] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

template <>
__constexpr__ Fraction<> w<2, 4>[4] = {{1, 4}, {1, 4}, {1, 4}, {1, 4}};

template <>
__constexpr__ int opp<2, 4>[4] = {1, 0, 3, 2};
}  // namespace latsetdata

// D2Q5
template <typename T>
struct D2Q5 : public Basic_Lattice_Set<2, 5> {
  using FloatType = T;
  static constexpr Vector<int, 2> c[q] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}};
  static constexpr T w[q] = {T(1) / T(3), T(1) / T(6), T(1) / T(6), T(1) / T(6),
                             T(1) / T(6)};
  static constexpr int opp[q] = {0, 2, 1, 4, 3};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
};

namespace latsetdata {
template <>
__constexpr__ Vector<int, 2> c<2, 5>[5] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}};

template <>
__constexpr__ Fraction<> w<2, 5>[5] = {{1, 3}, {1, 6}, {1, 6}, {1, 6}, {1, 6}};

template <>
__constexpr__ int opp<2, 5>[5] = {0, 2, 1, 4, 3};
}  // namespace latsetdata

// D2Q9
template <typename T>
struct D2Q9 : public Basic_Lattice_Set<2, 9> {
  using FloatType = T;
  static constexpr Vector<int, 2> c[q] = {{0, 0}, {1, 0},   {-1, 0}, {0, 1}, {0, -1},
                                          {1, 1}, {-1, -1}, {1, -1}, {-1, 1}};
  static constexpr T w[q] = {T(4) / T(9),  T(1) / T(9),  T(1) / T(9),
                             T(1) / T(9),  T(1) / T(9),  T(1) / T(36),
                             T(1) / T(36), T(1) / T(36), T(1) / T(36)};
  static constexpr int opp[q] = {0, 2, 1, 4, 3, 6, 5, 8, 7};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
};

namespace latsetdata {
template <>
__constexpr__ Vector<int, 2> c<2, 9>[9] = {{0, 0}, {1, 0},   {-1, 0}, {0, 1}, {0, -1},
                                           {1, 1}, {-1, -1}, {1, -1}, {-1, 1}};

template <>
__constexpr__ Fraction<> w<2, 9>[9] = {{4, 9},  {1, 9},  {1, 9},  {1, 9}, {1, 9},
                                       {1, 36}, {1, 36}, {1, 36}, {1, 36}};

template <>
__constexpr__ int opp<2, 9>[9] = {0, 2, 1, 4, 3, 6, 5, 8, 7};
}  // namespace latsetdata

// D3Q7
template <typename T>
struct D3Q7 : public Basic_Lattice_Set<3, 7> {
  using FloatType = T;
  static constexpr Vector<int, 3> c[q] = {{0, 0, 0},  {1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
                                          {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
  static constexpr T w[q] = {T(1) / T(4), T(1) / T(8), T(1) / T(8), T(1) / T(8),
                             T(1) / T(8), T(1) / T(8), T(1) / T(8)};
  static constexpr int opp[q] = {0, 2, 1, 4, 3, 6, 5};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
};

namespace latsetdata {
template <>
__constexpr__ Vector<int, 3> c<3, 7>[7] = {{0, 0, 0},  {1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
                                           {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};

template <>
__constexpr__ Fraction<> w<3, 7>[7] = {{1, 4}, {1, 8}, {1, 8}, {1, 8},
                                       {1, 8}, {1, 8}, {1, 8}};

template <>
__constexpr__ int opp<3, 7>[7] = {0, 2, 1, 4, 3, 6, 5};

}  // namespace latsetdata

// D3Q15
template <typename T>
struct D3Q15 : public Basic_Lattice_Set<3, 15> {
  using FloatType = T;
  static constexpr Vector<int, 3> c[q] = {
    {0, 0, 0},   {1, 0, 0},  {-1, 0, 0},  {0, 1, 0},    {0, -1, 0},
    {0, 0, 1},   {0, 0, -1}, {1, 1, 1},   {-1, -1, -1}, {1, 1, -1},
    {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1}, {-1, 1, 1},   {1, -1, -1}};
  static constexpr T w[q] = {T(2) / T(9),  T(1) / T(9),  T(1) / T(9),  T(1) / T(9),
                             T(1) / T(9),  T(1) / T(9),  T(1) / T(9),  T(1) / T(72),
                             T(1) / T(72), T(1) / T(72), T(1) / T(72), T(1) / T(72),
                             T(1) / T(72), T(1) / T(72), T(1) / T(72)};
  static constexpr int opp[q] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
};

namespace latsetdata {
template <>
__constexpr__ Vector<int, 3> c<3, 15>[15] = {
  {0, 0, 0},   {1, 0, 0},  {-1, 0, 0},  {0, 1, 0},    {0, -1, 0},
  {0, 0, 1},   {0, 0, -1}, {1, 1, 1},   {-1, -1, -1}, {1, 1, -1},
  {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1}, {-1, 1, 1},   {1, -1, -1}};

template <>
__constexpr__ Fraction<> w<3, 15>[15] = {{2, 9},  {1, 9},  {1, 9},  {1, 9},  {1, 9},
                                         {1, 9},  {1, 9},  {1, 72}, {1, 72}, {1, 72},
                                         {1, 72}, {1, 72}, {1, 72}, {1, 72}, {1, 72}};

template <>
__constexpr__ int opp<3, 15>[15] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};

}  // namespace latsetdata

// D3Q19
template <typename T>
struct D3Q19 : public Basic_Lattice_Set<3, 19> {
  using FloatType = T;
  static constexpr Vector<int, 3> c[q] = {
    {0, 0, 0},  {1, 0, 0},   {-1, 0, 0}, {0, 1, 0},   {0, -1, 0}, {0, 0, 1},   {0, 0, -1},
    {1, 1, 0},  {-1, -1, 0}, {1, 0, 1},  {-1, 0, -1}, {0, 1, 1},  {0, -1, -1}, {1, -1, 0},
    {-1, 1, 0}, {1, 0, -1},  {-1, 0, 1}, {0, 1, -1},  {0, -1, 1}};
  static constexpr T w[q] = {T(1) / T(3),  T(1) / T(18), T(1) / T(18), T(1) / T(18),
                             T(1) / T(18), T(1) / T(18), T(1) / T(18), T(1) / T(36),
                             T(1) / T(36), T(1) / T(36), T(1) / T(36), T(1) / T(36),
                             T(1) / T(36), T(1) / T(36), T(1) / T(36), T(1) / T(36),
                             T(1) / T(36), T(1) / T(36), T(1) / T(36)};
  static constexpr int opp[q] = {0, 2,  1,  4,  3,  6,  5,  8,  7, 10,
                                 9, 12, 11, 14, 13, 16, 15, 18, 17};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
};

namespace latsetdata {
template <>
__constexpr__ Vector<int, 3> c<3, 19>[19] = {
  {0, 0, 0},  {1, 0, 0},   {-1, 0, 0}, {0, 1, 0},   {0, -1, 0}, {0, 0, 1},   {0, 0, -1},
  {1, 1, 0},  {-1, -1, 0}, {1, 0, 1},  {-1, 0, -1}, {0, 1, 1},  {0, -1, -1}, {1, -1, 0},
  {-1, 1, 0}, {1, 0, -1},  {-1, 0, 1}, {0, 1, -1},  {0, -1, 1}};

template <>
__constexpr__ Fraction<> w<3, 19>[19] = {{1, 3},  {1, 18}, {1, 18}, {1, 18}, {1, 18},
                                         {1, 18}, {1, 18}, {1, 36}, {1, 36}, {1, 36},
                                         {1, 36}, {1, 36}, {1, 36}, {1, 36}, {1, 36},
                                         {1, 36}, {1, 36}, {1, 36}, {1, 36}};

template <>
__constexpr__ int opp<3, 19>[19] = {0, 2,  1,  4,  3,  6,  5,  8,  7, 10,
                                    9, 12, 11, 14, 13, 16, 15, 18, 17};

}  // namespace latsetdata

// D3Q27
template <typename T>
struct D3Q27 : public Basic_Lattice_Set<3, 27> {
  using FloatType = T;
  static constexpr Vector<int, 3> c[q] = {
    {0, 0, 0},  {1, 0, 0},    {-1, 0, 0}, {0, 1, 0},   {0, -1, 0}, {0, 0, 1},
    {0, 0, -1},  // 0-6
    {1, 1, 0},  {-1, -1, 0},  {1, 0, 1},  {-1, 0, -1}, {0, 1, 1},  {0, -1, -1},
    {1, -1, 0}, {-1, 1, 0},   {1, 0, -1}, {-1, 0, 1},  {0, 1, -1}, {0, -1, 1},  // 7-18
    {1, 1, 1},  {-1, -1, -1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1},
    {-1, 1, 1}, {1, -1, -1}};  // 19-26
  static constexpr T w[q] = {T(8) / T(27),  T(2) / T(27),  T(2) / T(27),  T(2) / T(27),
                             T(2) / T(27),  T(2) / T(27),  T(2) / T(27),  T(1) / T(54),
                             T(1) / T(54),  T(1) / T(54),  T(1) / T(54),  T(1) / T(54),
                             T(1) / T(54),  T(1) / T(54),  T(1) / T(54),  T(1) / T(54),
                             T(1) / T(54),  T(1) / T(54),  T(1) / T(54),  T(1) / T(216),
                             T(1) / T(216), T(1) / T(216), T(1) / T(216), T(1) / T(216),
                             T(1) / T(216), T(1) / T(216), T(1) / T(216)};
  static constexpr int opp[q] = {0,  2,  1,  4,  3,  6,  5,  8,  7,  10, 9,  12, 11, 14,
                                 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
};

namespace latsetdata {
template <>
__constexpr__ Vector<int, 3> c<3, 27>[27] = {
  {0, 0, 0},  {1, 0, 0},   {-1, 0, 0}, {0, 1, 0},   {0, -1, 0}, {0, 0, 1},   {0, 0, -1},
  {1, 1, 0},  {-1, -1, 0}, {1, 0, 1},  {-1, 0, -1}, {0, 1, 1},  {0, -1, -1}, {1, -1, 0},
  {-1, 1, 0}, {1, 0, -1},  {-1, 0, 1}, {0, 1, -1},  {0, -1, 1}, {1, 1, 1},   {-1, -1, -1},
  {1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1}, {-1, 1, 1}, {1, -1, -1}};

template <>
__constexpr__ Fraction<> w<3, 27>[27] = {
  {8, 27},  {2, 27},  {2, 27},  {2, 27},  {2, 27},  {2, 27},  {2, 27},
  {1, 54},  {1, 54},  {1, 54},  {1, 54},  {1, 54},  {1, 54},  {1, 54},
  {1, 54},  {1, 54},  {1, 54},  {1, 54},  {1, 54},  {1, 216}, {1, 216},
  {1, 216}, {1, 216}, {1, 216}, {1, 216}, {1, 216}, {1, 216}};

template <>
__constexpr__ int opp<3, 27>[27] = {0,  2,  1,  4,  3,  6,  5,  8,  7,
                                    10, 9,  12, 11, 14, 13, 16, 15, 18,
                                    17, 20, 19, 22, 21, 24, 23, 26, 25};

}  // namespace latsetdata