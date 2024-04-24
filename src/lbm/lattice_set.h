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
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 * 
 */

// lattice set for lbm and ca
// lattice_set.h

#pragma once

#include "data_struct/Vector.h"
#include "utils/util.h"

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
  static constexpr int c[q][d] = {{0}, {1}, {-1}};
  static constexpr T w[q] = {T(1) / T(3), T(1) / T(6), T(1) / T(6)};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
  static constexpr int opp[q] = {0, 2, 1};
};

// D2Q4
template <typename T>
struct D2Q4 : public Basic_Lattice_Set<2, 4> {
  static constexpr Vector<int, 2> c[q] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
  static constexpr T w[q] = {T(1) / T(4), T(1) / T(4), T(1) / T(4),
                             T(1) / T(4)};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
  static constexpr int opp[q] = {2, 3, 0, 1};
};

// D2Q5
template <typename T>
struct D2Q5 : public Basic_Lattice_Set<2, 5> {
  static constexpr Vector<int, 2> c[q] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}};
  static constexpr T w[q] = {T(1) / T(3), T(1) / T(6), T(1) / T(6), T(1) / T(6),
                             T(1) / T(6)};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
  static constexpr int opp[q] = {0, 3, 4, 1, 2};
};

// D2Q9
template <typename T>
struct D2Q9 : public Basic_Lattice_Set<2, 9> {
  static constexpr Vector<int, 2> c[q] = {{0, 0}, {1, 0},  {0, 1},   {-1, 0}, {0, -1},
                                  {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  static constexpr T w[q] = {T(4) / T(9),  T(1) / T(9),  T(1) / T(9),
                             T(1) / T(9),  T(1) / T(9),  T(1) / T(36),
                             T(1) / T(36), T(1) / T(36), T(1) / T(36)};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
  static constexpr int opp[q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
};

// D3Q7
template <typename T>
struct D3Q7 : public Basic_Lattice_Set<3, 7> {
  static constexpr Vector<int, 3> c[q] = {{0, 0, 0}, {1, 0, 0},  {-1, 0, 0},
                                          {0, 1, 0}, {0, -1, 0}, {0, 0, 1},
                                          {0, 0, -1}};
  static constexpr T w[q] = {T(1) / T(4), T(1) / T(8), T(1) / T(8), T(1) / T(8),
                             T(1) / T(8), T(1) / T(8), T(1) / T(8)};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
  static constexpr int opp[q] = {0, 2, 1, 4, 3, 6, 5};
};

// D3Q15
template <typename T>
struct D3Q15 : public Basic_Lattice_Set<3, 15> {
  static constexpr Vector<int, 3> c[q] = {
      {0, 0, 0},   {1, 0, 0},  {-1, 0, 0},  {0, 1, 0},    {0, -1, 0},
      {0, 0, 1},   {0, 0, -1}, {1, 1, 1},   {-1, -1, -1}, {1, 1, -1},
      {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1}, {-1, 1, 1},   {1, -1, -1}};
  static constexpr T w[q] = {
      T(2) / T(9),  T(1) / T(9),  T(1) / T(9),  T(1) / T(9),  T(1) / T(9),
      T(1) / T(9),  T(1) / T(9),  T(1) / T(72), T(1) / T(72), T(1) / T(72),
      T(1) / T(72), T(1) / T(72), T(1) / T(72), T(1) / T(72), T(1) / T(72)};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
  static constexpr int opp[q] = {0, 2,  1, 4,  3,  6,  5, 8,
                                 7, 10, 9, 12, 11, 14, 13};
};

// D3Q19
template <typename T>
struct D3Q19 : public Basic_Lattice_Set<3, 19> {
  static constexpr Vector<int, 3> c[q] = {
      {0, 0, 0},   {1, 0, 0},  {-1, 0, 0},  {0, 1, 0},   {0, -1, 0},
      {0, 0, 1},   {0, 0, -1}, {1, 1, 0},   {-1, -1, 0}, {1, 0, 1},
      {-1, 0, -1}, {0, 1, 1},  {0, -1, -1}, {1, -1, 0},  {-1, 1, 0},
      {1, 0, -1},  {-1, 0, 1}, {0, 1, -1},  {0, -1, 1}};
  static constexpr T w[q] = {
      T(1) / T(3),  T(1) / T(18), T(1) / T(18), T(1) / T(18), T(1) / T(18),
      T(1) / T(18), T(1) / T(18), T(1) / T(36), T(1) / T(36), T(1) / T(36),
      T(1) / T(36), T(1) / T(36), T(1) / T(36), T(1) / T(36), T(1) / T(36),
      T(1) / T(36), T(1) / T(36), T(1) / T(36), T(1) / T(36)};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
  static constexpr int opp[q] = {0, 2,  1,  4,  3,  6,  5,  8,  7, 10,
                                 9, 12, 11, 14, 13, 16, 15, 18, 17};
};

// D3Q27
template <typename T>
struct D3Q27 : public Basic_Lattice_Set<3, 27> {
  static constexpr Vector<int, 3> c[q] = {
      {0, 0, 0},  {1, 0, 0},    {-1, 0, 0}, {0, 1, 0},
      {0, -1, 0}, {0, 0, 1},    {0, 0, -1},  // 0-6
      {1, 1, 0},  {-1, -1, 0},  {1, 0, 1},  {-1, 0, -1},
      {0, 1, 1},  {0, -1, -1},  {1, -1, 0}, {-1, 1, 0},
      {1, 0, -1}, {-1, 0, 1},   {0, 1, -1}, {0, -1, 1},  // 7-18
      {1, 1, 1},  {-1, -1, -1}, {1, 1, -1}, {-1, -1, 1},
      {1, -1, 1}, {-1, 1, -1},  {-1, 1, 1}, {1, -1, -1}};  // 19-26
  static constexpr T w[q] = {
      T(8) / T(27),  T(2) / T(27),  T(2) / T(27),  T(2) / T(27),  T(2) / T(27),
      T(2) / T(27),  T(2) / T(27),  T(1) / T(54),  T(1) / T(54),  T(1) / T(54),
      T(1) / T(54),  T(1) / T(54),  T(1) / T(54),  T(1) / T(54),  T(1) / T(54),
      T(1) / T(54),  T(1) / T(54),  T(1) / T(54),  T(1) / T(54),  T(1) / T(216),
      T(1) / T(216), T(1) / T(216), T(1) / T(216), T(1) / T(216), T(1) / T(216),
      T(1) / T(216), T(1) / T(216)};
  static constexpr T cs2 = T(1) / T(3);
  static constexpr T InvCs2 = T(3);
  static constexpr T InvCs4 = T(9);
  static constexpr int opp[q] = {0,  2,  1,  4,  3,  6,  5,  8,  7,
                                 10, 9,  12, 11, 14, 13, 16, 15, 18,
                                 17, 20, 19, 22, 21, 24, 23, 26, 25};
};

