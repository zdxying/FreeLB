/* This file is part of FreeLB, modified from openLB and FluidX3D with the following copyright notice:
 *
 * // start of the original OpenLB's copyright notice
 * 
 * This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 * 
 * // end of the original OpenLB's copyright notice
 * 
 * FluidX3D: https://github.com/ProjectPhysX/FluidX3Ds
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

#include "data_struct/block_lattice.h"
#include "lbm/freeSurfaceBase.h"


namespace olbfs {

template <typename CELL>
typename CELL::FloatType getClampedVOF(CELL& cell) {
  using T = typename CELL::FloatType;
  return std::clamp(cell.template get<VOLUMEFRAC<T>>(), T{}, T{1});
}

template <typename CELL>
static bool hasNeighborType(CELL& cell, FSType fstype) {
  using LatSet = typename CELL::LatticeSet;
  for (int i = 1; i < LatSet::q; ++i) {
    if (util::isFlag(cell.template getField<STATE>().get(cell.getNeighborId(i)), fstype))
      return true;
  }
  return false;
}

template <typename CELL>
static bool hasNeighborFlag(CELL& cell, FSFlag fsflag) {
  using LatSet = typename CELL::LatticeSet;
  for (int i = 1; i < LatSet::q; ++i) {
    if (util::isFlag(cell.template getField<FLAG>().get(cell.getNeighborId(i)), fsflag))
      return true;
  }
  return false;
}

// Parker-Youngs normal
// Parker-Youngs weights, a corrected version:
// |c|^2  1  2  3
// weight 4  2  1
// openlb 1  2  4
template <typename LatSet>
constexpr std::array<int, LatSet::q> Parker_YoungsWeights() {
  return make_Array<int, LatSet::q>([&](unsigned int i) {
    // int weight = LatSet::d == 2 ? 4 : 8;
    // int weight = 8;
    // if (latset::c<LatSet>(i)[0] != 0) weight /= 2;
    // if (latset::c<LatSet>(i)[1] != 0) weight /= 2;
    // if constexpr (LatSet::d == 3) {
    //   if (latset::c<LatSet>(i)[2] != 0) weight /= 2;
    // }
    // return weight;
    
    int weight = 1;
    if (latset::c<LatSet>(i)[0] != 0) weight *= 2;
    if (latset::c<LatSet>(i)[1] != 0) weight *= 2;
    if constexpr (LatSet::d == 3) {
      if (latset::c<LatSet>(i)[2] != 0) weight *= 2;
    }
    return weight/2;
  });
}

// for D2Q9 and D3Q27 only
template <typename CELL>
void computeParker_YoungsNormal(
  CELL& cell, Vector<typename CELL::FloatType, CELL::LatticeSet::d>& normal) {
  using T = typename CELL::FloatType;
  using LatSet = std::conditional_t<CELL::LatticeSet::d == 2, D2Q9<T>, D3Q27<T>>;
  for (int i = 1; i < LatSet::q; ++i) {
    CELL celln = cell.getNeighbor(latset::c<LatSet>(i));
    const T clampedvof = getClampedVOF(celln);
    normal -= Parker_YoungsWeights<LatSet>()[i] * latset::c<LatSet>(i) * clampedvof;
  }
}

// openlb offset calculation
template <typename T>
T offsetHelper2D(T volume, const std::vector<T>& sorted_normal) {
  T d2 = volume * sorted_normal[1] + 0.5 * sorted_normal[0];
  if (d2 >= sorted_normal[0]) {
    return d2;
  }
  T d1 = std::sqrt(2. * sorted_normal[0] * sorted_normal[1] * volume);
  return d1;
}

// openlb
// A lot of magic numbers are happening here. Optimized algorithm taken from Moritz
// Lehmann
template <typename T>
T offsetHelper3D(T volume, const std::vector<T>& sorted_normal) {
  T sn0_plus_sn1 = sorted_normal[0] + sorted_normal[1];
  T sn0_times_sn1 = sorted_normal[0] * sorted_normal[1];
  T sn2_volume = sorted_normal[2] * volume;

  T min_sn0_plus_sn1_and_sn2 = std::min(sn0_plus_sn1, sorted_normal[2]);

  T d5 = sn2_volume + 0.5 * sn0_plus_sn1;
  if (d5 > min_sn0_plus_sn1_and_sn2 && d5 <= sorted_normal[2]) {
    return d5;
  }

  T d2 = 0.5 * sorted_normal[0] +
         0.28867513 * std::sqrt(std::max(0., 24. * sorted_normal[1] * sn2_volume -
                                               sorted_normal[0] * sorted_normal[0]));

  if (d2 > sorted_normal[0] && d2 <= sorted_normal[1]) {
    return d2;
  }
  // cubic root
  T d1 = std::cbrt(6.0 * sn0_times_sn1 * sn2_volume);
  if (d1 <= sorted_normal[0]) {
    return d1;
  }

  T x3 = 81.0 * sn0_times_sn1 * (sn0_plus_sn1 - 2. * sn2_volume);
  T y3 = std::sqrt(
    std::max(0., 23328. * sn0_times_sn1 * sn0_times_sn1 * sn0_times_sn1 - x3 * x3));
  T u3 = std::cbrt(x3 * x3 + y3 * y3);
  T d3 = sn0_plus_sn1 - (7.5595264 * sn0_times_sn1 + 0.26456684 * u3) *
                          (1. / std::sqrt(u3)) *
                          std::sin(0.5235988 - 0.3333334 * std::atan(y3 / x3));
  if (d3 > sorted_normal[1] && d3 <= min_sn0_plus_sn1_and_sn2) {
    return d3;
  }

  T t4 = 9. * std::pow(sn0_plus_sn1 + sorted_normal[2], 2) - 18.;
  T x4 =
    std::max(sn0_times_sn1 * sorted_normal[2] * (324. - 648. * volume), 1.1754944e-38);
  T y4 = std::sqrt(std::max(4. * t4 * t4 * t4 - x4 * x4, 0.));
  T u4 = std::cbrt(x4 * x4 + y4 * y4);
  T d4 = 0.5 * (sn0_plus_sn1 + sorted_normal[2]) -
         (0.20998684 * t4 + 0.13228342 * u4) * (1. / std::sqrt(u4)) *
           std::sin(0.5235988 - 0.3333334 * std::atan(y4 / x4));

  return d4;
}

// openlb
// A lot of magic numbers are happening here. Optimized algorithm taken from Moritz
// Lehmann
template <typename T, typename LatSet>
T offsetHelperOpt(T vol, const Vector<T, LatSet::d>& sn) {
  const T sn0_p_sn1 = sn[0] + sn[1];
  const T sn2_t_V = sn[2] * vol;

  if (sn0_p_sn1 <= 2. * sn2_t_V) {
    return sn2_t_V + 0.5 * sn0_p_sn1;
  }

  const T sq_sn0 = std::pow(sn[0], 2), sn1_6 = 6. * sn[1], v1 = sq_sn0 / sn1_6;

  if (v1 <= sn2_t_V && sn2_t_V < v1 + 0.5 * (sn[1] - sn[0])) {
    return 0.5 * (sn[0] + std::sqrt(sq_sn0 + 8.0 * sn[1] * (sn2_t_V - v1)));
  }

  const T v6 = sn[0] * sn1_6 * sn2_t_V;
  if (sn2_t_V < v1) {
    return std::cbrt(v6);
  }

  const T v3 = sn[2] < sn0_p_sn1 ? (std::pow(sn[2], 2) * (3. * sn0_p_sn1 - sn[2]) +
                                    sq_sn0 * (sn[0] - 3.0 * sn[2]) +
                                    std::pow(sn[1], 2) * (sn[1] - 3.0 * sn[2])) /
                                     (sn[0] * sn1_6)
                                 : 0.5 * sn0_p_sn1;

  const T sq_sn0_sq_sn1 = sq_sn0 + std::pow(sn[1], 2),
          v6_cb_sn0_sn1 = v6 - std::pow(sn[0], 3) - std::pow(sn[1], 3);

  const bool case34 = sn2_t_V < v3;
  const T a = case34 ? v6_cb_sn0_sn1 : 0.5 * (v6_cb_sn0_sn1 - std::pow(sn[2], 3));
  const T b = case34 ? sq_sn0_sq_sn1 : 0.5 * (sq_sn0_sq_sn1 + std::pow(sn[2], 2));
  const T c = case34 ? sn0_p_sn1 : 0.5;
  const T t = std::sqrt(std::pow(c, 2) - b);
  return c - 2.0 * t *
               std::sin(0.33333334 * std::asin((std::pow(c, 3) - 0.5 * a - 1.5 * b * c) /
                                               std::pow(t, 3)));
}

// openlb
template <typename T, typename LatSet>
T calculateCubeOffset(T volume, const Vector<T, LatSet::d>& normal) {
  std::vector<T> abs_normal(LatSet::d, T{});
  for (unsigned int i = 0; i < LatSet::d; i++) {
    abs_normal[i] = std::abs(normal[i]);
  }

  T volume_symmetry = 0.5 - std::abs(volume - 0.5);

  std::sort(abs_normal.begin(), abs_normal.end());

  if constexpr (LatSet::d == 2) {
    abs_normal[0] = std::max(normal[0], 1e-5);
  } else if (LatSet::d == 3) {
    abs_normal[0] = std::max(normal[0], 1e-12);
    abs_normal[1] = std::max(normal[1], 1e-12);
  }

  T d{};
  if constexpr (LatSet::d == 2) {
    d = offsetHelper2D<T>(volume_symmetry, abs_normal);
  } else if (LatSet::d == 3) {
    d = offsetHelper3D<T>(volume_symmetry, abs_normal);
  }

  T sorted_normal_acc = 0;
  for (unsigned int i = 0; i < LatSet::d; i++) {
    sorted_normal_acc += abs_normal[i];
  }

  return std::copysign(d - 0.5 * sorted_normal_acc, volume - 0.5);
}

// openlb Optimized version of calculateCubeOffset
template <typename T, typename LatSet>
T calculateCubeOffsetOpt(T volume, const Vector<T, LatSet::d>& normal) {
  Vector<T, LatSet::d> abs_normal = normal.getabs();

  T a_l1 = abs_normal[0] + abs_normal[1] + abs_normal[2];

  T volume_symmetry = 0.5 - std::abs(volume - 0.5);

  Vector<T, LatSet::d> sorted_normal{};
  sorted_normal[0] =
    std::min(std::min(abs_normal[0], abs_normal[1]), abs_normal[2]) / a_l1;
  sorted_normal[1] = 0.;
  sorted_normal[2] =
    std::max(std::max(abs_normal[0], abs_normal[1]), abs_normal[2]) / a_l1;

  sorted_normal[1] = std::max(1. - sorted_normal[0] - sorted_normal[2], 0.);

  T d = offsetHelperOpt<T, LatSet>(volume_symmetry, sorted_normal);

  return a_l1 * std::copysign(0.5 - d, volume - 0.5);
}

// openlb pivoted LU solver
template <typename T, size_t S>
std::array<T, S> solvePivotedLU(std::array<std::array<T, S>, S>& matrix,
                                const std::array<T, S>& b, size_t N) {
  std::array<T, S> x;
  std::array<T, S> pivots;
  for (size_t i = 0; i < S; ++i) {
    pivots[i] = i;
    x[i] = 0.;
  }

  N = std::min(N, S);

  for (size_t i = 0; i < N; ++i) {
    T max = 0.;
    size_t max_index = i;

    for (size_t j = i; j < N; ++j) {
      T abs = std::abs(matrix[pivots[j]][i]);
      if (abs > max) {
        max_index = j;
        max = abs;
      }
    }

    if (max_index != i) {
      size_t tmp_index = pivots[i];
      pivots[i] = pivots[max_index];
      pivots[max_index] = tmp_index;
    }

    for (size_t j = i + 1; j < N; ++j) {
      matrix[pivots[j]][i] /= matrix[pivots[i]][i];

      for (size_t k = i + 1; k < N; ++k) {
        matrix[pivots[j]][k] -= matrix[pivots[j]][i] * matrix[pivots[i]][k];
      }
    }
  }

  for (size_t i = 0; i < N; ++i) {
    x[i] = b[pivots[i]];

    for (size_t j = 0; j < i; ++j) {
      x[i] -= matrix[pivots[i]][j] * x[j];
    }
  }

  for (size_t i = N; i > 0; --i) {
    for (size_t j = i; j < N; ++j) {
      x[i - 1] -= matrix[pivots[i - 1]][j] * x[j];
    }

    x[i - 1] /= matrix[pivots[i - 1]][i - 1];
  }

  return x;
}


// openlb
template <typename CELL>
typename CELL::FloatType ComputeCurvature2D(CELL& cell) {
  using T = typename CELL::FloatType;
  // using LatSet = typename CELL::LatticeSet;
  using LatSet = D2Q9<T>;

  // Parker-Youngs Normal
  Vector<T, 2> normal{};
  computeParker_YoungsNormal(cell, normal);
  T norm = normal.getnorm();
  if (norm < 1e-6) return T{};
  normal /= norm;

  // Rotation matrix is
  // ( n1 | -n0 )
  // ( n0 |  n1 )

  // It is 2 because of the amount of fitting parameters. Not because of the dimension
  constexpr size_t S = 2;
  std::array<std::array<T, S>, S> lq_matrix;
  std::array<T, S> b_rhs;
  for (size_t i = 0; i < S; ++i) {
    for (size_t j = 0; j < S; ++j) {
      lq_matrix[i][j] = 0.;
    }
    b_rhs[i] = 0.;
  }

  // Offset for the plic correction
  T origin_offset{};
  T VOF = getClampedVOF(cell);
  origin_offset = calculateCubeOffset<T, LatSet>(VOF, normal);

  std::size_t healthy_interfaces = 0;

  for (int iPop = 1; iPop < LatSet::q; ++iPop) {
    CELL cellC = cell.getNeighbor(latset::c<LatSet>(iPop));
    if (!util::isFlag(cellC.template get<STATE>(), FSType::Interface) ||
        !hasNeighborType(cellC, FSType::Gas)) {
      continue;
    }

    ++healthy_interfaces;

    T fill_level = getClampedVOF(cellC);

    T cube_offset = calculateCubeOffset<T, LatSet>(fill_level, normal);

    T x_pos = latset::c<LatSet>(iPop)[0];
    T y_pos = latset::c<LatSet>(iPop)[1];

    // Rotation
    T rot_x_pos = x_pos * normal[1] - y_pos * normal[0];
    T rot_y_pos = x_pos * normal[0] + y_pos * normal[1] + (cube_offset - origin_offset);

    T rot_x_pos_2 = rot_x_pos * rot_x_pos;
    T rot_x_pos_3 = rot_x_pos_2 * rot_x_pos;
    T rot_x_pos_4 = rot_x_pos_3 * rot_x_pos;

    lq_matrix[1][1] += rot_x_pos_2;
    lq_matrix[1][0] += rot_x_pos_3;
    lq_matrix[0][0] += rot_x_pos_4;

    b_rhs[0] += rot_x_pos_2 * (rot_y_pos);
    b_rhs[1] += rot_x_pos * (rot_y_pos);
  }

  lq_matrix[0][1] = lq_matrix[1][0];

  // Thikonov regularization parameter
  T alpha{};
  for (size_t i = 0; i < LatSet::d; ++i) {
    lq_matrix[i][i] += alpha;
  }

  // It is 2 because of the fitting parameters. Not dependent on the dimension
  std::array<T, S> solved_fit =
    solvePivotedLU<T, S>(lq_matrix, b_rhs, healthy_interfaces);

  // signed curvature -> kappa = y'' / ( (1 + y'²)^(3/2) )
  T denom = std::sqrt(1. + solved_fit[1] * solved_fit[1]);
  denom = denom * denom * denom;
  T curvature = 2. * solved_fit[0] / denom;
  return std::max(-1., std::min(1., curvature));
}

// openlb
template <typename CELL>
typename CELL::FloatType ComputeCurvature3D(CELL& cell) {
  using T = typename CELL::FloatType;
  // using LatSet = typename CELL::LatticeSet;
  using LatSet = D3Q27<T>;
  // This is b_z
  // Parker-Youngs Normal
  Vector<T, 3> normal{};
  computeParker_YoungsNormal(cell, normal);
  T norm = normal.getnorm();
  if (norm < 1e-6) return T{};
  normal /= norm;

  std::array<T, 3> r_vec{0.56270900, 0.32704452, 0.75921047};
  /*
  std::array<T,DESCRIPTOR::d> r_vec{
    0.,0.,1.
  };
  */
  std::array<std::array<T, 3>, 3> rotation{
    {{{0., 0., 0.}},
     //{{normal[1], -normal[0], 0.}},
     {{normal[1] * r_vec[2] - normal[2] * r_vec[1],
       normal[2] * r_vec[0] - normal[0] * r_vec[2],
       normal[0] * r_vec[1] - normal[1] * r_vec[0]}},
     {{normal[0], normal[1], normal[2]}}}};

  // Cross product with (0,0,1) x normal
  // This is b_y

  // (normal[0], normal[1], normal[2])

  T cross_norm = 0.;
  for (size_t i = 0; i < LatSet::d; ++i) {
    cross_norm += rotation[1][i] * rotation[1][i];
  }

  // If too close too each other use the identity matrix
  if (cross_norm > 1e-6) {
    cross_norm = std::sqrt(cross_norm);
    for (size_t i = 0; i < LatSet::d; ++i) {
      rotation[1][i] /= cross_norm;
    }
  } else {
    rotation[1] = {{-normal[2], 0., normal[0]}};

    cross_norm = 0.;
    for (size_t i = 0; i < LatSet::d; ++i) {
      cross_norm += rotation[1][i] * rotation[1][i];
    }

    cross_norm = std::sqrt(cross_norm);

    for (size_t i = 0; i < LatSet::d; ++i) {
      rotation[1][i] /= cross_norm;
    }
  }

  // Cross product of ((0,0,1) x normal / | (0,0,1) x normal |) x normal
  // This is b_x
  rotation[0] = {{rotation[1][1] * normal[2] - rotation[1][2] * normal[1],
                  rotation[1][2] * normal[0] - rotation[1][0] * normal[2],
                  rotation[1][0] * normal[1] - rotation[1][1] * normal[0]}};

  // These three form a matrix and are entered into each row
  // ( b_x )
  // ( b_y )
  // ( b_z )

  constexpr size_t S = 5;
  std::array<std::array<T, S>, S> lq_matrix;
  std::array<T, S> b_rhs;
  for (size_t i = 0; i < S; ++i) {
    for (size_t j = 0; j < S; ++j) {
      lq_matrix[i][j] = 0.;
    }
    b_rhs[i] = 0.;
  }
  T origin_offset{};
  {
    T fill_level = getClampedVOF(cell);
    origin_offset = calculateCubeOffsetOpt<T, LatSet>(fill_level, normal);
  }

  size_t healthy_interfaces = 0;
  for (int iPop = 1; iPop < LatSet::q; iPop++) {
    auto cellC = cell.getNeighbor(latset::c<LatSet>(iPop));

    if (!util::isFlag(cellC.template get<STATE>(), FSType::Interface) ||
        !hasNeighborType(cellC, FSType::Gas)) {
      continue;
    }

    ++healthy_interfaces;

    T fill_level = getClampedVOF(cellC);

    T cube_offset = calculateCubeOffsetOpt<T, LatSet>(fill_level, normal);

    int i = latset::c<LatSet>(iPop)[0];
    int j = latset::c<LatSet>(iPop)[1];
    int k = latset::c<LatSet>(iPop)[2];

    std::array<T, 3> pos{static_cast<T>(i), static_cast<T>(j), static_cast<T>(k)};
    std::array<T, 3> r_pos{0., 0., cube_offset - origin_offset};

    for (size_t a = 0; a < LatSet::d; ++a) {
      for (size_t b = 0; b < LatSet::d; ++b) {
        r_pos[a] += rotation[a][b] * pos[b];
      }
    }

    T r_x_2 = r_pos[0] * r_pos[0];
    T r_x_3 = r_x_2 * r_pos[0];
    T r_x_4 = r_x_3 * r_pos[0];

    T r_y_2 = r_pos[1] * r_pos[1];
    T r_y_3 = r_y_2 * r_pos[1];
    T r_y_4 = r_y_3 * r_pos[1];

    T r_x_2_y_2 = r_x_2 * r_y_2;
    T r_x_3_y = r_x_3 * r_pos[1];
    T r_x_2_y = r_x_2 * r_pos[1];

    T r_x_y_3 = r_pos[0] * r_y_3;
    T r_x_y_2 = r_pos[0] * r_y_2;

    T r_x_y = r_pos[0] * r_pos[1];

    lq_matrix[0][0] += r_x_4;
    lq_matrix[1][1] += r_y_4;
    lq_matrix[2][2] += r_x_2_y_2;
    lq_matrix[3][3] += r_x_2;
    lq_matrix[4][4] += r_y_2;

    // skip [1][0] copy later from [2][2]
    lq_matrix[2][0] += r_x_3_y;
    lq_matrix[3][0] += r_x_3;
    lq_matrix[4][0] += r_x_2_y;

    lq_matrix[2][1] += r_x_y_3;
    lq_matrix[3][1] += r_x_y_2;
    lq_matrix[4][1] += r_y_3;

    // skip [3][2] copy from [4][0]
    // skip [4][2] copy from [3][1]

    lq_matrix[4][3] += r_x_y;

    b_rhs[0] += r_x_2 * r_pos[2];
    b_rhs[1] += r_y_2 * r_pos[2];
    b_rhs[2] += r_x_y * r_pos[2];
    b_rhs[3] += r_pos[0] * r_pos[2];
    b_rhs[4] += r_pos[1] * r_pos[2];
  }

  lq_matrix[1][0] = lq_matrix[2][2];
  lq_matrix[3][2] = lq_matrix[4][0];
  lq_matrix[4][2] = lq_matrix[3][1];

  for (size_t i = 0; i < S; ++i) {
    for (size_t j = i + 1; j < S; ++j) {
      lq_matrix[i][j] = lq_matrix[j][i];
    }
  }

  // Consider using Thikonov regularization?
  // T alpha = 1e-8;
  T alpha{};
  for (size_t i = 0; i < S; ++i) {
    lq_matrix[i][i] += alpha;
  }

  std::array<T, S> solved_fit =
    solvePivotedLU<T, S>(lq_matrix, b_rhs, healthy_interfaces);

  T denom = std::sqrt(1. + solved_fit[3] * solved_fit[3] + solved_fit[4] * solved_fit[4]);
  denom = denom * denom * denom;
  T curvature = ((1. + solved_fit[4] * solved_fit[4]) * solved_fit[0] +
                 (1. + solved_fit[3] * solved_fit[3]) * solved_fit[1] -
                 solved_fit[3] * solved_fit[4] * solved_fit[2]) /
                denom;

  return std::max(-1., std::min(1., curvature));
}


template <typename CELL>
typename CELL::FloatType ComputeCurvature(CELL& cell) {
  using LatSet = typename CELL::LatticeSet;
  if constexpr (LatSet::d == 2) {
    return ComputeCurvature2D(cell);
  } else if constexpr (LatSet::d == 3) {
    return ComputeCurvature3D(cell);
  }
}

}  // namespace olbfs


// FluidX3D's implementation of free surface model

namespace fx3dfs {

template <typename CELL>
typename CELL::FloatType getClampedVOF(CELL& cell) {
  using T = typename CELL::FloatType;
  return std::clamp(cell.template get<VOLUMEFRAC<T>>(), T{}, T{1});
}

template <typename CELL>
typename CELL::FloatType computeVOF(CELL& cell) {
  using T = typename CELL::FloatType;
  const T mass = cell.template get<MASS<T>>();
  const T rho = cell.template get<RHO<T>>();
  if (util::isFlag(cell.template get<STATE>(), FSType::Fluid)) {
    return T{1};
  } else if (util::isFlag(cell.template get<STATE>(), FSType::Interface)) {
    return rho > T{} ? std::clamp(mass / rho, T{}, T{1}) : T{0.5};
  } else {
    return T{};
  }
}
template <typename T>
T computeVOF(T mass, T rho) {
  return rho > T{} ? std::clamp(mass / rho, T{}, T{1}) : T{0.5};
}

template <typename CELL>
static bool hasNeighborType(CELL& cell, FSType fstype) {
  using LatSet = typename CELL::LatticeSet;
  for (int i = 1; i < LatSet::q; ++i) {
    if (util::isFlag(cell.template getField<STATE>().get(cell.getNeighborId(i)), fstype))
      return true;
  }
  return false;
}

// Parker-Youngs normal
// Parker-Youngs weights, a corrected version:
// D2Q9
// |c|^2  1  2
// e.g. (1,0) (1,1)
// weight 2  1
// D3Q27
// |c|^2  1  2  3
// e.g. (1,0,0) (1,1,0) (1,1,1)
// weight 4  2  1
// for D2Q9 and D3Q27 only
template <typename LatSet>
constexpr auto ParkerYoungsWeights() {
  if constexpr (LatSet::d == 2) {
    return make_Array<Vector<int, 2>, 9>([&](unsigned int i) {
      int weight = 4;
      if (latset::c<LatSet>(i)[0] != 0) weight /= 2;
      if (latset::c<LatSet>(i)[1] != 0) weight /= 2;
      return weight * latset::c<LatSet>(i);
    });
  } else if constexpr (LatSet::d == 3) {
    return make_Array<Vector<int, 3>, 27>([&](unsigned int i) {
      int weight = 8;
      if (latset::c<LatSet>(i)[0] != 0) weight /= 2;
      if (latset::c<LatSet>(i)[1] != 0) weight /= 2;
      if (latset::c<LatSet>(i)[2] != 0) weight /= 2;
      return weight * latset::c<LatSet>(i);
    });
  }
}

// for D2Q9 and D3Q27 only
template <typename CELL>
void computeParkerYoungsNormal(
  CELL& cell, Vector<typename CELL::FloatType, CELL::LatticeSet::d>& normal) {
  using T = typename CELL::FloatType;
  using LatSet = std::conditional_t<CELL::LatticeSet::d == 2, D2Q9<T>, D3Q27<T>>;
  for (int i = 1; i < LatSet::q; ++i) {
    CELL celln = cell.getNeighbor(latset::c<LatSet>(i));
    normal -= ParkerYoungsWeights<LatSet>()[i] * computeVOF(celln);
  }
}

template <typename T>
T PlicCube_Reduced(const T V, const T n1, const T n2, const T n3) {
  // optimized solution from SZ and Kawano, source:
  // https://doi.org/10.3390/computation10020021

  const T n12 = n1 + n2, n3V = n3 * V;
  if (n12 <= T{2} * n3V) return n3V + T{0.5} * n12;  // case (5)
  const T sqn1 = n1 * n1, n26 = T{6} * n2,
          v1 = sqn1 / n26;  // after case (5) check n2>0 is true
  if (v1 <= n3V && n3V < v1 + T{0.5} * (n2 - n1))
    return T{0.5} * (n1 + std::sqrt(sqn1 + T{8} * n2 * (n3V - v1)));  // case (2)
  const T V6 = n1 * n26 * n3V;
  if (n3V < v1) return std::cbrt(V6);  // case (1)
  const T v3 = n3 < n12 ? (n3 * n3 * (T{3} * n12 - n3) + sqn1 * (n1 - T{3} * n3) +
                           n2 * n2 * (n2 - T{3} * n3)) /
                            (n1 * n26)
                        : T{0.5} * n12;  // after case (2) check n1>0 is true
  const T sqn12 = sqn1 + n2 * n2, V6cbn12 = V6 - n1 * n1 * n1 - n2 * n2 * n2;
  const bool case34 = n3V < v3;  // true: case (3), false: case (4)
  const T a = case34 ? V6cbn12 : T{0.5} * (V6cbn12 - n3 * n3 * n3);
  const T b = case34 ? sqn12 : T{0.5} * (sqn12 + n3 * n3);
  const T c = case34 ? n12 : T{0.5};
  const T t = std::sqrt(c * c - b);
  return c -
         T{2} * t *
           std::sin(T{0.33333334} *
                    std::asin(((c * c * c) - T{0.5} * a - T{1.5} * b * c) / (t * t * t)));
}

template <typename T>
T PlicCube(const T ClampedVOF, const Vector<T, 3>& normal) {
  // unit cube - plane intersection:
  // ClampedVOF in [0,1], normal vector n -> plane offset d0

  // eliminate symmetry cases, normalize n using L1 norm
  const T ax = std::fabs(normal[0]);
  const T ay = std::fabs(normal[1]);
  const T az = std::fabs(normal[2]);
  const T sum = ax + ay + az;
  const T V = T{0.5} - std::fabs(ClampedVOF - T{0.5});
  const T n1 = std::min(std::min(ax, ay), az) / sum;
  const T n3 = std::max(std::max(ax, ay), az) / sum;
  const T n2 = std::fdim(T{1}, n1 + n3);        // ensure n2>=0 //max(x-y,0)
  const T d = PlicCube_Reduced(V, n1, n2, n3);  // calculate PLIC with reduced symmetry
  return sum * std::copysign(
                 T{0.5} - d,
                 ClampedVOF - T{0.5});  // rescale result and apply symmetry for V0>0.5
}

template <typename T>
void LU_Solve(T* M, T* x, T* b, const int N, const int Nsol) {
  // solves system of N linear equations M*x=b within dimensionality Nsol<=N
  for (int i = 0; i < Nsol; i++) {  // decompose M in M=L*U
    for (int j = i + 1; j < Nsol; j++) {
      M[N * j + i] /= M[N * i + i];
      for (int k = i + 1; k < Nsol; k++) M[N * j + k] -= M[N * j + i] * M[N * i + k];
    }
  }
  for (int i = 0; i < Nsol; i++) {  // find solution of L*y=b
    x[i] = b[i];
    for (int k = 0; k < i; k++) x[i] -= M[N * i + k] * x[k];
  }
  for (int i = Nsol - 1; i >= 0; i--) {  // find solution of U*x=y
    for (int k = i + 1; k < Nsol; k++) x[i] -= M[N * i + k] * x[k];
    x[i] /= M[N * i + i];
  }
}

// for D2Q9
template <typename CELL>
typename CELL::FloatType ComputeCurvature2D(CELL& cell) {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  const T VOF = computeVOF(cell);
  // Parker-Youngs Normal
  Vector<T, 2> normal{};
  computeParkerYoungsNormal(cell, normal);
  const T norm = normal.getnorm();
  if (norm < 1e-6) return T{};
  normal /= norm;

  // new coordinate system: bz is normal to surface, bx and by are tangent to surface
  const Vector<T, 3> by{normal[0], normal[1], T{}};
  const Vector<T, 3> bx = CrossProduct3(by, Vector<T, 3>{T{}, T{}, T{1}});
  // const Vector<T, 3> bx{normal[1], -normal[0], T{}};

  // number of neighboring interface points
  unsigned int number{};
  Vector<T, 2> p[6];  // number of neighboring interface points is less or equal than than
                      // 8 minus 1 gas and minus 1 fluid point = 6
  const T center_offset =
    PlicCube(VOF, by);  // calculate z-offset PLIC of center point only once
  for (unsigned int i = 1; i < 9; ++i) {  // iterate over neighbors
    CELL celln = cell.getNeighbor(latset::c<D2Q9<T>>(i));
    const T nVOF = computeVOF(celln);
    if (nVOF > T{} && nVOF < T{1}) {  // limit neighbors to interface nodes
      const Vector<int, 3> ei{
        latset::c<D2Q9<T>>(i)[0], latset::c<D2Q9<T>>(i)[1],
        0};  // assume neighbor normal vector is the same as center normal vector
      const T offset = PlicCube(nVOF, by) - center_offset;
      // do coordinate system transformation into (x, f(x)) and apply PLIC pffsets
      p[number++] = Vector<T, 2>{ei * bx, ei * by + offset};  
    }
  }

  T M[4] = {T{}, T{}, T{}, T{}}, x[2] = {T{}, T{}}, b[2] = {T{}, T{}};
  for (unsigned int i = 0; i < number;
       ++i) {  // f(x,y)=A*x2+H*x, x=(A,H), Q=(x2,x), M*x=b, M=Q*Q^T, b=Q*z
    const T x = p[i][0], y = p[i][1], x2 = x * x, x3 = x2 * x;
    /**/ M[0] += x2 * x2;
    M[1] += x3;
    b[0] += x2 * y;
    /*M[2]+=x3   ;*/ M[3] += x2;
    b[1] += x * y;
  }
  M[2] = M[1];  // use symmetry of matrix to save arithmetic operations
  if (number >= (unsigned int){2})
    LU_Solve(M, x, b, 2, 2);
  else
    LU_Solve(M, x, b, 2, std::min((unsigned int){2}, number));
  const T A = x[0], H = x[1];
  const T temp = T{1} / std::sqrt(H * H + T{1});
  const T K =
    T{2} * A *
    (temp * temp * temp);  // mean curvature of Monge patch (x, f(x)), note that curvature
                           // definition in 2D is different than 3D (additional factor 2)

  if (std::isnan(K)) {
    if (std::signbit(K)) {
      return T{-1};
    } else {
      return T{1};
    }
  }
  // return std::max(T{-1}, std::min(T{1}, K));
  return (K < T{-1}) ? T{-1} : (T{1} < K) ? T{1} : K;
}

template <typename CELL>
typename CELL::FloatType ComputeCurvature3D(CELL& cell) {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  const T VOF = computeVOF(cell);
  // Parker-Youngs Normal
  // new coordinate system: bz is normal to surface, bx and by are tangent to surface
  Vector<T, 3> bz{};
  computeParkerYoungsNormal(cell, bz);
  const T norm = bz.getnorm();
  if (norm < 1e-6) return T{};
  bz /= norm;
  const Vector<T, 3> rn{T{0.56270900}, T{0.32704452},
                        T{0.75921047}};  // random normalized vector that is just by
                                         // random chance not collinear with bz
  const Vector<T, 3> by =
    CrossProduct3(bz, rn).getnormalize();  // normalize() is necessary here because bz and rn
                                        // are not perpendicular
  const Vector<T, 3> bx = CrossProduct3(by, bz);

  // number of neighboring interface points
  unsigned int number{};
  Vector<T, 3> p[24];  // number of neighboring interface points is less or equal than
                       // than 26 minus 1 gas and minus 1 fluid point = 24
  const T center_offset =
    PlicCube(VOF, bz);  // calculate z-offset PLIC of center point only once
  for (unsigned int i = 1; i < 27; ++i) {  // iterate over neighbors
    CELL celln = cell.getNeighbor(latset::c<D3Q27<T>>(i));
    const T nVOF = computeVOF(celln);
    if (nVOF > T{} && nVOF < T{1}) {  // limit neighbors to interface nodes
      // assume neighbor normal vector is the same as center normal vector
      const T offset = PlicCube(nVOF, bz) - center_offset;
      p[number++] = Vector<T, 3>{latset::c<D3Q27<T>>(i) * bx, latset::c<D3Q27<T>>(i) * by,
                                 latset::c<D3Q27<T>>(i) * bz +
                                   offset};  // do coordinate system transformation into
                                             // (x, y, f(x,y)) and apply PLIC pffsets
    }
  }
  T M[25], x[5] = {T{}, T{}, T{}, T{}, T{}}, b[5] = {T{}, T{}, T{}, T{}, T{}};
  for (unsigned int i = 0; i < 25; ++i) M[i] = T{};
  for (unsigned int i = 0; i < number;
       ++i) {  // f(x,y)=A*x2+B*y2+C*x*y+H*x+I*y, x=(A,B,C,H,I), Q=(x2,y2,x*y,x,y), M*x=b,
               // M=Q*Q^T, b=Q*z
    const T x = p[i][0], y = p[i][1], z = p[i][2], x2 = x * x, y2 = y * y, x3 = x2 * x,
            y3 = y2 * y;
    /**/ M[0] += x2 * x2;
    M[1] += x2 * y2;
    M[2] += x3 * y;
    M[3] += x3;
    M[4] += x2 * y;
    b[0] += x2 * z;
    /*M[ 5]+=x2*y2;*/ M[6] += y2 * y2;
    M[7] += x * y3;
    M[8] += x * y2;
    M[9] += y3;
    b[1] += y2 * z;
    /*M[10]+=x3*y ; M[11]+=x *y3;*/ M[12] += x2 * y2;
    M[13] += x2 * y;
    M[14] += x * y2;
    b[2] += x * y * z;
    /*M[15]+=x3   ; M[16]+=x *y2; M[17]+=x2*y ;*/ M[18] += x2;
    M[19] += x * y;
    b[3] += x * z;
    /*M[20]+=x2*y ; M[21]+=   y3; M[22]+=x *y2; M[23]+=x *y ;*/ M[24] += y2;
    b[4] += y * z;
  }
  for (unsigned int i = 1; i < 5;
       ++i) {  // use symmetry of matrix to save arithmetic operations
    for (unsigned int j = 0; j < i; ++j) M[i * 5 + j] = M[j * 5 + i];
  }
  if (number >= (unsigned int){5})
    LU_Solve(M, x, b, 5, 5);
  else
    LU_Solve(M, x, b, 5, std::min((unsigned int){5}, number));
  const T A = x[0], B = x[1], C = x[2], H = x[3], I = x[4];
  const T temp = T{1} / std::sqrt(H * H + I * I + T{1});
  const T K = (A * (I * I + 1.0f) + B * (H * H + 1.0f) - C * H * I) *
              (temp * temp * temp);  // mean curvature of Monge patch (x, y, f(x, y))
  return std::max(-1., std::min(1., K));
}

template <typename CELL>
typename CELL::FloatType ComputeCurvature(CELL& cell) {
  using LatSet = typename CELL::LatticeSet;
  if constexpr (LatSet::d == 2) {
    return ComputeCurvature2D(cell);
  } else if constexpr (LatSet::d == 3) {
    return ComputeCurvature3D(cell);
  }
}

}  // namespace fx3dfs