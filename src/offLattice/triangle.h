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

#pragma once

#include "data_struct/Vector.h"
#include "utils/util.h"

namespace offlat {

// only for 3D case

template <typename T>
using Triangle = Vector<Vector<T, 3>, 3>;

template <typename T>
struct LinearInterp {
  // the first of 2 indices of the edge where one vertex is on
  std::size_t id1;
  // the second of 2 indices of the edge where one vertex is on
  std::size_t id2;
  // T mu = (isoValue - val1) / (val2 - val1);
  // intp = p1 + mu * (p2 - p1)
  T mu;

  LinearInterp() = default;
  LinearInterp(std::size_t id1, std::size_t id2, T mu) : id1(id1), id2(id2), mu(mu) {}
  
  inline T getIntp(const T val1, const T val2) const { return val1 + mu * (val2 - val1); }
  template <unsigned int D>
  inline Vector<T, D> getIntp(const Vector<T, D>& val1, const Vector<T, D>& val2) const {
    return val1 + mu * (val2 - val1);
  }
};

// TriangleIdx contains 6 points' indices of the 3 edges where the triangle's vertices are on
// and 3 values for linear interpolation
template <typename T>
using TriangleIdx = Vector<LinearInterp<T>, 3>;

template <typename T>
void TriangleNormal(Vector<T, 3>& normal, const Triangle<T>& triangle) {
  Vector<T, 3> e0 = triangle[1] - triangle[0];
  Vector<T, 3> e1 = triangle[2] - triangle[0];

  normal = CrossProduct3(e0, e1);
  const T norm = normal.getnorm();
  if (util::nearZero(norm)) {
    normal = Vector<T, 3>{};
  } else {
    normal /= norm;
  }
}
template <typename T>
Vector<T, 3> getTriangleNormal(const Triangle<T>& triangle) {
  Vector<T, 3> normal;
  TriangleNormal(normal, triangle);
  return normal;
}
template <typename T>
void TriangleNormal(Vector<T, 3>& normal, const Vector<T, 3>& v0, const Vector<T, 3>& v1,
                    const Vector<T, 3>& v2) {
  Vector<T, 3> e0 = v1 - v0;
  Vector<T, 3> e1 = v2 - v0;

  normal = CrossProduct3(e0, e1);
  const T norm = normal.getnorm();
  if (util::nearZero(norm)) {
    normal = Vector<T, 3>{};
  } else {
    normal /= norm;
  }
}
template <typename T>
Vector<T, 3> getTriangleNormal(const Vector<T, 3>& v0, const Vector<T, 3>& v1,
                               const Vector<T, 3>& v2) {
  Vector<T, 3> normal;
  TriangleNormal(normal, v0, v1, v2);
  return normal;
}

template <typename T>
T getTriangleArea(const Vector<T, 3>& v0, const Vector<T, 3>& v1,
                  const Vector<T, 3>& v2) {
  Vector<T, 3> e0 = v1 - v0;
  Vector<T, 3> e1 = v2 - v0;

  return T{0.5} * CrossProduct3(e0, e1).getnorm();
}
template <typename T>
T getTriangleArea(const Triangle<T>& triangle) {
  return getTriangleArea(triangle[0], triangle[1], triangle[2]);
}


}  // namespace offlat
