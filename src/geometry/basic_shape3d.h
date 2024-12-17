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

#include "geometry/basic_geometry.h"

template <typename T>
class Cylinder final : public AABB<T, 3> {
 protected:
  T _Radius;
  T _Height;
  // normal of the circle
  Vector<T, 3> _height;
  // centre of the bottom circle of the cylinder
  Vector<T, 3> _centre;

  // normalized _height
  Vector<T, 3> _heightn;

 public:
  Cylinder(T radius, const Vector<T, 3>& height, const Vector<T, 3>& centre)
      :  AABB<T, 3>(), _Radius(radius), _height(height), _centre(centre) {
    // get radius and height
    _Height = _height.getnorm();
    // get normalized _height, i.e. get the normal of the circle
    _heightn = _height.getnormalize();
    // get projection of square of the circle's diameter
    // centre of the another circle
    Vector<T, 3> centre2 = _centre + _height;
    T x_min = _centre[0] - _Radius * std::sqrt(1 - _heightn[0] * _heightn[0]);
    T x_max = centre2[0] + _Radius * std::sqrt(1 - _heightn[0] * _heightn[0]);
    T y_min = _centre[1] - _Radius * std::sqrt(1 - _heightn[1] * _heightn[1]);
    T y_max = centre2[1] + _Radius * std::sqrt(1 - _heightn[1] * _heightn[1]);
    T z_min = _centre[2] - _Radius * std::sqrt(1 - _heightn[2] * _heightn[2]);
    T z_max = centre2[2] + _Radius * std::sqrt(1 - _heightn[2] * _heightn[2]);
    // extension
    Vector<T, 3> extension(x_max - x_min, y_max - y_min, z_max - z_min);
    // get AABB centre
    Vector<T, 3> centre_ = _centre + _height / T(2);
    // set AABB
    AABB<T, 3>::_center = centre_;
    AABB<T, 3>::_extension = extension;
    AABB<T, 3>::_min = centre_ - extension / T(2);
    AABB<T, 3>::_max = centre_ + extension / T(2);
  }
  bool isInside(const Vector<T, 3>& pt) const override {
    constexpr T eps = std::numeric_limits<T>::epsilon();
    // check if pt is inside the Cylinder
    const Vector<T, 3> pt_ = pt - _centre;
    // projection of pt_ on _height
    const T proj_h = pt_ * _heightn;
    if (proj_h < 0 - eps || proj_h > _Height + eps) return false;
    const T pt_norm = pt_.getnorm();
    // get cos(theta) between pt_ and _height
    const T cos_theta = proj_h / pt_norm;
    // projection of pt_ on _Radius
    const T proj_r = pt_norm * std::sqrt(T(1) - cos_theta * cos_theta);
    if (proj_r > _Radius + eps) return false;
    return true;
  }
};