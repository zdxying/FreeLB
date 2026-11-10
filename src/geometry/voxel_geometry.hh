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

// voxel_geometry.hh

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "geometry/voxel_geometry.h"

namespace flbgeo {

template <typename T>
VoxelGeometry<T>::VoxelGeometry(int nx, int ny, int nz, const AABB<T, 3>& domain,
                                T voxelSize)
    : _nx(nx), _ny(ny), _nz(nz), _domain(domain), _voxelSize(voxelSize) {
  if (nx <= 0 || ny <= 0 || nz <= 0) {
    std::cerr << "[VoxelGeometry] grid size must be positive, got " << nx << "x" << ny
              << "x" << nz << std::endl;
    std::exit(1);
  }
  if (!(voxelSize > T(0))) {
    std::cerr << "[VoxelGeometry] voxelSize must be positive, got " << voxelSize
              << std::endl;
    std::exit(1);
  }

  const Vector<T, 3>& extent = domain.getExtension();
  for (int i = 0; i < 3; ++i) {
    const T n = static_cast<T>(i == 0 ? nx : (i == 1 ? ny : nz));
    const T expected = n * voxelSize;
    const T tol = T(8) * std::numeric_limits<T>::epsilon() *
                  std::max(T(1), std::max(extent[i], expected));
    if (!std::isfinite(extent[i]) || std::abs(extent[i] - expected) > tol) {
      std::cerr << "[VoxelGeometry] domain extent " << i << " = " << extent[i]
                << " does not match n * voxelSize = " << expected << std::endl;
      std::exit(1);
    }
  }

  const std::size_t total =
    static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
    static_cast<std::size_t>(nz);
  _role.assign(total, VoxelRole::Outside);
  _entity.assign(total, kNoEntity);
  _material.assign(total, std::uint16_t{0});
}

template <typename T>
std::size_t VoxelGeometry<T>::getIndex(const Vector<int, 3>& loc) const {
  return (static_cast<std::size_t>(loc[2]) * static_cast<std::size_t>(_ny) +
          static_cast<std::size_t>(loc[1])) *
           static_cast<std::size_t>(_nx) +
         static_cast<std::size_t>(loc[0]);
}

template <typename T>
Vector<int, 3> VoxelGeometry<T>::getLoc(std::size_t id) const {
  const std::size_t nx = static_cast<std::size_t>(_nx);
  const std::size_t ny = static_cast<std::size_t>(_ny);
  return Vector<int, 3>{static_cast<int>(id % nx),
                        static_cast<int>((id / nx) % ny),
                        static_cast<int>(id / (nx * ny))};
}

template <typename T>
Vector<T, 3> VoxelGeometry<T>::getPosition(std::size_t id) const {
  const Vector<int, 3> loc = getLoc(id);
  const Vector<T, 3>& min = _domain.getMin();
  return Vector<T, 3>{min[0] + (static_cast<T>(loc[0]) + T(0.5)) * _voxelSize,
                      min[1] + (static_cast<T>(loc[1]) + T(0.5)) * _voxelSize,
                      min[2] + (static_cast<T>(loc[2]) + T(0.5)) * _voxelSize};
}

template <typename T>
std::size_t VoxelGeometry<T>::findCellId(const Vector<T, 3>& pos) const {
  const Vector<T, 3>& min = _domain.getMin();
  Vector<int, 3> loc;
  for (int i = 0; i < 3; ++i) {
    const T rel = (pos[i] - min[i]) / _voxelSize;
    const int n = i == 0 ? _nx : (i == 1 ? _ny : _nz);
    if (!(rel >= T(0)) || !(rel < static_cast<T>(n))) return kInvalidCell;
    loc[i] = static_cast<int>(rel);
  }
  return getIndex(loc);
}

template <typename T>
std::size_t VoxelGeometry<T>::countRole(VoxelRole role) const {
  return static_cast<std::size_t>(std::count(_role.begin(), _role.end(), role));
}

template <typename T>
std::size_t VoxelGeometry<T>::countEntity(std::uint16_t entity) const {
  return static_cast<std::size_t>(
    std::count(_entity.begin(), _entity.end(), entity));
}

template <typename T, typename FieldM>
void writeVoxelFlags(const VoxelGeometry<T>& geo, FieldM& fields,
                     std::uint8_t domainFlag, std::uint8_t voidFlag) {
  fields.forEach([&](auto& blockField, std::size_t id) {
    const Vector<T, 3> pos = blockField.getBlock().getLoc_t(id);
    const std::size_t gid = geo.findCellId(pos);
    blockField.SetField(
      id, gid != kInvalidCell && geo.isDomain(gid) ? domainFlag : voidFlag);
  });
}

}  // namespace flbgeo
