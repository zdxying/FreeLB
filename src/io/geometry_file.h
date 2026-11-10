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

// geometry_file.h
// Binary persistence for a VoxelGeometry and the recipe that built it.  The
// format is native-endian and versioned; a mismatched magic, version, scalar
// type, or digest is a hard error because a partially trusted geometry would
// silently change the simulation.

#pragma once

#include <string>

#include "geometry/voxel_geometry.h"

namespace flbgeo {

// Geometry and its construction recipe are stored together so a case can be
// rebuilt and audited from one file.
template <typename T>
struct GeometryFile {
  VoxelGeometry<T> Geometry;
  GeometryRecipe<T> Recipe;

 private:
  static GeometryFile<T> readFromFile(const std::string& path);

  template <typename U>
  friend GeometryFile<U> readGeometryFile(const std::string& path);
};

template <typename T>
void writeGeometryFile(const std::string& path, const VoxelGeometry<T>& geometry,
                       const GeometryRecipe<T>& recipe);

template <typename T>
GeometryFile<T> readGeometryFile(const std::string& path);

}  // namespace flbgeo
