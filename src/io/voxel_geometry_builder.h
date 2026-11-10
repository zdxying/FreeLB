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

// voxel_geometry_builder.h
// Rasterises STL meshes and primitives into a VoxelGeometry.  This lives in the
// IO layer because it consumes StlReader; the geometry record itself does not.

#pragma once

#include <cstdint>
#include <vector>

#include "geometry/voxel_geometry.h"

template <typename T>
class StlReader;

namespace flbgeo {

template <typename T>
class VoxelGeometryBuilder {
 public:
  // bandWidth is the half-width of the narrow surface band in physical units.
  VoxelGeometryBuilder(int nx, int ny, int nz, const AABB<T, 3>& domain, T voxelSize,
                       T bandWidth = T(1.5));

  // Rasterises reader into the grid at once, so the builder never keeps a
  // pointer to a reader that may be destroyed before build().  Cells outside
  // this geometry's domain are clipped.
  // requireClosed rejects an open or non-manifold mesh before any cell changes.
  // The reader's own voxel size controls the inside/outside classification;
  // construct it no coarser than this grid.  Band distances and normals are
  // computed from the triangles, independent of that sampling.
  // The reader is non-const because StlReader exposes getTree()/getMesh()
  // through non-const accessors.  Its signedDistance()/evalSurfaceNormal() do
  // not instantiate against the current Vector API, so the equivalent distance
  // and normal are computed from the mesh here.
  void addStl(StlReader<T>& reader, std::uint16_t entity, std::uint16_t material,
              VoxelRole insideRole, bool requireClosed = true);

  // Cell centers inside box receive role, subject to the same write priority as
  // addStl().  A domain box is usually added before the STLs.
  void addBox(const AABB<T, 3>& box, std::uint16_t entity, std::uint16_t material,
              VoxelRole role);

  void addSphere(const Vector<T, 3>& center, T radius, std::uint16_t entity,
                 std::uint16_t material, VoxelRole role);

  // Records a patch request; cells are selected from the built band when build()
  // runs, so this may be called before or after addStl().  The selector is an
  // AABB in physical coordinates.  A request that selects no cell aborts: a
  // silently empty inlet/outlet is worse than a failed build.
  void addPatch(std::uint16_t entity, PatchKind kind, const AABB<T, 3>& selector);

  // Sorts the band ids, applies all pending patches, and returns the geometry.
  VoxelGeometry<T> build();

 private:
  struct PendingPatch {
    std::uint16_t Entity;
    PatchKind Kind;
    AABB<T, 3> Selector;
  };

  // Solid always overwrites; Domain only fills Outside; Outside never writes.
  void writeVoxel(std::size_t id, VoxelRole role, std::uint16_t entity,
                  std::uint16_t material);

  static bool isInside(StlReader<T>& reader, const Vector<T, 3>& pos);

  // Inclusive cell-index range covering box, clamped to the grid.  Callers still
  // test cell centers because the range is deliberately a little generous.
  void cellRange(const AABB<T, 3>& box, Vector<int, 3>& lo,
                 Vector<int, 3>& hi) const;

  static void checkClosed(StlReader<T>& reader, T voxelSize,
                          std::uint16_t entity);

  VoxelGeometry<T> _geometry;
  T _bandWidth;
  std::vector<PendingPatch> _pendingPatches;
};

}  // namespace flbgeo
