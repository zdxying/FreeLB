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

// voxel_geometry_builder.hh

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <limits>
#include <utility>
#include <vector>

// stlreader.hh uses MPI_RANK in its print paths without including the macro
// header itself; include it first so this header also compiles on its own.
#include "utils/mpi_macros.h"
#include "io/stlreader.hh"
#include "io/voxel_geometry_builder.h"

namespace flbgeo {

template <typename T>
VoxelGeometryBuilder<T>::VoxelGeometryBuilder(int nx, int ny, int nz,
                                              const AABB<T, 3>& domain, T voxelSize,
                                              T bandWidth)
    : _geometry(nx, ny, nz, domain, voxelSize), _bandWidth(bandWidth) {
  if (!(_bandWidth > T(0))) {
    std::cerr << "[VoxelGeometryBuilder] bandWidth must be positive, got " << _bandWidth
              << std::endl;
    std::exit(1);
  }
}

template <typename T>
void VoxelGeometryBuilder<T>::addStl(StlReader<T>& reader, std::uint16_t entity,
                                     std::uint16_t material, VoxelRole insideRole,
                                     bool requireClosed) {
  if (reader.getMesh().triangleSize() == 0) {
    std::cerr << "[VoxelGeometryBuilder::addStl] entity " << entity
              << ": STL mesh has no triangles" << std::endl;
    std::exit(1);
  }
  if (requireClosed) checkClosed(reader, _geometry.getVoxelSize(), entity);

  const AABB<T, 3> meshBox(reader.getMin(), reader.getMax());

  // Occupancy: every cell whose center the octree marks as inside.  The octree
  // leaf grid is centered on the mesh and is not aligned with this geometry's
  // grid, so the point query, not a leaf-index copy, is the contract.
  Vector<int, 3> lo, hi;
  cellRange(meshBox, lo, hi);
  for (int k = lo[2]; k <= hi[2]; ++k) {
    for (int j = lo[1]; j <= hi[1]; ++j) {
      for (int i = lo[0]; i <= hi[0]; ++i) {
        const std::size_t id = _geometry.getIndex(Vector<int, 3>{i, j, k});
        const Vector<T, 3> pos = _geometry.getPosition(id);
        if (isInside(reader, pos)) writeVoxel(id, insideRole, entity, material);
      }
    }
  }

  // Surface band: mark a superset of the cells within bandWidth of the surface.
  // A boundary leaf contains surface, so any cell center within bandWidth of the
  // surface is inside that leaf's box expanded by bandWidth.  This keeps the
  // O(triangles) exact distance/normal evaluation off the interior and far field.
  std::vector<Octree<T>*> leaves;
  reader.getTree()->getLeafs(leaves);
  std::vector<std::uint8_t> candidate(_geometry.size(), std::uint8_t{0});
  for (Octree<T>* leaf : leaves) {
    if (!leaf->getBoundaryNode()) continue;
    const T reach = leaf->getRadius() + _bandWidth;
    const Vector<T, 3>& center = leaf->getCenter();
    const AABB<T, 3> leafBox(center - Vector<T, 3>(reach),
                             center + Vector<T, 3>(reach));
    Vector<int, 3> llo, lhi;
    cellRange(leafBox, llo, lhi);
    for (int k = llo[2]; k <= lhi[2]; ++k) {
      for (int j = llo[1]; j <= lhi[1]; ++j) {
        for (int i = llo[0]; i <= lhi[0]; ++i) {
          const std::size_t id = _geometry.getIndex(Vector<int, 3>{i, j, k});
          if (leafBox.isInside(_geometry.getPosition(id))) {
            candidate[id] = std::uint8_t{1};
          }
        }
      }
    }
  }

  for (std::size_t id = 0; id < candidate.size(); ++id) {
    if (!candidate[id]) continue;
    const Vector<T, 3> pos = _geometry.getPosition(id);
    const bool inside = isInside(reader, pos);
    const T distance = signedDistanceAt(reader, pos, inside);
    if (std::abs(distance) > _bandWidth) continue;
    const Vector<T, 3> normal = surfaceNormalAt(reader, pos, distance, inside);

    auto it = _geometry._band.find(id);
    if (it == _geometry._band.end() ||
        std::abs(distance) < std::abs(it->second.Distance)) {
      _geometry._band[id] = SurfaceSample<T>{distance, normal, entity};
    }
  }
}

template <typename T>
void VoxelGeometryBuilder<T>::addBox(const AABB<T, 3>& box, std::uint16_t entity,
                                     std::uint16_t material, VoxelRole role) {
  Vector<int, 3> lo, hi;
  cellRange(box, lo, hi);
  for (int k = lo[2]; k <= hi[2]; ++k) {
    for (int j = lo[1]; j <= hi[1]; ++j) {
      for (int i = lo[0]; i <= hi[0]; ++i) {
        const std::size_t id = _geometry.getIndex(Vector<int, 3>{i, j, k});
        if (box.isInside(_geometry.getPosition(id))) {
          writeVoxel(id, role, entity, material);
        }
      }
    }
  }
}

template <typename T>
void VoxelGeometryBuilder<T>::addSphere(const Vector<T, 3>& center, T radius,
                                        std::uint16_t entity,
                                        std::uint16_t material, VoxelRole role) {
  if (!(radius > T(0))) {
    std::cerr << "[VoxelGeometryBuilder::addSphere] radius must be positive, got "
              << radius << std::endl;
    std::exit(1);
  }
  const Vector<T, 3> half(radius);
  const AABB<T, 3> sphereBox(center - half, center + half);

  Vector<int, 3> lo, hi;
  cellRange(sphereBox, lo, hi);
  const T radius2 = radius * radius;
  for (int k = lo[2]; k <= hi[2]; ++k) {
    for (int j = lo[1]; j <= hi[1]; ++j) {
      for (int i = lo[0]; i <= hi[0]; ++i) {
        const std::size_t id = _geometry.getIndex(Vector<int, 3>{i, j, k});
        const Vector<T, 3> offset = _geometry.getPosition(id) - center;
        if (offset.getnorm2() <= radius2) {
          writeVoxel(id, role, entity, material);
        }
      }
    }
  }
}

template <typename T>
void VoxelGeometryBuilder<T>::addPatch(std::uint16_t entity, PatchKind kind,
                                       const AABB<T, 3>& selector) {
  _pendingPatches.push_back(PendingPatch{entity, kind, selector});
}

template <typename T>
VoxelGeometry<T> VoxelGeometryBuilder<T>::build() {
  std::vector<std::size_t> ids;
  ids.reserve(_geometry._band.size());
  for (const auto& entry : _geometry._band) ids.push_back(entry.first);
  std::sort(ids.begin(), ids.end());
  _geometry._bandIds = std::move(ids);

  _geometry._patches.clear();
  _geometry._patches.reserve(_pendingPatches.size());
  for (const PendingPatch& pending : _pendingPatches) {
    BoundaryPatch patch;
    patch.Entity = pending.Entity;
    patch.Kind = pending.Kind;
    patch.Cells.reserve(16);
    for (std::size_t id : _geometry._bandIds) {
      if (!_geometry.isDomain(id)) continue;
      if (_geometry.getSurfaceSample(id).Entity != pending.Entity) continue;
      if (!pending.Selector.isInside(_geometry.getPosition(id))) continue;
      patch.Cells.push_back(id);
    }
    std::sort(patch.Cells.begin(), patch.Cells.end());
    patch.Cells.erase(std::unique(patch.Cells.begin(), patch.Cells.end()),
                      patch.Cells.end());

    if (patch.Cells.empty()) {
      std::cerr << "[VoxelGeometryBuilder::build] patch for entity " << pending.Entity
                << " selected no Domain cell in the surface band; check the selector"
                << std::endl;
      std::exit(1);
    }
    _geometry._patches.push_back(std::move(patch));
  }

  return _geometry;
}

template <typename T>
void VoxelGeometryBuilder<T>::writeVoxel(std::size_t id, VoxelRole role,
                                         std::uint16_t entity,
                                         std::uint16_t material) {
  switch (role) {
    case VoxelRole::Solid:
      // A solid always wins: an embedded obstacle must not be erased by a later
      // domain box, and its entity/material must survive.
      _geometry._role[id] = VoxelRole::Solid;
      _geometry._entity[id] = entity;
      _geometry._material[id] = material;
      break;
    case VoxelRole::Domain:
      if (_geometry._role[id] == VoxelRole::Outside) {
        _geometry._role[id] = VoxelRole::Domain;
        _geometry._entity[id] = entity;
        _geometry._material[id] = material;
      }
      break;
    case VoxelRole::Outside:
      break;
  }
}

template <typename T>
void VoxelGeometryBuilder<T>::cellRange(const AABB<T, 3>& box, Vector<int, 3>& lo,
                                        Vector<int, 3>& hi) const {
  const Vector<T, 3>& min = _geometry.getDomain().getMin();
  const T h = _geometry.getVoxelSize();
  const int n[3] = {_geometry.getNx(), _geometry.getNy(), _geometry.getNz()};
  for (int i = 0; i < 3; ++i) {
    const T a = (box.getMin()[i] - min[i]) / h;
    const T b = (box.getMax()[i] - min[i]) / h;
    lo[i] = std::max(0, static_cast<int>(std::floor(a)));
    hi[i] = std::min(n[i] - 1, static_cast<int>(std::ceil(b)));
  }
}

template <typename T>
bool VoxelGeometryBuilder<T>::isInside(StlReader<T>& reader,
                                       const Vector<T, 3>& pos) {
  Octree<T>* node = reader.getTree()->find(pos);
  return node != nullptr && node->isLeaf() && node->getInside();
}

template <typename T>
T VoxelGeometryBuilder<T>::signedDistanceAt(StlReader<T>& reader,
                                            const Vector<T, 3>& pos, bool inside) {
  T distance = std::numeric_limits<T>::max();
  for (const Triangle<T>& triangle : reader.getMesh().getTriangles()) {
    const Vector<T, 3> point = triangle.closestPtPointTriangle(pos);
    distance = std::min(distance, (point - pos).getnorm());
  }
  return inside ? -distance : distance;
}

template <typename T>
Vector<T, 3> VoxelGeometryBuilder<T>::surfaceNormalAt(StlReader<T>& reader,
                                                      const Vector<T, 3>& pos,
                                                      T distance,
                                                      bool inside) const {
  const std::vector<Triangle<T>>& triangles = reader.getMesh().getTriangles();
  T best2 = std::numeric_limits<T>::max();
  Vector<T, 3> closest;
  const Triangle<T>* nearest = nullptr;
  for (const Triangle<T>& triangle : triangles) {
    const Vector<T, 3> point = triangle.closestPtPointTriangle(pos);
    const T d2 = (point - pos).getnorm2();
    if (d2 < best2) {
      best2 = d2;
      closest = point;
      nearest = &triangle;
    }
  }

  constexpr T eps = std::numeric_limits<T>::epsilon();
  Vector<T, 3> normal = closest - pos;
  const T length = normal.getnorm();
  if (!(length > eps)) {
    // A point exactly on the surface has no direction to the closest point;
    // average the triangles that contain it, as findNormalOnSurface() does.
    Vector<T, 3> sum;
    std::size_t count = 0;
    for (const Triangle<T>& triangle : triangles) {
      if (triangle.IsInside(pos)) {
        sum = sum + triangle.getNormal();
        ++count;
      }
    }
    const T sumNorm = sum.getnorm();
    if (count != 0 && sumNorm > eps) {
      sum /= sumNorm;
      return sum;
    }
    if (nearest != nullptr) return nearest->getNormal();
    return Vector<T, 3>{};
  }

  normal /= length;
  // Near the surface the octree's leaf-scale classification can disagree with
  // the nearest-triangle side, so probe both sides at the reader's sampling
  // scale before choosing the sign.
  const T h = reader.getVoxelSize();
  if (std::abs(distance) < h) {
    const bool inner = isInside(reader, closest - normal * h);
    const bool outer = isInside(reader, closest + normal * h);
    if (!(inner && !outer)) normal = normal * T(-1);
  } else if (!inside) {
    normal = normal * T(-1);
  }
  return normal;
}

template <typename T>
void VoxelGeometryBuilder<T>::checkClosed(StlReader<T>& reader, T voxelSize,
                                          std::uint16_t entity) {
  const std::vector<Triangle<T>>& triangles = reader.getMesh().getTriangles();
  if (triangles.empty()) {
    std::cerr << "[VoxelGeometryBuilder::addStl] entity " << entity
              << ": STL mesh has no triangles" << std::endl;
    std::exit(1);
  }

  // Quantise before comparing so vertices that were written with different
  // round-off collapse to one key.  The step is tied to the voxel size because
  // features smaller than that cannot be represented by this grid anyway.
  const T step = std::max(T(1e-9), voxelSize * T(1e-6));
  std::map<std::array<long long, 3>, int> vertexIds;
  std::map<std::pair<int, int>, int> edgeCount;

  auto vertexId = [&](const Vector<T, 3>& v) {
    std::array<long long, 3> key;
    for (int i = 0; i < 3; ++i) {
      key[i] = std::llround(static_cast<double>(v[i]) / static_cast<double>(step));
    }
    auto found = vertexIds.find(key);
    if (found != vertexIds.end()) return found->second;
    const int id = static_cast<int>(vertexIds.size());
    vertexIds.emplace(key, id);
    return id;
  };

  for (const Triangle<T>& triangle : triangles) {
    int ids[3];
    for (int i = 0; i < 3; ++i) ids[i] = vertexId(triangle.vertex[i]);
    for (int i = 0; i < 3; ++i) {
      const int a = ids[i];
      const int b = ids[(i + 1) % 3];
      ++edgeCount[std::minmax(a, b)];
    }
  }

  std::size_t badEdges = 0;
  for (const auto& entry : edgeCount) {
    if (entry.second != 2) ++badEdges;
  }
  if (badEdges != 0) {
    std::cerr << "[VoxelGeometryBuilder::addStl] entity " << entity
              << ": mesh is not closed; " << badEdges << " of " << edgeCount.size()
              << " undirected edges are not shared by exactly two triangles"
              << std::endl;
    std::exit(1);
  }
}

}  // namespace flbgeo
