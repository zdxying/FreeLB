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

// voxel_geometry.h
// Immutable, domain-sized voxel record: role, entity, material, boundary patches,
// and a narrow surface band.  It does not own a reader and has no IO dependency;
// VoxelGeometryBuilder is the only writer.  writeVoxelFlags converts the record
// into the active mask of a block field manager.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "data_struct/Vector.h"
#include "geometry/basic_geometry.h"

namespace flbgeo {

// Outside cells are not part of the lattice; Domain cells are active fluid;
// Solid cells are the obstacle.  The enum value is stable because geometry
// files will persist it (B1b).
enum class VoxelRole : std::uint8_t { Outside = 0, Domain = 1, Solid = 2 };

// A patch is an overlay on cells of one entity; it never changes occupancy.
enum class PatchKind : std::uint8_t { Wall, Inlet, Outlet, Symmetry };

inline constexpr std::uint16_t kNoEntity = std::uint16_t{0xFFFF};
inline constexpr std::size_t kInvalidCell = static_cast<std::size_t>(-1);

struct BoundaryPatch {
  std::uint16_t Entity{kNoEntity};
  PatchKind Kind{PatchKind::Wall};
  // Sorted cell ids; a patch selects each cell at most once.
  std::vector<std::size_t> Cells;
};

// Identifies the primitive that produced one entry of a build recipe.  The
// value is stable because the geometry file persists it.
enum class GeometrySourceKind : std::uint8_t { Stl = 0, Box = 1, Sphere = 2 };

// One construction step, kept so a case can be rebuilt from its geometry file
// instead of re-reading the original mesh.  Only the fields selected by Kind
// carry meaning; the rest keep their default values.
template <typename T>
struct GeometrySource {
  GeometrySourceKind Kind{GeometrySourceKind::Box};
  std::string Path;
  AABB<T, 3> Box;
  Vector<T, 3> Center;
  T Radius{};
  std::uint16_t Entity{kNoEntity};
  std::uint16_t Material{0};
  VoxelRole Role{VoxelRole::Domain};
  bool RequireClosed{false};
};

// The inputs of VoxelGeometryBuilder, in call order.
template <typename T>
struct GeometryRecipe {
  std::vector<GeometrySource<T>> Sources;
  T BandWidth{};
};

// Signed distance is negative inside the solid and positive in the fluid.
// The normal points from the solid towards the fluid.
template <typename T>
struct SurfaceSample {
  T Distance{};
  Vector<T, 3> Normal{};
  std::uint16_t Entity{kNoEntity};
};

template <typename T>
class VoxelGeometryBuilder;

// Defined in io/geometry_file.h; only the reader needs the private arrays.
template <typename T>
struct GeometryFile;

template <typename T>
class VoxelGeometry {
 public:
  // The domain and the grid must describe the same box: domain.extent ==
  // n * voxelSize within rounding, otherwise getPosition() would silently leave
  // the record's own domain.
  VoxelGeometry(int nx, int ny, int nz, const AABB<T, 3>& domain, T voxelSize);

  int getNx() const { return _nx; }
  int getNy() const { return _ny; }
  int getNz() const { return _nz; }
  std::size_t size() const { return _role.size(); }
  const AABB<T, 3>& getDomain() const { return _domain; }
  T getVoxelSize() const { return _voxelSize; }

  // Unchecked fast path: loc must come from getLoc() or a bounds-checked point.
  std::size_t getIndex(const Vector<int, 3>& loc) const;
  Vector<int, 3> getLoc(std::size_t id) const;
  // Center of the cell, not its lower corner.
  Vector<T, 3> getPosition(std::size_t id) const;
  // kInvalidCell when pos is outside the domain; points on the upper face are
  // outside because the cell's half-open interval ends there.
  std::size_t findCellId(const Vector<T, 3>& pos) const;

  VoxelRole getRole(std::size_t id) const { return _role[id]; }
  bool isDomain(std::size_t id) const { return _role[id] == VoxelRole::Domain; }
  bool isSolid(std::size_t id) const { return _role[id] == VoxelRole::Solid; }
  bool isOutside(std::size_t id) const { return _role[id] == VoxelRole::Outside; }
  std::uint16_t getEntity(std::size_t id) const { return _entity[id]; }
  std::uint16_t getMaterial(std::size_t id) const { return _material[id]; }

  std::size_t countRole(VoxelRole role) const;
  // Counts all cells carrying the entity, Domain and Solid alike.
  std::size_t countEntity(std::uint16_t entity) const;

  std::size_t getPatchCount() const { return _patches.size(); }
  const BoundaryPatch& getPatch(std::size_t i) const { return _patches[i]; }
  const std::vector<std::size_t>& getPatchCells(std::size_t i) const {
    return _patches[i].Cells;
  }

  bool hasSurfaceSample(std::size_t id) const {
    return _band.find(id) != _band.end();
  }
  // Only call after hasSurfaceSample(); the map lookup is the contract.
  const SurfaceSample<T>& getSurfaceSample(std::size_t id) const {
    return _band.at(id);
  }
  std::size_t getBandSize() const { return _bandIds.size(); }
  // Sorted, so iteration order is deterministic across builds and platforms.
  const std::vector<std::size_t>& getBandCells() const { return _bandIds; }

 private:
  friend class VoxelGeometryBuilder<T>;
  friend struct GeometryFile<T>;

  int _nx;
  int _ny;
  int _nz;
  AABB<T, 3> _domain;
  T _voxelSize;

  std::vector<VoxelRole> _role;
  std::vector<std::uint16_t> _entity;
  std::vector<std::uint16_t> _material;
  std::vector<BoundaryPatch> _patches;
  std::unordered_map<std::size_t, SurfaceSample<T>> _band;
  std::vector<std::size_t> _bandIds;
};

// Fills a 3D block field manager from the record: Domain cells get domainFlag,
// Solid and Outside get voidFlag.  Positions outside the record's domain are
// void.  FieldM is a template parameter so geometry does not include the field
// layer.
template <typename T, typename FieldM>
void writeVoxelFlags(const VoxelGeometry<T>& geo, FieldM& fields,
                     std::uint8_t domainFlag, std::uint8_t voidFlag);

}  // namespace flbgeo
