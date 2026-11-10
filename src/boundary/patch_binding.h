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

// patch_binding.h
// Binds the geometry record's boundary patches to the application's task flags:
// each PatchKind maps to one flag bit and the patch cells are OR-ed into the
// flag field.  Task dispatch and the flag-driven boundary managers then select
// those cells without an AABB painted per case.
//
// A kind without a mapping is counted (unmappedPatches/unmappedCells), never
// silently skipped: the caller decides whether that is an error for its setup.

#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>

#include "geometry/voxel_geometry.h"
#include "utils/util.h"

namespace flbbc {

// one flag bit per patch kind; 0 means "this kind has no boundary condition"
struct PatchFlagMap {
  std::uint8_t wall = 0;
  std::uint8_t inlet = 0;
  std::uint8_t outlet = 0;
  std::uint8_t symmetry = 0;
};

struct PatchBindStats {
  std::size_t wall = 0;
  std::size_t inlet = 0;
  std::size_t outlet = 0;
  std::size_t symmetry = 0;
  // a patch whose kind has no flag in the map
  std::size_t unmappedPatches = 0;
  std::size_t unmappedCells = 0;
  // band cells that do not touch a non-Domain cell: they belong to the narrow
  // band but are not wall cells, so no BC flag is painted on them
  std::size_t interiorCells = 0;
};

// A record cell is a boundary cell when at least one axis neighbour is not
// Domain (solid, outside, or beyond the grid).  The narrow band is wider than
// one cell, so a patch also contains cells that are not adjacent to the wall;
// painting a wall task on those would put a wall inside the fluid.
template <typename T>
bool isBoundaryCell(const flbgeo::VoxelGeometry<T>& record, std::size_t id) {
  const Vector<int, 3> loc = record.getLoc(id);
  const int dims[3] = {record.getNx(), record.getNy(), record.getNz()};
  for (int axis = 0; axis < 3; ++axis) {
    for (int sign = -1; sign <= 1; sign += 2) {
      Vector<int, 3> nloc = loc;
      nloc[axis] += sign;
      if (nloc[axis] < 0 || nloc[axis] >= dims[axis]) return true;
      if (!record.isDomain(record.getIndex(nloc))) return true;
    }
  }
  return false;
}

// OR the patch flags into every domain cell that the record maps to a patch
// cell.  `byEntity` (optional) overrides the kind mapping for one entity: two
// Wall patches can then drive different tasks (e.g. a static and a moving
// wall), which the kind alone cannot express.  The field manager must still be
// in sync with the record (the same geometry used by writeVoxelFlags); the
// mapping is by position, so halos are painted too and only cells carrying
// domainFlag are touched.
template <typename T, typename FieldM>
PatchBindStats bindPatches(
  const flbgeo::VoxelGeometry<T>& record, FieldM& flags,
  const PatchFlagMap& map, std::uint8_t domainFlag,
  const std::unordered_map<std::uint16_t, std::uint8_t>* byEntity = nullptr,
  std::uint8_t skipFlag = 0) {
  PatchBindStats stats;
  std::unordered_map<std::size_t, std::uint8_t> bitByCell;  // record id -> bits
  for (std::size_t p = 0; p < record.getPatchCount(); ++p) {
    const flbgeo::BoundaryPatch& patch = record.getPatch(p);
    std::uint8_t bit = 0;
    std::size_t* count = nullptr;
    switch (patch.Kind) {
      case flbgeo::PatchKind::Wall:
        bit = map.wall;
        count = &stats.wall;
        break;
      case flbgeo::PatchKind::Inlet:
        bit = map.inlet;
        count = &stats.inlet;
        break;
      case flbgeo::PatchKind::Outlet:
        bit = map.outlet;
        count = &stats.outlet;
        break;
      case flbgeo::PatchKind::Symmetry:
        bit = map.symmetry;
        count = &stats.symmetry;
        break;
    }
    if (byEntity != nullptr) {
      const auto overrideIt = byEntity->find(patch.Entity);
      if (overrideIt != byEntity->end()) bit = overrideIt->second;
    }
    if (bit == 0) {
      ++stats.unmappedPatches;
      stats.unmappedCells += patch.Cells.size();
      continue;
    }
    for (std::size_t id : patch.Cells) {
      if (!isBoundaryCell(record, id)) {
        ++stats.interiorCells;
        continue;
      }
      ++(*count);
      std::uint8_t& bits = bitByCell[id];
      bits = static_cast<std::uint8_t>(bits | bit);
    }
  }
  if (bitByCell.empty()) return stats;

  flags.forEach([&](auto& bf, std::size_t id) {
    const std::size_t rec = record.findCellId(bf.getBlock().getLoc_t(id));
    if (rec == flbgeo::kInvalidCell) return;
    const auto it = bitByCell.find(rec);
    if (it == bitByCell.end()) return;
    const std::uint8_t flag = bf.getField(0)[id];
    if (!util::isFlag(flag, domainFlag)) return;
    if (skipFlag != 0 && util::isFlag(flag, skipFlag)) return;
    // replace, not OR: a boundary flag is its own task key, as in the
    // existing AABB-painted cases; skipFlag keeps e.g. a periodic face intact
    bf.SetField(id, it->second);
  });
  return stats;
}

}  // namespace flbbc
