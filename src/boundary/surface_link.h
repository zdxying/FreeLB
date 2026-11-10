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
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// surface_link.h
// One fluid-to-wall link of the narrow-band surface description.  The table is
// derived data: it is built from VoxelGeometry and the block flag field, and
// never written back.  Direction is the lattice direction from the fluid cell
// towards the wall; Fraction is the wall distance divided by the link length.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "data_struct/Vector.h"

namespace flbbc {

inline constexpr std::uint16_t kNoEntity = std::uint16_t{0xFFFF};

template <typename T>
struct SurfaceLink {
  // Local cell id in the block's extended lattice.
  std::size_t CellId{};
  // Lattice direction from the fluid cell towards the wall.
  unsigned int Direction{};
  // Wall position along the link in (0, 1].
  T Fraction{};
  std::uint16_t Entity{kNoEntity};
  // Zero for a static wall; the moving-wall branch fills this per entity.
  Vector<T, 3> WallVelocity{};
  // The fluid node behind the wall is outside the block or void; use the
  // one-sided q >= 1/2 branch instead of reading a second fluid node.
  bool Thin{false};
  // The band normal did not define a wall crossing on this link; the builder
  // fell back to the half-way position.
  bool Fallback{false};
};

// Per-block link lists.  A block owns one entry per boundary cell/direction;
// the wall itself may be in another block, but the void flag is in the halo.
template <typename T>
struct SurfaceLinkTable {
  std::vector<std::vector<SurfaceLink<T>>> Blocks;
};

}  // namespace flbbc
