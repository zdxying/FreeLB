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

// nbr_direction.h
// neighbor block direction relative to current block
//
// A geometry-level concept shared by geometry (BasicBlock::resize/getFace),
// boundary (periodic_boundary.h) and parallel (communicator.h).  Kept in a
// leaf header so neither geometry nor parallel headers have to depend on the other module for it.
#pragma once

// std::uint8_t
#include <cstdint>

enum NbrDirection : std::uint8_t {
  // not a neighbor
  NONE = 0,
  // x-negative
  XN = 1,
  // x-positive
  XP = 2,
  // y-negative
  YN = 4,
  // y-positive
  YP = 8,
  // z-negative
  ZN = 16,
  // z-positive
  ZP = 32
};

// an edge or corner neighbor carries one flag per crossed axis
constexpr NbrDirection operator|(NbrDirection a, NbrDirection b) {
  return static_cast<NbrDirection>(static_cast<std::uint8_t>(a) |
                                   static_cast<std::uint8_t>(b));
}
