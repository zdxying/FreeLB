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

// link_key.h

#pragma once

#include <cstdint>
#include <tuple>

namespace flbcomm {

// Total order of one link inside a coalesced message.  The sender and the
// receiver derive it from the geometry alone, so the two ref lists can be
// packed in the same order without exchanging it: a Dir link carries the
// receiver's face/corner/edge split, while the Full/All and cross-level
// exchanges have one link per block pair and leave Kind/Index at zero.
struct LinkKey {
  std::int32_t Receiver{};
  std::int32_t Sender{};
  std::uint16_t Kind{};
  std::uint16_t Index{};

  bool operator<(const LinkKey& other) const {
    return std::tie(Receiver, Sender, Kind, Index) <
           std::tie(other.Receiver, other.Sender, other.Kind, other.Index);
  }
};

}  // namespace flbcomm
