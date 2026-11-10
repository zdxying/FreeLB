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

// voxel_map.h
// map between active voxels and cells of the AABB, with owned/halo identity

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>

#ifdef _VOX_ENABLED

class VoxelMap {
 private:
  static constexpr std::int32_t Invalid = -1;

  // AABB cell number and active voxel number
  std::size_t AabbNum;
  std::size_t VoxNum;
  std::size_t OwnedNum;

  // AABB cell -> active slot (Invalid for void)
  std::int32_t* Aabb_to_Vox;
  // active slot -> AABB cell
  std::size_t* Vox_to_AABB;
  // active slot -> owned (1) or halo (0)
  std::uint8_t* Vox_Owned;

  bool _initialized;

  void release() {
    delete[] Aabb_to_Vox;
    delete[] Vox_to_AABB;
    delete[] Vox_Owned;
    Aabb_to_Vox = nullptr;
    Vox_to_AABB = nullptr;
    Vox_Owned = nullptr;
    AabbNum = 0;
    VoxNum = 0;
    OwnedNum = 0;
    _initialized = false;
  }

  void copyFrom(const VoxelMap& other) {
    AabbNum = other.AabbNum;
    VoxNum = other.VoxNum;
    OwnedNum = other.OwnedNum;
    _initialized = other._initialized;
    Aabb_to_Vox = nullptr;
    Vox_to_AABB = nullptr;
    Vox_Owned = nullptr;
    if (AabbNum > 0) {
      Aabb_to_Vox = new std::int32_t[AabbNum];
      std::copy(other.Aabb_to_Vox, other.Aabb_to_Vox + AabbNum, Aabb_to_Vox);
    }
    if (VoxNum > 0) {
      Vox_to_AABB = new std::size_t[VoxNum];
      Vox_Owned = new std::uint8_t[VoxNum];
      std::copy(other.Vox_to_AABB, other.Vox_to_AABB + VoxNum, Vox_to_AABB);
      std::copy(other.Vox_Owned, other.Vox_Owned + VoxNum, Vox_Owned);
    }
  }

 public:
  VoxelMap()
      : AabbNum(0), VoxNum(0), OwnedNum(0), Aabb_to_Vox(nullptr),
        Vox_to_AABB(nullptr), Vox_Owned(nullptr), _initialized(false) {}
  explicit VoxelMap(std::size_t voxnum)
      : AabbNum(0), VoxNum(0), OwnedNum(0), Aabb_to_Vox(nullptr),
        Vox_to_AABB(nullptr), Vox_Owned(nullptr), _initialized(false) {
    Init(voxnum);
  }

  VoxelMap(const VoxelMap& other)
      : AabbNum(0), VoxNum(0), OwnedNum(0), Aabb_to_Vox(nullptr),
        Vox_to_AABB(nullptr), Vox_Owned(nullptr), _initialized(false) {
    copyFrom(other);
  }

  VoxelMap& operator=(const VoxelMap& other) {
    if (this == &other) return *this;
    release();
    copyFrom(other);
    return *this;
  }

  VoxelMap(VoxelMap&& other) noexcept
      : AabbNum(other.AabbNum), VoxNum(other.VoxNum), OwnedNum(other.OwnedNum),
        Aabb_to_Vox(other.Aabb_to_Vox), Vox_to_AABB(other.Vox_to_AABB),
        Vox_Owned(other.Vox_Owned), _initialized(other._initialized) {
    other.AabbNum = 0;
    other.VoxNum = 0;
    other.OwnedNum = 0;
    other.Aabb_to_Vox = nullptr;
    other.Vox_to_AABB = nullptr;
    other.Vox_Owned = nullptr;
    other._initialized = false;
  }

  VoxelMap& operator=(VoxelMap&& other) noexcept {
    if (this == &other) return *this;
    release();
    AabbNum = other.AabbNum;
    VoxNum = other.VoxNum;
    OwnedNum = other.OwnedNum;
    Aabb_to_Vox = other.Aabb_to_Vox;
    Vox_to_AABB = other.Vox_to_AABB;
    Vox_Owned = other.Vox_Owned;
    _initialized = other._initialized;
    other.AabbNum = 0;
    other.VoxNum = 0;
    other.OwnedNum = 0;
    other.Aabb_to_Vox = nullptr;
    other.Vox_to_AABB = nullptr;
    other.Vox_Owned = nullptr;
    other._initialized = false;
    return *this;
  }

  ~VoxelMap() { release(); }

  // Rebuild both directions.  A built-but-empty map is a valid state.
  void Init(std::size_t aabbnum, std::size_t voxnum) {
    release();
    AabbNum = aabbnum;
    VoxNum = voxnum;
    if (AabbNum > 0) {
      Aabb_to_Vox = new std::int32_t[AabbNum];
      std::fill(Aabb_to_Vox, Aabb_to_Vox + AabbNum, Invalid);
    }
    if (VoxNum > 0) {
      Vox_to_AABB = new std::size_t[VoxNum]{};
      Vox_Owned = new std::uint8_t[VoxNum]{};
    }
    _initialized = true;
  }

  // Backward-compatible form for standalone users that only need the forward map.
  void Init(std::size_t voxnum) { Init(0, voxnum); }

  void reset() { release(); }

  bool isInitialized() const { return _initialized; }
  bool empty() const { return VoxNum == 0; }

  void setSlot(std::size_t slot, std::size_t aabbId, bool owned) {
    assert(slot < VoxNum);
    assert(AabbNum == 0 || aabbId < AabbNum);
    assert(aabbId <= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()));
    assert(AabbNum == 0 || Aabb_to_Vox[aabbId] == Invalid);
    Vox_to_AABB[slot] = aabbId;
    Vox_Owned[slot] = owned ? std::uint8_t{1} : std::uint8_t{0};
    if (AabbNum > 0) Aabb_to_Vox[aabbId] = static_cast<std::int32_t>(slot);
    if (owned) ++OwnedNum;
  }

  // get idx with a bounds check in debug builds
  std::size_t operator[](std::size_t VoxId) const {
    assert(VoxId < VoxNum);
    return Vox_to_AABB[VoxId];
  }

  bool isOwned(std::size_t slot) const {
    assert(slot < VoxNum);
    return Vox_Owned[slot] != 0;
  }

  bool isActive(std::size_t aabbId) const {
    assert(aabbId < AabbNum);
    return Aabb_to_Vox[aabbId] != Invalid;
  }

  std::size_t aabbToSlot(std::size_t aabbId) const {
    assert(aabbId < AabbNum);
    assert(Aabb_to_Vox[aabbId] != Invalid);
    return static_cast<std::size_t>(Aabb_to_Vox[aabbId]);
  }

  std::size_t* getVox_to_AABB() { return Vox_to_AABB; }
  const std::size_t* getVox_to_AABB() const { return Vox_to_AABB; }
  const std::uint8_t* getOwnedMask() const { return Vox_Owned; }

  std::size_t getAabbNum() const { return AabbNum; }
  std::size_t getVoxNum() const { return VoxNum; }
  std::size_t getOwnedNum() const { return OwnedNum; }
};

#endif
