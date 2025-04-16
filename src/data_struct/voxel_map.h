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
// map between valid voxel and cell of AABB
// this may be useful for complex geometry

#pragma once

#include <cstddef>
#include "utils/util.h"

#ifdef _VOX_ENABLED

class VoxelMap {
 private:
  // voxel number
  std::size_t VoxNum;

  // voxel to AABB cell: size = VoxNum
  std::size_t* Vox_to_AABB;

 public:
  VoxelMap() : VoxNum(0), Vox_to_AABB(nullptr) {}
  VoxelMap(std::size_t voxnum) : 
    VoxNum(voxnum), Vox_to_AABB(new std::size_t[voxnum]{}) {}

  ~VoxelMap(){ if (Vox_to_AABB) delete[] Vox_to_AABB; }

  void Init(std::size_t voxnum) {
    VoxNum = voxnum;
    Vox_to_AABB = new std::size_t[VoxNum]{};
  }

  // get idx without boundary check
  std::size_t operator[](std::size_t VoxId) const { return Vox_to_AABB[VoxId]; }

  std::size_t* getVox_to_AABB() { return Vox_to_AABB; }

  std::size_t getVoxNum() const { return VoxNum; }

};

#endif