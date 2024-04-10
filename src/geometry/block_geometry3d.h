/* This file is part of FreeLB
 * 
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 * 
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 * 
 */

// block_geometry3d.h

#pragma once

#include "geometry/geometry3d.h"

template <typename T>
class Block3D : public AABB<T, 3> {
 private:
  AABB<T,3> _OutAABB;
  // int aabb
  AABB<int, 3> _OutBlock;
  AABB<int, 3> _Block;
  // ref level
  std::uint8_t _level;
  // mesh data
  std::vector<Vector<T, 3>> _Voxels;
  // abstract tree
  std::vector<AbstractTree<T, 3>> _Trees;
  // mesh param
  int _Nx;
  int _Ny;
  int _Nz;
  int N;
  // it is recommended to set voxelSize = 1.
  T _voxelSize;
  // facilitate getting index in a cubic array
  // _idx (1, Nx, Nx * Ny)
  //      (i, j, k)
  // id = i + j * Nx + k * Nx * Ny
  Vector<int, 3> Projection;

  FlagField GeometryFlag;

  std::uint8_t _AABBflag;
  std::uint8_t _voidflag;

  // facilitate getting index in a cubic array
  // {0, 1, -1, Nx, -Nx, Nx * Ny, -Nx * Ny, 
  //  1 + Nx, -1 - Nx, 1 + Nx*Ny, -1 - Nx*Ny, Nx + Nx*Ny, -Nx - Nx*Ny,
  // 1 - Nx, -1 + Nx, 1 - Nx*Ny, -1 + Nx*Ny, Nx - Nx*Ny, -Nx + Nx*Ny,
  // 1 + Nx + Nx*Ny, -1 - Nx - Nx*Ny, 1 + Nx - Nx*Ny, -1 - Nx + Nx*Ny,
  // 1 - Nx + Nx*Ny, -1 + Nx - Nx*Ny, -1 + Nx + Nx*Ny, 1 - Nx - Nx*Ny}
  std::array<int, 27> Delta_Index{
      0,  1,  -1, _Nx, -_Nx, _Nx*_Ny, -_Nx*_Ny, 
      1+_Nx, -1-_Nx, 1+_Nx*_Ny, -1-_Nx*_Ny, _Nx+_Nx*_Ny, -_Nx-_Nx*_Ny, 
      1-_Nx, -1+_Nx, 1-_Nx*_Ny, -1+_Nx*_Ny, _Nx-_Nx*_Ny, -_Nx+_Nx*_Ny, 
      1+_Nx+_Nx*_Ny, -1-_Nx-_Nx*_Ny, 1+_Nx-_Nx*_Ny, -1-_Nx+_Nx*_Ny, 
      1-_Nx+_Nx*_Ny, -1+_Nx-_Nx*_Ny, -1+_Nx+_Nx*_Ny, 1-_Nx-_Nx*_Ny};
  
  Vector<T, 3> _MinCenter;

 public:
  // constructors
  // default aabbflag = 1, voidflag = 0
  Block3D(const AABB<T, 3>& AABBs, const Vector<int, 3>& idxAABBs,
          std::uint8_t level, T voxelSize = T(1),
          std::uint8_t AABBflag = std::uint8_t(2),
          std::uint8_t voidflag = std::uint8_t(1));
  Block3D(const AABB<T, 3>& AABBs, const Vector<int, 3>& idxAABBs,
          T voxelSize = T(1), std::uint8_t AABBflag = std::uint8_t(2),
          std::uint8_t voidflag = std::uint8_t(1));

  ~Block3D() = default;
  // read from AABBs
  void ReadAABBs(const AABB<T, 3>& AABBs, std::uint8_t AABBflag = std::uint8_t(2));

  // set flag in a region defined by AABBs
  void setFlag(const AABB<T, 3>& AABBs, std::uint8_t flag);
  void setFlag(const AABB<T, 3>& AABBs, std::uint8_t fromflag,
               std::uint8_t flag);
  // lambda function set flag
  // call: setFlag<LatSet>(AABBs, [](int id){func(id);});
  template <typename Func>
  void forEachVoxel(const AABB<T, 3>& AABBs, Func func);
  template <typename Func>
  void forEachVoxel(const AABB<T, 3>& AABBs, std::uint8_t fromflag, Func func);
  template <typename Func>
  void forEachVoxel(std::uint8_t fromflag, Func func);
  // print

  // get
  int getNx() const { return _Nx; }
  int getNy() const { return _Ny; }
  int getNz() const { return _Nz; }
  // return N = Nx * Ny * Nz(1)
  int getN() const { return N; }
  T getVoxelSize() const { return _voxelSize; }
  const Vector<int, 3>& getProjection() const { return Projection; }
  int getIndex(int i, int j, int k) const {     return i + j * Projection[1] + k * Projection[2]; }
  int getIndex(const Vector<int, 3>& loc) const {
    return loc[0] + loc[1] * Projection[1] + loc[2] * Projection[2];
  }
  int getGridx(const T loc) const {
    return int((loc - this->_min[0]) / _voxelSize);
  }
  int getGridy(const T loc) const {
    return int((loc - this->_min[1]) / _voxelSize);
  }
  int getGridz(const T loc) const { return 0; }

  template <typename LatSet>
  int getNeighborId(int id, int dir) const {
    return id + LatSet::c[dir][0] + LatSet::c[dir][1] * Projection[1] +
           LatSet::c[dir][2] * Projection[2];
  }
  int getNbrId(int id, int dir) const { return id + Delta_Index[dir]; }

  const Vector<T, 3>& getMin() const { return this->_min; }
  const Vector<T, 3>& getMax() const { return this->_max; }

  const AABB<T,3>& getAABB() const { return *this; }
  const AABB<T,3>& getOutAABB() const { return _OutAABB; }
  const AABB<int,3>& getBlock() const { return _Block; }
  const AABB<int,3>& getOutBlock() const { return _OutBlock; }
  
  std::vector<Vector<T, 3>>& getVoxels() { return _Voxels; }
  const std::vector<Vector<T, 2>>& getVoxels() const { return _Voxels; }
  Vector<T, 3>& getVoxel(int id) { return _Voxels[id]; }
  const Vector<T, 3>& getVoxel(int id) const { return _Voxels[id]; }
  FlagField& getGeoFlagField() { return GeometryFlag; }
  std::uint8_t getGeoFlag(int id) const { return GeometryFlag.get(id); }
  std::uint8_t getAABBflag() const { return _AABBflag; }
  std::uint8_t getVoidflag() const { return _voidflag; }
  template <typename LatSet>
  bool hasNeighborFlag(int id, std::uint8_t flag) const {
    for (int dir = 1; dir < LatSet::q; ++dir) {
      if (GeometryFlag.get(getNbrId(id, dir)) == flag) return true;
    }
    return false;
  }
  int findIndex(const Vector<T, 3>& loc) const {
    int x = int((loc[0] - this->_min[0]) / _voxelSize);
    int y = int((loc[1] - this->_min[1]) / _voxelSize);
    int z = int((loc[2] - this->_min[2]) / _voxelSize);
    return x + y * Projection[1] + z * Projection[2];
  }
  int findCellId(const Vector<T, 3>& pos) const {
    int i = int((pos[0] - this->_min[0]) / _voxelSize);
    int j = int((pos[1] - this->_min[1]) / _voxelSize);
    int k = int((pos[2] - this->_min[2]) / _voxelSize);
    i = i < 0 ? 0 : i;
    j = j < 0 ? 0 : j;
    k = k < 0 ? 0 : k;
    i = i > _Nx - 1 ? _Nx - 1 : i;
    j = j > _Ny - 1 ? _Ny - 1 : j;
    k = k > _Nz - 1 ? _Nz - 1 : k;
    return i + j * Projection[1] + k * Projection[2];
  }
};

template <typename T>
class BlockGeometry3D : public AABB<T, 3> {};

template <typename T>
class RefinedGeometry3D{
  
};