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

// geometry2d.h

#pragma once

#include "geometry/basic_geometry.h"
#include "data_struct/quadtree.h"

template <typename T>
class Circle final : public AABB<T, 2> {
 protected:
  T _Radius;

 public:
  Circle(T radius, const Vector<T, 2>& centre)
      : AABB<T, 2>(centre, 2 * radius, 2 * radius), _Radius(radius) {}
  bool isInside(const Vector<T, 2>& pt) const override {
    T eps = std::numeric_limits<T>::epsilon();
    return (pt - this->_center).getnorm2() <= _Radius * _Radius + eps;
  }
  template <typename S>
  bool IsInside(const Vector<S, 2>& pt) const {
    T eps = std::numeric_limits<decltype(T{} * S{})>::epsilon();
    return (pt - this->_center).getnorm2() <= _Radius * _Radius + eps;
  }
};


////////////////////////////////////////////////////////////
/// Geometry2D ////
///////////////////////////////////////////////////////////

template <typename T>
class Geometry2D : public AABB<T, 2> {
 private:
  // mesh data
  std::vector<Vector<T, 2>> _Voxels;
  // mesh param
  int _Nx;
  int _Ny;
  int _Nz = 1;
  int N;
  // it is recommended to set voxelSize = 1.
  T _voxelSize;
  // facilitate getting index in a cubic array
  // _idx (1, Nx)
  //      (i, j)
  // id = i + j * Nx
  Vector<int, 2> Projection;

  FlagField GeometryFlag;

  std::uint8_t _AABBflag;
  std::uint8_t _voidflag;

  // facilitate getting index
  // {0, 1, _Nx, -1, -_Nx, 1 + _Nx, -1 + _Nx, -1 - _Nx, 1 - _Nx}
  std::array<int, 9> Delta_Index{0,       1,        _Nx,      -1,     -_Nx,
                                 1 + _Nx, -1 + _Nx, -1 - _Nx, 1 - _Nx};
  Vector<T, 2> _MinCenter;

 public:
  // constructors
  // default aabbflag = 1, voidflag = 0
  Geometry2D(int Nx, int Ny, const AABB<T, 2>& AABBs, T voxelSize = T(1),
             const Vector<T, 2>& min = Vector<T, 2>{},
             std::uint8_t AABBflag = std::uint8_t(2),
             std::uint8_t voidflag = std::uint8_t(1));

  ~Geometry2D() = default;
  // read from AABBs
  void ReadAABBs(const AABB<T, 2>& AABBs, std::uint8_t AABBflag = std::uint8_t(2));

  // setup boundary
  template <typename LatSet>
  void SetupBoundary(std::uint8_t AABBflag = std::uint8_t(2),
                     std::uint8_t boundaryflag = std::uint8_t(4));

  // set flag in a region defined by AABBs
  void setFlag(const AABB<T, 2>& AABBs, std::uint8_t flag);
  void setFlag(const AABB<T, 2>& AABBs, std::uint8_t fromflag,
               std::uint8_t flag);
  // lambda function set flag
  // call: setFlag<LatSet>(AABBs, [](int id){func(id);});
  template <typename Func>
  void forEachVoxel(const AABB<T, 2>& AABBs, Func func);
  template <typename Func>
  void forEachVoxelint(const AABB<int, 2>& AABBs, Func func);
  template <typename Func>
  void forEachVoxel(const AABB<T, 2>& AABBs, std::uint8_t fromflag, Func func);
  template <typename Func>
  void forEachVoxel(std::uint8_t fromflag, Func func);
  // print
  void print() {
    std::cout << "[Geometry2D]: "
              << "\n"
              << "(Nx,Ny) = (" << _Nx << ", " << _Ny << "); "
              << "min = (" << this->_min[0] << ", " << this->_min[1] << "); "
              << "max = (" << this->_max[0] << ", " << this->_max[1] << ")"
              << std::endl;
  }

  // get
  int getNx() const { return _Nx; }
  int getNy() const { return _Ny; }
  int getNz() const { return _Nz; }
  // return N = Nx * Ny * Nz(1)
  int getVoxelsNum() const { return N; }
  T getVoxelSize() const { return _voxelSize; }
  const Vector<int, 2>& getProjection() const { return Projection; }
  int getIndex(int i, int j, int k) const { return i + j * _Nx; }
  int getIndex(const Vector<int, 2>& loc) const {
    return loc[0] + loc[1] * _Nx;
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
    return id + LatSet::c[dir] * Projection;
  }
  int getNbrId(int id, int dir) const { return id + Delta_Index[dir]; }

  const Vector<T, 2>& getMin() const { return this->_min; }
  const Vector<T, 2>& getMax() const { return this->_max; }

  std::vector<Vector<T, 2>>& getVoxels() { return _Voxels; }
  const std::vector<Vector<T, 2>>& getVoxels() const { return _Voxels; }
  Vector<T, 2>& getVoxel(int id) { return _Voxels[id]; }
  const Vector<T, 2>& getVoxel(int id) const { return _Voxels[id]; }
  FlagField& getGeoFlagField() { return GeometryFlag; }
  const FlagField& getGeoFlagField() const { return GeometryFlag; }
  
  std::uint8_t getGeoFlag(int id) const { return GeometryFlag.get(id); }
  std::uint8_t getAABBflag() const { return _AABBflag; }
  std::uint8_t getVoidflag() const { return _voidflag; }
  template <typename LatSet>
  bool hasNeighborFlag(int id, std::uint8_t flag) const {
    for (int dir = 1; dir < LatSet::q; dir++) {
      if (GeometryFlag.get(getNeighborId<LatSet>(id, dir)) == flag) return true;
    }
    return false;
  }
  int findIndex(const Vector<T, 2>& loc) const {
    int x = int((loc[0] - this->_min[0]) / _voxelSize);
    int y = int((loc[1] - this->_min[1]) / _voxelSize);
    return x + y * _Nx;
  }
  int findCellId(const Vector<T, 2>& pos) const {
    int i = int((pos[0] - this->_min[0]) / _voxelSize);
    int j = int((pos[1] - this->_min[1]) / _voxelSize);
    i = i < 0 ? 0 : i;
    j = j < 0 ? 0 : j;
    i = i > _Nx - 1 ? _Nx - 1 : i;
    j = j > _Ny - 1 ? _Ny - 1 : j;
    return i + j * _Nx;
  }
};


