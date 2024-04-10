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

// geometry3d.h

#pragma once

#include "geometry/basic_geometry.h"
#include "data_struct/field.h"
#include "io/stlreader.h"

template <typename T>
class Cylinder final : public AABB<T, 3> {
 protected:
  T _Radius;
  T _Height;
  // normal of the circle
  Vector<T, 3> _height;
  // centre of the bottom circle of the cylinder
  Vector<T, 3> _centre;

  // normalized _height
  Vector<T, 3> _heightn;

 public:
  Cylinder(T radius, const Vector<T, 3>& height, const Vector<T, 3>& centre)
      : _Radius(radius), _height(height), _centre(centre), AABB<T, 3>() {
    // get radius and height
    _Height = _height.getnorm();
    // get normalized _height, i.e. get the normal of the circle
    _heightn = _height.getnormalize();
    // get projection of square of the circle's diameter
    // centre of the another circle
    Vector<T, 3> centre2 = _centre + _height;
    T x_min = _centre[0] - _Radius * sqrt(1 - pow(_heightn[0], 2));
    T x_max = centre2[0] + _Radius * sqrt(1 - pow(_heightn[0], 2));
    T y_min = _centre[1] - _Radius * sqrt(1 - pow(_heightn[1], 2));
    T y_max = centre2[1] + _Radius * sqrt(1 - pow(_heightn[1], 2));
    T z_min = _centre[2] - _Radius * sqrt(1 - pow(_heightn[2], 2));
    T z_max = centre2[2] + _Radius * sqrt(1 - pow(_heightn[2], 2));
    // extension
    Vector<T, 3> extension(x_max - x_min, y_max - y_min, z_max - z_min);
    // get AABB centre
    Vector<T, 3> centre_ = _centre + _height / T(2);
    // set AABB
    AABB<T, 3>::_center = centre_;
    AABB<T, 3>::_extension = extension;
    AABB<T, 3>::_min = centre - extension / T(2);
    AABB<T, 3>::_max = centre + extension / T(2);
  }
  bool isInside(const Vector<T, 3>& pt) const override {
    T eps = std::numeric_limits<T>::epsilon();
    // check if pt is inside the Cylinder
    Vector<T, 3> pt_ = pt - _centre;
    // projection of pt_ on _height
    T proj_h = pt_ * _heightn;
    if (proj_h < 0 - eps || proj_h > _Height + eps) return false;
    T pt_norm = pt_.getnorm();
    // get cos(theta) between pt_ and _height
    // T cos_theta = proj_h / pt_norm;
    // projection of pt_ on _Radius
    T proj_r = pt_norm * sqrt(1 - pow(proj_h / pt_norm, 2));
    if (proj_r > _Radius + eps) return false;
    return true;
  }
};

////////////////////////
/// Geometry3D ////
////////////////////////

// in Geometry3D, voxel is initialized in each cell in a cubic array defined by
// the AABB regardless of whether the cell is inside the geometry or not

template <typename T>
class Geometry3D : public AABB<T, 3> {
 private:
  // mesh data
  std::vector<Vector<T, 3>> _Voxels;
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
  // neighbor cell direction
  /*
  neighbor rules
  1: self
  6: nearest neighbor <100>
  12: second nearest neighbor <110>
  8: third nearest neighbor <111>
  */
  // const std::vector<Vector<int, 3>> _NeighborDir = {
  //     {0, 0, 0},  {1, 0, 0},   {0, 1, 0},    {-1, 0, 0},
  //     {0, -1, 0}, {0, 0, 1},   {0, 0, -1},  // 0-6
  //     {1, 1, 0},  {-1, 1, 0},  {-1, -1, 0},  {1, -1, 0},
  //     {1, 0, 1},  {-1, 0, 1},  {-1, 0, -1},  {1, 0, -1},
  //     {0, 1, 1},  {0, -1, 1},  {0, -1, -1},  {0, 1, -1},  // 7-18
  //     {1, 1, 1},  {-1, 1, 1},  {-1, -1, 1},  {1, -1, 1},
  //     {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1}, {1, -1, -1}};  // 19-26

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

 public:
  // constructors
  // construct from stl file
  Geometry3D(const StlReader<T>& reader,
             std::uint8_t stlflag = std::uint8_t(1));
  // construct from basic shapes provided by freelb
  Geometry3D(int Nx, int Ny, int Nz, const AABB<T, 3>& AABBs,
             T voxelSize = T(1), const Vector<T, 3>& min = Vector<T, 3>{},
             std::uint8_t AABBflag = std::uint8_t(2),
             std::uint8_t voidflag = std::uint8_t(1));
  // construct from stl file and basic shapes provided by freelb
  // this constructor get AABBs from function argument
  // voxelsize is gained from stlreader
  Geometry3D(int Nx, int Ny, int Nz, const AABB<T, 3>& AABBs,
             const StlReader<T>& reader, std::uint8_t AABBflag,
             std::uint8_t stlflag, const Vector<T, 3>& min = Vector<T, 3>{});
  ~Geometry3D() = default;
  // read from octree
  void ReadOctree(Octree<T>* tree, std::uint8_t stlflag);
  // read from AABBs
  void ReadAABBs(const AABB<T, 3>& AABBs, std::uint8_t AABBflag);

  // setup boundary
  template <typename LatSet>
  void SetupBoundary(std::uint8_t AABBflag = std::uint8_t(2),
                     std::uint8_t boundaryflag = std::uint8_t(4));

  // set flag in a region defined by AABBs
  void setFlag(const AABB<T, 3>& AABBs, std::uint8_t flag);
  void setFlag(const AABB<T, 3>& AABBs, std::uint8_t fromflag,
               std::uint8_t flag);
  // set flag from stlreader in a region defined by AABBs
  void setFlag(const AABB<T, 3>& AABBs, const StlReader<T>& reader,
               std::uint8_t flag);
  void setFlag(const AABB<T, 3>& AABBs, const StlReader<T>& reader,
               std::uint8_t fromflag, std::uint8_t flag);

  // lambda function set flag
  // call: setFlag<LatSet>(AABBs, [](int id){func(id);});
  template <typename Func>
  void forEachVoxel(const AABB<T, 3>& AABBs, Func func); 
  template <typename Func>
  void forEachVoxel(const AABB<T, 3>& AABBs, std::uint8_t fromflag,
                      Func func);
  template <typename Func>
  void forEachVoxel(std::uint8_t fromflag, Func func);            
  // print
  void print() {
    std::cout << "[VoxelGeometry3D]: "
              << "\n"
              << "(Nx,Ny,Nz) = (" << _Nx << ", " << _Ny << ", " << _Nz << "); "
              << "min = (" << this->_min[0] << ", " << this->_min[1] << ", "
              << this->_min[2] << "); " << "max = (" << this->_max[0] << ", " << this->_max[1] << ", "
              << this->_max[2] << "); " <<std::endl;
  }

  // get
  int getNx() const { return _Nx; }
  int getNy() const { return _Ny; }
  int getNz() const { return _Nz; }
  // return N = Nx * Ny * Nz
  int getVoxelsNum() const { return N; }
  T getVoxelSize() const { return _voxelSize; }
  const Vector<int, 3>& getProjection() const { return Projection; }
  int getIndex(int i, int j, int k) const {
    return i + j * Projection[1] + k * Projection[2];
  }
  int getIndex(const Vector<int, 3>& loc) const { return loc[0] + loc[1] * Projection[1] + loc[2] * Projection[2]; }

  template <typename LatSet>
  int getNeighborId(int id, int dir) const {
    return id + LatSet::c[dir][0] + LatSet::c[dir][1] * Projection[1] +
           LatSet::c[dir][2] * Projection[2];
  }
  // {0, 1, -1, Nx, -Nx, Nx * Ny, -Nx * Ny, 
  //  1 + Nx, -1 - Nx, 1 + Nx*Ny, -1 - Nx*Ny, Nx + Nx*Ny, -Nx - Nx*Ny,
  // 1 - Nx, -1 + Nx, 1 - Nx*Ny, -1 + Nx*Ny, Nx - Nx*Ny, -Nx + Nx*Ny,
  // 1 + Nx + Nx*Ny, -1 - Nx - Nx*Ny, 1 + Nx - Nx*Ny, -1 - Nx + Nx*Ny,
  // 1 - Nx + Nx*Ny, -1 + Nx - Nx*Ny, -1 + Nx + Nx*Ny, 1 - Nx - Nx*Ny}
  // NOT SUITABLE FOR D3Q15!
  int getNbrId(int id, int dir) const { return id + Delta_Index[dir]; }
  std::array<int, 27> getDeltaIndex() { return Delta_Index; }
  
  const Vector<T, 3>& getMin() const { return this->_min; }
  const Vector<T, 3>& getMax() const { return this->_max; }

  std::vector<Vector<T, 3>>& getVoxels() { return _Voxels; }
  const std::vector<Vector<T, 3>>& getVoxels() const { return _Voxels; }
  Vector<T, 3>& getVoxel(int id) { return _Voxels[id]; }
  const Vector<T, 3>& getVoxel(int id) const { return _Voxels[id]; }
  FlagField& getGeoFlagField() { return GeometryFlag; }

  std::uint8_t getAABBflag() const { return _AABBflag; }
  std::uint8_t getVoidflag() const { return _voidflag; }
  template <typename LatSet>
  bool hasNeighborFlag(int id, std::uint8_t flag) const {
    for (int dir = 1; dir < LatSet::q; dir++) {
      if (GeometryFlag.get(getNeighborId(id, dir)) == flag)
        return true;
    }
    return false;
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
