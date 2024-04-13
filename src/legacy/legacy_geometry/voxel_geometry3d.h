

#pragma once

#include <fstream>
#include <sstream>

#include "geometry/basic_geometry.h"
#include "data_struct/field_struct.h"
#include "io/stlreader.h"
#include "legacy_lattice/legacy_voxel.h"

template <typename T>
class VoxelGeometry3D : public AABB<T, 3> {
 private:
  // mesh data
  std::vector<Voxel<T, 3>> _Voxels;
  // mesh param
  int _Nx;
  int _Ny;
  int _Nz;
  // it is recommended to set voxelSize = 1.
  T _voxelSize;
  // facilitate getting index in a cubic array
  // _idx (1, Nx, Nx * Ny)
  //      (i, j, k)
  // id = i + j * Nx + k * Nx * Ny
  Vector<int, 3> _idx;
  // index to find voxel in std::vector<Voxel<T, 3>> _Voxels
  // std::vector<int> _VoxelIdx;
  // a cubic array(_Nx * _Ny * _Nz), index of voxel in _Voxels,
  std::vector<int> _GlobalIdx;
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

 public:
  // constructors
  // construct from stl file
  VoxelGeometry3D(const StlReader<T>& reader);
  // construct from basic shapes provided by freelb
  VoxelGeometry3D(int Nx, int Ny, int Nz, const AABB<T, 3>& AABBs,
                  T voxelSize = T(1),
                  const Vector<T, 3>& min = Vector<T, 3>{T(0), T(0), T(0)});
  // construct from stl file and basic shapes provided by freelb
  // this constructor get AABBs from function argument
  // voxelsize is gained from stlreader
  VoxelGeometry3D(int Nx, int Ny, int Nz, const AABB<T, 3>& AABBs,
                  const StlReader<T>& reader, int flag,
                  const Vector<T, 3>& min = Vector<T, 3>{T(0), T(0), T(0)});
  // read from octree
  void ReadOctree(Octree<T>* tree);
  // read from AABBs
  void ReadAABBs(const AABB<T, 3>& AABBs);
  // set neighbor of each voxel
  // this will be called after constructing VoxelGeometry3D
  // if constructed from stl file and basic shapes provided by freelb
  // DO NOT call this function
  template <typename LatSet>
  void Setup();
  // set neighbor of each voxel
  // for obj constructed from stl file and basic shapes provided by freelb
  template <typename LatSet>
  void Setup(int AABBflag, int AABBBdflag, int stlflag, int stlBdflag);
  // write voxels to file
  template <bool boundaryinfo = false>
  void WriteStruPoints();
  void WriteStruGrid();
  // print
  void print() {
    std::cout << "[VoxelGeometry3D]: "
              << "\n"
              << "(Nx,Ny,Nz) = (" << _Nx << ", " << _Ny << ", " << _Nz << "); "
              << "min = (" << this->_min[0] << ", " << this->_min[1] << ", "
              << this->_min[2] << ")" << std::endl;
  }
  // set
  void setFlag(int id, int flag) { _Voxels[id].setFlag(flag); }
  // set flag in a region defined by AABBs
  void setFlag(const AABB<T, 3>& AABBs, int flag);
  void setFlag(const AABB<T, 3>& AABBs, int fromflag, int flag);
  // set flag from stlreader in a region defined by AABBs
  void setFlag(const AABB<T, 3>& AABBs, const StlReader<T>& reader, int flag);
  void setFlag(const AABB<T, 3>& AABBs, const StlReader<T>& reader,
               int fromflag, int flag);

  // get
  int getNx() const { return _Nx; }
  int getNy() const { return _Ny; }
  int getNz() const { return _Nz; }
  T getVoxelSize() const { return _voxelSize; }
  const Vector<int, 3>& getIdx() const { return _idx; }

  std::vector<int>& getGlobalIdx() { return _GlobalIdx; }
  const std::vector<int>& getGlobalIdx() const { return _GlobalIdx; }

  const Vector<T, 3>& getMin() const { return this->_min; }
  const Vector<T, 3>& getMax() const { return this->_max; }

  // return std::vector<Voxel<T, 3>>
  std::vector<Voxel<T, 3>>& getVoxels() { return _Voxels; }
  const std::vector<Voxel<T, 3>>& getVoxels() const { return _Voxels; }
  // return Voxel<T, 3> at _Voxels[id]
  Voxel<T, 3>& getVoxel(int id) { return _Voxels[id]; }
  const Voxel<T, 3>& getVoxel(int id) const { return _Voxels[id]; }
};
