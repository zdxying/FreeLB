// voxel_geometry2d.h
#pragma once

#include "geometry/basic_geometry.h"
#include "legacy_lattice/legacy_voxel.h"

// legacy, low efficiency voxelgeometry
template <typename T>
class VoxelGeometry2D : public AABB<T, 2> {
 private:
  // mesh data
  std::vector<Voxel<T, 2>> _Voxels;
  // mesh param
  int _Nx;
  int _Ny;
  // it is recommended to set voxelSize = 1.
  T _voxelSize;
  // facilitate getting index in a cubic array
  // _idx (1, Nx)
  //      (i, j)
  // id = i + j * Nx
  Vector<int, 2> _idx;
  // index to find voxel in std::vector<Voxel<T, 2>> _Voxels
  // std::vector<int> _VoxelIdx;
  // a rectangular array(_Nx * _Ny), index of voxel in _Voxels,
  std::vector<int> _GlobalIdx;

  int InnerFlag;
  int BoundaryFlag;
  int InterfaceFlag;

 public:
  // constructors
  // construct from basic shapes provided by freelb
  VoxelGeometry2D(int Nx, int Ny, const AABB<T, 2>& AABBs, T voxelSize = T(1),
                  const Vector<T, 2>& min = Vector<T, 2>{});
  // read from octree
  // void ReadOctree(Octree<T>* tree);
  // read from AABBs
  void ReadAABBs(const AABB<T, 2>& AABBs);

  // normal setup
  template <typename LatSet>
  void Setup(int innerflag = 0, int boundaryflag = 1);
  // set neighbor of each voxel
  // for obj constructed from stl file and basic shapes provided by freelb
  // template <typename LatSet>
  // void Setup(int fromflag, int toBdflag);
  // write voxels to file
  template <bool boundaryinfo = false>
  void WriteStruPoints();
  // print
  void print() {
    std::cout << "[VoxelGeometry2D]: "
              << "\n"
              << "(Nx,Ny) = (" << _Nx << ", " << _Ny << "); "
              << "min = (" << this->_min[0] << ", " << this->_min[1] << ")"
              << std::endl;
  }
  // set
  void setFlag(int id, int flag) { _Voxels[id].setFlag(flag); }
  // set flag in a region defined by AABBs
  void setFlag(const AABB<T, 2>& AABBs, int flag);
  void setFlag(const AABB<T, 2>& AABBs, int fromflag, int toflag);

  // -------------set filed flag outside of this class
  template <typename U = int>
  void setFieldFlag(std::vector<U>& field, const U& flag) {
    for (const Voxel<T, 2>& vox : _Voxels) field[vox.getId()] = flag;
  }
  template <typename U = int>
  void setFieldFlag(std::vector<U>& field, const U& flag, int voxflag) {
    for (const Voxel<T, 2>& vox : _Voxels) {
      if (vox.getFlag() == voxflag) field[vox.getId()] = flag;
    }
  }
  // set filed flag in a region defined by AABBs
  template <typename U = int>
  void setFieldFlag(const AABB<T, 2>& AABBs, std::vector<U>& field,
                    const U& flag);
  template <typename U = int>
  void setFieldFlag(const AABB<T, 2>& AABBs, std::vector<U>& field,
                    const U& fromflag, const U& toflag);
  template <typename U = int>
  void setFieldFlag(const AABB<T, 2>& AABBs, std::vector<U>& field,
                    const U& flag, int voxflag);

  // -------------end set filed flag outside of this class

  // get
  int getNx() const { return _Nx; }
  int getNy() const { return _Ny; }
  T getVoxelSize() const { return _voxelSize; }
  const Vector<int, 2>& getIdx() const { return _idx; }

  std::vector<int>& getGlobalIdx() { return _GlobalIdx; }
  const std::vector<int>& getGlobalIdx() const { return _GlobalIdx; }

  const Vector<T, 2>& getMin() const { return this->_min; }
  const Vector<T, 2>& getMax() const { return this->_max; }

  std::vector<Voxel<T, 2>>& getVoxels() { return _Voxels; }
  const std::vector<Voxel<T, 2>>& getVoxels() const { return _Voxels; }

  Voxel<T, 2>& getVoxel(int id) { return _Voxels[id]; }
  const Voxel<T, 2>& getVoxel(int id) const { return _Voxels[id]; }
};