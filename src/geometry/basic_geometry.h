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

// basic_geometry.h

#pragma once

#include <limits>

#include "data_struct/field_struct.h"
#include "data_struct/interpolation.h"
#include "parallel/communicator.h"

// axis-aligned bounding box (AABB) class
// all the voxels in the computational domain are inside the AABB
template <typename T, unsigned int D>
class AABB {
 protected:
  // max and min point of AABB
  Vector<T, D> _min;
  Vector<T, D> _max;
  // extension of AABB
  Vector<T, D> _extension;
  // center of AABB
  Vector<T, D> _center;

 public:
  // get
  const Vector<T, D>& getMin() const { return _min; }
  const Vector<T, D>& getMax() const { return _max; }
  const Vector<T, D>& getExtension() const { return _extension; }
  const Vector<T, D>& getCenter() const { return _center; }
  //
  AABB(const Vector<T, D>& min, const Vector<T, D>& max)
      : _min(min), _max(max), _extension(max - min), _center((_min + _max) / T(2)) {}
  AABB(const Vector<T, D>& centre, T extension0, T extension1, T extension2 = T(0))
      : _center(centre) {
    if constexpr (D == 2) {
      _extension = Vector<T, 2>{extension0, extension1};
      _min = _center - _extension / T(2);
      _max = _center + _extension / T(2);
    } else {
      _extension = Vector<T, 3>{extension0, extension1, extension2};
      _min = _center - _extension / T(2);
      _max = _center + _extension / T(2);
    }
  }
  AABB(const AABB<T, D>& aabb)
      : _min(aabb.getMin()), _max(aabb.getMax()), _extension(aabb.getExtension()),
        _center(aabb.getCenter()) {}
  AABB() : _min(T(0)), _max(T(0)), _extension(T(0)), _center(T(0)) {}

  // check if a point Vector<S, D> is inside the AABB, with additional template parameter
  // S
  template <typename S>
  bool IsInside(const Vector<S, D>& pt) const;
  // check if a point Vector<T, D> is inside the AABB
  virtual bool isInside(const Vector<T, D>& pt) const;
  // get extended AABB
  AABB<T, D> getExtended(T extension = T(1)) const;
  AABB<T, D> getExtended(const Vector<T, D>& extension) const;
  // get Resized AABB
  AABB<T, D> getResized(T factor) const;
  // get overlapped local cell index
  // designed to call on AABB<int, D>
  // AABBa.getOverlappedCellIdx(AABBb, cellIdx) will add a's local cell index to cellIdx
  // legacy function, use getCellIdx instead
  void getOverlappedCellIdx(const AABB<int, D>& aabb, std::vector<int>& cellIdx) const;
  // get local index of a point in the block without checking if is inside
  // call on AABB<int, D>
  std::size_t getLocIdx(const Vector<int, D>& pt) const;
  // get overlapped global cell index
  // designed to call on AABB<int, D>,
  // 2 AABBs is supposed to be described in the same coordinate system
  // legacy function, use getCellIdx instead
  void getOverlappedCellIdx(const AABB<int, D>& aabb, std::vector<int>& cellIdx,
                            const Vector<int, D>& shift) const;
  // get local index of a point in the block without checking if is inside
  // call on AABB<int, D>
  std::size_t getLocIdx(const Vector<int, D>& pt, const Vector<int, D>& shift) const;

  // called on AABB<int, D>
  void divide(int Nx, int Ny, std::vector<AABB<int, 2>>& subAABBs) const;
  void divide(int Nx, int Ny, int Nz, std::vector<AABB<int, 3>>& subAABBs) const ;
};

// get intersection of 2 AABBs without checking if is intersected
template <typename T, unsigned int D>
static AABB<T, D> getIntersection(const AABB<T, D>& aabb0, const AABB<T, D>& aabb1) {
  if constexpr (D == 2) {
    Vector<T, 2> min{std::max(aabb0.getMin()[0], aabb1.getMin()[0]),
                     std::max(aabb0.getMin()[1], aabb1.getMin()[1])};
    Vector<T, 2> max{std::min(aabb0.getMax()[0], aabb1.getMax()[0]),
                     std::min(aabb0.getMax()[1], aabb1.getMax()[1])};
    return AABB<T, 2>{min, max};
  } else if constexpr (D == 3) {
    Vector<T, 3> min{std::max(aabb0.getMin()[0], aabb1.getMin()[0]),
                     std::max(aabb0.getMin()[1], aabb1.getMin()[1]),
                     std::max(aabb0.getMin()[2], aabb1.getMin()[2])};
    Vector<T, 3> max{std::min(aabb0.getMax()[0], aabb1.getMax()[0]),
                     std::min(aabb0.getMax()[1], aabb1.getMax()[1]),
                     std::min(aabb0.getMax()[2], aabb1.getMax()[2])};
    return AABB<T, 3>{min, max};
  }
}

// contiguous blocks(share a common edge or face) are NOT overlapped
template <typename T, typename U, unsigned int D>
bool isOverlapped(const AABB<T, D>& aabb0, const AABB<U, D>& aabb1) {
  static_assert(D == 2 || D == 3, "Error: Dimension is not supported!");
  const T eps = std::numeric_limits<T>::epsilon();
  if constexpr (D == 2) {
    return ((aabb0.getMin()[0] + eps) < aabb1.getMax()[0] &&
            (aabb0.getMin()[1] + eps) < aabb1.getMax()[1] &&
            (aabb1.getMin()[0] + eps) < aabb0.getMax()[0] &&
            (aabb1.getMin()[1] + eps) < aabb0.getMax()[1]);
  } else if constexpr (D == 3) {
    return ((aabb0.getMin()[0] + eps) < aabb1.getMax()[0] &&
            (aabb0.getMin()[1] + eps) < aabb1.getMax()[1] &&
            (aabb0.getMin()[2] + eps) < aabb1.getMax()[2] &&
            (aabb1.getMin()[0] + eps) < aabb0.getMax()[0] &&
            (aabb1.getMin()[1] + eps) < aabb0.getMax()[1] &&
            (aabb1.getMin()[2] + eps) < aabb0.getMax()[2]);
  }
}


/*
geometry of a block structure
public AABB<T, D>: position of the block
  _min{x, y, z} _max{_min + voxsize*{Nx-1, Ny-1, Nz-1}}
AABB<int, D> IndexBlock: global index block
  _min{i, j, k} _max{i + Nx - 1, j + Ny - 1, k + Nz - 1}
Vector<int, D> Mesh: Mesh parameter of the block
  {Nx, Ny, Nz}
Projection: {1, Nx} or {1, Nx, Nx*Ny}
int N: N = Nx * Ny * Nz
T VoxelSize: voxel size
Vector<T, D> MinCenter: min cell center
  _min + Vector<T, D>{T(0.5) * voxsize}
int BlockId: block id
*/
template <typename T, unsigned int D>
class BasicBlock : public AABB<T, D> {
 protected:
  // _min{minx, miny, minz}
  // _max{minx + Nx - 1, miny + Ny - 1, minz + Nz - 1}
  // (IndexBlock.getExtension()+1) * VoxelSize = _extension
  AABB<int, D> IndexBlock;
  // Mesh{Nx, Ny, Nz}
  Vector<int, D> Mesh;
  // N = Nx * Ny * Nz
  std::size_t N;
  // Projection: {1, Nx} or {1, Nx, Nx*Ny}
  Vector<int, D> Projection;
  // voxel size
  T VoxelSize;
  // min cell center
  Vector<T, D> MinCenter;
  // for block structure
  int BlockId;
  // ref level
  std::uint8_t _level;

 public:
  BasicBlock() = default;
  BasicBlock(T voxsize, const AABB<T, D>& aabb, const AABB<int, D>& idxblock,
             int blockid = -1)
      : VoxelSize(voxsize), AABB<T, D>(aabb),
        MinCenter(aabb.getMin() + Vector<T, D>{T(0.5) * voxsize}), IndexBlock(idxblock),
        Mesh(idxblock.getExtension() + Vector<int, D>{1}), BlockId(blockid),
        _level(std::uint8_t(0)) {
    if constexpr (D == 2) {
      N = Mesh[0] * Mesh[1];
      Projection = Vector<int, 2>{1, Mesh[0]};
    } else if constexpr (D == 3) {
      N = Mesh[0] * Mesh[1] * Mesh[2];
      Projection = Vector<int, 3>{1, Mesh[0], Mesh[0] * Mesh[1]};
    }
  }
  // for refined block, voxelsize should satisfy: newvoxsize = parent_voxsize / 2;
  // indexblock is the GLOBAL index block
  BasicBlock(std::uint8_t level, T newvoxsize, int blockid, const AABB<T, D>& aabb,
             const AABB<int, D>& idxblock, const Vector<int, D>& mesh)
      : VoxelSize(newvoxsize), _level(level), AABB<T, D>(aabb), BlockId(blockid),
        MinCenter(aabb.getMin() + Vector<T, D>{T(0.5) * newvoxsize}),
        IndexBlock(idxblock), Mesh(mesh) {
    if constexpr (D == 2) {
      N = Mesh[0] * Mesh[1];
      Projection = Vector<int, 2>{1, Mesh[0]};
    } else if constexpr (D == 3) {
      N = Mesh[0] * Mesh[1] * Mesh[2];
      Projection = Vector<int, 3>{1, Mesh[0], Mesh[0] * Mesh[1]};
    }
  }
  // get
  const AABB<int, D>& getIdxBlock() const { return IndexBlock; }
  const AABB<T, D>& getAABB() const { return *this; }
  const Vector<int, D>& getProjection() const { return Projection; }
  const Vector<int, D>& getMesh() const { return Mesh; }
  Vector<int, D>* getMeshPtr() { return &Mesh; }
  const Vector<T, D>& getMinCenter() const { return MinCenter; }

  int getNx() const { return Mesh[0]; }
  int getNy() const { return Mesh[1]; }
  int getNz() const {
    if constexpr (D == 3)
      return Mesh[2];
    else if constexpr (D == 2)
      return 1;
  }
  // return total number of cells in the block
  std::size_t getN() const { return N; }

  std::uint8_t getLevel() const { return _level; }
  T getVoxelSize() const { return VoxelSize; }
  int getBlockId() const { return BlockId; }
  void setBlockId(int id) { BlockId = id; }

  // getExtBlock(Vector<int, D>{ext})
  BasicBlock<T, D> getExtBlock(int ext = 1) const;
  // get the extended block, similar to getExtended in AABB
  BasicBlock<T, D> getExtBlock(const Vector<int, D>& extension) const;
  // suggested to call on base block to get refined base block, deltalevel>0
  BasicBlock<T, D> getRefinedBlock(std::uint8_t deltalevel = std::uint8_t(1)) const;
  void refine(std::uint8_t deltalevel = std::uint8_t(1));
  // suggested to call on base block to get coarsened base block, _deltalevel>0
  BasicBlock<T, D> getCoasenedBlock(std::uint8_t _deltalevel = std::uint8_t(1)) const;
  void coarsen(std::uint8_t _deltalevel = std::uint8_t(1));

  // get cell center
  inline Vector<T, D> getVoxel(const Vector<int, D>& locidx) const {
    return MinCenter + (locidx * VoxelSize);
  }
  // get LOCAL index of a point in the block, locidx is the Global position
  std::size_t getIndex_t(const Vector<T, D>& locidx) const;
  // get LOCAL index of a point in the block, locidx is the LOCAL IDX position
  std::size_t getIndex(const Vector<int, D>& locidx) const;
  // get location index of a point in the block, id is the LOCAL index
  void getLoc(std::size_t id, Vector<int, D>& locidx) const;
  Vector<int, D> getLoc(std::size_t id) const;
  // get real position in the block
  void getLoc_t(std::size_t id, Vector<T, D>& loc) const;
  Vector<T, D> getLoc_t(std::size_t id) const;
  Vector<T, D> getLoc_t(Vector<int, D>& locidx) const;

  // get min and max LOCAL index(within the current BasicBlock) of a given AABB
  void getLocIdxRange(const AABB<T, D>& AABBs, Vector<int, D>& idx_min,
                      Vector<int, D>& idx_max) const;

  // lambda functions take LOCAL index: func(std::size_t idx)
  template <typename Func>
  void forEach(Func func);
  
  template <typename Func>
  void forEach(const AABB<T, D>& AABBs, Func func);

  // lambda functions take LOCAL index, for cells with specific flag
  template <typename ArrayType, typename Func>
  void forEach(const AABB<T, D>& AABBs, const ArrayType& flag,
               std::uint8_t fromflag, Func func);

  template <typename ArrayType, typename Func>
  void forEach(const ArrayType& flag, std::uint8_t fromflag, Func func);

  // get certain cells' indices(usually in an overlapped AABB) in the block
  // similar to void getOverlappedCellIdx(const AABB<int, D>& aabb, std::vector<int>&
  // cellIdx) const; but this function is designed to recieve an AABB<T, D> always call on
  // extended BasicBlock base could be the base blockaabb or the extended blockaabb
  // cellIdx will be cleared before adding new indices
  void getCellIdx(const AABB<T, D>& base, const AABB<T, D>& AABBs,
                  std::vector<std::size_t>& cellIdx) const;
};

// communication structure for block geometry
template <typename T, unsigned int D>
struct BlockComm {
  std::vector<std::size_t> SendCells;  // index in sendblock
  std::vector<std::size_t> RecvCells;  // index in recvblock
  Block<T, D>* SendBlock;

  BlockComm(Block<T, D>* block) : SendBlock(block) {}
  int getSendId() const { return SendBlock->getBlockId(); }
};

// communication structure for blocks of different (refine) levels
// heterogeneous communication/ interpolation
template <typename T, unsigned int D>
struct InterpBlockComm {
  // index in sendblock
  std::vector<InterpSource<D>> SendCells;
  // interp weight
  // std::vector<InterpWeight<T, D>> InterpWeights;
  // index in recvblock
  std::vector<std::size_t> RecvCells;
  // receive data from sendblock
  Block<T, D>* SendBlock;
  // predefined interp weight
  static constexpr std::array<InterpWeight<T, 2>, 4> InterpWeight2D{{
    {T(0.0625), T(0.1875), T(0.1875), T(0.5625)}, 
    {T(0.1875), T(0.0625), T(0.5625), T(0.1875)}, 
    {T(0.1875), T(0.5625), T(0.0625), T(0.1875)}, 
    {T(0.5625), T(0.1875), T(0.1875), T(0.0625)}}};
  static constexpr std::array<InterpWeight<T, 3>, 8> InterpWeight3D{
    {T{0.015625}, T{0.046875}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.140625}, T{0.421875}},
    {T{0.046875}, T{0.015625}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.421875}, T{0.140625}},
    {T{0.046875}, T{0.140625}, T{0.015625}, T{0.046875}, T{0.140625}, T{0.421875}, T{0.046875}, T{0.140625}},
    {T{0.140625}, T{0.046875}, T{0.046875}, T{0.015625}, T{0.421875}, T{0.140625}, T{0.140625}, T{0.046875}},
    {T{0.046875}, T{0.140625}, T{0.140625}, T{0.421875}, T{0.015625}, T{0.046875}, T{0.046875}, T{0.140625}},
    {T{0.140625}, T{0.046875}, T{0.421875}, T{0.140625}, T{0.046875}, T{0.015625}, T{0.140625}, T{0.046875}},
    {T{0.140625}, T{0.421875}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.015625}, T{0.046875}},
    {T{0.421875}, T{0.140625}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.046875}, T{0.015625}}
  };

  static constexpr auto getIntpWeight() {
    if constexpr (D == 2) {
      return InterpWeight2D;
    } else {
      return InterpWeight3D;
    }
  }

  InterpBlockComm(Block<T, D>* block) : SendBlock(block) {}

  int getSendId() const { return SendBlock->getBlockId(); }
  // return (D == 2) ? 4 : 8;
  static constexpr int getSourceNum() { return (D == 2) ? 4 : 8; }
  // return (D == 2) ? T(0.25) : T(0.125);
  static constexpr T getUniformWeight() { return (D == 2) ? T(0.25) : T(0.125); }
};


template <typename T>
static void DivideBlock2D(BasicBlock<T, 2>& block, int blocknum,
                          std::vector<AABB<int, 2>>& blocksvec) {
  DivideBlock2D<T>(block.getNx(), block.getNy(), blocknum, block.getIdxBlock(),
                   blocksvec);
}

template <typename T>
static void DivideBlock2D(BasicBlock<T, 2>& block, int blocknum,
                          std::vector<BasicBlock<T, 2>>& basicblocksvec) {
  std::vector<AABB<int, 2>> blocksvec;
  const AABB<int, 2>& IdxBlock = block.getIdxBlock();
  // divide block
  DivideBlock2D<T>(block.getNx(), block.getNy(), blocknum, IdxBlock, blocksvec);
  // construct basicblocks from aabb
  std::uint8_t level = block.getLevel();
  T voxsize = block.getVoxelSize();
  const Vector<T, 2>& Min = block.getMin();
  int ratio = static_cast<int>(std::pow(2, static_cast<int>(level)));
  int blockid = 0;
  for (const AABB<int, 2>& idxaabb : blocksvec) {
    // relative index location
    Vector<int, 2> idxmin = idxaabb.getMin() - IdxBlock.getMin();
    Vector<int, 2> idxmax = idxaabb.getMax() - IdxBlock.getMin();
    Vector<T, 2> min = Min + voxsize * ratio * idxmin;
    // remember to add Vector<int, 2>{1}
    Vector<T, 2> max = Min + voxsize * ratio * (idxmax + Vector<int, 2>{1});
    AABB<T, 2> aabb{min, max};
    // mesh
    Vector<int, 2> Mesh = (idxaabb.getExtension() + Vector<int, 2>{1}) * ratio;
    basicblocksvec.emplace_back(level, voxsize, blockid, aabb, idxaabb, Mesh);
    ++blockid;
  }
}

template <typename T>
static void DivideBlock2D(int NX, int NY, int blocknum, const AABB<int, 2>& idxAABBs,
                          std::vector<AABB<int, 2>>& blocksvec) {
  // int NX = _Nx - 2;
  // int NY = _Ny - 2;
  // inner AABB of the domain
  // AABB<int, 2> idxAABBs(Vector<int, 2>(1, 1), Vector<int, 2>(NX, NY));
  if (blocknum == 1) {
    blocksvec.push_back(idxAABBs);
    return;
  }
  int Xblocks = 0;
  int Yblocks = 0;
  T bestRatio = T(NX) / T(NY);
  T difRatio = std::fabs(bestRatio - 1) + 1;  // >= 1
  for (int i = 1; i <= blocknum; i++) {
    int j = blocknum / i;
    if (i * j <= blocknum) {
      if (std::fabs(bestRatio - T(i) / T(j)) <= difRatio) {
        difRatio = std::fabs(bestRatio - T(i) / T(j));
        Xblocks = i;
        Yblocks = j;
      }
    }
  }

  T ratio = T(Xblocks) / T(Yblocks);
  int rest = blocknum - Xblocks * Yblocks;

  // create aabbs
  if (rest == 0) {
    idxAABBs.divide(Xblocks, Yblocks, blocksvec);
    return;
  }

  if (ratio < bestRatio && (Yblocks - rest) >= 0) {
    // divide along y direction
    int rest_blocknum = Xblocks * (Yblocks - rest);
    T bestVolume = (T)NX * NY * rest_blocknum / (T)blocknum;
    int seg_Y = (int)(bestVolume / (T)NX);
    Vector<int, 2> min0 = idxAABBs.getMin();
    Vector<int, 2> max0 = Vector<int, 2>{idxAABBs.getMax()[0], min0[1] + seg_Y - 1};
    Vector<int, 2> min1 = Vector<int, 2>{min0[0], min0[1] + seg_Y};
    Vector<int, 2> max1 = idxAABBs.getMax();
    AABB<int, 2> Child_0(min0, max0);
    AABB<int, 2> Child_1(min1, max1);
    Child_0.divide(Xblocks, Yblocks - rest, blocksvec);
    Child_1.divide(Xblocks + 1, rest, blocksvec);
  } else {
    // divide along x direction
    int rest_blocknum = Yblocks * (Xblocks - rest);
    T bestVolume = (T)NX * NY * rest_blocknum / (T)blocknum;
    int seg_X = (int)(bestVolume / (T)NY + 0.9999);
    Vector<int, 2> min0 = idxAABBs.getMin();
    Vector<int, 2> max0 = Vector<int, 2>{min0[0] + seg_X - 1, idxAABBs.getMax()[1]};
    Vector<int, 2> min1 = Vector<int, 2>{min0[0] + seg_X, min0[1]};
    Vector<int, 2> max1 = idxAABBs.getMax();
    AABB<int, 2> Child_0(min0, max0);
    AABB<int, 2> Child_1(min1, max1);
    Child_0.divide(Xblocks - rest, Yblocks, blocksvec);
    Child_1.divide(rest, Yblocks + 1, blocksvec);
  }
}


template <typename T>
static void DivideBlock3D(BasicBlock<T, 3>& block, int blocknum,
                          std::vector<AABB<int, 3>>& blocksvec) {
  DivideBlock3D<T>(block.getNx(), block.getNy(), block.getNz(), blocknum, block.getIdxBlock(),
                   blocksvec);
}

template <typename T>
static void DivideBlock3D(BasicBlock<T, 3>& block, int blocknum,
                          std::vector<BasicBlock<T, 3>>& basicblocksvec) {
  std::vector<AABB<int, 3>> blocksvec;
  const AABB<int, 3>& IdxBlock = block.getIdxBlock();
  // divide block
  DivideBlock3D<T>(block.getNx(), block.getNy(), block.getNz(), blocknum, IdxBlock, blocksvec);
  // construct basicblocks from aabb
  std::uint8_t level = block.getLevel();
  T voxsize = block.getVoxelSize();
  const Vector<T, 3>& Min = block.getMin();
  int ratio = static_cast<int>(std::pow(2, static_cast<int>(level)));
  int blockid = 0;
  for (const AABB<int, 3>& idxaabb : blocksvec) {
    // relative index location
    Vector<int, 3> idxmin = idxaabb.getMin() - IdxBlock.getMin();
    Vector<int, 3> idxmax = idxaabb.getMax() - IdxBlock.getMin();
    Vector<T, 3> min = Min + voxsize * ratio * idxmin;
    // remember to add Vector<int, 2>{1}
    Vector<T, 3> max = Min + voxsize * ratio * (idxmax + Vector<int, 3>{1});
    AABB<T, 3> aabb{min, max};
    // mesh
    Vector<int, 3> Mesh = (idxaabb.getExtension() + Vector<int, 3>{1}) * ratio;
    basicblocksvec.emplace_back(level, voxsize, blockid, aabb, idxaabb, Mesh);
    ++blockid;
  }
}

template <typename T>
static void DivideBlock3D(int NX, int NY, int Nz, int blocknum, const AABB<int, 3>& idxAABBs,
                          std::vector<AABB<int, 3>>& blocksvec) {

}