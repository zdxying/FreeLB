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

// used in isOverlapped() 
#include <limits>

#include <vector>
// uint8_t size_t
#include <cstdint>
// pow()
#include <cmath>
// std::max(), std::min()
#include <algorithm>

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
  void divide(int Nx, int Ny, int Nz, std::vector<AABB<int, 3>>& subAABBs) const;
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
  void forEach(const AABB<T, D>& AABBs, const ArrayType& flag, std::uint8_t fromflag,
               Func func);

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
  // bestratio in the range of (0, +inf)
  T bestRatio = T(NX) / T(NY);
  // init difratio, >= 1, make sure in the first loop, dist_to_bestRatio <= difRatio
  T difRatio = std::fabs(bestRatio - 1) + 1;
  for (int i = 1; i <= blocknum; i++) {
    int j = blocknum / i;
    // find the Xblocks and Yblocks with the closest ratio to bestRatio
    T dist_to_bestRatio = std::fabs(bestRatio - T(i) / T(j));
    if (dist_to_bestRatio <= difRatio) {
      difRatio = dist_to_bestRatio;
      Xblocks = i;
      Yblocks = j;
    }
  }

  // exchange i and j and find the best ratio?
  // for (int j = 1; j <= blocknum; j++) {
  //   int i = blocknum / j;
  //   T dist_to_bestRatio = std::fabs(bestRatio - T(i) / T(j));
  //   if (dist_to_bestRatio <= difRatio) {
  //     difRatio = dist_to_bestRatio;
  //     Xblocks = i;
  //     Yblocks = j;
  //   }
  // }
  // end of exchange i and j and find the best ratio

  T ratio = T(Xblocks) / T(Yblocks);
  int rest = blocknum - Xblocks * Yblocks;
  // xy + rest = blocknum
  // = x(y - rest) + (x + 1)rest
  // = y(x - rest) + (y + 1)rest

  // create aabbs
  if (rest == 0) {
    idxAABBs.divide(Xblocks, Yblocks, blocksvec);
    return;
  }

  // (Yblocks - rest) may be negative, but (Xblocks - rest) is always positive
  if (ratio < bestRatio && (Yblocks - rest) >= 0) {
    // divide along y direction
    // x(y - rest) 
    int rest_blocknum = Xblocks * (Yblocks - rest);
    T bestVolume = (T)NX * NY * rest_blocknum / (T)blocknum;
    // int (x(y - rest)/blocknum) * NY
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
    // y(x - rest), note that (Xblocks - rest) >= 1
    int rest_blocknum = Yblocks * (Xblocks - rest);
    T bestVolume = (T)NX * NY * rest_blocknum / (T)blocknum;
    // int (y(x - rest)/blocknum + 0.9999) * NX
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
  DivideBlock3D<T>(block.getNx(), block.getNy(), block.getNz(), blocknum,
                   block.getIdxBlock(), blocksvec);
}

template <typename T>
static void DivideBlock3D(BasicBlock<T, 3>& block, int blocknum,
                          std::vector<BasicBlock<T, 3>>& basicblocksvec) {
  std::vector<AABB<int, 3>> blocksvec;
  const AABB<int, 3>& IdxBlock = block.getIdxBlock();
  // divide block
  DivideBlock3D<T>(block.getNx(), block.getNy(), block.getNz(), blocknum, IdxBlock,
                   blocksvec);
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
    // remember to add Vector<int, 3>{1}
    Vector<T, 3> max = Min + voxsize * ratio * (idxmax + Vector<int, 3>{1});
    AABB<T, 3> aabb{min, max};
    // mesh
    Vector<int, 3> Mesh = (idxaabb.getExtension() + Vector<int, 3>{1}) * ratio;
    basicblocksvec.emplace_back(level, voxsize, blockid, aabb, idxaabb, Mesh);
    ++blockid;
  }
}

template <typename T>
static void DivideBlock3D(int NX, int NY, int NZ, int blocknum,
                          const AABB<int, 3>& idxAABBs,
                          std::vector<AABB<int, 3>>& blocksvec) {
  if (blocknum == 1) {
    blocksvec.push_back(idxAABBs);
    return;
  }
  Vector<int, 3> AABBminIdx = idxAABBs.getMin();
  Vector<int, 3> AABBmaxIdx = idxAABBs.getMax();

  int iXX = 1;
  int iYY = 1;
  int iZZ = blocknum;
  int nX = NX / iXX;
  int bestIx = iXX;
  int nY = NY / iYY;
  int bestIy = iYY;
  int nZ = NZ / iZZ;
  int bestIz = iZZ;
  T bestRatio =
    ((T)(NX / iXX) / (T)(NY / iYY) - 1) * ((T)(NX / iXX) / (T)(NY / iYY) - 1) +
    ((T)(NY / iYY) / (T)(NZ / iZZ) - 1) * ((T)(NY / iYY) / (T)(NZ / iZZ) - 1) +
    ((T)(NZ / iZZ) / (T)(NX / iXX) - 1) * ((T)(NZ / iZZ) / (T)(NX / iXX) - 1);

  for (int iX = 1; iX <= blocknum; iX++) {
    for (int iY = 1; iY * iX <= blocknum; iY++) {
      for (int iZ = blocknum / (iX * iY); iZ * iY * iX <= blocknum; iZ++) {
        if ((iX + 1) * iY * iZ > blocknum && iX * (iY + 1) * iZ > blocknum) {
          T ratio =
            ((T)(NX / iX) / (T)(NY / iY) - 1) * ((T)(NX / iX) / (T)(NY / iY) - 1) +
            ((T)(NY / iY) / (T)(NZ / iZ) - 1) * ((T)(NY / iY) / (T)(NZ / iZ) - 1) +
            ((T)(NZ / iZ) / (T)(NX / iX) - 1) * ((T)(NZ / iZ) / (T)(NX / iX) - 1);
          if (ratio < bestRatio) {
            bestRatio = ratio;
            bestIx = iX;
            bestIy = iY;
            bestIz = iZ;
            nX = NX / iX;
            nY = NY / iY;
            nZ = NZ / iZ;
          }
        }
      }
    }
  }

  int rest = blocknum - bestIx * bestIy * bestIz;

  if (rest == 0) {
    idxAABBs.divide(bestIx, bestIy, bestIz, blocksvec);
    return;
  } else {

    // add in z than in y direction
    // 1
    if (nZ > nY && nZ > nX) {
      int restY = rest % bestIy;
      // split in two cuboid
      if (restY == 0) {
        int restX = rest / bestIy;

        std::vector<AABB<int, 2>> GeoHelperVec;
        AABB<int, 2> AABBXZ(Vector<int, 2>{AABBminIdx[0], AABBminIdx[2]},
                            Vector<int, 2>{AABBmaxIdx[0], AABBmaxIdx[2]});
        DivideBlock2D<T>(NX, NZ, bestIx * bestIz + restX, AABBXZ, GeoHelperVec);

        int miny_child = AABBminIdx[1];
        for (int iY = 0; iY < bestIy; iY++) {
          int Ny_child = (NY + bestIy - iY - 1) / bestIy;
          for (AABB<int, 2>& aabb : GeoHelperVec) {
            blocksvec.emplace_back(Vector<int, 3>{aabb.getMin()[0], miny_child, aabb.getMin()[1]},
                                   Vector<int, 3>{aabb.getMax()[0], miny_child + Ny_child - 1, aabb.getMax()[1]});
          }
          miny_child += Ny_child;
        }
        return;
      }

      // split in four cuboid

      int restX = rest / bestIy + 1;
      int Ny_child = 0;
      int miny_child = AABBminIdx[1];
      int splited_nY = (int)(NY * (T)((bestIx * bestIz + restX) * restY) / (T)blocknum);

      std::vector<AABB<int, 2>> GeoHelperVec0;
      AABB<int, 2> AABBXZ0(Vector<int, 2>{AABBminIdx[0], AABBminIdx[2]},
                            Vector<int, 2>{AABBmaxIdx[0], AABBmaxIdx[2]});
      DivideBlock2D<T>(NX, NZ, bestIx * bestIz + restX, AABBXZ0, GeoHelperVec0);

      for (int iY = 0; iY < restY; iY++) {
        Ny_child = (splited_nY + restY - iY - 1) / restY;
        for (AABB<int, 2>& aabb : GeoHelperVec0) {
          blocksvec.emplace_back(Vector<int, 3>{aabb.getMin()[0], miny_child, aabb.getMin()[1]},
                                 Vector<int, 3>{aabb.getMax()[0], miny_child + Ny_child - 1, aabb.getMax()[1]});
        }
        miny_child += Ny_child;
      }

      splited_nY = NY - splited_nY;
      restX = rest / bestIy;

      std::vector<AABB<int, 2>> GeoHelperVec1;
      AABB<int, 2> AABBXZ1(Vector<int, 2>{AABBminIdx[0], AABBminIdx[2]},
                            Vector<int, 2>{AABBmaxIdx[0], AABBmaxIdx[2]});
      DivideBlock2D<T>(NX, NZ, bestIx * bestIz + restX, AABBXZ1, GeoHelperVec1);

      miny_child = AABBminIdx[1];;
      for (int iY = 0; iY < bestIy - restY; iY++) {
        Ny_child = (splited_nY + bestIy - restY - iY - 1) / (bestIy - restY);
        for (AABB<int, 2>& aabb : GeoHelperVec1) {
          blocksvec.emplace_back(Vector<int, 3>{aabb.getMin()[0], miny_child, aabb.getMin()[1]},
                                 Vector<int, 3>{aabb.getMax()[0], miny_child + Ny_child - 1, aabb.getMax()[1]});
        }
        miny_child += Ny_child;
      }
      return;
    }
    // 1

    // add in x than in y direction
    // 2
    else if (nX > nY && nX > nZ) {
      int restY = rest % bestIy;
      // split in two cuboid
      if (restY == 0) {
        int restZ = rest / bestIy;

        std::vector<AABB<int, 2>> GeoHelperVec;
        AABB<int, 2> AABBXZ(Vector<int, 2>{AABBminIdx[0], AABBminIdx[2]},
                            Vector<int, 2>{AABBmaxIdx[0], AABBmaxIdx[2]});
        DivideBlock2D<T>(NX, NZ, bestIx * bestIz + restZ, AABBXZ, GeoHelperVec);

        int miny_child = AABBminIdx[1];
        for (int iY = 0; iY < bestIy; iY++) {
          int Ny_child = (NY + bestIy - iY - 1) / bestIy;
          for (AABB<int, 2>& aabb : GeoHelperVec) {
            blocksvec.emplace_back(Vector<int, 3>{aabb.getMin()[0], miny_child, aabb.getMin()[1]},
                                   Vector<int, 3>{aabb.getMax()[0], miny_child + Ny_child - 1, aabb.getMax()[1]});
          }
          miny_child += Ny_child;
        }
        return;
      }

      // split in four cuboid

      int restZ = rest / bestIy + 1;
      int Ny_child = 0;
      int miny_child = AABBminIdx[1];
      int splited_nY = (int)(NY * (T)((bestIx * bestIz + restZ) * restY) / (T)blocknum);

      std::vector<AABB<int, 2>> GeoHelperVec0;
      AABB<int, 2> AABBXZ0(Vector<int, 2>{AABBminIdx[0], AABBminIdx[2]},
                            Vector<int, 2>{AABBmaxIdx[0], AABBmaxIdx[2]});
      DivideBlock2D<T>(NX, NZ, bestIx * bestIz + restZ, AABBXZ0, GeoHelperVec0);

      for (int iY = 0; iY < restY; iY++) {
        Ny_child = (splited_nY + restY - iY - 1) / restY;
        for (AABB<int, 2>& aabb : GeoHelperVec0) {
          blocksvec.emplace_back(Vector<int, 3>{aabb.getMin()[0], miny_child, aabb.getMin()[1]},
                                 Vector<int, 3>{aabb.getMax()[0], miny_child + Ny_child - 1, aabb.getMax()[1]});
        }
        miny_child += Ny_child;
      }

      splited_nY = NY - splited_nY;
      restZ = rest / bestIy;

      std::vector<AABB<int, 2>> GeoHelperVec1;
      AABB<int, 2> AABBXZ1(Vector<int, 2>{AABBminIdx[0], AABBminIdx[2]},
                            Vector<int, 2>{AABBmaxIdx[0], AABBmaxIdx[2]});
      DivideBlock2D<T>(NX, NZ, bestIx * bestIz + restZ, AABBXZ1, GeoHelperVec1);

      miny_child = AABBminIdx[1];
      for (int iY = 0; iY < bestIy - restY; iY++) {
        Ny_child = (splited_nY + bestIy - restY - iY - 1) / (bestIy - restY);
        for (AABB<int, 2>& aabb : GeoHelperVec1) {
          blocksvec.emplace_back(Vector<int, 3>{aabb.getMin()[0], miny_child, aabb.getMin()[1]},
                                 Vector<int, 3>{aabb.getMax()[0], miny_child + Ny_child - 1, aabb.getMax()[1]});
        }
        miny_child += Ny_child;
      }
      return;
    }
    // 2

    // add in y than in x direction
    // 3
    else {
      int restX = rest % bestIx;
      // split in two cuboid
      if (restX == 0) {
        int restZ = rest / bestIx;
        
        std::vector<AABB<int, 2>> GeoHelperVec;
        AABB<int, 2> AABBZY(Vector<int, 2>{AABBminIdx[2], AABBminIdx[1]},
                            Vector<int, 2>{AABBmaxIdx[2], AABBmaxIdx[1]});
        DivideBlock2D<T>(NZ, NY, bestIz * bestIy + restZ, AABBZY, GeoHelperVec);

        int minx_child = AABBminIdx[0];
        for (int iX = 0; iX < bestIx; iX++) {
          int Nx_child = (NX + bestIx - iX - 1) / bestIx;
          for (AABB<int, 2>& aabb : GeoHelperVec) {
            blocksvec.emplace_back(Vector<int, 3>{minx_child, aabb.getMin()[1], aabb.getMin()[0]},
                                   Vector<int, 3>{minx_child + Nx_child - 1, aabb.getMax()[1], aabb.getMax()[0]});
          }
          minx_child += Nx_child;
        }
        return;
      }

      // split in four cuboid

      int restZ = rest / bestIx + 1;
      int Nx_child = 0;
      int minx_child = AABBminIdx[0];
      int splited_nX = (int)(NX * (T)((bestIz * bestIy + restZ) * restX) / (T)blocknum);

      std::vector<AABB<int, 2>> GeoHelperVec0;
      AABB<int, 2> AABBZY0(Vector<int, 2>{AABBminIdx[2], AABBminIdx[1]},
                          Vector<int, 2>{AABBmaxIdx[2], AABBmaxIdx[1]});
      DivideBlock2D<T>(NZ, NY, bestIz * bestIy + restZ, AABBZY0, GeoHelperVec0);

      for (int iX = 0; iX < restX; iX++) {
        Nx_child = (splited_nX + restX - iX - 1) / restX;
        for (AABB<int, 2>& aabb : GeoHelperVec0) {
          blocksvec.emplace_back(Vector<int, 3>{minx_child, aabb.getMin()[1], aabb.getMin()[0]},
                                 Vector<int, 3>{minx_child + Nx_child - 1, aabb.getMax()[1], aabb.getMax()[0]});
        }
        minx_child += Nx_child;
      }

      splited_nX = NX - splited_nX;
      restZ = rest / bestIx;

      std::vector<AABB<int, 2>> GeoHelperVec1;
      AABB<int, 2> AABBZY1(Vector<int, 2>{AABBminIdx[2], AABBminIdx[1]},
                          Vector<int, 2>{AABBmaxIdx[2], AABBmaxIdx[1]});
      DivideBlock2D<T>(NZ, NY, bestIz * bestIy + restZ, AABBZY1, GeoHelperVec1);

      minx_child = AABBminIdx[0];
      for (int iX = 0; iX < bestIx - restX; iX++) {
        Nx_child = (splited_nX + bestIx - restX - iX - 1) / (bestIx - restX);
        for (AABB<int, 2>& aabb : GeoHelperVec1) {
          blocksvec.emplace_back(Vector<int, 3>{minx_child, aabb.getMin()[1], aabb.getMin()[0]},
                                 Vector<int, 3>{minx_child + Nx_child - 1, aabb.getMax()[1], aabb.getMax()[0]});
        }
        minx_child += Nx_child;
      }
      return;
    }
  }
}


// calculate the Standard Deviation of the number of points in each block
template <typename T, unsigned int D>
T ComputeStdDev(const std::vector<BasicBlock<T, D>> &Blocks) {
  T mean{};
  for (const BasicBlock<T, D> &block : Blocks) {
    mean += block.getN();
  }
  mean /= Blocks.size();
  T stdDev{};
  for (const BasicBlock<T, D> &block : Blocks) {
    stdDev += std::pow((static_cast<T>(block.getN()) / mean - T(1)), 2);
  }
  stdDev = std::sqrt(stdDev / Blocks.size());
  return stdDev;
}


template <typename T, unsigned int D>
class BlockGeometryHelperBase : public BasicBlock<T, D> {
  protected:
  // base block
  BasicBlock<T, D> _BaseBlock;

  // rectanglar blocks
  std::vector<BasicBlock<T, D>> _BasicBlocks0;
  // store old geometry info for data transfer, blockcomms are not needed here
  // simply store basicblocks is ok since block could be constructed from basicblock
  std::vector<BasicBlock<T, D>> _BasicBlocks1;
  // exchange flag for _BasicBlocks0 and _BasicBlocks1
  bool _Exchanged;

  // for mpi
  std::vector<std::vector<int>> _BlockIndex0;
  std::vector<std::vector<int>> _BlockIndex1;
  std::vector<std::vector<BasicBlock<T, D>*>> _BlockIndexPtr0;
  std::vector<std::vector<BasicBlock<T, D>*>> _BlockIndexPtr1;
  bool _IndexExchanged;

  std::size_t getTotalBaseCellNum() {
    std::size_t sum = 0;
    for (BasicBlock<T, 2>& block : getAllBasicBlocks()) {
      sum += block.getN();
    }
    return sum;
  }

  // get all new basic blocks
  std::vector<BasicBlock<T, 2>>& getAllBasicBlocks() {
    if (_Exchanged)
      return _BasicBlocks1;
    else
      return _BasicBlocks0;
  }
  BasicBlock<T, 2>& getAllBasicBlock(int id) { return getAllBasicBlocks()[id]; }
  const BasicBlock<T, 2>& getAllBasicBlock(int id) const {
    return getAllBasicBlocks()[id];
  }
  // get all old basic blocks
  std::vector<BasicBlock<T, 2>>& getAllOldBasicBlocks() {
    if (_Exchanged)
      return _BasicBlocks0;
    else
      return _BasicBlocks1;
  }
  BasicBlock<T, 2>& getAllOldBasicBlock(int id) { return getAllOldBasicBlocks()[id]; }
  const BasicBlock<T, 2>& getAllOldBasicBlock(int id) const {
    return getAllOldBasicBlocks()[id];
  }


  std::vector<std::vector<int>>& getAllBlockIndices() {
    if (_IndexExchanged)
      return _BlockIndex1;
    else
      return _BlockIndex0;
  }
  std::vector<std::vector<int>>& getAllOldBlockIndices() {
    if (_IndexExchanged)
      return _BlockIndex0;
    else
      return _BlockIndex1;
  }


  std::vector<std::vector<BasicBlock<T, 2>*>>& getAllBlockIndexPtrs() {
    if (_IndexExchanged)
      return _BlockIndexPtr1;
    else
      return _BlockIndexPtr0;
  }
  // get basic blocks from std::vector<std::vector<BasicBlock<T, 2>*>> _BlockIndexPtr
  std::vector<BasicBlock<T, 2>*>& getBasicBlocks() {
    return getAllBlockIndexPtrs()[mpi().getRank()];
  }
  BasicBlock<T, 2>& getBasicBlock(int id) { return *(getBasicBlocks()[id]); }
  const BasicBlock<T, 2>& getBasicBlock(int id) const { return *(getBasicBlocks()[id]); }

  std::vector<std::vector<BasicBlock<T, 2>*>>& getAllOldBlockIndexPtrs() {
    if (_IndexExchanged)
      return _BlockIndexPtr0;
    else
      return _BlockIndexPtr1;
  }
  std::vector<BasicBlock<T, 2>*>& getOldBasicBlocks() {
    return getAllOldBlockIndexPtrs()[mpi().getRank()];
  }
  BasicBlock<T, 2>& getOldBasicBlock(int id) { return *(getOldBasicBlocks()[id]); }
  const BasicBlock<T, 2>& getOldBasicBlock(int id) const {
    return *(getOldBasicBlocks()[id]);
  }


};