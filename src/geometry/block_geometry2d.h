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

// block_geometry2d.h

#pragma once

#include "geometry/geometry2d.h"

// basic block 2d structure, BasicBlock stores the original AABB and index AABB(not
// extended)

template <typename T>
class Block2D : public BasicBlock<T, 2> {
 private:
  // base block
  BasicBlock<T, 2> _BaseBlock;

  // --- communication structure ---
  // neighbor block
  std::vector<Block2D<T>*> _Neighbors;
  // conmmunicate with same level block
  std::vector<BlockComm<T, 2>> Communicators;
  // average block comm, get from higher level block
  std::vector<InterpBlockComm<T, 2>> AverageComm;
  // interp block comm, get from lower level block
  std::vector<InterpBlockComm<T, 2>> InterpComm;

  // overlap
  int _overlap;

#ifdef MPI_ENABLED
  int _Rank;
  bool _NeedMPIComm;

  MPIBlockCommStru MPIComm;
  MPIInterpBlockCommStru<2> MPIInterpComm;
  // for flag comm
  MPIBlockBuffer<std::uint8_t> MPIBuffer;
#endif

 public:
  // construct directly from basicblock
  Block2D(const BasicBlock<T, 2>& baseblock, int olap = 1);

  // block2d for uniform block structure
  Block2D(const AABB<T, 2>& block, const AABB<int, 2>& idxblock, int blockid,
          T voxelSize = T(1), int olap = 1);
#ifdef MPI_ENABLED
  Block2D(int rank, int blockid, const AABB<T, 2>& block, const AABB<int, 2>& idxblock,
          T voxelSize = T(1), int olap = 1);
#endif

  ~Block2D() = default;

  const BasicBlock<T, 2>& getSelfBlock() const { return *this; }
  const BasicBlock<T, 2>& getBaseBlock() const { return _BaseBlock; }
  BasicBlock<T, 2>& getBaseBlock() { return _BaseBlock; }

  int getOverlap() const { return _overlap; }

  // setup boundary
  template <typename FieldType, typename datatype, typename LatSet>
  void SetupBoundary(const AABB<T, 2>& block, FieldType& field, datatype bdvalue);

  // refine block
  void Refine(std::uint8_t deltalevel = std::uint8_t(1));
  // coarsen block
  void Coarsen(std::uint8_t deltalevel = std::uint8_t(1));

  std::vector<Block2D<T>*>& getNeighbors() { return _Neighbors; }
  const Block2D<T>& getNeighbor(int id) const { return *_Neighbors[id]; }

  std::vector<BlockComm<T, 2>>& getCommunicators() { return Communicators; }
  BlockComm<T, 2>& getCommunicator(int id) { return Communicators[id]; }

  std::vector<InterpBlockComm<T, 2>>& getInterpBlockComm() { return InterpComm; }
  InterpBlockComm<T, 2>& getInterpBlockComm(int id) { return InterpComm[id]; }

  std::vector<InterpBlockComm<T, 2>>& getAverageBlockComm() { return AverageComm; }
  InterpBlockComm<T, 2>& getAverageBlockComm(int id) { return AverageComm[id]; }

#ifdef MPI_ENABLED
  int getRank() const { return _Rank; }
  bool getNeedMPIComm() const { return _NeedMPIComm; }

  MPIBlockCommStru& getMPIBlockComm() { return MPIComm; }
  const MPIBlockCommStru& getMPIBlockComm() const { return MPIComm; }

  MPIBlockBuffer<std::uint8_t>& getMPIBlockBuffer() { return MPIBuffer; }
  const MPIBlockBuffer<std::uint8_t>& getMPIBlockBuffer() const { return MPIBuffer; }

  // mpi communicate FlagField GeometryFlag, non blocking
  void GeoFlagMPIComm();
#endif
};

template <typename T>
class BlockGeometryHelper2D;

// block geometry
template <typename T>
class BlockGeometry2D : public BasicBlock<T, 2> {
 private:
  // base block
  BasicBlock<T, 2> _BaseBlock;
  // TODO: _BlockAABBs may be removed
  std::vector<AABB<int, 2>> _BlockAABBs;
  std::vector<Block2D<T>> _Blocks;

  // ext(overlap) of the whole domain
  int _overlap;
  // max level
  std::uint8_t _MaxLevel;

 public:
  BlockGeometry2D(int Nx, int Ny, int blocknum, const AABB<T, 2>& block,
                  T voxelSize = T(1), int overlap = 1, bool refine = false);
  BlockGeometry2D(BlockGeometryHelper2D<T>& GeoHelper);
  ~BlockGeometry2D() = default;

  void Init(BlockGeometryHelper2D<T>& GeoHelper);

  void UpdateMaxLevel();
  inline std::uint8_t getMaxLevel() const { return _MaxLevel; }

  int getBlockNum() const { return _Blocks.size(); }
  // get total number of cells, overlapped cells included
  int getTotalCellNum() const;
  // get total number of cells, overlapped cells not included
  int getBaseCellNum() const;

  // first divide blockgeometry into blocks, stored in _BlockAABBs
  void DivideBlocks(int blocknum);
  void CreateBlocks();
  void SetupNbrs();
  void InitCommunicators();
  // communicate GeometryFlag
  void GeoFlagComm();
  // Communicators with different level
  // Low level block(coarse) get info from High level block(fine) using average
  void InitAverComm(int highlevelovlap = 2);
  // experimental
  // void InitAverComm2();
  // High level block(fine) get info from Low level block(coarse) using interpolation
  void InitIntpComm();
  // High level block(fine) get info from Low level block(coarse) using extrapolation
  // only support one single layer of overlapped cells
  // void InitExtpComm();
  // init all commuicators
  void InitAllComm();

  const BasicBlock<T, 2>& getSelfBlock() const { return *this; }
  const BasicBlock<T, 2>& getBaseBlock() const { return _BaseBlock; }

  std::vector<Block2D<T>>& getBlocks() { return _Blocks; }
  const std::vector<Block2D<T>>& getBlocks() const { return _Blocks; }

  Block2D<T>& getBlock(int id) { return _Blocks[id]; }
  const Block2D<T>& getBlock(int id) const { return _Blocks[id]; }

  // get size of each block
  std::vector<std::size_t> getBlockSizes() const {
    std::vector<std::size_t> sizes;
    for (const Block2D<T>& block : _Blocks) sizes.push_back(block.getN());
    return sizes;
  }
};

enum BlockCellTag : std::uint8_t { None = 1, Refine = 2, Coarsen = 4, Solid = 8 };

// especially designed for refine blockgeometry, uniform blockgeometry does not need this
// all BasicBlocks here refer to the base block, overlaps will not be handled here
template <typename T>
class BlockGeometryHelper2D : public BasicBlock<T, 2> {
 private:
  // base block
  BasicBlock<T, 2> _BaseBlock;
  // resulting basic block Cells, each cell containing several points
  std::vector<BasicBlock<T, 2>> _BlockCells;
  // tags for BlockCells
  std::vector<BlockCellTag> _BlockCellTags;
  // rectanglar blocks
  std::vector<BasicBlock<T, 2>> _BasicBlocks0;
  // cell geometry info
  int CellsNx;
  int CellsNy;
  int CellsN;
  // block length
  int BlockLen;
  // extension of the whole domain
  int Ext;
  // max level
  std::uint8_t _MaxLevel;
  // store old geometry info for data transfer, blockcomms are not needed here
  // simply store basicblocks is ok since block could be constructed from basicblock
  std::vector<BasicBlock<T, 2>> _BasicBlocks1;
  // exchange flag for _BasicBlocks0 and _BasicBlocks1
  bool _Exchanged;

 public:
  // domain of Nx * Ny will be divided into (Nx/blocklen)*(Ny/blocklen) blocks
  BlockGeometryHelper2D(int Nx, int Ny, int blocklen, const AABB<T, 2>& AABBs,
                        T voxelSize = T(1), int ext = 1);
  ~BlockGeometryHelper2D() = default;

  // get
  int getExt() const { return Ext; }
  BasicBlock<T, 2>& getBaseBlock() { return _BaseBlock; }
  BasicBlock<T, 2>& getBlockCell(int id) { return _BlockCells[id]; }
  std::vector<BasicBlock<T, 2>>& getBlockCells() { return _BlockCells; }

  BasicBlock<T, 2>& getBasicBlock(int id) { return getBasicBlocks()[id]; }
  const BasicBlock<T, 2>& getBasicBlock(int id) const { return getBasicBlocks()[id]; }

  BasicBlock<T, 2>& getOldBasicBlock(int id) { return getOldBasicBlocks()[id]; }
  const BasicBlock<T, 2>& getOldBasicBlock(int id) const { return getOldBasicBlocks()[id]; }

  // get new basic blocks
  std::vector<BasicBlock<T, 2>>& getBasicBlocks() {
    if (_Exchanged)
      return _BasicBlocks1;
    else
      return _BasicBlocks0;
  }

  // get old basic blocks
  std::vector<BasicBlock<T, 2>>& getOldBasicBlocks() {
    if (_Exchanged)
      return _BasicBlocks0;
    else
      return _BasicBlocks1;
  }

  BlockCellTag& getBlockCellTag(int id) { return _BlockCellTags[id]; }
  std::vector<BlockCellTag>& getBlockCellTags() { return _BlockCellTags; }


  int getCellsNx() const { return CellsNx; }
  int getCellsNy() const { return CellsNy; }

  void CreateBlockCells();
  // create block from BlockCells, this should be called after refinement
  void CreateBlocks();
  // check refine cell status
  void PostRefine();

  void UpdateMaxLevel();
  std::uint8_t getMaxLevel() const { return _MaxLevel; }

  // lambda function for each block
  template <typename Func>
  void forEachBlockCell(Func func);

  // optimal procNum usually set to actual procNum, 
  // MaxProcNum usually set to x * optProcNum
  void AdaptiveOptimization(int OptProcNum, int MaxProcNum = -1, bool enforce = true);
  // void AdaptiveOptimization(std::vector<BasicBlock<T, 2>>& Blocks, int OptProcNum,
  //                           int MaxProcNum = -1, bool enforce = true);
  // optimize for parallel computing
  void Optimize(int ProcessNum, bool enforce = true);
  void Optimize(std::vector<BasicBlock<T, 2>>& Blocks, int ProcessNum,
                bool enforce = true);
  // calculate the Standard Deviation of the number of points in each block
  T ComputeStdDev() const;
  T ComputeStdDev(const std::vector<BasicBlock<T, 2>>& Blocks) const;
};

// helper class for MPI
template <typename T>
class BlockGeometryMPIHelper2D : public BasicBlock<T, 2> {
 private:
  // base block
  BasicBlock<T, 2> _BaseBlock;
  // blocks
  int _BlockNum;
  // for each process
  std::vector<AABB<int, 2>> _BlockAABBs;
  // base blocks
  std::vector<BasicBlock<T, 2>> _Blocks;

  int _overlap;

  std::vector<std::vector<int>> _BlockNbrRanks;

 public:
  BlockGeometryMPIHelper2D(int Nx, int Ny, int blocknum, const AABB<T, 2>& AABBs,
                           T voxelSize = T(1), int overlap = 1, bool refine = false);
  ~BlockGeometryMPIHelper2D() = default;

  // get
  int getBlockNum() const { return _BlockNum; }
  const BasicBlock<T, 2>& getBaseBlock() { return _BaseBlock; }
  const BasicBlock<T, 2>& getSelfBlock(int id) { return *this; }

  // get AABB<int, 2>& in std::vector<AABB<int, 2>>
  const AABB<int, 2>& getBlockAABB(int id) const { return _BlockAABBs[id]; }
  // get std::vector<AABB<int, 2>>&
  std::vector<AABB<int, 2>>& getBlockAABBs() { return _BlockAABBs; }

  // get BasicBlock<T, 2>& in std::vector<BasicBlock<T, 2>>
  const BasicBlock<T, 2>& getBlock(int id) const { return _Blocks[id]; }
  // get std::vector<BasicBlock<T, 2>>&
  std::vector<BasicBlock<T, 2>>& getBlocks() { return _Blocks; }

  // get std::vector<int>& int std::vector<std::vector<int>>
  const std::vector<int>& getBlockNbrRanks(int id) const { return _BlockNbrRanks[id]; }
  // get std::vector<std::vector<int>>&
  std::vector<std::vector<int>>& getBlockNbrRanks() { return _BlockNbrRanks; }

  // the following can be called on each process
  void DivideBlocks(int blocknum);
  void CreateBlocks();
  void SetupNbrs();
  // the following should be called on specific process
  void InitMPIBlockCommStru(MPIBlockCommStru& BlockComm);

  void InitMPIAverComm(int highlevelovlap = 2);

  void InitMPIIntpComm();

  void InitAllMPIComm();
  // void GeoFlagComm();
};

// block geometry manager for MPI
// divide blocks for MPI, each process has its own BlockGeometry which manages several
// blocks template <typename T> class BlockGeometryManager2D : public BasicBlock<T, 2> {
//  private:
//   // base block
//   BasicBlock<T, 2> BaseBlock;
//   // blocks
//   int ToatlBlockNum;
//   // mpi blocks num
//   int MPIBlockNum;
//   std::vector<AABB<int, 2>> Blocks;
//   std::vector<int> BlockNumsEachProcess;
//   // flag
//   std::uint8_t _AABBflag;
//   std::uint8_t _voidflag;
//  public:
//   BlockGeometryManager2D(int Nx, int Ny, int blocknum, const AABB<T, 2>& AABBs, T
//   voxelSize = T(1),
//                          std::uint8_t AABBflag = std::uint8_t(1),
//                          std::uint8_t voidflag = std::uint8_t(0));
//   ~BlockGeometryManager2D() = default;
//   void DivideBlocks(int processnum);
// };
