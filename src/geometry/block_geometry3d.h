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

#include "geometry/basic_geometry.h"
#include "io/stlreader.h"
#include "io/block_reader.h"

template <typename T>
class Block3D : public BasicBlock<T, 3> {
 private:
  // base block
  BasicBlock<T, 3> _BaseBlock;

  // --- communication structure ---
  Communicator _Comm;
  // neighbor block
  std::vector<Block3D<T>*> _Neighbors;

  // overlap
  int _overlap;

 public:
  // constructors
  // default aabbflag = 1, voidflag = 0
  Block3D(const BasicBlock<T, 3>& baseblock, int olap = 1);

  Block3D(const AABB<T, 3>& block, const AABB<int, 3>& idxblock, int blockid,
          T voxelSize = T{1}, int olap = 1);

  ~Block3D() = default;

  const BasicBlock<T, 3>& getSelfBlock() const { return *this; }
  const BasicBlock<T, 3>& getBaseBlock() const { return _BaseBlock; }
  BasicBlock<T, 3>& getBaseBlock() { return _BaseBlock; }

  int getOverlap() const { return _overlap; }
  
  // this is a generalized function to iterate over all cells in _BaseBlock
  // TODO: change similar code to this?
  template <typename Func>
  void forEachBaseCell(const Func& func);
  // clean lonely flags
  template <typename FieldType, typename LatSet>
  void CleanLonelyFlags(FieldType& field, std::uint8_t flag, std::uint8_t voidflag, 
    unsigned int lonelyth, bool& cleaned, std::size_t& count);

  // setup boundary
  template <typename FieldType, typename LatSet>
  void SetupBoundary(const AABB<T, 3>& block, FieldType& field,
                     typename FieldType::value_type bdvalue);
  template <typename FieldType, typename LatSet>
  void SetupBoundary(FieldType& field, typename FieldType::value_type fromvalue,
  typename FieldType::value_type voidvalue, typename FieldType::value_type bdvalue);
  
  template <typename FieldType>
  void ReadOctree(Octree<T>* tree, FieldType& field, typename FieldType::value_type stlflag);

  std::vector<Block3D<T>*>& getNeighbors() { return _Neighbors; }
  const Block3D<T>& getNeighbor(int id) const { return *_Neighbors[id]; }

  Communicator& getCommunicator() { return _Comm; }
  const Communicator& getCommunicator() const { return _Comm; }
};

template <typename T>
class BlockGeometryHelper3D;

template <typename T>
class BlockGeometry3D : public BasicBlock<T, 3> {
 private:
  // base block
  BasicBlock<T, 3> _BaseBlock;
  // TODO: _BlockAABBs may be removed
  std::vector<AABB<int, 3>> _BlockAABBs;
  // _BasicBlocks may be removed, now is used in some functions take a vector of BasicBlock
  std::vector<BasicBlock<T, 3>> _BasicBlocks;
  // info of blockaabbs and basicblocks is contained in _Blocks
  std::vector<Block3D<T>> _Blocks;
  // block id - element index in std::vector<Block3D<T>> _Blocks;
  std::unordered_map<int, std::size_t> _BlockIndexMap;

  // ext(overlap) of the whole domain
  int _overlap;
  // max level
  std::uint8_t _MaxLevel;

 public:
  // construct uniform blockgeometry
  BlockGeometry3D(int Nx, int Ny, int Nz, int blocknum, const AABB<T, 3>& block,
                  T voxelSize = T(1), int overlap = 1,
                  int blockXNum = 0, int blockYNum = 0, int blockZNum = 0);
  // construct uniform/ refined blockgeometry from GeoHelper
  BlockGeometry3D(BlockGeometryHelper3D<T>& GeoHelper, bool useHelperOlap = true);
  // construct blockgeometry from blockreader
  BlockGeometry3D(const BlockReader<T,3>& blockreader, bool useReaderOlap = true);
  // construct blockgeometry from stlreader, 
  // not recommended, use BlockGeometryHelper3D with StlReader instead, may be removed?
  BlockGeometry3D(const StlReader<T>& reader, int blocknum);
  ~BlockGeometry3D() = default;

  void PrintInfo() const;
  void Init(BlockGeometryHelper3D<T>& GeoHelper);

  const BasicBlock<T, 3>& getSelfBlock() const { return *this; }
  const BasicBlock<T, 3>& getBaseBlock() const { return _BaseBlock; }

  std::vector<Block3D<T>>& getBlocks() { return _Blocks; }
  const std::vector<Block3D<T>>& getBlocks() const { return _Blocks; }

  const std::vector<BasicBlock<T, 3>>& getBasicBlocks() const { return _BasicBlocks; }

  Block3D<T>& getBlock(int id) { return _Blocks[id]; }
  const Block3D<T>& getBlock(int id) const { return _Blocks[id]; }

  // mpi: return _BlockIndexMap[id];
  std::size_t findBlockIndex(int id) const { 
#ifdef MPI_ENABLED
    return _BlockIndexMap.at(id); 
#else
    return id;
#endif
  }
  bool hasBlock(int id) const {
    return _BlockIndexMap.find(id) != _BlockIndexMap.end();
  }

  inline std::uint8_t getMaxLevel() const { return _MaxLevel; }

  int getBlockNum() const { return _Blocks.size(); }
  // get total number of cells, overlapped cells included
  std::size_t getTotalCellNum() const;
  // get total number of cells, overlapped cells not included
  std::size_t getBaseCellNum() const;

  // first divide blockgeometry into blocks, stored in _BlockAABBs
  void DivideBlocks(int blocknum, int blockXNum = 0, int blockYNum = 0, int blockZNum = 0);
  void CreateBlocks(int blocknum, int blockXNum = 0, int blockYNum = 0, int blockZNum = 0);
  void SetupNbrs();

  // --- communication ---
  void InitComm();
  // Low level block(coarse) get info from High level block(fine) using average
  void InitAverComm();
  // High level block(fine) get info from Low level block(coarse) using interpolation
  void InitIntpComm();
  // init all commuicators
  void InitAllComm();

// --- mpi communication ---
#ifdef MPI_ENABLED
  // construct uniform/ refined blockgeometry from GeoHelper and a vector of BasicBlock<T, 2>
  // for mpi communication with direction info 
  BlockGeometry3D(BlockGeometryHelper3D<T>& GeoHelper, 
  std::vector<BasicBlock<T, 3>>& BasicBlocks);

  void InitMPIComm(BlockGeometryHelper3D<T>& GeoHelper);
  // Low level block(coarse) get info from High level block(fine) using average
  void InitMPIAverComm(BlockGeometryHelper3D<T>& GeoHelper);
  // High level block(fine) get info from Low level block(coarse) using interpolation
  void InitMPIIntpComm(BlockGeometryHelper3D<T>& GeoHelper);
  // init all commuicators
  void InitAllMPIComm(BlockGeometryHelper3D<T>& GeoHelper);

#endif

private:
  void AddtoSharedAverComm(std::vector<SharedComm> &Comms, const BasicBlock<T, 3>& baseblock_extx, Block3D<T>& block, Block3D<T> *nblock);
  void AddtoSharedIntpComm(std::vector<SharedComm> &Comms, Block3D<T>& block, Block3D<T> *nblock);
  void BuildBlockIndexMap();
};

template <typename T>
class BlockGeometryHelper3D : public BasicBlock<T, 3> {
  private:
  // base block
  BasicBlock<T, 3> _BaseBlock;
  // resulting basic block Cells, each cell containing several points
  std::vector<BasicBlock<T, 3>> _BlockCells;
  // tags for BlockCells
  std::vector<BlockCellTag> _BlockCellTags;
  // rectanglar blocks
  std::vector<BasicBlock<T, 3>> _BasicBlocks0;
  // store old geometry info for data transfer, blockcomms are not needed here
  // simply store basicblocks is ok since block could be constructed from basicblock
  std::vector<BasicBlock<T, 3>> _BasicBlocks1;
  // cell geometry info
  int CellsNx;
  int CellsNy;
  int CellsNz;
  int CellsN;
  // block length
  int BlockCellLen;
  // overlap of blocks
  int _Overlap;
  // extension of the whole domain
  int _Ext;
  // max level limit
  std::uint8_t _LevelLimit;
  // max level
  std::uint8_t _MaxLevel;
  // exchange flag for _BasicBlocks0 and _BasicBlocks1
  bool _Exchanged;
  // delta index for cell block
  std::array<int, 26> Delta_Cellidx;

  // for mpi
  std::vector<std::vector<int>> _BlockIndex0;
  std::vector<std::vector<int>> _BlockIndex1;
  std::vector<std::vector<BasicBlock<T, 3>*>> _BlockIndexPtr0;
  std::vector<std::vector<BasicBlock<T, 3>*>> _BlockIndexPtr1;
  bool _IndexExchanged;
#ifdef MPI_ENABLED
  // mpi neighbors: first: rank, second: blockid
  std::vector<std::vector<std::pair<int, int>>> _MPIBlockNbrs;
  // for mpi direction communication init
  std::unique_ptr<BlockGeometry3D<T>> _BlockGeometry3D;
#endif


 public:
  // domain of Nx * Ny will be divided into (Nx/blocklen)*(Ny/blocklen) blocks
  BlockGeometryHelper3D(int Nx, int Ny, int Nz, const AABB<T, 3>& AABBs, T voxelSize = T(1),
                        int blockcelllen = 10, int olap = 1, int ext = 0, 
                        std::uint8_t llimit = std::uint8_t(2));
  // olap is the num of overlapped cells between blocks
  // ext is the extension of the whole domain, used for creating solid cells after reading stl
  BlockGeometryHelper3D(const StlReader<T>& reader, int blockcelllen = 10, 
                        int olap = 1, int ext = 0, std::uint8_t llimit = std::uint8_t(2));
  // construct blockgeometryhelper from blockreader
  BlockGeometryHelper3D(const BlockReader<T,3>& blockreader, bool useReaderOlap = true);
  ~BlockGeometryHelper3D() = default;

  // get
  void UpdateMaxLevel();
  std::uint8_t getMaxLevel() const { return _MaxLevel; }
  std::uint8_t getLevelLimit() const { return _LevelLimit; }
  const std::array<int, 26>& getDeltaCellidx() const { return Delta_Cellidx; }

  int getOverlap() const { return _Overlap; }
  BasicBlock<T, 3>& getBaseBlock() { return _BaseBlock; }

  BasicBlock<T, 3>& getBlockCell(int id) { return _BlockCells[id]; }
  std::vector<BasicBlock<T, 3>>& getBlockCells() { return _BlockCells; }

  BlockCellTag& getBlockCellTag(int id) { return _BlockCellTags[id]; }
  std::vector<BlockCellTag>& getBlockCellTags() { return _BlockCellTags; }

  int getCellsNx() const { return CellsNx; }
  int getCellsNy() const { return CellsNy; }
  int getCellsNz() const { return CellsNz; }

  std::size_t getTotalBaseCellNum() {
    std::size_t sum = 0;
    for (BasicBlock<T, 3>& block : getAllBasicBlocks()) {
      sum += block.getN();
    }
    return sum;
  }

  // get all new basic blocks
  std::vector<BasicBlock<T, 3>>& getAllBasicBlocks() {
    if (_Exchanged) return _BasicBlocks1;
    else return _BasicBlocks0;
  }
  BasicBlock<T, 3>& getAllBasicBlock(std::size_t id) { return getAllBasicBlocks()[id]; }
  const BasicBlock<T, 3>& getAllBasicBlock(std::size_t id) const {
    return getAllBasicBlocks()[id];
  }
  // get all old basic blocks
  std::vector<BasicBlock<T, 3>>& getAllOldBasicBlocks() {
    if (_Exchanged) return _BasicBlocks0;
    else return _BasicBlocks1;
  }
  BasicBlock<T, 3>& getAllOldBasicBlock(std::size_t id) { return getAllOldBasicBlocks()[id]; }
  const BasicBlock<T, 3>& getAllOldBasicBlock(std::size_t id) const {
    return getAllOldBasicBlocks()[id];
  }


  std::vector<std::vector<int>>& getAllBlockIndices() {
    if (_IndexExchanged) return _BlockIndex1;
    else return _BlockIndex0;
  }
  std::vector<std::vector<int>>& getAllOldBlockIndices() {
    if (_IndexExchanged) return _BlockIndex0;
    else return _BlockIndex1;
  }


  std::vector<std::vector<BasicBlock<T, 3>*>>& getAllBlockIndexPtrs() {
    if (_IndexExchanged) return _BlockIndexPtr1;
    else return _BlockIndexPtr0;
  }
  // get basic blocks from std::vector<std::vector<BasicBlock<T, 3>*>> _BlockIndexPtr
  std::vector<BasicBlock<T, 3>*>& getBasicBlocks() {
    return getAllBlockIndexPtrs()[mpi().getRank()];
  }
  BasicBlock<T, 3>& getBasicBlock(int id) { return *(getBasicBlocks()[id]); }
  const BasicBlock<T, 3>& getBasicBlock(int id) const { return *(getBasicBlocks()[id]); }

  std::vector<std::vector<BasicBlock<T, 3>*>>& getAllOldBlockIndexPtrs() {
    if (_IndexExchanged) return _BlockIndexPtr0;
    else return _BlockIndexPtr1;
  }
  std::vector<BasicBlock<T, 3>*>& getOldBasicBlocks() {
    return getAllOldBlockIndexPtrs()[mpi().getRank()];
  }
  BasicBlock<T, 3>& getOldBasicBlock(int id) { return *(getOldBasicBlocks()[id]); }
  const BasicBlock<T, 3>& getOldBasicBlock(int id) const {
    return *(getOldBasicBlocks()[id]);
  }

#ifdef MPI_ENABLED
  std::vector<std::pair<int, int>>& getMPIBlockNbrs(int blockid) {
    return _MPIBlockNbrs[blockid];
  }
  const BlockGeometry3D<T>& getBlockGeometry3D() const { return *_BlockGeometry3D; }
#endif

  void CreateBlockCells();
  void TagBlockCells(const StlReader<T>& reader);
  // find if one cell of the block is inside the octree
  bool IsInside(Octree<T>* tree, const BasicBlock<T, 3> &block) const;
  // find if one cell of the block is outside the octree
  bool hasOutSideCell(Octree<T>* tree, const BasicBlock<T, 3> &block) const;
  // create block from BlockCells, this should be called after refinement
  void CreateBlocks(bool CreateFromInsideTag = false);
  // create blocks manually, if used, AdaptiveOptimization() is NOT necessary
  void CreateBlocks(int blockXNum, int blockYNum, int blockZNum);
  // tag neighbor refine cells
  void TagRefineLayer(std::vector<std::uint8_t>& refine, bool& refined);
  // check refine cell status
  void CheckRefine();
  // perform refine
  void Refine();

  // shrink created BasicBlocks to fit geometry held by octree
  // this should be called after CreateBlocks(bool CreateFromInsideTag);
  void ShrinkBasicBlocks(const StlReader<T>& reader);

  // lambda function for each cell block
  template <typename Func>
  void forEachBlockCell(const Func& func);

  // optimal procNum usually set to actual procNum,
  // MaxProcNum usually set to x * optProcNum
  void AdaptiveOptimization(int OptProcNum, int MaxProcNum = -1, bool enforce = true);
  // void AdaptiveOptimization(std::vector<BasicBlock<T, 3>>& Blocks, int OptProcNum,
  //                           int MaxProcNum = -1, bool enforce = true);
  // optimize for parallel computing
  void Optimize(int ProcessNum, bool enforce = true, bool info = false);
  void Optimize(std::vector<BasicBlock<T, 3>>& Blocks, int ProcessNum,
                bool enforce = true);
	
	// optimize block's geometry to make each block has similar number of cells
	// this should be the FINAL step before LoadBalancing() and should be called after ShrinkBasicBlocks()
	// LoadOptimization() will NOT create or delete blocks, it changes the position and size of blocks
	// to get the best load result
  void LoadOptimization(int maxiter = 1000, T tolstddev = T(0.05));
	
	// greedy algorithm to balance the load of each process
  void LoadBalancing(int ProcessNum = mpi().getSize());

#ifdef MPI_ENABLED
  // call this after LoadBalancing()
  void SetupMPINbrs();
  // construct BlockGeometry2D from BlockGeometryHelper2D's this pointer
  void InitBlockGeometry3D();
  // find which rank the block belongs to
  int whichRank(int blockid);
#endif

};
