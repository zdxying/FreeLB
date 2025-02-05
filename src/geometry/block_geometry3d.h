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
  std::vector<BasicBlock<T, 3>>& BasicBlocks, bool useHelperOlap = true);

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
class GeometryFlagField {
  private:
  // geometry info
  int _Nx;
  int _Ny;
  int _Nz;
  int _NxNy;
  int _Ext;
  // shift 1 in each direction
  // avoid out-of-bound access when checking neighbors
  int _Shift;
  // voxel size
  T _VoxelSize;
  // min voxel position with extension
  Vector<T, 3> _MinVoxel;
  // flag field
  std::uint8_t* _FlagField = nullptr;

  void readOctree(Octree<T>* tree, std::uint8_t stlflag = std::uint8_t{2}) {
    for (int z = 0; z < _Nz; ++z) {
      for (int y = 0; y < _Ny; ++y) {
        for (int x = 0; x < _Nx; ++x) {
          const Vector<T, 3> vox = _MinVoxel + Vector<int, 3>{x, y, z} * _VoxelSize;
          Octree<T>* node = tree->find(vox);
          if (node != nullptr) {
            if (node->isLeaf() && node->getInside())
              _FlagField[z * _NxNy + y * _Nx + x + _Shift] = stlflag;
          }
        }
      }
    }
  }

  public:
  std::uint8_t _VoidFlag = std::uint8_t{1};
  GeometryFlagField() = default;
  GeometryFlagField(const StlReader<T>& reader, int ext = 1, 
    std::uint8_t stlflag = std::uint8_t{2}, 
    std::uint8_t voidflag = std::uint8_t{1}) {
      Init(reader, ext, stlflag, voidflag);
    }
  ~GeometryFlagField() {
    if (_FlagField) delete[] _FlagField;
  }

  void Init(const StlReader<T>& reader, int ext = 1, 
    std::uint8_t stlflag = std::uint8_t{2}, 
    std::uint8_t voidflag = std::uint8_t{1}) {
    _Nx = int(std::ceil(reader.getMesh().getMax_Min()[0] / reader.getVoxelSize())) + 2 * ext;
    _Ny = int(std::ceil(reader.getMesh().getMax_Min()[1] / reader.getVoxelSize())) + 2 * ext;
    _Nz = int(std::ceil(reader.getMesh().getMax_Min()[2] / reader.getVoxelSize())) + 2 * ext;
    _NxNy = _Nx * _Ny;
    _Ext = ext;
    _Shift = _NxNy + _Nx + 1;
    _VoxelSize = reader.getVoxelSize();
    _MinVoxel = reader.getMesh().getMin() - ext * _VoxelSize / 2;
    _VoidFlag = voidflag;
    if (_FlagField) delete[] _FlagField;
    const std::size_t actaulsize = (_Nx + 2) * (_Ny + 2) * (_Nz + 2);
    _FlagField = new std::uint8_t[actaulsize];
    std::fill(_FlagField, _FlagField + actaulsize, voidflag);
    readOctree(reader.getTree(), stlflag);
  }

  std::uint8_t get(int x, int y, int z) const {
    // boundary check
    if (x < 0 || x >= _Nx || y < 0 || y >= _Ny || z < 0 || z >= _Nz) return _VoidFlag;
    return _FlagField[z * _NxNy + y * _Nx + x + _Shift];
  }
  std::uint8_t operator[](const Vector<T, 3>& pos) const {
    const int x = int(std::round((pos[0] - _MinVoxel[0]) / _VoxelSize));
    const int y = int(std::round((pos[1] - _MinVoxel[1]) / _VoxelSize));
    const int z = int(std::round((pos[2] - _MinVoxel[2]) / _VoxelSize));
    return get(x, y, z);
  }
  void set(int x, int y, int z, std::uint8_t flag) {
    // boundary check
    if (x < 0 || x >= _Nx || y < 0 || y >= _Ny || z < 0 || z >= _Nz) return;
    _FlagField[z * _NxNy + y * _Nx + x + _Shift] = flag;
  }

  template <typename LatSet>
  void SetupBoundary(std::uint8_t fromvalue, std::uint8_t voidvalue, std::uint8_t bdvalue) {
    // temp flag field store the transition flag
    GenericArray<bool> TransFlag(_Nx*_Ny*_Nz, false);
    // use 0 index, _FlagField is already extended, accessing neighbor cells is safe
    for (int z = 0; z < _Nz; ++z) {
      for (int y = 0; y < _Ny; ++y) {
        for (int x = 0; x < _Nx; ++x) {
          if (this->get(x, y, z) == fromvalue) {
            for (unsigned int i = 1; i < LatSet::q; ++i) {
              const int nx = x + latset::c<LatSet>(i)[0];
              const int ny = y + latset::c<LatSet>(i)[1];
              const int nz = z + latset::c<LatSet>(i)[2];
              if (this->get(nx, ny, nz) == voidvalue) {
                const std::size_t id = z * _NxNy + y * _Nx + x;
                TransFlag.set(id, true);
                break;
              }
            }
          }
        }
      }
    }
  
    for (int z = 0; z < _Nz; ++z) {
      for (int y = 0; y < _Ny; ++y) {
        for (int x = 0; x < _Nx; ++x) {
          const std::size_t id = z * _NxNy + y * _Nx + x;
          if (TransFlag[id]) this->set(x, y, z, bdvalue);
        }
      }
    }
  }
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
  // 2025/1/18: _Ext is used to create boundary flag from void flag
  // set _Ext to 1, then at the whole domain's boundary(the boundary of the overall AABB),
  // additional 1 layer of cells will be created and tagged with VoidFlag
  // reading octree will NOT set them to AABBFlag
  // however this is NOT perfect for complex geometry:
  // when creating blocks from inside flag, the blockcell containing only voidflag(find no cell in octree)
  // will be neglected, it is possible that some of the rest blockcells lost the additional 1 layer of cells created by _Ext
  // this could be partially solved by using AddVoidCellLayer() after CreateBlocks()
  // however when the resulting BlockGeometry has hollows inside, and the size of the hollow does not match the adjacent block
  // AddVoidCellLayer() could not solve the problem
  // we may design another work flow for complex geometry:
  // 1. read stl, create octree
  // 2. read octree and store all kinds of flags in a buffer flag field
  //    note that the buffer flag field will NOT affect the subsequent real FlagField
  // 2. tag block cells, NOT from octree, from buffer flag field instead
  // 3. create blocks from tagged block cells, and do optimization...
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

  // pointer to stlreader
  const StlReader<T>* _Reader = nullptr;
  // helper flag field
  GeometryFlagField<T> _FlagField;

 public:
  // domain of Nx * Ny will be divided into (Nx/blocklen)*(Ny/blocklen) blocks
  BlockGeometryHelper3D(int Nx, int Ny, int Nz, const AABB<T, 3>& AABBs, T voxelSize = T(1),
                        int blockcelllen = 10, int olap = 1, int ext = 0, bool useblockcell = true,
                        std::uint8_t llimit = std::uint8_t(2));
  // olap is the num of overlapped cells between blocks
  // ext is the extension of the whole domain, used for creating solid cells after reading stl
  BlockGeometryHelper3D(const StlReader<T>& reader, int blockcelllen = 10, 
                        int olap = 1, int ext = 0, bool useblockcell = true, 
                        std::uint8_t llimit = std::uint8_t(2));
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
    std::size_t sum{};
    for (BasicBlock<T, 3>& block : getAllBasicBlocks()) sum += block.getN();
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

  GeometryFlagField<T>& getFlagField() { return _FlagField; }

  void CreateBlockCells();
  // tag from stlreader
  void TagBlockCells(const StlReader<T>& reader);
  // tag from helper flag field
  void TagBlockCells(std::uint8_t voidflag = std::uint8_t{1});
  // create block from BlockCells, this should be called after refinement
  void CreateBlocks(bool CreateFromInsideTag = false, bool outputinfo = true);
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
  void RemoveUnusedCells(const StlReader<T>& reader, bool outputinfo = true);
  // remove cells using helper flag field
  void RemoveUnusedCells(std::uint8_t voidflag = std::uint8_t{1}, bool outputinfo = true);
  // this works not well for complex geometry, may be removed
  void AddVoidCellLayer(const StlReader<T>& reader, bool outputinfo = true);

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

  T IterateAndOptimizeImp(int ProcNum, int BlockCellLen);
  // this will iterate all int blockcell length within the input range
  // each trial will divide(and merge and do LoadOptimization) the domain into ProcNum blocks
  // and calculate the standard deviation of the number of cells in each block
  // the best blockcell length will be the one with the smallest standard deviation
  void IterateAndOptimize(int ProcNum, int MinBlockCellLen = 10, int MaxBlockCellLen = 100, bool stepinfo = true);

  // recursive coordinate bisection along longest axis based on _FlagField
  // geometry is divided into exactly ProcNum blocks, ProcNum should be power of 2 for simplicity
  // this will NOT use the _BlockCells
  void RCBOptimization(int ProcNum, bool verbose = false);

  // (triple) recursive coordinate bisection along z - y - x axis based on _FlagField
  // one z - y - x bisection makes 8 blocks of similar length ratio with the original block
  // to get the minimal communication overhead
  // geometry is divided into exactly ProcNum blocks, ProcNum should be power of 8 for best result
  // this will NOT use the _BlockCells
	void TRCBOptimization(int ProcNum, bool verbose = false);

	// optimize block's geometry to make each block has similar number of cells
	// this should be the FINAL step before LoadBalancing() and should be called after RemoveUnusedCells()
	// LoadOptimization() will NOT create or delete blocks, it changes the position and size of blocks
	// to get the best load result
  void LoadOptimization(int maxiter = 1000, T tolstddev = T(0.05), bool outputinfo = true);
  
  // init GeoHelper using a (different) blockcelllen, 
	// this function is valid only if GeoHelper is constructed from stlreader
  void Init(int blockcelllen);
	
	// greedy algorithm to balance the load of each process
  void LoadBalancing(int ProcessNum = mpi().getSize());

#ifdef MPI_ENABLED
  // call this after LoadBalancing()
  void SetupMPINbrs();
  // construct BlockGeometry2D from BlockGeometryHelper2D's this pointer
  void InitBlockGeometry3D(bool useHelperOlap = true);
  // find which rank the block belongs to
  int whichRank(int blockid);
#endif

};
