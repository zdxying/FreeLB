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

// block_geometry2d.hh

#pragma once

#include "block_geometry2d.h"
#include "geometry/block_geometry2d.h"


template <typename T>
Block2D<T>::Block2D(const BasicBlock<T, 2> &baseblock, int olap)
    : _BaseBlock(baseblock), _overlap(olap),
      BasicBlock<T, 2>(baseblock.getExtBlock(olap)) {}


template <typename T>
Block2D<T>::Block2D(const AABB<T, 2> &block, const AABB<int, 2> &idxblock, int blockid,
                    T voxelSize, int olap)
    : BasicBlock<T, 2>(voxelSize, block.getExtended(Vector<T, 2>{voxelSize}),
                       idxblock.getExtended(Vector<int, 2>{1}), blockid),
      _BaseBlock(voxelSize, block, idxblock, blockid), _overlap(olap) {
  // read from AABBs
  // ReadAABBs(AABBs, AABBflag);
}

#ifdef MPI_ENABLED

template <typename T>
Block2D<T>::Block2D(int rank, int blockid, const AABB<T, 2> &block,
                    const AABB<int, 2> &idxblock, T voxelSize, int olap)
    : BasicBlock<T, 2>(voxelSize, block.getExtended(Vector<T, 2>{voxelSize}),
                       idxblock.getExtended(Vector<int, 2>{1}), blockid),
      _BaseBlock(voxelSize, block, idxblock, blockid), _overlap(olap), _Rank(rank),
      _voidflag(voidflag) {}

#endif

template <typename T>
template <typename FieldType, typename datatype, typename LatSet>
void Block2D<T>::SetupBoundary(const AABB<T, 2> &block, FieldType &field,
                               datatype bdvalue) {
  // temp flag field store the transition flag
  GenericArray<bool> TransFlag(BasicBlock<T, 2>::N, false);

  for (int y = _overlap; y < BasicBlock<T, 2>::Mesh[1] - _overlap; ++y) {
    for (int x = _overlap; x < BasicBlock<T, 2>::Mesh[0] - _overlap; ++x) {
      const Vector<int, 2> locidx{x, y};
      const Vector<T, 2> vox = BasicBlock<T, 2>::getVoxel(locidx);
      for (int i = 1; i < LatSet::q; ++i) {
        const Vector<T, 2> nvox = vox + LatSet::c[i] * BasicBlock<T, 2>::VoxelSize;
        if (!block.isInside(nvox)) {
          TransFlag.set(BasicBlock<T, 2>::getIndex(locidx), true);
          break;
        }
      }
    }
  }

  for (std::size_t id = 0; id < BasicBlock<T, 2>::N; ++id) {
    if (TransFlag[id]) field.SetField(id, bdvalue);
  }
#ifdef MPI_ENABLED
  GeoFlagMPIComm();
#endif
}

template <typename T>
void Block2D<T>::Refine(std::uint8_t deltalevel) {
  _BaseBlock = _BaseBlock.getRefinedBlock(deltalevel);
  _overlap = 2;
  BasicBlock<T, 2>::operator=(_BaseBlock.getExtBlock(_overlap));
}

template <typename T>
void Block2D<T>::Coarsen(std::uint8_t deltalevel) {
  _BaseBlock = _BaseBlock.getCoasenedBlock(deltalevel);
  if (_BaseBlock.getLevel() == std::uint8_t(0)) {
    _overlap = 1;
  }
  BasicBlock<T, 2>::operator=(_BaseBlock.getExtBlock(_overlap));
}

#ifdef MPI_ENABLED

template <typename T>
void Block2D<T>::GeoFlagMPIComm() {
  Mpi().barrier();
  // send data
  // send buffer
  // add to send buffer
  for (int i = 0; i < MPIComm.Senders.size(); ++i) {
    std::vector<std::uint8_t> &buffer = MPIBuffer.SendBuffers[i];
    const std::vector<std::size_t> &sendcells = MPIComm.Senders[i].SendCells;
    for (std::size_t id = 0; id < sendcells.size(); ++id) {
      buffer[id] = GeometryFlag.get(sendcells[id]);
    }
  }
  // non-blocking send
  std::vector<MPI_Request> requests;
  for (int i = 0; i < MPIComm.Senders.size(); ++i) {
    std::vector<std::uint8_t> &buffer = MPIBuffer.SendBuffers[i];
    MPI_Request request;
    Mpi().iSend(buffer.data(), buffer.size(), MPIComm.Senders[i].RecvRank, &request, 0,
                MPI_COMM_WORLD);
    // MPI_ISend()
    requests.push_back(request);
  }
  // non-blocking recv
  for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
    std::vector<std::uint8_t> &buffer = MPIBuffer.RecvBuffers[i];
    MPI_Request request;
    Mpi().iRecv(buffer.data(), buffer.size(), MPIComm.Recvers[i].SendRank, &request, 0,
                MPI_COMM_WORLD);
    requests.push_back(request);
  }
  // wait for all requests
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  // set flag based on recv buffer
  for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
    const std::vector<std::uint8_t> &buffer = MPIBuffer.RecvBuffers[i];
    const std::vector<int> &recvcells = MPIComm.Recvers[i].RecvCells;
    for (int j = 0; j < buffer.size(); ++j) {
      GeometryFlag.SetField(recvcells[j], buffer[j]);
    }
  }
  Mpi().barrier();
}

#endif

// -----------blockgeometry2d----------------


template <typename T>
BlockGeometry2D<T>::BlockGeometry2D(int Nx, int Ny, int blocknum, const AABB<T, 2> &block,
                                    T voxelSize, int overlap, bool refine)
    : BasicBlock<T, 2>(voxelSize, block.getExtended(Vector<T, 2>{voxelSize}),
                       AABB<int, 2>(Vector<int, 2>{0}, Vector<int, 2>{Nx + 1, Ny + 1})),
      _BaseBlock(voxelSize, block,
                 AABB<int, 2>(Vector<int, 2>{1}, Vector<int, 2>{Nx, Ny})),
      _overlap(overlap) {
  DivideBlocks(blocknum);
  CreateBlocks();
  SetupNbrs();
  if (!refine) {
    InitCommunicators();
  }
}

template <typename T>
BlockGeometry2D<T>::BlockGeometry2D(BlockGeometryHelper2D<T> &GeoHelper)
    : _BaseBlock(GeoHelper.getBaseBlock()), BasicBlock<T, 2>(GeoHelper),
      _overlap(GeoHelper.getExt()) {
  // create blocks from GeoHelper
  for (BasicBlock<T, 2> &baseblock : GeoHelper.getBasicBlocks()) {
    int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
    _Blocks.emplace_back(baseblock, overlap);
  }
  SetupNbrs();
  InitAllComm();
  UpdateMaxLevel();
}

template <typename T>
void BlockGeometry2D<T>::Init(BlockGeometryHelper2D<T> &GeoHelper) {
  ;
  _Blocks.clear();
  // create blocks from GeoHelper
  for (BasicBlock<T, 2> &baseblock : GeoHelper.getBasicBlocks()) {
    int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
    _Blocks.emplace_back(baseblock, overlap);
  }
  SetupNbrs();
  InitAllComm();
  UpdateMaxLevel();
}

template <typename T>
void BlockGeometry2D<T>::UpdateMaxLevel() {
  _MaxLevel = std::uint8_t(0);
  for (const Block2D<T> &block : _Blocks)
    _MaxLevel = std::max(_MaxLevel, block.getLevel());
}

template <typename T>
std::size_t BlockGeometry2D<T>::getTotalCellNum() const {
  std::size_t sum = 0;
  for (const Block2D<T> &block : _Blocks) {
    sum += block.getN();
  }
  return sum;
}

template <typename T>
std::size_t BlockGeometry2D<T>::getBaseCellNum() const {
  std::size_t sum = 0;
  for (const Block2D<T> &block : _Blocks) {
    sum += block.getBaseBlock().getN();
  }
  return sum;
}

template <typename T>
void BlockGeometry2D<T>::DivideBlocks(int blocknum) {
  _BlockAABBs.clear();
  _BlockAABBs.reserve(blocknum);
  DivideBlock2D(_BaseBlock, blocknum, _BlockAABBs);
}

template <typename T>
void BlockGeometry2D<T>::CreateBlocks() {
  _Blocks.clear();
  _Blocks.reserve(_BlockAABBs.size());
  // create blocks
  int blockid = 0;
  for (const AABB<int, 2> &blockaabb : _BlockAABBs) {
    Vector<T, 2> MIN =
      blockaabb.getMin() * _BaseBlock.getVoxelSize() + BasicBlock<T, 2>::_min;
    Vector<T, 2> MAX =
      (blockaabb.getMax() + Vector<T, 2>{T(1)}) * _BaseBlock.getVoxelSize() +
      BasicBlock<T, 2>::_min;
    AABB<T, 2> aabb(MIN, MAX);
    _Blocks.emplace_back(aabb, blockaabb, blockid, _BaseBlock.getVoxelSize(), _overlap);
    blockid++;
  }
}

template <typename T>
void BlockGeometry2D<T>::SetupNbrs() {
  for (Block2D<T> &block : _Blocks) {
    int id = block.getBlockId();
    std::vector<Block2D<T> *> &nbrsvec = block.getNeighbors();
    nbrsvec.clear();
    for (Block2D<T> &blockn : _Blocks) {
      int idn = blockn.getBlockId();
      if (id != idn) {
        if (isOverlapped(block, blockn.getBaseBlock())) {
          nbrsvec.push_back(&blockn);
        }
      }
    }
  }
}

template <typename T>
void BlockGeometry2D<T>::InitCommunicators() {
  for (int i = 0; i < _Blocks.size(); ++i) {
    // getCellIdx() called from ext block
    Block2D<T> &block = _Blocks[i];
    // get the first layer of overlapped cells(counted from inside to outside)
    BasicBlock<T, 2> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
    std::uint8_t blocklevel = block.getLevel();
    std::vector<BlockComm<T, 2>> &Communicators = block.getCommunicators();
    Communicators.clear();
    for (Block2D<T> *nblock : block.getNeighbors()) {
      // check if 2 blocks are of the same level
      if (nblock->getLevel() == blocklevel) {
        Communicators.emplace_back(nblock);
        BlockComm<T, 2> &comm = Communicators.back();
        // blocks of the same level only communicate with the first layer of overlapped
        // cells
        block.getCellIdx(baseblock_ext1, nblock->getBaseBlock(), comm.RecvCells);
        nblock->getCellIdx(nblock->getBaseBlock(), baseblock_ext1, comm.SendCells);
      }
    }
  }
}

template <typename T>
void BlockGeometry2D<T>::InitAverComm(int highlevelovlap) {
  for (int i = 0; i < _Blocks.size(); ++i) {
    // getCellIdx() called from ext block
    Block2D<T> &block = _Blocks[i];
    // get the first layer of overlapped cells(counted from inside to outside)
    BasicBlock<T, 2> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
    std::uint8_t blocklevel = block.getLevel();
    std::vector<InterpBlockComm<T, 2>> &Communicators = block.getAverageBlockComm();
    Communicators.clear();
    for (Block2D<T> *nblock : block.getNeighbors()) {
      // find block of blocklevel+1
      if (nblock->getLevel() == blocklevel + 1) {
        Communicators.emplace_back(nblock);
        InterpBlockComm<T, 2> &comm = Communicators.back();
        // init recv cells
        // low level block only recv from higher using the first layer of overlapped cells
        block.getCellIdx(baseblock_ext1, nblock->getBaseBlock(), comm.RecvCells);
        // get coarsened basicblock of nblock
        BasicBlock<T, 2> Cnbaseblock = nblock->getBaseBlock().getCoasenedBlock();
        // buffer container to store virtual cell indices
        // TODO: more efficient implementation without using virtual cells?
        std::vector<std::size_t> VSendCells;
        // get virtual (coarse) cell indices on coarsened basicblock of nblock
        Cnbaseblock.getCellIdx(baseblock_ext1, Cnbaseblock, VSendCells);
        // get real cell indices from virtual cell indices
        for (std::size_t id : VSendCells) {
          // virtual (coarse) cell location
          Vector<int, 2> Clocidx = Cnbaseblock.getLoc(id);
          // real cell location
          Vector<int, 2> locidx0 = Clocidx * 2 + Vector<int, 2>{highlevelovlap};
          std::size_t id0 = nblock->getIndex(locidx0);
          Vector<int, 2> locidx1 = locidx0 + Vector<int, 2>{1, 0};
          std::size_t id1 = nblock->getIndex(locidx1);
          Vector<int, 2> locidx2 = locidx0 + Vector<int, 2>{0, 1};
          std::size_t id2 = nblock->getIndex(locidx2);
          Vector<int, 2> locidx3 = locidx0 + Vector<int, 2>{1, 1};
          std::size_t id3 = nblock->getIndex(locidx3);
          comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
        }
      } else if (nblock->getLevel() > blocklevel + 1) {
        std::cerr << "[BlockGeometry2D<T>::InitAverComm] Error: block level difference "
                     "larger than 1"
                  << std::endl;
      }
    }
  }
}

template <typename T>
void BlockGeometry2D<T>::InitIntpComm() {
  for (int i = 0; i < _Blocks.size(); ++i) {
    // getCellIdx() called from ext block
    Block2D<T> &block = _Blocks[i];
    std::uint8_t blocklevel = block.getLevel();
    std::vector<InterpBlockComm<T, 2>> &Communicators = block.getInterpBlockComm();
    Communicators.clear();
    for (Block2D<T> *nblock : block.getNeighbors()) {
      // find block of blocklevel-1
      if (nblock->getLevel() == blocklevel - 1) {
        Communicators.emplace_back(nblock);
        InterpBlockComm<T, 2> &comm = Communicators.back();
        // init recv cells
        // high level block recv from lower using 2^(diff_level) layers of overlapped
        // cells
        block.getCellIdx(block, nblock->getBaseBlock(), comm.RecvCells);
        // get refined basicblock of nblock
        BasicBlock<T, 2> Rnbaseblock = nblock->getBaseBlock().getRefinedBlock();
        // buffer container to store virtual cell indices
        std::vector<std::size_t> VSendCells;
        // get virtual (refined) cell indices on refined basicblock of nblock
        Rnbaseblock.getCellIdx(block, Rnbaseblock, VSendCells);
        // get real cell indices from virtual cell indices
        T halfvoxsize = nblock->getVoxelSize() / T(2);
        for (std::size_t id : VSendCells) {
          // refined cell location
          Vector<T, 2> Rloc = Rnbaseblock.getLoc_t(id);
          // get dist of Rloc and MinCenter of nblock
          Vector<T, 2> dist = Rloc - nblock->getMinCenter();
          // find nearest real cell location
          int x = static_cast<int>(dist[0] / nblock->getVoxelSize());
          int y = static_cast<int>(dist[1] / nblock->getVoxelSize());
          T remx = dist[0] - x * nblock->getVoxelSize();
          T remy = dist[1] - y * nblock->getVoxelSize();
          int scheme = 0;
          if (remx < halfvoxsize) {
            scheme += 1;
          }
          if (remy < halfvoxsize) {
            scheme += 2;
          }
          /*
          2 3
          0 1
          */
          std::size_t id0 = nblock->getIndex(Vector<int, 2>{x, y});
          std::size_t id1 = nblock->getIndex(Vector<int, 2>{x + 1, y});
          std::size_t id2 = nblock->getIndex(Vector<int, 2>{x, y + 1});
          std::size_t id3 = nblock->getIndex(Vector<int, 2>{x + 1, y + 1});
          comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
          switch (scheme) {
            case 0: {
              comm.InterpWeights.emplace_back(
                InterpWeight<T, 2>{T(0.0625), T(0.1875), T(0.1875), T(0.5625)});
            } break;
            case 1: {
              comm.InterpWeights.emplace_back(
                InterpWeight<T, 2>{T(0.1875), T(0.0625), T(0.5625), T(0.1875)});
            } break;
            case 2: {
              comm.InterpWeights.emplace_back(
                InterpWeight<T, 2>{T(0.1875), T(0.5625), T(0.0625), T(0.1875)});
            } break;
            case 3: {
              comm.InterpWeights.emplace_back(
                InterpWeight<T, 2>{T(0.5625), T(0.1875), T(0.1875), T(0.0625)});
            } break;
          }
        }
      } else if (nblock->getLevel() < blocklevel - 1) {
        std::cerr << "[BlockGeometry2D<T>::InitIntpComm] Error: block level difference "
                     "larger than 1"
                  << std::endl;
      }
    }
  }
}

template <typename T>
void BlockGeometry2D<T>::InitAllComm() {
  InitCommunicators();
  InitAverComm();
  InitIntpComm();
}

template <typename T>
void BlockGeometry2D<T>::GeoFlagComm() {
  // comm of same level
  for (Block2D<T> &block : _Blocks) {
    FlagField &GeoFlags = block.getGeoFlagField();
    for (BlockComm<T, 2> &comm : block.getCommunicators()) {
      Block2D<T> *nblock = comm.SendBlock;
      FlagField &nGeoFlags = nblock->getGeoFlagField();
      int size = comm.RecvCells.size();
      for (int i = 0; i < size; ++i) {
        std::size_t idrecv = comm.RecvCells[i];
        std::size_t idsend = comm.SendCells[i];
        GeoFlags.SetField(idrecv, nGeoFlags.get(idsend));
      }
    }
  }
}

// BlockGeometryHelper2D

template <typename T>
BlockGeometryHelper2D<T>::BlockGeometryHelper2D(int Nx, int Ny, int blocklen,
                                                const AABB<T, 2> &AABBs, T voxelSize,
                                                std::uint8_t llimit, int ext)
    : BasicBlock<T, 2>(voxelSize, AABBs.getExtended(Vector<T, 2>{voxelSize * ext}),
                       AABB<int, 2>(Vector<int, 2>{0}, Vector<int, 2>{Nx + 1, Ny + 1})),
      _BaseBlock(voxelSize, AABBs,
                 AABB<int, 2>(Vector<int, 2>{1}, Vector<int, 2>{Nx, Ny})),
      BlockLen(blocklen), Ext(ext), _MaxLevel(std::uint8_t(0)), _Exchanged(true),
      _LevelLimit(llimit) {
  if (BlockLen < 4) {
    std::cerr << "BlockGeometryHelper2D<T>, BlockLen < 4" << std::endl;
  }
  CellsNx = _BaseBlock.getNx() / BlockLen;
  CellsNy = _BaseBlock.getNy() / BlockLen;
  CellsN = CellsNx * CellsNy;

  Delta_Cellidx = {-CellsNx - 1, -CellsNx, -CellsNx + 1, -1, 1,
                   CellsNx - 1,  CellsNx,  CellsNx + 1};

  CreateBlockCells();
}

template <typename T>
void BlockGeometryHelper2D<T>::CreateBlockCells() {
  // buffer vector to store cell aabbs
  std::vector<AABB<int, 2>> AABBCells;
  AABBCells.reserve(CellsN);
  // divide base block into cell aabbs
  _BaseBlock.getIdxBlock().divide(CellsNx, CellsNy, AABBCells);
  // create cell blocks from cell aabbs
  _BlockCells.reserve(CellsN);
  int blockid = 0;
  T voxsize = BasicBlock<T, 2>::getVoxelSize();
  Vector<T, 2> _Min = BasicBlock<T, 2>::_min;
  for (const AABB<int, 2> &aabbcell : AABBCells) {
    Vector<T, 2> MIN = aabbcell.getMin() * voxsize + _Min;
    Vector<T, 2> MAX = (aabbcell.getMax() + Vector<T, 2>{T(1)}) * voxsize + _Min;
    AABB<T, 2> aabb(MIN, MAX);
    _BlockCells.emplace_back(voxsize, aabb, aabbcell, blockid);
    _BlockCellTags.emplace_back(BlockCellTag::none);
    ++blockid;
  }
}

template <typename T>
void BlockGeometryHelper2D<T>::UpdateMaxLevel() {
  _MaxLevel = std::uint8_t(0);
  for (const BasicBlock<T, 2> &block : _BlockCells) {
    if (block.getLevel() > _MaxLevel) _MaxLevel = block.getLevel();
  }
}

template <typename T>
void BlockGeometryHelper2D<T>::CreateBlocks() {
  // create new blocks on relatively older blocks
  std::vector<BasicBlock<T, 2>> &BasicBlocks = getOldBasicBlocks();
  BasicBlocks.clear();
  // now old blocks become new blocks
  _Exchanged = !_Exchanged;
  // get max refine level
  // int maxlevel = 0;
  // for (const BasicBlock<T, 2>& block : _BlockCells) {
  //   maxlevel = std::max(maxlevel, static_cast<int>(block.getLevel()));
  // }
  // buffer array visited
  bool *visited = new bool[CellsN];
  std::fill(visited, visited + CellsN, false);
  // create blocks
  int blockid = 0;
  for (int j = 0; j < CellsNy; ++j) {
    for (int i = 0; i < CellsNx; ++i) {
      std::size_t id = i + j * CellsNx;
      if (visited[id]) {
        continue;
      }
      // get level
      std::uint8_t level = _BlockCells[id].getLevel();
      T voxsize = _BlockCells[id].getVoxelSize();
      Vector<int, 2> NewMesh = _BlockCells[id].getMesh();
      // expand block along x
      int Nx = 1;
      std::size_t tempid = id + Nx;
      while (i + Nx < CellsNx && _BlockCells[tempid].getLevel() == level &&
             !visited[tempid]) {
        NewMesh[0] += _BlockCells[tempid].getNx();
        ++Nx;
        ++tempid;
      }
      // expand block along y
      int Ny = 1;
      std::size_t startid = id + Ny * CellsNx;
      while (j + Ny < CellsNy) {
        startid = id + Ny * CellsNx;
        tempid = startid;
        for (int k = 0; k < Nx; ++k) {
          if (_BlockCells[tempid].getLevel() == level && !visited[tempid]) {
            ++tempid;
          } else {
            goto end_y_expansion;
          }
        }
        ++Ny;
        NewMesh[1] += _BlockCells[startid].getNy();
      }
    end_y_expansion:
      // create block
      Vector<int, 2> Ext = BlockLen * Vector<int, 2>{Nx, Ny};
      Vector<int, 2> min = _BlockCells[id].getIdxBlock().getMin();
      Vector<int, 2> max = min + Ext - Vector<int, 2>{1};
      AABB<int, 2> idxblock(min, max);
      Vector<T, 2> MIN = _BlockCells[id].getMin();
      int ratio = static_cast<int>(std::pow(2, static_cast<int>(level)));
      Vector<T, 2> MAX = Ext * voxsize * ratio + MIN;
      AABB<T, 2> block(MIN, MAX);
      BasicBlocks.emplace_back(level, voxsize, blockid, block, idxblock, NewMesh);
      blockid++;
      // set visited
      for (int jj = 0; jj < Ny; ++jj) {
        startid = id + jj * CellsNx;
        for (int ii = 0; ii < Nx; ++ii) {
          visited[startid + ii] = true;
        }
      }
    }
  }
  // delete buffer
  delete[] visited;
  // update max level
  UpdateMaxLevel();
}

template <typename T>
template <typename Func>
void BlockGeometryHelper2D<T>::forEachBlockCell(Func func) {
  for (BasicBlock<T, 2> &block : _BlockCells) {
    func(block);
  }
}

template <typename T>
void BlockGeometryHelper2D<T>::AdaptiveOptimization(int OptProcNum, int MaxProcNum,
                                                    bool enforce) {
  if (MaxProcNum == -1) {
    MaxProcNum = 2 * OptProcNum;
  }
  // get new basic blocks
  std::vector<BasicBlock<T, 2>> &BasicBlocks = getBasicBlocks();
  // vector store stddev of each optimization scheme
  std::vector<T> StdDevs;
  int minProcNum = std::min(OptProcNum, static_cast<int>(BasicBlocks.size()));
  int maxProcNum = std::max(MaxProcNum, static_cast<int>(BasicBlocks.size()));
  for (int i = minProcNum; i <= maxProcNum; ++i) {
    // buffer vector copied from BasicBlocks
    std::vector<BasicBlock<T, 2>> Blocks = BasicBlocks;
    Optimize(Blocks, i, enforce);
    StdDevs.push_back(ComputeStdDev(Blocks));
  }
  // find shcemes with minimum stddev
  std::vector<int> bestSchemesvec;
  T minStdDev = StdDevs[0];
  bestSchemesvec.push_back(minProcNum);
  for (int i = 1; i < StdDevs.size(); ++i) {
    if (StdDevs[i] < minStdDev) {
      minStdDev = StdDevs[i];
      bestSchemesvec.clear();
      bestSchemesvec.push_back(i + minProcNum);
    } else if (StdDevs[i] == minStdDev) {
      bestSchemesvec.push_back(i + minProcNum);
    }
  }
  // find the best shceme, which is closest to OptProcNum
  int bestScheme = bestSchemesvec[0];
  int minDiff = std::abs(bestScheme - OptProcNum);
  for (int scheme : bestSchemesvec) {
    int diff = std::abs(scheme - OptProcNum);
    if (diff < minDiff) {
      minDiff = diff;
      bestScheme = scheme;
    }
  }
  // apply the best scheme
  Optimize(bestScheme, enforce);
  std::cout << "Optimization result: " << BasicBlocks.size()
            << " Blocks with stdDev: " << minStdDev << std::endl;
}

template <typename T>
void BlockGeometryHelper2D<T>::Optimize(int ProcessNum, bool enforce) {
  Optimize(getBasicBlocks(), ProcessNum, enforce);
}

template <typename T>
void BlockGeometryHelper2D<T>::Optimize(std::vector<BasicBlock<T, 2>> &Blocks,
                                        int ProcessNum, bool enforce) {
  // get total number of points
  std::size_t Total = 0;
  for (const BasicBlock<T, 2> &block : Blocks) {
    Total += block.getN();
  }
  // get number of points per process
  std::size_t NumPerProcess = Total / ProcessNum;

  // divide large blocks
  // T threshold = static_cast<T>(ProcessNum) / size;
  // buffer vector to store new blocks
  std::vector<BasicBlock<T, 2>> NewBasicBlocks;
  // iterate through all blocks
  typename std::vector<BasicBlock<T, 2>>::iterator it = Blocks.begin();

  if (enforce && Blocks.size() < ProcessNum) {
    // sort blocks by number of points, large blocks first
    std::sort(Blocks.begin(), Blocks.end(),
              [](const BasicBlock<T, 2> &a, const BasicBlock<T, 2> &b) {
                return a.getN() > b.getN();
              });
    // make size = ProcessNum
    int count = Blocks.size();
    while (count < ProcessNum && it != Blocks.end()) {
      T ratio = std::round(static_cast<T>(it->getN()) / NumPerProcess);
      int part = static_cast<int>(ratio);
      if (ratio < 2) {
        part = 2;
      }
      // divide block
      DivideBlock2D(*it, part, NewBasicBlocks);
      // remove block
      it = Blocks.erase(it);
      count += part - 1;
    }
  } else {
    while (it != Blocks.end()) {
      T ratio = std::round(static_cast<T>(it->getN()) / NumPerProcess);
      if (ratio >= 2) {
        // divide block
        DivideBlock2D(*it, static_cast<int>(ratio), NewBasicBlocks);
        // remove block
        it = Blocks.erase(it);
      } else {
        ++it;
      }
    }
  }
  // insert new blocks
  Blocks.insert(Blocks.end(), NewBasicBlocks.begin(), NewBasicBlocks.end());
  // update block id
  for (int i = 0; i < Blocks.size(); ++i) {
    Blocks[i].setBlockId(i);
  }
}

template <typename T>
T BlockGeometryHelper2D<T>::ComputeStdDev() const {
  return ComputeStdDev(getBasicBlocks());
}

template <typename T>
T BlockGeometryHelper2D<T>::ComputeStdDev(
  const std::vector<BasicBlock<T, 2>> &Blocks) const {
  T mean = 0;
  for (const BasicBlock<T, 2> &block : Blocks) {
    mean += block.getN();
  }
  mean /= Blocks.size();
  T stdDev = 0;
  for (const BasicBlock<T, 2> &block : Blocks) {
    stdDev += std::pow((static_cast<T>(block.getN()) / mean - T(1)), 2);
  }
  stdDev = std::sqrt(stdDev / Blocks.size());
  return stdDev;
}

template <typename T>
void BlockGeometryHelper2D<T>::TagRefineLayer(std::vector<bool> &refine, bool &refined) {
  UpdateMaxLevel();
  // refine one additional layer if has neighbor with lower level
  // tag cells to be refined
  for (std::uint8_t level = _MaxLevel; level > std::uint8_t(0); --level) {
    for (int celly = 1; celly < CellsNy - 1; ++celly) {
      for (int cellx = 1; cellx < CellsNx - 1; ++cellx) {
        int cellid = celly * CellsNx + cellx;
        if (refine[cellid]) {
          if (_BlockCells[cellid].getLevel() == level) {
            for (int delta : Delta_Cellidx) {
              int ncellid = cellid + delta;
              if (_BlockCells[ncellid].getLevel() < level &&
                  util::isFlag(_BlockCellTags[ncellid], BlockCellTag::none)) {
                _BlockCellTags[ncellid] = BlockCellTag::refine;
                refined = true;
              }
            }
          }
        }
      }
    }
  }
}

template <typename T>
void BlockGeometryHelper2D<T>::CheckRefine() {
  Refine();
  UpdateMaxLevel();
  // refine one additional layer if has neighbor with 2 lower level
  // tag cells to be refined
  for (std::uint8_t level = _MaxLevel; level > std::uint8_t(1); --level) {
    for (int celly = 1; celly < CellsNy - 1; ++celly) {
      for (int cellx = 1; cellx < CellsNx - 1; ++cellx) {
        int cellid = celly * CellsNx + cellx;
        if (_BlockCells[cellid].getLevel() == level) {
          for (int delta : Delta_Cellidx) {
            int ncellid = cellid + delta;
            if (_BlockCells[ncellid].getLevel() < (level - std::uint8_t(1)) &&
                util::isFlag(_BlockCellTags[ncellid], BlockCellTag::none)) {
              _BlockCellTags[ncellid] = BlockCellTag::refine;
            }
          }
        }
      }
    }
    Refine();
  }
}

template <typename T>
void BlockGeometryHelper2D<T>::Refine() {
  for (int i = 0; i < CellsN; ++i) {
    if (util::isFlag(_BlockCellTags[i], BlockCellTag::refine)) {
      _BlockCells[i].refine();
      // reset tag
      _BlockCellTags[i] = BlockCellTag::none;
    }
  }
}

// helper class for MPI

template <typename T>
BlockGeometryMPIHelper2D<T>::BlockGeometryMPIHelper2D(int Nx, int Ny, int blocknum,
                                                      const AABB<T, 2> &block,
                                                      T voxelSize, int overlap,
                                                      bool refine)
    : BasicBlock<T, 2>(voxelSize, block.getExtended(Vector<T, 2>{voxelSize}),
                       AABB<int, 2>(Vector<int, 2>{0}, Vector<int, 2>{Nx + 1, Ny + 1})),
      _BaseBlock(voxelSize, block,
                 AABB<int, 2>(Vector<int, 2>{1}, Vector<int, 2>{Nx, Ny})),
      _BlockNum(blocknum) {
  DivideBlocks(_BlockNum);
  CreateBlocks();
  SetupNbrs();
}

template <typename T>
void BlockGeometryMPIHelper2D<T>::DivideBlocks(int blocknum) {
  _BlockAABBs.clear();
  _BlockAABBs.reserve(blocknum);
  DivideBlock2D(_BaseBlock, blocknum, _BlockAABBs);
}

template <typename T>
void BlockGeometryMPIHelper2D<T>::CreateBlocks() {
  // get std::vector<BasicBlock<T, 2>>
  _Blocks.clear();
  _Blocks.reserve(_BlockNum);
  int blockid = 0;
  for (const AABB<int, 2> &block : _BlockAABBs) {
    Vector<T, 2> MIN =
      block.getMin() * _BaseBlock.getVoxelSize() + BasicBlock<T, 2>::_min;
    Vector<T, 2> MAX = (block.getMax() + Vector<T, 2>{T(1)}) * _BaseBlock.getVoxelSize() +
                       BasicBlock<T, 2>::_min;
    AABB<T, 2> aabb(MIN, MAX);
    _Blocks.emplace_back(_BaseBlock.getVoxelSize(), aabb, block, blockid);
    ++blockid;
  }
}

template <typename T>
void BlockGeometryMPIHelper2D<T>::SetupNbrs() {
  _BlockNbrRanks.resize(_BlockNum, std::vector<int>{});
  for (const BasicBlock<T, 2> &block : _Blocks) {
    int selfrank = block.getBlockId();
    AABB<int, 2> ExtIdxblock = block.getIdxBlock().getExtended(Vector<int, 2>{1});
    for (const BasicBlock<T, 2> &blockn : _Blocks) {
      int rankn = blockn.getBlockId();
      if (rankn != selfrank) {
        if (isOverlapped(ExtIdxblock, blockn.getIdxBlock())) {
          _BlockNbrRanks[selfrank].push_back(rankn);
        }
      }
    }
  }
}

template <typename T>
void BlockGeometryMPIHelper2D<T>::InitMPIBlockCommStru(MPIBlockCommStru &BlockComm) {
  int rank = Mpi().getRank();
  // senders
  // get basicblock to call getCellIdx
  const BasicBlock<T, 2> ExtBlock = _Blocks[rank].getExtBlock();
  // ext block to init recvers
  const AABB<T, 2> ExtAABB = ExtBlock.getAABB();
  // base block to init senders
  const AABB<T, 2> BaseAABB = _Blocks[rank].getAABB();
  // const AABB<int, 2> Idxblock = _Blocks[rank].getIdxBlock();
  // get neighbor ranks
  std::vector<int> &nbrRanks = _BlockNbrRanks[rank];
  for (int nbrRank : nbrRanks) {
    BlockComm.Senders.emplace_back(nbrRank, nbrRank);
    MPIBlockSendStru &send = BlockComm.Senders.back();
    const AABB<T, 2> ExtAABBn = _Blocks[nbrRank].getExtBlock().getAABB();
    ExtBlock.getCellIdx(BaseAABB, ExtAABBn, send.SendCells);
  }
  // recvers
  // get basicblock
  const AABB<int, 2> ExtIdxblock =
    _Blocks[rank].getIdxBlock().getExtended(Vector<int, 2>{1});
  for (int nbrRank : nbrRanks) {
    BlockComm.Recvers.emplace_back(nbrRank, nbrRank);
    MPIBlockRecvStru &recv = BlockComm.Recvers.back();
    const AABB<T, 2> BaseAABBn = _Blocks[nbrRank].getAABB();
    ExtBlock.getCellIdx(BaseAABBn, ExtAABB, recv.RecvCells);
  }
}

// template <typename T>
// void BlockGeometry2D<T>::InitAverComm2() {
//   for (int i = 0; i < _Blocks.size(); ++i) {
//     // getCellIdx() called from ext block
//     Block2D<T> &block = _Blocks[i];
//     BasicBlock<T, 2> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
//     std::uint8_t blocklevel = block.getLevel();
//     std::vector<InterpBlockComm<T, 2>> &Communicators = block.getAverageBlockComm();
//     for (Block2D<T> *nblock : block.getNeighbors()) {
//       // find block of blocklevel+1
//       if (nblock->getLevel() == blocklevel + 1) {
//         Communicators.emplace_back(nblock);
//         InterpBlockComm<T, 2> &comm = Communicators.back();
//         // init recv cells
//         // low level block only recv from higher using the first layer of overlapped
//         cells block.getCellIdx(baseblock_ext1, nblock->getBaseBlock(), comm.RecvCells);
//         // buffer container to store cell indices
//         std::vector<int> Sends;
//         // get virtual (coarse) cell indices on coarsened basicblock of nblock
//         nblock->getCellIdx(baseblock_ext1, nblock->getBaseBlock(), Sends);
//         // check
//         if (Sends.size() != 4 * comm.RecvCells.size()) {
//           std::cerr << "Error: mismatch of SendCells and RecvCells" << std::endl;
//         }
//         // get scheme
//         Vector<int, 2> loc0 = nblock->getLoc(Sends[0]);
//         Vector<int, 2> loc2 = nblock->getLoc(Sends[2]);
//         Vector<int, 2> Align_dir = loc2 - loc0;
//         if (Align_dir[0] == 0) {
//           // vertical
//           for (int id = 0; id < Sends.size(); id += 4) {
//             comm.SendCells.emplace_back(
//               InterpSource<2>{Sends[id], Sends[id + 1], Sends[id + 2], Sends[id + 3]});
//           }
//         } else if (Align_dir[1] == 0) {
//           // horizontal
//           int halfsize = Sends.size() / 2;
//           for (int id = 0; id < halfsize; id += 2) {
//             comm.SendCells.emplace_back(InterpSource<2>{
//               Sends[id], Sends[id + 1], Sends[id + halfsize], Sends[id + halfsize +
//               1]});
//           }
//         } else {
//           std::cerr << "Error: non-aligned cells" << std::endl;
//         }
//       } else if (nblock->getLevel() > blocklevel + 1) {
//         std::cerr << "Error: block level difference larger than 1" << std::endl;
//       }
//     }
//   }
// }

// template <typename T>
// void BlockGeometry2D<T>::InitExtpComm() {
//   for (int i = 0; i < _Blocks.size(); ++i) {
//     // getCellIdx() called from ext block
//     Block2D<T> &block = _Blocks[i];
//     std::uint8_t blocklevel = block.getLevel();
//     std::vector<InterpBlockComm<T, 2>> &Communicators = block.getInterpBlockComm();
//     for (Block2D<T> *nblock : block.getNeighbors()) {
//       // find block of blocklevel-1
//       if (nblock->getLevel() == blocklevel - 1) {
//         Communicators.emplace_back(nblock);
//         InterpBlockComm<T, 2> &comm = Communicators.back();
//         // init recv cells
//         block.getCellIdx(block.getAABB(), nblock->getBaseBlock().getAABB(),
//         comm.RecvCells);
//         // get refined basicblock of nblock
//         BasicBlock<T, 2> Rnbaseblock = nblock->getBaseBlock().getRefinedBlock();
//         // buffer container to store virtual cell indices
//         std::vector<int> VSendCells;
//         // get virtual (refined) cell indices on refined basicblock of nblock
//         Rnbaseblock.getCellIdx(block.getAABB(), Rnbaseblock.getAABB(), VSendCells);
//         // get real cell indices from virtual cell indices
//         /*
//              2
//             ---
//            3\ \1
//             ---
//              0
//         */
//         // get virtual cells' alignment direction
//         Vector<int, 2> Loc_id0 = Rnbaseblock.getLoc(VSendCells[0]);
//         Vector<int, 2> Loc_id1 = Rnbaseblock.getLoc(VSendCells[1]);
//         Vector<int, 2> Align_dir = Loc_id1 - Loc_id0;
//         // real cell location of cell 0
//         Vector<T, 2> Rloc = Rnbaseblock.getLoc_t(VSendCells[0]);
//         // get dist of Rloc and MinCenter of nblock
//         Vector<T, 2> dist = Rloc - nblock->getMinCenter();
//         // find nearest real cell location
//         int x = static_cast<int>(dist[0] / nblock->getVoxelSize());
//         int y = static_cast<int>(dist[1] / nblock->getVoxelSize());
//         T remx = dist[0] - x * nblock->getVoxelSize();
//         T remy = dist[1] - y * nblock->getVoxelSize();
//         T halfvoxsize = nblock->getVoxelSize() / T(2);
//         int scheme = 0;
//         if (Align_dir[0] == 0) {
//           // scheme 1 or 3
//           if (remx > halfvoxsize) {
//             scheme = 3;
//           } else if (remx < halfvoxsize) {
//             scheme = 1;
//           }
//         } else if (Align_dir[1] == 0) {
//           // scheme 0 or 2
//           if (remy > halfvoxsize) {
//             scheme = 0;
//           } else if (remy < halfvoxsize) {
//             scheme = 2;
//           }
//         } else {
//           std::cerr << "Error: InitExtpComm, Align_dir error" << std::endl;
//         }
//         switch (scheme) {
//           case 0: {
//             // first cell
//             Vector<T, 2> Rloc = Rnbaseblock.getLoc_t(VSendCells[0]);
//             Vector<T, 2> dist = Rloc - nblock->getMinCenter();
//             int x = static_cast<int>(dist[0] / nblock->getVoxelSize());
//             int y = static_cast<int>(dist[1] / nblock->getVoxelSize());
//             int basex = x + 1;
//             int basey = y + 1;
//             int id0 = nblock->getIndex(Vector<int, 2>{basex, basey});
//             int id1 = nblock->getIndex(Vector<int, 2>{basex + 1, basey});
//             int id2 = nblock->getIndex(Vector<int, 2>{basex, basey + 1});
//             int id3 = nblock->getIndex(Vector<int, 2>{basex + 1, basey + 1});
//             comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(1.25), T(0), T(0),
//             T(-0.25)}); for (int j = 0; j < VSendCells.size() - 1; j += 2) {
//               comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(0.9375), T(0.3125), T(-0.1875), T(-0.0625)});
//               comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(0.3125), T(0.9375), T(-0.0625), T(-0.1875)});
//               // update ids
//               id0 = id1;
//               ++id1;
//               id2 = id3;
//               ++id3;
//             }
//             // last cell
//             --id2;
//             comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(1.25), T(0), T(-0.25),
//             T(0)});
//           } break;
//           case 1: {
//             // first cell
//             Vector<T, 2> Rloc = Rnbaseblock.getLoc_t(VSendCells[0]);
//             Vector<T, 2> dist = Rloc - nblock->getMinCenter();
//             int x = static_cast<int>(dist[0] / nblock->getVoxelSize());
//             int y = static_cast<int>(dist[1] / nblock->getVoxelSize());
//             int basex = x - 1;
//             int basey = y + 1;
//             int id0 = nblock->getIndex(Vector<int, 2>{basex, basey});
//             int id1 = nblock->getIndex(Vector<int, 2>{basex + 1, basey});
//             int id2 = nblock->getIndex(Vector<int, 2>{basex, basey + 1});
//             int id3 = nblock->getIndex(Vector<int, 2>{basex + 1, basey + 1});
//             comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(0), T(1.25), T(-0.25),
//             T(0)}); for (int j = 0; j < VSendCells.size() - 1; j += 2) {
//               comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(-0.1875), T(0.9375), T(-0.0625), T(0.3125)});
//               comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(-0.0625), T(0.3125), T(-0.1875), T(0.9375)});
//               // update ids
//               id0 = id2;
//               id1 = id3;
//               id2 += nblock->getNx();
//               id3 += nblock->getNx();
//             }
//             // last cell
//             id0 -= nblock->getNx();
//             comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(-0.25), T(1.25), T(0),
//             T(0)});
//           } break;
//           case 2: {
//             // first cell
//             Vector<T, 2> Rloc = Rnbaseblock.getLoc_t(VSendCells[0]);
//             Vector<T, 2> dist = Rloc - nblock->getMinCenter();
//             int x = static_cast<int>(dist[0] / nblock->getVoxelSize());
//             int y = static_cast<int>(dist[1] / nblock->getVoxelSize());
//             int basex = x + 1;
//             int basey = y - 1;
//             int id0 = nblock->getIndex(Vector<int, 2>{basex, basey});
//             int id1 = nblock->getIndex(Vector<int, 2>{basex + 1, basey});
//             int id2 = nblock->getIndex(Vector<int, 2>{basex, basey + 1});
//             int id3 = nblock->getIndex(Vector<int, 2>{basex + 1, basey + 1});
//             comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(0), T(-0.25), T(1.25),
//             T(0)}); for (int j = 0; j < VSendCells.size() - 1; j += 2) {
//               comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(-0.1875), T(-0.0625), T(0.9375), T(0.3125)});
//               comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(-0.0625), T(-0.1875), T(0.3125), T(0.9375)});
//               // update ids
//               id0 = id1;
//               ++id1;
//               id2 = id3;
//               ++id3;
//             }
//             // last cell
//             --id0;
//             comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(-0.25), T(0), T(1.25),
//             T(0)});
//           } break;
//           case 3: {
//             // first cell
//             Vector<T, 2> Rloc = Rnbaseblock.getLoc_t(VSendCells[0]);
//             Vector<T, 2> dist = Rloc - nblock->getMinCenter();
//             int x = static_cast<int>(dist[0] / nblock->getVoxelSize());
//             int y = static_cast<int>(dist[1] / nblock->getVoxelSize());
//             int basex = x + 1;
//             int basey = y + 1;
//             int id0 = nblock->getIndex(Vector<int, 2>{basex, basey});
//             int id1 = nblock->getIndex(Vector<int, 2>{basex + 1, basey});
//             int id2 = nblock->getIndex(Vector<int, 2>{basex, basey + 1});
//             int id3 = nblock->getIndex(Vector<int, 2>{basex + 1, basey + 1});
//             comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(1.25), T(0), T(-0.25),
//             T(0)}); for (int j = 0; j < VSendCells.size() - 1; j += 2) {
//               comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(0.9375), T(-0.1875), T(0.3125), T(-0.0625)});
//               comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(0.3125), T(-0.0625), T(0.9375), T(-0.1875)});
//               // update ids
//               id0 = id2;
//               id1 = id3;
//               id2 += nblock->getNx();
//               id3 += nblock->getNx();
//             }
//             // last cell
//             id1 -= nblock->getNx();
//             id2 -= nblock->getNx();
//             comm.SendCells.emplace_back(InterpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(0), T(-0.25), T(1.25),
//             T(0)});
//           }
//         }
//       }
//     }
//   }
// }