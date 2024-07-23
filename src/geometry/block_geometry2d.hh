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

#include "geometry/block_geometry2d.h"

template <typename T>
Block2D<T>::Block2D(const BasicBlock<T, 2> &baseblock, int olap)
    : BasicBlock<T, 2>(baseblock.getExtBlock(olap)), 
      _BaseBlock(baseblock), _overlap(olap){}


template <typename T>
Block2D<T>::Block2D(const AABB<T, 2> &block, const AABB<int, 2> &idxblock, int blockid,
                    T voxelSize, int olap)
    : BasicBlock<T, 2>(voxelSize, block.getExtended(Vector<T, 2>{voxelSize}),
                       idxblock.getExtended(Vector<int, 2>{1}), blockid),
      _BaseBlock(voxelSize, block, idxblock, blockid), _overlap(olap) {
  // read from AABBs
  // ReadAABBs(AABBs, AABBflag);
}

template <typename T>
template <typename FieldType, typename LatSet>
void Block2D<T>::SetupBoundary(const AABB<T, 2> &block, FieldType &field,
                               typename FieldType::value_type bdvalue) {
  // temp flag field store the transition flag
  GenericArray<bool> TransFlag(BasicBlock<T, 2>::N, false);

  for (int y = _overlap; y < BasicBlock<T, 2>::Mesh[1] - _overlap; ++y) {
    for (int x = _overlap; x < BasicBlock<T, 2>::Mesh[0] - _overlap; ++x) {
      const Vector<int, 2> locidx{x, y};
      const Vector<T, 2> vox = BasicBlock<T, 2>::getVoxel(locidx);
      for (unsigned int i = 1; i < LatSet::q; ++i) {
        const Vector<T, 2> nvox = vox + latset::c<LatSet>(i) * BasicBlock<T, 2>::VoxelSize;
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
}

// -----------blockgeometry2d----------------


template <typename T>
BlockGeometry2D<T>::BlockGeometry2D(int Nx, int Ny, int blocknum, const AABB<T, 2> &block,
                                    T voxelSize, int overlap)
    : BasicBlock<T, 2>(voxelSize, block.getExtended(Vector<T, 2>{voxelSize}),
                       AABB<int, 2>(Vector<int, 2>{0}, Vector<int, 2>{Nx + 1, Ny + 1})),
      _BaseBlock(voxelSize, block,
                 AABB<int, 2>(Vector<int, 2>{1}, Vector<int, 2>{Nx, Ny})),
      _overlap(overlap), _MaxLevel(std::uint8_t(0)) {
  CreateBlocks(blocknum);
  SetupNbrs();
  InitComm();
#ifndef MPI_ENABLED
  PrintInfo();
#endif
}

template <typename T>
BlockGeometry2D<T>::BlockGeometry2D(BlockGeometryHelper2D<T> &GeoHelper)
    : BasicBlock<T, 2>(GeoHelper), _BaseBlock(GeoHelper.getBaseBlock()), 
      _overlap(GeoHelper.getExt()), _MaxLevel(GeoHelper.getMaxLevel()) {
  // create blocks from GeoHelper
  for (BasicBlock<T, 2> *baseblock : GeoHelper.getBasicBlocks()) {
    int overlap = (baseblock->getLevel() != std::uint8_t(0)) ? 2 : 1;
    _Blocks.emplace_back(*baseblock, overlap);
  }
  SetupNbrs();
  InitAllComm();
#ifdef MPI_ENABLED
  InitAllMPIComm(GeoHelper);
#else
  PrintInfo();
#endif
}

template <typename T>
void BlockGeometry2D<T>::PrintInfo() const {
  std::cout << "[BlockGeometry2D]: "
            << "Total Cell Num: " << getTotalCellNum() << std::endl;
}

template <typename T>
void BlockGeometry2D<T>::Init(BlockGeometryHelper2D<T> &GeoHelper) {
  ;
  _Blocks.clear();
  // create blocks from GeoHelper
  for (BasicBlock<T, 2> *baseblock : GeoHelper.getBasicBlocks()) {
    int overlap = (baseblock->getLevel() != std::uint8_t(0)) ? 2 : 1;
    _Blocks.emplace_back(*baseblock, overlap);
  }
  SetupNbrs();
  InitAllComm();
#ifdef MPI_ENABLED
  InitAllMPIComm(GeoHelper);
#else
  PrintInfo();
#endif
}

template <typename T>
std::size_t BlockGeometry2D<T>::getTotalCellNum() const {
  std::size_t sum = 0;
  for (const Block2D<T> &block : _Blocks) sum += block.getN();
  return sum;
}

template <typename T>
std::size_t BlockGeometry2D<T>::getBaseCellNum() const {
  std::size_t sum = 0;
  for (const Block2D<T> &block : _Blocks) sum += block.getBaseBlock().getN();
  return sum;
}

template <typename T>
void BlockGeometry2D<T>::DivideBlocks(int blocknum) {
  _BlockAABBs.clear();
  _BlockAABBs.reserve(blocknum);
  DivideBlock2D(_BaseBlock, blocknum, _BlockAABBs);
}

template <typename T>
void BlockGeometry2D<T>::CreateBlocks(int blocknum) {
  DivideBlocks(blocknum);

  _Blocks.clear();
  _Blocks.reserve(_BlockAABBs.size());
  _BasicBlocks.clear();
  _BasicBlocks.reserve(_BlockAABBs.size());
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
    _BasicBlocks.emplace_back(_BaseBlock.getVoxelSize(), aabb, blockaabb, blockid);
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
void BlockGeometry2D<T>::InitComm() {
  for (Block2D<T> &block : _Blocks) {
    std::vector<BlockComm<T, 2>> &Communicators = block.getCommunicators();
    Communicators.clear();
    // get the first layer of overlapped cells(counted from inside to outside)
    BasicBlock<T, 2> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
    std::uint8_t blocklevel = block.getLevel();
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
void BlockGeometry2D<T>::InitAverComm() {
  for (Block2D<T> &block : _Blocks) {
    std::vector<IntpBlockComm<T, 2>> &Communicators = block.getAverageBlockComm();
    Communicators.clear();
    // get the first layer of overlapped cells(counted from inside to outside)
    BasicBlock<T, 2> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
    std::uint8_t blocklevel = block.getLevel();
    for (Block2D<T> *nblock : block.getNeighbors()) {
      // find block of blocklevel+1
      if (nblock->getLevel() == blocklevel + 1) {
        Communicators.emplace_back(nblock);
        IntpBlockComm<T, 2> &comm = Communicators.back();
        // vox size
        const T Cvoxsize = block.getVoxelSize();
        const T Fvoxsize = nblock->getVoxelSize();
        // get intersection
        const AABB<T, 2> intsec = getIntersection(baseblock_ext1, nblock->getBaseBlock());
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        // start index of intsec in block
        Vector<T, 2> startC = intsec.getMin() - block.getMin();
        int startCx = static_cast<int>(std::round(startC[0] / Cvoxsize));
        int startCy = static_cast<int>(std::round(startC[1] / Cvoxsize));
        // start index of intsec in nblock
        Vector<T, 2> startF = intsec.getMin() - nblock->getMin();
        int startFx = static_cast<int>(std::round(startF[0] / Fvoxsize));
        int startFy = static_cast<int>(std::round(startF[1] / Fvoxsize));

        for (int iy = 0; iy < CNy; ++iy) {
          for (int ix = 0; ix < CNx; ++ix) {
            std::size_t Cid = (iy + startCy) * block.getNx() + ix + startCx;
            std::size_t Fid0 = (iy * 2 + startFy) * nblock->getNx() + ix * 2 + startFx;
            std::size_t Fid1 = Fid0 + 1;
            std::size_t Fid2 = Fid0 + nblock->getNx();
            std::size_t Fid3 = Fid2 + 1;
            comm.RecvCells.push_back(Cid);
            comm.SendCells.emplace_back(IntpSource<2>{Fid0, Fid1, Fid2, Fid3});
          }
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
  for (Block2D<T> &block : _Blocks) {
    std::uint8_t blocklevel = block.getLevel();
    std::vector<IntpBlockComm<T, 2>> &Communicators = block.getIntpBlockComm();
    Communicators.clear();
    for (Block2D<T> *nblock : block.getNeighbors()) {
      // find block of blocklevel-1
      if (nblock->getLevel() == blocklevel - 1) {
        Communicators.emplace_back(nblock);
        IntpBlockComm<T, 2> &comm = Communicators.back();
        // vox size
        const T Cvoxsize = nblock->getVoxelSize();
        const T Fvoxsize = block.getVoxelSize();
        // get intersection
        const AABB<T, 2> intsec = getIntersection(block, nblock->getBaseBlock());
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        // get start index of intsec in nblock
        Vector<T, 2> startC = intsec.getMin() - nblock->getMin();
        // shift 1 voxel to left bottom for interpolation
        int startCx = static_cast<int>(std::round(startC[0] / Cvoxsize)) - 1;
        int startCy = static_cast<int>(std::round(startC[1] / Cvoxsize)) - 1;
        // start index of intsec in FBlock
        Vector<T, 2> startF = intsec.getMin() - block.getMin();
        int startFx = static_cast<int>(std::round(startF[0] / Fvoxsize));
        int startFy = static_cast<int>(std::round(startF[1] / Fvoxsize));
        for (int iy = 0; iy < CNy; ++iy) {
          for (int ix = 0; ix < CNx; ++ix) {
            std::size_t Cid0 = (iy + startCy) * nblock->getNx() + ix + startCx;
            std::size_t Cid1 = Cid0 + 1;
            std::size_t Cid2 = Cid0 + nblock->getNx();
            std::size_t Cid3 = Cid2 + 1;
            std::size_t Fid = (iy * 2 + startFy) * block.getNx() + ix * 2 + startFx;

            // shift 1 voxel upward(+y direction)
            std::size_t Cid0_ = Cid0 + nblock->getNx();
            std::size_t Cid1_ = Cid1 + nblock->getNx();
            std::size_t Cid2_ = Cid2 + nblock->getNx();
            std::size_t Cid3_ = Cid3 + nblock->getNx();
            std::size_t Fid_ = Fid + block.getNx();

            // 0
            comm.RecvCells.push_back(Fid);
            comm.SendCells.emplace_back(IntpSource<2>{Cid0, Cid1, Cid2, Cid3});
            // comm.InterpWeights.emplace_back(
            //   InterpWeight<T, 2>{T(0.0625), T(0.1875), T(0.1875), T(0.5625)});

            // 1, shift along +x direction
            Cid0 += 1;
            Cid1 += 1;
            Cid2 += 1;
            Cid3 += 1;
            Fid += 1;
            comm.RecvCells.push_back(Fid);
            comm.SendCells.emplace_back(IntpSource<2>{Cid0, Cid1, Cid2, Cid3});
            // comm.InterpWeights.emplace_back(
            //   InterpWeight<T, 2>{T(0.1875), T(0.0625), T(0.5625), T(0.1875)});

            // 2
            comm.RecvCells.push_back(Fid_);
            comm.SendCells.emplace_back(IntpSource<2>{Cid0_, Cid1_, Cid2_, Cid3_});
            // comm.InterpWeights.emplace_back(
            //   InterpWeight<T, 2>{T(0.1875), T(0.5625), T(0.0625), T(0.1875)});

            // 3, shift along +x direction
            Cid0_ += 1;
            Cid1_ += 1;
            Cid2_ += 1;
            Cid3_ += 1;
            Fid_ += 1;
            comm.RecvCells.push_back(Fid_);
            comm.SendCells.emplace_back(IntpSource<2>{Cid0_, Cid1_, Cid2_, Cid3_});
            // comm.InterpWeights.emplace_back(
            //   InterpWeight<T, 2>{T(0.5625), T(0.1875), T(0.1875), T(0.0625)});
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
  InitComm();
  InitAverComm();
  InitIntpComm();
}

#ifdef MPI_ENABLED

template <typename T>
void BlockGeometry2D<T>::InitMPIComm(BlockGeometryHelper2D<T> &GeoHelper) {
  for (Block2D<T> &block : _Blocks) {
    MPIBlockComm &MPIComm = block.getMPIBlockComm();
    MPIComm.clear();
    std::uint8_t blocklevel = block.getLevel();
    // base block of block
    const BasicBlock<T, 2> &baseblock = block.getBaseBlock();
    // get the first layer of overlapped cells(counted from inside to outside)
    // int overlap = (blocklevel != std::uint8_t(0)) ? 2 : 1;
    BasicBlock<T, 2> baseblock_ext1 = baseblock.getExtBlock(1);
    // find neighbors
    std::vector<std::pair<int, int>> &nbrs =
      GeoHelper.getMPIBlockNbrs(block.getBlockId());
    for (const std::pair<int, int> &nbr : nbrs) {
      // check if 2 blocks are of the same level
      const BasicBlock<T, 2> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel) {
        BasicBlock<T, 2> nbaseblock_ext1 = nbaseblock.getExtBlock(1);
        // init sender
        MPIComm.Senders.emplace_back(nbr.first, nbaseblock.getBlockId());
        MPIBlockSendStru &sender = MPIComm.Senders.back();
        block.getCellIdx(baseblock, nbaseblock_ext1, sender.SendCells);
        // init receiver
        MPIComm.Recvers.emplace_back(nbr.first, nbaseblock.getBlockId());
        MPIBlockRecvStru &recver = MPIComm.Recvers.back();
        block.getCellIdx(nbaseblock, baseblock_ext1, recver.RecvCells);

        block._NeedMPIComm = true;
      }
    }
  }
}

template <typename T>
void BlockGeometry2D<T>::InitMPIAverComm(BlockGeometryHelper2D<T> &GeoHelper) {
  for (Block2D<T> &block : _Blocks) {
    MPIIntpBlockComm<T, 2> &MPIComm = block.getMPIAverBlockComm();
    MPIComm.clear();
    std::uint8_t blocklevel = block.getLevel();
    const BasicBlock<T, 2> &baseblock = block.getBaseBlock();
    // first layer of overlapped cells(counted from inside to outside)
    BasicBlock<T, 2> baseblock_ext1 = baseblock.getExtBlock(1);
    std::vector<std::pair<int, int>> &nbrs =
      GeoHelper.getMPIBlockNbrs(block.getBlockId());
    for (const std::pair<int, int> &nbr : nbrs) {
      const BasicBlock<T, 2> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel - 1) {
        BasicBlock<T, 2> nbaseblock_ext1 = nbaseblock.getExtBlock(1);
        // init sender
        MPIComm.Senders.emplace_back(nbr.first, nbaseblock.getBlockId());
        MPIIntpBlockSendStru<T, 2> &sender = MPIComm.Senders.back();
        // vox size
        const T Cvoxsize = nbaseblock.getVoxelSize();
        const T Fvoxsize = block.getVoxelSize();
        // init sender
        // get intersection
        const AABB<T, 2> intsec = getIntersection(baseblock, nbaseblock_ext1);
        // use coarse grid size here, convient for calculating fine cell index
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        // start index of intsec in block
        Vector<T, 2> startF = intsec.getMin() - block.getMin();
        int startFx = static_cast<int>(std::round(startF[0] / Fvoxsize));
        int startFy = static_cast<int>(std::round(startF[1] / Fvoxsize));
        for (int iy = 0; iy < CNy; ++iy) {
          for (int ix = 0; ix < CNx; ++ix) {
            std::size_t Fid0 = (iy * 2 + startFy) * block.getNx() + ix * 2 + startFx;
            std::size_t Fid1 = Fid0 + 1;
            std::size_t Fid2 = Fid0 + block.getNx();
            std::size_t Fid3 = Fid2 + 1;
            sender.SendCells.emplace_back(IntpSource<2>{Fid0, Fid1, Fid2, Fid3});
          }
        }
        block._NeedMPIComm = true;
      } else if (nbaseblock.getLevel() == blocklevel + 1) {
        // init receiver
        MPIComm.Recvers.emplace_back(nbr.first, nbaseblock.getBlockId());
        MPIBlockRecvStru &recver = MPIComm.Recvers.back();
        // vox size
        const T Cvoxsize = block.getVoxelSize();
        // const T Fvoxsize = nbaseblock.getVoxelSize();
        // get intersection
        const AABB<T, 2> intsec = getIntersection(baseblock_ext1, nbaseblock);
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        // start index of intsec in block
        Vector<T, 2> startC = intsec.getMin() - block.getMin();
        int startCx = static_cast<int>(std::round(startC[0] / Cvoxsize));
        int startCy = static_cast<int>(std::round(startC[1] / Cvoxsize));
        for (int iy = 0; iy < CNy; ++iy) {
          for (int ix = 0; ix < CNx; ++ix) {
            std::size_t Cid = (iy + startCy) * block.getNx() + ix + startCx;
            recver.RecvCells.push_back(Cid);
          }
        }
        block._NeedMPIComm = true;
      } else if (nbaseblock.getLevel() > blocklevel + 1 ||
                 nbaseblock.getLevel() < blocklevel - 1) {
        std::cerr << "[BlockGeometry2D<T>::InitMPIAverComm] Error: block level "
                     "difference larger than 1"
                  << std::endl;
      }
    }
  }
}

template <typename T>
void BlockGeometry2D<T>::InitMPIIntpComm(BlockGeometryHelper2D<T> &GeoHelper) {
  for (Block2D<T> &block : _Blocks) {
    MPIIntpBlockComm<T, 2> &MPIComm = block.getMPIIntpBlockComm();
    MPIComm.clear();
    std::uint8_t blocklevel = block.getLevel();
    const BasicBlock<T, 2> &baseblock = block.getBaseBlock();
    // 2 layers of overlapped cells(counted from inside to outside)
    // BasicBlock<T, 2> baseblock_ext2 = baseblock.getExtBlock(2);
    std::vector<std::pair<int, int>> &nbrs =
      GeoHelper.getMPIBlockNbrs(block.getBlockId());
    for (const std::pair<int, int> &nbr : nbrs) {
      const BasicBlock<T, 2> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel + 1) {
        BasicBlock<T, 2> nbaseblock_ext2 = nbaseblock.getExtBlock(2);
        // init sender
        MPIComm.Senders.emplace_back(nbr.first, nbaseblock.getBlockId());
        MPIIntpBlockSendStru<T, 2> &sender = MPIComm.Senders.back();
        // vox size
        const T Cvoxsize = block.getVoxelSize();
        // const T Fvoxsize = nbaseblock.getVoxelSize();
        // init sender
        // get intersection
        const AABB<T, 2> intsec = getIntersection(baseblock, nbaseblock_ext2);
        // use coarse grid size here, convient for calculating coarse cell index
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        // start index of intsec in block
        Vector<T, 2> startC = intsec.getMin() - block.getMin();
        // shift 1 voxel to left bottom for interpolation
        int startCx = static_cast<int>(std::round(startC[0] / Cvoxsize)) - 1;
        int startCy = static_cast<int>(std::round(startC[1] / Cvoxsize)) - 1;
        for (int iy = 0; iy < CNy; ++iy) {
          for (int ix = 0; ix < CNx; ++ix) {
            std::size_t Cid0 = (iy + startCy) * block.getNx() + ix + startCx;
            std::size_t Cid1 = Cid0 + 1;
            std::size_t Cid2 = Cid0 + block.getNx();
            std::size_t Cid3 = Cid2 + 1;

            std::size_t Cid0_ = Cid0 + block.getNx();
            std::size_t Cid1_ = Cid1 + block.getNx();
            std::size_t Cid2_ = Cid2 + block.getNx();
            std::size_t Cid3_ = Cid3 + block.getNx();

            // 0
            sender.SendCells.emplace_back(IntpSource<2>{Cid0, Cid1, Cid2, Cid3});

            // 1
            Cid0 += 1;
            Cid1 += 1;
            Cid2 += 1;
            Cid3 += 1;
            sender.SendCells.emplace_back(IntpSource<2>{Cid0, Cid1, Cid2, Cid3});

            // 2
            sender.SendCells.emplace_back(IntpSource<2>{Cid0_, Cid1_, Cid2_, Cid3_});

            // 3
            Cid0_ += 1;
            Cid1_ += 1;
            Cid2_ += 1;
            Cid3_ += 1;
            sender.SendCells.emplace_back(IntpSource<2>{Cid0_, Cid1_, Cid2_, Cid3_});
          }
        }
        block._NeedMPIComm = true;
      } else if (nbaseblock.getLevel() == blocklevel - 1) {
        // init receiver
        MPIComm.Recvers.emplace_back(nbr.first, nbaseblock.getBlockId());
        MPIBlockRecvStru &recver = MPIComm.Recvers.back();
        // vox size
        const T Cvoxsize = nbaseblock.getVoxelSize();
        const T Fvoxsize = block.getVoxelSize();
        // get intersection
        const AABB<T, 2> intsec = getIntersection(block, nbaseblock);
        // use coarse grid size here, convient for calculating fine cell index
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        // start index of intsec in block
        Vector<T, 2> startF = intsec.getMin() - block.getMin();
        int startFx = static_cast<int>(std::round(startF[0] / Fvoxsize));
        int startFy = static_cast<int>(std::round(startF[1] / Fvoxsize));
        for (int iy = 0; iy < CNy; ++iy) {
          for (int ix = 0; ix < CNx; ++ix) {
            std::size_t Fid0 = (iy * 2 + startFy) * block.getNx() + ix * 2 + startFx;
            std::size_t Fid1 = Fid0 + 1;
            std::size_t Fid2 = Fid0 + block.getNx();
            std::size_t Fid3 = Fid2 + 1;
            recver.RecvCells.push_back(Fid0);
            recver.RecvCells.push_back(Fid1);
            recver.RecvCells.push_back(Fid2);
            recver.RecvCells.push_back(Fid3);
          }
        }
        block._NeedMPIComm = true;
      } else if (nbaseblock.getLevel() > blocklevel + 1 ||
                 nbaseblock.getLevel() < blocklevel - 1) {
        std::cerr << "[BlockGeometry2D<T>::InitMPIIntpComm] Error: block level "
                     "difference larger than 1"
                  << std::endl;
      }
    }
  }
}

template <typename T>
void BlockGeometry2D<T>::InitAllMPIComm(BlockGeometryHelper2D<T> &GeoHelper) {
  InitMPIComm(GeoHelper);
  InitMPIAverComm(GeoHelper);
  InitMPIIntpComm(GeoHelper);
}

#endif

// BlockGeometryHelper2D

template <typename T>
BlockGeometryHelper2D<T>::BlockGeometryHelper2D(int Nx, int Ny, const AABB<T, 2> &AABBs,
                                                T voxelSize, int blockcelllen,
                                                std::uint8_t llimit, int ext)
    : BasicBlock<T, 2>(voxelSize, AABBs.getExtended(Vector<T, 2>{voxelSize * ext}),
                       AABB<int, 2>(Vector<int, 2>{0}, Vector<int, 2>{Nx + 1, Ny + 1})),
      _BaseBlock(voxelSize, AABBs,
                 AABB<int, 2>(Vector<int, 2>{1}, Vector<int, 2>{Nx, Ny})),
      BlockCellLen(blockcelllen), Ext(ext), _LevelLimit(llimit), _MaxLevel(std::uint8_t(0)), 
      _Exchanged(true), _IndexExchanged(true) {
  if (BlockCellLen < 4) {
    std::cerr << "BlockGeometryHelper2D<T>, BlockCellLen < 4" << std::endl;
  }
  CellsNx = _BaseBlock.getNx() / BlockCellLen;
  CellsNy = _BaseBlock.getNy() / BlockCellLen;
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
  std::vector<BasicBlock<T, 2>> &BasicBlocks = getAllOldBasicBlocks();
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
      Vector<int, 2> Ext = BlockCellLen * Vector<int, 2>{Nx, Ny};
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
  // check
  if (OptProcNum < 1) {
    std::cerr << "[BlockGeometryHelper2D<T>::AdaptiveOptimization]: OptProcNum < 1"
              << std::endl;
  }
  if (MaxProcNum == -1) {
    MaxProcNum = 2 * OptProcNum;
  }
  if (MaxProcNum < 1) {
    std::cerr << "[BlockGeometryHelper2D<T>::AdaptiveOptimization]: MaxProcNum < 1"
              << std::endl;
  }
  // get new basic blocks
  std::vector<BasicBlock<T, 2>> &BasicBlocks = getAllBasicBlocks();
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
  for (std::size_t i = 1; i < StdDevs.size(); ++i) {
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
  MPI_RANK(0)
  std::cout << "Optimization result: " << BasicBlocks.size()
            << " Blocks with stdDev: " << minStdDev << std::endl;
}

template <typename T>
void BlockGeometryHelper2D<T>::Optimize(int ProcessNum, bool enforce) {
  Optimize(getAllBasicBlocks(), ProcessNum, enforce);
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

  if (enforce && Blocks.size() < static_cast<std::size_t>(ProcessNum)) {
    // sort blocks by number of points, large blocks first
    std::sort(Blocks.begin(), Blocks.end(),
              [](const BasicBlock<T, 2> &a, const BasicBlock<T, 2> &b) {
                return a.getN() > b.getN();
              });
    // make size = ProcessNum
    std::size_t count = Blocks.size();
    while (count < static_cast<std::size_t>(ProcessNum) && it != Blocks.end()) {
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
  for (std::size_t i = 0; i < Blocks.size(); ++i) {
    Blocks[i].setBlockId(i);
  }
}

template <typename T>
void BlockGeometryHelper2D<T>::LoadBalancing(int ProcessNum) {
  std::vector<std::vector<int>> &BlockIndex = getAllOldBlockIndices();
  BlockIndex.clear();
  std::vector<std::vector<BasicBlock<T, 2> *>> &BlockIndexPtr = getAllOldBlockIndexPtrs();
  BlockIndexPtr.clear();
  _IndexExchanged = !_IndexExchanged;
#ifdef MPI_ENABLED
  BlockIndex.resize(ProcessNum);
  // map block's cell num to block id
  std::vector<std::pair<std::size_t, int>> BlockId_CellNum;
  for (BasicBlock<T, 2> &block : getAllBasicBlocks()) {
    BlockId_CellNum.emplace_back(std::make_pair(block.getN(), block.getBlockId()));
  }
  // sort blocks by number of cells in descending order
  std::sort(BlockId_CellNum.begin(), BlockId_CellNum.end(),
            [](const std::pair<std::size_t, int> &a,
               const std::pair<std::size_t, int> &b) { return a.first > b.first; });
  // divide blocks into ProcessNum parts, each part has number of cells as close as
  // possible a greedy algorithm
  for (std::pair<std::size_t, int> &pair : BlockId_CellNum) {
    int blockid = pair.second;
    // in std::vector<std::vector<int>> BlockIndex, find the part with smallest sum
    std::size_t minsum = std::numeric_limits<std::size_t>::max();
    int minidx = 0;
    for (int i = 0; i < ProcessNum; ++i) {
      std::size_t sum = 0;
      for (int id : BlockIndex[i]) {
        sum += getAllBasicBlock(static_cast<std::size_t>(id)).getN();
      }
      if (sum < minsum) {
        minsum = sum;
        minidx = i;
      }
    }
    // add blockid to vector with minidx
    BlockIndex[minidx].push_back(blockid);
  }
#else
  BlockIndex.resize(1);
  for (BasicBlock<T, 2> &block : getAllBasicBlocks()) {
    BlockIndex[0].push_back(block.getBlockId());
  }
#endif
  BlockIndexPtr.resize(BlockIndex.size());
  for (std::size_t i = 0; i < BlockIndex.size(); ++i) {
    for (int id : BlockIndex[i]) {
      BlockIndexPtr[i].push_back(&getAllBasicBlock(static_cast<std::size_t>(id)));
    }
  }

#ifdef MPI_ENABLED
  SetupMPINbrs();
  // print info
  MPI_RANK(0)
  std::cout << "[LoadBalancing result]: " << "\n";
  std::cout << "Rank:     ";
  for (std::size_t i = 0; i < BlockIndex.size(); ++i) {
    std::cout << i << " | " ;
  }
  std::cout << std::endl;
  std::cout << "BlockNum: ";
  for (std::size_t i = 0; i < BlockIndex.size(); ++i) {
    std::cout << BlockIndex[i].size() << " | ";
  }
  std::cout << std::endl;
#endif
}

#ifdef MPI_ENABLED

template <typename T>
void BlockGeometryHelper2D<T>::SetupMPINbrs() {
  _MPIBlockNbrs.clear();
  _MPIBlockNbrs.resize(getAllBasicBlocks().size(), std::vector<std::pair<int, int>>{});
  for (std::size_t iRank = 0; iRank < getAllBlockIndices().size(); ++iRank) {
    for (int blockid : getAllBlockIndices()[iRank]) {
      std::vector<std::pair<int, int>> &nbrsvec = _MPIBlockNbrs[blockid];
      const BasicBlock<T, 2> block = getAllBasicBlock(static_cast<std::size_t>(blockid)).getExtBlock(1);
      // for all blocks in all ranks except iRank
      for (std::size_t nRank = 0; nRank < getAllBlockIndices().size(); ++nRank) {
        if (nRank != iRank) {
          for (int nblockid : getAllBlockIndices()[nRank]) {
            const BasicBlock<T, 2> nblock = getAllBasicBlock(static_cast<std::size_t>(nblockid)).getExtBlock(1);
            if (isOverlapped(block, nblock))
              nbrsvec.push_back(std::make_pair(nRank, nblockid));
          }
        }
      }
    }
  }
}

#endif

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
                  util::isFlag(static_cast<CellTagType>(_BlockCellTags[ncellid]), static_cast<CellTagType>(BlockCellTag::none))) {
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
                util::isFlag(static_cast<CellTagType>(_BlockCellTags[ncellid]), static_cast<CellTagType>(BlockCellTag::none))) {
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
    if (util::isFlag(static_cast<CellTagType>(_BlockCellTags[i]), static_cast<CellTagType>(BlockCellTag::refine))) {
      _BlockCells[i].refine();
      // reset tag
      _BlockCellTags[i] = BlockCellTag::none;
    }
  }
}

// template <typename T>
// void BlockGeometry2D<T>::InitIntpComm() {
//   for (int i = 0; i < _Blocks.size(); ++i) {
//     // getCellIdx() called from ext block
//     Block2D<T> &block = _Blocks[i];
//     std::uint8_t blocklevel = block.getLevel();
//     std::vector<IntpBlockComm<T, 2>> &Communicators = block.getIntpBlockComm();
//     Communicators.clear();
//     for (Block2D<T> *nblock : block.getNeighbors()) {
//       // find block of blocklevel-1
//       if (nblock->getLevel() == blocklevel - 1) {
//         Communicators.emplace_back(nblock);
//         IntpBlockComm<T, 2> &comm = Communicators.back();
//         // init recv cells
//         // high level block recv from lower using 2^(diff_level) layers of overlapped
//         // cells
//         block.getCellIdx(block, nblock->getBaseBlock(), comm.RecvCells);
//         // get refined basicblock of nblock
//         BasicBlock<T, 2> Rnbaseblock = nblock->getBaseBlock().getRefinedBlock();
//         // buffer container to store virtual cell indices
//         std::vector<std::size_t> VSendCells;
//         // get virtual (refined) cell indices on refined basicblock of nblock
//         Rnbaseblock.getCellIdx(block, Rnbaseblock, VSendCells);
//         // get real cell indices from virtual cell indices
//         T halfvoxsize = nblock->getVoxelSize() / T(2);
//         for (std::size_t id : VSendCells) {
//           // refined cell location
//           Vector<T, 2> Rloc = Rnbaseblock.getLoc_t(id);
//           // get dist of Rloc and MinCenter of nblock
//           Vector<T, 2> dist = Rloc - nblock->getMinCenter();
//           // find nearest real cell location
//           int x = static_cast<int>(dist[0] / nblock->getVoxelSize());
//           int y = static_cast<int>(dist[1] / nblock->getVoxelSize());
//           T remx = dist[0] - x * nblock->getVoxelSize();
//           T remy = dist[1] - y * nblock->getVoxelSize();
//           int scheme = 0;
//           if (remx < halfvoxsize) {
//             scheme += 1;
//           }
//           if (remy < halfvoxsize) {
//             scheme += 2;
//           }
//           /*
//           2 3
//           0 1
//           */
//           std::size_t id0 = nblock->getIndex(Vector<int, 2>{x, y});
//           std::size_t id1 = nblock->getIndex(Vector<int, 2>{x + 1, y});
//           std::size_t id2 = nblock->getIndex(Vector<int, 2>{x, y + 1});
//           std::size_t id3 = nblock->getIndex(Vector<int, 2>{x + 1, y + 1});
//           comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//           switch (scheme) {
//             case 0: {
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(0.0625), T(0.1875), T(0.1875), T(0.5625)});
//             } break;
//             case 1: {
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(0.1875), T(0.0625), T(0.5625), T(0.1875)});
//             } break;
//             case 2: {
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(0.1875), T(0.5625), T(0.0625), T(0.1875)});
//             } break;
//             case 3: {
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(0.5625), T(0.1875), T(0.1875), T(0.0625)});
//             } break;
//           }
//         }
//       } else if (nblock->getLevel() < blocklevel - 1) {
//         std::cerr << "[BlockGeometry2D<T>::InitIntpComm] Error: block level difference
//         "
//                      "larger than 1"
//                   << std::endl;
//       }
//     }
//   }
// }

// template <typename T>
// void BlockGeometry2D<T>::InitAverComm(int highlevelovlap) {
//   for (int i = 0; i < _Blocks.size(); ++i) {
//     // getCellIdx() called from ext block
//     Block2D<T> &block = _Blocks[i];
//     // get the first layer of overlapped cells(counted from inside to outside)
//     BasicBlock<T, 2> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
//     std::uint8_t blocklevel = block.getLevel();
//     std::vector<IntpBlockComm<T, 2>> &Communicators = block.getAverageBlockComm();
//     Communicators.clear();
//     for (Block2D<T> *nblock : block.getNeighbors()) {
//       // find block of blocklevel+1
//       if (nblock->getLevel() == blocklevel + 1) {
//         Communicators.emplace_back(nblock);
//         IntpBlockComm<T, 2> &comm = Communicators.back();
//         // init recv cells
//         // low level block only recv from higher using the first layer of overlapped
//         cells block.getCellIdx(baseblock_ext1, nblock->getBaseBlock(), comm.RecvCells);
//         // get coarsened basicblock of nblock
//         BasicBlock<T, 2> Cnbaseblock = nblock->getBaseBlock().getCoasenedBlock();
//         // buffer container to store virtual cell indices
//         // TODO: more efficient implementation without using virtual cells?
//         std::vector<std::size_t> VSendCells;
//         // get virtual (coarse) cell indices on coarsened basicblock of nblock
//         Cnbaseblock.getCellIdx(baseblock_ext1, Cnbaseblock, VSendCells);
//         // get real cell indices from virtual cell indices
//         for (std::size_t id : VSendCells) {
//           // virtual (coarse) cell location
//           Vector<int, 2> Clocidx = Cnbaseblock.getLoc(id);
//           // real cell location
//           Vector<int, 2> locidx0 = Clocidx * 2 + Vector<int, 2>{highlevelovlap};
//           std::size_t id0 = nblock->getIndex(locidx0);
//           Vector<int, 2> locidx1 = locidx0 + Vector<int, 2>{1, 0};
//           std::size_t id1 = nblock->getIndex(locidx1);
//           Vector<int, 2> locidx2 = locidx0 + Vector<int, 2>{0, 1};
//           std::size_t id2 = nblock->getIndex(locidx2);
//           Vector<int, 2> locidx3 = locidx0 + Vector<int, 2>{1, 1};
//           std::size_t id3 = nblock->getIndex(locidx3);
//           comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//         }
//       } else if (nblock->getLevel() > blocklevel + 1) {
//         std::cerr << "[BlockGeometry2D<T>::InitAverComm] Error: block level difference
//         "
//                      "larger than 1"
//                   << std::endl;
//       }
//     }
//   }
// }

// template <typename T>
// void BlockGeometry2D<T>::InitAverComm2() {
//   for (int i = 0; i < _Blocks.size(); ++i) {
//     // getCellIdx() called from ext block
//     Block2D<T> &block = _Blocks[i];
//     BasicBlock<T, 2> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
//     std::uint8_t blocklevel = block.getLevel();
//     std::vector<IntpBlockComm<T, 2>> &Communicators = block.getAverageBlockComm();
//     for (Block2D<T> *nblock : block.getNeighbors()) {
//       // find block of blocklevel+1
//       if (nblock->getLevel() == blocklevel + 1) {
//         Communicators.emplace_back(nblock);
//         IntpBlockComm<T, 2> &comm = Communicators.back();
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
//               IntpSource<2>{Sends[id], Sends[id + 1], Sends[id + 2], Sends[id + 3]});
//           }
//         } else if (Align_dir[1] == 0) {
//           // horizontal
//           int halfsize = Sends.size() / 2;
//           for (int id = 0; id < halfsize; id += 2) {
//             comm.SendCells.emplace_back(IntpSource<2>{
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
//     std::vector<IntpBlockComm<T, 2>> &Communicators = block.getIntpBlockComm();
//     for (Block2D<T> *nblock : block.getNeighbors()) {
//       // find block of blocklevel-1
//       if (nblock->getLevel() == blocklevel - 1) {
//         Communicators.emplace_back(nblock);
//         IntpBlockComm<T, 2> &comm = Communicators.back();
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
//             comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(1.25), T(0), T(0),
//             T(-0.25)}); for (int j = 0; j < VSendCells.size() - 1; j += 2) {
//               comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(0.9375), T(0.3125), T(-0.1875), T(-0.0625)});
//               comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
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
//             comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
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
//             comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(0), T(1.25), T(-0.25),
//             T(0)}); for (int j = 0; j < VSendCells.size() - 1; j += 2) {
//               comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(-0.1875), T(0.9375), T(-0.0625), T(0.3125)});
//               comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
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
//             comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
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
//             comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(0), T(-0.25), T(1.25),
//             T(0)}); for (int j = 0; j < VSendCells.size() - 1; j += 2) {
//               comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(-0.1875), T(-0.0625), T(0.9375), T(0.3125)});
//               comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
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
//             comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
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
//             comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(1.25), T(0), T(-0.25),
//             T(0)}); for (int j = 0; j < VSendCells.size() - 1; j += 2) {
//               comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//               comm.InterpWeights.emplace_back(
//                 InterpWeight<T, 2>{T(0.9375), T(-0.1875), T(0.3125), T(-0.0625)});
//               comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
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
//             comm.SendCells.emplace_back(IntpSource<2>{id0, id1, id2, id3});
//             comm.InterpWeights.emplace_back(InterpWeight<T, 2>{T(0), T(-0.25), T(1.25),
//             T(0)});
//           }
//         }
//       }
//     }
//   }
// }