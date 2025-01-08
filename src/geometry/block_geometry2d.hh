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
    : BasicBlock<T, 2>(voxelSize, block.getExtended(Vector<T, 2>{voxelSize*olap}),
                       idxblock.getExtended(Vector<int, 2>{olap}), blockid),
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
  const int overlap = _overlap;
  for (int y = overlap; y < BasicBlock<T, 2>::Mesh[1] - overlap; ++y) {
    for (int x = overlap; x < BasicBlock<T, 2>::Mesh[0] - overlap; ++x) {
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

template <typename T>
template <typename FieldType, typename LatSet>
void Block2D<T>::SetupBoundary(FieldType &field, typename FieldType::value_type fromvalue,
typename FieldType::value_type voidvalue, typename FieldType::value_type bdvalue) {
  // temp flag field store the transition flag
  GenericArray<bool> TransFlag(BasicBlock<T, 2>::N, false);
  const int overlap = _overlap;
  for (int y = overlap; y < BasicBlock<T, 2>::Mesh[1] - overlap; ++y) {
    for (int x = overlap; x < BasicBlock<T, 2>::Mesh[0] - overlap; ++x) {
      const std::size_t idx = x + y * BasicBlock<T, 2>::Mesh[0];
      if (field.get(idx) == fromvalue) {
        for (unsigned int i = 1; i < LatSet::q; ++i) {
          const std::size_t nbridx = idx + latset::c<LatSet>(i) * BasicBlock<T, 2>::Projection;
          if (field.get(nbridx) == voidvalue) {
            TransFlag.set(idx, true);
            break;
          }
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
                                    T voxelSize, int overlap, int blockXNum, int blockYNum)
    : BasicBlock<T, 2>(voxelSize, block.getExtended(Vector<T, 2>{voxelSize*overlap}),
                       AABB<int, 2>(Vector<int, 2>{0}, Vector<int, 2>{Nx - 1 + 2*overlap, Ny - 1 + 2*overlap})),
      _BaseBlock(voxelSize, block,
                 AABB<int, 2>(Vector<int, 2>{overlap}, Vector<int, 2>{Nx - 1 + overlap, Ny - 1 + overlap})),
      _overlap(overlap), _MaxLevel(std::uint8_t(0)) {
  CreateBlocks(blocknum, blockXNum, blockYNum);
  BuildBlockIndexMap();
  SetupNbrs();
  InitComm();
  PrintInfo();
}

template <typename T>
BlockGeometry2D<T>::BlockGeometry2D(BlockGeometryHelper2D<T> &GeoHelper, bool useHelperOlap)
    : BasicBlock<T, 2>(GeoHelper), _BaseBlock(GeoHelper.getBaseBlock()), 
      _overlap(GeoHelper.getOverlap()), _MaxLevel(GeoHelper.getMaxLevel()) {
  // create blocks from GeoHelper
  for (BasicBlock<T, 2> *baseblock : GeoHelper.getBasicBlocks()) {
    int overlap{};
    if (useHelperOlap) {
      overlap = _overlap;
    } else {
      overlap = (baseblock->getLevel() != std::uint8_t(0)) ? 2 : 1;
    }
    _Blocks.emplace_back(*baseblock, overlap);
  }
  BuildBlockIndexMap();
  SetupNbrs();
  InitAllComm();
#ifdef MPI_ENABLED
  GeoHelper.InitBlockGeometry2D();
  InitAllMPIComm(GeoHelper);
#endif
  PrintInfo();
}

template <typename T>
BlockGeometry2D<T>::BlockGeometry2D(const BlockReader<T,2>& blockreader, bool useReaderOlap) 
    : BasicBlock<T, 2>(blockreader.getBasicBlock()), _BaseBlock(blockreader.getBaseBlock()), 
      _overlap(1), _MaxLevel(blockreader.getMaxLevel()) {
  // create blocks from Block Reader
  int iblock{};
  for (const BasicBlock<T, 2> &baseblock : blockreader.getBlocks()) {
    int overlap{};
    if (useReaderOlap) {
      overlap = blockreader.getOverlaps()[iblock];
    } else {
      overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
    }
    _Blocks.emplace_back(baseblock, overlap);
    ++iblock;
  }
  BuildBlockIndexMap();
  SetupNbrs();
  InitAllComm();
  PrintInfo();
}

template <typename T>
void BlockGeometry2D<T>::PrintInfo() const {
  std::size_t cellnum = getTotalCellNum();
  // DO NOT use MPI_RANK(0) before getTotalCellNum() cause we have: mpi().barrier();
  MPI_RANK(0)
  std::cout << "[BlockGeometry2D]: "
            << "Total Cell Num: " << cellnum << std::endl;
#ifndef MPI_ENABLED
  std::cout << "BlockNum: " << _Blocks.size() << " with StdDev: "
            << ComputeBlockNStdDev(_BasicBlocks) << std::endl;
#endif
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
  BuildBlockIndexMap();
  SetupNbrs();
  InitAllComm();
#ifdef MPI_ENABLED
  // GeoHelper.InitBlockGeometry2D();
  InitAllMPIComm(GeoHelper);
#endif
  PrintInfo();
}

template <typename T>
std::size_t BlockGeometry2D<T>::getTotalCellNum() const {
  std::size_t sum{};
  for (const Block2D<T> &block : _Blocks) sum += block.getN();
#ifdef MPI_ENABLED
  std::size_t Result{};
  mpi().barrier();
	mpi().reduce(sum, Result, MPI_SUM);
  return Result;
#else
  return sum;
#endif
}

template <typename T>
std::size_t BlockGeometry2D<T>::getBaseCellNum() const {
  std::size_t sum{};
  for (const Block2D<T> &block : _Blocks) sum += block.getBaseBlock().getN();
#ifdef MPI_ENABLED
  std::size_t Result{};
  mpi().barrier();
	mpi().reduce(sum, Result, MPI_SUM);
  return Result;
#else
  return sum;
#endif
}

template <typename T>
void BlockGeometry2D<T>::DivideBlocks(int blocknum, int blockXNum, int blockYNum) {
  _BlockAABBs.clear();
  _BlockAABBs.reserve(blocknum);
  if (blockXNum == 0 || blockYNum == 0) {
    // default scheme
    DivideBlock2D(_BaseBlock, blocknum, _BlockAABBs);
  } else {
    // manually divide
    DivideBlock2D(_BaseBlock, blockXNum, blockYNum, _BlockAABBs);
  }
}

template <typename T>
void BlockGeometry2D<T>::CreateBlocks(int blocknum, int blockXNum, int blockYNum) {
  DivideBlocks(blocknum, blockXNum, blockYNum);

  _Blocks.clear();
  _Blocks.reserve(_BlockAABBs.size());
  _BasicBlocks.clear();
  _BasicBlocks.reserve(_BlockAABBs.size());
  // create blocks
  int blockid = 0;
  const T voxsize = _BaseBlock.getVoxelSize();
  for (const AABB<int, 2> &blockaabb : _BlockAABBs) {
    Vector<T, 2> MIN =
      (blockaabb.getMin() - Vector<int,2>{_overlap}) * voxsize + _BaseBlock.getMin();
    Vector<T, 2> MAX =
      (blockaabb.getMax() - Vector<int,2>{_overlap-1}) * voxsize + _BaseBlock.getMin();
    AABB<T, 2> aabb(MIN, MAX);
    _Blocks.emplace_back(aabb, blockaabb, blockid, voxsize, _overlap);
    _BasicBlocks.emplace_back(voxsize, aabb, blockaabb, blockid);
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
        if (isOverlapped(block.getBaseBlock().getExtBlock(1), blockn.getBaseBlock())) {
          nbrsvec.push_back(&blockn);
        }
      }
    }
  }
}

template <typename T>
void BlockGeometry2D<T>::InitComm() {
  int Tag{};
  for (Block2D<T> &block : _Blocks) {
    // (inner) overlapped cell communicators
    std::vector<SharedComm> &Comms = block.getCommunicator().Comm.Comms;
    // all overlapped cell communicators
    std::vector<SharedComm> &AllComms = block.getCommunicator().AllComm.Comms;
    Comms.clear();
    AllComms.clear();
    // get block with overlap 1, for Comms
    BasicBlock<T, 2> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
    // get block with overlap = _overlap, for AllComms
    BasicBlock<T, 2> baseblock_exto = block.getSelfBlock();

    std::uint8_t blocklevel = block.getLevel();
    for (Block2D<T> *nblock : block.getNeighbors()) {
      // check if 2 blocks are of the same level
      if (nblock->getLevel() == blocklevel) {
        // ------ add to Comms
        // get overlapped cells
        std::vector<std::size_t> Recvs;
        std::vector<std::size_t> Sends;
        block.getCellIdx(baseblock_ext1, nblock->getBaseBlock(), Recvs);
        nblock->getCellIdx(nblock->getBaseBlock(), baseblock_ext1, Sends);
        // exclude corner cells
        std::vector<std::size_t> CornerRecvs;
        std::vector<std::size_t> CornerSends;
        block.ExcludeCornerIdx(Recvs, Sends, CornerRecvs, CornerSends);
        // avoid empty communicator
        if (Recvs.size() > 0) {
          Comms.emplace_back(nblock->getBlockId(), Tag);
          ++Tag;
          SharedComm &comm = Comms.back();
          comm.setRecvSendIdx(Recvs, Sends);
          comm.Direction = getEdgeNbrDirection<2>(block.whichEdge(comm.SendRecvCells[1]));
        }
        // add corner cells to communicators
        if (CornerRecvs.size() > 0) {
          // add corner cells to communicators and find direction
          for (std::size_t i = 0; i < CornerRecvs.size(); ++i) {
            Comms.emplace_back(nblock->getBlockId(), Tag);
            ++Tag;
            SharedComm &commcorner = Comms.back();
            commcorner.SendRecvCells = {CornerSends[i], CornerRecvs[i]};
            commcorner.Direction = getCornerNbrDirection<2>(block.whichCorner(CornerRecvs[i]));
          }
        }

        // ------ add to AllComms
        // here we will NOT consider direction in ALLComms for now
        std::vector<std::size_t> AllRecvs;
        std::vector<std::size_t> AllSends;
        block.getCellIdx(baseblock_exto, nblock->getBaseBlock(), AllRecvs);
        nblock->getCellIdx(nblock->getBaseBlock(), baseblock_exto, AllSends);
        // add cell index to communicator
        AllComms.emplace_back(nblock->getBlockId());
        SharedComm &allcomm = AllComms.back();
        allcomm.setRecvSendIdx(AllRecvs, AllSends);
      }
    }
  }
}

template <typename T>
void BlockGeometry2D<T>::InitAverComm() {
  for (Block2D<T> &block : _Blocks) {
    // (inner) overlapped cell communicators
    std::vector<SharedComm> &Comms = block.getCommunicator().Comm.AverComm;
    // all overlapped cell communicators
    std::vector<SharedComm> &AllComms = block.getCommunicator().AllComm.AverComm;
    Comms.clear();
    AllComms.clear();
    // get block with overlap 1, for Comms
    BasicBlock<T, 2> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
    // get block with overlap = _overlap, for AllComms
    BasicBlock<T, 2> baseblock_exto = block.getSelfBlock();

    std::uint8_t blocklevel = block.getLevel();
    for (Block2D<T> *nblock : block.getNeighbors()) {
      // find block of blocklevel+1
      if (nblock->getLevel() == blocklevel + 1) {
        // ------ add to Comms
        AddtoSharedAverComm(Comms, baseblock_ext1, block, nblock);
        // ------ add to AllComms
        AddtoSharedAverComm(AllComms, baseblock_exto, block, nblock);
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
    // (inner) overlapped cell communicators
    std::vector<SharedComm> &Comms = block.getCommunicator().Comm.IntpComm;
    // all overlapped cell communicators
    std::vector<SharedComm> &AllComms = block.getCommunicator().AllComm.IntpComm;
    Comms.clear();
    AllComms.clear();
    // get block with overlap = _overlap, for both Comms and AllComms
    // BasicBlock<T, 2> baseblock_exto = block.getSelfBlock();
    std::uint8_t blocklevel = block.getLevel();
    for (Block2D<T> *nblock : block.getNeighbors()) {
      // find block of blocklevel-1
      if (nblock->getLevel() == blocklevel - 1) {
        // ------ add to Comms
        AddtoSharedIntpComm(Comms, block, nblock);
        // ------ add to AllComms, for now it is the same as Comms
        AddtoSharedIntpComm(AllComms, block, nblock);
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
BlockGeometry2D<T>::BlockGeometry2D(BlockGeometryHelper2D<T> &GeoHelper, 
std::vector<BasicBlock<T, 2>>& BasicBlocks)
    : BasicBlock<T, 2>(GeoHelper), _BaseBlock(GeoHelper.getBaseBlock()), 
      _overlap(GeoHelper.getOverlap()), _MaxLevel(GeoHelper.getMaxLevel()) {
  // create blocks from GeoHelper
  for (BasicBlock<T, 2>& baseblock : BasicBlocks) {
    int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
    _Blocks.emplace_back(baseblock, overlap);
  }
  BuildBlockIndexMap();
  SetupNbrs();
  // this uses all blocks of the same level to init shared communicators
  // for mpi communication with direction info 
  InitComm();
}

template <typename T>
void BlockGeometry2D<T>::InitMPIComm(BlockGeometryHelper2D<T> &GeoHelper) {
  // mpi efficient send/recv using direction info for pop communication
  const BlockGeometry2D<T>& HelperBlockGeometry = GeoHelper.getBlockGeometry2D();
  for (Block2D<T> &block : _Blocks) {
    // efficient MPI Recvs/Sends using direction info for pop communication
    std::vector<DistributedComm>& Recvs = block.getCommunicator().DirRecvs;
    std::vector<DistributedComm>& Sends = block.getCommunicator().DirSends;
    Recvs.clear();
    Sends.clear();
    // find hblock in HelperBlockGeometry with:
    // the same blockid as block or
    // the same sendblockid in hblock's communicator
    for (const Block2D<T> &hblock : HelperBlockGeometry.getBlocks()) {
      // get all normal communicators of hblock
      const std::vector<SharedComm>& hComms = hblock.getCommunicator().Comm.Comms;

      if (hblock.getBlockId() == block.getBlockId()) {
        // find if nbr block of hblock is NOT in _Blocks(not in this rank)
        bool UseMPIComm = false;
        std::vector<const SharedComm*> SendhComms;
        for(const SharedComm& hComm : hComms) {
          if (!hasBlock(hComm.SendBlockId)) {
            UseMPIComm = true;
            SendhComms.push_back(&hComm);
          }
        }
        if (UseMPIComm) {
          // init MPI communicators
          for (const SharedComm* SendhComm : SendhComms) {
            Recvs.emplace_back(GeoHelper.whichRank(SendhComm->SendBlockId), SendhComm->SendBlockId, SendhComm->Tag);
            DistributedComm& mpirecv = Recvs.back();
            mpirecv.Direction = SendhComm->Direction;
            SendhComm->getRecvvector(mpirecv.Cells);
          }
        }
      } else {
        // find if hblock is NOT in _Blocks(not in this rank)
        if (!hasBlock(hblock.getBlockId())) {
          // find if block is the target of hblock's communicator
          for(const SharedComm& hComm : hComms) {
            if (hComm.SendBlockId == block.getBlockId()) {
              Sends.emplace_back(GeoHelper.whichRank(hblock.getBlockId()), hblock.getBlockId(), hComm.Tag);
              DistributedComm& mpisend = Sends.back();
              mpisend.Direction = hComm.Direction;
              hComm.getSendvector(mpisend.Cells);
            }
          }
        }
      }
    }
  }

  for (Block2D<T> &block : _Blocks) {
    // (inner) overlapped cell communicators
    DistributedCommSet &MPIComm = block.getCommunicator().MPIComm;
    // all overlapped cell communicators
    DistributedCommSet &AllMPIComm = block.getCommunicator().AllMPIComm;
    MPIComm.Recvs.clear();
    MPIComm.Sends.clear();
    AllMPIComm.Recvs.clear();
    AllMPIComm.Sends.clear();
    // base block of block
    const BasicBlock<T, 2> &baseblock = block.getBaseBlock();
    // get block with overlap 1, for Comms
    BasicBlock<T, 2> baseblock_ext1 = baseblock.getExtBlock(1);
    // get block with overlap = _overlap, for AllComms
    BasicBlock<T, 2> baseblock_exto = block.getSelfBlock();

    std::uint8_t blocklevel = block.getLevel();
    std::vector<std::pair<int, int>> &nbrs = GeoHelper.getMPIBlockNbrs(block.getBlockId());

    for (const std::pair<int, int> &nbr : nbrs) {
      // check if 2 blocks are of the same level
      const BasicBlock<T, 2> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel) {
        BasicBlock<T, 2> nbaseblock_ext1 = nbaseblock.getExtBlock(1);
        BasicBlock<T, 2> nbaseblock_exto = nbaseblock.getExtBlock(block.getOverlap());
        // ------ add to Comms
        // init receiver
        MPIComm.Recvs.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm& mpirecv = MPIComm.Recvs.back();
        // get overlapped recv/send cells
        block.getCellIdx(baseblock_ext1, nbaseblock, mpirecv.Cells);

        // init sender
        MPIComm.Sends.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm& mpisend = MPIComm.Sends.back();
        // get overlapped recv/send cells
        block.getCellIdx(baseblock, nbaseblock_ext1, mpisend.Cells);

        // ------ add to AllComms
        // init receiver
        AllMPIComm.Recvs.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm& allmpirecv = AllMPIComm.Recvs.back();
        // get overlapped recv/send cells
        block.getCellIdx(baseblock_exto, nbaseblock, allmpirecv.Cells);
        // init sender
        AllMPIComm.Sends.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm& allmpisend = AllMPIComm.Sends.back();
        // get overlapped recv/send cells
        block.getCellIdx(baseblock, nbaseblock_exto, allmpisend.Cells);

        block.getCommunicator()._NeedMPIComm = true;
      }
    }
  }
  // // add direction info for mpi send
  // mpi().barrier();
  // std::vector<std::vector<std::uint8_t>> SendBuffers(
  //   _Blocks.size(), std::vector<std::uint8_t>{});
  // std::vector<std::vector<std::uint8_t>> RecvBuffers(
  //   _Blocks.size(), std::vector<std::uint8_t>{});
  // std::size_t iblock{};
  // // --- send data --- we send recv direction here
  // std::vector<MPI_Request> SendRequests;
  // for (Block2D<T> &block : _Blocks) {
  //   const std::vector<DistributedComm> &Sends = block.getCommunicator().MPIComm.Recvs;
  //   std::vector<std::uint8_t>& SendBuffer = SendBuffers[iblock];
  //   SendBuffer.resize(Sends.size(), std::uint8_t{});
  //   for (std::size_t i = 0; i < Sends.size(); ++i) {
  //     const DistributedComm& mpisend = Sends[i];
  //     std::uint8_t& buffer = SendBuffer[i];
  //     buffer = static_cast<std::uint8_t>(mpisend.Direction);
  //     // non-blocking send
  //     MPI_Request request;
  //     mpi().iSend(&buffer, 1, mpisend.TargetRank, &request, mpisend.TargetBlockId);
  //     SendRequests.push_back(request);
  //   }
  //   ++iblock;
  // }
  // // --- receive data --- we receive send direction here
  // iblock = 0;
  // std::vector<MPI_Request> RecvRequests;
  // for (Block2D<T> &block : _Blocks) {
  //   const std::vector<DistributedComm> &Recvs = block.getCommunicator().MPIComm.Sends;
  //   std::vector<std::uint8_t>& RecvBuffer = RecvBuffers[iblock];
  //   RecvBuffer.resize(Recvs.size(), std::uint8_t{});
  //   for (std::size_t i = 0; i < Recvs.size(); ++i) {
  //     const DistributedComm& mpirecv = Recvs[i];
  //     std::uint8_t& buffer = RecvBuffer[i];
  //     // non-blocking recv
  //     MPI_Request request;
  //     mpi().iRecv(&buffer, 1, mpirecv.TargetRank, &request, block.getBlockId());
  //     RecvRequests.push_back(request);
  //   }
  //   ++iblock;
  // }
  // // wait for all send and recv requests to complete
  // MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // // --- wait and set data ---
  // iblock = 0;
  // for (Block2D<T> &block : _Blocks) {
  //   std::vector<DistributedComm> &Recvs = block.getCommunicator().MPIComm.Sends;
  //   const std::vector<std::uint8_t>& RecvBuffer = RecvBuffers[iblock];
  //   for (std::size_t i = 0; i < Recvs.size(); ++i) {
  //     DistributedComm& mpirecv = Recvs[i];
  //     mpirecv.Direction = static_cast<NbrDirection>(RecvBuffer[i]);
  //   }
  //   ++iblock;
  // }
  // mpi().barrier();
}

template <typename T>
void BlockGeometry2D<T>::InitMPIAverComm(BlockGeometryHelper2D<T> &GeoHelper) {
  for (Block2D<T> &block : _Blocks) {
    // (inner) overlapped cell communicators
    DistributedCommSet &MPIComm = block.getCommunicator().MPIComm;
    // all overlapped cell communicators
    DistributedCommSet &AllMPIComm = block.getCommunicator().AllMPIComm;
    MPIComm.AverRecvs.clear();
    MPIComm.AverSends.clear();
    AllMPIComm.AverRecvs.clear();
    AllMPIComm.AverSends.clear();

    // base block of block
    const BasicBlock<T, 2> &baseblock = block.getBaseBlock();
    // get block with overlap 1, for Comms
    BasicBlock<T, 2> baseblock_ext1 = baseblock.getExtBlock(1);
    // get block with overlap = _overlap, for AllComms
    BasicBlock<T, 2> baseblock_exto = block.getSelfBlock();

    std::uint8_t blocklevel = block.getLevel();
    std::vector<std::pair<int, int>> &nbrs = GeoHelper.getMPIBlockNbrs(block.getBlockId());

    for (const std::pair<int, int> &nbr : nbrs) {
      const BasicBlock<T, 2> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel + 1) {
        // init receiver
        MPIComm.AverRecvs.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &recver = MPIComm.AverRecvs.back();
        block.getCellIdx(baseblock_ext1, nbaseblock, recver.Cells);

        AllMPIComm.AverRecvs.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &Allrecver = AllMPIComm.AverRecvs.back();
        block.getCellIdx(baseblock_exto, nbaseblock, Allrecver.Cells);

        block.getCommunicator()._NeedMPIComm = true;
      } else if (nbaseblock.getLevel() == blocklevel - 1) {
        // init sender
        // virtual coarse block
        BasicBlock<T, 2> Cblock = block.getCoasenedBlock();

        BasicBlock<T, 2> nbaseblock_ext1 = nbaseblock.getExtBlock(1);
        int noverlap = nbaseblock.getLevel() == std::uint8_t(0) ? 1 : 2;
        BasicBlock<T, 2> nbaseblock_exto = nbaseblock.getExtBlock(noverlap);

        MPIComm.AverSends.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &sender = MPIComm.AverSends.back();
        std::vector<std::size_t> VSends;
        Cblock.getCellIdx(baseblock, nbaseblock_ext1, VSends);
        sender.Cells.reserve(VSends.size()*4);
        for (std::size_t vid : VSends) {
          std::vector<std::size_t> idxs;
          Cblock.getRefinedCellIdx(vid, idxs);
          sender.Cells.insert(sender.Cells.end(), idxs.begin(), idxs.end());          
        }

        AllMPIComm.AverSends.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &allsender = AllMPIComm.AverSends.back();
        std::vector<std::size_t> AllVSends;
        Cblock.getCellIdx(baseblock, nbaseblock_exto, AllVSends);
        allsender.Cells.reserve(AllVSends.size()*4);
        for (std::size_t vid : AllVSends) {
          std::vector<std::size_t> idxs;
          Cblock.getRefinedCellIdx(vid, idxs);
          allsender.Cells.insert(allsender.Cells.end(), idxs.begin(), idxs.end());          
        }

        block.getCommunicator()._NeedMPIComm = true;
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
    // (inner) overlapped cell communicators
    DistributedCommSet &MPIComm = block.getCommunicator().MPIComm;
    // all overlapped cell communicators
    DistributedCommSet &AllMPIComm = block.getCommunicator().AllMPIComm;
    MPIComm.IntpRecvs.clear();
    MPIComm.IntpSends.clear();
    AllMPIComm.IntpRecvs.clear();
    AllMPIComm.IntpSends.clear();

    // base block of block
    const BasicBlock<T, 2> &baseblock = block.getBaseBlock();
    // get block with overlap 1, for Comms
    // BasicBlock<T, 2> baseblock_ext1 = baseblock.getExtBlock(1);
    // get block with overlap = _overlap, for AllComms
    // BasicBlock<T, 2> baseblock_exto = block.getSelfBlock();

    std::uint8_t blocklevel = block.getLevel();
    std::vector<std::pair<int, int>> &nbrs = GeoHelper.getMPIBlockNbrs(block.getBlockId());

    for (const std::pair<int, int> &nbr : nbrs) {
      const BasicBlock<T, 2> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel - 1) {
        // init receiver
        // virtual coarse block
        BasicBlock<T, 2> Cblock = block.getCoasenedBlock();

        MPIComm.IntpRecvs.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &recver = MPIComm.IntpRecvs.back();
        AllMPIComm.IntpRecvs.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &allrecver = AllMPIComm.IntpRecvs.back();

        std::vector<std::size_t> VRecvs;
        Cblock.getCellIdx(block, nbaseblock, VRecvs);

        recver.Cells.reserve(VRecvs.size()*4);
        allrecver.Cells.reserve(VRecvs.size()*4);
        for (std::size_t vid : VRecvs) {
          std::vector<std::size_t> idxs;
          Cblock.getRefinedCellIdx(vid, idxs);
          recver.Cells.insert(recver.Cells.end(), idxs.begin(), idxs.end());
          allrecver.Cells.insert(allrecver.Cells.end(), idxs.begin(), idxs.end());        
        }

        block.getCommunicator()._NeedMPIComm = true;
      } else if (nbaseblock.getLevel() == blocklevel + 1) {
        // init sender
        BasicBlock<T, 2> nbaseblock_ext2 = nbaseblock.getExtBlock(2);
        
        MPIComm.IntpSends.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &sender = MPIComm.IntpSends.back();
        AllMPIComm.IntpSends.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &allsender = AllMPIComm.IntpSends.back();

        // vox size
        const T Cvoxsize = block.getVoxelSize();
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

        std::vector<std::size_t>& Sends = sender.Cells;
        std::vector<std::size_t>& AllSends = allsender.Cells;
        Sends.reserve(CNx * CNy * 16);
        AllSends.reserve(CNx * CNy * 16);

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
            Sends.insert(Sends.end(), {Cid0, Cid1, Cid2, Cid3});
            AllSends.insert(AllSends.end(), {Cid0, Cid1, Cid2, Cid3});

            // 1
            Cid0 += 1;
            Cid1 += 1;
            Cid2 += 1;
            Cid3 += 1;
            Sends.insert(Sends.end(), {Cid0, Cid1, Cid2, Cid3});
            AllSends.insert(AllSends.end(), {Cid0, Cid1, Cid2, Cid3});

            // 2
            Sends.insert(Sends.end(), {Cid0_, Cid1_, Cid2_, Cid3_});
            AllSends.insert(AllSends.end(), {Cid0_, Cid1_, Cid2_, Cid3_});

            // 3
            Cid0_ += 1;
            Cid1_ += 1;
            Cid2_ += 1;
            Cid3_ += 1;
            Sends.insert(Sends.end(), {Cid0_, Cid1_, Cid2_, Cid3_});
            AllSends.insert(AllSends.end(), {Cid0_, Cid1_, Cid2_, Cid3_});
          }
        }
        block.getCommunicator()._NeedMPIComm = true;
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

template <typename T>
void BlockGeometry2D<T>::AddtoSharedAverComm(std::vector<SharedComm> &Comms, 
const BasicBlock<T, 2>& baseblock_extx, Block2D<T>& block, Block2D<T> *nblock){
  // ------ add to Comms
  Comms.emplace_back(nblock->getBlockId());
  SharedComm &comm = Comms.back();
  // virtual coarse nblock
  BasicBlock<T, 2> nCblock = nblock->getCoasenedBlock();
  // get overlapped coarse cells
  std::vector<std::size_t> Recvs;
  std::vector<std::size_t> VSends;
  block.getCellIdx(baseblock_extx, nblock->getBaseBlock(), Recvs);
  nCblock.getCellIdx(nblock->getBaseBlock(), baseblock_extx, VSends);
  // get real send cell index
  std::vector<std::size_t> RSends;
  RSends.reserve(VSends.size()*4);
  for(std::size_t vid : VSends){
    std::vector<std::size_t> idxs;
    nCblock.getRefinedCellIdx(vid, idxs);
    RSends.insert(RSends.end(), idxs.begin(), idxs.end());
  }
  comm.setRecvSendIdx(Recvs, RSends);
}

template <typename T>
void BlockGeometry2D<T>::AddtoSharedIntpComm(std::vector<SharedComm> &Comms, Block2D<T>& block, Block2D<T> *nblock){
  // ------ add to Comms
  Comms.emplace_back(nblock->getBlockId());
  SharedComm &comm = Comms.back();
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

  std::vector<std::size_t>& SendRecvs = comm.SendRecvCells;
  // (4+1) * 4
  SendRecvs.reserve(CNx*CNy*20);

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
      SendRecvs.insert(SendRecvs.end(), {Cid0, Cid1, Cid2, Cid3});
      SendRecvs.push_back(Fid);
      // {T(0.0625), T(0.1875), T(0.1875), T(0.5625)});

      // 1, shift along +x direction
      Cid0 += 1;
      Cid1 += 1;
      Cid2 += 1;
      Cid3 += 1;
      Fid += 1;
      SendRecvs.insert(SendRecvs.end(), {Cid0, Cid1, Cid2, Cid3});
      SendRecvs.push_back(Fid);
      // {T(0.1875), T(0.0625), T(0.5625), T(0.1875)});

      // 2
      SendRecvs.insert(SendRecvs.end(), {Cid0_, Cid1_, Cid2_, Cid3_});
      SendRecvs.push_back(Fid_);
      // {T(0.1875), T(0.5625), T(0.0625), T(0.1875)});

      // 3, shift along +x direction
      Cid0_ += 1;
      Cid1_ += 1;
      Cid2_ += 1;
      Cid3_ += 1;
      Fid_ += 1;
      SendRecvs.insert(SendRecvs.end(), {Cid0_, Cid1_, Cid2_, Cid3_});
      SendRecvs.push_back(Fid_);
      // {T(0.5625), T(0.1875), T(0.1875), T(0.0625)});
    }
  }
}

template <typename T>
void BlockGeometry2D<T>::BuildBlockIndexMap() {
  _BlockIndexMap.clear();
  std::size_t count{};
  for (const Block2D<T> &block : _Blocks) {
    _BlockIndexMap[block.getBlockId()] = count;
    ++count;
  }
}

// BlockGeometryHelper2D

template <typename T>
BlockGeometryHelper2D<T>::BlockGeometryHelper2D(int Nx, int Ny, const AABB<T, 2> &AABBs,
                                                T voxelSize, int blockcelllen,
                                                std::uint8_t llimit, int overlap)
    : BasicBlock<T, 2>(voxelSize, AABBs.getExtended(Vector<T, 2>{voxelSize * overlap}),
                       AABB<int, 2>(Vector<int, 2>{0}, Vector<int, 2>{Nx - 1 + 2*overlap, Ny - 1 + 2*overlap})),
      _BaseBlock(voxelSize, AABBs,
                 AABB<int, 2>(Vector<int, 2>{overlap}, Vector<int, 2>{Nx - 1 + overlap, Ny - 1 + overlap})),
      BlockCellLen(blockcelllen), _Overlap(overlap), _LevelLimit(llimit), _MaxLevel(std::uint8_t(0)), 
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
BlockGeometryHelper2D<T>::BlockGeometryHelper2D(const BlockReader<T,2>& blockreader, bool useReaderOlap) 
    : BasicBlock<T, 2>(blockreader.getBasicBlock()), _BaseBlock(blockreader.getBaseBlock()), 
      BlockCellLen(0), _Overlap(1), _LevelLimit(std::uint8_t{}), _MaxLevel(blockreader.getMaxLevel()),
      _Exchanged(true), _IndexExchanged(true) {
  // create new blocks on relatively older blocks
  std::vector<BasicBlock<T, 2>> &BasicBlocks = getAllOldBasicBlocks();
  BasicBlocks.clear();
  // now old blocks become new blocks
  _Exchanged = !_Exchanged;
  // create blocks from Block Reader
  int iblock{};
  for (const BasicBlock<T, 2> &baseblock : blockreader.getBlocks()) {
    int overlap{};
    if (useReaderOlap) {
      overlap = blockreader.getOverlaps()[iblock];
    } else {
      overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
    }
    BasicBlocks.emplace_back(baseblock, overlap);
    ++iblock;
  }
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
  // collect info
  std::vector<std::size_t> BlockCellNumVec;
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

      BlockCellNumVec.push_back(Nx * Ny);
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

  // print info
  // find min and max block cell num
  std::size_t min = BlockCellNumVec[0];
  std::size_t max = BlockCellNumVec[0];
  for (std::size_t x : BlockCellNumVec) {
    min = x < min ? x : min;
    max = x > max ? x : max;
  }
  T aver = T(CellsN) / BasicBlocks.size();

  MPI_RANK(0)
  std::cout << "[BlockGeometryHelper2D<T>::CreateBlocks]: " << BasicBlocks.size() << " Blocks created, with: \n"
            << "  Min BlockCell Num: " << min << ", Max BlockCell Num: " << max << ", Average Block Cell Num: " << aver
            << std::endl;
}

template <typename T>
void BlockGeometryHelper2D<T>::CreateBlocks(int blockXNum, int blockYNum) {
  // create new blocks on relatively older blocks
  std::vector<BasicBlock<T, 2>> &BasicBlocks = getAllOldBasicBlocks();
  BasicBlocks.clear();
  // now old blocks become new blocks
  _Exchanged = !_Exchanged;
  // divide blocks manually
  DivideBlock2D(_BaseBlock, blockXNum, blockYNum, BasicBlocks);
}

template <typename T>
template <typename Func>
void BlockGeometryHelper2D<T>::forEachBlockCell(const Func& func) {
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
    StdDevs.push_back(ComputeBlockNStdDev(Blocks));
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
  // possible: a greedy algorithm
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
            // like SetupNbrs() in BlockGeometry2D, we use baseblock here, NO need to use extblock
            const BasicBlock<T, 2> nbaseblock = getAllBasicBlock(static_cast<std::size_t>(nblockid));
            if (isOverlapped(block, nbaseblock))
              nbrsvec.push_back(std::make_pair(nRank, nblockid));
          }
        }
      }
    }
  }
}

template <typename T>
void BlockGeometryHelper2D<T>::InitBlockGeometry2D() {
  _BlockGeometry2D = std::make_unique<BlockGeometry2D<T>>(*this, getAllBasicBlocks());
}

template <typename T>
int BlockGeometryHelper2D<T>::whichRank(int blockid) {
  for (std::size_t iRank = 0; iRank < getAllBlockIndices().size(); ++iRank) {
    for (int id : getAllBlockIndices()[iRank]) {
      if (id == blockid) return static_cast<int>(iRank);
    }
  }
  return -1;
}

#endif

template <typename T>
void BlockGeometryHelper2D<T>::TagRefineLayer(std::vector<std::uint8_t> &refine, bool &refined) {
  UpdateMaxLevel();
  // refine one additional layer if has neighbor with lower level
  // tag cells to be refined
  for (std::uint8_t level = _MaxLevel; level > std::uint8_t(0); --level) {
    for (int celly = 1; celly < CellsNy - 1; ++celly) {
      for (int cellx = 1; cellx < CellsNx - 1; ++cellx) {
        int cellid = celly * CellsNx + cellx;
        if (static_cast<bool>(refine[cellid])) {
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