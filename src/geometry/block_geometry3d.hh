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

// block_geometry3d.hh

#include "geometry/block_geometry3d.h"

template <typename T>
Block3D<T>::Block3D(const BasicBlock<T, 3> &baseblock, int olap)
    : BasicBlock<T, 3>(baseblock.getExtBlock(olap)), 
    _BaseBlock(baseblock), _overlap(olap) {}


template <typename T>
Block3D<T>::Block3D(const AABB<T, 3> &block, const AABB<int, 3> &idxblock, int blockid,
                    T voxelSize, int olap)
    : BasicBlock<T, 3>(voxelSize, block.getExtended(Vector<T, 3>{voxelSize*olap}),
                       idxblock.getExtended(Vector<int, 3>{olap}), blockid),
      _BaseBlock(voxelSize, block, idxblock, blockid), _overlap(olap) {
  // read from AABBs
  // ReadAABBs(AABBs, AABBflag);
}

template <typename T>
template <typename FieldType, typename LatSet>
void Block3D<T>::CleanLonelyFlags(FieldType& field, std::uint8_t flag, std::uint8_t voidflag, 
  unsigned int lonelyth, bool& cleaned, std::size_t& count) {
  GenericArray<bool> TransFlag(BasicBlock<T, 3>::N, false);
  const int overlap = _overlap;
  for (int z = overlap; z < BasicBlock<T, 3>::Mesh[2] - overlap; ++z) {
    for (int y = overlap; y < BasicBlock<T, 3>::Mesh[1] - overlap; ++y) {
      std::size_t id = BasicBlock<T, 3>::getIndex(Vector<int, 3>{overlap, y, z});
      for (int x = overlap; x < BasicBlock<T, 3>::Mesh[0] - overlap; ++x) {
        if (util::isFlag(field.get(id), flag)) {
          unsigned int lonely{};
          for (unsigned int i = 1; i < LatSet::q; ++i) {
            // get neighbor cell index
            const std::size_t nbridx = id + latset::c<LatSet>(i) * BasicBlock<T, 3>::getProjection();
            if (util::isFlag(field.get(nbridx), flag)) ++lonely;
          }
          if (lonely < lonelyth) {
            TransFlag[id] = true;
            cleaned = true;
            ++count;
          }
        }
        ++id;
      }
    }
  }
  
  if (cleaned) {
    for (std::size_t id = 0; id < BasicBlock<T, 3>::N; ++id) {
      if (TransFlag[id]) field.SetField(id, voidflag);
    }
  }
}

template <typename T>
template <typename FieldType, typename LatSet>
void Block3D<T>::SetupBoundary(const AABB<T, 3> &block, FieldType &field,
                               typename FieldType::value_type bdvalue) {
  // temp flag field store the transition flag
  GenericArray<bool> TransFlag(BasicBlock<T, 3>::N, false);
  const int overlap = _overlap;
  for (int z = overlap; z < BasicBlock<T, 3>::Mesh[2] - overlap; ++z) {
    for (int y = overlap; y < BasicBlock<T, 3>::Mesh[1] - overlap; ++y) {
      for (int x = overlap; x < BasicBlock<T, 3>::Mesh[0] - overlap; ++x) {
        const Vector<int, 3> locidx{x, y, z};
        const Vector<T, 3> vox = BasicBlock<T, 3>::getVoxel(locidx);
        for (unsigned int i = 1; i < LatSet::q; ++i) {
          const Vector<T, 3> nvox = vox + latset::c<LatSet>(i) * BasicBlock<T, 3>::VoxelSize;
          if (!block.isInside(nvox)) {
            TransFlag.set(BasicBlock<T, 3>::getIndex(locidx), true);
            break;
          }
        }
      }
    }
  }

  for (std::size_t id = 0; id < BasicBlock<T, 3>::N; ++id) {
    if (TransFlag[id]) field.SetField(id, bdvalue);
  }
}

template <typename T>
template <typename FieldType, typename LatSet>
void Block3D<T>::SetupBoundary(FieldType &field, typename FieldType::value_type fromvalue,
typename FieldType::value_type voidvalue, typename FieldType::value_type bdvalue) {
  // temp flag field store the transition flag
  GenericArray<bool> TransFlag(BasicBlock<T, 3>::N, false);
  const int overlap = _overlap;
  for (int z = overlap; z < BasicBlock<T, 3>::Mesh[2] - overlap; ++z) {
    for (int y = overlap; y < BasicBlock<T, 3>::Mesh[1] - overlap; ++y) {
      for (int x = overlap; x < BasicBlock<T, 3>::Mesh[0] - overlap; ++x) {
        const std::size_t idx = x + y * BasicBlock<T, 3>::Projection[1] + z * BasicBlock<T, 3>::Projection[2];
        if (field.get(idx) == fromvalue) {
          for (unsigned int i = 1; i < LatSet::q; ++i) {
            const std::size_t nbridx = idx + latset::c<LatSet>(i) * BasicBlock<T, 3>::Projection;
            if (field.get(nbridx) == voidvalue) {
              TransFlag.set(idx, true);
              break;
            }
          }
        }
      }
    }
  }

  for (std::size_t id = 0; id < BasicBlock<T, 3>::N; ++id) {
    if (TransFlag[id]) field.SetField(id, bdvalue);
  }
}

template <typename T>
template <typename FieldType>
void Block3D<T>::ReadOctree(Octree<T>* tree, FieldType& field, typename FieldType::value_type stlflag) {
  for (int z = 0; z < BasicBlock<T, 3>::Mesh[2]; ++z) {
    for (int y = 0; y < BasicBlock<T, 3>::Mesh[1]; ++y) {
      for (int x = 0; x < BasicBlock<T, 3>::Mesh[0]; ++x) {
        const Vector<int, 3> locidx{x, y, z};
        const Vector<T, 3> vox = BasicBlock<T, 3>::getVoxel(locidx);
        // get the node containing the voxel
        Octree<T>* node = tree->find(vox);
        if (node != nullptr) {
          // check if it is a [leaf] node and if it is [inside]
          if (node->isLeaf() && node->getInside())
            field.SetField(BasicBlock<T, 3>::getIndex(locidx), stlflag);
        }
      }
    }
  }
}

// -----------blockgeometry3d----------------


template <typename T>
BlockGeometry3D<T>::BlockGeometry3D(int Nx, int Ny, int Nz, int blocknum,
                                    const AABB<T, 3> &block, T voxelSize, int overlap,
                                    int blockXNum, int blockYNum, int blockZNum)
    : BasicBlock<T, 3>(
        voxelSize, block.getExtended(Vector<T, 3>{voxelSize*overlap}),
        AABB<int, 3>(Vector<int, 3>{0}, Vector<int, 3>{Nx - 1 + 2*overlap, Ny - 1 + 2*overlap, Nz - 1 + 2*overlap})),
      _BaseBlock(voxelSize, block,
                 AABB<int, 3>(Vector<int, 3>{overlap}, Vector<int, 3>{Nx - 1 + overlap, Ny - 1 + overlap, Nz - 1 + overlap})),
      _overlap(overlap), _MaxLevel(std::uint8_t(0)) {
  CreateBlocks(blocknum, blockXNum, blockYNum, blockZNum);
  BuildBlockIndexMap();
  SetupNbrs();
  InitComm();
  PrintInfo();
}

template <typename T>
BlockGeometry3D<T>::BlockGeometry3D(BlockGeometryHelper3D<T> &GeoHelper, bool useHelperOlap)
    : BasicBlock<T, 3>(GeoHelper), _BaseBlock(GeoHelper.getBaseBlock()), 
      _overlap(GeoHelper.getOverlap()), _MaxLevel(GeoHelper.getMaxLevel()) {
  // create blocks from GeoHelper
  for (BasicBlock<T, 3> *baseblock : GeoHelper.getBasicBlocks()) {
    int overlap{};
    if (useHelperOlap) {
      overlap = _overlap;
    } else {
      overlap = (baseblock->getLevel() != std::uint8_t(0)) ? 2 : 1;
    }
    _Blocks.emplace_back(*baseblock, overlap);
    _BasicBlocks.emplace_back(*baseblock);
  }
  BuildBlockIndexMap();
  SetupNbrs();
  InitAllComm();
#ifdef MPI_ENABLED
  GeoHelper.InitBlockGeometry3D(useHelperOlap);
  InitAllMPIComm(GeoHelper);
#endif
  PrintInfo();
}

template <typename T>
BlockGeometry3D<T>::BlockGeometry3D(const BlockReader<T,3>& blockreader, bool useReaderOlap) 
    : BasicBlock<T, 3>(blockreader.getBasicBlock()), _BaseBlock(blockreader.getBaseBlock()), 
      _overlap(1), _MaxLevel(blockreader.getMaxLevel()) {
  // create blocks from Block Reader
  int iblock{};
  for (const BasicBlock<T, 3> &baseblock : blockreader.getBlocks()) {
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
BlockGeometry3D<T>::BlockGeometry3D(const StlReader<T>& reader, int blocknum) 
    : BasicBlock<T, 3>(reader.getVoxelSize(), 
        AABB<T, 3>(reader.getMesh().getMin() - reader.getVoxelSize(), reader.getMesh().getMax() + reader.getVoxelSize()),
        AABB<int, 3>(Vector<int, 3>{0}, 
                     Vector<int, 3>{int(std::ceil(reader.getMesh().getMax_Min()[0] / reader.getVoxelSize())) + 1, 
                                    int(std::ceil(reader.getMesh().getMax_Min()[1] / reader.getVoxelSize())) + 1,
                                    int(std::ceil(reader.getMesh().getMax_Min()[2] / reader.getVoxelSize())) + 1})),
      _BaseBlock(reader.getVoxelSize(), 
        AABB<T, 3>(reader.getMesh().getMin(), reader.getMesh().getMax()),
        AABB<int, 3>(Vector<int, 3>{1}, 
                     Vector<int, 3>{int(std::ceil(reader.getMesh().getMax_Min()[0] / reader.getVoxelSize())), 
                                    int(std::ceil(reader.getMesh().getMax_Min()[1] / reader.getVoxelSize())),
                                    int(std::ceil(reader.getMesh().getMax_Min()[2] / reader.getVoxelSize()))})), 
      _overlap(1), _MaxLevel(std::uint8_t(0)) {
  CreateBlocks(blocknum);
  BuildBlockIndexMap();
  SetupNbrs();
  InitComm();
  PrintInfo();
}

template <typename T>
void BlockGeometry3D<T>::PrintInfo() const {
  std::size_t total_cellnum = getTotalCellNum();
  std::size_t base_cellnum = getBaseCellNum();
  // DO NOT use MPI_RANK(0) before getTotalCellNum() cause we have: mpi().barrier();
  MPI_RANK(0)
  std::cout << "[BlockGeometry3D]:\n" 
            << "  Total CellNum: " << total_cellnum << " \n"
            << "  Base CellNum: " << base_cellnum << std::endl;
#ifndef MPI_ENABLED
  std::cout << "  BlockNum: " << _Blocks.size() << " with StdDev: "
            << ComputeBlockNStdDev(_BasicBlocks) << std::endl;
#endif
}

template <typename T>
void BlockGeometry3D<T>::Init(BlockGeometryHelper3D<T> &GeoHelper) {
  ;
  _Blocks.clear();
  _BasicBlocks.clear();
  // create blocks from GeoHelper
  for (BasicBlock<T, 3> *baseblock : GeoHelper.getBasicBlocks()) {
    int overlap = (baseblock->getLevel() != std::uint8_t(0)) ? 2 : 1;
    _Blocks.emplace_back(*baseblock, overlap);
    _BasicBlocks.emplace_back(*baseblock);
  }
  BuildBlockIndexMap();
  SetupNbrs();
  InitAllComm();
#ifdef MPI_ENABLED
// GeoHelper.InitBlockGeometry3D();
  InitAllMPIComm(GeoHelper);
#endif
  PrintInfo();
}


template <typename T>
template <typename FieldType>
void BlockGeometry3D<T>::InitVoxelMap(FieldType& FieldM, std::uint8_t voidflag) {
  #ifdef _VOX_ENABLED
  std::size_t TotalRemoved{};
  int iblock{};
  for (Block3D<T>& block : _Blocks) {
    auto& field = FieldM.getBlockField(iblock).getFieldType().getField(0);
    std::size_t N = block.getN();
    std::size_t voxNum{};
    std::vector<std::size_t> vox_to_aabb;
    vox_to_aabb.reserve(N);

    for (std::size_t id = 0; id < N; ++id) {
      if (field[id] != voidflag) {
        vox_to_aabb.emplace_back(id);
        ++voxNum;
      }
    }

    block.getVoxelMap().Init(voxNum);
    std::size_t* vox_aabb = block.getVoxelMap().getVox_to_AABB();
    for (std::size_t idx = 0; idx < voxNum; ++idx) {
      vox_aabb[idx] = vox_to_aabb[idx];
    }

    // remove void cells from communicator
    // remove void cells from communicator
    Communicator& COMM = block.getCommunicator();
    // comm
    std::vector<SharedComm>& comms = COMM.Comm.Comms;
    for (SharedComm& comm : comms) {
      // because indices are ordered as: Scell0, Rcell0; Scell1, Rcell1, ...
      // use iterator seems complex, we create a new vector instead
      std::vector<std::size_t> newSendRecvs;
      newSendRecvs.reserve(comm.SendRecvCells.size());
      for (std::size_t i = 0; i < comm.SendRecvCells.size(); i += 2) {
        if (field[comm.SendRecvCells[i+1]] != voidflag) {
          newSendRecvs.emplace_back(comm.SendRecvCells[i]);
          newSendRecvs.emplace_back(comm.SendRecvCells[i+1]);
        }
      }
      if (newSendRecvs.size() != comm.SendRecvCells.size()) {
        comm.SendRecvCells = newSendRecvs;
        comm.SendRecvCells.shrink_to_fit();
      }
    }
    // all comm
    std::vector<SharedComm>& allcomms = COMM.AllComm.Comms;
    for (SharedComm& comm : allcomms) {
      // because indices are ordered as: Scell0, Rcell0; Scell1, Rcell1, ...
      // use iterator seems complex, we create a new vector instead
      std::vector<std::size_t> newSendRecvs;
      newSendRecvs.reserve(comm.SendRecvCells.size());
      for (std::size_t i = 0; i < comm.SendRecvCells.size(); i += 2) {
        if (field[comm.SendRecvCells[i+1]] != voidflag) {
          newSendRecvs.emplace_back(comm.SendRecvCells[i]);
          newSendRecvs.emplace_back(comm.SendRecvCells[i+1]);
        }
      }
      if (newSendRecvs.size() != comm.SendRecvCells.size()) {
        TotalRemoved += (comm.SendRecvCells.size() - newSendRecvs.size()) / 2;
        comm.SendRecvCells = newSendRecvs;
        comm.SendRecvCells.shrink_to_fit();
      }
    }

  #ifdef MPI_ENABLED

    std::vector<DistributedComm>& Recvs = COMM.MPIComm.Recvs;
    for (DistributedComm& comm : Recvs) {
      std::vector<std::size_t> newcells;
      newcells.reserve(comm.Cells.size());
      for (std::size_t id : comm.Cells) {
        if (field[id] != voidflag) newcells.emplace_back(id);
      }
      if (newcells.size() != comm.Cells.size()) {
        comm.Cells = newcells;
        comm.Cells.shrink_to_fit();
      }
    }

    std::vector<DistributedComm>& Sends = COMM.MPIComm.Sends;
    for (DistributedComm& comm : Sends) {
      std::vector<std::size_t> newcells;
      newcells.reserve(comm.Cells.size());
      for (std::size_t id : comm.Cells) {
        if (field[id] != voidflag) newcells.emplace_back(id);
      }
      if (newcells.size() != comm.Cells.size()) {
        comm.Cells = newcells;
        comm.Cells.shrink_to_fit();
      }
    }

    std::vector<DistributedComm>& allRecvs = COMM.AllMPIComm.Recvs;
    for (DistributedComm& comm : allRecvs) {
      std::vector<std::size_t> newcells;
      newcells.reserve(comm.Cells.size());
      for (std::size_t id : comm.Cells) {
        if (field[id] != voidflag) newcells.emplace_back(id);
      }
      if (newcells.size() != comm.Cells.size()) {
        TotalRemoved += comm.Cells.size() - newcells.size();
        comm.Cells = newcells;
        comm.Cells.shrink_to_fit();
      }
    }

    std::vector<DistributedComm>& allSends = COMM.AllMPIComm.Sends;
    for (DistributedComm& comm : allSends) {
      std::vector<std::size_t> newcells;
      newcells.reserve(comm.Cells.size());
      for (std::size_t id : comm.Cells) {
        if (field[id] != voidflag) newcells.emplace_back(id);
      }
      if (newcells.size() != comm.Cells.size()) {
        comm.Cells = newcells;
        comm.Cells.shrink_to_fit();
      }
    }

    std::vector<DistributedComm>& dirRecvs = COMM.DirRecvs;
    for (DistributedComm& comm : dirRecvs) {
      std::vector<std::size_t> newcells;
      newcells.reserve(comm.Cells.size());
      for (std::size_t id : comm.Cells) {
        if (field[id] != voidflag) newcells.emplace_back(id);
      }
      if (newcells.size() != comm.Cells.size()) {
        comm.Cells = newcells;
        comm.Cells.shrink_to_fit();
      }
    }

    std::vector<DistributedComm>& dirSends = COMM.DirSends;
    for (DistributedComm& comm : dirSends) {
      std::vector<std::size_t> newcells;
      newcells.reserve(comm.Cells.size());
      for (std::size_t id : comm.Cells) {
        if (field[id] != voidflag) newcells.emplace_back(id);
      }
      if (newcells.size() != comm.Cells.size()) {
        comm.Cells = newcells;
        comm.Cells.shrink_to_fit();
      }
    }

  #endif

    ++iblock;
  }

  #ifdef MPI_ENABLED
		mpi().reduceAndBcast(TotalRemoved, MPI_MAX);
  #endif

  MPI_RANK(0)
  std::cout << "[BlockGeometry::InitVoxelMap]:\n"
              << "  Removed: " << TotalRemoved << " cell pairs in Communicator" << std::endl;
  
  #endif
}


template <typename T>
std::size_t BlockGeometry3D<T>::getTotalCellNum() const {
  std::size_t sum{};
  for (const Block3D<T> &block : _Blocks) sum += block.getN();
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
std::size_t BlockGeometry3D<T>::getBaseCellNum() const {
  std::size_t sum{};
  for (const Block3D<T> &block : _Blocks) sum += block.getBaseBlock().getN();
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
void BlockGeometry3D<T>::DivideBlocks(int blocknum, int blockXNum, int blockYNum, int blockZNum) {
  _BlockAABBs.clear();
  _BlockAABBs.reserve(blocknum);
  if (blockXNum == 0 || blockYNum == 0 || blockZNum == 0) {
    // default scheme
    DivideBlock3D(_BaseBlock, blocknum, _BlockAABBs);
  } else {
    // manually divide
    DivideBlock3D(_BaseBlock, blockXNum, blockYNum, blockZNum, _BlockAABBs);
  }
}

template <typename T>
void BlockGeometry3D<T>::CreateBlocks(int blocknum, int blockXNum, int blockYNum, int blockZNum) {
  DivideBlocks(blocknum);

  _Blocks.clear();
  _Blocks.reserve(_BlockAABBs.size());
  _BasicBlocks.clear();
  _BasicBlocks.reserve(_BlockAABBs.size());
  // create blocks
  int blockid = 0;
  const T voxsize = _BaseBlock.getVoxelSize();
  for (const AABB<int, 3> &blockaabb : _BlockAABBs) {
    Vector<T, 3> MIN =
      (blockaabb.getMin() - Vector<int,3>{_overlap}) * voxsize + _BaseBlock.getMin();
    Vector<T, 3> MAX =
      (blockaabb.getMax() - Vector<int,3>{_overlap-1}) * voxsize + _BaseBlock.getMin();
    AABB<T, 3> aabb(MIN, MAX);
    _Blocks.emplace_back(aabb, blockaabb, blockid, voxsize, _overlap);
    _BasicBlocks.emplace_back(voxsize, aabb, blockaabb, blockid);
    blockid++;
  }
}

template <typename T>
void BlockGeometry3D<T>::SetupNbrs() {
  for (Block3D<T> &block : _Blocks) {
    int id = block.getBlockId();
    std::vector<Block3D<T> *> &nbrsvec = block.getNeighbors();
    nbrsvec.clear();
    for (Block3D<T> &blockn : _Blocks) {
      int idn = blockn.getBlockId();
      if (id != idn) {
        if (isOverlapped(block.getBaseBlock().getExtBlock(1), blockn.getBaseBlock())) {
          nbrsvec.push_back(&blockn);
        }
      }
    }
  }
}

// comms

template <typename T>
void BlockGeometry3D<T>::InitComm() {
  int Tag{};
  for (Block3D<T> &block : _Blocks) {
    // (inner) overlapped cell communicators
    std::vector<SharedComm> &Comms = block.getCommunicator().Comm.Comms;
    // all overlapped cell communicators
    std::vector<SharedComm> &AllComms = block.getCommunicator().AllComm.Comms;
    Comms.clear();
    AllComms.clear();
    // get block with overlap 1, for Comms
    BasicBlock<T, 3> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
    // get block with overlap = _overlap, for AllComms
    BasicBlock<T, 3> baseblock_exto = block.getSelfBlock();

    std::uint8_t blocklevel = block.getLevel();
    for (Block3D<T> *nblock : block.getNeighbors()) {
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
        // exclude edge cells
        std::vector<std::vector<std::size_t>> EdgeRecvs;
        std::vector<std::vector<std::size_t>> EdgeSends;
        block.ExcludeEdgeIdx(Recvs, Sends, EdgeRecvs, EdgeSends);
        // avoid empty communicator
        if (Recvs.size() > 0) {
          Comms.emplace_back(nblock->getBlockId(), Tag);
          ++Tag;
          SharedComm &comm = Comms.back();
          comm.setRecvSendIdx(Recvs, Sends);
          comm.Direction = getFaceNbrDirection(block.whichFace(comm.SendRecvCells[1]));
        }
        if (CornerRecvs.size() > 0) {
          // add corner cells to communicators and find direction
          for (std::size_t i = 0; i < CornerRecvs.size(); ++i) {
            Comms.emplace_back(nblock->getBlockId(), Tag);
            ++Tag;
            SharedComm &commcorner = Comms.back();
            commcorner.SendRecvCells = {CornerSends[i], CornerRecvs[i]};
            commcorner.Direction = getCornerNbrDirection<3>(block.whichCorner(CornerRecvs[i]));
          }
        }
        // add edge cells to communicators and find direction
        for (std::size_t i = 0; i < EdgeRecvs.size(); ++i) {
          if (EdgeRecvs[i].size() > 0) {
            Comms.emplace_back(nblock->getBlockId(), Tag);
            ++Tag;
            SharedComm &commedge = Comms.back();
            commedge.setRecvSendIdx(EdgeRecvs[i], EdgeSends[i]);
            commedge.Direction = getEdgeNbrDirection<3>(block.whichEdge(EdgeRecvs[i][0]));
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
void BlockGeometry3D<T>::InitAverComm() {
  for (Block3D<T> &block : _Blocks) {
    // (inner) overlapped cell communicators
    std::vector<SharedComm> &Comms = block.getCommunicator().Comm.AverComm;
    // all overlapped cell communicators
    std::vector<SharedComm> &AllComms = block.getCommunicator().AllComm.AverComm;
    Comms.clear();
    AllComms.clear();
    // get block with overlap 1, for Comms
    BasicBlock<T, 3> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
    // get block with overlap = _overlap, for AllComms
    BasicBlock<T, 3> baseblock_exto = block.getSelfBlock();

    std::uint8_t blocklevel = block.getLevel();
    for (Block3D<T> *nblock : block.getNeighbors()) {
      // find block of blocklevel+1
      if (nblock->getLevel() == blocklevel + 1) {
        // ------ add to Comms
        AddtoSharedAverComm(Comms, baseblock_ext1, block, nblock);
        // ------ add to AllComms
        AddtoSharedAverComm(AllComms, baseblock_exto, block, nblock);
      } else if (nblock->getLevel() > blocklevel + 1) {
        std::cerr << "[BlockGeometry3D<T>::InitAverComm] Error: block level difference "
                     "larger than 1"
                  << std::endl;
      }
    }
  }
}

template <typename T>
void BlockGeometry3D<T>::InitIntpComm() {
  for (Block3D<T> &block : _Blocks) {
    // (inner) overlapped cell communicators
    std::vector<SharedComm> &Comms = block.getCommunicator().Comm.IntpComm;
    // all overlapped cell communicators
    std::vector<SharedComm> &AllComms = block.getCommunicator().AllComm.IntpComm;
    Comms.clear();
    AllComms.clear();
    // get block with overlap = _overlap, for both Comms and AllComms
    std::uint8_t blocklevel = block.getLevel();
    for (Block3D<T> *nblock : block.getNeighbors()) {
      // find block of blocklevel-1
      if (nblock->getLevel() == blocklevel - 1) {
        // ------ add to Comms
        AddtoSharedIntpComm(Comms, block, nblock);
        // ------ add to AllComms, for now it is the same as Comms
        AddtoSharedIntpComm(AllComms, block, nblock);
      } else if (nblock->getLevel() < blocklevel - 1) {
        std::cerr << "[BlockGeometry3D<T>::InitIntpComm] Error: block level difference "
                     "larger than 1"
                  << std::endl;
      }
    }
  }
}
template <typename T>
void BlockGeometry3D<T>::InitAllComm() {
  InitComm();
  InitAverComm();
  InitIntpComm();
}

#ifdef MPI_ENABLED

template <typename T>
BlockGeometry3D<T>::BlockGeometry3D(BlockGeometryHelper3D<T> &GeoHelper, 
std::vector<BasicBlock<T, 3>>& BasicBlocks, bool useHelperOlap)
    : BasicBlock<T, 3>(GeoHelper), _BaseBlock(GeoHelper.getBaseBlock()), 
      _overlap(GeoHelper.getOverlap()), _MaxLevel(GeoHelper.getMaxLevel()) {
  // create blocks from GeoHelper
  for (BasicBlock<T, 3>& baseblock : BasicBlocks) {
    int overlap{};
    if (useHelperOlap) {
      overlap = _overlap;
    } else {
      overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
    }
    _Blocks.emplace_back(baseblock, overlap);
    _BasicBlocks.emplace_back(baseblock);
  }
  BuildBlockIndexMap();
  SetupNbrs();
  // this uses all blocks of the same level to init shared communicators
  // for mpi communication with direction info 
  InitComm();
}

template <typename T>
void BlockGeometry3D<T>::InitMPIComm(BlockGeometryHelper3D<T> &GeoHelper) {
  // mpi efficient send/recv using direction info for pop communication
  const BlockGeometry3D<T>& HelperBlockGeometry = GeoHelper.getBlockGeometry3D();
  for (Block3D<T> &block : _Blocks) {
    // efficient MPI Recvs/Sends using direction info for pop communication
    std::vector<DistributedComm>& Recvs = block.getCommunicator().DirRecvs;
    std::vector<DistributedComm>& Sends = block.getCommunicator().DirSends;
    Recvs.clear();
    Sends.clear();
    // find hblock in HelperBlockGeometry with:
    // the same blockid as block or
    // the same sendblockid in hblock's communicator
    for (const Block3D<T> &hblock : HelperBlockGeometry.getBlocks()) {
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

for (Block3D<T> &block : _Blocks) {
    // (inner) overlapped cell communicators
    DistributedCommSet &MPIComm = block.getCommunicator().MPIComm;
    // all overlapped cell communicators
    DistributedCommSet &AllMPIComm = block.getCommunicator().AllMPIComm;
    MPIComm.Recvs.clear();
    MPIComm.Sends.clear();
    AllMPIComm.Recvs.clear();
    AllMPIComm.Sends.clear();
    // base block of block
    const BasicBlock<T, 3> &baseblock = block.getBaseBlock();
    // get block with overlap 1, for Comms
    BasicBlock<T, 3> baseblock_ext1 = baseblock.getExtBlock(1);
    // get block with overlap = _overlap, for AllComms
    BasicBlock<T, 3> baseblock_exto = block.getSelfBlock();

    std::uint8_t blocklevel = block.getLevel();
    std::vector<std::pair<int, int>> &nbrs = GeoHelper.getMPIBlockNbrs(block.getBlockId());

    for (const std::pair<int, int> &nbr : nbrs) {
      // check if 2 blocks are of the same level
      const BasicBlock<T, 3> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel) {
        BasicBlock<T, 3> nbaseblock_ext1 = nbaseblock.getExtBlock(1);
        BasicBlock<T, 3> nbaseblock_exto = nbaseblock.getExtBlock(block.getOverlap());
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
  // for (Block3D<T> &block : _Blocks) {
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
  // for (Block3D<T> &block : _Blocks) {
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
  // for (Block3D<T> &block : _Blocks) {
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
void BlockGeometry3D<T>::InitMPIAverComm(BlockGeometryHelper3D<T> &GeoHelper) {
  for (Block3D<T> &block : _Blocks) {
    // (inner) overlapped cell communicators
    DistributedCommSet &MPIComm = block.getCommunicator().MPIComm;
    // all overlapped cell communicators
    DistributedCommSet &AllMPIComm = block.getCommunicator().AllMPIComm;
    MPIComm.AverRecvs.clear();
    MPIComm.AverSends.clear();
    AllMPIComm.AverRecvs.clear();
    AllMPIComm.AverSends.clear();

    // base block of block
    const BasicBlock<T, 3> &baseblock = block.getBaseBlock();
    // get block with overlap 1, for Comms
    BasicBlock<T, 3> baseblock_ext1 = baseblock.getExtBlock(1);
    // get block with overlap = _overlap, for AllComms
    BasicBlock<T, 3> baseblock_exto = block.getSelfBlock();

    std::uint8_t blocklevel = block.getLevel();
    std::vector<std::pair<int, int>> &nbrs = GeoHelper.getMPIBlockNbrs(block.getBlockId());

    for (const std::pair<int, int> &nbr : nbrs) {
      const BasicBlock<T, 3> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
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
        BasicBlock<T, 3> Cblock = block.getCoasenedBlock();

        BasicBlock<T, 3> nbaseblock_ext1 = nbaseblock.getExtBlock(1);
        int noverlap = nbaseblock.getLevel() == std::uint8_t(0) ? 1 : 2;
        BasicBlock<T, 3> nbaseblock_exto = nbaseblock.getExtBlock(noverlap);

        MPIComm.AverSends.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &sender = MPIComm.AverSends.back();
        std::vector<std::size_t> VSends;
        Cblock.getCellIdx(baseblock, nbaseblock_ext1, VSends);
        sender.Cells.reserve(VSends.size()*8);
        for (std::size_t vid : VSends) {
          std::vector<std::size_t> idxs;
          Cblock.getRefinedCellIdx(vid, idxs);
          sender.Cells.insert(sender.Cells.end(), idxs.begin(), idxs.end());          
        }

        AllMPIComm.AverSends.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &allsender = AllMPIComm.AverSends.back();
        std::vector<std::size_t> AllVSends;
        Cblock.getCellIdx(baseblock, nbaseblock_exto, AllVSends);
        allsender.Cells.reserve(AllVSends.size()*8);
        for (std::size_t vid : AllVSends) {
          std::vector<std::size_t> idxs;
          Cblock.getRefinedCellIdx(vid, idxs);
          allsender.Cells.insert(allsender.Cells.end(), idxs.begin(), idxs.end());          
        }

        block.getCommunicator()._NeedMPIComm = true;
      } else if (nbaseblock.getLevel() > blocklevel + 1 ||
                 nbaseblock.getLevel() < blocklevel - 1) {
        std::cerr << "[BlockGeometry3D<T>::InitMPIAverComm] Error: block level "
                     "difference larger than 1"
                  << std::endl;
      }
    }
  }
}

template <typename T>
void BlockGeometry3D<T>::InitMPIIntpComm(BlockGeometryHelper3D<T> &GeoHelper) {
  for (Block3D<T> &block : _Blocks) {
    // (inner) overlapped cell communicators
    DistributedCommSet &MPIComm = block.getCommunicator().MPIComm;
    // all overlapped cell communicators
    DistributedCommSet &AllMPIComm = block.getCommunicator().AllMPIComm;
    MPIComm.IntpRecvs.clear();
    MPIComm.IntpSends.clear();
    AllMPIComm.IntpRecvs.clear();
    AllMPIComm.IntpSends.clear();

    // base block of block
    const BasicBlock<T, 3> &baseblock = block.getBaseBlock();
    // get block with overlap 1, for Comms
    // BasicBlock<T, 3> baseblock_ext1 = baseblock.getExtBlock(1);
    // get block with overlap = _overlap, for AllComms
    // BasicBlock<T, 3> baseblock_exto = block.getSelfBlock();

    std::uint8_t blocklevel = block.getLevel();
    std::vector<std::pair<int, int>> &nbrs = GeoHelper.getMPIBlockNbrs(block.getBlockId());

    for (const std::pair<int, int> &nbr : nbrs) {
      const BasicBlock<T, 3> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel - 1) {
        // init receiver
        // virtual coarse block
        BasicBlock<T, 3> Cblock = block.getCoasenedBlock();

        MPIComm.IntpRecvs.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &recver = MPIComm.IntpRecvs.back();
        AllMPIComm.IntpRecvs.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &allrecver = AllMPIComm.IntpRecvs.back();

        std::vector<std::size_t> VRecvs;
        Cblock.getCellIdx(block, nbaseblock, VRecvs);

        recver.Cells.reserve(VRecvs.size()*8);
        allrecver.Cells.reserve(VRecvs.size()*8);
        for (std::size_t vid : VRecvs) {
          std::vector<std::size_t> idxs;
          Cblock.getRefinedCellIdx(vid, idxs);
          recver.Cells.insert(recver.Cells.end(), idxs.begin(), idxs.end());
          allrecver.Cells.insert(allrecver.Cells.end(), idxs.begin(), idxs.end());        
        }

        block.getCommunicator()._NeedMPIComm = true;
      } else if (nbaseblock.getLevel() == blocklevel + 1) {
        // init sender
        BasicBlock<T, 3> nbaseblock_ext2 = nbaseblock.getExtBlock(2);
        
        MPIComm.IntpSends.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &sender = MPIComm.IntpSends.back();
        AllMPIComm.IntpSends.emplace_back(nbr.first, nbaseblock.getBlockId());
        DistributedComm &allsender = AllMPIComm.IntpSends.back();

        // vox size
        const T Cvoxsize = block.getVoxelSize();
        // get intersection
        const AABB<T, 3> intsec = getIntersection(baseblock, nbaseblock_ext2);
        // use coarse grid size here, convient for calculating coarse cell index
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        int CNz = static_cast<int>(std::round(intsec.getExtension()[2] / Cvoxsize));
        // start index of intsec in block
        Vector<T, 3> startC = intsec.getMin() - block.getMin();
        // shift 1 voxel to left bottom for interpolation
        int startCx = static_cast<int>(std::round(startC[0] / Cvoxsize)) - 1;
        int startCy = static_cast<int>(std::round(startC[1] / Cvoxsize)) - 1;
        int startCz = static_cast<int>(std::round(startC[2] / Cvoxsize)) - 1;

        std::vector<std::size_t>& Sends = sender.Cells;
        std::vector<std::size_t>& AllSends = allsender.Cells;

        Sends.reserve(CNx * CNy * CNz * 64);
        AllSends.reserve(CNx * CNy * CNz * 64);

        std::size_t XY = block.getNx() * block.getNy();
        for (int iz = 0; iz < CNz; ++iz) {
          for (int iy = 0; iy < CNy; ++iy) {
            for (int ix = 0; ix < CNx; ++ix) {
              // original
              std::size_t Cid0 = (iz + startCz) * XY + (iy + startCy) * block.getNx() + ix + startCx;
              std::size_t Cid1 = Cid0 + 1;
              std::size_t Cid2 = Cid0 + block.getNx();
              std::size_t Cid3 = Cid2 + 1;
              std::size_t Cid4 = Cid0 + XY;
              std::size_t Cid5 = Cid4 + 1;
              std::size_t Cid6 = Cid4 + block.getNx();
              std::size_t Cid7 = Cid6 + 1;

              // shift 1 voxel along +x direction
              std::size_t Cid0_x = Cid0 + 1;
              std::size_t Cid1_x = Cid1 + 1;
              std::size_t Cid2_x = Cid2 + 1;
              std::size_t Cid3_x = Cid3 + 1;
              std::size_t Cid4_x = Cid4 + 1;
              std::size_t Cid5_x = Cid5 + 1;
              std::size_t Cid6_x = Cid6 + 1;
              std::size_t Cid7_x = Cid7 + 1;

              // shift 1 voxel along +y direction
              std::size_t Cid0_y = Cid0 + block.getNx();
              std::size_t Cid1_y = Cid1 + block.getNx();
              std::size_t Cid2_y = Cid2 + block.getNx();
              std::size_t Cid3_y = Cid3 + block.getNx();
              std::size_t Cid4_y = Cid4 + block.getNx();
              std::size_t Cid5_y = Cid5 + block.getNx();
              std::size_t Cid6_y = Cid6 + block.getNx();
              std::size_t Cid7_y = Cid7 + block.getNx();

              // shift 1 voxel along +z direction
              std::size_t Cid0_z = Cid0 + XY;
              std::size_t Cid1_z = Cid1 + XY;
              std::size_t Cid2_z = Cid2 + XY;
              std::size_t Cid3_z = Cid3 + XY;
              std::size_t Cid4_z = Cid4 + XY;
              std::size_t Cid5_z = Cid5 + XY;
              std::size_t Cid6_z = Cid6 + XY;
              std::size_t Cid7_z = Cid7 + XY;

              // 0
              Sends.insert(Sends.end(), {Cid0, Cid1, Cid2, Cid3, Cid4, Cid5, Cid6, Cid7});
              AllSends.insert(AllSends.end(), {Cid0, Cid1, Cid2, Cid3, Cid4, Cid5, Cid6, Cid7});

              // 1
              Sends.insert(Sends.end(), {Cid0_x, Cid1_x, Cid2_x, Cid3_x, Cid4_x, Cid5_x, Cid6_x, Cid7_x});
              AllSends.insert(AllSends.end(), {Cid0_x, Cid1_x, Cid2_x, Cid3_x, Cid4_x, Cid5_x, Cid6_x, Cid7_x});

              // 2
              Sends.insert(Sends.end(), {Cid0_y, Cid1_y, Cid2_y, Cid3_y, Cid4_y, Cid5_y, Cid6_y, Cid7_y});
              AllSends.insert(AllSends.end(), {Cid0_y, Cid1_y, Cid2_y, Cid3_y, Cid4_y, Cid5_y, Cid6_y, Cid7_y});

              // 3
              std::size_t Cid0_xy = Cid0_y + 1;
              std::size_t Cid1_xy = Cid1_y + 1;
              std::size_t Cid2_xy = Cid2_y + 1;
              std::size_t Cid3_xy = Cid3_y + 1;
              std::size_t Cid4_xy = Cid4_y + 1;
              std::size_t Cid5_xy = Cid5_y + 1;
              std::size_t Cid6_xy = Cid6_y + 1;
              std::size_t Cid7_xy = Cid7_y + 1;

              Sends.insert(Sends.end(), {Cid0_xy, Cid1_xy, Cid2_xy, Cid3_xy, Cid4_xy, Cid5_xy, Cid6_xy, Cid7_xy});
              AllSends.insert(AllSends.end(), {Cid0_xy, Cid1_xy, Cid2_xy, Cid3_xy, Cid4_xy, Cid5_xy, Cid6_xy, Cid7_xy});

              // 4
              Sends.insert(Sends.end(), {Cid0_z, Cid1_z, Cid2_z, Cid3_z, Cid4_z, Cid5_z, Cid6_z, Cid7_z});
              AllSends.insert(AllSends.end(), {Cid0_z, Cid1_z, Cid2_z, Cid3_z, Cid4_z, Cid5_z, Cid6_z, Cid7_z});

              // 5
              std::size_t Cid0_xz = Cid0_z + 1;
              std::size_t Cid1_xz = Cid1_z + 1;
              std::size_t Cid2_xz = Cid2_z + 1;
              std::size_t Cid3_xz = Cid3_z + 1;
              std::size_t Cid4_xz = Cid4_z + 1;
              std::size_t Cid5_xz = Cid5_z + 1;
              std::size_t Cid6_xz = Cid6_z + 1;
              std::size_t Cid7_xz = Cid7_z + 1;

              Sends.insert(Sends.end(), {Cid0_xz, Cid1_xz, Cid2_xz, Cid3_xz, Cid4_xz, Cid5_xz, Cid6_xz, Cid7_xz});
              AllSends.insert(AllSends.end(), {Cid0_xz, Cid1_xz, Cid2_xz, Cid3_xz, Cid4_xz, Cid5_xz, Cid6_xz, Cid7_xz});

              // 6
              std::size_t Cid0_yz = Cid0_z + block.getNx();
              std::size_t Cid1_yz = Cid1_z + block.getNx();
              std::size_t Cid2_yz = Cid2_z + block.getNx();
              std::size_t Cid3_yz = Cid3_z + block.getNx();
              std::size_t Cid4_yz = Cid4_z + block.getNx();
              std::size_t Cid5_yz = Cid5_z + block.getNx();
              std::size_t Cid6_yz = Cid6_z + block.getNx();
              std::size_t Cid7_yz = Cid7_z + block.getNx();

              Sends.insert(Sends.end(), {Cid0_yz, Cid1_yz, Cid2_yz, Cid3_yz, Cid4_yz, Cid5_yz, Cid6_yz, Cid7_yz});
              AllSends.insert(AllSends.end(), {Cid0_yz, Cid1_yz, Cid2_yz, Cid3_yz, Cid4_yz, Cid5_yz, Cid6_yz, Cid7_yz});

              // 7
              std::size_t Cid0_xyz = Cid0_yz + 1;
              std::size_t Cid1_xyz = Cid1_yz + 1;
              std::size_t Cid2_xyz = Cid2_yz + 1;
              std::size_t Cid3_xyz = Cid3_yz + 1;
              std::size_t Cid4_xyz = Cid4_yz + 1;
              std::size_t Cid5_xyz = Cid5_yz + 1;
              std::size_t Cid6_xyz = Cid6_yz + 1;
              std::size_t Cid7_xyz = Cid7_yz + 1;

              Sends.insert(Sends.end(), {
                Cid0_xyz, Cid1_xyz, Cid2_xyz, Cid3_xyz, Cid4_xyz, Cid5_xyz, Cid6_xyz, Cid7_xyz});
              AllSends.insert(AllSends.end(), {
                Cid0_xyz, Cid1_xyz, Cid2_xyz, Cid3_xyz, Cid4_xyz, Cid5_xyz, Cid6_xyz, Cid7_xyz});
              
            }
          }
        }
        block.getCommunicator()._NeedMPIComm = true;
      } else if (nbaseblock.getLevel() > blocklevel + 1 ||
                 nbaseblock.getLevel() < blocklevel - 1) {
        std::cerr << "[BlockGeometry3D<T>::InitMPIIntpComm] Error: block level "
                     "difference larger than 1"
                  << std::endl;
      }
    }
  }
}

template <typename T>
void BlockGeometry3D<T>::InitAllMPIComm(BlockGeometryHelper3D<T> &GeoHelper) {
  InitMPIComm(GeoHelper);
  InitMPIAverComm(GeoHelper);
  InitMPIIntpComm(GeoHelper);
}

#endif

template <typename T>
void BlockGeometry3D<T>::AddtoSharedAverComm(std::vector<SharedComm> &Comms, 
const BasicBlock<T, 3>& baseblock_extx, Block3D<T>& block, Block3D<T> *nblock){
  // ------ add to Comms
  Comms.emplace_back(nblock->getBlockId());
  SharedComm &comm = Comms.back();
  // virtual coarse nblock
  BasicBlock<T, 3> nCblock = nblock->getCoasenedBlock();
  // get overlapped coarse cells
  std::vector<std::size_t> Recvs;
  std::vector<std::size_t> VSends;
  block.getCellIdx(baseblock_extx, nblock->getBaseBlock(), Recvs);
  nCblock.getCellIdx(nblock->getBaseBlock(), baseblock_extx, VSends);
  // get real send cell index
  std::vector<std::size_t> RSends;
  RSends.reserve(VSends.size()*8);
  for(std::size_t vid : VSends){
    std::vector<std::size_t> idxs;
    nCblock.getRefinedCellIdx(vid, idxs);
    RSends.insert(RSends.end(), idxs.begin(), idxs.end());
  }
  comm.setRecvSendIdx(Recvs, RSends);
}

template <typename T>
void BlockGeometry3D<T>::AddtoSharedIntpComm(std::vector<SharedComm> &Comms, Block3D<T>& block, Block3D<T> *nblock){
  // ------ add to Comms
  Comms.emplace_back(nblock->getBlockId());
  SharedComm &comm = Comms.back();
  // vox size
  const T Cvoxsize = nblock->getVoxelSize();
  const T Fvoxsize = block.getVoxelSize();
  // get intersection
  const AABB<T, 3> intsec = getIntersection(block, nblock->getBaseBlock());
  int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
  int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
  int CNz = static_cast<int>(std::round(intsec.getExtension()[2] / Cvoxsize));
  // get start index of intsec in nblock
  Vector<T, 3> startC = intsec.getMin() - nblock->getMin();
  // shift 1 voxel to left bottom for interpolation
  int startCx = static_cast<int>(std::round(startC[0] / Cvoxsize)) - 1;
  int startCy = static_cast<int>(std::round(startC[1] / Cvoxsize)) - 1;
  int startCz = static_cast<int>(std::round(startC[2] / Cvoxsize)) - 1;
  // start index of intsec in FBlock
  Vector<T, 3> startF = intsec.getMin() - block.getMin();
  int startFx = static_cast<int>(std::round(startF[0] / Fvoxsize));
  int startFy = static_cast<int>(std::round(startF[1] / Fvoxsize));
  int startFz = static_cast<int>(std::round(startF[2] / Fvoxsize));

  std::vector<std::size_t>& SendRecvs = comm.SendRecvCells;
  // (8 + 1) * 8
  SendRecvs.reserve(CNx*CNy*CNz*72);

  std::size_t XY = block.getNx() * block.getNy();
  std::size_t nXY = nblock->getNx() * nblock->getNy();
  for (int iz = 0; iz < CNz; ++iz) {
    for (int iy = 0; iy < CNy; ++iy) {
      for (int ix = 0; ix < CNx; ++ix) {
        // original
        std::size_t Cid0 = (iz + startCz) * nXY + (iy + startCy) * nblock->getNx() + ix + startCx;
        std::size_t Cid1 = Cid0 + 1;
        std::size_t Cid2 = Cid0 + nblock->getNx();
        std::size_t Cid3 = Cid2 + 1;
        std::size_t Cid4 = Cid0 + nXY;
        std::size_t Cid5 = Cid4 + 1;
        std::size_t Cid6 = Cid4 + nblock->getNx();
        std::size_t Cid7 = Cid6 + 1;
        std::size_t Fid = (iz * 2 + startFz) * XY + (iy * 2 + startFy) * block.getNx() + ix * 2 + startFx;

        // shift 1 voxel along +x direction
        std::size_t Cid0_x = Cid0 + 1;
        std::size_t Cid1_x = Cid1 + 1;
        std::size_t Cid2_x = Cid2 + 1;
        std::size_t Cid3_x = Cid3 + 1;
        std::size_t Cid4_x = Cid4 + 1;
        std::size_t Cid5_x = Cid5 + 1;
        std::size_t Cid6_x = Cid6 + 1;
        std::size_t Cid7_x = Cid7 + 1;
        std::size_t Fid_x = Fid + 1;

        // shift 1 voxel along +y direction
        std::size_t Cid0_y = Cid0 + nblock->getNx();
        std::size_t Cid1_y = Cid1 + nblock->getNx();
        std::size_t Cid2_y = Cid2 + nblock->getNx();
        std::size_t Cid3_y = Cid3 + nblock->getNx();
        std::size_t Cid4_y = Cid4 + nblock->getNx();
        std::size_t Cid5_y = Cid5 + nblock->getNx();
        std::size_t Cid6_y = Cid6 + nblock->getNx();
        std::size_t Cid7_y = Cid7 + nblock->getNx();
        std::size_t Fid_y = Fid + block.getNx();

        // shift 1 voxel along +z direction
        std::size_t Cid0_z = Cid0 + nXY;
        std::size_t Cid1_z = Cid1 + nXY;
        std::size_t Cid2_z = Cid2 + nXY;
        std::size_t Cid3_z = Cid3 + nXY;
        std::size_t Cid4_z = Cid4 + nXY;
        std::size_t Cid5_z = Cid5 + nXY;
        std::size_t Cid6_z = Cid6 + nXY;
        std::size_t Cid7_z = Cid7 + nXY;
        std::size_t Fid_z = Fid + XY;

        // 0
        SendRecvs.insert(SendRecvs.end(), {Cid0, Cid1, Cid2, Cid3, Cid4, Cid5, Cid6, Cid7});
        SendRecvs.push_back(Fid);

        // 1
        SendRecvs.insert(SendRecvs.end(), {Cid0_x, Cid1_x, Cid2_x, Cid3_x, Cid4_x, Cid5_x, Cid6_x, Cid7_x});
        SendRecvs.push_back(Fid_x);

        // 2
        SendRecvs.insert(SendRecvs.end(), {Cid0_y, Cid1_y, Cid2_y, Cid3_y, Cid4_y, Cid5_y, Cid6_y, Cid7_y});
        SendRecvs.push_back(Fid_y);

        // 3
        std::size_t Cid0_xy = Cid0_y + 1;
        std::size_t Cid1_xy = Cid1_y + 1;
        std::size_t Cid2_xy = Cid2_y + 1;
        std::size_t Cid3_xy = Cid3_y + 1;
        std::size_t Cid4_xy = Cid4_y + 1;
        std::size_t Cid5_xy = Cid5_y + 1;
        std::size_t Cid6_xy = Cid6_y + 1;
        std::size_t Cid7_xy = Cid7_y + 1;
        std::size_t Fid_xy = Fid_y + 1;

        SendRecvs.insert(SendRecvs.end(), {Cid0_xy, Cid1_xy, Cid2_xy, Cid3_xy, Cid4_xy, Cid5_xy, Cid6_xy, Cid7_xy});
        SendRecvs.push_back(Fid_xy);

        // 4
        SendRecvs.insert(SendRecvs.end(), {Cid0_z, Cid1_z, Cid2_z, Cid3_z, Cid4_z, Cid5_z, Cid6_z, Cid7_z});
        SendRecvs.push_back(Fid_z);

        // 5
        std::size_t Cid0_xz = Cid0_z + 1;
        std::size_t Cid1_xz = Cid1_z + 1;
        std::size_t Cid2_xz = Cid2_z + 1;
        std::size_t Cid3_xz = Cid3_z + 1;
        std::size_t Cid4_xz = Cid4_z + 1;
        std::size_t Cid5_xz = Cid5_z + 1;
        std::size_t Cid6_xz = Cid6_z + 1;
        std::size_t Cid7_xz = Cid7_z + 1;
        std::size_t Fid_xz = Fid_z + XY;

        SendRecvs.insert(SendRecvs.end(), {Cid0_xz, Cid1_xz, Cid2_xz, Cid3_xz, Cid4_xz, Cid5_xz, Cid6_xz, Cid7_xz});
        SendRecvs.push_back(Fid_xz);

        // 6
        std::size_t Cid0_yz = Cid0_z + nblock->getNx();
        std::size_t Cid1_yz = Cid1_z + nblock->getNx();
        std::size_t Cid2_yz = Cid2_z + nblock->getNx();
        std::size_t Cid3_yz = Cid3_z + nblock->getNx();
        std::size_t Cid4_yz = Cid4_z + nblock->getNx();
        std::size_t Cid5_yz = Cid5_z + nblock->getNx();
        std::size_t Cid6_yz = Cid6_z + nblock->getNx();
        std::size_t Cid7_yz = Cid7_z + nblock->getNx();
        std::size_t Fid_yz = Fid_z + XY;

        SendRecvs.insert(SendRecvs.end(), {Cid0_yz, Cid1_yz, Cid2_yz, Cid3_yz, Cid4_yz, Cid5_yz, Cid6_yz, Cid7_yz});
        SendRecvs.push_back(Fid_yz);

        // 7
        std::size_t Cid0_xyz = Cid0_yz + 1;
        std::size_t Cid1_xyz = Cid1_yz + 1;
        std::size_t Cid2_xyz = Cid2_yz + 1;
        std::size_t Cid3_xyz = Cid3_yz + 1;
        std::size_t Cid4_xyz = Cid4_yz + 1;
        std::size_t Cid5_xyz = Cid5_yz + 1;
        std::size_t Cid6_xyz = Cid6_yz + 1;
        std::size_t Cid7_xyz = Cid7_yz + 1;
        std::size_t Fid_xyz = Fid_yz + 1;

        SendRecvs.insert(SendRecvs.end(), {
        Cid0_xyz, Cid1_xyz, Cid2_xyz, Cid3_xyz, Cid4_xyz, Cid5_xyz, Cid6_xyz, Cid7_xyz});
        SendRecvs.push_back(Fid_xyz);
      }
    }
  }
}

template <typename T>
void BlockGeometry3D<T>::BuildBlockIndexMap() {
  _BlockIndexMap.clear();
  std::size_t count{};
  for (const Block3D<T> &block : _Blocks) {
    _BlockIndexMap[block.getBlockId()] = count;
    ++count;
  }
}


// --- helper functions ---
template <typename T>
bool IsInside(Octree<T>* tree, const BasicBlock<T, 3> &block) {
  for (int z = 0; z < block.getNz(); ++z) {
    for (int y = 0; y < block.getNy(); ++y) {
      for (int x = 0; x < block.getNx(); ++x) {
        const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{x, y, z});
        // get the node containing the voxel
        Octree<T>* node = tree->find(vox);
        if (node != nullptr) {
          // check if it is a [leaf] node and if it is [inside]
          if (node->isLeaf() && node->getInside()) {
            return true;
          }
        }
      }
    }
  }
  return false;
}
template <typename T>
bool hasOutSideCell(Octree<T>* tree, const BasicBlock<T, 3> &block) {
  for (int z = 0; z < block.getNz(); ++z) {
    for (int y = 0; y < block.getNy(); ++y) {
      for (int x = 0; x < block.getNx(); ++x) {
        const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{x, y, z});
        // get the node containing the voxel
        Octree<T>* node = tree->find(vox);
        if (node == nullptr) {
          return true;
        } else if (node->isLeaf() && !node->getInside()) {
          return true;
        }
      }
    }
  }
  return false;
}


// BlockGeometryHelper3D
#include "lbm/lattice_set.h"

template <typename T>
BlockGeometryHelper3D<T>::BlockGeometryHelper3D(int Nx, int Ny, int Nz,
                                                const AABB<T, 3> &AABBs, T voxelSize,
                                                int blockcelllen, int olap, int ext, bool useblockcell,
                                                std::uint8_t llimit)
  : BasicBlock<T, 3>(
      voxelSize, AABBs.getExtended(Vector<T, 3>{voxelSize * olap}),
      AABB<int, 3>(Vector<int, 3>{0}, Vector<int, 3>{Nx - 1 + 2*olap, Ny - 1 + 2*olap, Nz - 1 + 2*olap})),
    _BaseBlock(voxelSize, AABBs,
                AABB<int, 3>(Vector<int, 3>{olap}, Vector<int, 3>{Nx - 1 + olap, Ny - 1 + olap, Nz - 1 + olap})),
    BlockCellLen(blockcelllen), _Overlap(olap), _Ext(ext), _LevelLimit(llimit), _MaxLevel(std::uint8_t(0)), 
    _Exchanged(true), _IndexExchanged(true) {
  if (BlockCellLen < 4) {
    std::cerr << "BlockGeometryHelper3D<T>, BlockCellLen < 4" << std::endl;
  }

  if (useblockcell){
    CellsNx = std::ceil(T(_BaseBlock.getNx()) / T(BlockCellLen));
    CellsNy = std::ceil(T(_BaseBlock.getNy()) / T(BlockCellLen));
    CellsNz = std::ceil(T(_BaseBlock.getNz()) / T(BlockCellLen));
    CellsN = CellsNx * CellsNy * CellsNz;

    Vector<int, 3> Projection{1, CellsNx, CellsNx * CellsNy};

    Delta_Cellidx = make_Array<int, D3Q27<T>::q - 1>(
      [&](int i) { return D3Q27<T>::c[i + 1] * Projection; });

    CreateBlockCells();
  } else {
    std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllOldBasicBlocks();
    BasicBlocks.clear();
    _Exchanged = !_Exchanged;
    BasicBlocks.push_back(_BaseBlock);
  }
}

template <typename T>
BlockGeometryHelper3D<T>::BlockGeometryHelper3D(const StlReader<T>& reader, int blockcelllen, 
  int olap, int ext, bool useblockcell, std::uint8_t llimit)
  : BasicBlock<T, 3>(
      reader.getVoxelSize(), 
      reader.getMesh().getAABB().getExtended((olap+ext)*reader.getVoxelSize()),
      AABB<int, 3>(
        Vector<int, 3>{}, 
        Vector<int, 3>{
          int(std::ceil(reader.getMesh().getMax_Min()[0] / reader.getVoxelSize())) + 2*(olap+ext) - 1, 
          int(std::ceil(reader.getMesh().getMax_Min()[1] / reader.getVoxelSize())) + 2*(olap+ext) - 1,
          int(std::ceil(reader.getMesh().getMax_Min()[2] / reader.getVoxelSize())) + 2*(olap+ext) - 1})),
    _BaseBlock(
      reader.getVoxelSize(), 
      reader.getMesh().getAABB().getExtended(ext*reader.getVoxelSize()),
      AABB<int, 3>(
        Vector<int, 3>{olap}, 
        Vector<int, 3>{
          int(std::ceil(reader.getMesh().getMax_Min()[0] / reader.getVoxelSize())) + (olap+2*ext) - 1, 
          int(std::ceil(reader.getMesh().getMax_Min()[1] / reader.getVoxelSize())) + (olap+2*ext) - 1,
          int(std::ceil(reader.getMesh().getMax_Min()[2] / reader.getVoxelSize())) + (olap+2*ext) - 1})),
    BlockCellLen(blockcelllen), _Overlap(olap), _Ext(ext), _LevelLimit(llimit), _MaxLevel(std::uint8_t(0)), 
    _Exchanged(true), _IndexExchanged(true), _Reader(&reader) {
  if (BlockCellLen < 4) {
    std::cerr << "BlockGeometryHelper3D<T>, BlockCellLen < 4" << std::endl;
  }

  if (useblockcell){
    CellsNx = std::ceil(T(_BaseBlock.getNx()) / T(BlockCellLen));
    CellsNy = std::ceil(T(_BaseBlock.getNy()) / T(BlockCellLen));
    CellsNz = std::ceil(T(_BaseBlock.getNz()) / T(BlockCellLen));
    CellsN = CellsNx * CellsNy * CellsNz;

    // correct the mesh size by CellsNx, CellsNy, CellsNz
    int NewNx = CellsNx * BlockCellLen;
    int NewNy = CellsNy * BlockCellLen;
    int NewNz = CellsNz * BlockCellLen;
    T MeshSizeX = NewNx * _BaseBlock.getVoxelSize();
    T MeshSizeY = NewNy * _BaseBlock.getVoxelSize();
    T MeshSizeZ = NewNz * _BaseBlock.getVoxelSize();

    // ext is already included in the _BaseBlock, CellsNx, CellsNy, CellsNz
    static_cast<BasicBlock<T, 3>&>(*this) = BasicBlock<T, 3>(
      reader.getVoxelSize(), 
      AABB<T, 3>(reader.getMesh().getMin() - (olap+ext)*reader.getVoxelSize(), 
                reader.getMesh().getMin() - (olap+ext)*reader.getVoxelSize() + Vector<T, 3>{MeshSizeX, MeshSizeY, MeshSizeZ} + 2*olap*reader.getVoxelSize()),
      AABB<int, 3>(Vector<int, 3>{}, 
                  Vector<int, 3>{NewNx + 2*olap-1, NewNy + 2*olap-1, NewNz + 2*olap-1}));
    _BaseBlock = BasicBlock<T, 3>(
      reader.getVoxelSize(), 
      AABB<T, 3>(reader.getMesh().getMin() - ext*reader.getVoxelSize(), 
                reader.getMesh().getMin() - ext*reader.getVoxelSize() + Vector<T, 3>{MeshSizeX, MeshSizeY, MeshSizeZ}),
      AABB<int, 3>(Vector<int, 3>{olap}, 
                  Vector<int, 3>{NewNx + olap-1, NewNy + olap-1, NewNz + olap-1}));

    // end correct

    Vector<int, 3> Projection{1, CellsNx, CellsNx * CellsNy};

    Delta_Cellidx = make_Array<int, D3Q27<T>::q - 1>(
      [&](int i) { return D3Q27<T>::c[i + 1] * Projection; });

    CreateBlockCells();
    TagBlockCells(reader);
  } else {
    // even if we do not use block cell, BasicBlock<T, 3> and _BaseBlock should be corrected
    // to make mesh size = integer multiple of voxel size
    const int Nx = _BaseBlock.getNx();
    const int Ny = _BaseBlock.getNy();
    const int Nz = _BaseBlock.getNz();
    const T voxsize = _BaseBlock.getVoxelSize();
    const T meshNx = Nx * voxsize;
    const T meshNy = Ny * voxsize;
    const T meshNz = Nz * voxsize;

    static_cast<BasicBlock<T, 3>&>(*this) = BasicBlock<T, 3>(
      voxsize, 
      AABB<T, 3>(reader.getMesh().getMin() - (olap+ext)*voxsize, 
                reader.getMesh().getMin() - (olap+ext)*voxsize + Vector<T, 3>{meshNx, meshNy, meshNz} + 2*olap*voxsize),
      AABB<int, 3>(Vector<int, 3>{}, 
                  Vector<int, 3>{Nx + 2*olap-1, Ny + 2*olap-1, Nz + 2*olap-1}));
    
    _BaseBlock = BasicBlock<T, 3>(
      voxsize, 
      AABB<T, 3>(reader.getMesh().getMin() - ext*voxsize, 
                reader.getMesh().getMin() - ext*voxsize + Vector<T, 3>{meshNx, meshNy, meshNz}),
      AABB<int, 3>(Vector<int, 3>{olap}, 
                  Vector<int, 3>{Nx + olap-1, Ny + olap-1, Nz + olap-1}));

    std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllOldBasicBlocks();
    BasicBlocks.clear();
    _Exchanged = !_Exchanged;
    BasicBlocks.push_back(_BaseBlock);
  }
}

template <typename T>
BlockGeometryHelper3D<T>::BlockGeometryHelper3D(const BlockReader<T,3>& blockreader, bool useReaderOlap) 
    : BasicBlock<T, 3>(blockreader.getBasicBlock()), _BaseBlock(blockreader.getBaseBlock()), 
      BlockCellLen(0), _Overlap(1), _Ext(0), _LevelLimit(std::uint8_t{}), _MaxLevel(blockreader.getMaxLevel()),
      _Exchanged(true), _IndexExchanged(true) {
  // create new blocks on relatively older blocks
  std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllOldBasicBlocks();
  BasicBlocks.clear();
  // now old blocks become new blocks
  _Exchanged = !_Exchanged;
  // create blocks from Block Reader
  int iblock{};
  for (const BasicBlock<T, 3> &baseblock : blockreader.getBlocks()) {
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
void BlockGeometryHelper3D<T>::CreateBlockCells() {
  _BlockCells.clear();
  _BlockCellTags.clear();
  // buffer vector to store cell aabbs
  std::vector<AABB<int, 3>> AABBCells;
  AABBCells.reserve(CellsN);
  // divide base block into cell aabbs
  _BaseBlock.getIdxBlock().divide(CellsNx, CellsNy, CellsNz, AABBCells);
  // create cell blocks from cell aabbs
  _BlockCells.reserve(CellsN);
  int blockid = 0;
  const T voxsize = BasicBlock<T, 3>::getVoxelSize();
  Vector<T, 3> _Min = BasicBlock<T, 3>::_min;
  for (const AABB<int, 3> &aabbcell : AABBCells) {
    Vector<T, 3> MIN = aabbcell.getMin() * voxsize + _Min;
    Vector<T, 3> MAX = (aabbcell.getMax() + Vector<T, 3>{T(1)}) * voxsize + _Min;
    AABB<T, 3> aabb(MIN, MAX);
    _BlockCells.emplace_back(voxsize, aabb, aabbcell, blockid);
    _BlockCellTags.emplace_back(BlockCellTag::none);
    ++blockid;
  }
}

template <typename T>
void BlockGeometryHelper3D<T>::TagBlockCells(const StlReader<T>& reader) {
  Octree<T>* tree = reader.getTree();
  std::size_t i{};
  for (const BasicBlock<T, 3> &blockcell : _BlockCells) {
    if (IsInside(tree, blockcell)) _BlockCellTags[i] = BlockCellTag::inside;
    ++i;
  }
}
template <typename T>
void BlockGeometryHelper3D<T>::TagBlockCells(std::uint8_t voidflag) {
  std::size_t i{};
  for (const BasicBlock<T, 3> &blockcell : _BlockCells) {
    
    bool has_stl = false;
    for (int z = 0; z < blockcell.getNz(); ++z) {
      for (int y = 0; y < blockcell.getNy(); ++y) {
        for (int x = 0; x < blockcell.getNx(); ++x) {
          const Vector<T, 3> vox = blockcell.getVoxel(Vector<int, 3>{x, y, z});
          if (_FlagField[vox] != voidflag) {
            has_stl = true;
            goto check_has_stl;
          }
        }
      }
    }
    check_has_stl:
    if (has_stl) {
      _BlockCellTags[i] = BlockCellTag::inside;
    }
    ++i;
  }
}


template <typename T>
void BlockGeometryHelper3D<T>::UpdateMaxLevel() {
  _MaxLevel = std::uint8_t(0);
  for (const BasicBlock<T, 3> &block : _BlockCells) {
    if (block.getLevel() > _MaxLevel) _MaxLevel = block.getLevel();
  }
}

template <typename T>
void BlockGeometryHelper3D<T>::CreateBlocks(bool CreateFromInsideTag, bool outputinfo) {
  // collect info
  std::vector<std::size_t> BlockCellNumVec;

  // create new blocks on relatively older blocks
  std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllOldBasicBlocks();
  BasicBlocks.clear();
  // now old blocks become new blocks
  _Exchanged = !_Exchanged;
  // get max refine level
  // int maxlevel = 0;
  // for (const BasicBlock<T, 3>& block : _BlockCells) {
  //   maxlevel = std::max(maxlevel, static_cast<int>(block.getLevel()));
  // }
  // buffer array visited
  bool *visited = new bool[CellsN];
  std::fill(visited, visited + CellsN, false);
  // create blocks
  int blockid = 0;
  int XY = CellsNx * CellsNy;
  for (int k = 0; k < CellsNz; ++k) {
    for (int j = 0; j < CellsNy; ++j) {
      for (int i = 0; i < CellsNx; ++i) {
        std::size_t id = i + j * CellsNx + k * XY;
        if (CreateFromInsideTag && !util::isFlag(_BlockCellTags[id], BlockCellTag::inside)) {
          continue;
        }
        if (visited[id]) {
          continue;
        }
        // get level
        std::uint8_t level = _BlockCells[id].getLevel();
        T voxsize = _BlockCells[id].getVoxelSize();
        Vector<int, 3> NewMesh = _BlockCells[id].getMesh();

        // expand block along x
        int Nx = 1;
        while (i + Nx < CellsNx) {
          std::size_t tempid = id + Nx;
          if (_BlockCells[tempid].getLevel() != level || visited[tempid]) {
            break;
          }
          if (CreateFromInsideTag && !util::isFlag(_BlockCellTags[tempid], BlockCellTag::inside)) {
            break;
          }
          NewMesh[0] += _BlockCells[tempid].getNx();
          ++Nx;
        }

        // expand block along y
        int Ny = 1;
        while (j + Ny < CellsNy) {
          std::size_t startid = id + Ny * CellsNx;
          for (int iNx = 0; iNx < Nx; ++iNx) {
            std::size_t tempid = startid + iNx;
            if (_BlockCells[tempid].getLevel() != level || visited[tempid]) {
              goto end_y_expansion;
            }
            if (CreateFromInsideTag && !util::isFlag(_BlockCellTags[tempid], BlockCellTag::inside)) {
              goto end_y_expansion;
            }
          }
          NewMesh[1] += _BlockCells[startid].getNy();
          ++Ny;
        }
      end_y_expansion:
        // expand block along z
        int Nz = 1;
        while (k + Nz < CellsNz) {
          for (int jNy = 0; jNy < Ny; ++jNy) {
            std::size_t startid = id + jNy * CellsNx + Nz * XY;
            for (int iNx = 0; iNx < Nx; ++iNx) {
              std::size_t tempid = startid + iNx;
              if (_BlockCells[tempid].getLevel() != level || visited[tempid]) {
                goto end_z_expansion;
              }
              if (CreateFromInsideTag && !util::isFlag(_BlockCellTags[tempid], BlockCellTag::inside)) {
                goto end_z_expansion;
              }
            }
          }
          NewMesh[2] += _BlockCells[id + Nz * XY].getNz();
          ++Nz;
        }
      end_z_expansion:

        BlockCellNumVec.push_back(Nx * Ny * Nz);
        // create block
        Vector<int, 3> Ext = BlockCellLen * Vector<int, 3>{Nx, Ny, Nz};
        Vector<int, 3> min = _BlockCells[id].getIdxBlock().getMin();
        Vector<int, 3> max = min + Ext - Vector<int, 3>{1};
        AABB<int, 3> idxblock(min, max);
        Vector<T, 3> MIN = _BlockCells[id].getMin();
        int ratio = static_cast<int>(std::pow(2, static_cast<int>(level)));
        Vector<T, 3> MAX = Ext * voxsize * ratio + MIN;
        AABB<T, 3> block(MIN, MAX);
        BasicBlocks.emplace_back(level, voxsize, blockid, block, idxblock, NewMesh);
        blockid++;
        // set visited
        for (int kk = 0; kk < Nz; ++kk) {
          for (int jj = 0; jj < Ny; ++jj) {
            std::size_t tempid = id + jj * CellsNx + kk * XY;
            for (int ii = 0; ii < Nx; ++ii) {
              visited[tempid + ii] = true;
            }
          }
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

  if (!outputinfo) return;
  MPI_RANK(0)
  std::cout << "[BlockGeometryHelper3D<T>::CreateBlocks]:\n" 
            << "  Created " << BasicBlocks.size() << " Blocks\n"
            << "  [BlockCell Num]: \n" 
            << "    min: " << min << ", max: " << max << ", average: " << aver << std::endl;
}

template <typename T>
void BlockGeometryHelper3D<T>::CreateBlocks(int blockXNum, int blockYNum, int blockZNum) {
  // create new blocks on relatively older blocks
  std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllOldBasicBlocks();
  BasicBlocks.clear();
  // now old blocks become new blocks
  _Exchanged = !_Exchanged;
  // divide blocks manually
  DivideBlock3D(_BaseBlock, blockXNum, blockYNum, blockZNum, BasicBlocks);
}

template <typename T>
template <typename Func>
void BlockGeometryHelper3D<T>::forEachBlockCell(const Func& func) {
  for (BasicBlock<T, 3> &block : _BlockCells) {
    func(block);
  }
}

template <typename T>
void BlockGeometryHelper3D<T>::AdaptiveOptimization(int OptProcNum, int MaxProcNum,
                                                    bool enforce) {
  // check
  if (OptProcNum < 1) {
    std::cerr << "[BlockGeometryHelper3D<T>::AdaptiveOptimization]: OptProcNum < 1"
              << std::endl;
  }
  if (MaxProcNum == -1) {
    MaxProcNum = 2 * OptProcNum;
  }
  if (MaxProcNum < 1) {
    std::cerr << "[BlockGeometryHelper3D<T>::AdaptiveOptimization]: MaxProcNum < 1"
              << std::endl;
  }
  // get new basic blocks
  std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllBasicBlocks();
  // vector store stddev of each optimization scheme
  std::vector<T> StdDevs;
  int minProcNum = std::min(OptProcNum, static_cast<int>(BasicBlocks.size()));
  int maxProcNum = std::max(MaxProcNum, static_cast<int>(BasicBlocks.size()));
  for (int i = minProcNum; i <= maxProcNum; ++i) {
    // buffer vector copied from BasicBlocks
    std::vector<BasicBlock<T, 3>> Blocks = BasicBlocks;
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
  std::cout << "[BlockGeometryHelper3D<T>::AdaptiveOptimization]:\n"
            << "  BlockNum: " << BasicBlocks.size() << " with stdDev: " << minStdDev << std::endl;
}

template <typename T>
void BlockGeometryHelper3D<T>::Optimize(int ProcessNum, bool enforce, bool info) {
  Optimize(getAllBasicBlocks(), ProcessNum, enforce);
  if (info) {
    MPI_RANK(0)
    std::cout << "[BlockGeometryHelper3D<T>::Optimize]:\n"
              << "  BlockNum: " << getAllBasicBlocks().size() 
              << " with stdDev: " << ComputeBlockNStdDev(getAllBasicBlocks()) << std::endl;
  }
}

template <typename T>
void BlockGeometryHelper3D<T>::Optimize(std::vector<BasicBlock<T, 3>> &Blocks,
                                        int ProcessNum, bool enforce) {
  // get total number of points
  std::size_t Total = 0;
  for (const BasicBlock<T, 3> &block : Blocks) {
    Total += block.getN();
  }
  // get number of points per process
  std::size_t NumPerProcess = Total / ProcessNum;

  // divide large blocks
  // T threshold = static_cast<T>(ProcessNum) / size;
  // buffer vector to store new blocks
  std::vector<BasicBlock<T, 3>> NewBasicBlocks;
  // iterate through all blocks
  typename std::vector<BasicBlock<T, 3>>::iterator it = Blocks.begin();

  // [enforce] will try to divide blocks into ProcessNum parts
  // but will NOT divide blocks with number of cells less than NumPerProcess
  if (enforce && Blocks.size() < static_cast<std::size_t>(ProcessNum)) {
    // sort blocks in descending order based on the number of cells
    std::sort(Blocks.begin(), Blocks.end(),
              [](const BasicBlock<T, 3> &a, const BasicBlock<T, 3> &b) {
                return a.getN() > b.getN();
              });
    // make size = ProcessNum
    std::size_t count = Blocks.size();
    while (count < static_cast<std::size_t>(ProcessNum) && it != Blocks.end()) {
      T ratio = static_cast<T>(it->getN()) / NumPerProcess;
      if (ratio <= T{1}) { break; }
      int part = static_cast<int>(ratio);
      if (ratio <= T{2}) { part = 2; }
      if ((count + part - 1) > static_cast<std::size_t>(ProcessNum)) {
        part = ProcessNum - count + 1;
      }
      if (part == 1) { break; }
      // divide block
      DivideBlock3D(*it, part, NewBasicBlocks);
      // remove block
      it = Blocks.erase(it);
      count += part - 1;
    }
  } else {
    while (it != Blocks.end()) {
      T ratio = std::round(static_cast<T>(it->getN()) / NumPerProcess);
      if (ratio >= T{2}) {
        // divide block
        DivideBlock3D(*it, static_cast<int>(ratio), NewBasicBlocks);
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

// only for movable neighbor blocks with face overlap
struct FaceNbrInfo {
  // neighbor block id
  std::size_t NbrBlockId;
  // neighbor direction
  NbrDirection Direction;
  // is face movable
  bool Movable;

  FaceNbrInfo(std::size_t id, NbrDirection dir) : NbrBlockId(id), Direction(dir), Movable(false) {}
  FaceNbrInfo(std::size_t id, NbrDirection dir, bool movable) : NbrBlockId(id), Direction(dir), Movable(movable) {}
};

template <typename T>
void buildFaceNbrInfosImpl(const std::vector<BasicBlock<T, 3>> &Blocks,
                       std::vector<std::vector<FaceNbrInfo>> &NbrInfos) {
  if (Blocks.size() == 1) return;
  NbrInfos.clear();
  NbrInfos.resize(Blocks.size(), std::vector<FaceNbrInfo>{});

  // build NbrInfos
  for (std::size_t iblock = 0; iblock < Blocks.size(); ++iblock) {
    std::vector<FaceNbrInfo> &NbrInfo = NbrInfos[iblock];
    NbrInfo.clear();

    const BasicBlock<T, 3> &block = Blocks[iblock];
    for (std::size_t jblock = 0; jblock < Blocks.size(); ++jblock) {
      if (iblock == jblock) continue;

      const BasicBlock<T, 3> &nblock = Blocks[jblock];
      if (block.whichNbrType(nblock) == 3) {
        // has face neighbor
        NbrDirection dir = getFaceNbrDirection(block.whichFace(nblock));
        // check if the face is movable
        AABB<T, 3> face = block.getFace(dir);
        AABB<T, 3> nface = nblock.getFace(getOpposite(dir));
        if (util::Vectorfpequal(face.getMin(), nface.getMin()) &&
            util::Vectorfpequal(face.getMax(), nface.getMax())) {
          // if 2 blocks sharing the same face
          NbrInfo.emplace_back(jblock, dir, true);
        } else {
          NbrInfo.emplace_back(jblock, dir, false);
        }
      }
    }
    // check NbrInfo
    // 2025/1/16: in some cases, NbrInfo of one block may be empty
    // we do not know how to handle this
    // for now each time we use NbrInfo, we check if it is empty
    // if (NbrInfo.empty()) {
    //   IF_MPI_RANK(0){
    //     std::cerr << "[buildFaceNbrInfos<T>]: NbrInfo is empty" << std::endl;
    //   }
    // }
  }
}
template <typename T>
void buildFaceNbrInfos(std::vector<BasicBlock<T, 3>> &Blocks,
                       std::vector<std::vector<FaceNbrInfo>> &NbrInfos) {
  if (Blocks.size() == 1) return;
  // build NbrInfos
  buildFaceNbrInfosImpl(Blocks, NbrInfos);
  // 2025/1/17: empty NbrInfo may be caused by RemoveUnusedCells()
  // block without face neighbor is bad block
  // we could remove the block if NbrInfo is empty
  // check NbrInfos of each block and remove the block if NbrInfo is empty
  bool removed = true;
  while (removed) {
    removed = false;
    typename std::vector<BasicBlock<T, 3>>::iterator it = Blocks.begin();
    while (it != Blocks.end()) {
      if (NbrInfos[it->getBlockId()].empty()) {
        it = Blocks.erase(it);
        removed = true;
      } else {
        ++it;
      }
    }
    // if removed, rebuild block index and NbrInfos
    if (removed) {
      for (std::size_t i = 0; i < Blocks.size(); ++i) Blocks[i].setBlockId(i);
      buildFaceNbrInfosImpl(Blocks, NbrInfos);
    }
  }
}

// this will try to make the number of cells in the base block as close to the bestN as possible
template <typename T>
void AdjustAdjacentBlockSize(BasicBlock<T, 3>& base, BasicBlock<T, 3>& nbr, std::size_t bestN, NbrDirection direction) {
  const int baseNx = base.getNx();
  const int baseNy = base.getNy();
  const int baseNz = base.getNz();
  // delta > 0, base: shrink, nbr: expand
  // delta < 0, base: expand, nbr: shrink
  const int delta = base.getN() - bestN;
  if (delta == 0) return;

  const int nbrNx = nbr.getNx();
  const int nbrNy = nbr.getNy();
  const int nbrNz = nbr.getNz();

  int shift{};
  if (util::isFlag(direction, NbrDirection::XN | NbrDirection::XP)) {
    shift = static_cast<int>(std::round(T(delta) / (baseNy * baseNz)));
    if (shift == 0) return;
    if (shift < 0 && nbrNx + shift < 2) return;
  } else if (util::isFlag(direction, NbrDirection::YN | NbrDirection::YP)) {
    shift = static_cast<int>(std::round(T(delta) / (baseNx * baseNz)));
    if (shift == 0) return;
    if (shift < 0 && nbrNy + shift < 2) return;
  } else if (util::isFlag(direction, NbrDirection::ZN | NbrDirection::ZP)) {
    shift = static_cast<int>(std::round(T(delta) / (baseNx * baseNy)));
    if (shift == 0) return;
    if (shift < 0 && nbrNz + shift < 2) return;
  }

  base.resize(-shift, direction);
  nbr.resize(shift, getOpposite(direction));
}

// this will try to make the number of cells in the input 2 blocks as close as possible
template <typename T>
void AdjustAdjacentBlockSize(BasicBlock<T, 3>& base, BasicBlock<T, 3>& nbr, NbrDirection direction) {
  const int baseNx = base.getNx();
  const int baseNy = base.getNy();
  const int baseNz = base.getNz();
  // delta > 0, base: shrink, nbr: expand
  // delta < 0, base: expand, nbr: shrink
  const int delta = base.getN() + nbr.getN() - 2 * base.getN();
  if (delta == 0) return;

  const int nbrNx = nbr.getNx();
  const int nbrNy = nbr.getNy();
  const int nbrNz = nbr.getNz();

  int shift{};
  if (util::isFlag(direction, NbrDirection::XN | NbrDirection::XP)) {
    shift = static_cast<int>(std::round(T(delta) / (baseNy * baseNz)));
    if (shift == 0) return;
    if (shift > 0 && baseNx - shift < 2) return;
    if (shift < 0 && nbrNx + shift < 2) return;
  } else if (util::isFlag(direction, NbrDirection::YN | NbrDirection::YP)) {
    shift = static_cast<int>(std::round(T(delta) / (baseNx * baseNz)));
    if (shift == 0) return;
    if (shift > 0 && baseNy - shift < 2) return;
    if (shift < 0 && nbrNy + shift < 2) return;
  } else if (util::isFlag(direction, NbrDirection::ZN | NbrDirection::ZP)) {
    shift = static_cast<int>(std::round(T(delta) / (baseNx * baseNy)));
    if (shift == 0) return;
    if (shift > 0 && baseNz - shift < 2) return;
    if (shift < 0 && nbrNz + shift < 2) return;
  }

  base.resize(-shift, direction);
  nbr.resize(shift, getOpposite(direction));
}

// this will adjust the block size of the input 2 blocks in [dir] by [step](usually small integer)
template <typename T>
void AdjustAdjacentBlockSize_byStep(BasicBlock<T, 3>& base, BasicBlock<T, 3>& nbr, std::size_t bestN,
  NbrDirection direction, int step = 1) {

  const int baseNx = base.getNx();
  const int baseNy = base.getNy();
  const int baseNz = base.getNz();
  // delta > 0, base: shrink, nbr: expand
  // delta < 0, base: expand, nbr: shrink
  const int delta = base.getN() - bestN;
  if (delta == 0) return;
  if (delta < 0) step = -step;

  const int nbrNx = nbr.getNx();
  const int nbrNy = nbr.getNy();
  const int nbrNz = nbr.getNz();

  if (util::isFlag(direction, NbrDirection::XN | NbrDirection::XP)) {
    if (step > 0 && baseNx - step < 2) return;
    if (step < 0 && nbrNx + step < 2) return;
  } else if (util::isFlag(direction, NbrDirection::YN | NbrDirection::YP)) {
    if (step > 0 && baseNy - step < 2) return;
    if (step < 0 && nbrNy + step < 2) return;
  } else if (util::isFlag(direction, NbrDirection::ZN | NbrDirection::ZP)) {
    if (step > 0 && baseNz - step < 2) return;
    if (step < 0 && nbrNz + step < 2) return;
  }

  base.resize(-step, direction);
  nbr.resize(step, getOpposite(direction));
}

// Blocks[id] will be merged into Blocks[NbrId]
// note that MergeBlock() will NOT update NbrInfos
template <typename T>
bool MergeBlock(std::vector<BasicBlock<T, 3>> &Blocks, 
  std::vector<std::vector<FaceNbrInfo>>& NbrInfos, std::vector<FaceNbrInfo>& IdNbrInfo, std::size_t id, std::size_t NbrId) {
  // find which face is shared by Blocks[id] and Blocks[NbrId]
  NbrDirection dir{};
  for (const FaceNbrInfo &info : NbrInfos[NbrId]) {
    if (info.NbrBlockId == id) {
      dir = info.Direction;
      break;
    }
  }
  // check if the merge operation could be performed
  BasicBlock<T, 3> &smallblock = Blocks[id];
  BasicBlock<T, 3> &nblock = Blocks[NbrId];
  AABB<T, 3> smallface = smallblock.getFace(getOpposite(dir));
  AABB<T, 3> nface = nblock.getFace(dir);
  if (isInside(smallface, nface)) {
    // this is a test block buffer to check if the resize is valid
    BasicBlock<T, 3> nblock_buffer = nblock;
    // adjust block size along dir
    if (util::isFlag(dir, NbrDirection::XN | NbrDirection::XP)) {
      nblock_buffer.resize(Blocks[id].getNx(), dir);
    } else if (util::isFlag(dir, NbrDirection::YN | NbrDirection::YP)) {
      nblock_buffer.resize(Blocks[id].getNy(), dir);
    } else if (util::isFlag(dir, NbrDirection::ZN | NbrDirection::ZP)) {
      nblock_buffer.resize(Blocks[id].getNz(), dir);
    }
    // check if the resize is valid: use overlapped criteria
    bool valid = true;
    for (const BasicBlock<T, 3> &block : Blocks) {
      if (block.getBlockId() == int(id) || block.getBlockId() == int(NbrId)) continue;
      if (isOverlapped(nblock_buffer, block)) valid = false;
    }
    if (valid) {
      // perform resize
      nblock = nblock_buffer;
      return true;
    }
  }
  // if not inside or resize is invalid:
  // erase Blocks[NbrId] from nbr list of Blocks[id]
  typename std::vector<FaceNbrInfo>::iterator it = IdNbrInfo.begin();
  while (it != IdNbrInfo.end()) {
    if (it->NbrBlockId == NbrId) it = IdNbrInfo.erase(it);
    else ++it;
  }
  if (IdNbrInfo.empty()) return false;
  // find another neighbor block
  std::size_t NbrIdx = IdNbrInfo[0].NbrBlockId;
  std::size_t minN = Blocks[NbrIdx].getN();
  for (const FaceNbrInfo &info : IdNbrInfo) {
    std::size_t nid = info.NbrBlockId;
    if (Blocks[nid].getN() < minN) {
      NbrIdx = nid;
      minN = Blocks[nid].getN();
    }
  }
  // recursively merge Blocks[id] into Blocks[NbrIdx]
  return MergeBlock(Blocks, NbrInfos, IdNbrInfo, id, NbrIdx);
}

template <typename T>
void MergeSmallBlocks(std::vector<BasicBlock<T, 3>> &Blocks, 
  std::vector<std::vector<FaceNbrInfo>>& NbrInfos, T threshold = T{0.25}) {
  // get total number of cells
  std::size_t Total{};
  for (const BasicBlock<T, 3> &block : Blocks) {
    Total += block.getN();
  }
  // get number of cells per process
  std::size_t AverageN = Total / Blocks.size();
  std::size_t thN = static_cast<std::size_t>(AverageN * threshold);
  // iterate through all blocks and find small blocks, N < AverageN * threshold
  // smallblocks will be erased and merged into neighbor blocks with the least cells
  typename std::vector<BasicBlock<T, 3>>::iterator it = Blocks.begin();
  while (it != Blocks.end()) {
    if (it->getN() <= thN) {
      // find the neighbor block with the least cells
      std::size_t id = it->getBlockId();
      // check empty
      // if (NbrInfos[id].empty()) {
      //   // we can't perform merge operation
      //   ++it;
      //   continue;
      // }
      std::size_t NbrId = NbrInfos[id][0].NbrBlockId;
      std::size_t minN = Blocks[NbrId].getN();
      for (const FaceNbrInfo &info : NbrInfos[id]) {
        std::size_t nid = info.NbrBlockId;
        if (Blocks[nid].getN() < minN) {
          NbrId = nid;
          minN = Blocks[nid].getN();
        }
      }
      // merge small block into the neighbor block
      std::vector<FaceNbrInfo> IdNbrInfo = NbrInfos[id];
      bool merged = MergeBlock(Blocks, NbrInfos, IdNbrInfo, id, NbrId);
      if (merged) {
        // remove small block
        it = Blocks.erase(it);
        // update block id
        for (std::size_t i = 0; i < Blocks.size(); ++i) {
          Blocks[i].setBlockId(i);
        }
        // update NbrInfos
        buildFaceNbrInfos(Blocks, NbrInfos);
      } else {
        ++it;
      }
    } else {
      ++it;
    }
  }
}

template <typename T>
void ForcedMergeSmallBlocks(std::vector<BasicBlock<T, 3>> &Blocks, 
  std::vector<std::vector<FaceNbrInfo>>& NbrInfos, int TargetBlockNum) {
  if (Blocks.size() <= static_cast<std::size_t>(TargetBlockNum)) return;
  // sort blocks in ascending order based on the number of cells
  std::sort(Blocks.begin(), Blocks.end(),
            [](const BasicBlock<T, 3> &a, const BasicBlock<T, 3> &b) {
              return a.getN() < b.getN();
            });
  // DO NOT Forget to update block id and NbrInfos
  for (std::size_t i = 0; i < Blocks.size(); ++i) {
    Blocks[i].setBlockId(i);
  }
  buildFaceNbrInfos(Blocks, NbrInfos);
  // merge small blocks
  typename std::vector<BasicBlock<T, 3>>::iterator it = Blocks.begin();
  while (Blocks.size() > static_cast<std::size_t>(TargetBlockNum) && it != Blocks.end()) {
    std::size_t id = it->getBlockId();
    // check empty
    // if (NbrInfos[id].empty()) {
    //   // we can't perform merge operation
    //   ++it;
    //   continue;
    // }
    std::size_t NbrId = NbrInfos[id][0].NbrBlockId;
    std::size_t minN = Blocks[NbrId].getN();
    for (const FaceNbrInfo &info : NbrInfos[id]) {
      std::size_t nid = info.NbrBlockId;
      if (Blocks[nid].getN() < minN) {
        NbrId = nid;
        minN = Blocks[nid].getN();
      }
    }
    // merge small block into the neighbor block
    std::vector<FaceNbrInfo> IdNbrInfo = NbrInfos[id];
    bool merged = MergeBlock(Blocks, NbrInfos, IdNbrInfo, id, NbrId);
    if (merged) {
      // remove small block
      it = Blocks.erase(it);
      // update block id
      for (std::size_t i = 0; i < Blocks.size(); ++i) {
        Blocks[i].setBlockId(i);
      }
      // update NbrInfos
      buildFaceNbrInfos(Blocks, NbrInfos);
    } else {
      ++it;
    }
  }
}


template <typename T>
T BlockGeometryHelper3D<T>::IterateAndOptimizeImp(int ProcNum, int BlockCellLen) {
  // ------ create blocks and remove unused cells ------
  // 1. initialize BlockGeometryHelper3D with BlockCellLen
  //    blockcells will be tagged with _FlagField instead of _Reader
  //    make sure that _FlagField is initialized before calling this function
  Init(BlockCellLen);
  // 2. create blocks
  CreateBlocks(true, false);
  // if (_Ext != 0) AddVoidCellLayer(*_Reader, false);
  // 3. remove unused cells using _FlagField
  RemoveUnusedCells(std::uint8_t{1}, false);
  //    build face neighbor infos for current BlockCellLen
  std::vector<std::vector<FaceNbrInfo>> NbrInfos;
  std::vector<BasicBlock<T, 3>> &Blocks = getAllBasicBlocks();
  
  // ------ divide-merge loop ------
  T maxRatio = T{2};
  T minRatio = T{0.1};
  int iter{};
  while (maxRatio > T{1.2} && minRatio < T{0.5} && iter < 100) {
    ++iter;
  // 4. divide blocks into ProcNum parts, but will NOT enforce exactly ProcNum parts
  //    even if the Blocks.size() > ProcNum, when one block has number of cells more than twice of the average
  //    it will be divided into x parts, x = std::round(block.getN() / average)
    Optimize(ProcNum, false, false);
  //    update block id before building NbrInfos
    for (std::size_t i = 0; i < Blocks.size(); ++i) Blocks[i].setBlockId(i);
    buildFaceNbrInfos(Blocks, NbrInfos);
  // 5. merge small blocks, use default threshold = 0.125
  //    this step mainly aims to eleminate small blocks created in step 2: CreateBlocks()
  //    if Blocks.size() > ProcNum, ForcedMergeSmallBlocks() will enforce exactly ProcNum parts
  //    otherwise, MergeSmallBlocks() will not enforce exactly ProcNum parts, instead, this will be done in step 6
    if (Blocks.size() > static_cast<std::size_t>(ProcNum)) 
      ForcedMergeSmallBlocks(Blocks, NbrInfos, ProcNum);
    else
      MergeSmallBlocks(Blocks, NbrInfos);
  // 6. divide blocks again, this time enforce exactly ProcNum parts
    Optimize(ProcNum, true, false);
  //    and remove unused cells using _FlagField
    RemoveUnusedCells(std::uint8_t{1}, false);
  //    get max and min ratio
    if (Blocks.size() != static_cast<std::size_t>(ProcNum)) continue;
    std::size_t Total{};
    for (const BasicBlock<T, 3> &block : Blocks) Total += block.getN();
    std::size_t AverageN = Total / ProcNum;
    for (const BasicBlock<T, 3> &block : Blocks) {
      T ratio = static_cast<T>(block.getN()) / AverageN;
      maxRatio = std::max(maxRatio, ratio);
      minRatio = std::min(minRatio, ratio);
    }
  }

  // ------ load optimization ------
  // 7. do load optimization, this will go beyond the limit of blockcells 
  //    and try to make the number of cells in each block as close as possible
  LoadOptimization(500, 0.01, false);
  // 8. get stdDev
  return ComputeBlockNStdDev(Blocks);
}

template <typename T>
void BlockGeometryHelper3D<T>::IterateAndOptimize(int ProcNum, int MinBlockCellLen, int MaxBlockCellLen, bool stepinfo) {
  T stdDev = std::numeric_limits<T>::max();
  int bestBlockCellLen = MinBlockCellLen;
  IF_MPI_RANK(0) {
    std::cout << "[BlockGeometryHelper3D<T>::IterateAndOptimize]:\n"
              << "  min BlockCellLen: " << MinBlockCellLen << ", max BlockCellLen: " << MaxBlockCellLen << std::endl;
    if (stepinfo) {
      std::cout << "  BlockCellLen | StdDev" << std::endl;
    }
  }
  // iterate all possible blockcell length
  for (int BlockCellLen = MinBlockCellLen; BlockCellLen < MaxBlockCellLen + 1; ++BlockCellLen) {
    // 8. get stdDev
    T newStdDev = IterateAndOptimizeImp(ProcNum, BlockCellLen);
    // 9. update, for now we force the number of blocks to be exactly ProcNum
    if (newStdDev < stdDev && getAllBasicBlocks().size() == static_cast<std::size_t>(ProcNum)) {
      stdDev = newStdDev;
      bestBlockCellLen = BlockCellLen;
    }
    // print info
    if (stepinfo) {
      MPI_RANK(0)
      std::cout << "  " << std::setw(12) << BlockCellLen << " | " << newStdDev << std::endl;
    }
  }
  // perform the best optimization scheme
  IterateAndOptimizeImp(ProcNum, bestBlockCellLen);
  // print info
  MPI_RANK(0)
  std::size_t MaxBlockCellNum{};
  std::size_t MinBlockCellNum = std::numeric_limits<std::size_t>::max();
  std::size_t TotalCellNum{};
  for (const BasicBlock<T, 3> &block : getAllBasicBlocks()) {
    MaxBlockCellNum = std::max(MaxBlockCellNum, block.getN());
    MinBlockCellNum = std::min(MinBlockCellNum, block.getN());
    TotalCellNum += block.getN();
  }
  std::size_t AverageCellNum = static_cast<std::size_t>(TotalCellNum / getAllBasicBlocks().size());

  std::size_t OverlappedCellNum = ComputeOverlappedCellNum(getAllBasicBlocks());
  std::cout << "Final Optimization Result:\n"
            << "  using BlockCellLen:  " << bestBlockCellLen << " \n"
            << "  Total CellNum:       " << TotalCellNum << " \n"
            << "  Overlapped CellNum:  " << OverlappedCellNum << " \n"
            << "  Average CellNum:     " << AverageCellNum << " \n"
            << "  Max Block's CellNum: " << MaxBlockCellNum << " \n"
            << "  Min Block's CellNum: " << MinBlockCellNum << " \n"
            << "  StdDev:              " << stdDev << std::endl;
}

template <typename T>
T ComputeBlockNStdDev(const std::vector<BasicBlock<T, 3>> &Blocks, const GeometryFlagField<T>& flagF) {
  std::vector<std::size_t> Nvecs;
  T mean{};
  for (const BasicBlock<T, 3> &block : Blocks) {
    std::size_t N{};
    for (int z = 0; z < block.getNz(); ++z) {
      for (int y = 0; y < block.getNy(); ++y) {
        for (int x = 0; x < block.getNx(); ++x) {
          const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{x, y, z});
          if (flagF[vox] != flagF._VoidFlag) ++N;
        }
      }
    }
    Nvecs.push_back(N);
    mean += static_cast<T>(N);
  }
  mean /= Blocks.size();
  T stdDev{};
  for (std::size_t N : Nvecs) {
    stdDev += std::pow((static_cast<T>(N) / mean - T(1)), 2);
  }
  stdDev = std::sqrt(stdDev / Blocks.size());
  return stdDev;
}

template <typename T>
void buildSAT(int dir, const BasicBlock<T, 3>& block,
  const GeometryFlagField<T>& flagF, std::vector<std::size_t> &SAT) {
  SAT.clear();
  const int Nx = block.getNx();
  const int Ny = block.getNy();
  const int Nz = block.getNz();
  const std::uint8_t voidflag = flagF._VoidFlag;
  if (dir == 0) {
    SAT.reserve(Nx);
    for (int x = 0; x < Nx; ++x) {
      std::size_t sum{};
      for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
          const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{x, y, z});
          if (flagF[vox] != voidflag) ++sum; 
        }
      }
      SAT.push_back(sum);
    }
  } else if (dir == 1) {
    SAT.reserve(Ny);
    for (int y = 0; y < Ny; ++y) {
      std::size_t sum{};
      for (int z = 0; z < Nz; ++z) {
        for (int x = 0; x < Nx; ++x) {
          const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{x, y, z});
          if (flagF[vox] != voidflag) ++sum; 
        }
      }
      SAT.push_back(sum);
    }
  } else if (dir == 2) {
    SAT.reserve(Nz);
    for (int z = 0; z < Nz; ++z) {
      std::size_t sum{};
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{x, y, z});
          if (flagF[vox] != voidflag) ++sum; 
        }
      }
      SAT.push_back(sum);
    }
  }
}

template <typename T>
void cobisectImpl(BasicBlock<T, 3>& originBlock, std::vector<BasicBlock<T, 3>>& divideBlocks, 
  const GeometryFlagField<T>& flagF) {
  // 1. find the longest edge: 0 -> x, 1 -> y, 2 -> z
  int longestEdge = 2;
  NbrDirection dir = NbrDirection::ZN;
  int maxEdgeLen = originBlock.getNz();
  if (originBlock.getNy() > maxEdgeLen) {
    longestEdge = 1;
    dir = NbrDirection::YN;
    maxEdgeLen = originBlock.getNy();
  }
  if (originBlock.getNx() > maxEdgeLen) {
    longestEdge = 0;
    dir = NbrDirection::XN;
    maxEdgeLen = originBlock.getNx();
  }
  // 2. build the SAT of the longest edge from originBlock
  std::vector<std::size_t> SAT;
  buildSAT(longestEdge, originBlock, flagF, SAT);
  // 3. divide the block along the longest edge into 2 parts
  //    get the total number of cells and half total number of cells
  std::size_t total{};
  for (std::size_t i : SAT) total += i;
  std::size_t halftotal = total / 2;
  //   find the index where the sum of SAT is larger than halftotal
  std::size_t sum{};
  int idx{};
  for (std::size_t i = 0; i < SAT.size(); ++i) {
    sum += SAT[i];
    ++idx;
    if (sum > halftotal) break;
  }
  //    find which one is closer to halftotal
  if (idx > 0) {
    if (sum - halftotal > halftotal - sum + SAT[idx - 1]) --idx;
  }
  //    perform division
  int restidx = maxEdgeLen - idx;
  BasicBlock<T, 3> divideBlock1 = originBlock;
  BasicBlock<T, 3> divideBlock2 = originBlock;
  divideBlock1.resize(-restidx, getOpposite(dir));
  divideBlock2.resize(-idx, dir);
  //    add divideBlocks
  divideBlocks.push_back(divideBlock1);
  divideBlocks.push_back(divideBlock2);
}

template <typename T>
void cobisectImpl(int Edge, BasicBlock<T, 3>& originBlock, std::vector<BasicBlock<T, 3>>& divideBlocks, 
  const GeometryFlagField<T>& flagF) {
  // 1. find which edge
  NbrDirection dir = NbrDirection::NONE;
  int maxEdgeLen = originBlock.getMesh()[Edge];
  if (Edge == 0) {
    dir = NbrDirection::XN;
  } else if (Edge == 1) {
    dir = NbrDirection::YN;
  } else if (Edge == 2) {
    dir = NbrDirection::ZN;
  } else {
    IF_MPI_RANK(0) {std::cerr << "[cobisectImpl]: Edge must be 0, 1 or 2" << std::endl;}
    exit(1);
  }
  // 2. build the SAT of the longest edge from originBlock
  std::vector<std::size_t> SAT;
  buildSAT(Edge, originBlock, flagF, SAT);
  // 3. divide the block along the longest edge into 2 parts
  //    get the total number of cells and half total number of cells
  std::size_t total{};
  for (std::size_t i : SAT) total += i;
  std::size_t halftotal = total / 2;
  //   find the index where the sum of SAT is larger than halftotal
  std::size_t sum{};
  int idx{};
  for (std::size_t i = 0; i < SAT.size(); ++i) {
    sum += SAT[i];
    ++idx;
    if (sum > halftotal) break;
  }
  //    find which one is closer to halftotal
  if (idx > 0) {
    if (sum - halftotal > halftotal - sum + SAT[idx - 1]) --idx;
  }
  //    perform division
  int restidx = maxEdgeLen - idx;
  BasicBlock<T, 3> divideBlock1 = originBlock;
  BasicBlock<T, 3> divideBlock2 = originBlock;
  divideBlock1.resize(-restidx, getOpposite(dir));
  divideBlock2.resize(-idx, dir);
  //    add divideBlocks
  divideBlocks.push_back(divideBlock1);
  divideBlocks.push_back(divideBlock2);
}

template <typename T>
void cobisect(std::vector<BasicBlock<T, 3>>& originBlocks, std::vector<BasicBlock<T, 3>>& divideBlocks, 
  const GeometryFlagField<T>& flagF) {
  for(BasicBlock<T, 3>& block : originBlocks) {
    cobisectImpl(block, divideBlocks, flagF);
  }
}

template <typename T>
void cobisect(int Edge, std::vector<BasicBlock<T, 3>>& originBlocks, std::vector<BasicBlock<T, 3>>& divideBlocks, 
  const GeometryFlagField<T>& flagF) {
  for(BasicBlock<T, 3>& block : originBlocks) {
    cobisectImpl(Edge, block, divideBlocks, flagF);
  }
}

template <typename T>
void BlockGeometryHelper3D<T>::RCBOptimization(int ProcNum, bool verbose) {

  if (ProcNum == 1) {
    std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllBasicBlocks();
    BasicBlocks.clear();
    BasicBlocks.push_back(_BaseBlock);
    BasicBlocks[0].setBlockId(0);
    RemoveUnusedCells(std::uint8_t{1}, true);
    return;
  }
  // 1. check if the ProcNum is power of 2 
  //    and decide number of iterations
  int ProcNumx = 1;
  int iter{};
  for (;ProcNumx < ProcNum;) {
    ++iter;
    ProcNumx *= 2;
  }
  if (ProcNumx != ProcNum) {
    IF_MPI_RANK(0) {
      std::cerr << "[BlockGeometryHelper3D<T>::RCBOptimization]: ProcNum must be power of 2" << std::endl;
    }
    exit(1);
  }
  
  std::vector<BasicBlock<T, 3>> originBlocks;
  originBlocks.push_back(_BaseBlock);
  std::vector<BasicBlock<T, 3>> divideBlocks;
  // 2. recursively divide blocks
  for (int i = 0; i < iter; ++i) {
    divideBlocks.clear();
    cobisect(originBlocks, divideBlocks, _FlagField);
    //  update originBlocks and divideBlocks
    originBlocks.clear();
    originBlocks = divideBlocks;
  }
  // 3. add divideBlocks to Blocks
  std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllBasicBlocks();
  BasicBlocks.clear();
  BasicBlocks = divideBlocks;

  // 4. update block id
  for (std::size_t i = 0; i < BasicBlocks.size(); ++i) {
    BasicBlocks[i].setBlockId(i);
  }
  // 5. remove unused cells
  RemoveUnusedCells(std::uint8_t{1}, true);

  // 6. output info
  T stddev = ComputeBlockNStdDev(BasicBlocks, _FlagField);
  IF_MPI_RANK(0) {
    std::cout << "[BlockGeometryHelper3D<T>::RCBOptimization]:\n"
              << "  BlockNum: " << BasicBlocks.size() << "  with StdDev: " << stddev << std::endl;
  }
}

template <typename T>
void BlockGeometryHelper3D<T>::TRCBOptimization(int ProcNum, bool verbose) {
  // 0. check if the ProcNum is power of 2 
  //    and decide number of iterations
  int ProcNumx = 1;
  int iter{};
  for (;ProcNumx < ProcNum;) {
    ++iter;
    ProcNumx *= 2;
  }
  if (ProcNumx != ProcNum) {
    IF_MPI_RANK(0) {
      std::cerr << "[BlockGeometryHelper3D<T>::RCBOptimization]: ProcNum must be power of 2" << std::endl;
    }
    exit(1);
  }
  int titer = iter / 3;
  int riter = iter - titer * 3;
  
  std::vector<BasicBlock<T, 3>> originBlocks;
  originBlocks.push_back(_BaseBlock);
  std::vector<BasicBlock<T, 3>> divideBlocks;
  //    recursively divide blocks
  for (int i = 0; i < titer; ++i) {
    for (int j = 2; j >= 0; --j) {
      divideBlocks.clear();
      cobisect(j, originBlocks, divideBlocks, _FlagField);
      //  update originBlocks and divideBlocks
      originBlocks.clear();
      originBlocks = divideBlocks;
    }
  }
  for (int i = 0; i < riter; ++i) {
    divideBlocks.clear();
    cobisect(originBlocks, divideBlocks, _FlagField);
    //  update originBlocks and divideBlocks
    originBlocks.clear();
    originBlocks = divideBlocks;
  }
  // 4. add divideBlocks to Blocks
  std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllBasicBlocks();
  BasicBlocks.clear();
  BasicBlocks = divideBlocks;

  // 5. remove unused cells
  RemoveUnusedCells(std::uint8_t{1}, false);
  //    update block id
  for (std::size_t i = 0; i < BasicBlocks.size(); ++i) {
    BasicBlocks[i].setBlockId(i);
  }

  // 6. output info
  T stddev = ComputeBlockNStdDev(BasicBlocks, _FlagField);
  IF_MPI_RANK(0) {
    std::cout << "[BlockGeometryHelper3D<T>::RCBOptimization]:\n"
              << "  BlockNum: " << BasicBlocks.size() << "  with StdDev: " << stddev << std::endl;
  }
}


template <typename T>
void BlockGeometryHelper3D<T>::LoadOptimization(int maxiter, T tolstddev, bool outputinfo) {
  std::vector<BasicBlock<T, 3>> &originBlocks = getAllBasicBlocks();
  if (originBlocks.size() == 1) return;
  // copy Blocks to a buffer vector
  std::vector<BasicBlock<T, 3>> Blocks = originBlocks;

  // old info
  // total number of cells
  std::size_t OldTotalCellNum{};
  for (const BasicBlock<T, 3> &block : Blocks) {
    OldTotalCellNum += block.getN();
  }
  // std dev
  T OldStdDev = ComputeBlockNStdDev(Blocks);
  // overlapped cell number
  std::size_t OldOverlappedCellNum = ComputeOverlappedCellNum(Blocks);  
  // end old info

  // update block id, make sure block id is consistent with the index in Blocks
  for (std::size_t i = 0; i < Blocks.size(); ++i) {
    Blocks[i].setBlockId(i);
  }
  // get nbr info for each BasicBlocks, each block has a vector of FaceNbrInfo
  std::vector<std::vector<FaceNbrInfo>> NbrInfos;
  
  buildFaceNbrInfos(Blocks, NbrInfos);

  // best number of cells per block
  const std::size_t BestN = OldTotalCellNum / Blocks.size();

  // get extension of the whole domain
  Vector<T, 3> min = Blocks[0].getMin();
  Vector<T, 3> max = Blocks[0].getMax();
  for (const BasicBlock<T, 3> &block : Blocks) {
    for (int i = 0; i < 3; ++i) {
      min[i] = std::min(min[i], block.getMin()[i]);
      max[i] = std::max(max[i], block.getMax()[i]);
    }
  }
  Vector<T, 3> ext = max - min;

  // get number of cells in each direction
  const T voxsize = Blocks[0].getVoxelSize();
  int TotalNx = static_cast<int>(std::round(ext[0] / voxsize));
  int TotalNy = static_cast<int>(std::round(ext[1] / voxsize));
  int TotalNz = static_cast<int>(std::round(ext[2] / voxsize));

  // let TotalNx: TotalNy: TotalNz = 1: a: b
  T a = static_cast<T>(TotalNy) / static_cast<T>(TotalNx);
  T b = static_cast<T>(TotalNz) / static_cast<T>(TotalNx);
  // BestNx: BestNy: BestNz = 1: a: b
  // BestNx*BestNy*BestNz = BestN
  // a*b*BestNx^3 = BestN
  int BestNx = static_cast<int>(std::round(std::pow(BestN / (a * b), 1.0 / 3.0)));
  int BestNy = static_cast<int>(std::round(a * BestNx));
  int BestNz = static_cast<int>(std::round(b * BestNx));
  int BestNxNy = BestNx * BestNy;
  int BestNxNz = BestNx * BestNz;
  int BestNyNz = BestNy * BestNz;
  // BestNx*BestNy*BestNz may not be exactly equal to BestN
  // adjust to make them as close as possible
  int diff = BestNx * BestNy * BestNz - BestN;
  while (diff != 0) {
    if (diff > 0) {
      if (diff >= BestNxNy) {
        BestNz -= 1;
      } else if (diff >= BestNxNz) {
        BestNy -= 1;
      } else if (diff >= BestNyNz) {
        BestNx -= 1;
      } else {
        // BestNx*BestNy*BestNz is already close to BestN
        // find if they could be closer
        int minface = std::min(std::min(BestNxNy, BestNxNz), BestNyNz);
        int newdiff = diff - minface;
        if (std::abs(newdiff) < diff) {
          // use newdiff
          if (BestNxNy <= BestNxNz && BestNxNy <= BestNyNz) {
            BestNz -= 1;
          } else if (BestNxNz <= BestNxNy && BestNxNz <= BestNyNz) {
            BestNy -= 1;
          } else {
            BestNx -= 1;
          }
        } else {
          // this is already the best
          break;
        }
      }
    } else if (diff < 0) {
      int _diff = std::abs(diff);
      if (_diff >= BestNxNy) {
        BestNz += 1;
      } else if (_diff >= BestNxNz) {
        BestNy += 1;
      } else if (_diff >= BestNyNz) {
        BestNx += 1;
      } else {
        // BestNx*BestNy*BestNz is already close to BestN
        // find if they could be closer
        int minface = std::min(std::min(BestNxNy, BestNxNz), BestNyNz);
        int newdiff = diff + minface;
        if (newdiff < _diff) {
          // use newdiff
          if (BestNxNy <= BestNxNz && BestNxNy <= BestNyNz) {
            BestNz += 1;
          } else if (BestNxNz <= BestNxNy && BestNxNz <= BestNyNz) {
            BestNy += 1;
          } else {
            BestNx += 1;
          }
        } else {
          // this is already the best
          break;
        }
      }
    }
    diff = BestNx * BestNy * BestNz - BestN;
  }

  // perform optimization
  int iter{};
  T min_stddev = OldStdDev;
  // only one layer of overlapped cells are considered
  std::size_t overlappedCellNum = OldOverlappedCellNum;

  // we have 2 optimization guidelines:
  // 1. size first: shrink or expand blocks to their cell number be as close as possible to BestN
  //    which is good for load balancing
  // 2. scale first: shrink or expand blocks to make their Nx, Ny, Nz be as close as possible to BestNx, BestNy, BestNz
  //    which is good for efficient communication(least number of overlapped cells)
  // we will try both guidelines and find the shceme with the minimum stddev 

  // we used overlappedCellNum here, and in the future we may design an algorithm to 
  // include both stddev and overlappedCellNum to choose the best scheme

  // 1. size first optimization
  while (iter < maxiter && min_stddev > tolstddev) {

    for (BasicBlock<T, 3> &block : Blocks) {
      // get block id
      std::size_t id = block.getBlockId();
      std::size_t N = block.getN();
      // int Nx = block.getNx();
      // int Ny = block.getNy();
      // int Nz = block.getNz();
      // get block's neighbors
      std::vector<FaceNbrInfo> &NbrInfo = NbrInfos[id];
      if (NbrInfo.empty()) continue;
      // get all neighbor types and determine movable faces/ resizable directions
      // movable faces: XN, XP| YN, YP| ZN, ZP
      // resizable directions: X| Y| Z
      // get all neighbor directions
      // std::uint8_t AllNbrDir{};
      // for (const FaceNbrInfo &info : NbrInfo) {
      //   AllNbrDir |= static_cast<std::uint8_t>(info.Direction);
      // }

      BasicBlock<T, 3> *nblockptr = nullptr;
      // shrink or expand
      if (N > BestN) {
        // find neighbors with the minimum number of cells
        std::size_t min_nblockN = Blocks[NbrInfo[0].NbrBlockId].getN();
        FaceNbrInfo* minN_info = &NbrInfo[0];
        for (FaceNbrInfo &info : NbrInfo) {
          if (!info.Movable) continue;
          if (Blocks[info.NbrBlockId].getN() < min_nblockN) {
            min_nblockN = Blocks[info.NbrBlockId].getN();
            minN_info = &info;
          }
        }
        if (!minN_info->Movable) continue;
        nblockptr = &Blocks[minN_info->NbrBlockId];
        // base and neighbors are already large enough
        // if (min_nblockN >= BestN) continue;
        // shrink base and expand neighbors
        // AdjustAdjacentBlockSize(block, *nblockptr, BestN, minN_info->Direction);
        // AdjustAdjacentBlockSize(block, *nblockptr, minN_info->Direction);
        AdjustAdjacentBlockSize_byStep(block, *nblockptr, BestN, minN_info->Direction);
      } else if (N < BestN) {
        // find neighbors with the maximum number of cells
        std::size_t max_nblockN = Blocks[NbrInfo[0].NbrBlockId].getN();
        FaceNbrInfo* maxN_info = &NbrInfo[0];
        for (FaceNbrInfo &info : NbrInfo) {
          if (!info.Movable) continue;
          if (Blocks[info.NbrBlockId].getN() > max_nblockN) {
            max_nblockN = Blocks[info.NbrBlockId].getN();
            maxN_info = &info;
          }
        }
        if (!maxN_info->Movable) continue;
        nblockptr = &Blocks[maxN_info->NbrBlockId];
        // base and neighbors are already small enough
        // if (max_nblockN <= BestN) continue;
        // expand base and shrink neighbors
        // AdjustAdjacentBlockSize(block, *nblockptr, BestN, maxN_info->Direction);
        // AdjustAdjacentBlockSize(block, *nblockptr, maxN_info->Direction);
        AdjustAdjacentBlockSize_byStep(block, *nblockptr, BestN, maxN_info->Direction);
      } else if (N == BestN) {
        // do nothing
        continue;
      }
      // end shrink or expand
      // nbr info may be changed, update it
      buildFaceNbrInfos(Blocks, NbrInfos);

    }

    ++iter;

    T stddev = ComputeBlockNStdDev(Blocks);
    if (stddev < min_stddev) {
      min_stddev = stddev;
      overlappedCellNum = ComputeOverlappedCellNum(Blocks);
      // copy Blocks to originBlocks
      originBlocks = Blocks;
    }
    
  }

  // 2. scale first optimization


  // print info
  if (!outputinfo) return;
  MPI_RANK(0)
  std::size_t TotalCellNum{};
  for (const BasicBlock<T, 3> &block : originBlocks) {
    TotalCellNum += block.getN();
  }
  std::cout << "[BlockGeometryHelper3D<T>::LoadOptimization]:\n"
            << "  After              " << iter << " iterations\n"
            << "  Total CellNum:     " << OldTotalCellNum << " -> " << TotalCellNum << "\n"
            << "  StdDev:            " << OldStdDev << " -> " << min_stddev << "\n"
            << "  OverlappedCellNum: " << OldOverlappedCellNum << " -> " << overlappedCellNum << std::endl;
}


template <typename T>
void BlockGeometryHelper3D<T>::Init(int blockcelllen) {
  if (blockcelllen < 4) {
    std::cerr << "BlockGeometryHelper3D<T>, BlockCellLen < 4" << std::endl;
  }
  // clear all vectors
  _BasicBlocks0.clear();
  _BasicBlocks1.clear();
  static_cast<BasicBlock<T, 3>&>(*this) = BasicBlock<T, 3>(
    _Reader->getVoxelSize(), 
    _Reader->getMesh().getAABB().getExtended((_Overlap+_Ext) * _Reader->getVoxelSize()),
    AABB<int, 3>(
      Vector<int, 3>{}, 
      Vector<int, 3>{
        int(std::ceil(_Reader->getMesh().getMax_Min()[0] / _Reader->getVoxelSize())) + 2*(_Overlap+_Ext) - 1, 
        int(std::ceil(_Reader->getMesh().getMax_Min()[1] / _Reader->getVoxelSize())) + 2*(_Overlap+_Ext) - 1,
        int(std::ceil(_Reader->getMesh().getMax_Min()[2] / _Reader->getVoxelSize())) + 2*(_Overlap+_Ext) - 1}));
  _BaseBlock = BasicBlock<T, 3>(
    _Reader->getVoxelSize(), 
    _Reader->getMesh().getAABB().getExtended(_Ext*_Reader->getVoxelSize()),
    AABB<int, 3>(
      Vector<int, 3>{_Overlap}, 
      Vector<int, 3>{
        int(std::ceil(_Reader->getMesh().getMax_Min()[0] / _Reader->getVoxelSize())) + (_Overlap+2*_Ext) - 1, 
        int(std::ceil(_Reader->getMesh().getMax_Min()[1] / _Reader->getVoxelSize())) + (_Overlap+2*_Ext) - 1,
        int(std::ceil(_Reader->getMesh().getMax_Min()[2] / _Reader->getVoxelSize())) + (_Overlap+2*_Ext) - 1}));
  
                
  BlockCellLen = blockcelllen;
  _Exchanged = true;
  _IndexExchanged = true;
  
  CellsNx = std::ceil(T(_BaseBlock.getNx()) / T(BlockCellLen));
  CellsNy = std::ceil(T(_BaseBlock.getNy()) / T(BlockCellLen));
  CellsNz = std::ceil(T(_BaseBlock.getNz()) / T(BlockCellLen));
  CellsN = CellsNx * CellsNy * CellsNz;

  // correct the mesh size by CellsNx, CellsNy, CellsNz
  int NewNx = CellsNx * BlockCellLen;
  int NewNy = CellsNy * BlockCellLen;
  int NewNz = CellsNz * BlockCellLen;
  T MeshSizeX = NewNx * _BaseBlock.getVoxelSize();
  T MeshSizeY = NewNy * _BaseBlock.getVoxelSize();
  T MeshSizeZ = NewNz * _BaseBlock.getVoxelSize();

  // _Ext is already included in the _BaseBlock, CellsNx, CellsNy, CellsNz
  static_cast<BasicBlock<T, 3>&>(*this) = BasicBlock<T, 3>(
    _Reader->getVoxelSize(), 
    AABB<T, 3>(_Reader->getMesh().getMin() - (_Overlap+_Ext)*_Reader->getVoxelSize(), 
               _Reader->getMesh().getMin() - (_Overlap+_Ext)*_Reader->getVoxelSize() + Vector<T, 3>{MeshSizeX, MeshSizeY, MeshSizeZ} + 2*_Overlap*_Reader->getVoxelSize()),
    AABB<int, 3>(Vector<int, 3>{}, 
                 Vector<int, 3>{NewNx + 2*_Overlap-1, NewNy + 2*_Overlap-1, NewNz + 2*_Overlap-1}));
  _BaseBlock = BasicBlock<T, 3>(
    _Reader->getVoxelSize(), 
    AABB<T, 3>(_Reader->getMesh().getMin() - _Ext*_Reader->getVoxelSize(), 
               _Reader->getMesh().getMin() - _Ext*_Reader->getVoxelSize() + Vector<T, 3>{MeshSizeX, MeshSizeY, MeshSizeZ}),
    AABB<int, 3>(Vector<int, 3>{_Overlap}, 
                 Vector<int, 3>{NewNx + _Overlap-1, NewNy + _Overlap-1, NewNz + _Overlap-1}));

  // end correct

  Vector<int, 3> Projection{1, CellsNx, CellsNx * CellsNy};

  Delta_Cellidx = make_Array<int, D3Q27<T>::q - 1>(
    [&](int i) { return D3Q27<T>::c[i + 1] * Projection; });

  CreateBlockCells();
  TagBlockCells();
}


template <typename T>
void BlockGeometryHelper3D<T>::LoadBalancing(int ProcessNum) {
  std::vector<std::vector<int>> &BlockIndex = getAllOldBlockIndices();
  BlockIndex.clear();
  std::vector<std::vector<BasicBlock<T, 3> *>> &BlockIndexPtr = getAllOldBlockIndexPtrs();
  BlockIndexPtr.clear();
  _IndexExchanged = !_IndexExchanged;
#ifdef MPI_ENABLED
  BlockIndex.resize(ProcessNum);
  // map block's cell num to block id
  std::vector<std::pair<std::size_t, int>> BlockId_CellNum;
  for (BasicBlock<T, 3> &block : getAllBasicBlocks()) {
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
  for (BasicBlock<T, 3> &block : getAllBasicBlocks()) {
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
    std::cout << i << " | ";
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
void BlockGeometryHelper3D<T>::SetupMPINbrs() {
  _MPIBlockNbrs.clear();
  _MPIBlockNbrs.resize(getAllBasicBlocks().size(), std::vector<std::pair<int, int>>{});
  for (std::size_t iRank = 0; iRank < getAllBlockIndices().size(); ++iRank) {
    for (int blockid : getAllBlockIndices()[iRank]) {
      std::vector<std::pair<int, int>> &nbrsvec = _MPIBlockNbrs[blockid];
      const BasicBlock<T, 3> block = getAllBasicBlock(static_cast<std::size_t>(blockid)).getExtBlock(1);
      // for all blocks in all ranks except iRank
      for (std::size_t nRank = 0; nRank < getAllBlockIndices().size(); ++nRank) {
        if (nRank != iRank) {
          for (int nblockid : getAllBlockIndices()[nRank]) {
            // like SetupNbrs() in BlockGeometry3D, we use baseblock here, NO need to use extblock
            const BasicBlock<T, 3> nbaseblock = getAllBasicBlock(static_cast<std::size_t>(nblockid));
            if (isOverlapped(block, nbaseblock))
              nbrsvec.push_back(std::make_pair(nRank, nblockid));
          }
        }
      }
    }
  }
}

template <typename T>
void BlockGeometryHelper3D<T>::InitBlockGeometry3D(bool useHelperOlap) {
  _BlockGeometry3D = std::make_unique<BlockGeometry3D<T>>(*this, getAllBasicBlocks(), useHelperOlap);
}

template <typename T>
int BlockGeometryHelper3D<T>::whichRank(int blockid) {
  for (std::size_t iRank = 0; iRank < getAllBlockIndices().size(); ++iRank) {
    for (int id : getAllBlockIndices()[iRank]) {
      if (id == blockid) return static_cast<int>(iRank);
    }
  }
  return -1;
}

#endif

template <typename T>
void BlockGeometryHelper3D<T>::TagRefineLayer(std::vector<std::uint8_t> &refine, bool &refined) {
  UpdateMaxLevel();
  // refine one additional layer if has neighbor with lower level
  // tag cells to be refined
  int XY = CellsNx * CellsNy;
  for (std::uint8_t level = _MaxLevel; level > std::uint8_t(0); --level) {
    for (int cellz = 1; cellz < CellsNz - 1; ++cellz) {
      for (int celly = 1; celly < CellsNy - 1; ++celly) {
        for (int cellx = 1; cellx < CellsNx - 1; ++cellx) {
          int cellid = celly * CellsNx + cellx + cellz * XY;
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
}

template <typename T>
void BlockGeometryHelper3D<T>::CheckRefine() {
  Refine();
  UpdateMaxLevel();
  // refine one additional layer if has neighbor with 2 lower level
  // tag cells to be refined
  int XY = CellsNx * CellsNy;
  for (std::uint8_t level = _MaxLevel; level > std::uint8_t(1); --level) {
    for (int cellz = 1; cellz < CellsNz - 1; ++cellz) {
      for (int celly = 1; celly < CellsNy - 1; ++celly) {
        for (int cellx = 1; cellx < CellsNx - 1; ++cellx) {
          int cellid = celly * CellsNx + cellx + cellz * XY;
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
    }
    Refine();
  }
}

template <typename T>
void BlockGeometryHelper3D<T>::Refine() {
  for (int i = 0; i < CellsN; ++i) {
    if (util::isFlag(static_cast<CellTagType>(_BlockCellTags[i]), static_cast<CellTagType>(BlockCellTag::refine))) {
      _BlockCells[i].refine();
      // reset tag
      _BlockCellTags[i] = BlockCellTag::none;
    }
  }
}

template <typename T>
void BlockGeometryHelper3D<T>::RemoveUnusedCells(const StlReader<T>& reader, bool outputinfo) {

  Octree<T>* tree = reader.getTree();
  std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllBasicBlocks();

  std::size_t Total{};
  for (const BasicBlock<T, 3> &block : BasicBlocks) {
    Total += block.getN();
  }

  // statistics
  std::size_t totalremoved{};

  for (BasicBlock<T, 3> &block : BasicBlocks) {
    // if contains cell outside of the geometry, test if it could be shrinked
    if (hasOutSideCell(tree, block)) {
      // 1. XN face, shrink along +x direction
      int xnlayer{};
      for (int ix = 0; ix < block.getNx() - _Overlap; ++ix) {
        // find if all cells on YZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int iy = 0; iy < block.getNy(); ++iy) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto xncheckalloutside;
              }
            }
          }
        }
        xncheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++xnlayer;
        }
      }
      xnlayer -= _Ext;
      if (xnlayer > 0) {
        totalremoved += xnlayer * block.getNy() * block.getNz();
        block.resize(-xnlayer, NbrDirection::XN);
      }

      // 2. XP face, shrink along -x direction
      int xplayer{};
      for (int ix = block.getNx() - 1; ix >= _Overlap; --ix) {
        // find if all cells on YZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int iy = 0; iy < block.getNy(); ++iy) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto xpcheckalloutside;
              }
            }
          }
        }
        xpcheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++xplayer;
        }
      }
      xplayer -= _Ext;
      if (xplayer > 0) {
        totalremoved += xplayer * block.getNy() * block.getNz();
        block.resize(-xplayer, NbrDirection::XP);
      }

      // 3. YN face, shrink along +y direction
      int ynlayer{};
      for (int iy = 0; iy < block.getNy() - _Overlap; ++iy) {
        // find if all cells on XZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto yncheckalloutside;
              }
            }
          }
        }
        yncheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++ynlayer;
        }
      }
      ynlayer -= _Ext;
      if (ynlayer > 0) {
        totalremoved += ynlayer * block.getNx() * block.getNz();
        block.resize(-ynlayer, NbrDirection::YN);
      }

      // 4. YP face, shrink along -y direction
      int yplayer{};
      for (int iy = block.getNy() - 1; iy >= _Overlap; --iy) {
        // find if all cells on XZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto ypcheckalloutside;
              }
            }
          }
        }
        ypcheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++yplayer;
        }
      }
      yplayer -= _Ext;
      if (yplayer > 0) {
        totalremoved += yplayer * block.getNx() * block.getNz();
        block.resize(-yplayer, NbrDirection::YP);
      }

      // 5. ZN face, shrink along +z direction
      int znlayer{};
      for (int iz = 0; iz < block.getNz() - _Overlap; ++iz) {
        // find if all cells on XY plane are OUTSIDE
        bool alloutside = true;
        for (int iy = 0; iy < block.getNy(); ++iy) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto zncheckalloutside;
              }
            }
          }
        }
        zncheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++znlayer;
        }
      }
      znlayer -= _Ext;
      if (znlayer > 0) {
        totalremoved += znlayer * block.getNx() * block.getNy();
        block.resize(-znlayer, NbrDirection::ZN);
      }

      // 6. ZP face, shrink along -z direction
      int zplayer{};
      for (int iz = block.getNz() - 1; iz >= _Overlap; --iz) {
        // find if all cells on XY plane are OUTSIDE
        bool alloutside = true;
        for (int iy = 0; iy < block.getNy(); ++iy) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto zpcheckalloutside;
              }
            }
          }
        }
        zpcheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++zplayer;
        }
      }
      zplayer -= _Ext;
      if (zplayer > 0) {
        totalremoved += zplayer * block.getNx() * block.getNy();
        block.resize(-zplayer, NbrDirection::ZP);
      }

    }
  }

  std::size_t TotalR{};
  for (const BasicBlock<T, 3> &block : BasicBlocks) {
    TotalR += block.getN();
  }

  // print statistics
  if (!outputinfo) return;
  MPI_RANK(0)
  std::cout << "[BlockGeometryHelper3D<T>::RemoveUnusedCells]: \n" 
  << "  Removed: " << totalremoved << " cells\n"
  << "  Total CellNum: " << Total << " -> " << TotalR << " cells\n";
}

template <typename T>
void BlockGeometryHelper3D<T>::RemoveUnusedCells(std::uint8_t voidflag, bool outputinfo) {

  std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllBasicBlocks();

  std::size_t Total{};
  for (const BasicBlock<T, 3> &block : BasicBlocks) {
    Total += block.getN();
  }

  // statistics
  std::size_t totalremoved{};

  for (BasicBlock<T, 3> &block : BasicBlocks) {

    // if contains cell outside of the geometry, test if it could be shrinked
    bool has_void = false;
    for (int z = 0; z < block.getNz(); ++z) {
      for (int y = 0; y < block.getNy(); ++y) {
        for (int x = 0; x < block.getNx(); ++x) {
          const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{x, y, z});
          if (_FlagField[vox] == voidflag) {
            has_void = true;
            goto check_has_void;
          }
        }
      }
    }
    check_has_void:
    if (has_void) {
      // 1. XN face, shrink along +x direction
      int xnlayer{};
      for (int ix = 0; ix < block.getNx() - _Overlap; ++ix) {
        // find if all cells on YZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int iy = 0; iy < block.getNy(); ++iy) {
            const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{ix, iy, iz});
            if (_FlagField[vox] != voidflag) {
              alloutside = false;
              goto xncheckalloutside;
            }
          }
        }
        xncheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++xnlayer;
        }
      }
      xnlayer -= _Ext;
      if (xnlayer > 0) {
        totalremoved += xnlayer * block.getNy() * block.getNz();
        block.resize(-xnlayer, NbrDirection::XN);
      }

      // 2. XP face, shrink along -x direction
      int xplayer{};
      for (int ix = block.getNx() - 1; ix >= _Overlap; --ix) {
        // find if all cells on YZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int iy = 0; iy < block.getNy(); ++iy) {
            const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{ix, iy, iz});
            if (_FlagField[vox] != voidflag) {
              alloutside = false;
              goto xpcheckalloutside;
            }
          }
        }
        xpcheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++xplayer;
        }
      }
      xplayer -= _Ext;
      if (xplayer > 0) {
        totalremoved += xplayer * block.getNy() * block.getNz();
        block.resize(-xplayer, NbrDirection::XP);
      }

      // 3. YN face, shrink along +y direction
      int ynlayer{};
      for (int iy = 0; iy < block.getNy() - _Overlap; ++iy) {
        // find if all cells on XZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{ix, iy, iz});
            if (_FlagField[vox] != voidflag) {
              alloutside = false;
              goto yncheckalloutside;
            }
          }
        }
        yncheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++ynlayer;
        }
      }
      ynlayer -= _Ext;
      if (ynlayer > 0) {
        totalremoved += ynlayer * block.getNx() * block.getNz();
        block.resize(-ynlayer, NbrDirection::YN);
      }

      // 4. YP face, shrink along -y direction
      int yplayer{};
      for (int iy = block.getNy() - 1; iy >= _Overlap; --iy) {
        // find if all cells on XZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{ix, iy, iz});
            if (_FlagField[vox] != voidflag) {
              alloutside = false;
              goto ypcheckalloutside;
            }
          }
        }
        ypcheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++yplayer;
        }
      }
      yplayer -= _Ext;
      if (yplayer > 0) {
        totalremoved += yplayer * block.getNx() * block.getNz();
        block.resize(-yplayer, NbrDirection::YP);
      }

      // 5. ZN face, shrink along +z direction
      int znlayer{};
      for (int iz = 0; iz < block.getNz() - _Overlap; ++iz) {
        // find if all cells on XY plane are OUTSIDE
        bool alloutside = true;
        for (int iy = 0; iy < block.getNy(); ++iy) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{ix, iy, iz});
            if (_FlagField[vox] != voidflag) {
              alloutside = false;
              goto zncheckalloutside;
            }
          }
        }
        zncheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++znlayer;
        }
      }
      znlayer -= _Ext;
      if (znlayer > 0) {
        totalremoved += znlayer * block.getNx() * block.getNy();
        block.resize(-znlayer, NbrDirection::ZN);
      }

      // 6. ZP face, shrink along -z direction
      int zplayer{};
      for (int iz = block.getNz() - 1; iz >= _Overlap; --iz) {
        // find if all cells on XY plane are OUTSIDE
        bool alloutside = true;
        for (int iy = 0; iy < block.getNy(); ++iy) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<T, 3> vox = block.getVoxel(Vector<int, 3>{ix, iy, iz});
            if (_FlagField[vox] != voidflag) {
              alloutside = false;
              goto zpcheckalloutside;
            }
          }
        }
        zpcheckalloutside:
        if (!alloutside) {
          break;
        } else {
          ++zplayer;
        }
      }
      zplayer -= _Ext;
      if (zplayer > 0) {
        totalremoved += zplayer * block.getNx() * block.getNy();
        block.resize(-zplayer, NbrDirection::ZP);
      }
    }

  }

  std::size_t TotalR{};
  for (const BasicBlock<T, 3> &block : BasicBlocks) {
    TotalR += block.getN();
  }

  // print statistics
  if (!outputinfo) return;
  MPI_RANK(0)
  std::cout << "[BlockGeometryHelper3D<T>::RemoveUnusedCells]: \n" 
  << "  Removed: " << totalremoved << " cells\n"
  << "  Total CellNum: " << Total << " -> " << TotalR << " cells\n";
}

template <typename T>
void BlockGeometryHelper3D<T>::AddVoidCellLayer(const StlReader<T>& reader, bool outputinfo) {
  if (_Ext == 0) return;

  Octree<T>* tree = reader.getTree();
  std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllBasicBlocks();
  if (BasicBlocks.size() == 1) return;
  // build face neighbor infos
  std::vector<std::vector<FaceNbrInfo>> NbrInfos;
  // update block id
  for (std::size_t i = 0; i < BasicBlocks.size(); ++i) BasicBlocks[i].setBlockId(i);
  buildFaceNbrInfos(BasicBlocks, NbrInfos);

  std::size_t Total{};
  for (const BasicBlock<T, 3> &block : BasicBlocks) {
    Total += block.getN();
  }

  // statistics
  std::size_t totaladded{};

  for (BasicBlock<T, 3> &block : BasicBlocks) {
    const std::vector<FaceNbrInfo> &NbrInfo = NbrInfos[block.getBlockId()];

    std::vector<NbrDirection> dirs{NbrDirection::XN, NbrDirection::XP, 
      NbrDirection::YN, NbrDirection::YP, NbrDirection::ZN, NbrDirection::ZP};
    // remove direction that has neighbor
    for (const FaceNbrInfo &info : NbrInfo) {
      dirs.erase(std::remove(dirs.begin(), dirs.end(), info.Direction), dirs.end());
    }
    if (dirs.empty()) continue;

    for (NbrDirection dir : dirs) {
      // 1. XN face, expand along -x direction
      if (dir == NbrDirection::XN) {
        int ix = 0;
        // find if all cells on YZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int iy = 0; iy < block.getNy(); ++iy) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto xncheckalloutside;
              }
            }
          }
        }
        xncheckalloutside:
        if (!alloutside) {
          // has inside cell, expand one layer
          block.resize(1, NbrDirection::XN);
          totaladded += block.getNy() * block.getNz();
        }
      }
      // 2. XP face, expand along +x direction
      if (dir == NbrDirection::XP) {
        int ix = block.getNx() - 1;
        // find if all cells on YZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int iy = 0; iy < block.getNy(); ++iy) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto xpcheckalloutside;
              }
            }
          }
        }
        xpcheckalloutside:
        if (!alloutside) {
          // has inside cell, expand one layer
          block.resize(1, NbrDirection::XP);
          totaladded += block.getNy() * block.getNz();
        }
      }
      // 3. YN face, expand along -y direction
      if (dir == NbrDirection::YN) {
        int iy = 0;
        // find if all cells on XZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto yncheckalloutside;
              }
            }
          }
        }
        yncheckalloutside:
        if (!alloutside) {
          // has inside cell, expand one layer
          block.resize(1, NbrDirection::YN);
          totaladded += block.getNx() * block.getNz();
        }
      }
      // 4. YP face, expand along +y direction
      if (dir == NbrDirection::YP) {
        int iy = block.getNy() - 1;
        // find if all cells on XZ plane are OUTSIDE
        bool alloutside = true;
        for (int iz = 0; iz < block.getNz(); ++iz) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto ypcheckalloutside;
              }
            }
          }
        }
        ypcheckalloutside:
        if (!alloutside) {
          // has inside cell, expand one layer
          block.resize(1, NbrDirection::YP);
          totaladded += block.getNx() * block.getNz();
        }
      }
      // 5. ZN face, expand along -z direction
      if (dir == NbrDirection::ZN) {
        int iz = 0;
        // find if all cells on XY plane are OUTSIDE
        bool alloutside = true;
        for (int iy = 0; iy < block.getNy(); ++iy) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto zncheckalloutside;
              }
            }
          }
        }
        zncheckalloutside:
        if (!alloutside) {
          // has inside cell, expand one layer
          block.resize(1, NbrDirection::ZN);
          totaladded += block.getNx() * block.getNy();
        }
      }
      // 6. ZP face, expand along +z direction
      if (dir == NbrDirection::ZP) {
        int iz = block.getNz() - 1;
        // find if all cells on XY plane are OUTSIDE
        bool alloutside = true;
        for (int iy = 0; iy < block.getNy(); ++iy) {
          for (int ix = 0; ix < block.getNx(); ++ix) {
            const Vector<int, 3> locidx{ix, iy, iz};
            const Vector<T, 3> vox = block.getVoxel(locidx);
            Octree<T>* node = tree->find(vox);
            if (node != nullptr) {
              if (node->isLeaf() && node->getInside()) {
                alloutside = false;
                goto zpcheckalloutside;
              }
            }
          }
        }
        zpcheckalloutside:
        if (!alloutside) {
          // has inside cell, expand one layer
          block.resize(1, NbrDirection::ZP);
          totaladded += block.getNx() * block.getNy();
        }
      }
    }

  }

  std::size_t TotalR{};
  for (const BasicBlock<T, 3> &block : BasicBlocks) {
    TotalR += block.getN();
  }

  // print statistics
  if (!outputinfo) return;
  MPI_RANK(0)
  std::cout << "[BlockGeometryHelper3D<T>::AddVoidCellLayer]: \n" 
  << "  Added: " << totaladded << " cells\n"
  << "  Total CellNum: " << Total << " -> " << TotalR << " cells\n";
}