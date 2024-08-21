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
    _BaseBlock(baseblock), _overlap(olap)
       {}


template <typename T>
Block3D<T>::Block3D(const AABB<T, 3> &block, const AABB<int, 3> &idxblock, int blockid,
                    T voxelSize, int olap)
    : BasicBlock<T, 3>(voxelSize, block.getExtended(Vector<T, 3>{voxelSize}),
                       idxblock.getExtended(Vector<int, 3>{1}), blockid),
      _BaseBlock(voxelSize, block, idxblock, blockid), _overlap(olap) {
  // read from AABBs
  // ReadAABBs(AABBs, AABBflag);
}

template <typename T>
template <typename FieldType, typename LatSet>
void Block3D<T>::SetupBoundary(const AABB<T, 3> &block, FieldType &field,
                               typename FieldType::value_type bdvalue) {
  // temp flag field store the transition flag
  GenericArray<bool> TransFlag(BasicBlock<T, 3>::N, false);
  for (int z = _overlap; z < BasicBlock<T, 3>::Mesh[2] - _overlap; ++z) {
    for (int y = _overlap; y < BasicBlock<T, 3>::Mesh[1] - _overlap; ++y) {
      for (int x = _overlap; x < BasicBlock<T, 3>::Mesh[0] - _overlap; ++x) {
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

// -----------blockgeometry3d----------------


template <typename T>
BlockGeometry3D<T>::BlockGeometry3D(int Nx, int Ny, int Nz, int blocknum,
                                    const AABB<T, 3> &block, T voxelSize, int overlap)
    : BasicBlock<T, 3>(
        voxelSize, block.getExtended(Vector<T, 3>{voxelSize}),
        AABB<int, 3>(Vector<int, 3>{0}, Vector<int, 3>{Nx + 1, Ny + 1, Nz + 1})),
      _BaseBlock(voxelSize, block,
                 AABB<int, 3>(Vector<int, 3>{1}, Vector<int, 3>{Nx, Ny, Nz})),
      _overlap(overlap), _MaxLevel(std::uint8_t(0)) {
  CreateBlocks(blocknum);
  SetupNbrs();
  InitComm();
#ifndef MPI_ENABLED
  PrintInfo();
#endif
}

template <typename T>
BlockGeometry3D<T>::BlockGeometry3D(BlockGeometryHelper3D<T> &GeoHelper)
    : BasicBlock<T, 3>(GeoHelper), _BaseBlock(GeoHelper.getBaseBlock()), 
      _overlap(GeoHelper.getExt()), _MaxLevel(GeoHelper.getMaxLevel()) {
  // create blocks from GeoHelper
  for (BasicBlock<T, 3> *baseblock : GeoHelper.getBasicBlocks()) {
    int overlap = (baseblock->getLevel() != std::uint8_t(0)) ? 2 : 1;
    _Blocks.emplace_back(*baseblock, overlap);
    _BasicBlocks.emplace_back(*baseblock);
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
BlockGeometry3D<T>::BlockGeometry3D(const BlockReader3D<T>& blockreader) 
    : BasicBlock<T, 3>(blockreader.getBasicBlock()), _BaseBlock(blockreader.getBaseBlock()), 
      _overlap(1), _MaxLevel(blockreader.getMaxLevel()) {
  // create blocks from Block Reader
  for (const BasicBlock<T, 3> &baseblock : blockreader.getBlocks()) {
    int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
    _Blocks.emplace_back(baseblock, overlap);
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
void BlockGeometry3D<T>::PrintInfo() const {
  std::cout << "[BlockGeometry3D]: " << "Total Cell Num: " << getTotalCellNum()
            << std::endl;
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
  SetupNbrs();
  InitAllComm();
#ifdef MPI_ENABLED
  InitAllMPIComm(GeoHelper);
#else
  PrintInfo();
#endif
}

template <typename T>
std::size_t BlockGeometry3D<T>::getTotalCellNum() const {
  std::size_t sum = 0;
  for (const Block3D<T> &block : _Blocks) sum += block.getN();
  return sum;
}

template <typename T>
std::size_t BlockGeometry3D<T>::getBaseCellNum() const {
  std::size_t sum = 0;
  for (const Block3D<T> &block : _Blocks) sum += block.getBaseBlock().getN();
  return sum;
}

template <typename T>
void BlockGeometry3D<T>::DivideBlocks(int blocknum) {
  _BlockAABBs.clear();
  _BlockAABBs.reserve(blocknum);
  DivideBlock3D(_BaseBlock, blocknum, _BlockAABBs);
}

template <typename T>
void BlockGeometry3D<T>::CreateBlocks(int blocknum) {
  DivideBlocks(blocknum);

  _Blocks.clear();
  _Blocks.reserve(_BlockAABBs.size());
  _BasicBlocks.clear();
  _BasicBlocks.reserve(_BlockAABBs.size());
  // create blocks
  int blockid = 0;
  for (const AABB<int, 3> &blockaabb : _BlockAABBs) {
    Vector<T, 3> MIN =
      blockaabb.getMin() * _BaseBlock.getVoxelSize() + BasicBlock<T, 3>::_min;
    Vector<T, 3> MAX =
      (blockaabb.getMax() + Vector<T, 3>{T(1)}) * _BaseBlock.getVoxelSize() +
      BasicBlock<T, 3>::_min;
    AABB<T, 3> aabb(MIN, MAX);
    _Blocks.emplace_back(aabb, blockaabb, blockid, _BaseBlock.getVoxelSize(), _overlap);
    _BasicBlocks.emplace_back(_BaseBlock.getVoxelSize(), aabb, blockaabb, blockid);
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
        if (isOverlapped(block, blockn.getBaseBlock())) {
          nbrsvec.push_back(&blockn);
        }
      }
    }
  }
}

// comms

template <typename T>
void BlockGeometry3D<T>::InitComm() {
  for (Block3D<T> &block : _Blocks) {
    std::vector<BlockComm<T, 3>> &Communicators = block.getCommunicators();
    Communicators.clear();
    // get the first layer of overlapped cells(counted from inside to outside)
    BasicBlock<T, 3> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
    std::uint8_t blocklevel = block.getLevel();
    for (Block3D<T> *nblock : block.getNeighbors()) {
      // check if 2 blocks are of the same level
      if (nblock->getLevel() == blocklevel) {
        Communicators.emplace_back(nblock);
        BlockComm<T, 3> &comm = Communicators.back();
        // blocks of the same level only communicate with the first layer of overlapped
        // cells
        block.getCellIdx(baseblock_ext1, nblock->getBaseBlock(), comm.RecvCells);
        nblock->getCellIdx(nblock->getBaseBlock(), baseblock_ext1, comm.SendCells);
      }
    }
  }
}

template <typename T>
void BlockGeometry3D<T>::InitAverComm() {
  for (Block3D<T> &block : _Blocks) {
    std::vector<IntpBlockComm<T, 3>> &Communicators = block.getAverageBlockComm();
    Communicators.clear();
    // get the first layer of overlapped cells(counted from inside to outside)
    BasicBlock<T, 3> baseblock_ext1 = block.getBaseBlock().getExtBlock(1);
    std::uint8_t blocklevel = block.getLevel();
    std::size_t XY = block.getNx() * block.getNy();
    for (Block3D<T> *nblock : block.getNeighbors()) {
      // find block of blocklevel+1
      if (nblock->getLevel() == blocklevel + 1) {
        Communicators.emplace_back(nblock);
        IntpBlockComm<T, 3> &comm = Communicators.back();
        // vox size
        const T Cvoxsize = block.getVoxelSize();
        const T Fvoxsize = nblock->getVoxelSize();
        // get intersection
        const AABB<T, 3> intsec = getIntersection(baseblock_ext1, nblock->getBaseBlock());
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        int CNz = static_cast<int>(std::round(intsec.getExtension()[2] / Cvoxsize));
        // start index of intsec in block
        Vector<T, 3> startC = intsec.getMin() - block.getMin();
        int startCx = static_cast<int>(std::round(startC[0] / Cvoxsize));
        int startCy = static_cast<int>(std::round(startC[1] / Cvoxsize));
        int startCz = static_cast<int>(std::round(startC[2] / Cvoxsize));
        // start index of intsec in nblock
        Vector<T, 3> startF = intsec.getMin() - nblock->getMin();
        int startFx = static_cast<int>(std::round(startF[0] / Fvoxsize));
        int startFy = static_cast<int>(std::round(startF[1] / Fvoxsize));
        int startFz = static_cast<int>(std::round(startF[2] / Fvoxsize));

        std::size_t nXY = nblock->getNx() * nblock->getNy();
        for (int iz = 0; iz < CNz; ++iz) {
          for (int iy = 0; iy < CNy; ++iy) {
            for (int ix = 0; ix < CNx; ++ix) {
              std::size_t Cid =
                (iz + startCz) * XY + (iy + startCy) * block.getNx() + ix + startCx;
              std::size_t Fid0 = (iz * 2 + startFz) * nXY +
                                 (iy * 2 + startFy) * nblock->getNx() + ix * 2 + startFx;
              std::size_t Fid1 = Fid0 + 1;
              std::size_t Fid2 = Fid0 + nblock->getNx();
              std::size_t Fid3 = Fid2 + 1;

              std::size_t Fid4 = Fid0 + nXY;
              std::size_t Fid5 = Fid4 + 1;
              std::size_t Fid6 = Fid4 + nblock->getNx();
              std::size_t Fid7 = Fid6 + 1;

              comm.RecvCells.push_back(Cid);
              comm.SendCells.emplace_back(
                IntpSource<3>{Fid0, Fid1, Fid2, Fid3, Fid4, Fid5, Fid6, Fid7});
            }
          }
        }
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
    std::uint8_t blocklevel = block.getLevel();
    std::vector<IntpBlockComm<T, 3>> &Communicators = block.getIntpBlockComm();
    Communicators.clear();
    std::size_t XY = block.getNx() * block.getNy();
    for (Block3D<T> *nblock : block.getNeighbors()) {
      // find block of blocklevel-1
      if (nblock->getLevel() == blocklevel - 1) {
        Communicators.emplace_back(nblock);
        IntpBlockComm<T, 3> &comm = Communicators.back();
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
        // int startFz = static_cast<int>(std::round(startF[2] / Fvoxsize));

        std::size_t nXY = nblock->getNx() * nblock->getNy();
        for (int iz = 0; iz < CNz; ++iz) {
          for (int iy = 0; iy < CNy; ++iy) {
            for (int ix = 0; ix < CNx; ++ix) {
              // original
              std::size_t Cid0 =
                (iz + startCz) * nXY + (iy + startCy) * nblock->getNx() + ix + startCx;
              std::size_t Cid1 = Cid0 + 1;
              std::size_t Cid2 = Cid0 + nblock->getNx();
              std::size_t Cid3 = Cid2 + 1;
              std::size_t Cid4 = Cid0 + nXY;
              std::size_t Cid5 = Cid4 + 1;
              std::size_t Cid6 = Cid4 + nblock->getNx();
              std::size_t Cid7 = Cid6 + 1;
              std::size_t Fid = (iy * 2 + startFy) * block.getNx() + ix * 2 + startFx;

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
              comm.RecvCells.push_back(Fid);
              comm.SendCells.emplace_back(
                IntpSource<3>{Cid0, Cid1, Cid2, Cid3, Cid4, Cid5, Cid6, Cid7});

              // 1
              comm.RecvCells.push_back(Fid_x);
              comm.SendCells.emplace_back(IntpSource<3>{
                Cid0_x, Cid1_x, Cid2_x, Cid3_x, Cid4_x, Cid5_x, Cid6_x, Cid7_x});

              // 2
              comm.RecvCells.push_back(Fid_y);
              comm.SendCells.emplace_back(IntpSource<3>{
                Cid0_y, Cid1_y, Cid2_y, Cid3_y, Cid4_y, Cid5_y, Cid6_y, Cid7_y});

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

              comm.RecvCells.push_back(Fid_xy);
              comm.SendCells.emplace_back(IntpSource<3>{
                Cid0_xy, Cid1_xy, Cid2_xy, Cid3_xy, Cid4_xy, Cid5_xy, Cid6_xy, Cid7_xy});

              // 4
              comm.RecvCells.push_back(Fid_z);
              comm.SendCells.emplace_back(IntpSource<3>{
                Cid0_z, Cid1_z, Cid2_z, Cid3_z, Cid4_z, Cid5_z, Cid6_z, Cid7_z});

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

              comm.RecvCells.push_back(Fid_xz);
              comm.SendCells.emplace_back(IntpSource<3>{
                Cid0_xz, Cid1_xz, Cid2_xz, Cid3_xz, Cid4_xz, Cid5_xz, Cid6_xz, Cid7_xz});

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

              comm.RecvCells.push_back(Fid_yz);
              comm.SendCells.emplace_back(IntpSource<3>{
                Cid0_yz, Cid1_yz, Cid2_yz, Cid3_yz, Cid4_yz, Cid5_yz, Cid6_yz, Cid7_yz});

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

              comm.RecvCells.push_back(Fid_xyz);
              comm.SendCells.emplace_back(IntpSource<3>{Cid0_xyz, Cid1_xyz, Cid2_xyz,
                                                          Cid3_xyz, Cid4_xyz, Cid5_xyz,
                                                          Cid6_xyz, Cid7_xyz});
            }
          }
        }
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
void BlockGeometry3D<T>::InitMPIComm(BlockGeometryHelper3D<T> &GeoHelper) {
  for (Block3D<T> &block : _Blocks) {
    MPIBlockComm &MPIComm = block.getMPIBlockComm();
    MPIComm.clear();
    std::uint8_t blocklevel = block.getLevel();
    // base block of block
    const BasicBlock<T, 3> &baseblock = block.getBaseBlock();
    // get the first layer of overlapped cells(counted from inside to outside)
    // int overlap = (blocklevel != std::uint8_t(0)) ? 2 : 1;
    BasicBlock<T, 3> baseblock_ext1 = baseblock.getExtBlock(1);
    // find neighbors
    std::vector<std::pair<int, int>> &nbrs =
      GeoHelper.getMPIBlockNbrs(block.getBlockId());
    for (const std::pair<int, int> &nbr : nbrs) {
      // check if 2 blocks are of the same level
      const BasicBlock<T, 3> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel) {
        BasicBlock<T, 3> nbaseblock_ext1 = nbaseblock.getExtBlock(1);
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
void BlockGeometry3D<T>::InitMPIAverComm(BlockGeometryHelper3D<T> &GeoHelper) {
  for (Block3D<T> &block : _Blocks) {
    MPIIntpBlockComm<T, 3> &MPIComm = block.getMPIAverBlockComm();
    MPIComm.clear();
    std::uint8_t blocklevel = block.getLevel();
    const BasicBlock<T, 3> &baseblock = block.getBaseBlock();
    // first layer of overlapped cells(counted from inside to outside)
    BasicBlock<T, 3> baseblock_ext1 = baseblock.getExtBlock(1);
    std::vector<std::pair<int, int>> &nbrs =
      GeoHelper.getMPIBlockNbrs(block.getBlockId());
    std::size_t XY = block.getNx() * block.getNy();
    for (const std::pair<int, int> &nbr : nbrs) {
      const BasicBlock<T, 3> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel - 1) {
        BasicBlock<T, 3> nbaseblock_ext1 = nbaseblock.getExtBlock(1);
        // init sender
        MPIComm.Senders.emplace_back(nbr.first, nbaseblock.getBlockId());
        MPIIntpBlockSendStru<T, 3> &sender = MPIComm.Senders.back();
        // vox size
        const T Cvoxsize = nbaseblock.getVoxelSize();
        const T Fvoxsize = block.getVoxelSize();
        // init sender
        // get intersection
        const AABB<T, 3> intsec = getIntersection(baseblock, nbaseblock_ext1);
        // use coarse grid size here, convient for calculating fine cell index
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        int CNz = static_cast<int>(std::round(intsec.getExtension()[2] / Cvoxsize));
        // start index of intsec in block
        Vector<T, 3> startF = intsec.getMin() - block.getMin();
        int startFx = static_cast<int>(std::round(startF[0] / Fvoxsize));
        int startFy = static_cast<int>(std::round(startF[1] / Fvoxsize));
        int startFz = static_cast<int>(std::round(startF[2] / Fvoxsize));

        for (int iz = 0; iz < CNz; ++iz) {
          for (int iy = 0; iy < CNy; ++iy) {
            for (int ix = 0; ix < CNx; ++ix) {
              std::size_t Fid0 = (iz * 2 + startFz) * XY +
                                 (iy * 2 + startFy) * block.getNx() + ix * 2 + startFx;
              std::size_t Fid1 = Fid0 + 1;
              std::size_t Fid2 = Fid0 + block.getNx();
              std::size_t Fid3 = Fid2 + 1;

              std::size_t Fid4 = Fid0 + XY;
              std::size_t Fid5 = Fid4 + 1;
              std::size_t Fid6 = Fid4 + block.getNx();
              std::size_t Fid7 = Fid6 + 1;

              sender.SendCells.emplace_back(
                IntpSource<3>{Fid0, Fid1, Fid2, Fid3, Fid4, Fid5, Fid6, Fid7});
            }
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
        const AABB<T, 3> intsec = getIntersection(baseblock_ext1, nbaseblock);
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        int CNz = static_cast<int>(std::round(intsec.getExtension()[2] / Cvoxsize));
        // start index of intsec in block
        Vector<T, 3> startC = intsec.getMin() - block.getMin();
        int startCx = static_cast<int>(std::round(startC[0] / Cvoxsize));
        int startCy = static_cast<int>(std::round(startC[1] / Cvoxsize));
        int startCz = static_cast<int>(std::round(startC[2] / Cvoxsize));

        for (int iz = 0; iz < CNz; ++iz) {
          for (int iy = 0; iy < CNy; ++iy) {
            for (int ix = 0; ix < CNx; ++ix) {
              std::size_t Cid =
                (iz + startCz) * XY + (iy + startCy) * block.getNx() + ix + startCx;
              recver.RecvCells.push_back(Cid);
            }
          }
        }
        block._NeedMPIComm = true;
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
    MPIIntpBlockComm<T, 3> &MPIComm = block.getMPIIntpBlockComm();
    MPIComm.clear();
    std::uint8_t blocklevel = block.getLevel();
    const BasicBlock<T, 3> &baseblock = block.getBaseBlock();
    // 2 layers of overlapped cells(counted from inside to outside)
    // BasicBlock<T, 3> baseblock_ext2 = baseblock.getExtBlock(2);
    std::vector<std::pair<int, int>> &nbrs =
      GeoHelper.getMPIBlockNbrs(block.getBlockId());
    std::size_t XY = block.getNx() * block.getNy();
    for (const std::pair<int, int> &nbr : nbrs) {
      const BasicBlock<T, 3> &nbaseblock = GeoHelper.getAllBasicBlock(static_cast<std::size_t>(nbr.second));
      if (nbaseblock.getLevel() == blocklevel + 1) {
        BasicBlock<T, 3> nbaseblock_ext2 = nbaseblock.getExtBlock(2);
        // init sender
        MPIComm.Senders.emplace_back(nbr.first, nbaseblock.getBlockId());
        MPIIntpBlockSendStru<T, 3> &sender = MPIComm.Senders.back();
        // vox size
        const T Cvoxsize = block.getVoxelSize();
        // const T Fvoxsize = nbaseblock.getVoxelSize();
        // init sender
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

        for (int iz = 0; iz < CNz; ++iz) {
          for (int iy = 0; iy < CNy; ++iy) {
            for (int ix = 0; ix < CNx; ++ix) {
              // original
              std::size_t Cid0 =
                (iz + startCz) * XY + (iy + startCy) * block.getNx() + ix + startCx;
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
              sender.SendCells.emplace_back(
                IntpSource<3>{Cid0, Cid1, Cid2, Cid3, Cid4, Cid5, Cid6, Cid7});

              // 1
              sender.SendCells.emplace_back(IntpSource<3>{
                Cid0_x, Cid1_x, Cid2_x, Cid3_x, Cid4_x, Cid5_x, Cid6_x, Cid7_x});

              // 2
              sender.SendCells.emplace_back(IntpSource<3>{
                Cid0_y, Cid1_y, Cid2_y, Cid3_y, Cid4_y, Cid5_y, Cid6_y, Cid7_y});

              // 3
              std::size_t Cid0_xy = Cid0_y + 1;
              std::size_t Cid1_xy = Cid1_y + 1;
              std::size_t Cid2_xy = Cid2_y + 1;
              std::size_t Cid3_xy = Cid3_y + 1;
              std::size_t Cid4_xy = Cid4_y + 1;
              std::size_t Cid5_xy = Cid5_y + 1;
              std::size_t Cid6_xy = Cid6_y + 1;
              std::size_t Cid7_xy = Cid7_y + 1;

              sender.SendCells.emplace_back(IntpSource<3>{
                Cid0_xy, Cid1_xy, Cid2_xy, Cid3_xy, Cid4_xy, Cid5_xy, Cid6_xy, Cid7_xy});

              // 4
              sender.SendCells.emplace_back(IntpSource<3>{
                Cid0_z, Cid1_z, Cid2_z, Cid3_z, Cid4_z, Cid5_z, Cid6_z, Cid7_z});

              // 5
              std::size_t Cid0_xz = Cid0_z + 1;
              std::size_t Cid1_xz = Cid1_z + 1;
              std::size_t Cid2_xz = Cid2_z + 1;
              std::size_t Cid3_xz = Cid3_z + 1;
              std::size_t Cid4_xz = Cid4_z + 1;
              std::size_t Cid5_xz = Cid5_z + 1;
              std::size_t Cid6_xz = Cid6_z + 1;
              std::size_t Cid7_xz = Cid7_z + 1;

              sender.SendCells.emplace_back(IntpSource<3>{
                Cid0_xz, Cid1_xz, Cid2_xz, Cid3_xz, Cid4_xz, Cid5_xz, Cid6_xz, Cid7_xz});

              // 6
              std::size_t Cid0_yz = Cid0_z + block.getNx();
              std::size_t Cid1_yz = Cid1_z + block.getNx();
              std::size_t Cid2_yz = Cid2_z + block.getNx();
              std::size_t Cid3_yz = Cid3_z + block.getNx();
              std::size_t Cid4_yz = Cid4_z + block.getNx();
              std::size_t Cid5_yz = Cid5_z + block.getNx();
              std::size_t Cid6_yz = Cid6_z + block.getNx();
              std::size_t Cid7_yz = Cid7_z + block.getNx();

              sender.SendCells.emplace_back(IntpSource<3>{
                Cid0_yz, Cid1_yz, Cid2_yz, Cid3_yz, Cid4_yz, Cid5_yz, Cid6_yz, Cid7_yz});

              // 7
              std::size_t Cid0_xyz = Cid0_yz + 1;
              std::size_t Cid1_xyz = Cid1_yz + 1;
              std::size_t Cid2_xyz = Cid2_yz + 1;
              std::size_t Cid3_xyz = Cid3_yz + 1;
              std::size_t Cid4_xyz = Cid4_yz + 1;
              std::size_t Cid5_xyz = Cid5_yz + 1;
              std::size_t Cid6_xyz = Cid6_yz + 1;
              std::size_t Cid7_xyz = Cid7_yz + 1;

              sender.SendCells.emplace_back(IntpSource<3>{Cid0_xyz, Cid1_xyz, Cid2_xyz,
                                                            Cid3_xyz, Cid4_xyz, Cid5_xyz,
                                                            Cid6_xyz, Cid7_xyz});
            }
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
        const AABB<T, 3> intsec = getIntersection(block, nbaseblock);
        // use coarse grid size here, convient for calculating fine cell index
        int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
        int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
        int CNz = static_cast<int>(std::round(intsec.getExtension()[2] / Cvoxsize));
        // start index of intsec in block
        Vector<T, 3> startF = intsec.getMin() - block.getMin();
        int startFx = static_cast<int>(std::round(startF[0] / Fvoxsize));
        int startFy = static_cast<int>(std::round(startF[1] / Fvoxsize));
        int startFz = static_cast<int>(std::round(startF[2] / Fvoxsize));

        for (int iz = 0; iz < CNz; ++iz) {
          for (int iy = 0; iy < CNy; ++iy) {
            for (int ix = 0; ix < CNx; ++ix) {
              std::size_t Fid0 = (iz * 2 + startFz) * XY +
                                 (iy * 2 + startFy) * block.getNx() + ix * 2 + startFx;
              std::size_t Fid1 = Fid0 + 1;
              std::size_t Fid2 = Fid0 + block.getNx();
              std::size_t Fid3 = Fid2 + 1;
              std::size_t Fid4 = Fid0 + XY;
              std::size_t Fid5 = Fid4 + 1;
              std::size_t Fid6 = Fid4 + block.getNx();
              std::size_t Fid7 = Fid6 + 1;
              recver.RecvCells.push_back(Fid0);
              recver.RecvCells.push_back(Fid1);
              recver.RecvCells.push_back(Fid2);
              recver.RecvCells.push_back(Fid3);
              recver.RecvCells.push_back(Fid4);
              recver.RecvCells.push_back(Fid5);
              recver.RecvCells.push_back(Fid6);
              recver.RecvCells.push_back(Fid7);
            }
          }
        }
        block._NeedMPIComm = true;
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


// BlockGeometryHelper3D
#include "lbm/lattice_set.h"

template <typename T>
BlockGeometryHelper3D<T>::BlockGeometryHelper3D(int Nx, int Ny, int Nz,
                                                const AABB<T, 3> &AABBs, T voxelSize,
                                                int blockcelllen, std::uint8_t llimit,
                                                int ext)
    : BasicBlock<T, 3>(
        voxelSize, AABBs.getExtended(Vector<T, 3>{voxelSize * ext}),
        AABB<int, 3>(Vector<int, 3>{0}, Vector<int, 3>{Nx + 1, Ny + 1, Nz + 1})),
      _BaseBlock(voxelSize, AABBs,
                 AABB<int, 3>(Vector<int, 3>{1}, Vector<int, 3>{Nx, Ny, Nz})),
      BlockCellLen(blockcelllen), Ext(ext), _LevelLimit(llimit), _MaxLevel(std::uint8_t(0)), 
      _Exchanged(true), _IndexExchanged(true) {
  if (BlockCellLen < 4) {
    std::cerr << "BlockGeometryHelper3D<T>, BlockCellLen < 4" << std::endl;
  }
  CellsNx = _BaseBlock.getNx() / BlockCellLen;
  CellsNy = _BaseBlock.getNy() / BlockCellLen;
  CellsNz = _BaseBlock.getNz() / BlockCellLen;
  CellsN = CellsNx * CellsNy * CellsNz;

  Vector<int, 3> Projection{1, CellsNx, CellsNx * CellsNy};

  Delta_Cellidx = make_Array<int, D3Q27<T>::q - 1>(
    [&](int i) { return D3Q27<T>::c[i + 1] * Projection; });

  CreateBlockCells();
}

template <typename T>
void BlockGeometryHelper3D<T>::CreateBlockCells() {
  // buffer vector to store cell aabbs
  std::vector<AABB<int, 3>> AABBCells;
  AABBCells.reserve(CellsN);
  // divide base block into cell aabbs
  _BaseBlock.getIdxBlock().divide(CellsNx, CellsNy, CellsNz, AABBCells);
  // create cell blocks from cell aabbs
  _BlockCells.reserve(CellsN);
  int blockid = 0;
  T voxsize = BasicBlock<T, 3>::getVoxelSize();
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
void BlockGeometryHelper3D<T>::UpdateMaxLevel() {
  _MaxLevel = std::uint8_t(0);
  for (const BasicBlock<T, 3> &block : _BlockCells) {
    if (block.getLevel() > _MaxLevel) _MaxLevel = block.getLevel();
  }
}

template <typename T>
void BlockGeometryHelper3D<T>::CreateBlocks() {
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
            }
          }
          NewMesh[2] += _BlockCells[id + Nz * XY].getNz();
          ++Nz;
        }
      end_z_expansion:

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
}

template <typename T>
void BlockGeometryHelper3D<T>::CreateBlocks(int blocknum) {
  std::vector<AABB<int, 3>> BlockAABBs;
  BlockAABBs.reserve(blocknum);
  DivideBlock3D(_BaseBlock, blocknum, BlockAABBs);

  // create new blocks on relatively older blocks
  std::vector<BasicBlock<T, 3>> &BasicBlocks = getAllOldBasicBlocks();
  // now old blocks become new blocks
  _Exchanged = !_Exchanged;

  BasicBlocks.clear();
  BasicBlocks.reserve(BlockAABBs.size());
  // create blocks
  int blockid = 0;
  for (const AABB<int, 3> &blockaabb : BlockAABBs) {
    Vector<T, 3> MIN =
      blockaabb.getMin() * _BaseBlock.getVoxelSize() + BasicBlock<T, 3>::_min;
    Vector<T, 3> MAX =
      (blockaabb.getMax() + Vector<T, 3>{T(1)}) * _BaseBlock.getVoxelSize() +
      BasicBlock<T, 3>::_min;
    AABB<T, 3> aabb(MIN, MAX);
    BasicBlocks.emplace_back(_BaseBlock.getVoxelSize(), aabb, blockaabb, blockid);
    blockid++;
  }
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
void BlockGeometryHelper3D<T>::Optimize(int ProcessNum, bool enforce) {
  Optimize(getAllBasicBlocks(), ProcessNum, enforce);
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

  if (enforce && Blocks.size() < static_cast<std::size_t>(ProcessNum)) {
    // sort blocks by number of points, large blocks first
    std::sort(Blocks.begin(), Blocks.end(),
              [](const BasicBlock<T, 3> &a, const BasicBlock<T, 3> &b) {
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
      DivideBlock3D(*it, part, NewBasicBlocks);
      // remove block
      it = Blocks.erase(it);
      count += part - 1;
    }
  } else {
    while (it != Blocks.end()) {
      T ratio = std::round(static_cast<T>(it->getN()) / NumPerProcess);
      if (ratio >= 2) {
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
            const BasicBlock<T, 3> nblock = getAllBasicBlock(static_cast<std::size_t>(nblockid)).getExtBlock(1);
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