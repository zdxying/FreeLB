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

// field_struct.hh

#pragma once

#include "data_struct/field_struct.h"


// BlockField

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockField<FieldType, FloatType, Dim>::BlockField(Block<FloatType, Dim>& block)
    : FieldType(block.getN()), _Block(block) {
  constructInDevice();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockField<FieldType, FloatType, Dim>::BlockField(
  Block<FloatType, Dim>& block, datatype initvalue)
    : FieldType(block.getN(), initvalue), _Block(block) {
  constructInDevice();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockField<FieldType, FloatType, Dim>::BlockField(const BlockField& blockF)
    : FieldType(blockF), _Block(blockF._Block) {
  constructInDevice();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockField<FieldType, FloatType, Dim>::BlockField(BlockField&& blockF) noexcept
    : FieldType(std::move(blockF)), _Block(blockF._Block) {
  constructInDevice();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockField<FieldType, FloatType, Dim>::~BlockField() {
#ifdef __CUDACC__
  if (dev_BlockField) cuda_free(dev_BlockField);
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockField<FieldType, FloatType, Dim>::constructInDevice() {
#ifdef __CUDACC__
  dev_BlockField = cuda_malloc<cudev_BlockFieldType>(1);
  // temp host object
  cudev_BlockFieldType temp(this->get_devptr());
  // copy to device
  host_to_device(dev_BlockField, &temp, 1);
#endif
}

#ifdef MPI_ENABLED

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockField<FieldType, FloatType, Dim>::mpiNormalSend(
  std::vector<MPI_Request>& SendRequests, std::vector<std::vector<datatype>>& SendBuffers,
  const std::vector<DistributedComm>& MPISends) {
  // add to send buffer
  for (std::size_t i = 0; i < MPISends.size(); ++i) {
    const DistributedComm& comm = MPISends[i];
    std::vector<datatype>& buffer = SendBuffers[i];
    buffer.resize(comm.Cells.size() * array_dim);
    const std::vector<std::size_t>& sends = comm.Cells;
    std::size_t bufidx{};
    for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
      const auto& Array = FieldType::getField(iArr);
      for (std::size_t id : sends) {
        buffer[bufidx] = Array[id];
        ++bufidx;
      }
    }
    // non-blocking send
    MPI_Request request;
    mpi().iSend(
      buffer.data(), buffer.size(), comm.TargetRank, &request, comm.TargetBlockId);
    SendRequests.push_back(request);
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockField<FieldType, FloatType, Dim>::mpiAverSend(
  std::vector<MPI_Request>& SendRequests, std::vector<std::vector<datatype>>& SendBuffers,
  const std::vector<DistributedComm>& MPISends) {
  static constexpr std::size_t SendSize = Dim == 2 ? 4 : 8;
  // add to send buffer
  for (std::size_t i = 0; i < MPISends.size(); ++i) {
    const DistributedComm& comm = MPISends[i];
    std::vector<datatype>& buffer = SendBuffers[i];
    buffer.resize(comm.Cells.size() / SendSize * array_dim);
    const std::vector<std::size_t>& sends = comm.Cells;
    std::size_t bufidx{};
    for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
      const auto& Array = FieldType::getField(iArr);
      for (std::size_t id = 0; id < sends.size(); id += SendSize) {
        buffer[bufidx] = getAverage<FloatType, Dim>(Array, id, sends);
        ++bufidx;
      }
    }
    // non-blocking send
    MPI_Request request;
    mpi().iSend(
      buffer.data(), buffer.size(), comm.TargetRank, &request, comm.TargetBlockId);
    SendRequests.push_back(request);
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockField<FieldType, FloatType, Dim>::mpiIntpSend(
  std::vector<MPI_Request>& SendRequests, std::vector<std::vector<datatype>>& SendBuffers,
  const std::vector<DistributedComm>& MPISends) {
  static constexpr std::size_t SendSize = Dim == 2 ? 4 : 8;
  // add to send buffer
  for (int i = 0; i < MPISends.size(); ++i) {
    const DistributedComm& comm = MPISends[i];
    std::vector<datatype>& buffer = SendBuffers[i];
    buffer.resize(comm.Cells.size() / SendSize * array_dim);
    const std::vector<std::size_t>& sends = comm.Cells;
    std::size_t bufidx{};
    for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
      const auto& Array = FieldType::getField(iArr);
      for (std::size_t id = 0; id < sends.size();) {
        Interpolation<FloatType, Dim>(buffer, bufidx, Array, sends, id);
      }
    }
    // non-blocking send
    MPI_Request request;
    mpi().iSend(
      buffer.data(), buffer.size(), comm.TargetRank, &request, comm.TargetBlockId);
    SendRequests.push_back(request);
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockField<FieldType, FloatType, Dim>::mpiRecv(
  std::vector<MPI_Request>& RecvRequests, std::vector<std::vector<datatype>>& RecvBuffers,
  const std::vector<DistributedComm>& MPIRecvs) {
  // non-blocking recv
  for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
    const DistributedComm& comm = MPIRecvs[i];
    std::vector<datatype>& buffer = RecvBuffers[i];
    buffer.resize(comm.Cells.size() * array_dim);
    MPI_Request request;
    mpi().iRecv(
      buffer.data(), buffer.size(), comm.TargetRank, &request, _Block.getBlockId());
    RecvRequests.push_back(request);
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockField<FieldType, FloatType, Dim>::mpiSet(int& reqidx,
  std::vector<MPI_Request>& RecvRequests,
  const std::vector<std::vector<datatype>>& RecvBuffers,
  const std::vector<DistributedComm>& MPIRecvs) {
  // wait and set field data
  for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
    MPI_Wait(&RecvRequests[i + reqidx], MPI_STATUS_IGNORE);
    const DistributedComm& comm = MPIRecvs[i];
    const std::vector<datatype>& buffer = RecvBuffers[i];
    const std::vector<std::size_t>& recvs = comm.Cells;
    std::size_t bufidx{};
    for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
      auto& Array = FieldType::getField(iArr);
      for (std::size_t id : recvs) {
        Array[id] = buffer[bufidx];
        ++bufidx;
      }
    }
  }
  reqidx += MPIRecvs.size();
}

#endif


// BlockFieldManager

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockFieldManager<FieldType, FloatType, Dim>::BlockFieldManager(
  BlockGeometry<FloatType, Dim>& blockgeometry)
    : _BlockGeo(blockgeometry) {
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    _Fields.emplace_back(block);
  }
  InitDeviceData();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockFieldManager<FieldType, FloatType, Dim>::BlockFieldManager(
  BlockGeometry<FloatType, Dim>& blockgeometry, datatype initvalue)
    : _BlockGeo(blockgeometry) {
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    _Fields.emplace_back(block, initvalue);
  }
  InitDeviceData();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockFieldManager<FieldType, FloatType, Dim>::BlockFieldManager(
  const BlockFieldManager& blockFManager)
    : _Fields(blockFManager._Fields), _BlockGeo(blockFManager._BlockGeo) {
#ifdef __CUDACC__
  std::size_t FieldNum = _Fields.size();
  dev_Fields = cuda_malloc<cudev_BlockFieldType*>(FieldNum);
  device_to_device(dev_Fields, blockFManager.dev_Fields, FieldNum);
  constructInDevice();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockFieldManager<FieldType, FloatType, Dim>::BlockFieldManager(
  BlockFieldManager&& blockFManager) noexcept
    : _Fields(std::move(blockFManager._Fields)), _BlockGeo(blockFManager._BlockGeo) {
#ifdef __CUDACC__
  dev_Fields = blockFManager.dev_Fields;
  constructInDevice();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::InitDeviceData() {
#ifdef __CUDACC__
  std::size_t FieldNum = _Fields.size();
  dev_Fields = cuda_malloc<cudev_BlockFieldType*>(FieldNum);
  cudev_BlockFieldType* host_Data[FieldNum];
  for (std::size_t i = 0; i < FieldNum; ++i) host_Data[i] = _Fields[i].get_devObj();
  host_to_device(dev_Fields, host_Data, FieldNum);
  constructInDevice();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockField<FieldType, FloatType, Dim>*
BlockFieldManager<FieldType, FloatType, Dim>::findBlockField(int blockid) {
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getBlockId() == blockid) return &blockF;
  }
  std::cerr << "[BlockFieldManager]: can't find blockF with blockid " << blockid
            << std::endl;
  exit(1);
}


#ifdef __CUDACC__

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::constructInDevice() {
  dev_BlockFieldManager = cuda_malloc<cudev_BlockFieldManagerType>(1);
  // temp host object
  cudev_BlockFieldManagerType temp(dev_Fields);
  // copy to device
  host_to_device(dev_BlockFieldManager, &temp, 1);
}

#endif


template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::InitAndComm() {
  Init();
  AllCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::InitAndComm(datatype initvalue) {
  Init(initvalue);
  AllCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::InitAndComm(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper) {
  Init(GeoHelper);
  AllCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::InitAndComm(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue) {
  Init(GeoHelper, initvalue);
  AllCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::Init() {
  _Fields.clear();
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    _Fields.emplace_back(block);
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::Init(datatype initvalue) {
  _Fields.clear();
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    _Fields.emplace_back(block, initvalue);
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::NonFieldInit(datatype initvalue) {
  _Fields.clear();
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    _Fields.emplace_back(block, initvalue);
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::Init(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper) {
  std::vector<BlockField<FieldType, FloatType, Dim>> NewFields;
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    NewFields.emplace_back(block);
  }
  // data transfer
  FieldDataTransfer(GeoHelper, NewFields);
  _Fields.swap(NewFields);
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::Init(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue) {
  std::vector<BlockField<FieldType, FloatType, Dim>> NewFields;
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    NewFields.emplace_back(block, initvalue);
  }
  // data transfer
  FieldDataTransfer(GeoHelper, NewFields);
  _Fields.swap(NewFields);
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename FlagFieldType, typename Func>
void BlockFieldManager<FieldType, FloatType, Dim>::InitCopy(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue,
  const BlockFieldManager<FlagFieldType, FloatType, Dim>& FlagFManager, std::uint8_t flag,
  const Func& func) {
  std::vector<BlockField<FieldType, FloatType, Dim>> NewFields;
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    NewFields.emplace_back(block, initvalue);
  }
  // func
  // forEach(FlagFManager, flag, func); but this is operate on NewFields
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    auto& field = NewFields[iblock];
    const auto& flagarr = FlagFManager.getBlockField(iblock).getField(0);
    blockgeo.forEach(flagarr, flag, [&field, &func](std::size_t id) { func(field, id); });
    ++iblock;
  }
  // data transfer
  FieldDataCopyTransfer(GeoHelper, NewFields);
  _Fields.swap(NewFields);
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::FieldDataTransfer(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper,
  std::vector<BlockField<FieldType, FloatType, Dim>>& NewFields) {
  if constexpr (!FieldType::isField) {
    return;
  }
  // copy from old field to new field
#pragma omp parallel for num_threads(Thread_Num)
  for (std::size_t inewblock = 0; inewblock < NewFields.size(); ++inewblock) {
    BlockField<FieldType, FloatType, Dim>& NewField = NewFields[inewblock];
    // newbaseblock could be accessed either from GeoHelper or _BlockGeo
    const BasicBlock<FloatType, Dim>& newbaseblock =
      GeoHelper.getAllBasicBlock(inewblock);
    const BasicBlock<FloatType, Dim>& newblock = NewField.getBlock();
    std::uint8_t Level = newblock.getLevel();
    // find overlapped old block field
    for (std::size_t iblock = 0; iblock < GeoHelper.getAllOldBasicBlocks().size();
         ++iblock) {
      // oldbaseblock could be accessed only from GeoHelper
      const BasicBlock<FloatType, Dim>& baseblock = GeoHelper.getAllOldBasicBlock(iblock);
      if (isOverlapped(newbaseblock, baseblock)) {
        const BlockField<FieldType, FloatType, Dim>& Field = _Fields[iblock];
        int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
        const BasicBlock<FloatType, Dim>& block = baseblock.getExtBlock(overlap);
        if (Level == block.getLevel()) {
          // copy
          FieldCopy2D(Field, NewField, block, baseblock, newblock, newbaseblock);
        } else if (Level > block.getLevel()) {
          // interp
          FieldInterpolation2D(Field, NewField, block, baseblock, newblock, newbaseblock);
        } else if (Level < block.getLevel()) {
          // average
          FieldAverage2D(Field, NewField, block, baseblock, newblock, newbaseblock);
        }
      }
    }
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::FieldDataCopyTransfer(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper,
  std::vector<BlockField<FieldType, FloatType, Dim>>& NewFields) {
  // copy from old field to new field
#pragma omp parallel for num_threads(Thread_Num)
  for (std::size_t inewblock = 0; inewblock < NewFields.size(); ++inewblock) {
    BlockField<FieldType, FloatType, Dim>& NewField = NewFields[inewblock];
    // newbaseblock could be accessed either from GeoHelper or _BlockGeo
    const BasicBlock<FloatType, Dim>& newbaseblock =
      GeoHelper.getAllBasicBlock(inewblock);
    const BasicBlock<FloatType, Dim>& newblock = NewField.getBlock();
    std::uint8_t Level = newblock.getLevel();
    // find overlapped old block field
    for (std::size_t iblock = 0; iblock < GeoHelper.getAllOldBasicBlocks().size();
         ++iblock) {
      // oldbaseblock could be accessed only from GeoHelper
      const BasicBlock<FloatType, Dim>& baseblock = GeoHelper.getAllOldBasicBlock(iblock);
      if (isOverlapped(newbaseblock, baseblock)) {
        const BlockField<FieldType, FloatType, Dim>& Field = _Fields[iblock];
        int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
        const BasicBlock<FloatType, Dim>& block = baseblock.getExtBlock(overlap);
        if (Level == block.getLevel()) {
          // copy
          FieldCopy2D(Field, NewField, block, baseblock, newblock, newbaseblock);
        }
      }
    }
  }
}


template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename LatSet>
void BlockFieldManager<FieldType, FloatType, Dim>::SetupBoundary(
  const AABB<FloatType, Dim>& block, datatype bdvalue) {
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    auto& field = _Fields[iblock];
    blockgeo.template SetupBoundary<FieldType, LatSet>(block, field, bdvalue);
    ++iblock;
  }
  AllNormalCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename LatSet>
void BlockFieldManager<FieldType, FloatType, Dim>::SetupBoundary(
  datatype fromvalue, datatype voidvalue, datatype bdvalue) {
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    auto& field = _Fields[iblock];
    blockgeo.template SetupBoundary<FieldType, LatSet>(
      field, fromvalue, voidvalue, bdvalue);
    ++iblock;
  }
  AllNormalCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::ReadOctree(
  Octree<FloatType>* tree, datatype stlflag) {
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    auto& field = _Fields[iblock];
    blockgeo.template ReadOctree<FieldType>(tree, field, stlflag);
    ++iblock;
  }
  AllNormalCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename LatSet>
void BlockFieldManager<FieldType, FloatType, Dim>::CleanLonelyFlags(
  std::uint8_t flag, std::uint8_t voidflag, unsigned int lonelyth, bool recursive) {
  static_assert(std::is_same<datatype, std::uint8_t>::value,
    "CleanLonelyFlags only works for uint8_t");
  bool cleaned = false;
  // statistic info
  std::size_t total{};
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    auto& field = _Fields[iblock];
    blockgeo.template CleanLonelyFlags<FieldType, LatSet>(
      field, flag, voidflag, lonelyth, cleaned, total);
    ++iblock;
  }
  AllNormalCommunicate();

  if (cleaned) {
    MPI_RANK(0)
    std::cout << "[BlockFieldManager::CleanLonelyFlags]: " << total
              << " lonely cells cleaned" << std::endl;

    if (recursive) CleanLonelyFlags<LatSet>(flag, voidflag, lonelyth, recursive);
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename Func>
void BlockFieldManager<FieldType, FloatType, Dim>::forEachField(const Func& func) {
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    func(blockF);
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::InitValue(datatype initvalue) {
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    blockF.Init(initvalue);
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename Func>
void BlockFieldManager<FieldType, FloatType, Dim>::forEach(const Func& func) {
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    BlockField<FieldType, FloatType, Dim>& blockfield = _Fields[iblock];
    blockgeo.forEach([&blockfield, &func](std::size_t id) { func(blockfield, id); });
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename Func>
void BlockFieldManager<FieldType, FloatType, Dim>::forEachInner(const Func& func) {
  for (const auto& blockfield : this->getBlockFields()) {
    const auto& field = blockfield.getFieldType().getField(0);
    const auto& blockxd = blockfield.getBlock();
    const int overlap = blockxd.getOverlap();
    const int Nx = blockxd.getNx();
    const int Ny = blockxd.getNy();
    const int Nz = blockxd.getNz();

    if constexpr (Dim == 2) {
      for (int j = overlap; j < Ny - overlap; ++j) {
        for (int i = overlap; i < Nx - overlap; ++i) {
          std::size_t id = j * Nx + i;
          func(field[id]);
        }
      }
    } else if constexpr (Dim == 3) {
      std::size_t NxNy = Nx * Ny;
      for (int k = overlap; k < Nz - overlap; ++k) {
        for (int j = overlap; j < Ny - overlap; ++j) {
          for (int i = overlap; i < Nx - overlap; ++i) {
            std::size_t id = k * NxNy + j * Nx + i;
            func(field[id]);
          }
        }
      }
    }
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename Func, typename Func1>
void BlockFieldManager<FieldType, FloatType, Dim>::forEach_TransFlag(
  const Func& func, const Func1& func1) {
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    GenericArray<bool> TransFlag(blockgeo.getN(), false);
    BlockField<FieldType, FloatType, Dim>& blockfield = _Fields[iblock];
    for (std::size_t id = 0; id < blockgeo.getN(); ++id) {
      TransFlag[id] = func(blockfield, id);
    }
    for (std::size_t id = 0; id < blockgeo.getN(); ++id) {
      if (TransFlag[id]) {
        func1(blockfield, id);
      }
    }
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename Func>
void BlockFieldManager<FieldType, FloatType, Dim>::forEach(
  const AABB<FloatType, Dim>& AABBs, const Func& func) {
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    auto& field = _Fields[iblock];
    blockgeo.forEach(AABBs, [&field, &func](std::size_t id) { func(field, id); });
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename FlagFieldType, typename Func>
void BlockFieldManager<FieldType, FloatType, Dim>::forEach(
  const AABB<FloatType, Dim>& AABBs,
  const BlockFieldManager<FlagFieldType, FloatType, Dim>& FlagFManager, std::uint8_t flag,
  const Func& func) {
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    auto& field = _Fields[iblock];
    const auto& flagarr = FlagFManager.getBlockField(iblock).getField(0);
    blockgeo.forEach(
      AABBs, flagarr, flag, [&field, &func](std::size_t id) { func(field, id); });
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename FlagFieldType, typename Func>
void BlockFieldManager<FieldType, FloatType, Dim>::forEach(
  const BlockFieldManager<FlagFieldType, FloatType, Dim>& FlagFManager, std::uint8_t flag,
  const Func& func) {
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    auto& field = _Fields[iblock];
    const auto& flagarr = FlagFManager.getBlockField(iblock).getField(0);
    blockgeo.forEach(flagarr, flag, [&field, &func](std::size_t id) { func(field, id); });
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename FieldTypeX, typename Func>
void BlockFieldManager<FieldType, FloatType, Dim>::forEach(
  const BlockFieldManager<FieldTypeX, FloatType, Dim>& BlockFM, const Func& func) {
  int iblock = 0;
  for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
    auto& thisbfield = _Fields[iblock];
    const auto& bfield = BlockFM.getBlockField(iblock);
    blockgeo.forEach(
      [&thisbfield, &bfield, &func](std::size_t id) { func(thisbfield, bfield, id); });
    ++iblock;
  }
}


#ifdef MPI_ENABLED

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPINormalCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends =
        blockF.getBlock().getCommunicator().MPIComm.Sends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
      blockF.mpiNormalSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.Recvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.Recvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPIAverCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends =
        blockF.getBlock().getCommunicator().MPIComm.AverSends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
      blockF.mpiAverSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.AverRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.AverRecvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPIIntpCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends =
        blockF.getBlock().getCommunicator().MPIComm.IntpSends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
      blockF.mpiIntpSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.IntpRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.IntpRecvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}


template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPINormalCommunicate(
  std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends =
        blockF.getBlock().getCommunicator().MPIComm.Sends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
      blockF.mpiNormalSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.Recvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.Recvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPIAverCommunicate(
  std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data --- send data to lower level, deLevel+1
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel()) + 1;
    if (blockF.getBlock().getLevel() != std::uint8_t(0)) {
      if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
          blockF.getBlock().getCommunicator()._NeedMPIComm) {
        const std::vector<DistributedComm>& Sends =
          blockF.getBlock().getCommunicator().MPIComm.AverSends;
        SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
        blockF.mpiAverSend(SendRequests, SendBuffers[iblock], Sends);
      }
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.AverRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.AverRecvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPIIntpCommunicate(
  std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data --- send data to higher level, deLevel-1
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel()) - 1;
    if (deLevel != -1) {
      if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
          blockF.getBlock().getCommunicator()._NeedMPIComm) {
        const std::vector<DistributedComm>& Sends =
          blockF.getBlock().getCommunicator().MPIComm.IntpSends;
        SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
        blockF.mpiIntpSend(SendRequests, SendBuffers[iblock], Sends);
      }
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.IntpRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().MPIComm.IntpRecvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}


template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPINormalCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends =
        blockF.getBlock().getCommunicator().AllMPIComm.Sends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
      blockF.mpiNormalSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.Recvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.Recvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPIAverCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends =
        blockF.getBlock().getCommunicator().AllMPIComm.AverSends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
      blockF.mpiAverSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.AverRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.AverRecvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPIIntpCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends =
        blockF.getBlock().getCommunicator().AllMPIComm.IntpSends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
      blockF.mpiIntpSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.IntpRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    if (blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.IntpRecvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}


template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPINormalCommunicate(
  std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends =
        blockF.getBlock().getCommunicator().AllMPIComm.Sends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
      blockF.mpiNormalSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.Recvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.Recvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPIAverCommunicate(
  std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends =
        blockF.getBlock().getCommunicator().AllMPIComm.AverSends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
      blockF.mpiAverSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.AverRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.AverRecvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPIIntpCommunicate(
  std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<datatype>>> SendBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::vector<std::vector<std::vector<datatype>>> RecvBuffers(
    _Fields.size(), std::vector<std::vector<datatype>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends =
        blockF.getBlock().getCommunicator().AllMPIComm.IntpSends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<datatype>{});
      blockF.mpiIntpSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.IntpRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<datatype>{});
      blockF.mpiRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) &&
        blockF.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs =
        blockF.getBlock().getCommunicator().AllMPIComm.IntpRecvs;
      blockF.mpiSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}
#endif


template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::NormalCommunicate() {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    // normal communicate
    const Block<FloatType, Dim>& block = blockF.getBlock();
    const Communicator& communicator = block.getCommunicator();
    const std::vector<SharedComm>& comms = communicator.Comm.Comms;

    for (const SharedComm& comm : comms) {
      const BlockField<FieldType, FloatType, Dim>& nblockF =
        _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
      const std::size_t size = comm.SendRecvCells.size();
      for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
        const auto& nArray = nblockF.getField(iArr);
        auto& Array = blockF.getField(iArr);
        for (std::size_t i = 0; i < size; i += 2) {
          Array[comm.SendRecvCells[i + 1]] = nArray[comm.SendRecvCells[i]];
        }
      }
    }
  }

#ifdef MPI_ENABLED
  MPINormalCommunicate();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AverCommunicate() {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
  static constexpr std::size_t SendRecvPairSize = Dim == 2 ? 5 : 9;
  static constexpr std::size_t RecvOffset = Dim == 2 ? 4 : 8;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    // average communicate
    const Block<FloatType, Dim>& block = blockF.getBlock();
    const Communicator& communicator = block.getCommunicator();
    const std::vector<SharedComm>& comms = communicator.Comm.AverComm;

    for (const SharedComm& comm : comms) {
      const BlockField<FieldType, FloatType, Dim>& nblockF =
        _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
      const std::size_t size = comm.SendRecvCells.size();
      for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
        const auto& nArray = nblockF.getField(iArr);
        auto& Array = blockF.getField(iArr);
        for (std::size_t i = 0; i < size; i += SendRecvPairSize) {
          Array[comm.SendRecvCells[i + RecvOffset]] =
            getAverage<FloatType, Dim>(nArray, i, comm.SendRecvCells);
        }
      }
    }
  }

#ifdef MPI_ENABLED
  MPIAverCommunicate();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::IntpCommunicate() {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
  static constexpr std::size_t SendRecvPairSize = Dim == 2 ? 5 : 9;
  static constexpr std::size_t RecvOffset = Dim == 2 ? 4 : 8;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    // interp communicate
    const Block<FloatType, Dim>& block = blockF.getBlock();
    const Communicator& communicator = block.getCommunicator();
    const std::vector<SharedComm>& comms = communicator.Comm.IntpComm;

    for (const SharedComm& comm : comms) {
      const BlockField<FieldType, FloatType, Dim>& nblockF =
        _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
      const std::size_t size = comm.SendRecvCells.size();
      for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
        const auto& nArray = nblockF.getField(iArr);
        auto& Array = blockF.getField(iArr);
        for (std::size_t id = 0; id < size;) {
          Interpolation<FloatType, Dim>(Array, nArray, comm.SendRecvCells, id);
        }
      }
    }
  }

#ifdef MPI_ENABLED
  MPIIntpCommunicate();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::Communicate() {
  NormalCommunicate();
  AverCommunicate();
  IntpCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::NormalCommunicate(std::int64_t count) {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      // normal communicate
      const Block<FloatType, Dim>& block = blockF.getBlock();
      const Communicator& communicator = block.getCommunicator();
      const std::vector<SharedComm>& comms = communicator.Comm.Comms;

      for (const SharedComm& comm : comms) {
        const BlockField<FieldType, FloatType, Dim>& nblockF =
          _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
        const std::size_t size = comm.SendRecvCells.size();
        for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
          const auto& nArray = nblockF.getField(iArr);
          auto& Array = blockF.getField(iArr);
          for (std::size_t i = 0; i < size; i += 2) {
            Array[comm.SendRecvCells[i + 1]] = nArray[comm.SendRecvCells[i]];
          }
        }
      }
    }
  }

#ifdef MPI_ENABLED
  MPINormalCommunicate(count);
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AverCommunicate(std::int64_t count) {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
  static constexpr std::size_t SendRecvPairSize = Dim == 2 ? 5 : 9;
  static constexpr std::size_t RecvOffset = Dim == 2 ? 4 : 8;
  static constexpr FloatType weight = Dim == 2 ? FloatType(0.25) : FloatType(0.125);
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      // average communicate
      const Block<FloatType, Dim>& block = blockF.getBlock();
      const Communicator& communicator = block.getCommunicator();
      const std::vector<SharedComm>& comms = communicator.Comm.AverComm;

      for (const SharedComm& comm : comms) {
        const BlockField<FieldType, FloatType, Dim>& nblockF =
          _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
        const std::size_t size = comm.SendRecvCells.size();
        for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
          const auto& nArray = nblockF.getField(iArr);
          auto& Array = blockF.getField(iArr);
          for (std::size_t i = 0; i < size; i += SendRecvPairSize) {
            Array[comm.SendRecvCells[i + RecvOffset]] =
              getAverage<FloatType, Dim>(nArray, i, comm.SendRecvCells);
          }
        }
      }
    }
  }

#ifdef MPI_ENABLED
  MPIAverCommunicate(count);
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::IntpCommunicate(std::int64_t count) {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
  static constexpr std::size_t SendRecvPairSize = Dim == 2 ? 5 : 9;
  static constexpr std::size_t RecvOffset = Dim == 2 ? 4 : 8;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      // interp communicate
      const Block<FloatType, Dim>& block = blockF.getBlock();
      const Communicator& communicator = block.getCommunicator();
      const std::vector<SharedComm>& comms = communicator.Comm.IntpComm;

      for (const SharedComm& comm : comms) {
        const BlockField<FieldType, FloatType, Dim>& nblockF =
          _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
        const std::size_t size = comm.SendRecvCells.size();
        for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
          const auto& nArray = nblockF.getField(iArr);
          auto& Array = blockF.getField(iArr);
          for (std::size_t id = 0; id < size;) {
            Interpolation<FloatType, Dim>(Array, nArray, comm.SendRecvCells, id);
          }
        }
      }
    }
  }

#ifdef MPI_ENABLED
  MPIIntpCommunicate(count);
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::Communicate(std::int64_t count) {
  NormalCommunicate(count);
  AverCommunicate(count);
  IntpCommunicate(count);
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllNormalCommunicate() {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    // normal communicate
    const Block<FloatType, Dim>& block = blockF.getBlock();
    const Communicator& communicator = block.getCommunicator();
    const std::vector<SharedComm>& comms = communicator.AllComm.Comms;

    for (const SharedComm& comm : comms) {
      const BlockField<FieldType, FloatType, Dim>& nblockF =
        _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
      const std::size_t size = comm.SendRecvCells.size();
      for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
        const auto& nArray = nblockF.getField(iArr);
        auto& Array = blockF.getField(iArr);
        for (std::size_t i = 0; i < size; i += 2) {
          Array[comm.SendRecvCells[i + 1]] = nArray[comm.SendRecvCells[i]];
        }
      }
    }
  }

#ifdef MPI_ENABLED
  AllMPINormalCommunicate();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllAverCommunicate() {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
  static constexpr std::size_t SendRecvPairSize = Dim == 2 ? 5 : 9;
  static constexpr std::size_t RecvOffset = Dim == 2 ? 4 : 8;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    // average communicate
    const Block<FloatType, Dim>& block = blockF.getBlock();
    const Communicator& communicator = block.getCommunicator();
    const std::vector<SharedComm>& comms = communicator.AllComm.AverComm;

    for (const SharedComm& comm : comms) {
      const BlockField<FieldType, FloatType, Dim>& nblockF =
        _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
      const std::size_t size = comm.SendRecvCells.size();
      for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
        const auto& nArray = nblockF.getField(iArr);
        auto& Array = blockF.getField(iArr);
        for (std::size_t i = 0; i < size; i += SendRecvPairSize) {
          Array[comm.SendRecvCells[i + RecvOffset]] =
            getAverage<FloatType, Dim>(nArray, i, comm.SendRecvCells);
        }
      }
    }
  }

#ifdef MPI_ENABLED
  AllMPIAverCommunicate();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllIntpCommunicate() {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
  static constexpr std::size_t SendRecvPairSize = Dim == 2 ? 5 : 9;
  static constexpr std::size_t RecvOffset = Dim == 2 ? 4 : 8;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    // interp communicate
    const Block<FloatType, Dim>& block = blockF.getBlock();
    const Communicator& communicator = block.getCommunicator();
    const std::vector<SharedComm>& comms = communicator.AllComm.IntpComm;

    for (const SharedComm& comm : comms) {
      const BlockField<FieldType, FloatType, Dim>& nblockF =
        _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
      const std::size_t size = comm.SendRecvCells.size();
      for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
        const auto& nArray = nblockF.getField(iArr);
        auto& Array = blockF.getField(iArr);
        for (std::size_t id = 0; id < size;) {
          Interpolation<FloatType, Dim>(Array, nArray, comm.SendRecvCells, id);
        }
      }
    }
  }

#ifdef MPI_ENABLED
  AllMPIIntpCommunicate();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllCommunicate() {
  AllNormalCommunicate();
  AllAverCommunicate();
  AllIntpCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllNormalCommunicate(
  std::int64_t count) {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      // normal communicate
      const Block<FloatType, Dim>& block = blockF.getBlock();
      const Communicator& communicator = block.getCommunicator();
      const std::vector<SharedComm>& comms = communicator.AllComm.Comms;

      for (const SharedComm& comm : comms) {
        const BlockField<FieldType, FloatType, Dim>& nblockF =
          _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
        const std::size_t size = comm.SendRecvCells.size();
        for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
          const auto& nArray = nblockF.getField(iArr);
          auto& Array = blockF.getField(iArr);
          for (std::size_t i = 0; i < size; i += 2) {
            Array[comm.SendRecvCells[i + 1]] = nArray[comm.SendRecvCells[i]];
          }
        }
      }
    }
  }

#ifdef MPI_ENABLED
  AllMPINormalCommunicate(count);
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllAverCommunicate(
  std::int64_t count) {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
  static constexpr std::size_t SendRecvPairSize = Dim == 2 ? 5 : 9;
  static constexpr std::size_t RecvOffset = Dim == 2 ? 4 : 8;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      // average communicate
      const Block<FloatType, Dim>& block = blockF.getBlock();
      const Communicator& communicator = block.getCommunicator();
      const std::vector<SharedComm>& comms = communicator.AllComm.AverComm;

      for (const SharedComm& comm : comms) {
        const BlockField<FieldType, FloatType, Dim>& nblockF =
          _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
        const std::size_t size = comm.SendRecvCells.size();
        for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
          const auto& nArray = nblockF.getField(iArr);
          auto& Array = blockF.getField(iArr);
          for (std::size_t i = 0; i < size; i += SendRecvPairSize) {
            Array[comm.SendRecvCells[i + RecvOffset]] =
              getAverage<FloatType, Dim>(nArray, i, comm.SendRecvCells);
          }
        }
      }
    }
  }

#ifdef MPI_ENABLED
  AllMPIAverCommunicate(count);
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllIntpCommunicate(
  std::int64_t count) {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
    const int deLevel =
      static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      // interp communicate
      const Block<FloatType, Dim>& block = blockF.getBlock();
      const Communicator& communicator = block.getCommunicator();
      const std::vector<SharedComm>& comms = communicator.AllComm.IntpComm;

      for (const SharedComm& comm : comms) {
        const BlockField<FieldType, FloatType, Dim>& nblockF =
          _Fields[_BlockGeo.findBlockIndex(comm.SendBlockId)];
        const std::size_t size = comm.SendRecvCells.size();
        for (unsigned int iArr = 0; iArr < ArrayDim; ++iArr) {
          const auto& nArray = nblockF.getField(iArr);
          auto& Array = blockF.getField(iArr);
          for (std::size_t id = 0; id < size;) {
            Interpolation<FloatType, Dim>(Array, nArray, comm.SendRecvCells, id);
          }
        }
      }
    }
  }

#ifdef MPI_ENABLED
  AllMPIIntpCommunicate(count);
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllCommunicate(std::int64_t count) {
  AllNormalCommunicate(count);
  AllAverCommunicate(count);
  AllIntpCommunicate(count);
}
