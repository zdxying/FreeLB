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
// shared plan post/wait/unpack runner
#include "parallel/plan_runner.h"


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
template <typename buffer_type>
void BlockField<FieldType, FloatType, Dim>::packSend(
  std::size_t commIdx, buffer_type& buffer, std::size_t offset, bool all) {
  const std::vector<DistributedComm>& sends =
    all ? _Block.getCommunicator().AllMPIComm.Sends : _Block.getCommunicator().MPIComm.Sends;
  const std::vector<std::size_t>& cells = sends[commIdx].Cells;
  for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
    const auto& Array = FieldType::getField(iArr);
    for (std::size_t id : cells) {
      buffer[offset] = Array[id];
      ++offset;
    }
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename buffer_type>
void BlockField<FieldType, FloatType, Dim>::unpackRecv(
  std::size_t commIdx, const buffer_type& buffer, std::size_t offset, bool all) {
  const std::vector<DistributedComm>& recvs =
    all ? _Block.getCommunicator().AllMPIComm.Recvs : _Block.getCommunicator().MPIComm.Recvs;
  const std::vector<std::size_t>& cells = recvs[commIdx].Cells;
  for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
    auto& Array = FieldType::getField(iArr);
    for (std::size_t id : cells) {
      Array[id] = buffer[offset];
      ++offset;
    }
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <bool Coarse, typename buffer_type>
void BlockField<FieldType, FloatType, Dim>::packCrossSend(
  const DistributedComm& comm, buffer_type& buffer, std::size_t offset) {
  // Aver packs the 2^Dim-cell average of each coarse cell, Intp the
  // interpolated value of each refined cell; the layout is the same.
  static constexpr std::size_t SendSize = flbcomm::AverSrcNum(Dim);
  const std::vector<std::size_t>& sends = comm.Cells;
  std::size_t bufidx = offset;
  for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
    const auto& Array = FieldType::getField(iArr);
    if constexpr (Coarse) {
      for (std::size_t id = 0; id < sends.size(); id += SendSize) {
        buffer[bufidx] = getAverage<FloatType, Dim>(Array, id, sends);
        ++bufidx;
      }
    } else {
      for (std::size_t id = 0; id < sends.size();) {
        Interpolation<FloatType, Dim>(buffer, bufidx, Array, sends, id);
      }
    }
  }
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename buffer_type>
void BlockField<FieldType, FloatType, Dim>::unpackCrossRecv(
  const DistributedComm& comm, const buffer_type& buffer, std::size_t offset) {
  const std::vector<std::size_t>& recvs = comm.Cells;
  std::size_t bufidx = offset;
  for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
    auto& Array = FieldType::getField(iArr);
    for (std::size_t id : recvs) {
      Array[id] = buffer[bufidx];
      ++bufidx;
    }
  }
}

#endif


// BlockFieldManager

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockFieldManager<FieldType, FloatType, Dim>::BlockFieldManager(
  BlockGeometry<FloatType, Dim>& blockgeometry)
    : _BlockGeo(blockgeometry), _tagBase(flbcomm::allocateFieldTagBase()) {
  _Fields.reserve(_BlockGeo.getBlocks().size());
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    _Fields.emplace_back(block);
  }
  InitDeviceData();
  _fieldsGeneration = _BlockGeo.getGeometryGeneration();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockFieldManager<FieldType, FloatType, Dim>::BlockFieldManager(
  BlockGeometry<FloatType, Dim>& blockgeometry, datatype initvalue)
    : _BlockGeo(blockgeometry), _tagBase(flbcomm::allocateFieldTagBase()) {
  _Fields.reserve(_BlockGeo.getBlocks().size());
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    _Fields.emplace_back(block, initvalue);
  }
  InitDeviceData();
  _fieldsGeneration = _BlockGeo.getGeometryGeneration();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockFieldManager<FieldType, FloatType, Dim>::BlockFieldManager(
  const BlockFieldManager& blockFManager)
    : _Fields(blockFManager._Fields),
      _fieldsGeneration(blockFManager._fieldsGeneration),
      _BlockGeo(blockFManager._BlockGeo), _tagBase(blockFManager._tagBase) {
  InitDeviceData();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
BlockFieldManager<FieldType, FloatType, Dim>::BlockFieldManager(
  BlockFieldManager&& blockFManager) noexcept
    : _Fields(std::move(blockFManager._Fields)),
      _fieldsGeneration(blockFManager._fieldsGeneration),
      _BlockGeo(blockFManager._BlockGeo), _tagBase(blockFManager._tagBase) {
#ifdef __CUDACC__
  // The source keeps its own table/mirror; its destructor frees them.  The
  // moved-to manager builds a fresh table from the block fields it now owns.
  InitDeviceData();
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
#ifdef MPI_ENABLED
  resetCommPlans();
#endif
  if (canReuseFields()) {
    resetFields(datatype{});
    return;
  }
  _Fields.clear();
  _Fields.reserve(_BlockGeo.getBlocks().size());
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    _Fields.emplace_back(block);
  }
  _fieldsGeneration = _BlockGeo.getGeometryGeneration();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::Init(datatype initvalue) {
#ifdef MPI_ENABLED
  resetCommPlans();
#endif
  if (canReuseFields()) {
    resetFields(initvalue);
    return;
  }
  _Fields.clear();
  _Fields.reserve(_BlockGeo.getBlocks().size());
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    _Fields.emplace_back(block, initvalue);
  }
  _fieldsGeneration = _BlockGeo.getGeometryGeneration();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::NonFieldInit(datatype initvalue) {
#ifdef MPI_ENABLED
  resetCommPlans();
#endif
  if (canReuseFields()) {
    resetFields(initvalue);
    return;
  }
  _Fields.clear();
  _Fields.reserve(_BlockGeo.getBlocks().size());
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    _Fields.emplace_back(block, initvalue);
  }
  _fieldsGeneration = _BlockGeo.getGeometryGeneration();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::Init(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper) {
  // The old set must stay alive until the transfer finishes, so peak storage
  // is old + new; this is inherent to an arbitrary block remap.  reserve()
  // removes the extra vector-growth copies on top of it.
  std::vector<BlockField<FieldType, FloatType, Dim>> NewFields;
  NewFields.reserve(_BlockGeo.getBlocks().size());
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    NewFields.emplace_back(block);
  }
  // data transfer
  FieldDataTransfer(GeoHelper, NewFields);
  _Fields.swap(NewFields);
  _fieldsGeneration = _BlockGeo.getGeometryGeneration();
#ifdef MPI_ENABLED
  resetCommPlans();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::Init(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue) {
  // See the note in the overload above: old and new coexist during transfer.
  std::vector<BlockField<FieldType, FloatType, Dim>> NewFields;
  NewFields.reserve(_BlockGeo.getBlocks().size());
  for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
    NewFields.emplace_back(block, initvalue);
  }
  // data transfer
  FieldDataTransfer(GeoHelper, NewFields);
  _Fields.swap(NewFields);
  _fieldsGeneration = _BlockGeo.getGeometryGeneration();
#ifdef MPI_ENABLED
  resetCommPlans();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
template <typename FlagFieldType, typename Func>
void BlockFieldManager<FieldType, FloatType, Dim>::InitCopy(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue,
  const BlockFieldManager<FlagFieldType, FloatType, Dim>& FlagFManager, std::uint8_t flag,
  const Func& func) {
  std::vector<BlockField<FieldType, FloatType, Dim>> NewFields;
  NewFields.reserve(_BlockGeo.getBlocks().size());
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
  _fieldsGeneration = _BlockGeo.getGeometryGeneration();
#ifdef MPI_ENABLED
  resetCommPlans();
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::FieldDataTransfer(
  BlockGeometryHelper<FloatType, Dim>& GeoHelper,
  std::vector<BlockField<FieldType, FloatType, Dim>>& NewFields) {
  static_assert(Dim == 2, "FieldDataTransfer only support 2D");
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

// Build the coalesced plan for the plain normal (all == false) or all
// (all == true) exchange.  Same ordering rule as the lattice
// full exchange: the links are keyed by (receiver block, sender block), which
// is the same pair on both ends and unique per link.  The tag comes from this
// manager's reserved range, so the plan can be in flight while another field
// or a lattice exchange is pending.
template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::initFieldPlan(bool all,
                                                                std::size_t depth) {
  auto commsOf = [&](std::size_t bi, bool send) -> const std::vector<DistributedComm>& {
    const auto& comm = _Fields[bi].getBlock().getCommunicator();
    if (all) return send ? comm.AllMPIComm.Sends : comm.AllMPIComm.Recvs;
    return send ? comm.MPIComm.Sends : comm.MPIComm.Recvs;
  };
  // Field exchanges use their own tag range so that a field exchange may be
  // in flight while a lattice exchange is posted (the hidden step): both are
  // then pending on MPI_COMM_WORLD, and MPI matches on (communicator, source,
  // tag) alone.  Each manager owns two slots; normal and all use +0/+1, so two
  // field managers can also be in flight at the same time.
  flbcomm::CommPlan<buffer_type>& plan = all ? _allPlans[depth] : _normalPlans[depth];
  flbcomm::buildPlan<buffer_type>(
    _Fields.size(), static_cast<std::size_t>(_BlockGeo.getMaxLevel()), depth,
    _BlockGeo.getGeometryGeneration(),
    [&](std::size_t bi) { return _Fields[bi].getBlock().getLevel(); },
    [&](std::size_t bi, bool send, std::vector<flbcomm::CommRef>& refs) {
      const std::vector<DistributedComm>& comms = commsOf(bi, send);
      const int self = _Fields[bi].getBlock().getBlockId();
      for (std::size_t ci = 0; ci < comms.size(); ++ci) {
        const int peerBlock = comms[ci].TargetBlockId;
        const int tag = _tagBase + (all ? 1 : 0);
        const int receiver = send ? peerBlock : self;
        const int sender = send ? self : peerBlock;
        refs.push_back(flbcomm::CommRef{
          bi, ci, 0, comms[ci].TargetRank, tag,
          flbcomm::LinkKey{receiver, sender, 0, 0},
          comms[ci].Cells.size() * FieldType::array_dim});
      }
    },
    plan);
}

// Same lazy per-depth lookup as the lattice: grow the vector up to `depth` on
// first use and build that entry once.
template <typename FieldType, typename FloatType, unsigned int Dim>
flbcomm::CommPlan<typename BlockFieldManager<FieldType, FloatType, Dim>::buffer_type>&
BlockFieldManager<FieldType, FloatType, Dim>::planFor(bool all, std::size_t depth) {
  std::vector<flbcomm::CommPlan<buffer_type>>& plans = all ? _allPlans : _normalPlans;
  const std::size_t maxDepth = static_cast<std::size_t>(_BlockGeo.getMaxLevel());
  if (depth > maxDepth) depth = maxDepth;
  if (plans.size() <= depth) plans.resize(depth + 1);
  flbcomm::CommPlan<buffer_type>& plan = plans[depth];
  if (!plan.Ready) initFieldPlan(all, depth);
  checkPlanGeneration(plan);
  return plan;
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::postFieldPlan(
  flbcomm::CommPlan<buffer_type>& plan, bool all,
  std::vector<MPI_Request>& sendRequests, std::vector<MPI_Request>& recvRequests) {
  flbcomm::postPlan(
    plan, sendRequests, recvRequests,
    [&](const flbcomm::CommRef& r, buffer_type& buffer, std::size_t offset) {
      _Fields[r.Block].packSend(r.Comm, buffer, offset, all);
    });
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::finishFieldPlan(
  flbcomm::CommPlan<buffer_type>& plan, bool all,
  std::vector<MPI_Request>& sendRequests, std::vector<MPI_Request>& recvRequests) {
  flbcomm::waitSends(sendRequests);
  flbcomm::finishRecvs(
    plan, recvRequests,
    [&](const flbcomm::CommRef& r, const buffer_type& buffer, std::size_t offset) {
      _Fields[r.Block].unpackRecv(r.Comm, buffer, offset, all);
    });
  sendRequests.clear();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::runFieldPlan(
  flbcomm::CommPlan<buffer_type>& plan, bool all) {
  std::vector<MPI_Request> sendRequests, recvRequests;
  postFieldPlan(plan, all, sendRequests, recvRequests);
  finishFieldPlan(plan, all, sendRequests, recvRequests);
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPINormalCommunicate(
  std::int64_t count) {
  const std::size_t depth =
    flbcomm::activeDepth(count, static_cast<std::size_t>(_BlockGeo.getMaxLevel()));
  runFieldPlan(planFor(false, depth), false);
}




// Aver (fine -> coarse) and Intp (coarse -> fine) differ only in the conversion
// and the direction of the level filter.  The plan shape matches the lattice:
// one message per peer, refs keyed by (receiver block, sender block), a fixed
// tag, and independent send/receive level filters.
template <typename FieldType, typename FloatType, unsigned int Dim>
const std::vector<DistributedComm>&
BlockFieldManager<FieldType, FloatType, Dim>::crossComms(flbcomm::CrossMode mode,
                                                         std::size_t block, bool send) {
  const auto& comm = _Fields[block].getBlock().getCommunicator();
  const bool aver = flbcomm::crossIsAver(mode);
  const bool all = flbcomm::crossIsAll(mode);
  if (aver) {
    return all ? (send ? comm.AllMPIComm.AverSends : comm.AllMPIComm.AverRecvs)
               : (send ? comm.MPIComm.AverSends : comm.MPIComm.AverRecvs);
  }
  return all ? (send ? comm.AllMPIComm.IntpSends : comm.AllMPIComm.IntpRecvs)
             : (send ? comm.MPIComm.IntpSends : comm.MPIComm.IntpRecvs);
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::initCrossPlan(flbcomm::CrossMode mode,
                                                                std::size_t depth) {
  const bool aver = flbcomm::crossIsAver(mode);
  const bool all = flbcomm::crossIsAll(mode);
  std::vector<flbcomm::CommPlan<buffer_type>>& plans =
    aver ? (all ? _allAverPlans : _averPlans)
         : (all ? _allIntpPlans : _intpPlans);
  flbcomm::CommPlan<buffer_type>& plan = plans[depth];
  const std::size_t maxLevel = static_cast<std::size_t>(_BlockGeo.getMaxLevel());
  const int tag = _tagBase + flbcomm::fieldTagSlot(mode);
  flbcomm::buildPlan<buffer_type>(
    _Fields.size(), _BlockGeo.getGeometryGeneration(),
    [&](std::size_t bi, bool send) {
      const int deLevel = static_cast<int>(maxLevel - _Fields[bi].getBlock().getLevel());
      if (aver) {
        return send ? (deLevel + 1 <= static_cast<int>(depth) &&
                       _Fields[bi].getBlock().getLevel() != 0)
                    : (deLevel <= static_cast<int>(depth));
      }
      return send ? (deLevel - 1 <= static_cast<int>(depth) && deLevel != 0)
                  : (deLevel <= static_cast<int>(depth));
    },
    [&](std::size_t bi, bool send, std::vector<flbcomm::CommRef>& refs) {
      const std::vector<DistributedComm>& comms = crossComms(mode, bi, send);
      const int self = _Fields[bi].getBlock().getBlockId();
      for (std::size_t ci = 0; ci < comms.size(); ++ci) {
        const int peerBlock = comms[ci].TargetBlockId;
        const int receiver = send ? peerBlock : self;
        const int sender = send ? self : peerBlock;
        refs.push_back(flbcomm::CommRef{
          bi, ci, 0, comms[ci].TargetRank, tag,
          flbcomm::LinkKey{receiver, sender, 0, 0},
          send ? _Fields[bi].crossSendSize(comms[ci])
               : _Fields[bi].crossRecvSize(comms[ci])});
      }
    },
    plan);
}

template <typename FieldType, typename FloatType, unsigned int Dim>
flbcomm::CommPlan<typename BlockFieldManager<FieldType, FloatType, Dim>::buffer_type>&
BlockFieldManager<FieldType, FloatType, Dim>::crossPlanFor(flbcomm::CrossMode mode,
                                                           std::size_t depth) {
  const bool aver = flbcomm::crossIsAver(mode);
  const bool all = flbcomm::crossIsAll(mode);
  std::vector<flbcomm::CommPlan<buffer_type>>& plans =
    aver ? (all ? _allAverPlans : _averPlans)
         : (all ? _allIntpPlans : _intpPlans);
  const std::size_t maxDepth = static_cast<std::size_t>(_BlockGeo.getMaxLevel());
  if (depth > maxDepth) depth = maxDepth;
  if (plans.size() <= depth) plans.resize(depth + 1);
  flbcomm::CommPlan<buffer_type>& plan = plans[depth];
  if (!plan.Ready) initCrossPlan(mode, depth);
  checkPlanGeneration(plan);
  return plan;
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::runCrossPlan(flbcomm::CrossMode mode,
                                                                std::size_t depth) {
  flbcomm::CommPlan<buffer_type>& plan = crossPlanFor(mode, depth);
  const bool aver = flbcomm::crossIsAver(mode);
  flbcomm::runPlan(
    plan,
    [&](const flbcomm::CommRef& r, buffer_type& buffer, std::size_t offset) {
      const auto& comms = crossComms(mode, r.Block, true);
      if (aver) {
        _Fields[r.Block].template packCrossSend<true>(comms[r.Comm], buffer, offset);
      } else {
        _Fields[r.Block].template packCrossSend<false>(comms[r.Comm], buffer, offset);
      }
    },
    [&](const flbcomm::CommRef& r, const buffer_type& buffer, std::size_t offset) {
      const auto& comms = crossComms(mode, r.Block, false);
      _Fields[r.Block].unpackCrossRecv(comms[r.Comm], buffer, offset);
    });
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPIAverCommunicate() {
  runCrossPlan(flbcomm::CrossMode::Aver,
               static_cast<std::size_t>(_BlockGeo.getMaxLevel()));
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPIAverCommunicate(
  std::int64_t count) {
  runCrossPlan(flbcomm::CrossMode::Aver,
               flbcomm::activeDepth(count, static_cast<std::size_t>(_BlockGeo.getMaxLevel())));
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPIIntpCommunicate() {
  runCrossPlan(flbcomm::CrossMode::Intp,
               static_cast<std::size_t>(_BlockGeo.getMaxLevel()));
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::MPIIntpCommunicate(
  std::int64_t count) {
  runCrossPlan(flbcomm::CrossMode::Intp,
               flbcomm::activeDepth(count, static_cast<std::size_t>(_BlockGeo.getMaxLevel())));
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPIAverCommunicate() {
  runCrossPlan(flbcomm::CrossMode::AllAver,
               static_cast<std::size_t>(_BlockGeo.getMaxLevel()));
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPIIntpCommunicate() {
  runCrossPlan(flbcomm::CrossMode::AllIntp,
               static_cast<std::size_t>(_BlockGeo.getMaxLevel()));
}


template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPINormalCommunicate(
  std::int64_t count) {
  const std::size_t depth =
    flbcomm::activeDepth(count, static_cast<std::size_t>(_BlockGeo.getMaxLevel()));
  runFieldPlan(planFor(true, depth), true);
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPIAverCommunicate(
  std::int64_t count) {
  runCrossPlan(flbcomm::CrossMode::AllAver,
               flbcomm::activeDepth(count, static_cast<std::size_t>(_BlockGeo.getMaxLevel())));
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllMPIIntpCommunicate(
  std::int64_t count) {
  runCrossPlan(flbcomm::CrossMode::AllIntp,
               flbcomm::activeDepth(count, static_cast<std::size_t>(_BlockGeo.getMaxLevel())));
}
#endif


template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::PostNormalCommunicate() {
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
  assert(!_normalPosted);
  postFieldPlan(planFor(false, static_cast<std::size_t>(_BlockGeo.getMaxLevel())), false,
                _normalSendRequests, _normalRecvRequests);
  _normalPosted = true;
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::FinishNormalCommunicate() {
#ifdef MPI_ENABLED
  assert(_normalPosted);
  finishFieldPlan(planFor(false, static_cast<std::size_t>(_BlockGeo.getMaxLevel())), false,
                  _normalSendRequests, _normalRecvRequests);
  _normalPosted = false;
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::NormalCommunicate() {
  PostNormalCommunicate();
  FinishNormalCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AverCommunicate() {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
  static constexpr std::size_t SendRecvPairSize = flbcomm::AverPairStride(Dim);
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(Dim);
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
          flbcomm::checkCellIndex(comm.SendRecvCells[i + RecvOffset], blockF.getBlock().getN(),
                                  block.getBlockId(), comm.SendBlockId,
                                  "BlockFieldManager::AverCommunicate");
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
  static constexpr std::size_t SendRecvPairSize = flbcomm::AverPairStride(Dim);
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(Dim);
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
          // Interpolation() writes through recvs[recvidx++], i.e. the slot right
          // after the sources; check it before the call, the same slot the Aver
          // path indexes explicitly.
          flbcomm::checkCellIndex(comm.SendRecvCells[id + RecvOffset], block.getN(),
                                  block.getBlockId(), comm.SendBlockId,
                                  "BlockFieldManager::IntpCommunicate");
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
  static constexpr std::size_t SendRecvPairSize = flbcomm::AverPairStride(Dim);
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(Dim);
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
            flbcomm::checkCellIndex(comm.SendRecvCells[i + RecvOffset],
                                    blockF.getBlock().getN(), blockF.getBlock().getBlockId(),
                                    comm.SendBlockId,
                                    "BlockFieldManager::AverCommunicate(count)");
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
  static constexpr std::size_t SendRecvPairSize = flbcomm::AverPairStride(Dim);
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(Dim);
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
            flbcomm::checkCellIndex(comm.SendRecvCells[id + RecvOffset], block.getN(),
                                    block.getBlockId(), comm.SendBlockId,
                                    "BlockFieldManager::IntpCommunicate(count)");
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
void BlockFieldManager<FieldType, FloatType, Dim>::PostAllNormalCommunicate() {
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
  assert(!_allPosted);
  postFieldPlan(planFor(true, static_cast<std::size_t>(_BlockGeo.getMaxLevel())), true,
                _allSendRequests, _allRecvRequests);
  _allPosted = true;
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::FinishAllNormalCommunicate() {
#ifdef MPI_ENABLED
  assert(_allPosted);
  finishFieldPlan(planFor(true, static_cast<std::size_t>(_BlockGeo.getMaxLevel())), true,
                  _allSendRequests, _allRecvRequests);
  _allPosted = false;
#endif
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllNormalCommunicate() {
  PostAllNormalCommunicate();
  FinishAllNormalCommunicate();
}

template <typename FieldType, typename FloatType, unsigned int Dim>
void BlockFieldManager<FieldType, FloatType, Dim>::AllAverCommunicate() {
  static constexpr unsigned int ArrayDim = FieldType::array_dim;
  static constexpr std::size_t SendRecvPairSize = flbcomm::AverPairStride(Dim);
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(Dim);
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
          flbcomm::checkCellIndex(comm.SendRecvCells[i + RecvOffset], blockF.getBlock().getN(),
                                  blockF.getBlock().getBlockId(), comm.SendBlockId,
                                  "BlockFieldManager::AllAverCommunicate");
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
  static constexpr std::size_t SendRecvPairSize = flbcomm::AverPairStride(Dim);
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(Dim);
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
          flbcomm::checkCellIndex(comm.SendRecvCells[id + RecvOffset], block.getN(),
                                  block.getBlockId(), comm.SendBlockId,
                                  "BlockFieldManager::AllIntpCommunicate");
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
  static constexpr std::size_t SendRecvPairSize = flbcomm::AverPairStride(Dim);
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(Dim);
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
            flbcomm::checkCellIndex(comm.SendRecvCells[i + RecvOffset],
                                    blockF.getBlock().getN(), blockF.getBlock().getBlockId(),
                                    comm.SendBlockId,
                                    "BlockFieldManager::AllAverCommunicate(count)");
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
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(Dim);
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
            flbcomm::checkCellIndex(comm.SendRecvCells[id + RecvOffset], block.getN(),
                                    block.getBlockId(), comm.SendBlockId,
                                    "BlockFieldManager::AllIntpCommunicate(count)");
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
