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

// field_struct.h

#pragma once

#include "data_struct/field.h"
#include "utils/alias.h"

template <typename FieldType, typename dataType>
class BlockFieldStruct {
 private:
  // block field data
  std::vector<FieldType> _Data;
  // block field size
  std::vector<std::size_t> _Sizes;

 public:
  BlockFieldStruct(std::vector<std::size_t> size) : _Sizes(size) {
    for (int i = 0; i < _Sizes.size(); ++i) _Data.emplace_back(_Sizes[i]);
  }
  BlockFieldStruct(std::vector<std::size_t> size, dataType initvalue) : _Sizes(size) {
    for (int i = 0; i < _Sizes.size(); ++i) _Data.emplace_back(_Sizes[i], initvalue);
  }
  std::size_t getBlockNum() const { return _Sizes.size(); }
  FieldType& getBlockField(int i) { return _Data[i]; }
  const FieldType& getBlockField(int i) const { return _Data[i]; }
  void Erase(int i) {
    _Data.erase(_Data.begin() + i);
    _Sizes.erase(_Sizes.begin() + i);
  }
  void Pushback(std::size_t size) {
    _Sizes.push_back(size);
    _Data.emplace_back(size);
  }
  void Pushback(std::size_t size, dataType initvalue) {
    _Sizes.push_back(size);
    _Data.emplace_back(size, initvalue);
  }
  void Clear() {
    _Data.clear();
    _Sizes.clear();
  }
};

template <typename T, unsigned int D>
struct BlockComm;

template <typename T, unsigned int D>
struct InterpBlockComm;

template <typename FieldType, typename FloatType, unsigned int Dim>
class BlockField;

template <typename FieldType, typename FloatType, unsigned int Dim>
struct BlockFieldComm {
  BlockField<FieldType, FloatType, Dim>* BlockF;
  BlockComm<FloatType, Dim>* Comm;

  BlockFieldComm(BlockField<FieldType, FloatType, Dim>* blockF,
                 BlockComm<FloatType, Dim>* comm)
      : BlockF(blockF), Comm(comm) {}

  const std::vector<std::size_t>& getSends() const { return Comm->SendCells; }
  const std::size_t& getSend(std::size_t i) const { return (Comm->SendCells)[i]; }

  const std::vector<std::size_t>& getRecvs() const { return Comm->RecvCells; }
  const std::size_t& getRecv(std::size_t i) const { return (Comm->RecvCells)[i]; }
};

template <typename FieldType, typename FloatType, unsigned int Dim>
struct InterpBlockFieldComm {
  BlockField<FieldType, FloatType, Dim>* BlockF;
  InterpBlockComm<FloatType, Dim>* Comm;

  InterpBlockFieldComm(BlockField<FieldType, FloatType, Dim>* blockF,
                       InterpBlockComm<FloatType, Dim>* comm)
      : BlockF(blockF), Comm(comm) {}

  const std::vector<InterpSource<Dim>>& getSends() const { return Comm->SendCells; }
  const InterpSource<Dim>& getSend(std::size_t i) const { return (Comm->SendCells)[i]; }

  const std::vector<std::size_t>& getRecvs() const { return Comm->RecvCells; }
  const std::size_t& getRecv(std::size_t i) const { return (Comm->RecvCells)[i]; }

  const std::vector<InterpWeight<FloatType, Dim>>& getWeights() const {
    return Comm->InterpWeights;
  }
  const InterpWeight<FloatType, Dim>& getWeight(std::size_t i) const {
    return (Comm->InterpWeights)[i];
  }
};

// BlockField with communication structure(in Block2D)
template <typename FieldType, typename FloatType, unsigned int Dim>
class BlockField {
 private:
  // block field
  FieldType _Field;
  // block(geometry) structure of the field
  BasicBlock<FloatType, Dim>& _Block;

  std::vector<BlockFieldComm<FieldType, FloatType, Dim>> Comms;
  // average block comm, get from higher level block
  std::vector<InterpBlockFieldComm<FieldType, FloatType, Dim>> AverComms;
  // interp block comm, get from lower level block
  std::vector<InterpBlockFieldComm<FieldType, FloatType, Dim>> InterpComms;

 public:
  BlockField() = default;
  BlockField(BasicBlock<FloatType, Dim>& block) : _Block(block), _Field(block.getN()) {}
  template <typename datatype>
  BlockField(BasicBlock<FloatType, Dim>& block, datatype initvalue)
      : _Block(block), _Field(block.getN(), initvalue) {}
  ~BlockField() = default;

  FieldType& getField() { return _Field; }
  const FieldType& getField() const { return _Field; }

  BasicBlock<FloatType, Dim>& getBlock() { return _Block; }
  const BasicBlock<FloatType, Dim>& getBlock() const { return _Block; }

  std::vector<BlockFieldComm<FieldType, FloatType, Dim>>& getComms() { return Comms; }
  std::vector<InterpBlockFieldComm<FieldType, FloatType, Dim>>& getAverComms() {
    return AverComms;
  }
  std::vector<InterpBlockFieldComm<FieldType, FloatType, Dim>>& getInterpComms() {
    return InterpComms;
  }

  void normalcommunicate() {
    for (BlockFieldComm<FieldType, FloatType, Dim>& comm : Comms) {
      BlockField<FieldType, FloatType, Dim>* nblockF = comm.BlockF;
      std::size_t size = comm.getRecvs().size();
      for (int iArr = 0; iArr < _Field.Size(); ++iArr) {
        const auto& nArray = nblockF->getField().getField(iArr);
        auto& Array = _Field.getField(iArr);
        for (std::size_t id = 0; id < size; ++id) {
          Array.set(comm.getRecv(id), nArray[comm.getSend(id)]);
        }
      }
    }
  }
  void avercommunicate() {
    for (InterpBlockFieldComm<FieldType, FloatType, Dim>& comm : AverComms) {
      BlockField<FieldType, FloatType, Dim>* nblockF = comm.BlockF;
      std::size_t size = comm.getRecvs().size();
      for (int iArr = 0; iArr < _Field.Size(); ++iArr) {
        const auto& nArray = nblockF->getField().getField(iArr);
        auto& Array = _Field.getField(iArr);
        for (std::size_t id = 0; id < size; ++id) {
          const InterpSource<Dim>& sends = comm.getSend(id);
          auto sum = nArray[sends[0]];
          for (int i = 1; i < comm.Comm->getSourceNum(); ++i) {
            sum += nArray[sends[i]];
          }
          Array.set(comm.getRecv(id), sum * comm.Comm->getUniformWeight());
        }
      }
    }
  }
  void interpcommunicate() {
    for (InterpBlockFieldComm<FieldType, FloatType, Dim>& comm : InterpComms) {
      BlockField<FieldType, FloatType, Dim>* nblockF = comm.BlockF;
      std::size_t size = comm.getRecvs().size();
      for (int iArr = 0; iArr < _Field.Size(); ++iArr) {
        const auto& nArray = nblockF->getField().getField(iArr);
        auto& Array = _Field.getField(iArr);
        for (std::size_t id = 0; id < size; ++id) {
          const InterpSource<Dim>& sends = comm.getSend(id);
          const InterpWeight<FloatType, Dim>& weights = comm.getWeight(id);
          auto value = nArray[sends[0]] * weights[0];
          for (int i = 1; i < comm.Comm->getSourceNum(); ++i) {
            value += nArray[sends[i]] * weights[i];
          }
          Array.set(comm.getRecv(id), value);
        }
      }
    }
  }
};

// BlockField manager with communication enabled
template <typename FieldType, typename FloatType, unsigned int Dim>
class BlockFieldManager {
 private:
  // block fields
  std::vector<BlockField<FieldType, FloatType, Dim>> _Fields;
  // blockstructure
  // obj of _BlockGeo should not be destroyed
  BlockGeometry<FloatType, Dim>& _BlockGeo;

 public:
  BlockFieldManager(BlockGeometry<FloatType, Dim>& blockgeometry)
      : _BlockGeo(blockgeometry) {
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      _Fields.emplace_back(block);
    }
    InitComm();
  }
  template <typename datatype>
  BlockFieldManager(BlockGeometry<FloatType, Dim>& blockgeometry, datatype initvalue)
      : _BlockGeo(blockgeometry) {
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      _Fields.emplace_back(block, initvalue);
    }
    InitComm();
  }
  ~BlockFieldManager() = default;


  void InitAndComm() {
    Init();
    CommunicateAll();
  }
  template <typename datatype>
  void InitAndComm(datatype initvalue) {
    Init<datatype>(initvalue);
    CommunicateAll();
  }
  void InitAndComm(BlockGeometryHelper<FloatType, Dim>& GeoHelper) {
    Init(GeoHelper);
    CommunicateAll();
  }
  template <typename datatype>
  void InitAndComm(BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue) {
    Init(GeoHelper, initvalue);
    CommunicateAll();
  }
  // pure init
  void Init() {
    _Fields.clear();
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      _Fields.emplace_back(block);
    }
    InitComm();
  }
  // pure init with initvalue
  template <typename datatype>
  void Init(datatype initvalue) {
    _Fields.clear();
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      _Fields.emplace_back(block, initvalue);
    }
    InitComm();
  }
  // init with initvalue, this assumes that the BlockGeo is already initialized
  void Init(BlockGeometryHelper<FloatType, Dim>& GeoHelper) {
    std::vector<BlockField<FieldType, FloatType, Dim>> NewFields;
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      NewFields.emplace_back(block);
    }
    // data transfer
    FieldDataTransfer(GeoHelper, NewFields);
    _Fields.swap(NewFields);
    InitComm();
  }
  // init with initvalue and data transfer
  template <typename datatype>
  void Init(BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue) {
    std::vector<BlockField<FieldType, FloatType, Dim>> NewFields;
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      NewFields.emplace_back(block, initvalue);
    }
    // data transfer
    FieldDataTransfer(GeoHelper, NewFields);
    _Fields.swap(NewFields);
    InitComm();
  }
  // init with initvalue and data transfer
  template <typename datatype, typename flagtype, typename Func>
  void InitCopy(
    BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue,
    const BlockFieldManager<ScalerField<flagtype>, FloatType, Dim>& FlagFManager,
    std::uint8_t flag, Func func) {
    std::vector<BlockField<FieldType, FloatType, Dim>> NewFields;
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      NewFields.emplace_back(block);
    }
    // func
    // forEach(FlagFManager, flag, func); but this is operate on NewFields
    int iblock = 0;
    for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
      FieldType& field = NewFields[iblock].getField();
      const GenericArray<flagtype>& flagarr =
        FlagFManager.getBlockField(iblock).getField().getField(0);
      blockgeo.forEach(flagarr, flag,
                       [&field, &func](std::size_t id) { func(field, id); });
      ++iblock;
    }
    // data transfer
    FieldDataCopyTransfer(GeoHelper, NewFields);
    _Fields.swap(NewFields);
    InitComm();
  }

  BlockGeometry<FloatType, Dim>& getGeo() { return _BlockGeo; }
  const BlockGeometry<FloatType, Dim>& getGeo() const { return _BlockGeo; }

  BlockField<FieldType, FloatType, Dim>& getBlockField(int i) { return _Fields[i]; }
  const BlockField<FieldType, FloatType, Dim>& getBlockField(int i) const {
    return _Fields[i];
  }

  std::vector<BlockField<FieldType, FloatType, Dim>>& getBlockFields() { return _Fields; }
  const std::vector<BlockField<FieldType, FloatType, Dim>>& getBlockFields() const {
    return _Fields;
  }

  void InitComm() {
    int blockid = 0;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      Block<FloatType, Dim>& block = _BlockGeo.getBlock(blockid);
      ++blockid;
      // normal communication
      std::vector<BlockFieldComm<FieldType, FloatType, Dim>>& Fcomms = blockF.getComms();
      Fcomms.clear();
      for (BlockComm<FloatType, Dim>& comm : block.getCommunicators()) {
        Fcomms.emplace_back(&_Fields[comm.getSendId()], &comm);
      }
      // average communication
      std::vector<InterpBlockFieldComm<FieldType, FloatType, Dim>>& Avercomms =
        blockF.getAverComms();
      Avercomms.clear();
      for (InterpBlockComm<FloatType, Dim>& comm : block.getAverageBlockComm()) {
        Avercomms.emplace_back(&_Fields[comm.getSendId()], &comm);
      }
      // interp communication
      std::vector<InterpBlockFieldComm<FieldType, FloatType, Dim>>& Interpcomms =
        blockF.getInterpComms();
      Interpcomms.clear();
      for (InterpBlockComm<FloatType, Dim>& comm : block.getInterpBlockComm()) {
        Interpcomms.emplace_back(&_Fields[comm.getSendId()], &comm);
      }
    }
  }

  // construct new field with Geohelper and copy data from old field, then swap
  void FieldDataTransfer(BlockGeometryHelper<FloatType, Dim>& GeoHelper,
                         std::vector<BlockField<FieldType, FloatType, Dim>>& NewFields) {
    // copy from old field to new field
#pragma omp parallel for num_threads(Thread_Num)
    for (int inewblock = 0; inewblock < NewFields.size(); ++inewblock) {
      BlockField<FieldType, FloatType, Dim>& NewField = NewFields[inewblock];
      // newbaseblock could be accessed either from GeoHelper or _BlockGeo
      const BasicBlock<FloatType, Dim>& newbaseblock = GeoHelper.getBasicBlock(inewblock);
      const BasicBlock<FloatType, Dim>& newblock = NewField.getBlock();
      std::uint8_t Level = newblock.getLevel();
      // find overlapped old block field
      for (int iblock = 0; iblock < GeoHelper.getOldBasicBlocks().size(); ++iblock) {
        // oldbaseblock could be accessed only from GeoHelper
        const BasicBlock<FloatType, Dim>& baseblock = GeoHelper.getOldBasicBlock(iblock);
        if (isOverlapped(newbaseblock, baseblock)) {
          const BlockField<FieldType, FloatType, Dim>& Field = _Fields[iblock];
          int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
          const BasicBlock<FloatType, Dim>& block = baseblock.getExtBlock(overlap);
          if (Level == block.getLevel()) {
            // copy
            FieldCopy2D(Field.getField(), NewField.getField(), block, baseblock, newblock,
                        newbaseblock);
          } else if (Level > block.getLevel()) {
            // interp
            FieldInterpolation2D(Field.getField(), NewField.getField(), block, baseblock,
                                 newblock, newbaseblock);
          } else if (Level < block.getLevel()) {
            // average
            FieldAverage2D(Field.getField(), NewField.getField(), block, baseblock,
                           newblock, newbaseblock);
          }
        }
      }
    }
  }

  // construct new field with Geohelper and transfer with average, then swap
  void FieldDataCopyTransfer(
    BlockGeometryHelper<FloatType, Dim>& GeoHelper,
    std::vector<BlockField<FieldType, FloatType, Dim>>& NewFields) {
    // copy from old field to new field
#pragma omp parallel for num_threads(Thread_Num)
    for (int inewblock = 0; inewblock < NewFields.size(); ++inewblock) {
      BlockField<FieldType, FloatType, Dim>& NewField = NewFields[inewblock];
      // newbaseblock could be accessed either from GeoHelper or _BlockGeo
      const BasicBlock<FloatType, Dim>& newbaseblock = GeoHelper.getBasicBlock(inewblock);
      const BasicBlock<FloatType, Dim>& newblock = NewField.getBlock();
      std::uint8_t Level = newblock.getLevel();
      // find overlapped old block field
      for (int iblock = 0; iblock < GeoHelper.getOldBasicBlocks().size(); ++iblock) {
        // oldbaseblock could be accessed only from GeoHelper
        const BasicBlock<FloatType, Dim>& baseblock = GeoHelper.getOldBasicBlock(iblock);
        if (isOverlapped(newbaseblock, baseblock)) {
          const BlockField<FieldType, FloatType, Dim>& Field = _Fields[iblock];
          int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
          const BasicBlock<FloatType, Dim>& block = baseblock.getExtBlock(overlap);
          if (Level == block.getLevel()) {
            // copy
            FieldCopy2D(Field.getField(), NewField.getField(), block, baseblock, newblock,
                        newbaseblock);
          }
        }
      }
    }
  }

  template <typename LatSet, typename datatype>
  void SetupBoundary(const AABB<FloatType, 2>& block, datatype bdvalue) {
    int iblock = 0;
    for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
      FieldType& field = _Fields[iblock].getField();
      blockgeo.template SetupBoundary<FieldType, datatype, LatSet>(block, field, bdvalue);
      ++iblock;
    }
    NormalCommunicate();
  }

  // call forEach(AABBs, [&](FieldType& field, std::size_t id){});
  template <typename Func>
  void forEach(const AABB<FloatType, Dim>& AABBs, Func func) {
    int iblock = 0;
    for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
      FieldType& field = _Fields[iblock].getField();
      blockgeo.forEach(AABBs, [&field, &func](std::size_t id) { func(field, id); });
      ++iblock;
    }
  }

  // call forEach(AABBs, FlagFManager, flag, [&](FieldType& field, std::size_t id){});
  template <typename flagtype, typename Func>
  void forEach(
    const AABB<FloatType, Dim>& AABBs,
    const BlockFieldManager<ScalerField<flagtype>, FloatType, Dim>& FlagFManager,
    std::uint8_t flag, Func func) {
    int iblock = 0;
    for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
      FieldType& field = _Fields[iblock].getField();
      const GenericArray<flagtype>& flagarr =
        FlagFManager.getBlockField(iblock).getField().getField(0);
      blockgeo.forEach(AABBs, flagarr, flag,
                       [&field, &func](std::size_t id) { func(field, id); });
      ++iblock;
    }
  }

  // call forEach(FlagFManager, flag, [&](FieldType& field, std::size_t id){});
  template <typename flagtype, typename Func>
  void forEach(
    const BlockFieldManager<ScalerField<flagtype>, FloatType, Dim>& FlagFManager,
    std::uint8_t flag, Func func) {
    int iblock = 0;
    for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
      FieldType& field = _Fields[iblock].getField();
      const GenericArray<flagtype>& flagarr =
        FlagFManager.getBlockField(iblock).getField().getField(0);
      blockgeo.forEach(flagarr, flag,
                       [&field, &func](std::size_t id) { func(field, id); });
      ++iblock;
    }
  }

  // communication
  void NormalCommunicate(std::int64_t count) {
#pragma omp parallel for num_threads(Thread_Num)
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if (count % (static_cast<int>(pow(2, deLevel))) == 0) {
        blockF.normalcommunicate();
      }
    }
  }

  // average communicate, do not use this function for pop field communication
  void AverCommunicate(std::int64_t count) {
#pragma omp parallel for num_threads(Thread_Num)
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if (count % (static_cast<int>(pow(2, deLevel))) == 0) {
        blockF.avercommunicate();
      }
    }
  }

  // interp communicate, do not use this function for pop field communication
  void InterpCommunicate(std::int64_t count) {
#pragma omp parallel for num_threads(Thread_Num)
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if (count % (static_cast<int>(pow(2, deLevel))) == 0) {
        blockF.interpcommunicate();
      }
    }
  }

  void CommunicateAll(std::int64_t count) {
    NormalCommunicate(count);
    AverCommunicate(count);
    InterpCommunicate(count);
  }

  void NormalCommunicate() {
#pragma omp parallel for num_threads(Thread_Num)
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      blockF.normalcommunicate();
    }
  }
  void AverCommunicate() {
#pragma omp parallel for num_threads(Thread_Num)
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      blockF.avercommunicate();
    }
  }
  void InterpCommunicate() {
#pragma omp parallel for num_threads(Thread_Num)
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      blockF.interpcommunicate();
    }
  }
  void CommunicateAll() {
    NormalCommunicate();
    AverCommunicate();
    InterpCommunicate();
  }
};
