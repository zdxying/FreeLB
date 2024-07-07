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
#include "parallel/communicator.h"
#include "utils/alias.h"

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

  // const std::vector<InterpWeight<FloatType, Dim>>& getWeights() const {
  //   return Comm->InterpWeights;
  // }
};

// BlockField with communication structure(in Block2D)
template <typename FieldType, typename FloatType, unsigned int Dim>
class BlockField : public FieldType {
 public:
  using datatype = typename FieldType::value_type;
  static constexpr unsigned int array_dim = FieldType::array_dim;
  static constexpr bool isField = FieldType::isField;

 private:
  // block(geometry) structure of the field
  Block<FloatType, Dim>& _Block;

  std::vector<BlockFieldComm<FieldType, FloatType, Dim>> Comms;
  // average block comm, get from higher level block
  std::vector<InterpBlockFieldComm<FieldType, FloatType, Dim>> AverComms;
  // interp block comm, get from lower level block
  std::vector<InterpBlockFieldComm<FieldType, FloatType, Dim>> InterpComms;

#ifdef MPI_ENABLED
  // MPI buffer
  MPIBlockBuffer<datatype> MPIBuffer;
  MPIBlockBuffer<datatype> MPIAverBuffer;
  MPIBlockBuffer<datatype> MPIInterpBuffer;
#endif

 public:
  BlockField() = default;
  BlockField(Block<FloatType, Dim>& block) : _Block(block), FieldType(block.getN()) {
    MPIBufferInit();
  }
  BlockField(Block<FloatType, Dim>& block, datatype initvalue)
      : _Block(block), FieldType(block.getN(), initvalue) {
    MPIBufferInit();
  }
  // copy constructor
  BlockField(const BlockField& blockF) : _Block(blockF._Block), FieldType(blockF) {
    MPIBufferInit();
  }
  // move constructor
  BlockField(BlockField&& blockF) noexcept
      : _Block(blockF._Block), FieldType(std::move(blockF)) {
    MPIBufferInit();
  }
  // copy assignment operator
  BlockField& operator=(const BlockField& blockF) {
    if (this != &blockF) {
      _Block = blockF._Block;
      FieldType::operator=(blockF);
      MPIBufferInit();
    }
    return *this;
  }
  // move assignment operator
  BlockField& operator=(BlockField&& blockF) noexcept {
    if (this != &blockF) {
      _Block = blockF._Block;
      FieldType::operator=(std::move(blockF));
      MPIBufferInit();
    }
    return *this;
  }


  ~BlockField() = default;

  void MPIBufferInit() {
#ifdef MPI_ENABLED
    MPIBlockBufferInit(_Block.getMPIBlockComm(), MPIBuffer, array_dim);
    MPIBlockBufferInit(_Block.getMPIAverBlockComm(), MPIAverBuffer, array_dim);
    MPIBlockBufferInit(_Block.getMPIInterpBlockComm(), MPIInterpBuffer, array_dim);
#endif
  }

  Block<FloatType, Dim>& getBlock() { return _Block; }
  const Block<FloatType, Dim>& getBlock() const { return _Block; }

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
      for (int iArr = 0; iArr < nblockF->Size(); ++iArr) {
        const auto& nArray = nblockF->getField(iArr);
        auto& Array = this->getField(iArr);
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
      for (int iArr = 0; iArr < FieldType::Size(); ++iArr) {
        const auto& nArray = nblockF->getField(iArr);
        auto& Array = FieldType::getField(iArr);
        for (std::size_t id = 0; id < size; ++id) {
          const InterpSource<Dim>& sends = comm.getSend(id);
          auto sum = nArray[sends[0]];
          if constexpr (Dim == 2) {
            sum += nArray[sends[1]] + nArray[sends[2]] + nArray[sends[3]];
          } else if constexpr (Dim == 3) {
            sum += nArray[sends[1]] + nArray[sends[2]] + nArray[sends[3]] +
                   nArray[sends[4]] + nArray[sends[5]] + nArray[sends[6]] +
                   nArray[sends[7]];
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
      for (int iArr = 0; iArr < FieldType::Size(); ++iArr) {
        const auto& nArray = nblockF->getField(iArr);
        auto& Array = FieldType::getField(iArr);
        const auto& sends = comm.getSends();
        const auto& recvs = comm.getRecvs();
        std::size_t idx = 0;
        for (std::size_t id = 0; id < size;) {
          Interpolation<FloatType, Dim>(Array, nArray, sends, recvs, id, idx);
            // {
            //   const InterpSource<Dim>& sends = comm.getSend(id);
            //   auto value = nArray[sends[0]] * comm.Comm->getIntpWeight()[0][0] +
            //                nArray[sends[1]] * comm.Comm->getIntpWeight()[0][1] +
            //                nArray[sends[2]] * comm.Comm->getIntpWeight()[0][2] +
            //                nArray[sends[3]] * comm.Comm->getIntpWeight()[0][3];
            //   Array.set(comm.getRecv(id), value);
            // }
            // {
            //   const InterpSource<Dim>& sends = comm.getSend(id + 1);
            //   auto value = nArray[sends[0]] * comm.Comm->getIntpWeight()[1][0] +
            //                nArray[sends[1]] * comm.Comm->getIntpWeight()[1][1] +
            //                nArray[sends[2]] * comm.Comm->getIntpWeight()[1][2] +
            //                nArray[sends[3]] * comm.Comm->getIntpWeight()[1][3];
            //   Array.set(comm.getRecv(id + 1), value);
            // }
            // {
            //   const InterpSource<Dim>& sends = comm.getSend(id + 2);
            //   auto value = nArray[sends[0]] * comm.Comm->getIntpWeight()[2][0] +
            //                nArray[sends[1]] * comm.Comm->getIntpWeight()[2][1] +
            //                nArray[sends[2]] * comm.Comm->getIntpWeight()[2][2] +
            //                nArray[sends[3]] * comm.Comm->getIntpWeight()[2][3];
            //   Array.set(comm.getRecv(id + 2), value);
            // }
            // {
            //   const InterpSource<Dim>& sends = comm.getSend(id + 3);
            //   auto value = nArray[sends[0]] * comm.Comm->getIntpWeight()[3][0] +
            //                nArray[sends[1]] * comm.Comm->getIntpWeight()[3][1] +
            //                nArray[sends[2]] * comm.Comm->getIntpWeight()[3][2] +
            //                nArray[sends[3]] * comm.Comm->getIntpWeight()[3][3];
            //   Array.set(comm.getRecv(id + 3), value);
            // }
        }
      }
    }
  }
#ifdef MPI_ENABLED
  MPIBlockBuffer<datatype>& getMPIBuffer() { return MPIBuffer; }
  MPIBlockBuffer<datatype>& getMPIAverBuffer() { return MPIAverBuffer; }
  MPIBlockBuffer<datatype>& getMPIInterpBuffer() { return MPIInterpBuffer; }
  // send data for normal communication
  void mpiNormalSend(std::vector<MPI_Request>& SendRequests) {
    const MPIBlockComm& MPIComm = _Block.getMPIBlockComm();
    // add to send buffer
    for (int i = 0; i < MPIComm.Senders.size(); ++i) {
      std::vector<datatype>& buffer = MPIBuffer.SendBuffers[i];
      const std::vector<std::size_t>& sends = MPIComm.Senders[i].SendCells;
      std::size_t bufidx = 0;
      for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
        const auto& Array = FieldType::getField(iArr);
        for (std::size_t id : sends) {
          buffer[bufidx] = Array[id];
          ++bufidx;
        }
      }
    }
    // non-blocking send
    for (int i = 0; i < MPIComm.Senders.size(); ++i) {
      MPI_Request request;
      std::vector<datatype>& buffer = MPIBuffer.SendBuffers[i];
      mpi().iSend(buffer.data(), buffer.size(), MPIComm.Senders[i].RecvRank, &request,
                  MPIComm.Senders[i].RecvBlockid);
      SendRequests.push_back(request);
    }
  }
  // recv data for normal communication
  void mpiNormalRecv(std::vector<MPI_Request>& RecvRequests) {
    const MPIBlockComm& MPIComm = _Block.getMPIBlockComm();
    // non-blocking recv
    for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
      MPI_Request request;
      std::vector<datatype>& buffer = MPIBuffer.RecvBuffers[i];
      mpi().iRecv(buffer.data(), buffer.size(), MPIComm.Recvers[i].SendRank, &request,
                  _Block.getBlockId());
      RecvRequests.push_back(request);
    }
  }
  // set data for normal communication
  void mpiNormalSet(int& reqidx, std::vector<MPI_Request>& RecvRequests) {
    const MPIBlockComm& MPIComm = _Block.getMPIBlockComm();
    // wait and set field data
    for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
      MPI_Wait(&RecvRequests[i + reqidx], MPI_STATUS_IGNORE);
      const std::vector<datatype>& buffer = MPIBuffer.RecvBuffers[i];
      const std::vector<std::size_t>& recvs = MPIComm.Recvers[i].RecvCells;
      std::size_t bufidx = 0;
      for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
        auto& Array = FieldType::getField(iArr);
        for (std::size_t id : recvs) {
          Array[id] = buffer[bufidx];
          ++bufidx;
        }
      }
    }
    reqidx += MPIComm.Recvers.size();
  }

  // send data for average communication
  void mpiAverSend(std::vector<MPI_Request>& SendRequests) {
    const MPIInterpBlockComm<FloatType, Dim>& MPIComm = _Block.getMPIAverBlockComm();
    // add to send buffer
    for (int i = 0; i < MPIComm.Senders.size(); ++i) {
      std::vector<datatype>& buffer = MPIAverBuffer.SendBuffers[i];
      const std::vector<InterpSource<Dim>>& sendcells = MPIComm.Senders[i].SendCells;
      constexpr FloatType weight = InterpBlockComm<FloatType, Dim>::getUniformWeight();
      std::size_t bufidx = 0;
      for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
        const auto& Array = FieldType::getField(iArr);
        for (const InterpSource<Dim>& sends : sendcells) {
          buffer[bufidx] = getAverage<FloatType, Dim>(Array, sends);
          ++bufidx;
        }
      }
    }
    // non-blocking send
    for (int i = 0; i < MPIComm.Senders.size(); ++i) {
      MPI_Request request;
      std::vector<datatype>& buffer = MPIAverBuffer.SendBuffers[i];
      mpi().iSend(buffer.data(), buffer.size(), MPIComm.Senders[i].RecvRank, &request,
                  MPIComm.Senders[i].RecvBlockid);
      SendRequests.push_back(request);
    }
  }
  // recv data for average communication
  void mpiAverRecv(std::vector<MPI_Request>& RecvRequests) {
    const MPIInterpBlockComm<FloatType, Dim>& MPIComm = _Block.getMPIAverBlockComm();
    // non-blocking recv
    for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
      MPI_Request request;
      std::vector<datatype>& buffer = MPIAverBuffer.RecvBuffers[i];
      mpi().iRecv(buffer.data(), buffer.size(), MPIComm.Recvers[i].SendRank, &request,
                  _Block.getBlockId());
      RecvRequests.push_back(request);
    }
  }
  // set data for average communication
  void mpiAverSet(int& reqidx, std::vector<MPI_Request>& RecvRequests) {
    const MPIInterpBlockComm<FloatType, Dim>& MPIComm = _Block.getMPIAverBlockComm();
    // wait and set field data
    for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
      MPI_Wait(&RecvRequests[i + reqidx], MPI_STATUS_IGNORE);
      const std::vector<datatype>& buffer = MPIAverBuffer.RecvBuffers[i];
      const std::vector<std::size_t>& recvcells = MPIComm.Recvers[i].RecvCells;
      std::size_t bufidx = 0;
      for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
        auto& Array = FieldType::getField(iArr);
        for (std::size_t id : recvcells) {
          Array[id] = buffer[bufidx];
          ++bufidx;
        }
      }
    }
    reqidx += MPIComm.Recvers.size();
  }

  // send data for interp communication
  void mpiIntpSend(std::vector<MPI_Request>& SendRequests) {
    const MPIInterpBlockComm<FloatType, Dim>& MPIComm = _Block.getMPIInterpBlockComm();
    // add to send buffer
    for (int i = 0; i < MPIComm.Senders.size(); ++i) {
      std::vector<datatype>& buffer = MPIInterpBuffer.SendBuffers[i];
      const std::vector<InterpSource<Dim>>& sendcells = MPIComm.Senders[i].SendCells;
      const std::size_t size = sendcells.size();
      std::size_t bufidx = 0;
      for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
        const auto& Array = FieldType::getField(iArr);
        for (std::size_t i = 0; i < size;) {
          Interpolation<FloatType, Dim>(Array, sendcells, i, buffer, bufidx);
        }
      }
    }
    // non-blocking send
    for (int i = 0; i < MPIComm.Senders.size(); ++i) {
      MPI_Request request;
      std::vector<datatype>& buffer = MPIInterpBuffer.SendBuffers[i];
      mpi().iSend(buffer.data(), buffer.size(), MPIComm.Senders[i].RecvRank, &request,
                  MPIComm.Senders[i].RecvBlockid);
      SendRequests.push_back(request);
    }
  }
  // recv data for interp communication
  void mpiIntpRecv(std::vector<MPI_Request>& RecvRequests) {
    const MPIInterpBlockComm<FloatType, Dim>& MPIComm = _Block.getMPIInterpBlockComm();
    // non-blocking recv
    for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
      MPI_Request request;
      std::vector<datatype>& buffer = MPIInterpBuffer.RecvBuffers[i];
      mpi().iRecv(buffer.data(), buffer.size(), MPIComm.Recvers[i].SendRank, &request,
                  _Block.getBlockId());
      RecvRequests.push_back(request);
    }
  }
  // set data for interp communication
  void mpiIntpSet(int& reqidx, std::vector<MPI_Request>& RecvRequests) {
    const MPIInterpBlockComm<FloatType, Dim>& MPIComm = _Block.getMPIInterpBlockComm();
    // wait and set field data
    for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
      MPI_Wait(&RecvRequests[i + reqidx], MPI_STATUS_IGNORE);
      const std::vector<datatype>& buffer = MPIInterpBuffer.RecvBuffers[i];
      const std::vector<std::size_t>& recvcells = MPIComm.Recvers[i].RecvCells;
      std::size_t bufidx = 0;
      for (unsigned int iArr = 0; iArr < array_dim; ++iArr) {
        auto& Array = FieldType::getField(iArr);
        for (std::size_t id : recvcells) {
          Array[id] = buffer[bufidx];
          ++bufidx;
        }
      }
    }
    reqidx += MPIComm.Recvers.size();
  }

#endif
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
  using datatype = typename FieldType::value_type;
  using array_type = typename FieldType::array_type;
  using field_type = FieldType;
  using float_type = FloatType;
  static constexpr unsigned int dim = Dim;
  static constexpr bool isField = FieldType::isField;

  BlockFieldManager(BlockGeometry<FloatType, Dim>& blockgeometry)
      : _BlockGeo(blockgeometry) {
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      _Fields.emplace_back(block);
    }
    InitComm();
  }
  BlockFieldManager(BlockGeometry<FloatType, Dim>& blockgeometry, datatype initvalue)
      : _BlockGeo(blockgeometry) {
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      _Fields.emplace_back(block, initvalue);
    }
    InitComm();
  }

  // copy constructor
  BlockFieldManager(const BlockFieldManager& blockFManager)
      : _Fields(blockFManager._Fields), _BlockGeo(blockFManager._BlockGeo) {}
  // move constructor
  BlockFieldManager(BlockFieldManager&& blockFManager) noexcept
      : _Fields(std::move(blockFManager._Fields)), _BlockGeo(blockFManager._BlockGeo) {}

  // copy assignment operator
  BlockFieldManager& operator=(const BlockFieldManager& blockFManager) {
    if (this != &blockFManager) {
      _Fields = blockFManager._Fields;
      _BlockGeo = blockFManager._BlockGeo;
    }
    return *this;
  }
  // move assignment operator
  BlockFieldManager& operator=(BlockFieldManager&& blockFManager) noexcept {
    if (this != &blockFManager) {
      _Fields = std::move(blockFManager._Fields);
      _BlockGeo = blockFManager._BlockGeo;
    }
    return *this;
  }


  ~BlockFieldManager() = default;


  void InitAndComm() {
    Init();
    CommunicateAll();
  }
  void InitAndComm(datatype initvalue) {
    Init(initvalue);
    CommunicateAll();
  }
  void InitAndComm(BlockGeometryHelper<FloatType, Dim>& GeoHelper) {
    Init(GeoHelper);
    CommunicateAll();
  }
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
  void Init(datatype initvalue) {
    _Fields.clear();
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      _Fields.emplace_back(block, initvalue);
    }
    InitComm();
  }
  void NonFieldInit(datatype initvalue) {
    _Fields.clear();
    for (Block<FloatType, Dim>& block : _BlockGeo.getBlocks()) {
      _Fields.emplace_back(block, initvalue);
    }
  }
  // this assumes that the BlockGeo is already initialized
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
  template <typename FlagFieldType, typename Func>
  void InitCopy(BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue,
                const BlockFieldManager<FlagFieldType, FloatType, Dim>& FlagFManager,
                std::uint8_t flag, Func func) {
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

  // get block field with block id
  BlockField<FieldType, FloatType, Dim>* findBlockField(int blockid) {
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      if (blockF.getBlock().getBlockId() == blockid) return &blockF;
    }
    std::cerr << "[BlockFieldManager]: can't find blockF with blockid " << blockid
              << std::endl;
    exit(1);
  }

  void InitComm() {
    if constexpr (!FieldType::isField) {
      return;
    }
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      Block<FloatType, Dim>& block = blockF.getBlock();
      // normal communication
      std::vector<BlockFieldComm<FieldType, FloatType, Dim>>& Fcomms = blockF.getComms();
      Fcomms.clear();
      for (BlockComm<FloatType, Dim>& comm : block.getCommunicators()) {
        Fcomms.emplace_back(findBlockField(comm.getSendId()), &comm);
      }
      // average communication
      std::vector<InterpBlockFieldComm<FieldType, FloatType, Dim>>& Avercomms =
        blockF.getAverComms();
      Avercomms.clear();
      for (InterpBlockComm<FloatType, Dim>& comm : block.getAverageBlockComm()) {
        Avercomms.emplace_back(findBlockField(comm.getSendId()), &comm);
      }
      // interp communication
      std::vector<InterpBlockFieldComm<FieldType, FloatType, Dim>>& Interpcomms =
        blockF.getInterpComms();
      Interpcomms.clear();
      for (InterpBlockComm<FloatType, Dim>& comm : block.getInterpBlockComm()) {
        Interpcomms.emplace_back(findBlockField(comm.getSendId()), &comm);
      }
    }
  }

  // construct new field with Geohelper and copy data from old field, then swap
  void FieldDataTransfer(BlockGeometryHelper<FloatType, Dim>& GeoHelper,
                         std::vector<BlockField<FieldType, FloatType, Dim>>& NewFields) {
    if constexpr (!FieldType::isField) {
      return;
    }
    // copy from old field to new field
#pragma omp parallel for num_threads(Thread_Num)
    for (int inewblock = 0; inewblock < NewFields.size(); ++inewblock) {
      BlockField<FieldType, FloatType, Dim>& NewField = NewFields[inewblock];
      // newbaseblock could be accessed either from GeoHelper or _BlockGeo
      const BasicBlock<FloatType, Dim>& newbaseblock =
        GeoHelper.getAllBasicBlock(inewblock);
      const BasicBlock<FloatType, Dim>& newblock = NewField.getBlock();
      std::uint8_t Level = newblock.getLevel();
      // find overlapped old block field
      for (int iblock = 0; iblock < GeoHelper.getAllOldBasicBlocks().size(); ++iblock) {
        // oldbaseblock could be accessed only from GeoHelper
        const BasicBlock<FloatType, Dim>& baseblock =
          GeoHelper.getAllOldBasicBlock(iblock);
        if (isOverlapped(newbaseblock, baseblock)) {
          const BlockField<FieldType, FloatType, Dim>& Field = _Fields[iblock];
          int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
          const BasicBlock<FloatType, Dim>& block = baseblock.getExtBlock(overlap);
          if (Level == block.getLevel()) {
            // copy
            FieldCopy2D(Field, NewField, block, baseblock, newblock, newbaseblock);
          } else if (Level > block.getLevel()) {
            // interp
            FieldInterpolation2D(Field, NewField, block, baseblock, newblock,
                                 newbaseblock);
          } else if (Level < block.getLevel()) {
            // average
            FieldAverage2D(Field, NewField, block, baseblock, newblock, newbaseblock);
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
      const BasicBlock<FloatType, Dim>& newbaseblock =
        GeoHelper.getAllBasicBlock(inewblock);
      const BasicBlock<FloatType, Dim>& newblock = NewField.getBlock();
      std::uint8_t Level = newblock.getLevel();
      // find overlapped old block field
      for (int iblock = 0; iblock < GeoHelper.getAllOldBasicBlocks().size(); ++iblock) {
        // oldbaseblock could be accessed only from GeoHelper
        const BasicBlock<FloatType, Dim>& baseblock =
          GeoHelper.getAllOldBasicBlock(iblock);
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

  template <typename LatSet>
  void SetupBoundary(const AABB<FloatType, Dim>& block, datatype bdvalue) {
    int iblock = 0;
    for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
      auto& field = _Fields[iblock];
      blockgeo.template SetupBoundary<FieldType, LatSet>(block, field, bdvalue);
      ++iblock;
    }
    NormalCommunicate();
#ifdef MPI_ENABLED
    MPINormalCommunicate();
#endif
  }

  template <typename Func>
  void forEachField(Func func) {
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      func(blockF);
    }
  }

  void InitValue(datatype initvalue) {
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      blockF.Init(initvalue);
    }
  }

  template <typename Func>
  void forEach(Func func) {
    int iblock = 0;
    for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
      BlockField<FieldType, FloatType, Dim>& blockfield = _Fields[iblock];
      blockgeo.forEach([&blockfield, &func](std::size_t id) { func(blockfield, id); });
      ++iblock;
    }
  }

  // call forEach(AABBs, [&](FieldType& field, std::size_t id){});
  template <typename Func>
  void forEach(const AABB<FloatType, Dim>& AABBs, Func func) {
    int iblock = 0;
    for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
      auto& field = _Fields[iblock];
      blockgeo.forEach(AABBs, [&field, &func](std::size_t id) { func(field, id); });
      ++iblock;
    }
  }

  // call forEach(AABBs, FlagFManager, flag, [&](FieldType& field, std::size_t id){});
  template <typename FlagFieldType, typename Func>
  void forEach(const AABB<FloatType, Dim>& AABBs,
               const BlockFieldManager<FlagFieldType, FloatType, Dim>& FlagFManager,
               std::uint8_t flag, Func func) {
    int iblock = 0;
    for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
      auto& field = _Fields[iblock];
      const auto& flagarr = FlagFManager.getBlockField(iblock).getField(0);
      blockgeo.forEach(AABBs, flagarr, flag,
                       [&field, &func](std::size_t id) { func(field, id); });
      ++iblock;
    }
  }

  // call forEach(FlagFManager, flag, [&](FieldType& field, std::size_t id){});
  template <typename FlagFieldType, typename Func>
  void forEach(const BlockFieldManager<FlagFieldType, FloatType, Dim>& FlagFManager,
               std::uint8_t flag, Func func) {
    int iblock = 0;
    for (Block<FloatType, Dim>& blockgeo : _BlockGeo.getBlocks()) {
      auto& field = _Fields[iblock];
      const auto& flagarr = FlagFManager.getBlockField(iblock).getField(0);
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
#ifdef MPI_ENABLED
    MPINormalCommunicate(count);
#endif
    AverCommunicate(count);
#ifdef MPI_ENABLED
    MPIAverCommunicate(count);
#endif
    InterpCommunicate(count);
#ifdef MPI_ENABLED
    MPIInterpCommunicate(count);
#endif
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
#ifdef MPI_ENABLED
    MPINormalCommunicate();
#endif
    AverCommunicate();
#ifdef MPI_ENABLED
    MPIAverCommunicate();
#endif
    InterpCommunicate();
#ifdef MPI_ENABLED
    MPIInterpCommunicate();
#endif
  }

#ifdef MPI_ENABLED

  void MPINormalCommunicate(std::int64_t count) {
    mpi().barrier();
    // --- send data ---
    std::vector<MPI_Request> SendRequests;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if ((count % (static_cast<int>(pow(2, deLevel))) == 0) &&
          blockF.getBlock()._NeedMPIComm) {
        blockF.mpiNormalSend(SendRequests);
      }
    }
    // --- receive data ---
    std::vector<MPI_Request> RecvRequests;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if ((count % (static_cast<int>(pow(2, deLevel))) == 0) &&
          blockF.getBlock()._NeedMPIComm) {
        blockF.mpiNormalRecv(RecvRequests);
      }
    }
    // wait for all send requests to complete
    MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
    // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
    // --- wait and set field data ---
    int reqidx = 0;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if ((count % (static_cast<int>(pow(2, deLevel))) == 0) &&
          blockF.getBlock()._NeedMPIComm) {
        blockF.mpiNormalSet(reqidx, RecvRequests);
      }
    }
  }

  // send data to lower level, deLevel+1
  void MPIAverSend(std::int64_t count, std::vector<MPI_Request>& SendRequests) {
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if (blockF.getBlock().getLevel() != std::uint8_t(0)) {
        if ((count % (static_cast<int>(pow(2, deLevel + 1))) == 0) &&
            blockF.getBlock()._NeedMPIComm) {
          blockF.mpiAverSend(SendRequests);
        }
      }
    }
  }
  void MPIAverRecv(std::int64_t count, std::vector<MPI_Request>& RecvRequests) {
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if ((count % (static_cast<int>(pow(2, deLevel))) == 0) &&
          blockF.getBlock()._NeedMPIComm) {
        blockF.mpiAverRecv(RecvRequests);
      }
    }
  }
  void MPIAverSet(std::int64_t count, std::vector<MPI_Request>& RecvRequests) {
    int reqidx = 0;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if ((count % (static_cast<int>(pow(2, deLevel))) == 0) &&
          blockF.getBlock()._NeedMPIComm) {
        blockF.mpiAverSet(reqidx, RecvRequests);
      }
    }
  }

  void MPIAverCommunicate(std::int64_t count) {
    mpi().barrier();
    // --- send data ---
    std::vector<MPI_Request> SendRequests;
    MPIAverSend(count, SendRequests);
    // --- receive data ---
    std::vector<MPI_Request> RecvRequests;
    MPIAverRecv(count, RecvRequests);
    // wait for all send requests to complete
    MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
    // --- wait and set field data ---
    MPIAverSet(count, RecvRequests);
  }
  // send data to higher level, deLevel-1
  void MPIInterpSend(std::int64_t count, std::vector<MPI_Request>& SendRequests) {
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel()) - 1;
      if (deLevel != -1) {
        if ((count % (static_cast<int>(pow(2, deLevel))) == 0) &&
            blockF.getBlock()._NeedMPIComm) {
          blockF.mpiIntpSend(SendRequests);
        }
      }
    }
  }
  void MPIInterpRecv(std::int64_t count, std::vector<MPI_Request>& RecvRequests) {
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if ((count % (static_cast<int>(pow(2, deLevel))) == 0) &&
          blockF.getBlock()._NeedMPIComm) {
        blockF.mpiIntpRecv(RecvRequests);
      }
    }
  }
  void MPIInterpSet(std::int64_t count, std::vector<MPI_Request>& RecvRequests) {
    int reqidx = 0;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      const int deLevel =
        static_cast<int>(_BlockGeo.getMaxLevel() - blockF.getBlock().getLevel());
      if ((count % (static_cast<int>(pow(2, deLevel))) == 0) &&
          blockF.getBlock()._NeedMPIComm) {
        blockF.mpiIntpSet(reqidx, RecvRequests);
      }
    }
  }
  void MPIInterpCommunicate(std::int64_t count) {
    mpi().barrier();
    // --- send data ---
    std::vector<MPI_Request> SendRequests;
    MPIInterpSend(count, SendRequests);
    // --- receive data ---
    std::vector<MPI_Request> RecvRequests;
    MPIInterpRecv(count, RecvRequests);
    // wait for all send requests to complete
    MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
    // --- wait and set field data ---
    MPIInterpSet(count, RecvRequests);
  }

  void MPINormalCommunicate() {
    mpi().barrier();
    // --- send data ---
    std::vector<MPI_Request> SendRequests;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      if (blockF.getBlock()._NeedMPIComm) {
        blockF.mpiNormalSend(SendRequests);
      }
    }
    // --- receive data ---
    std::vector<MPI_Request> RecvRequests;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      if (blockF.getBlock()._NeedMPIComm) {
        blockF.mpiNormalRecv(RecvRequests);
      }
    }
    // wait for all send requests to complete
    MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
    // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
    // --- wait and set field data ---
    int reqidx = 0;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      if (blockF.getBlock()._NeedMPIComm) {
        blockF.mpiNormalSet(reqidx, RecvRequests);
      }
    }
  }
  void MPIAverCommunicate() {
    mpi().barrier();
    // --- send data ---
    std::vector<MPI_Request> SendRequests;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      if (blockF.getBlock()._NeedMPIComm) {
        blockF.mpiAverSend(SendRequests);
      }
    }
    // --- receive data ---
    std::vector<MPI_Request> RecvRequests;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      if (blockF.getBlock()._NeedMPIComm) {
        blockF.mpiAverRecv(RecvRequests);
      }
    }
    // wait for all send requests to complete
    MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
    // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
    // --- wait and set field data ---
    int reqidx = 0;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      if (blockF.getBlock()._NeedMPIComm) {
        blockF.mpiAverSet(reqidx, RecvRequests);
      }
    }
  }
  void MPIInterpCommunicate() {
    mpi().barrier();
    // --- send data ---
    std::vector<MPI_Request> SendRequests;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      if (blockF.getBlock()._NeedMPIComm) {
        blockF.mpiIntpSend(SendRequests);
      }
    }
    // --- receive data ---
    std::vector<MPI_Request> RecvRequests;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      if (blockF.getBlock()._NeedMPIComm) {
        blockF.mpiIntpRecv(RecvRequests);
      }
    }
    // wait for all send requests to complete
    MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
    // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
    // --- wait and set field data ---
    int reqidx = 0;
    for (BlockField<FieldType, FloatType, Dim>& blockF : _Fields) {
      if (blockF.getBlock()._NeedMPIComm) {
        blockF.mpiIntpSet(reqidx, RecvRequests);
      }
    }
  }
#endif
};


template <typename T, typename LatSet, typename Pack>
class BlockFieldManagerCollection;

template <typename T, typename LatSet, typename... Fields>
class BlockFieldManagerCollection<T, LatSet, TypePack<Fields...>> {
 public:
  static constexpr std::size_t FieldNum = sizeof...(Fields);

  BlockFieldManagerCollection() : fields(std::make_tuple(Fields()...)) {}
  template <typename BLOCKGEOMETRY>
  BlockFieldManagerCollection(BLOCKGEOMETRY& blockgeo)
      : fields(std::make_tuple(BlockFieldManager<Fields, T, LatSet::d>(blockgeo)...)) {}

  template <typename BLOCKGEOMETRY, typename... InitValues>
  BlockFieldManagerCollection(BLOCKGEOMETRY& blockgeo,
                              const std::tuple<InitValues...>& initvalues)
      : fields(make_fields(blockgeo, initvalues, std::make_index_sequence<FieldNum>())) {}

  ~BlockFieldManagerCollection() = default;

  template <typename FieldType>
  auto& getField() {
    return std::get<BlockFieldManager<FieldType, T, LatSet::d>>(fields);
  }
  template <typename FieldType>
  const auto& getField() const {
    return std::get<BlockFieldManager<FieldType, T, LatSet::d>>(fields);
  }

  auto& getFields() { return fields; }
  template <std::size_t... Is>
  auto get_ith(std::index_sequence<Is...>, int i) {
    return std::make_tuple(&(std::get<Is>(fields).getBlockField(i))...);
  }
  auto get_ith(int i) { return get_ith(std::make_index_sequence<FieldNum>(), i); }

  // for each field, call func
  template <typename Func>
  void forEachField(Func&& func) {
    forEachFieldImpl(std::forward<Func>(func), std::make_index_sequence<FieldNum>());
  }
  template <typename Func, std::size_t... Is>
  void forEachFieldImpl(Func&& func, std::index_sequence<Is...>) {
    (func(std::get<Is>(fields), std::integral_constant<std::size_t, Is>{}), ...);
  }

 private:
  std::tuple<BlockFieldManager<Fields, T, LatSet::d>...> fields;

  template <typename BLOCKGEOMETRY, typename... InitValues, std::size_t... Is>
  auto make_fields(BLOCKGEOMETRY& blockgeo, const std::tuple<InitValues...>& initvalues,
                   std::index_sequence<Is...>) {
    return std::make_tuple(
      BlockFieldManager<Fields, T, LatSet::d>(blockgeo, std::get<Is>(initvalues))...);
  }
};

template <typename T, typename LatSet, typename Pack>
class BlockFieldManagerPtrCollection;

template <typename T, typename LatSet, typename... Fields>
class BlockFieldManagerPtrCollection<T, LatSet, TypePack<Fields...>> {
 public:
  static constexpr std::size_t FieldNum = sizeof...(Fields);

  BlockFieldManagerPtrCollection() : fields() {}
  template <typename First, typename... Rest>
  BlockFieldManagerPtrCollection(First* first, Rest*... rest) : fields(first, rest...) {}
  ~BlockFieldManagerPtrCollection() = default;

  void Init(Fields&... fieldrefs) { fields = std::make_tuple(&fieldrefs...); }

  template <typename FieldType>
  static constexpr unsigned int is_at_index() {
    return is_same_at_index<FieldType, Fields...>::value;
  }
  // not real 'add', just assign the valid field pointer to an existing null pointer
  template <typename FieldType>
  void addField(BlockFieldManager<FieldType, T, LatSet::d>& field) {
    // std::get<is_same_at_index<FieldType, Fields...>::value>(fields) = &field;
    std::get<BlockFieldManager<FieldType, T, LatSet::d>*>(fields) = &field;
  }

  template <typename FieldType>
  auto& getField() {
    return *(std::get<BlockFieldManager<FieldType, T, LatSet::d>*>(fields));
  }
  template <typename FieldType>
  const auto& getField() const {
    return *(std::get<BlockFieldManager<FieldType, T, LatSet::d>*>(fields));
  }

  auto& getFields() { return fields; }

  template <std::size_t... Is>
  auto get_ith(std::index_sequence<Is...>, int i) {
    // return std::make_tuple(&(std::get<Is>(fields)->getBlockField(i))...);
    return std::make_tuple(
      (std::get<Is>(fields) ? &(std::get<Is>(fields)->getBlockField(i)) : nullptr)...);
  }
  auto get_ith(int i) { return get_ith(std::make_index_sequence<FieldNum>(), i); }

 private:
  std::tuple<BlockFieldManager<Fields, T, LatSet::d>*...> fields;
};

template <typename T, typename LatSet, typename Pack>
struct ExtractFieldPtrs {};

template <typename T, typename LatSet, typename... Fields, typename... FieldPtrs>
struct ExtractFieldPtrs<T, LatSet,
                        TypePack<TypePack<Fields...>, TypePack<FieldPtrs...>>> {
  static auto getFieldPtrTuple(
    int i, BlockFieldManagerCollection<T, LatSet, TypePack<Fields...>>& BFMCol,
    BlockFieldManagerPtrCollection<T, LatSet, TypePack<FieldPtrs...>>& BFMPtrCol) {
    return std::tuple_cat(BFMCol.get_ith(i), BFMPtrCol.get_ith(i));
  }
};

template <typename T, typename LatSet, typename... Fields>
struct ExtractFieldPtrs<T, LatSet, TypePack<Fields...>> {
  static auto getFieldPtrTuple(
    int i, BlockFieldManagerCollection<T, LatSet, TypePack<Fields...>>& BFMCol) {
    return BFMCol.get_ith(i);
  }
  static auto getFieldPtrTuple(
    int i, BlockFieldManagerCollection<T, LatSet, TypePack<Fields...>>& BFMCol,
    BlockFieldManagerPtrCollection<T, LatSet, TypePack<>>& BFMPtrCol) {
    return BFMCol.get_ith(i);
  }
};


// vector of Genericvector
template <typename T>
class GenericvectorManager{
  private:
  // generic vectors
  std::vector<Genericvector<T>> _vectors;

 public:
  using datatype = T;
  using array_type = Genericvector<T>;

  GenericvectorManager() = default;
  GenericvectorManager(std::size_t size) : _vectors(size) {}
  GenericvectorManager(std::size_t size, T initvalue) : _vectors(size, initvalue) {}
  template <typename FLAGFM, typename FLAGTYPE>
  GenericvectorManager(std::size_t size, FLAGFM& FlagFM, FLAGTYPE Flag, std::string name = "GenericvectorManager"): _vectors(size) {
    Init(FlagFM, Flag);
    // get vector info
    std::size_t sumsize{};
    for(auto& vec : _vectors){
      sumsize += vec.size();
    }
    std::cout <<"[" << name <<"Num of indices]: " << sumsize << std::endl;
  }
  // copy constructor
  GenericvectorManager(const GenericvectorManager& vecManager) : _vectors(vecManager._vectors) {}
  // move constructor
  GenericvectorManager(GenericvectorManager&& vecManager) noexcept
      : _vectors(std::move(vecManager._vectors)) {}
  // copy assignment operator
  GenericvectorManager& operator=(const GenericvectorManager& vecManager) {
    if (this != &vecManager) {
      _vectors = vecManager._vectors;
    }
    return *this;
  }
  // move assignment operator
  GenericvectorManager& operator=(GenericvectorManager&& vecManager) noexcept {
    if (this != &vecManager) {
      _vectors = std::move(vecManager._vectors);
    }
    return *this;
  }

  ~GenericvectorManager() = default;

  void Init(std::size_t size) { _vectors.resize(size); }

  // init using flags
  template <typename FLAGFM, typename FLAGTYPE>
  void Init(FLAGFM& FlagFM, FLAGTYPE Flag){
    FlagFM.forEachField([&](auto& blockfield){
    auto& block = blockfield.getBlock();
    auto& VecIds = getvector(block.getBlockId()).getvector();
    block.forEach([&](std::size_t id){
      if (util::isFlag(blockfield.get(id), Flag)) {
        VecIds.push_back(id);
      }});
  });
  }
  
  Genericvector<T>& getvector(int i) { return _vectors[i]; }
  const Genericvector<T>& getvector(int i) const { return _vectors[i]; }

  std::vector<Genericvector<T>>& getallvectors() { return _vectors; }
  const std::vector<Genericvector<T>>& getallvectors() const { return _vectors; }

  std::vector<std::vector<T>*> getallvectorptrs() const {
    std::vector<std::vector<T>*> vecptrs;
    for (auto& vec : _vectors) {
      vecptrs.push_back(&vec.getvector());
    }
    return vecptrs;
  }
};