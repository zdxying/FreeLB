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

// communicator.h

#pragma once

#include <vector>

#include "parallel/mpi_manager.h"
#include "utils/alias.h"

// communication structure for block geometry
template <typename T, unsigned int D>
struct BlockComm {
  std::vector<std::size_t> SendCells;  // index in sendblock
  std::vector<std::size_t> RecvCells;  // index in recvblock
  Block<T, D>* SendBlock;

  BlockComm(Block<T, D>* block) : SendBlock(block) {}
  int getSendId() const { return SendBlock->getBlockId(); }
};

// communication structure for blocks of different (refine) levels
// heterogeneous communication/ interpolation
template <typename T, unsigned int D>
struct InterpBlockComm {
  // index in sendblock
  std::vector<InterpSource<D>> SendCells;
  // interp weight
  // std::vector<InterpWeight<T, D>> InterpWeights;
  // index in recvblock
  std::vector<std::size_t> RecvCells;
  // receive data from sendblock
  Block<T, D>* SendBlock;
  // predefined interp weight
  static constexpr std::array<InterpWeight<T, 2>, 4> InterpWeight2D{
    {{T(0.0625), T(0.1875), T(0.1875), T(0.5625)},
     {T(0.1875), T(0.0625), T(0.5625), T(0.1875)},
     {T(0.1875), T(0.5625), T(0.0625), T(0.1875)},
     {T(0.5625), T(0.1875), T(0.1875), T(0.0625)}}};
  static constexpr std::array<InterpWeight<T, 3>, 8> InterpWeight3D{
   {{T{0.015625}, T{0.046875}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.140625}, T{0.421875}},
    {T{0.046875}, T{0.015625}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.421875}, T{0.140625}},
    {T{0.046875}, T{0.140625}, T{0.015625}, T{0.046875}, T{0.140625}, T{0.421875}, T{0.046875}, T{0.140625}},
    {T{0.140625}, T{0.046875}, T{0.046875}, T{0.015625}, T{0.421875}, T{0.140625}, T{0.140625}, T{0.046875}},
    {T{0.046875}, T{0.140625}, T{0.140625}, T{0.421875}, T{0.015625}, T{0.046875}, T{0.046875}, T{0.140625}},
    {T{0.140625}, T{0.046875}, T{0.421875}, T{0.140625}, T{0.046875}, T{0.015625}, T{0.140625}, T{0.046875}},
    {T{0.140625}, T{0.421875}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.015625}, T{0.046875}},
    {T{0.421875}, T{0.140625}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.046875}, T{0.015625}}}};

  static constexpr auto getIntpWeight() {
    if constexpr (D == 2) {
      return InterpWeight2D;
    } else {
      return InterpWeight3D;
    }
  }

  InterpBlockComm(Block<T, D>* block) : SendBlock(block) {}

  int getSendId() const { return SendBlock->getBlockId(); }
  // return (D == 2) ? 4 : 8;
  static constexpr int getSourceNum() { return (D == 2) ? 4 : 8; }
  // return (D == 2) ? T(0.25) : T(0.125);
  static constexpr T getUniformWeight() { return (D == 2) ? T(0.25) : T(0.125); }
};

struct MPIBlockSendStru {
  std::vector<std::size_t> SendCells;
  // receive rank
  int RecvRank;
  // receive block id
  int RecvBlockid;

  MPIBlockSendStru(int rank, int blockid) : RecvRank(rank), RecvBlockid(blockid) {}
};

struct MPIBlockRecvStru {
  std::vector<std::size_t> RecvCells;
  // send rank
  int SendRank;
  // receive block id
  int SendBlockid;

  MPIBlockRecvStru(int rank, int blockid) : SendRank(rank), SendBlockid(blockid) {}
};

// a collection of MPIBlockSendStru and MPIBlockRecvStru
// std::vector<MPIBlockSendStru<T, D>> Senders;
// std::vector<MPIBlockRecvStru<T, D>> Recvers;
struct MPIBlockComm {
  std::vector<MPIBlockSendStru> Senders;
  std::vector<MPIBlockRecvStru> Recvers;
  void clear() {
    Senders.clear();
    Recvers.clear();
  }
};

template <typename T, unsigned int Dim>
struct MPIInterpBlockSendStru {
  std::vector<InterpSource<Dim>> SendCells;
  int RecvRank;
  int RecvBlockid;
  // predefined interp weight
  static constexpr std::array<InterpWeight<T, 2>, 4> InterpWeight2D{
    {{T(0.0625), T(0.1875), T(0.1875), T(0.5625)},
     {T(0.1875), T(0.0625), T(0.5625), T(0.1875)},
     {T(0.1875), T(0.5625), T(0.0625), T(0.1875)},
     {T(0.5625), T(0.1875), T(0.1875), T(0.0625)}}};
  static constexpr std::array<InterpWeight<T, 3>, 8> InterpWeight3D{
    {T{0.015625}, T{0.046875}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.140625}, T{0.421875}},
    {T{0.046875}, T{0.015625}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.421875}, T{0.140625}},
    {T{0.046875}, T{0.140625}, T{0.015625}, T{0.046875}, T{0.140625}, T{0.421875}, T{0.046875}, T{0.140625}},
    {T{0.140625}, T{0.046875}, T{0.046875}, T{0.015625}, T{0.421875}, T{0.140625}, T{0.140625}, T{0.046875}},
    {T{0.046875}, T{0.140625}, T{0.140625}, T{0.421875}, T{0.015625}, T{0.046875}, T{0.046875}, T{0.140625}},
    {T{0.140625}, T{0.046875}, T{0.421875}, T{0.140625}, T{0.046875}, T{0.015625}, T{0.140625}, T{0.046875}},
    {T{0.140625}, T{0.421875}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.015625}, T{0.046875}},
    {T{0.421875}, T{0.140625}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.046875}, T{0.015625}}};

  static constexpr auto getIntpWeight() {
    return (Dim == 2) ? InterpWeight2D : InterpWeight3D;
  }
  static constexpr int getSourceNum() { return (Dim == 2) ? 4 : 8; }
  static constexpr T getUniformWeight() {
    return (Dim == 2) ? T(0.25) : T(0.125);
  }

  MPIInterpBlockSendStru(int rank, int blockid) : RecvRank(rank), RecvBlockid(blockid) {}
};

template <typename FloatType, unsigned int Dim>
struct MPIInterpBlockComm {
  std::vector<MPIInterpBlockSendStru<FloatType, Dim>> Senders;
  std::vector<MPIBlockRecvStru> Recvers;
  void clear() {
    Senders.clear();
    Recvers.clear();
  }
};

// a collection of buffers
// std::vector<std::vector<T>> SendBuffers;
// std::vector<std::vector<T>> RecvBuffers;
template <typename T>
struct MPIBlockBuffer {
  std::vector<std::vector<T>> SendBuffers;
  std::vector<std::vector<T>> RecvBuffers;
};

template <typename T>
void MPIBlockBufferInit(const MPIBlockComm& MPIComm, MPIBlockBuffer<T>& MPIBuffer,
                        unsigned int Size = 1) {
  MPIBuffer.SendBuffers.resize(MPIComm.Senders.size(), std::vector<T>{});
  for (std::size_t i = 0; i < MPIComm.Senders.size(); ++i) {
    MPIBuffer.SendBuffers[i].resize(Size * MPIComm.Senders[i].SendCells.size(), T{});
  }
  MPIBuffer.RecvBuffers.resize(MPIComm.Recvers.size(), std::vector<T>{});
  for (std::size_t i = 0; i < MPIComm.Recvers.size(); ++i) {
    MPIBuffer.RecvBuffers[i].resize(Size * MPIComm.Recvers[i].RecvCells.size(), T{});
  }
}

template <typename FloatType, unsigned int Dim, typename T>
void MPIBlockBufferInit(const MPIInterpBlockComm<FloatType,Dim>& MPIComm, MPIBlockBuffer<T>& MPIBuffer,
                        unsigned int Size = 1) {
  MPIBuffer.SendBuffers.resize(MPIComm.Senders.size(), std::vector<T>{});
  for (std::size_t i = 0; i < MPIComm.Senders.size(); ++i) {
    MPIBuffer.SendBuffers[i].resize(Size * MPIComm.Senders[i].SendCells.size(), T{});
  }
  MPIBuffer.RecvBuffers.resize(MPIComm.Recvers.size(), std::vector<T>{});
  for (std::size_t i = 0; i < MPIComm.Recvers.size(); ++i) {
    MPIBuffer.RecvBuffers[i].resize(Size * MPIComm.Recvers[i].RecvCells.size(), T{});
  }
}

// add data to send buffer
// scalar data
template <typename ArrayType, typename T>
void addtoBuffer(std::vector<T>& buffer, const ArrayType& arr,
                 const std::vector<int>& index) {
  for (int i : index) {
    buffer.push_back(arr[i]);
  }
}
// pop data
template <typename T, unsigned int q>
void addtoBuffer(std::vector<T>& buffer, const PopulationField<T, q>& poparr,
                 const std::vector<int>& index) {
  for (int k = 0; k < q; ++k) {
    const CyclicArray<T>& arrk = poparr.getField(k);
    for (int i : index) {
      buffer.push_back(arrk[i]);
    }
  }
}
