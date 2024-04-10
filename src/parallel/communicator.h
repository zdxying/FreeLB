/* This file is part of FreeLB
 * 
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 * 
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 * 
 */

// communicator.h

#pragma once

#include <vector>

#include "parallel/mpi_manager.h"
#include "utils/alias.h"

struct MPIBlockSendStru {
  std::vector<std::size_t> SendCells;
  // send data buffer
  // std::vector<T> SendBuffer;
  // receive rank
  int RecvRank;
  // receive block id
  int RecvBlockid;

  MPIBlockSendStru(int rank, int blockid) : RecvRank(rank), RecvBlockid(blockid) {}
};


struct MPIBlockRecvStru {
  std::vector<std::size_t> RecvCells;
  // recv data buffer
  // std::vector<T> RecvBuffer;
  // send rank
  int SendRank;
  // receive block id
  int SendBlockid;

  MPIBlockRecvStru(int rank, int blockid) : SendRank(rank), SendBlockid(blockid) {}
};

// a collection of MPIBlockSendStru and MPIBlockRecvStru
// std::vector<MPIBlockSendStru<T, D>> Senders;
// std::vector<MPIBlockRecvStru<T, D>> Recvers;
struct MPIBlockCommStru {
  std::vector<MPIBlockSendStru> Senders;
  std::vector<MPIBlockRecvStru> Recvers;
};

template <unsigned int D>
struct MPIBlockInterpSendStru {
  std::vector<InterpSource<D>> SendCells;
  int RecvRank;
  int RecvBlockid;
  // constructor
  MPIBlockInterpSendStru(int rank, int blockid) : RecvRank(rank), RecvBlockid(blockid) {}
};

template <unsigned int D>
struct MPIInterpBlockCommStru {
  std::vector<MPIBlockInterpSendStru<D>> Senders;
  std::vector<MPIBlockRecvStru> Recvers;
};

// a collection of buffers
// std::vector<std::vector<T>> SendBuffers;
// std::vector<std::vector<T>> RecvBuffers;
template <typename T>
struct MPIBlockBuffer {
  std::vector<std::vector<T>> SendBuffers;
  std::vector<std::vector<T>> RecvBuffers;
};

namespace mpi {
// helper to init MPIBlockBuffer
template <typename T>
void MPIBlockBufferInit(const MPIBlockCommStru& MPIComm, MPIBlockBuffer<T>& MPIBuffer,
                        int Size = 1) {
  MPIBuffer.SendBuffers.resize(MPIComm.Senders.size(), std::vector<T>{});
  for (int i = 0; i < MPIComm.Senders.size(); ++i) {
    MPIBuffer.SendBuffers[i].resize(Size * MPIComm.Senders[i].SendCells.size(), T(0));
  }
  MPIBuffer.RecvBuffers.resize(MPIComm.Recvers.size(), std::vector<T>{});
  for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
    MPIBuffer.RecvBuffers[i].resize(Size * MPIComm.Recvers[i].RecvCells.size(), T(0));
  }
  Mpi().barrier();
}
// add data to send buffer
// scaler data
template <typename ArrayType, typename T>
void addtoBuffer(std::vector<T>& buffer, const ArrayType& arr, const std::vector<int>& index) {
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


}  // namespace mpi
