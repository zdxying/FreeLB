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

#ifdef __CUDACC__
#include "utils/cuda_device.h"
#endif


namespace cudev{

#ifdef __CUDACC__

template <typename T, unsigned int D>
struct BlockComm {
  std::size_t* SendCells;  // index in sendblock
  std::size_t* RecvCells;  // index in recvblock

  __any__ BlockComm(std::size_t* sendcells, std::size_t* recvcells)
      : SendCells(sendcells), RecvCells(recvcells) {}
};

// communication structure for blocks of different (refine) levels
// heterogeneous communication/ interpolation
template <typename T, unsigned int D>
struct IntpBlockComm {
  // index in sendblock, unfold to 1D array
  std::size_t* SendCells;
  // index in recvblock
  std::size_t* RecvCells;
  // predefined interp weight
  static constexpr T InterpWeight2D[4][4] = {
    {T(0.0625), T(0.1875), T(0.1875), T(0.5625)},
    {T(0.1875), T(0.0625), T(0.5625), T(0.1875)},
    {T(0.1875), T(0.5625), T(0.0625), T(0.1875)},
    {T(0.5625), T(0.1875), T(0.1875), T(0.0625)}
  };
  static constexpr T InterpWeight3D[8][8] = {
    {T(0.015625), T(0.046875), T(0.046875), T(0.140625), T(0.046875), T(0.140625), T(0.140625), T(0.421875)},
    {T(0.046875), T(0.015625), T(0.140625), T(0.046875), T(0.140625), T(0.046875), T(0.421875), T(0.140625)},
    {T(0.046875), T(0.140625), T(0.015625), T(0.046875), T(0.140625), T(0.421875), T(0.046875), T(0.140625)},
    {T(0.140625), T(0.046875), T(0.046875), T(0.015625), T(0.421875), T(0.140625), T(0.140625), T(0.046875)},
    {T(0.046875), T(0.140625), T(0.140625), T(0.421875), T(0.015625), T(0.046875), T(0.046875), T(0.140625)},
    {T(0.140625), T(0.046875), T(0.421875), T(0.140625), T(0.046875), T(0.015625), T(0.140625), T(0.046875)},
    {T(0.140625), T(0.421875), T(0.046875), T(0.140625), T(0.046875), T(0.140625), T(0.015625), T(0.046875)},
    {T(0.421875), T(0.140625), T(0.140625), T(0.046875), T(0.140625), T(0.046875), T(0.046875), T(0.015625)}
  };

  __device__ static constexpr auto getIntpWeight() {
    if constexpr (D == 2) {
      return InterpWeight2D;
    } else {
      return InterpWeight3D;
    }
  }

  __any__ IntpBlockComm(std::size_t* sendcells, std::size_t* recvcells)
      : SendCells(sendcells), RecvCells(recvcells) {}

  // return (D == 2) ? 4 : 8;
  __device__ static constexpr int getSourceNum() { return (D == 2) ? 4 : 8; }
  // return (D == 2) ? T(0.25) : T(0.125);
  __device__ static constexpr T getUniformWeight() { return (D == 2) ? T(0.25) : T(0.125); }
};

#endif

} // namespace cudev

// neighbor block direction relative to current block
enum NbrDirection : std::uint8_t {
  // not a neighbor
  NONE = 0,
  // x-negative
  XN = 1,
  // x-positive
  XP = 2,
  // y-negative
  YN = 4,
  // y-positive
  YP = 8,
  // z-negative
  ZN = 16,
  // z-positive
  ZP = 32
};

template <unsigned int D>
NbrDirection getCornerNbrDirection(int corner) {
  if (corner == -1) {
    return NbrDirection::NONE;
  }
  std::uint8_t direction{};
  if constexpr (D == 2) {
    if (corner == 0) {
      direction |= NbrDirection::XN | NbrDirection::YN;
    } else if (corner == 1) {
      direction |= NbrDirection::XP | NbrDirection::YN;
    } else if (corner == 2) {
      direction |= NbrDirection::XN | NbrDirection::YP;
    } else if (corner == 3) {
      direction |= NbrDirection::XP | NbrDirection::YP;
    }
  } else if constexpr (D == 3) {
    if (corner == 0) {
      direction |= NbrDirection::XN | NbrDirection::YN | NbrDirection::ZN;
    } else if (corner == 1) {
      direction |= NbrDirection::XP | NbrDirection::YN | NbrDirection::ZN;
    } else if (corner == 2) {
      direction |= NbrDirection::XN | NbrDirection::YP | NbrDirection::ZN;
    } else if (corner == 3) {
      direction |= NbrDirection::XP | NbrDirection::YP | NbrDirection::ZN;
    } else if (corner == 4) {
      direction |= NbrDirection::XN | NbrDirection::YN | NbrDirection::ZP;
    } else if (corner == 5) {
      direction |= NbrDirection::XP | NbrDirection::YN | NbrDirection::ZP;
    } else if (corner == 6) {
      direction |= NbrDirection::XN | NbrDirection::YP | NbrDirection::ZP;
    } else if (corner == 7) {
      direction |= NbrDirection::XP | NbrDirection::YP | NbrDirection::ZP;
    }
  }
  return static_cast<NbrDirection>(direction);
}

template <unsigned int D>
NbrDirection getEdgeNbrDirection(int edge) {
  if (edge == -1) {
    return NbrDirection::NONE;
  }
  std::uint8_t direction{};
  if constexpr (D == 2) {
    if (edge == 0) {
      direction |= NbrDirection::XN;
    } else if (edge == 1) {
      direction |= NbrDirection::XP;
    } else if (edge == 2) {
      direction |= NbrDirection::YN;
    } else if (edge == 3) {
      direction |= NbrDirection::YP;
    }
  } else if constexpr (D == 3) {
    if (edge == 0) {
      direction |= NbrDirection::XN | NbrDirection::YN;
    } else if (edge == 1) {
      direction |= NbrDirection::XP | NbrDirection::YN;
    } else if (edge == 2) {
      direction |= NbrDirection::XN | NbrDirection::YP;
    } else if (edge == 3) {
      direction |= NbrDirection::XP | NbrDirection::YP;
    } else if (edge == 4) {
      direction |= NbrDirection::XN | NbrDirection::ZN;
    } else if (edge == 5) {
      direction |= NbrDirection::XP | NbrDirection::ZN;
    } else if (edge == 6) {
      direction |= NbrDirection::XN | NbrDirection::ZP;
    } else if (edge == 7) {
      direction |= NbrDirection::XP | NbrDirection::ZP;
    } else if (edge == 8) {
      direction |= NbrDirection::YN | NbrDirection::ZN;
    } else if (edge == 9) {
      direction |= NbrDirection::YP | NbrDirection::ZN;
    } else if (edge == 10) {
      direction |= NbrDirection::YN | NbrDirection::ZP;
    } else if (edge == 11) {
      direction |= NbrDirection::YP | NbrDirection::ZP;
    }
  }
  return static_cast<NbrDirection>(direction);
}


NbrDirection getFaceNbrDirection(int face) {
  if (face == -1) {
    return NbrDirection::NONE;
  }
  std::uint8_t direction{};

  if (face == 0) {
    direction |= NbrDirection::XN;
  } else if (face == 1) {
    direction |= NbrDirection::XP;
  } else if (face == 2) {
    direction |= NbrDirection::YN;
  } else if (face == 3) {
    direction |= NbrDirection::YP;
  } else if (face == 4) {
    direction |= NbrDirection::ZN;
  } else if (face == 5) {
    direction |= NbrDirection::ZP;
  }
  return static_cast<NbrDirection>(direction);
}

// input two AABB should have no intersection
// not correct for now
template <typename T, unsigned int D>
NbrDirection getNbrDirection(const AABB<T, D>& refAABB, const AABB<T, D>& nbrAABB) {
  std::uint8_t direction{};
  // NbrDirection direction_ = NbrDirection::NONE;

  const Vector<T, D>& refcenter = refAABB.getCenter();
  const Vector<T, D>& nbrcenter = nbrAABB.getCenter();
  const Vector<T, D>& refextension = refAABB.getExtension();
  const Vector<T, D>& nbrextension = nbrAABB.getExtension();
  const Vector<T, D> diff = nbrcenter - refcenter;
  const Vector<T, D> absdiff = diff.getabs();
  const Vector<T, D> ext = (refextension + nbrextension) / T{2};
  
  constexpr T epsilon = std::numeric_limits<T>::epsilon();
  // check if is neighbor
  for (unsigned int i = 0; i < D; ++i) {
    if (absdiff[i] > ext[i] + epsilon) {
      return NbrDirection::NONE;
    }
  }

  if (diff[0] < -epsilon) {
    direction |= NbrDirection::XN;
  } else if (diff[0] > epsilon){
    direction |= NbrDirection::XP;
  }
  if (diff[1] < -epsilon) {
    direction |= NbrDirection::YN;
  } else if (diff[1] > epsilon){
    direction |= NbrDirection::YP;
  }

  if constexpr (D == 3) {
    if (diff[2] < -epsilon) {
      direction |= NbrDirection::ZN;
    } else if (diff[2] > epsilon){
      direction |= NbrDirection::ZP;
    }
  }

  return static_cast<NbrDirection>(direction);
}

std::vector<int> getCornerIdxFromNbr(std::vector<std::size_t>& cornerId, const std::vector<std::size_t>& allcornerIdx, const std::vector<NbrDirection>& nbrdirs) {
  std::vector<int> corners{};
  const int size = allcornerIdx.size();
  cornerId.clear();
  // 2d
  constexpr std::uint8_t xnyn = NbrDirection::XN | NbrDirection::YN;
  constexpr std::uint8_t xpyn = NbrDirection::XP | NbrDirection::YN;
  constexpr std::uint8_t xnyp = NbrDirection::XN | NbrDirection::YP;
  constexpr std::uint8_t xpyp = NbrDirection::XP | NbrDirection::YP;

  // 3d
  constexpr std::uint8_t xnynzn = NbrDirection::XN | NbrDirection::YN | NbrDirection::ZN;
  constexpr std::uint8_t xpynzn = NbrDirection::XP | NbrDirection::YN | NbrDirection::ZN;
  constexpr std::uint8_t xnypzn = NbrDirection::XN | NbrDirection::YP | NbrDirection::ZN;
  constexpr std::uint8_t xpypzn = NbrDirection::XP | NbrDirection::YP | NbrDirection::ZN;

  constexpr std::uint8_t xnynzp = NbrDirection::XN | NbrDirection::YN | NbrDirection::ZP;
  constexpr std::uint8_t xpynzp = NbrDirection::XP | NbrDirection::YN | NbrDirection::ZP;
  constexpr std::uint8_t xnypzp = NbrDirection::XN | NbrDirection::YP | NbrDirection::ZP;
  constexpr std::uint8_t xpypzp = NbrDirection::XP | NbrDirection::YP | NbrDirection::ZP;

  // collect all nbr directions
  std::uint8_t flag{};
  for (const auto& dir : nbrdirs) flag |= dir;

  if (size == 4) {
    // 2d
    if (util::isFlag(flag, xnyn)) {cornerId.push_back(allcornerIdx[0]); corners.push_back(0);}
    if (util::isFlag(flag, xpyn)) {cornerId.push_back(allcornerIdx[1]); corners.push_back(1);}
    if (util::isFlag(flag, xnyp)) {cornerId.push_back(allcornerIdx[2]); corners.push_back(2);}
    if (util::isFlag(flag, xpyp)) {cornerId.push_back(allcornerIdx[3]); corners.push_back(3);}
  } else if (size == 8) {
    // 3d
    if (util::isFlag(flag, xnynzn)) {cornerId.push_back(allcornerIdx[0]); corners.push_back(0);}
    if (util::isFlag(flag, xpynzn)) {cornerId.push_back(allcornerIdx[1]); corners.push_back(1);}
    if (util::isFlag(flag, xnypzn)) {cornerId.push_back(allcornerIdx[2]); corners.push_back(2);}
    if (util::isFlag(flag, xpypzn)) {cornerId.push_back(allcornerIdx[3]); corners.push_back(3);}
    if (util::isFlag(flag, xnynzp)) {cornerId.push_back(allcornerIdx[4]); corners.push_back(4);}
    if (util::isFlag(flag, xpynzp)) {cornerId.push_back(allcornerIdx[5]); corners.push_back(5);}
    if (util::isFlag(flag, xnypzp)) {cornerId.push_back(allcornerIdx[6]); corners.push_back(6);}
    if (util::isFlag(flag, xpypzp)) {cornerId.push_back(allcornerIdx[7]); corners.push_back(7);}
  }
  return corners;
}

// communication structure for block geometry
template <typename T, unsigned int D>
struct BlockComm {
  std::vector<std::size_t> SendCells;  // index in sendblock
  std::vector<std::size_t> RecvCells;  // index in recvblock
  // neighbor block
  Block<T, D>* SendBlock;
  // neighbor direction
  NbrDirection Direction;

  #ifdef __CUDACC__
  std::size_t* dev_SendCells;
  std::size_t* dev_RecvCells;
  cudev::BlockComm<T, D>* dev_BlockComm;
  #endif

  BlockComm(Block<T, D>* block) : SendBlock(block), Direction(NbrDirection::NONE) {}
  int getSendId() const { return SendBlock->getBlockId(); }
};

// communication structure for blocks of different (refine) levels
// heterogeneous communication/ interpolation
template <typename T, unsigned int D>
struct IntpBlockComm {
  // index in sendblock
  std::vector<IntpSource<D>> SendCells;
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

  // neighbor direction
  NbrDirection Direction;

  static constexpr auto getIntpWeight() {
    if constexpr (D == 2) {
      return InterpWeight2D;
    } else {
      return InterpWeight3D;
    }
  }

  #ifdef __CUDACC__
  std::size_t* dev_SendCells;
  std::size_t* dev_RecvCells;
  cudev::IntpBlockComm<T, D>* dev_BlockComm;
  #endif

  IntpBlockComm(Block<T, D>* block) : SendBlock(block), Direction(NbrDirection::NONE) {}

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
struct MPIIntpBlockSendStru {
  std::vector<IntpSource<Dim>> SendCells;
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

  MPIIntpBlockSendStru(int rank, int blockid) : RecvRank(rank), RecvBlockid(blockid) {}
};

template <typename FloatType, unsigned int Dim>
struct MPIIntpBlockComm {
  std::vector<MPIIntpBlockSendStru<FloatType, Dim>> Senders;
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
void MPIBlockBufferInit(const MPIIntpBlockComm<FloatType,Dim>& MPIComm, MPIBlockBuffer<T>& MPIBuffer,
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
