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

// block_lattice.hh

#pragma once

#include "data_struct/block_lattice.h"

template <typename T, typename LatSet, typename TypePack>
template <typename... FIELDPTRS>
BlockLattice<T, LatSet, TypePack>::BlockLattice(Block<T, LatSet::d>& block,
                                                AbstractConverter<T>& conv,
                                                std::tuple<FIELDPTRS...> fieldptrs)
    : BlockLatticeBase<T, LatSet, TypePack>(block, fieldptrs),
      Omega(RefineConverter<T>::getOmegaF(conv.getOMEGA(), block.getLevel())) {
  _Omega = T{1} - Omega;
  fOmega = T{1} - T{0.5} * Omega;

  std::size_t i{};
  CommDirection.resize(block.getCommunicator().Comm.Comms.size(), std::vector<unsigned int>{});
  for (SharedComm& comm : block.getCommunicator().Comm.Comms) {
    getCommPopDir<LatSet>(comm.Direction, CommDirection[i]);
    ++i;
  }
#ifdef MPI_ENABLED
  i = 0;
  RecvDirection.resize(block.getCommunicator().DirRecvs.size(), std::vector<unsigned int>{});
  for (DistributedComm& comm : block.getCommunicator().DirRecvs) {
    getCommPopDir<LatSet>(comm.Direction, RecvDirection[i]);
    ++i;
  }
  i = 0;
  SendDirection.resize(block.getCommunicator().DirSends.size(), std::vector<unsigned int>{});
  for (DistributedComm& comm : block.getCommunicator().DirSends) {
    getCommPopDir<LatSet>(comm.Direction, SendDirection[i]);
    ++i;
  }
#endif

#ifdef __CUDACC__
  InitDeviceData();
#endif
}


template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::Stream() {
  for (unsigned int i = 1; i < LatSet::q; ++i) {
    // this->template getField<POP<T, LatSet::q>>().getField(i).rotate(this->Delta_Index[i]);
    this->template getField<POP<T, LatSet::q>>().getField(i).rotate();
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename ArrayType>
void BlockLattice<T, LatSet, TypePack>::ApplyCellDynamics(const ArrayType& flagarr) {
  Cell<T, LatSet, TypePack> cell(0, *this);
  const std::size_t voxNum = this->getVoxNum();
  #ifdef _VOX_ENABLED
  const VoxelMap& map = this->getBlock().getVoxelMap();
  #endif
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
  for (std::size_t id = 0; id < voxNum; ++id) {
    #ifdef _VOX_ENABLED
    cell.setId(map[id]);
    CELLDYNAMICS::Execute(flagarr[map[id]], cell);
    #else
    cell.setId(id);
    CELLDYNAMICS::Execute(flagarr[id], cell);
    #endif
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLattice<T, LatSet, TypePack>::ApplyCellDynamics() {
  Cell<T, LatSet, TypePack> cell(0, *this);
  const std::size_t voxNum = this->getVoxNum();
  #ifdef _VOX_ENABLED
  const VoxelMap& map = this->getBlock().getVoxelMap();
  #endif
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
  for (std::size_t id = 0; id < voxNum; ++id) {
    #ifdef _VOX_ENABLED
    cell.setId(map[id]);
    #else
    cell.setId(id);
    CELLDYNAMICS::apply(cell);
    #endif
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLattice<T, LatSet, TypePack>::ApplyCellDynamics(const Genericvector<std::size_t>& Idx) {
  Cell<T, LatSet, TypePack> cell(0, *this);
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
  for (std::size_t id : Idx.getvector()) {
    cell.setId(id);
    CELLDYNAMICS::apply(cell);
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename DYNAMICS, typename elementType>
void BlockLattice<T, LatSet, TypePack>::ApplyDynamics(const Genericvector<elementType>& Idx) {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
  for (auto& element : Idx.getvector()) {
    DYNAMICS::apply(element, *this);
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename ArrayType>
void BlockLattice<T, LatSet, TypePack>::ApplyInnerCellDynamics(const ArrayType& flagarr) {
  Cell<T, LatSet, TypePack> cell(0, *this);
  if constexpr (LatSet::d == 2) {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
    for (int j = this->getOverlap(); j < this->getNy() - this->getOverlap(); ++j) {
      std::size_t id = j * this->getNx() + this->getOverlap();
      for (int i = this->getOverlap(); i < this->getNx() - this->getOverlap(); ++i) {
        cell.setId(id);
        CELLDYNAMICS::Execute(flagarr[id], cell);
        ++id;
      }
    }
  } else if constexpr (LatSet::d == 3) {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
    for (int k = this->getOverlap(); k < this->getNz() - this->getOverlap(); ++k) {
      for (int j = this->getOverlap(); j < this->getNy() - this->getOverlap(); ++j) {
        std::size_t id = k * this->getProjection()[2] + j * this->getProjection()[1] +
                         this->getOverlap();
        for (int i = this->getOverlap(); i < this->getNx() - this->getOverlap(); ++i) {
          cell.setId(id);
          CELLDYNAMICS::Execute(flagarr[id], cell);
          ++id;
        }
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLattice<T, LatSet, TypePack>::ApplyInnerCellDynamics() {
  Cell<T, LatSet, TypePack> cell(0, *this);
  if constexpr (LatSet::d == 2) {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
    for (int j = this->getOverlap(); j < this->getNy() - this->getOverlap(); ++j) {
      std::size_t id = j * this->getNx() + this->getOverlap();
      for (int i = this->getOverlap(); i < this->getNx() - this->getOverlap(); ++i) {
        cell.setId(id);
        CELLDYNAMICS::apply(cell);
        ++id;
      }
    }
  } else if constexpr (LatSet::d == 3) {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
    for (int k = this->getOverlap(); k < this->getNz() - this->getOverlap(); ++k) {
      for (int j = this->getOverlap(); j < this->getNy() - this->getOverlap(); ++j) {
        std::size_t id = k * this->getProjection()[2] + j * this->getProjection()[1] +
                         this->getOverlap();
        for (int i = this->getOverlap(); i < this->getNx() - this->getOverlap(); ++i) {
          cell.setId(id);
          CELLDYNAMICS::apply(cell);
          ++id;
        }
      }
    }
  }
}


#ifdef __CUDACC__

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::CuDevStream() {
  CuDevStreamKernel<<<LatSet::q, 1>>>(dev_BlockLat);
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename ArrayType>
void BlockLattice<T, LatSet, TypePack>::CuDevApplyCellDynamics(ArrayType& flagarr) {
  const unsigned int blockSize = THREADS_PER_BLOCK;
  const unsigned int blockNum = (this->getN() + blockSize - 1) / blockSize;
  CuDevApplyCellDynamicsKernel<T, LatSet, cudev_TypePack, CELLDYNAMICS, typename ArrayType::cudev_array_type><<<blockNum, blockSize>>>(dev_BlockLat, flagarr.get_devObj());
  // CuDevApplyCellDynamicsKernel<T, LatSet, cudev_TypePack, CELLDYNAMICS, typename ArrayType::cudev_array_type><<<blockNum, blockSize>>>(dev_BlockLat, flagarr.get_devObj(), this->getN());
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLattice<T, LatSet, TypePack>::CuDevApplyCellDynamics() {
  const unsigned int blockSize = THREADS_PER_BLOCK;
  const unsigned int blockNum = (this->getN() + blockSize - 1) / blockSize;
  CuDevApplyCellDynamicsKernel<T, LatSet, cudev_TypePack, CELLDYNAMICS><<<blockNum, blockSize>>>(dev_BlockLat);
  // CuDevApplyCellDynamicsKernel<T, LatSet, cudev_TypePack, CELLDYNAMICS><<<blockNum, blockSize>>>(dev_BlockLat, this->getN());
}

#endif


template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::EnableToleranceRho(T rhores) {
  RhoRes = rhores;
  RhoOld.Resize(this->getN());
  for (std::size_t i = 0; i < this->getN(); ++i)
    RhoOld[i] = this->template getField<GenericRho>().get(i);
}

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::EnableToleranceU(T ures) {
  URes = ures;
  UOld.Resize(this->getN());
  for (std::size_t i = 0; i < this->getN(); ++i)
    UOld[i] = this->template getField<VELOCITY<T, LatSet::d>>().get(i);
}

template <typename T, typename LatSet, typename TypePack>
T BlockLattice<T, LatSet, TypePack>::getToleranceRho() {
  T maxres{};
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : maxres)
#endif
  for (std::size_t i = 0; i < this->getN(); ++i) {
    T res = std::abs(this->template getField<GenericRho>().get(i) - RhoOld[i]);
    maxres = res > maxres ? res : maxres;
    RhoOld[i] = this->template getField<GenericRho>().get(i);
  }
  return maxres;
}

template <typename T, typename LatSet, typename TypePack>
T BlockLattice<T, LatSet, TypePack>::getToleranceU() {
  T maxres{};
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : maxres)
#endif
  for (std::size_t i = 0; i < this->getN(); ++i) {
    T res0 = std::abs(this->template getField<VELOCITY<T, LatSet::d>>().get(i)[0] - UOld[i][0]);
    T res1 = std::abs(this->template getField<VELOCITY<T, LatSet::d>>().get(i)[1] - UOld[i][1]);
    if constexpr (LatSet::d == 3) {
      T res2 = std::abs(this->template getField<VELOCITY<T, LatSet::d>>().get(i)[2] - UOld[i][2]);
      res1 = res1 > res2 ? res1 : res2;
    }
    T res = res0 > res1 ? res0 : res1;
    maxres = res > maxres ? res : maxres;
    // set UOld
    UOld[i][0] = this->template getField<VELOCITY<T, LatSet::d>>().get(i)[0];
    UOld[i][1] = this->template getField<VELOCITY<T, LatSet::d>>().get(i)[1];
    if constexpr (LatSet::d == 3)
      UOld[i][2] = this->template getField<VELOCITY<T, LatSet::d>>().get(i)[2];
  }
  return maxres;
}

template <typename T, typename LatSet, typename TypePack>
T BlockLattice<T, LatSet, TypePack>::getTolRho(int shift) {
  T maxres{};
  if constexpr (LatSet::d == 2) {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : maxres)
#endif
    for (int j = shift; j < this->getNy() - shift; ++j) {
      for (int i = shift; i < this->getNx() - shift; ++i) {
        std::size_t id = j * this->getNx() + i;
        T res = std::abs(this->template getField<GenericRho>().get(id) - RhoOld[id]);
        maxres = res > maxres ? res : maxres;
        RhoOld[id] = this->template getField<GenericRho>().get(id);
      }
    }
  } else if constexpr (LatSet::d == 3) {
  std::size_t NxNy = this->getNx() * this->getNy();
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : maxres)
#endif
    for (int k = shift; k < this->getNz() - shift; ++k) {
      for (int j = shift; j < this->getNy() - shift; ++j) {
        for (int i = shift; i < this->getNx() - shift; ++i) {
          std::size_t id = k * NxNy + j * this->getNx() + i;
          T res = std::abs(this->template getField<GenericRho>().get(id) - RhoOld[id]);
          maxres = res > maxres ? res : maxres;
          RhoOld[id] = this->template getField<GenericRho>().get(id);
        }
      }
    }
  }
  return maxres;
}

template <typename T, typename LatSet, typename TypePack>
T BlockLattice<T, LatSet, TypePack>::getTolU(int shift) {
  T maxres{};
  if constexpr (LatSet::d == 2) {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : maxres) 
#endif
    for (int j = shift; j < this->getNy() - shift; ++j) {
      for (int i = shift; i < this->getNx() - shift; ++i) {
        std::size_t id = j * this->getNx() + i;
        T res0 = std::abs(this->template getField<VELOCITY<T, LatSet::d>>().get(id)[0] - UOld[id][0]);
        T res1 = std::abs(this->template getField<VELOCITY<T, LatSet::d>>().get(id)[1] - UOld[id][1]);
        T res = res0 > res1 ? res0 : res1;
        maxres = res > maxres ? res : maxres;
        // set UOld
        UOld[id][0] = this->template getField<VELOCITY<T, LatSet::d>>().get(id)[0];
        UOld[id][1] = this->template getField<VELOCITY<T, LatSet::d>>().get(id)[1];
      }
    }
  } else if constexpr (LatSet::d == 3) {
    std::size_t NxNy = this->getNx() * this->getNy();
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : maxres)
#endif
    for (int k = shift; k < this->getNz() - shift; ++k) {
      for (int j = shift; j < this->getNy() - shift; ++j) {
        for (int i = shift; i < this->getNx() - shift; ++i) {
          std::size_t id = k * NxNy + j * this->getNx() + i;
          T res0 = std::abs(this->template getField<VELOCITY<T, LatSet::d>>().get(id)[0] - UOld[id][0]);
          T res1 = std::abs(this->template getField<VELOCITY<T, LatSet::d>>().get(id)[1] - UOld[id][1]);
          T res2 = std::abs(this->template getField<VELOCITY<T, LatSet::d>>().get(id)[2] - UOld[id][2]);
          res1 = res1 > res2 ? res1 : res2;
          T res = res0 > res1 ? res0 : res1;
          maxres = res > maxres ? res : maxres;
          // set UOld
          UOld[id][0] = this->template getField<VELOCITY<T, LatSet::d>>().get(id)[0];
          UOld[id][1] = this->template getField<VELOCITY<T, LatSet::d>>().get(id)[1];
          UOld[id][2] = this->template getField<VELOCITY<T, LatSet::d>>().get(id)[2];
        }
      }
    }
  }
  return maxres;
}


#ifdef MPI_ENABLED

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::mpiNormalSend(std::vector<MPI_Request>& SendRequests, 
std::vector<std::vector<T>>& SendBuffers, const std::vector<DistributedComm>& MPISends) {
  // add to send buffer
  for (std::size_t i = 0; i < MPISends.size(); ++i) {
    const DistributedComm& comm = MPISends[i];
    std::vector<T>& buffer = SendBuffers[i];
    buffer.resize(comm.Cells.size() * SendDirection[i].size());
    const std::vector<std::size_t>& sends = comm.Cells;
    std::size_t bufidx{};
    for (unsigned int iArr : SendDirection[i]) {
      const auto& Array = this->template getField<POP<T, LatSet::q>>().getField(iArr);
      for (std::size_t id : sends) {
        buffer[bufidx] = Array[id];
        ++bufidx;
      }
    }
    // non-blocking send
    MPI_Request request;
    mpi().iSend(buffer.data(), buffer.size(), comm.TargetRank, &request, comm.Tag);
    SendRequests.push_back(request);
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::mpiNormalFullSend(std::vector<MPI_Request>& SendRequests, 
std::vector<std::vector<T>>& SendBuffers, const std::vector<DistributedComm>& MPISends) {
  // add to send buffer
  for (std::size_t i = 0; i < MPISends.size(); ++i) {
    const DistributedComm& comm = MPISends[i];
    std::vector<T>& buffer = SendBuffers[i];
    buffer.resize(comm.Cells.size() * LatSet::q);
    const std::vector<std::size_t>& sends = comm.Cells;
    std::size_t bufidx{};
    for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
      const auto& Array = this->template getField<POP<T, LatSet::q>>().getField(iArr);
      for (std::size_t id : sends) {
        buffer[bufidx] = Array[id];
        ++bufidx;
      }
    }
    // non-blocking send
    MPI_Request request;
    mpi().iSend(buffer.data(), buffer.size(), comm.TargetRank, &request, comm.TargetBlockId);
    SendRequests.push_back(request);
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::mpiAverSend(std::vector<MPI_Request>& SendRequests,
std::vector<std::vector<T>>& SendBuffers, const std::vector<DistributedComm>& MPISends) {
  static constexpr std::size_t SendSize = LatSet::d == 2 ? 4 : 8;

  const T OmegaF = this->getOmega();
  const GenericArray<T>& RhoArr = this->template getField<GenericRho>().getField(0);
  const GenericArray<Vector<T, LatSet::d>>& UArr = this->template getField<VELOCITY<T, LatSet::d>>().getField(0);

  // add to send buffer
  for (std::size_t i = 0; i < MPISends.size(); ++i) {
    const DistributedComm& comm = MPISends[i];
    std::vector<T>& buffer = SendBuffers[i];
    buffer.resize(comm.Cells.size() / SendSize * LatSet::q);
    const std::vector<std::size_t>& sends = comm.Cells;
    const std::size_t size = sends.size();
    std::size_t bufidx{};

    for (std::size_t i = 0; i < size; i += SendSize) {
      T averRho = getAverage<T, LatSet::d>(RhoArr, i, sends);
      Vector<T, LatSet::d> averU = getAverage<T, LatSet::d>(UArr, i, sends);

      std::array<T, LatSet::q> feq{};
      equilibrium::SecondOrder<CellType>::apply(feq, averRho, averU);

      for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
        T averPop = getAverage<T, LatSet::d>(this->template getField<POP<T, LatSet::q>>().getField(iArr), i, sends);
        buffer[bufidx] = RefineConverter<T>::getPopC(averPop, feq[iArr], OmegaF);
        ++bufidx;
      }
    }
    // non-blocking send
    MPI_Request request;
    mpi().iSend(buffer.data(), buffer.size(), comm.TargetRank, &request, comm.TargetBlockId);
    SendRequests.push_back(request);
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::mpiIntpSend(std::vector<MPI_Request>& SendRequests, 
std::vector<std::vector<T>>& SendBuffers, const std::vector<DistributedComm>& MPISends) {
  static constexpr std::size_t SendSize = LatSet::d == 2 ? 4 : 8;

  const T OmegaC = this->getOmega();
  const GenericArray<T>& RhoArr = this->template getField<GenericRho>().getField(0);
  const GenericArray<Vector<T, LatSet::d>>& UArr = this->template getField<VELOCITY<T, LatSet::d>>().getField(0);

  // add to send buffer
  for (std::size_t i = 0; i < MPISends.size(); ++i) {
    const DistributedComm& comm = MPISends[i];
    std::vector<T>& buffer = SendBuffers[i];
    buffer.resize(comm.Cells.size() / SendSize * LatSet::q);
    const std::vector<std::size_t>& sends = comm.Cells;
    const std::size_t size = sends.size();
    std::size_t bufidx{};
    for (std::size_t i = 0; i < size;) {
      if constexpr (LatSet::d == 2) {
        popIntp<0>(OmegaC, sends, RhoArr, UArr, i            , buffer, bufidx); 
        popIntp<1>(OmegaC, sends, RhoArr, UArr, i += SendSize, buffer, bufidx); 
        popIntp<2>(OmegaC, sends, RhoArr, UArr, i += SendSize, buffer, bufidx); 
        popIntp<3>(OmegaC, sends, RhoArr, UArr, i += SendSize, buffer, bufidx); 
        i += SendSize;
      } else if constexpr (LatSet::d == 3) {
        popIntp<0>(OmegaC, sends, RhoArr, UArr, i            , buffer, bufidx); 
        popIntp<1>(OmegaC, sends, RhoArr, UArr, i += SendSize, buffer, bufidx); 
        popIntp<2>(OmegaC, sends, RhoArr, UArr, i += SendSize, buffer, bufidx); 
        popIntp<3>(OmegaC, sends, RhoArr, UArr, i += SendSize, buffer, bufidx); 
        popIntp<4>(OmegaC, sends, RhoArr, UArr, i += SendSize, buffer, bufidx); 
        popIntp<5>(OmegaC, sends, RhoArr, UArr, i += SendSize, buffer, bufidx); 
        popIntp<6>(OmegaC, sends, RhoArr, UArr, i += SendSize, buffer, bufidx); 
        popIntp<7>(OmegaC, sends, RhoArr, UArr, i += SendSize, buffer, bufidx); 
        i += SendSize;
      }
    }
    // non-blocking send
    MPI_Request request;
    mpi().iSend(buffer.data(), buffer.size(), comm.TargetRank, &request, comm.TargetBlockId);
    SendRequests.push_back(request);
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::mpiNormalRecv(std::vector<MPI_Request>& RecvRequests, 
std::vector<std::vector<T>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs) {
  // non-blocking recv
  for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
    const DistributedComm& comm = MPIRecvs[i];
    std::vector<T>& buffer = RecvBuffers[i];
    buffer.resize(comm.Cells.size() * RecvDirection[i].size());
    MPI_Request request;
    mpi().iRecv(buffer.data(), buffer.size(), comm.TargetRank, &request, comm.Tag);
    RecvRequests.push_back(request);
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::mpiFullRecv(std::vector<MPI_Request>& RecvRequests, 
std::vector<std::vector<T>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs) {
  // non-blocking recv
  for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
    const DistributedComm& comm = MPIRecvs[i];
    std::vector<T>& buffer = RecvBuffers[i];
    buffer.resize(comm.Cells.size() * LatSet::q);
    MPI_Request request;
    mpi().iRecv(buffer.data(), buffer.size(), comm.TargetRank, &request, this->BlockGeo.getBlockId());
    RecvRequests.push_back(request);
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::mpiNormalSet(int& reqidx, std::vector<MPI_Request>& RecvRequests,
const std::vector<std::vector<T>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs) {
  // wait and set field data
  for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
    MPI_Wait(&RecvRequests[i + reqidx], MPI_STATUS_IGNORE);
    const DistributedComm& comm = MPIRecvs[i];
    const std::vector<T>& buffer = RecvBuffers[i];
    const std::vector<std::size_t>& recvs = comm.Cells;
    std::size_t bufidx{};
    for (unsigned int iArr : RecvDirection[i]) {
      auto& Array = this->template getField<POP<T, LatSet::q>>().getField(iArr);
      for (std::size_t id : recvs) {
        Array[id] = buffer[bufidx];
        ++bufidx;
      }
    }
  }
  reqidx += MPIRecvs.size();
}

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::mpiNormalFullSet(int& reqidx, std::vector<MPI_Request>& RecvRequests,
const std::vector<std::vector<T>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs) {
  // wait and set field data
  for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
    MPI_Wait(&RecvRequests[i + reqidx], MPI_STATUS_IGNORE);
    const DistributedComm& comm = MPIRecvs[i];
    const std::vector<T>& buffer = RecvBuffers[i];
    const std::vector<std::size_t>& recvs = comm.Cells;
    std::size_t bufidx{};
    for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
      auto& Array = this->template getField<POP<T, LatSet::q>>().getField(iArr);
      for (std::size_t id : recvs) {
        Array[id] = buffer[bufidx];
        ++bufidx;
      }
    }
  }
  reqidx += MPIRecvs.size();
}

template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::mpiAverIntpSet(int& reqidx, std::vector<MPI_Request>& RecvRequests,
const std::vector<std::vector<T>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs) {
  // wait and set field data
  for (std::size_t i = 0; i < MPIRecvs.size(); ++i) {
    MPI_Wait(&RecvRequests[i + reqidx], MPI_STATUS_IGNORE);
    const DistributedComm& comm = MPIRecvs[i];
    const std::vector<T>& buffer = RecvBuffers[i];
    const std::vector<std::size_t>& recvs = comm.Cells;
    std::size_t bufidx{};
    for (std::size_t id : recvs) {
      std::array<T*, LatSet::q> CellPop = this->template getField<POP<T, LatSet::q>>().getArray(id);
      for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
        *(CellPop[iArr]) = buffer[bufidx];
        ++bufidx;
      }
    }
  }
  reqidx += MPIRecvs.size();
}


template <typename T, typename LatSet, typename TypePack>
template <unsigned int D>
void BlockLattice<T, LatSet, TypePack>::popIntp(T OmegaC, const std::vector<std::size_t>& sends, 
const GenericArray<T>& RhoArr, const GenericArray<Vector<T, LatSet::d>>& UArr, std::size_t i, 
std::vector<T>& buffer, std::size_t& bufidx) {
  
  T averRho = getInterpolation<D, T, LatSet::d>(RhoArr, i, sends);
  Vector<T, LatSet::d> averU = getInterpolation<D, T, LatSet::d>(UArr, i, sends);

  std::array<T, LatSet::q> feq{};
  equilibrium::SecondOrder<CellType>::apply(feq, averRho, averU);

  for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
    T averPop = getInterpolation<D, T, LatSet::d>(this->template getField<POP<T, LatSet::q>>().getField(iArr), i, sends);
    buffer[bufidx] = RefineConverter<T>::getPopF(averPop, feq[iArr], OmegaC);
    ++bufidx;
  }
}

#endif

// BlockLatticeManager

template <typename T, typename LatSet, typename TypePack>
template <typename... FIELDPTRTYPES>
BlockLatticeManager<T, LatSet, TypePack>::BlockLatticeManager(
  BlockGeometry<T, LatSet::d>& blockgeo, AbstractConverter<T>& conv,
  FIELDPTRTYPES*... fieldptrs)
    : BlockLatticeManagerBase<T, LatSet, TypePack>(blockgeo, fieldptrs...), Conv(conv) {
  // init rho and pop
  if constexpr (this->template hasField<GenericRho>()) {
    this->template getField<GenericRho>().Init(Conv.getLatRhoInit());
  }

  // if constexpr (this->template hasField<POP<T, LatSet::q>>()) {
  //   this->template getField<POP<T, LatSet::q>>().forEachField([&](auto& field) {
  //     for (unsigned int i = 0; i < LatSet::q; ++i) {
  //       field.getField(i).Init(Conv.getLatRhoInit() * latset::w<LatSet>(i));
  //     }
  //   });
  // }

  // init block lattices
  for (std::size_t i = 0; i < this->BlockGeo.getBlocks().size(); ++i) {
    BlockLats.emplace_back(this->BlockGeo.getBlock(i), Conv,
                           ExtractFieldPtrs<T, LatSet, TypePack>::getFieldPtrTuple(
                             i, this->Fields, this->FieldPtrs));
  }
  if constexpr (this->template hasField<POP<T, LatSet::q>>()) {
    ForEachBlockLattice([&](auto& blocklat) {
      auto& field = blocklat.template getField<POP<T, LatSet::q>>();
      for (unsigned int i = 0; i < LatSet::q; ++i) {
        field.getField(i).Init(Conv.getLatRhoInit() * latset::w<LatSet>(i), blocklat.getDelta_Index()[i]);
      }
    });
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename INITVALUEPACK, typename... FIELDPTRTYPES>
BlockLatticeManager<T, LatSet, TypePack>::BlockLatticeManager(
  BlockGeometry<T, LatSet::d>& blockgeo, INITVALUEPACK& initvalues,
  AbstractConverter<T>& conv, FIELDPTRTYPES*... fieldptrs)
    : BlockLatticeManagerBase<T, LatSet, TypePack>(blockgeo, initvalues, fieldptrs...),
      Conv(conv) {
  // init pop
  // if constexpr (this->template hasField<POP<T, LatSet::q>>()) {
  //   this->template getField<POP<T, LatSet::q>>().forEachField([&](auto& field) {
  //     for (unsigned int i = 0; i < LatSet::q; ++i) {
  //       field.getField(i).Init(Conv.getLatRhoInit() * latset::w<LatSet>(i));
  //     }
  //   });
  // }

  // init block lattices
  for (std::size_t i = 0; i < this->BlockGeo.getBlocks().size(); ++i) {
    BlockLats.emplace_back(this->BlockGeo.getBlock(i), Conv,
                           ExtractFieldPtrs<T, LatSet, TypePack>::getFieldPtrTuple(
                             i, this->Fields, this->FieldPtrs));
  }
  if constexpr (this->template hasField<POP<T, LatSet::q>>()) {
    ForEachBlockLattice([&](auto& blocklat) {
      auto& field = blocklat.template getField<POP<T, LatSet::q>>();
      for (unsigned int i = 0; i < LatSet::q; ++i) {
        field.getField(i).Init(Conv.getLatRhoInit() * latset::w<LatSet>(i), blocklat.getDelta_Index()[i]);
      }
    });
  }
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::Init() {
  BlockLats.clear();
  for (std::size_t i = 0; i < this->BlockGeo.getBlocks().size(); ++i) {
    BlockLats.emplace_back(this->BlockGeo.getBlock(i), Conv,
                           ExtractFieldPtrs<T, LatSet, TypePack>::getFieldPtrTuple(
                             i, this->Fields, this->FieldPtrs));
  }
  if constexpr (this->template hasField<POP<T, LatSet::q>>()) {
    ForEachBlockLattice([&](auto& blocklat) {
      auto& field = blocklat.template getField<POP<T, LatSet::q>>();
      for (unsigned int i = 0; i < LatSet::q; ++i) {
        field.getField(i).setOffset(blocklat.getDelta_Index()[i]);
      }
    });
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename... InitValues>
void BlockLatticeManager<T, LatSet, TypePack>::Init(
  BlockGeometryHelper<T, LatSet::d>& GeoHelper,
  const std::tuple<InitValues...>& initvalues) {
  this->Fields.forEachField(
    // [&]<std::size_t index>(auto& field, std::integral_constant<std::size_t, index>) {
    [&](auto& field, auto int_const) {
      if constexpr (field.isField) {
        // field.InitAndComm(GeoHelper, std::get<index>(initvalues));
        field.InitAndComm(GeoHelper, std::get<int_const.value>(initvalues));
        // field.InitAndComm(GeoHelper, initvalues.template get<int_const.value>());
      } else {
        // field.NonFieldInit(std::get<index>(initvalues));
        field.NonFieldInit(std::get<int_const.value>(initvalues));
        // field.NonFieldInit(initvalues.template get<int_const.value>());
      }
    });
  BlockLats.clear();
  for (std::size_t i = 0; i < this->BlockGeo.getBlocks().size(); ++i) {
    BlockLats.emplace_back(this->BlockGeo.getBlock(i), Conv,
                           ExtractFieldPtrs<T, LatSet, TypePack>::getFieldPtrTuple(
                             i, this->Fields, this->FieldPtrs));
  }
  if constexpr (this->template hasField<POP<T, LatSet::q>>()) {
    ForEachBlockLattice([&](auto& blocklat) {
      auto& field = blocklat.template getField<POP<T, LatSet::q>>();
      for (unsigned int i = 0; i < LatSet::q; ++i) {
        field.getField(i).setOffset(blocklat.getDelta_Index()[i]);
      }
    });
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::Stream(std::int64_t count) {
  for (auto& BLat : BlockLats) {
    if (count % (static_cast<int>(std::pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.Stream();
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::Stream() {
  for (auto& BLat : BlockLats) {BLat.Stream();}
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename FieldType>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyCellDynamics(
  std::int64_t count, const BlockFieldManager<FieldType, T, LatSet::d>& BFM) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      BlockLats[i]
        .template ApplyCellDynamics<CELLDYNAMICS, typename FieldType::array_type>(
          BFM.getBlockField(i).getField(0));
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename FieldType>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyCellDynamics(
  const BlockFieldManager<FieldType, T, LatSet::d>& BFM) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
// #pragma omp parallel for proc_bind(spread) num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
      BlockLats[i].template ApplyCellDynamics<CELLDYNAMICS, typename FieldType::array_type>(
          BFM.getBlockField(i).getField(0));
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyCellDynamics(std::int64_t count) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      BlockLats[i].template ApplyCellDynamics<CELLDYNAMICS>();
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyCellDynamics() {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
      BlockLats[i].template ApplyCellDynamics<CELLDYNAMICS>();
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyCellDynamics(
  std::int64_t count, const GenericvectorManager<std::size_t>& blockids) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      BlockLats[i].template ApplyCellDynamics<CELLDYNAMICS>(blockids.getvector(i));
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyCellDynamics(
  const GenericvectorManager<std::size_t>& blockids) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
      BlockLats[i].template ApplyCellDynamics<CELLDYNAMICS>(blockids.getvector(i));
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename DYNAMICS, typename elementType>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyDynamics(
  std::int64_t count, const GenericvectorManager<elementType>& blockids) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      BlockLats[i].template ApplyDynamics<DYNAMICS, elementType>(blockids.getvector(i));
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename DYNAMICS, typename elementType>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyDynamics(
  const GenericvectorManager<elementType>& blockids) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
      BlockLats[i].template ApplyDynamics<DYNAMICS, elementType>(blockids.getvector(i));
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename FieldType>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyInnerCellDynamics(
  std::int64_t count, const BlockFieldManager<FieldType, T, LatSet::d>& BFM) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      BlockLats[i]
        .template ApplyInnerCellDynamics<CELLDYNAMICS, typename FieldType::array_type>(
          BFM.getBlockField(i).getField(0));
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename FieldType>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyInnerCellDynamics(
const BlockFieldManager<FieldType, T, LatSet::d>& BFM) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
      BlockLats[i].template ApplyInnerCellDynamics<CELLDYNAMICS, typename FieldType::array_type>(
          BFM.getBlockField(i).getField(0));
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyInnerCellDynamics(
  std::int64_t count) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
      BlockLats[i].template ApplyInnerCellDynamics<CELLDYNAMICS>();
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyInnerCellDynamics() {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t i = 0; i < BlockLats.size(); ++i) {
      BlockLats[i].template ApplyInnerCellDynamics<CELLDYNAMICS>();
  }
}

#ifdef __CUDACC__

// for now only one block is supported

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::CuDevStream(){
  BlockLats[0].CuDevStream();
  // for (auto& BLat : BlockLats) {
  //   BLat.CuDevStream();
  // }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename FieldType>
void BlockLatticeManager<T, LatSet, TypePack>::CuDevApplyCellDynamics(BlockFieldManager<FieldType, T, LatSet::d>& BFM){
  BlockLats[0].template CuDevApplyCellDynamics<CELLDYNAMICS, typename FieldType::array_type>(
          BFM.getBlockField(0).getField(0));
  // for (std::size_t i = 0; i < BlockLats.size(); ++i) {
  //     BlockLats[i].template CuDevApplyCellDynamics<CELLDYNAMICS, typename FieldType::array_type>(
  //         BFM.getBlockField(i).getField(0));
  // }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLatticeManager<T, LatSet, TypePack>::CuDevApplyCellDynamics(){
  BlockLats[0].template CuDevApplyCellDynamics<CELLDYNAMICS>();
  // for (std::size_t i = 0; i < BlockLats.size(); ++i) {
  //   BlockLats[i].template CuDevApplyCellDynamics<CELLDYNAMICS>();
  // }
}

#endif

#ifdef MPI_ENABLED

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPINormalCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<T>>> SendBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::vector<std::vector<std::vector<T>>> RecvBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends = blat.getBlock().getCommunicator().DirSends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<T>{});
      blat.mpiNormalSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().DirRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<T>{});
      blat.mpiNormalRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().DirRecvs;
      blat.mpiNormalSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPINormalFullCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<T>>> SendBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::vector<std::vector<std::vector<T>>> RecvBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends = blat.getBlock().getCommunicator().MPIComm.Sends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<T>{});
      blat.mpiNormalFullSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.Recvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<T>{});
      blat.mpiFullRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.Recvs;
      blat.mpiNormalFullSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPINormalAllCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<T>>> SendBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::vector<std::vector<std::vector<T>>> RecvBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends = blat.getBlock().getCommunicator().AllMPIComm.Sends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<T>{});
      blat.mpiNormalFullSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().AllMPIComm.Recvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<T>{});
      blat.mpiFullRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().AllMPIComm.Recvs;
      blat.mpiNormalFullSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPIAverageCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<T>>> SendBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::vector<std::vector<std::vector<T>>> RecvBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends = blat.getBlock().getCommunicator().MPIComm.AverSends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<T>{});
      blat.mpiAverSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.AverRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<T>{});
      blat.mpiFullRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.AverRecvs;
      blat.mpiAverIntpSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPIInterpolateCommunicate() {
  mpi().barrier();
  std::vector<std::vector<std::vector<T>>> SendBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::vector<std::vector<std::vector<T>>> RecvBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends = blat.getBlock().getCommunicator().MPIComm.IntpSends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<T>{});
      blat.mpiIntpSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.IntpRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<T>{});
      blat.mpiFullRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BLOCKLATTICE& blat : BlockLats) {
    if (blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.IntpRecvs;
      blat.mpiAverIntpSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPINormalCommunicate(std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<T>>> SendBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::vector<std::vector<std::vector<T>>> RecvBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends = blat.getBlock().getCommunicator().DirSends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<T>{});
      blat.mpiNormalSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().DirRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<T>{});
      blat.mpiNormalRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().DirRecvs;
      blat.mpiNormalSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPINormalFullCommunicate(std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<T>>> SendBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::vector<std::vector<std::vector<T>>> RecvBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends = blat.getBlock().getCommunicator().MPIComm.Sends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<T>{});
      blat.mpiNormalFullSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.Recvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<T>{});
      blat.mpiFullRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.Recvs;
      blat.mpiNormalFullSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPINormalAllCommunicate(std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<T>>> SendBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::vector<std::vector<std::vector<T>>> RecvBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::size_t iblock{};
  // --- send data ---
  std::vector<MPI_Request> SendRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Sends = blat.getBlock().getCommunicator().AllMPIComm.Sends;
      SendBuffers[iblock].resize(Sends.size(), std::vector<T>{});
      blat.mpiNormalFullSend(SendRequests, SendBuffers[iblock], Sends);
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().AllMPIComm.Recvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<T>{});
      blat.mpiFullRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().AllMPIComm.Recvs;
      blat.mpiNormalFullSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPIAverageCommunicate(std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<T>>> SendBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::vector<std::vector<std::vector<T>>> RecvBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::size_t iblock{};
  // --- send data --- send data to lower level, deLevel+1
  std::vector<MPI_Request> SendRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel()) + 1;
    if (blat.getLevel() != std::uint8_t(0)) {
      if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
        const std::vector<DistributedComm>& Sends = blat.getBlock().getCommunicator().MPIComm.AverSends;
        SendBuffers[iblock].resize(Sends.size(), std::vector<T>{});
        blat.mpiAverSend(SendRequests, SendBuffers[iblock], Sends);
      }
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.AverRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<T>{});
      blat.mpiFullRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.AverRecvs;
      blat.mpiAverIntpSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPIInterpolateCommunicate(std::int64_t count) {
  mpi().barrier();
  std::vector<std::vector<std::vector<T>>> SendBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::vector<std::vector<std::vector<T>>> RecvBuffers(
    BlockLats.size(), std::vector<std::vector<T>>{});
  std::size_t iblock{};
  // --- send data --- send data to higher level, deLevel-1
  std::vector<MPI_Request> SendRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel()) - 1;
    if (deLevel != -1) {
      if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
        const std::vector<DistributedComm>& Sends = blat.getBlock().getCommunicator().MPIComm.IntpSends;
        SendBuffers[iblock].resize(Sends.size(), std::vector<T>{});
        blat.mpiIntpSend(SendRequests, SendBuffers[iblock], Sends);
      }
    }
    ++iblock;
  }
  // --- receive data ---
  iblock = 0;
  std::vector<MPI_Request> RecvRequests;
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.IntpRecvs;
      RecvBuffers[iblock].resize(Recvs.size(), std::vector<T>{});
      blat.mpiFullRecv(RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
  // wait for all send requests to complete
  MPI_Waitall(SendRequests.size(), SendRequests.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequests.size(), RecvRequests.data(), MPI_STATUSES_IGNORE);
  // --- wait and set field data ---
  iblock = 0;
  int reqidx{};
  for (BLOCKLATTICE& blat : BlockLats) {
    const int deLevel = static_cast<int>(getMaxLevel() - blat.getLevel());
    if ((count % (static_cast<int>(std::pow(2, deLevel))) == 0) && blat.getBlock().getCommunicator()._NeedMPIComm) {
      const std::vector<DistributedComm>& Recvs = blat.getBlock().getCommunicator().MPIComm.IntpRecvs;
      blat.mpiAverIntpSet(reqidx, RecvRequests, RecvBuffers[iblock], Recvs);
    }
    ++iblock;
  }
}

#endif


// Although the new communication scheme of populations is efficient, it did not communicate all directions,
// when you use other post-stream processors(after stream and before communication) in cells to be communicated, 
// if they might use invalid pop in one direction to get pop in anohter direction which will not be communicated
// e.g. the half-way bounce-back boundary condition, 
// then you should use Lattice.getField<POP<T, LatSet::q>>().Communicate() instead
template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalCommunicate() {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (BLOCKLATTICE& BLat : BlockLats) {
    normalCommunicate(BLat);
  }
#ifdef MPI_ENABLED
  MPINormalCommunicate();
#endif
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalCommunicate(std::int64_t count) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (auto& BLat : BlockLats) {
    if (count % (static_cast<int>(std::pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0){
      normalCommunicate(BLat);
    }
  }
#ifdef MPI_ENABLED
  MPINormalCommunicate(count);
#endif
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalFullCommunicate() {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (BLOCKLATTICE& BLat : BlockLats) {
    normalFullCommunicate(BLat);
  }
#ifdef MPI_ENABLED
  MPINormalFullCommunicate();
#endif
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalFullCommunicate(std::int64_t count) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (auto& BLat : BlockLats) {
    if (count % (static_cast<int>(std::pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0){
      normalFullCommunicate(BLat);
    }
  }
#ifdef MPI_ENABLED
  MPINormalFullCommunicate(count);
#endif
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalAllCommunicate() {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (BLOCKLATTICE& BLat : BlockLats) {
    normalAllCommunicate(BLat);
  }
#ifdef MPI_ENABLED
  MPINormalAllCommunicate();
#endif
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalAllCommunicate(std::int64_t count) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (auto& BLat : BlockLats) {
    if (count % (static_cast<int>(std::pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0){
      normalAllCommunicate(BLat);
    }
  }
#ifdef MPI_ENABLED
  MPINormalAllCommunicate(count);
#endif
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::AverageCommunicate() {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (auto& BLat : BlockLats) {
    averCommunicate(BLat, BLat.getBlock().getCommunicator().Comm.AverComm);
  }
#ifdef MPI_ENABLED
  MPIAverageCommunicate();
#endif
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::AverageCommunicate(std::int64_t count) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (auto& BLat : BlockLats) {
    if (count % (static_cast<int>(std::pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0){
      averCommunicate(BLat, BLat.getBlock().getCommunicator().Comm.AverComm);
    }
  }
#ifdef MPI_ENABLED
  MPIAverageCommunicate(count);
#endif
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::InterpolateCommunicate() {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (auto& BLat : BlockLats) {
    intpCommunicate(BLat, BLat.getBlock().getCommunicator().Comm.IntpComm);
  }
#ifdef MPI_ENABLED
  MPIInterpolateCommunicate();
#endif
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::InterpolateCommunicate(std::int64_t count) {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (auto& BLat : BlockLats) {
    if (count % (static_cast<int>(std::pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0) {
      intpCommunicate(BLat, BLat.getBlock().getCommunicator().Comm.IntpComm);
    }
  }
#ifdef MPI_ENABLED
  MPIInterpolateCommunicate(count);
#endif
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::Communicate() {
  NormalCommunicate();
  AverageCommunicate();
  InterpolateCommunicate();
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::Communicate(std::int64_t count) {
  NormalCommunicate(count);
  AverageCommunicate(count);
  InterpolateCommunicate(count);
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::FullCommunicate() {
  NormalFullCommunicate();
  AverageCommunicate();
  InterpolateCommunicate();
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::FullCommunicate(std::int64_t count) {
  NormalFullCommunicate(count);
  AverageCommunicate(count);
  InterpolateCommunicate(count);
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::AllCommunicate() {
  NormalAllCommunicate();
  AverageCommunicate();
  InterpolateCommunicate();
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::AllCommunicate(std::int64_t count) {
  NormalAllCommunicate(count);
  AverageCommunicate(count);
  InterpolateCommunicate(count);
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::EnableToleranceRho(T rhores) {
  for (auto& BLat : BlockLats) BLat.EnableToleranceRho(rhores);
}
template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::EnableToleranceU(T ures) {
  for (auto& BLat : BlockLats) BLat.EnableToleranceU(ures);
}

template <typename T, typename LatSet, typename TypePack>
T BlockLatticeManager<T, LatSet, TypePack>::getToleranceRho(int shift) {
  T maxres = T(0);
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) reduction(max : maxres)
#endif
  for (auto& BLat : BlockLats) {
    T temp;
    if (shift == 0) {
      temp = BLat.getToleranceRho();
    } else if (shift == -1) {
      int autoshift = BLat.getLevel() == std::uint8_t(0) ? 1 : 2;
      temp = BLat.getTolRho(autoshift);
    } else {
      temp = BLat.getTolRho(shift);
    }
    maxres = std::max(temp, maxres);
  }
#ifdef MPI_ENABLED
  mpi().barrier();
  mpi().reduceAndBcast(maxres, MPI_MAX);
#endif
  return maxres;
}

template <typename T, typename LatSet, typename TypePack>
T BlockLatticeManager<T, LatSet, TypePack>::getToleranceU(int shift) {
  T maxres = T(0);
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) reduction(max : maxres)
#endif
  for (auto& BLat : BlockLats) {
    T temp;
    if (shift == 0) {
      temp = BLat.getToleranceU();
    } else if (shift == -1) {
      int autoshift = BLat.getLevel() == std::uint8_t(0) ? 1 : 2;
      temp = BLat.getTolU(autoshift);
    } else {
      temp = BLat.getTolU(shift);
    }
    maxres = std::max(temp, maxres);
  }
#ifdef MPI_ENABLED
  mpi().barrier();
  mpi().reduceAndBcast(maxres, MPI_MAX);
#endif
  return maxres;
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::normalCommunicate(BLOCKLATTICE& BLat) {
  const Communicator& communicator = BLat.getBlock().getCommunicator();

  std::size_t icomm{};
  for (const SharedComm& comm : communicator.Comm.Comms) {
    const BLOCKLATTICE& nBLat = BlockLats[this->BlockGeo.findBlockIndex(comm.SendBlockId)];
    const std::size_t size = comm.SendRecvCells.size();
    for (unsigned int k : BLat.getCommDirection()[icomm]){
      const auto& nPopsk = nBLat.template getField<POP<T, LatSet::q>>().getField(k);
      auto& Popsk = BLat.template getField<POP<T, LatSet::q>>().getField(k);
      for (std::size_t i = 0; i < size; i+=2) {
        Popsk.set(comm.SendRecvCells[i+1], nPopsk[comm.SendRecvCells[i]]);
      }
    }
    ++icomm;
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::normalFullCommunicate(BLOCKLATTICE& BLat) {
  const Communicator& communicator = BLat.getBlock().getCommunicator();

  for (const SharedComm& comm : communicator.Comm.Comms) {
    const BLOCKLATTICE& nBLat = BlockLats[this->BlockGeo.findBlockIndex(comm.SendBlockId)];
    const std::size_t size = comm.SendRecvCells.size();

    for (unsigned int k = 0; k < LatSet::q; ++k){
      const auto& nPopsk = nBLat.template getField<POP<T, LatSet::q>>().getField(k);
      auto& Popsk = BLat.template getField<POP<T, LatSet::q>>().getField(k);
      for (std::size_t i = 0; i < size; i+=2) {
        Popsk.set(comm.SendRecvCells[i+1], nPopsk[comm.SendRecvCells[i]]);
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::normalAllCommunicate(BLOCKLATTICE& BLat) {
  const Communicator& communicator = BLat.getBlock().getCommunicator();

  for (const SharedComm& comm : communicator.AllComm.Comms) {
    const BLOCKLATTICE& nBLat = BlockLats[this->BlockGeo.findBlockIndex(comm.SendBlockId)];
    const std::size_t size = comm.SendRecvCells.size();

    for (unsigned int k = 0; k < LatSet::q; ++k){
      const auto& nPopsk = nBLat.template getField<POP<T, LatSet::q>>().getField(k);
      auto& Popsk = BLat.template getField<POP<T, LatSet::q>>().getField(k);
      for (std::size_t i = 0; i < size; i+=2) {
        Popsk.set(comm.SendRecvCells[i+1], nPopsk[comm.SendRecvCells[i]]);
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::averCommunicate(BLOCKLATTICE& BLat, const std::vector<SharedComm>& comms) {
  static constexpr std::size_t RecvOffset = LatSet::d == 2 ? 4 : 8;
  static constexpr std::size_t SendRecvOffset = RecvOffset + 1;

  for (const SharedComm& comm : comms) {
    const BLOCKLATTICE& nBLat = BlockLats[this->BlockGeo.findBlockIndex(comm.SendBlockId)];
    const std::size_t size = comm.SendRecvCells.size();
    const T OmegaF = nBLat.getOmega();
    const GenericArray<T>& nRhoF = nBLat.template getField<GenericRho>().getField(0);
    const GenericArray<Vector<T, LatSet::d>>& nUF =
      nBLat.template getField<VELOCITY<T, LatSet::d>>().getField(0);

    for (std::size_t i = 0; i < size; i += SendRecvOffset) {
      // get averaged(interpolated) rho and u
      T averRho = getAverage<T, LatSet::d>(nRhoF, i, comm.SendRecvCells);
      Vector<T, LatSet::d> averU = getAverage<T, LatSet::d>(nUF, i, comm.SendRecvCells);
      // get feq, peparing for pop conversion
      std::array<T, LatSet::q> feq{};
      Equilibrium<T, LatSet>::SecondOrder(feq, averU, averRho);

      std::array<T*, LatSet::q> CellPop =
        BLat.template getField<POP<T, LatSet::q>>().getArray(comm.SendRecvCells[i+RecvOffset]);
      for (unsigned int k = 0; k < LatSet::q; ++k) {
        T averpop = getAverage<T, LatSet::d>(
          nBLat.template getField<POP<T, LatSet::q>>().getField(k), i, comm.SendRecvCells);
        // convert from fine to coarse
        *(CellPop[k]) = RefineConverter<T>::getPopC(averpop, feq[k], OmegaF);
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::intpCommunicate(BLOCKLATTICE& BLat, const std::vector<SharedComm>& comms) {
  static constexpr std::size_t RecvOffset = LatSet::d == 2 ? 4 : 8;
  static constexpr std::size_t SendRecvOffset = RecvOffset + 1;

  for (const SharedComm& comm : comms) {
    const BLOCKLATTICE& nBLat = BlockLats[this->BlockGeo.findBlockIndex(comm.SendBlockId)];
    const std::size_t size = comm.SendRecvCells.size();
    const T OmegaC = nBLat.getOmega();
    const GenericArray<T>& nRhoF = nBLat.template getField<GenericRho>().getField(0);
    const GenericArray<Vector<T, LatSet::d>>& nUF =
      nBLat.template getField<VELOCITY<T, LatSet::d>>().getField(0);

    for (std::size_t i = 0; i < size;) {
      if constexpr (LatSet::d == 2) {
        popIntpCommunicate<0>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i                  );
        popIntpCommunicate<1>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i += SendRecvOffset);
        popIntpCommunicate<2>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i += SendRecvOffset);
        popIntpCommunicate<3>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i += SendRecvOffset);
        i += SendRecvOffset;
      } else if constexpr (LatSet::d == 3) {
        popIntpCommunicate<0>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i                  );
        popIntpCommunicate<1>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i += SendRecvOffset);
        popIntpCommunicate<2>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i += SendRecvOffset);
        popIntpCommunicate<3>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i += SendRecvOffset);
        popIntpCommunicate<4>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i += SendRecvOffset);
        popIntpCommunicate<5>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i += SendRecvOffset);
        popIntpCommunicate<6>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i += SendRecvOffset);
        popIntpCommunicate<7>(BLat, nBLat, OmegaC, comm, nRhoF, nUF, i += SendRecvOffset);
        i += SendRecvOffset;
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
template <unsigned int D>
void BlockLatticeManager<T, LatSet, TypePack>::popIntpCommunicate(BLOCKLATTICE& BLat, const BLOCKLATTICE& nBLat, T OmegaC, const SharedComm& comm, 
const GenericArray<T>& nRhoF, const GenericArray<Vector<T, LatSet::d>>& nUF, std::size_t i) {
  static constexpr std::size_t RecvOffset = LatSet::d == 2 ? 4 : 8;
  // static constexpr std::size_t SendRecvOffset = RecvOffset + 1;

  T averRho = getInterpolation<D, T, LatSet::d>(nRhoF, i, comm.SendRecvCells);
  Vector<T, LatSet::d> averU = getInterpolation<D, T, LatSet::d>(nUF, i, comm.SendRecvCells);

  std::array<T, LatSet::q> feq{};
  Equilibrium<T, LatSet>::SecondOrder(feq, averU, averRho);

  std::array<T*, LatSet::q> CellPop =
    BLat.template getField<POP<T, LatSet::q>>().getArray(comm.SendRecvCells[i+RecvOffset]);
  for (unsigned int k = 0; k < LatSet::q; ++k) {
    T averpop = getInterpolation<D, T, LatSet::d>(
      nBLat.template getField<POP<T, LatSet::q>>().getField(k), i, comm.SendRecvCells);
    // convert from coarse to fine
    *(CellPop[k]) = RefineConverter<T>::getPopF(averpop, feq[k], OmegaC);
  }
}

// DynamicBlockLatticeHelper2D

#include "utils/fdm_solver.h"

template <typename T, typename LatSet, typename TypePack>
void DynamicBlockLatticeHelper2D<T, LatSet, TypePack>::ComputeGradNorm2() {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (std::size_t icell = 0; icell < BlockGeoHelper.getBlockCells().size(); ++icell) {
    // rho field
    GenericArray<T>& GradNorm2 = _GradNorm2s[icell].getField(0);

    const BasicBlock<T, 2>& cellblock = BlockGeoHelper.getBlockCell(icell);
    T voxsize = cellblock.getVoxelSize();
    // find corresponding block
    for (int iblock = 0; iblock < BlockGeo.getBlockNum(); ++iblock) {
      const BasicBlock<T, 2>& block = BlockGeo.getBlock(iblock);
      const BasicBlock<T, 2>& baseblock = BlockGeo.getBlock(iblock).getBaseBlock();
      if (isOverlapped(cellblock, baseblock)) {
        const ScalarField<T>& RhoF =
          BlockLatMan.getBlockLat(iblock).template getField<GenericRho>();
        FDM2D<T> FDM(block.getNx(), block.getNy(), RhoF.getField(0));
        const AABB<T, 2> intsec = getIntersection(cellblock, baseblock);
        int Nx = static_cast<int>(std::round(intsec.getExtension()[0] / voxsize));
        int Ny = static_cast<int>(std::round(intsec.getExtension()[1] / voxsize));
        // cell start
        Vector<T, 2> cellstart = intsec.getMin() - cellblock.getMin();
        int cellstartx = static_cast<int>(std::round(cellstart[0] / voxsize));
        int cellstarty = static_cast<int>(std::round(cellstart[1] / voxsize));
        // block start
        Vector<T, 2> blockstart = intsec.getMin() - block.getMin();
        int blockstartx = static_cast<int>(std::round(blockstart[0] / voxsize));
        int blockstarty = static_cast<int>(std::round(blockstart[1] / voxsize));

        for (int iy = 1; iy < Ny - 1; ++iy) {
          for (int ix = 1; ix < Nx - 1; ++ix) {
            std::size_t idcell = (iy + cellstarty) * cellblock.getNx() + ix + cellstartx;
            std::size_t idblock = (iy + blockstarty) * block.getNx() + ix + blockstartx;
            GradNorm2[idcell] = FDM.gradnorm2(idblock);
          }
        }
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
void DynamicBlockLatticeHelper2D<T, LatSet, TypePack>::UpdateMaxGradNorm2() {
  for (std::size_t i = 0; i < _GradNorm2s.size(); ++i) {
    GenericArray<T>& GradNorm2A = _GradNorm2s[i].getField(0);
    // get max of gradnorm2
    T maxGradNorm2 = T(0);
    for (std::size_t id = 0; id < GradNorm2A.size(); ++id) {
      if (GradNorm2A[id] > maxGradNorm2) maxGradNorm2 = GradNorm2A[id];
    }
    _MaxGradNorm2s[i] = maxGradNorm2;
  }
}

template <typename T, typename LatSet, typename TypePack>
bool DynamicBlockLatticeHelper2D<T, LatSet, TypePack>::WillRefineOrCoarsen() {
  _GradNorm2s.clear();
  // init gradnorm2
  for (BasicBlock<T, 2>& block : BlockGeoHelper.getBlockCells()) {
    _GradNorm2s.emplace_back(block.getN(), T(0));
  }

  ComputeGradNorm2();
  UpdateMaxGradNorm2();
  // statistics
  int refineNum = 0;
  int coarsenNum = 0;
  for (int i = 0; i < _GradNorm2s.size(); ++i) {
    BasicBlock<T, 2>& block = BlockGeoHelper.getBlockCell(i);
    int level = static_cast<int>(block.getLevel());
    if (level < BlockGeoHelper.getLevelLimit()) {
      if (_MaxGradNorm2s[i] > _RefineTholds[level]) {
        block.refine();
        ++refineNum;
      }
    } else if (level > 0) {
      if (_MaxGradNorm2s[i] < _CoarsenTholds[level - 1]) {
        block.coarsen();
        ++coarsenNum;
      }
    }
  }

  if (refineNum > 0 || coarsenNum > 0) {
    std::cout << "Block Cell RefineNum: " << refineNum << " CoarsenNum: " << coarsenNum
              << std::endl;
    return true;
  } else {
    return false;
  }
}

template <typename T, typename LatSet, typename TypePack>
void DynamicBlockLatticeHelper2D<T, LatSet, TypePack>::GeoRefine(int OptProcNum,
                                                                 int MaxProcNum,
                                                                 bool enforce) {
  // post refine
  BlockGeoHelper.CheckRefine();
  // update BasicBlocks in GeoHelper
  BlockGeoHelper.CreateBlocks();
  BlockGeoHelper.AdaptiveOptimization(OptProcNum, MaxProcNum, enforce);
  BlockGeoHelper.LoadBalancing();
}


template <typename T, typename LatSet, typename TypePack>
void DynamicBlockLatticeHelper2D<T, LatSet, TypePack>::PopFieldInit() {
  // rho field data transfer could be done here
  // BlockLatMan.template getField<GenericRho>().InitAndComm(
  // BlockGeoHelper, BlockLatMan.getConverter().getLatRhoInit());
  // pop field data, this assumes that other fields like rho and velocity have been
  // transfered
  // BlockLatMan.template getField<POP<T, LatSet::q>>().Init(BlockGeoHelper);
  // now pops field data is transferred, but conversion of distribution function
  // between blocks of different refinement levels is not done yet
  for (std::size_t inewblock = 0;
       inewblock < BlockLatMan.template getField<POP<T, LatSet::q>>().getBlockFields().size();
       ++inewblock) {
    auto& PopsField =
      BlockLatMan.template getField<POP<T, LatSet::q>>().getBlockField(inewblock);
    const auto& RhoField =
      BlockLatMan.template getField<GenericRho>().getBlockField(inewblock);
    const auto& VELOCITYF =
      BlockLatMan.template getField<VELOCITY<T, LatSet::d>>().getBlockField(inewblock);

    const BasicBlock<T, 2>& newblock = BlockGeo.getBlock(inewblock);
    const BasicBlock<T, 2>& newbaseblock = BlockGeo.getBlock(inewblock).getBaseBlock();
    std::uint8_t Level = newblock.getLevel();
    // find overlapped old block field
    for (int iblock = 0; iblock < BlockGeo.getBlockNum(); ++iblock) {
      const BasicBlock<T, 2>& baseblock = BlockGeoHelper.getAllOldBasicBlock(static_cast<std::size_t>(iblock));
      int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
      const BasicBlock<T, 2>& block = baseblock.getExtBlock(overlap);
      if (isOverlapped(newbaseblock, baseblock)) {
        // get omega
        T omega = BlockLatMan.getBlockLat(iblock).getOmega();
        if (Level > block.getLevel()) {
          PopConversionCoarseToFine(RhoField, VELOCITYF, PopsField, baseblock, newblock,
                                    newbaseblock, omega);
        } else if (Level < block.getLevel()) {
          PopConversionFineToCoarse(RhoField, VELOCITYF, PopsField, baseblock, newblock,
                                    newbaseblock, omega);
        }
      }
    }
  }
  BlockLatMan.template getField<POP<T, LatSet::q>>().Communicate();
}

template <typename T, typename LatSet, typename TypePack>
void DynamicBlockLatticeHelper2D<T, LatSet, TypePack>::PopConversionFineToCoarse(
  const ScalarField<T>& RhoF, const VectorFieldAOS<T, LatSet::d>& VelocityF,
  PopulationField<T, LatSet::q>& PopsF, const BasicBlock<T, 2>& FBaseBlock,
  const BasicBlock<T, 2>& CBlock, const BasicBlock<T, 2>& CBaseBlock, T OmegaF) {
  // get intersection
  const AABB<T, 2> intsec = getIntersection(CBaseBlock, FBaseBlock);
  int Nx = intsec.getExtension()[0] / CBlock.getVoxelSize();
  Nx = Nx == 0 ? 1 : Nx;
  int Ny = intsec.getExtension()[1] / CBlock.getVoxelSize();
  Ny = Ny == 0 ? 1 : Ny;
  // get start index of intsec in CBlock
  Vector<T, 2> startC = intsec.getMin() - CBlock.getMin();
  int startx = static_cast<int>(startC[0] / CBlock.getVoxelSize());
  int starty = static_cast<int>(startC[1] / CBlock.getVoxelSize());

  const GenericArray<T>& RhoArr = RhoF.getField(0);
  const GenericArray<Vector<T, 2>>& VelocityArr = VelocityF.getField(0);

  for (int iy = starty; iy < starty + Ny; ++iy) {
    for (int ix = startx; ix < startx + Nx; ++ix) {
      std::size_t id = iy * CBlock.getNx() + ix;
      std::array<T*, LatSet::q> cell = PopsF.getArray(id);
      // get feq, peparing for pop conversion
      std::array<T, LatSet::q> feq{};
      Equilibrium<T, LatSet>::SecondOrder(feq, VelocityArr[id], RhoArr[id]);
      // convert from fine to coarse
      for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
        T popC = RefineConverter<T>::getPopC(*(cell[iArr]), feq[iArr], OmegaF);
        *(cell[iArr]) = popC;
      }
    }
  }
}


template <typename T, typename LatSet, typename TypePack>
void DynamicBlockLatticeHelper2D<T, LatSet, TypePack>::PopConversionCoarseToFine(
  const ScalarField<T>& RhoF, const VectorFieldAOS<T, LatSet::d>& VelocityF,
  PopulationField<T, LatSet::q>& PopsF, const BasicBlock<T, 2>& CBaseBlock,
  const BasicBlock<T, 2>& FBlock, const BasicBlock<T, 2>& FBaseBlock, T OmegaC) {
  // get intersection
  const AABB<T, 2> intsec = getIntersection(CBaseBlock, FBaseBlock);
  int Nx = intsec.getExtension()[0] / FBlock.getVoxelSize();
  Nx = Nx == 0 ? 1 : Nx;
  int Ny = intsec.getExtension()[1] / FBlock.getVoxelSize();
  Ny = Ny == 0 ? 1 : Ny;
  // get start index of intsec in FBlock
  Vector<T, 2> startF = intsec.getMin() - FBlock.getMin();
  int startx = static_cast<int>(startF[0] / FBlock.getVoxelSize());
  int starty = static_cast<int>(startF[1] / FBlock.getVoxelSize());

  const GenericArray<T>& RhoArr = RhoF.getField(0);
  const GenericArray<Vector<T, 2>>& VelocityArr = VelocityF.getField(0);

  for (int iy = starty; iy < starty + Ny; ++iy) {
    for (int ix = startx; ix < startx + Nx; ++ix) {
      std::size_t id = iy * FBlock.getNx() + ix;
      std::array<T*, LatSet::q> cell = PopsF.getArray(id);
      // get feq, peparing for pop conversion
      std::array<T, LatSet::q> feq{};
      Equilibrium<T, LatSet>::SecondOrder(feq, VelocityArr[id], RhoArr[id]);
      // convert from fine to coarse
      for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
        T popF = RefineConverter<T>::getPopF(*(cell[iArr]), feq[iArr], OmegaC);
        *(cell[iArr]) = popF;
      }
    }
  }
}


// template <typename T, typename LatSet>
// void DynamicBlockLatticeHelper2D<T, LatSet, TypePack>::UpdateGradNorm2F() {
//   // min centre of the whole domain
//   Vector<T, 2> GMinCentre = BlockGeo.getBaseBlock().getMinCenter();
//   int GNx = BlockGeo.getBaseBlock().getNx();
//   GenericArray<T>& GradNorm2F = _GradNorm2F.getField(0);
// #pragma omp parallel for num_threads(Thread_Num) schedule(static)
//   for (int i = 0; i < BlockGeoHelper.getBlockCells().size(); ++i) {
//     // rho field
//     GenericArray<T>& GradNorm2 = _GradNorm2s[i].getField(0);
//     BasicBlock<T, 2>& block = BlockGeoHelper.getBlockCell(i);
//     Vector<T, 2> minCentre = block.getMinCenter();
//     // get ext
//     Vector<T, 2> Ext = minCentre - GMinCentre;
//     int x = static_cast<int>(std::round(Ext[0] / block.getVoxelSize()));
//     int y = static_cast<int>(std::round(Ext[1] / block.getVoxelSize()));
//     // copy data from _GradNorm2s[i] to _GradNorm2F
//     std::size_t id = 0;
//     for (int iy = y; iy < y + block.getNy(); ++iy) {
//       // for (int ix = x; ix < x + block.getNx(); ++ix) {
//       std::size_t gidstart = GNx * iy + x;
//       std::copy(GradNorm2.getdataPtr(id), GradNorm2.getdataPtr(id + block.getNx()),
//                 GradNorm2F.getdataPtr(gidstart));
//       id += block.getNx();
//     }
//   }
// }
