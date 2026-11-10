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

#include <algorithm>
#include <utility>

#include "data_struct/block_lattice.h"
// Equilibrium / equilibrium::SecondOrder used below (MPI refinement path)
#include "lbm/equilibrium.h"
// shared plan post/wait/unpack runner (lattice, field and hidden steps)
#include "parallel/plan_runner.h"

#ifdef MPI_ENABLED
template <typename T, typename LatSet, typename TypePack>
void BlockLattice<T, LatSet, TypePack>::initHiddenCommunication() {
    if constexpr (LatSet::d == 2) {
      // add to _collideFirst: the inner ring adjacent to the overlap
      std::array<int, 2> bdj{this->getOverlap(), this->getNy() - this->getOverlap()-1};
      std::array<int, 2> bdi{this->getOverlap(), this->getNx() - this->getOverlap()-1};
      for (int j = this->getOverlap(); j < this->getNy() - this->getOverlap(); ++j) {
        for (int i = this->getOverlap(); i < this->getNx() - this->getOverlap(); ++i) {
          std::size_t id = j * this->getNx() + i;
          if (j == bdj[0] || j == bdj[1] || i == bdi[0] || i == bdi[1]) {
            _collideFirst.push_back(id);
          }
        }
      }
  } else if constexpr (LatSet::d == 3) {
    // add to _collideFirst: the inner ring adjacent to the overlap
    std::array<int, 2> bdk{this->getOverlap(), this->getNz() - this->getOverlap()-1};
    std::array<int, 2> bdj{this->getOverlap(), this->getNy() - this->getOverlap()-1};
    std::array<int, 2> bdi{this->getOverlap(), this->getNx() - this->getOverlap()-1};
    for (int k = this->getOverlap(); k < this->getNz() - this->getOverlap(); ++k) {
      for (int j = this->getOverlap(); j < this->getNy() - this->getOverlap(); ++j) {
        for (int i = this->getOverlap(); i < this->getNx() - this->getOverlap(); ++i) {
          std::size_t id = k * this->getProjection()[2] + j * this->getProjection()[1] + i;
          if (k == bdk[0] || k == bdk[1] || j == bdj[0] || j == bdj[1] || i == bdi[0] || i == bdi[1]) {
            _collideFirst.push_back(id);
          }
        }
      }
    }
  }

  // add to _collideLast: every cell that is not in the inner ring (the overlap
  // included).
  //
  // The overlap cells next to a neighbour are refreshed by the communicate
  // below, but the communicate only writes the directions the neighbour needs;
  // the other directions keep whatever is in the cell.  The plain path has just
  // collided them (it collides every cell), so these must be collided here too,
  // otherwise the two paths differ in the overlap, and that difference reaches
  // the fluid through the array rotation and the boundary conditions'
  // getPrevious() (which reads a neighbour cell).  In other words: both paths
  // must collide exactly the same set of cells; the only difference may be
  // whether the collision or the halo refresh comes first.
  std::vector<bool> isFirst(this->getN(), false);
  for (std::size_t id : _collideFirst) isFirst[id] = true;
  for (std::size_t id = 0; id < this->getN(); ++id) {
    if (!isFirst[id]) _collideLast.push_back(id);
  }
}
#endif

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
  initHiddenCommunication();
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
  if (map.isInitialized()) {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
    for (std::size_t id = 0; id < voxNum; ++id) {
      const std::size_t cellId = map[id];
      cell.setId(cellId);
      CELLDYNAMICS::Execute(flagarr[cellId], cell);
    }
    return;
  }
#endif
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
  for (std::size_t id = 0; id < voxNum; ++id) {
    cell.setId(id);
    CELLDYNAMICS::Execute(flagarr[id], cell);
  }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS>
void BlockLattice<T, LatSet, TypePack>::ApplyCellDynamics() {
  Cell<T, LatSet, TypePack> cell(0, *this);
  const std::size_t voxNum = this->getVoxNum();
#ifdef _VOX_ENABLED
  const VoxelMap& map = this->getBlock().getVoxelMap();
  if (map.isInitialized()) {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
    for (std::size_t id = 0; id < voxNum; ++id) {
      cell.setId(map[id]);
      CELLDYNAMICS::apply(cell);
    }
    return;
  }
#endif
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
  for (std::size_t id = 0; id < voxNum; ++id) {
    cell.setId(id);
    CELLDYNAMICS::apply(cell);
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

#ifdef MPI_ENABLED
// First half of the hidden communication: collide the cells adjacent to the
// overlap, then post this block's non-blocking MPI sends and receives.
//
// Only the cells adjacent to the overlap are collided here, and only the
// population directions the neighbour actually needs are exchanged -- that is
// what hides the MPI latency behind the collision of the interior cells.
template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename ArrayType>
void BlockLattice<T, LatSet, TypePack>::ApplyCellDynamicsWithHiddenCommFirst(
  const ArrayType& flagarr) {
  Cell<T, LatSet, TypePack> cell(0, *this);
#ifdef _HIDDEN_PROFILE
  const double _t0 = MPI_Wtime();
#endif
  // apply first cell dynamics
  for (std::size_t id : _collideFirst) {
    cell.setId(id);
    CELLDYNAMICS::Execute(flagarr[id], cell);
  }
#ifdef _HIDDEN_PROFILE
  g_hiddenCommProfile.collideFirst += MPI_Wtime() - _t0;
#endif
}

// Pack one communicator into the manager's per-rank send buffer.  The manager
// computes the offset and guarantees the buffer outlives the MPI_Isend.
template <typename T, typename LatSet, typename TypePack>
template <typename buffer_type>
void BlockLattice<T, LatSet, TypePack>::packNormalSend(
  std::size_t commIdx, buffer_type& buffer, std::size_t offset) {
  const std::vector<DistributedComm>& sends = this->getBlock().getCommunicator().DirSends;
  const std::vector<std::size_t>& cells = sends[commIdx].Cells;
  for (unsigned int iArr : SendDirection[commIdx]) {
    const auto& Array = this->template getField<POP<T, LatSet::q>>().getField(iArr);
    for (std::size_t id : cells) {
      buffer[offset] = Array[id];
      ++offset;
    }
  }
}

// Unpack one communicator from the manager's per-rank receive buffer.
template <typename T, typename LatSet, typename TypePack>
template <typename buffer_type>
void BlockLattice<T, LatSet, TypePack>::unpackNormalRecv(
  std::size_t commIdx, const buffer_type& buffer, std::size_t offset) {
  const std::vector<DistributedComm>& recvs = this->getBlock().getCommunicator().DirRecvs;
  const std::vector<std::size_t>& cells = recvs[commIdx].Cells;
  for (unsigned int iArr : RecvDirection[commIdx]) {
    auto& Array = this->template getField<POP<T, LatSet::q>>().getField(iArr);
    for (std::size_t id : cells) {
      Array[id] = buffer[offset];
      ++offset;
    }
  }
}

// Pack one full-direction communicator (all q populations) into the manager's
// per-rank send buffer.  all selects the wider AllMPIComm cell set.
template <typename T, typename LatSet, typename TypePack>
template <typename buffer_type>
void BlockLattice<T, LatSet, TypePack>::packFullSend(
  std::size_t commIdx, buffer_type& buffer, std::size_t offset, bool all) {
  const std::vector<DistributedComm>& sends =
    all ? this->getBlock().getCommunicator().AllMPIComm.Sends
        : this->getBlock().getCommunicator().MPIComm.Sends;
  const std::vector<std::size_t>& cells = sends[commIdx].Cells;
  for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
    const auto& Array = this->template getField<POP<T, LatSet::q>>().getField(iArr);
    for (std::size_t id : cells) {
      buffer[offset] = Array[id];
      ++offset;
    }
  }
}

// Unpack one full-direction communicator from the manager's per-rank receive
// buffer.
template <typename T, typename LatSet, typename TypePack>
template <typename buffer_type>
void BlockLattice<T, LatSet, TypePack>::unpackFullRecv(
  std::size_t commIdx, const buffer_type& buffer, std::size_t offset, bool all) {
  const std::vector<DistributedComm>& recvs =
    all ? this->getBlock().getCommunicator().AllMPIComm.Recvs
        : this->getBlock().getCommunicator().MPIComm.Recvs;
  const std::vector<std::size_t>& cells = recvs[commIdx].Cells;
  for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
    auto& Array = this->template getField<POP<T, LatSet::q>>().getField(iArr);
    for (std::size_t id : cells) {
      Array[id] = buffer[offset];
      ++offset;
    }
  }
}

// Second half: collide the remaining interior cells.  The overlap is still
// stale here -- exactly as in a normal communicate, where every cell is collided
// before the halos are refreshed -- so the result does not depend on where the
// halo refresh happens inside the step.
template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename ArrayType>
void BlockLattice<T, LatSet, TypePack>::ApplyCellDynamicsWithHiddenCommLast(
  const ArrayType& flagarr) {
  Cell<T, LatSet, TypePack> cell(0, *this);
#ifdef _HIDDEN_PROFILE
  const double _t0 = MPI_Wtime();
#endif
  for (std::size_t id : _collideLast) {
    cell.setId(id);
    CELLDYNAMICS::Execute(flagarr[id], cell);
  }
#ifdef _HIDDEN_PROFILE
  g_hiddenCommProfile.collideLast += MPI_Wtime() - _t0;
#endif
}
#endif

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
template <typename CELLDYNAMICS, typename CELLTYPE, typename ArrayType>
void BlockLattice<T, LatSet, TypePack>::CuDevApplyCellDynamics(ArrayType& flagarr) {
  // the launch geometry is a property of the cell's POP storage strategy
  const unsigned int blockSize = CELLTYPE::block_size;
  const unsigned int blockNum = (this->getN() + blockSize - 1) / blockSize;
  CuDevApplyCellDynamicsKernel<T, LatSet, cudev_TypePack, CELLTYPE, CELLDYNAMICS, typename ArrayType::cudev_array_type><<<blockNum, blockSize>>>(dev_BlockLat, flagarr.get_devObj(), this->getN());
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename CELLTYPE>
void BlockLattice<T, LatSet, TypePack>::CuDevApplyCellDynamics() {
  // the launch geometry is a property of the cell's POP storage strategy
  const unsigned int blockSize = CELLTYPE::block_size;
  const unsigned int blockNum = (this->getN() + blockSize - 1) / blockSize;
  CuDevApplyCellDynamicsKernel<T, LatSet, cudev_TypePack, CELLTYPE, CELLDYNAMICS><<<blockNum, blockSize>>>(dev_BlockLat, this->getN());
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
template <typename buffer_type>
void BlockLattice<T, LatSet, TypePack>::packAverSend(
  const DistributedComm& comm, buffer_type& buffer, std::size_t offset) {
  static constexpr std::size_t SendSize = flbcomm::AverSrcNum(LatSet::d);

  const T OmegaF = this->getOmega();
  const GenericArray<T>& RhoArr = this->template getField<GenericRho>().getField(0);
  const GenericArray<Vector<T, LatSet::d>>& UArr = this->template getField<VELOCITY<T, LatSet::d>>().getField(0);

  const std::vector<std::size_t>& sends = comm.Cells;
  const std::size_t size = sends.size();
  std::size_t bufidx = offset;

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
}

template <typename T, typename LatSet, typename TypePack>
template <typename buffer_type>
void BlockLattice<T, LatSet, TypePack>::packIntpSend(
  const DistributedComm& comm, buffer_type& buffer, std::size_t offset) {
  static constexpr std::size_t SendSize = flbcomm::AverSrcNum(LatSet::d);

  const T OmegaC = this->getOmega();
  const GenericArray<T>& RhoArr = this->template getField<GenericRho>().getField(0);
  const GenericArray<Vector<T, LatSet::d>>& UArr = this->template getField<VELOCITY<T, LatSet::d>>().getField(0);

  const std::vector<std::size_t>& sends = comm.Cells;
  const std::size_t size = sends.size();
  std::size_t bufidx = offset;
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
}

template <typename T, typename LatSet, typename TypePack>
template <typename buffer_type>
void BlockLattice<T, LatSet, TypePack>::unpackCrossRecv(
  const DistributedComm& comm, const buffer_type& buffer, std::size_t offset) {
  std::size_t bufidx = offset;
  for (std::size_t id : comm.Cells) {
    flbcomm::checkCellIndex(id, this->getN(), this->getBlock().getBlockId(),
                            comm.TargetBlockId, "unpackCrossRecv");
    std::array<T*, LatSet::q> CellPop = this->template getField<POP<T, LatSet::q>>().getArray(id);
    for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
      *(CellPop[iArr]) = buffer[bufidx];
      ++bufidx;
    }
  }
}


template <typename T, typename LatSet, typename TypePack>
template <unsigned int D, typename buffer_type>
void BlockLattice<T, LatSet, TypePack>::popIntp(T OmegaC, const std::vector<std::size_t>& sends, 
const GenericArray<T>& RhoArr, const GenericArray<Vector<T, LatSet::d>>& UArr, std::size_t i, 
buffer_type& buffer, std::size_t& bufidx) {
  
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
#ifdef MPI_ENABLED
  // The coalesced-communication plans refer to blocks and communicators that
  // are about to be rebuilt.  Every index of every depth goes, not just the
  // unfiltered one: a stale filtered plan reads block slots that no longer
  // hold what it was built against.
  _dirPlans.clear();
  _fullPlans.clear();
  _allPlans.clear();
  _averPlans.clear();
  _intpPlans.clear();
#endif
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
#ifdef MPI_ENABLED
  // Same as the other overload: every depth, not just the unfiltered plan.
  _dirPlans.clear();
  _fullPlans.clear();
  _allPlans.clear();
  _averPlans.clear();
  _intpPlans.clear();
#endif
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
template <typename CELLDYNAMICS, typename CELLTYPE, typename FieldType>
void BlockLatticeManager<T, LatSet, TypePack>::CuDevApplyCellDynamics(BlockFieldManager<FieldType, T, LatSet::d>& BFM){
  BlockLats[0].template CuDevApplyCellDynamics<CELLDYNAMICS, CELLTYPE, typename FieldType::array_type>(
          BFM.getBlockField(0).getField(0));
  // for (std::size_t i = 0; i < BlockLats.size(); ++i) {
  //     BlockLats[i].template CuDevApplyCellDynamics<CELLDYNAMICS, CELLTYPE, typename FieldType::array_type>(
  //         BFM.getBlockField(i).getField(0));
  // }
}

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, typename CELLTYPE>
void BlockLatticeManager<T, LatSet, TypePack>::CuDevApplyCellDynamics(){
  BlockLats[0].template CuDevApplyCellDynamics<CELLDYNAMICS, CELLTYPE>();
  // for (std::size_t i = 0; i < BlockLats.size(); ++i) {
  //   BlockLats[i].template CuDevApplyCellDynamics<CELLDYNAMICS, CELLTYPE>();
  // }
}

#endif

#ifdef MPI_ENABLED

template <typename T, typename LatSet, typename TypePack>
flbcomm::CommPlan<typename BlockLatticeManager<T, LatSet, TypePack>::buffer_type>&
BlockLatticeManager<T, LatSet, TypePack>::planFor(flbcomm::CommMode Mode,
                                                  std::size_t depth) {
  std::vector<flbcomm::CommPlan<buffer_type>>& plans =
    Mode == flbcomm::CommMode::Partial ? _dirPlans
    : Mode == flbcomm::CommMode::Full  ? _fullPlans
                                       : _allPlans;
  const std::size_t maxDepth = static_cast<std::size_t>(this->getMaxLevel());
  if (depth > maxDepth) depth = maxDepth;
  if (plans.size() <= depth) plans.resize(depth + 1);
  flbcomm::CommPlan<buffer_type>& plan = plans[depth];
  if (!plan.Ready) {
    if (Mode == flbcomm::CommMode::Partial) {
      initDirPlan(depth);
    } else {
      initFullPlan(Mode == flbcomm::CommMode::All, depth);
    }
  }
  checkPlanGeneration(plan);
  return plan;
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPIAverageCommunicate() {
  runCrossPlan(flbcomm::CrossMode::Aver,
               static_cast<std::size_t>(getMaxLevel()));
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPIInterpolateCommunicate() {
  runCrossPlan(flbcomm::CrossMode::Intp,
               static_cast<std::size_t>(getMaxLevel()));
}


template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPIAverageCommunicate(std::int64_t count) {
  runCrossPlan(flbcomm::CrossMode::Aver,
               flbcomm::activeDepth(count, static_cast<std::size_t>(getMaxLevel())));
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::MPIInterpolateCommunicate(std::int64_t count) {
  runCrossPlan(flbcomm::CrossMode::Intp,
               flbcomm::activeDepth(count, static_cast<std::size_t>(getMaxLevel())));
}

#endif

#ifdef MPI_ENABLED
// Build the coalesced plan for the plain Partial exchange: one message per peer,
// refs sorted by (target rank, Order).  Dir links carry an Order derived from
// the receiver's block ids and link kind/index; both ends derive it from the
// same receiver-based enumeration.  The group tag is a constant per exchange,
// so it no longer depends on a global link counter.
template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::initDirPlan(std::size_t depth) {
  flbcomm::CommPlan<buffer_type>& plan = _dirPlans[depth];
  flbcomm::buildPlan<buffer_type>(
    BlockLats.size(), static_cast<std::size_t>(getMaxLevel()), depth,
    this->BlockGeo.getGeometryGeneration(),
    [&](std::size_t bi) { return BlockLats[bi].getLevel(); },
    [&](std::size_t bi, bool send, std::vector<flbcomm::CommRef>& refs) {
      const std::vector<DistributedComm>& comms =
        send ? BlockLats[bi].getBlock().getCommunicator().DirSends
             : BlockLats[bi].getBlock().getCommunicator().DirRecvs;
      for (std::size_t ci = 0; ci < comms.size(); ++ci) {
        refs.push_back(flbcomm::CommRef{
          bi, ci, 0, comms[ci].TargetRank, flbcomm::DirTag, comms[ci].Order,
          send ? BlockLats[bi].normalSendSize(ci)
               : BlockLats[bi].normalRecvSize(ci)});
      }
    },
    plan);
}

// Build the coalesced plan for NormalFullCommunicate() (all == false) or
// NormalAllCommunicate() (all == true).  One message per peer; the refs are
// keyed by (receiver block, sender block), which is the same pair on both
// sides and unique per link.  The group tag is a constant per exchange (FullTag
// or AllTag), so it does not depend on block ids or a link counter.
template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::initFullPlan(bool all, std::size_t depth) {
  flbcomm::CommPlan<buffer_type>& plan = all ? _allPlans[depth] : _fullPlans[depth];
  flbcomm::buildPlan<buffer_type>(
    BlockLats.size(), static_cast<std::size_t>(getMaxLevel()), depth,
    this->BlockGeo.getGeometryGeneration(),
    [&](std::size_t bi) { return BlockLats[bi].getLevel(); },
    [&](std::size_t bi, bool send, std::vector<flbcomm::CommRef>& refs) {
      const std::vector<DistributedComm>& comms =
        all ? (send ? BlockLats[bi].getBlock().getCommunicator().AllMPIComm.Sends
                    : BlockLats[bi].getBlock().getCommunicator().AllMPIComm.Recvs)
            : (send ? BlockLats[bi].getBlock().getCommunicator().MPIComm.Sends
                    : BlockLats[bi].getBlock().getCommunicator().MPIComm.Recvs);
      const int self = BlockLats[bi].getBlock().getBlockId();
      const int tag = all ? flbcomm::AllTag : flbcomm::FullTag;
      for (std::size_t ci = 0; ci < comms.size(); ++ci) {
        const int peerBlock = comms[ci].TargetBlockId;
        const int receiver = send ? peerBlock : self;
        const int sender = send ? self : peerBlock;
        refs.push_back(flbcomm::CommRef{
          bi, ci, 0, comms[ci].TargetRank, tag,
          flbcomm::LinkKey{receiver, sender, 0, 0},
          send ? BlockLats[bi].fullSendSize(ci, all)
               : BlockLats[bi].fullRecvSize(ci, all)});
      }
    },
    plan);
}

// Cross-level plans.  One message per peer; the refs are keyed by the
// (receiver block, sender block) pair and each exchange uses a fixed tag.  The
// send and receive sides are filtered independently because Aver sends from a
// block one level finer than its neighbour and Intp the other way round.
template <typename T, typename LatSet, typename TypePack>
const std::vector<DistributedComm>&
BlockLatticeManager<T, LatSet, TypePack>::crossComms(flbcomm::CrossMode mode,
                                                     std::size_t block, bool send) {
  const auto& comm = BlockLats[block].getBlock().getCommunicator();
  if (flbcomm::crossIsAver(mode)) {
    return send ? comm.MPIComm.AverSends : comm.MPIComm.AverRecvs;
  }
  return send ? comm.MPIComm.IntpSends : comm.MPIComm.IntpRecvs;
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::initCrossPlan(flbcomm::CrossMode mode,
                                                            std::size_t depth) {
  const bool aver = flbcomm::crossIsAver(mode);
  std::vector<flbcomm::CommPlan<buffer_type>>& plans = aver ? _averPlans : _intpPlans;
  flbcomm::CommPlan<buffer_type>& plan = plans[depth];
  const std::size_t maxLevel = static_cast<std::size_t>(getMaxLevel());
  const int tag = flbcomm::crossTag(mode);
  flbcomm::buildPlan<buffer_type>(
    BlockLats.size(),
    this->BlockGeo.getGeometryGeneration(),
    [&](std::size_t bi, bool send) {
      const int deLevel = static_cast<int>(maxLevel - BlockLats[bi].getLevel());
      if (aver) {
        return send ? (deLevel + 1 <= static_cast<int>(depth) &&
                       BlockLats[bi].getLevel() != 0)
                    : (deLevel <= static_cast<int>(depth));
      }
      return send ? (deLevel - 1 <= static_cast<int>(depth) && deLevel != 0)
                  : (deLevel <= static_cast<int>(depth));
    },
    [&](std::size_t bi, bool send, std::vector<flbcomm::CommRef>& refs) {
      const std::vector<DistributedComm>& comms = crossComms(mode, bi, send);
      const int self = BlockLats[bi].getBlock().getBlockId();
      for (std::size_t ci = 0; ci < comms.size(); ++ci) {
        const int peerBlock = comms[ci].TargetBlockId;
        const int receiver = send ? peerBlock : self;
        const int sender = send ? self : peerBlock;
        refs.push_back(flbcomm::CommRef{
          bi, ci, 0, comms[ci].TargetRank, tag,
          flbcomm::LinkKey{receiver, sender, 0, 0},
          send ? BlockLats[bi].crossSendSize(comms[ci])
               : BlockLats[bi].crossRecvSize(comms[ci])});
      }
    },
    plan);
}

template <typename T, typename LatSet, typename TypePack>
flbcomm::CommPlan<typename BlockLatticeManager<T, LatSet, TypePack>::buffer_type>&
BlockLatticeManager<T, LatSet, TypePack>::crossPlanFor(flbcomm::CrossMode mode,
                                                       std::size_t depth) {
  const bool aver = flbcomm::crossIsAver(mode);
  std::vector<flbcomm::CommPlan<buffer_type>>& plans = aver ? _averPlans : _intpPlans;
  const std::size_t maxDepth = static_cast<std::size_t>(getMaxLevel());
  if (depth > maxDepth) depth = maxDepth;
  if (plans.size() <= depth) plans.resize(depth + 1);
  flbcomm::CommPlan<buffer_type>& plan = plans[depth];
  if (!plan.Ready) initCrossPlan(mode, depth);
  checkPlanGeneration(plan);
  return plan;
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::runCrossPlan(flbcomm::CrossMode mode,
                                                            std::size_t depth) {
  flbcomm::CommPlan<buffer_type>& plan = crossPlanFor(mode, depth);
  const bool aver = flbcomm::crossIsAver(mode);
  flbcomm::runPlan(
    plan,
    [&](const flbcomm::CommRef& r, buffer_type& buffer, std::size_t offset) {
      const auto& comms = crossComms(mode, r.Block, true);
      if (aver) {
        BlockLats[r.Block].packAverSend(comms[r.Comm], buffer, offset);
      } else {
        BlockLats[r.Block].packIntpSend(comms[r.Comm], buffer, offset);
      }
    },
    [&](const flbcomm::CommRef& r, const buffer_type& buffer, std::size_t offset) {
      const auto& comms = crossComms(mode, r.Block, false);
      BlockLats[r.Block].unpackCrossRecv(comms[r.Comm], buffer, offset);
    });
}

// Post a plan: pack every ref into its peer's buffer, one message per peer,
// then wait and unpack.  Mode only selects the packers -- the schedule is the
// same for the three exchanges, which differ in which cells and populations
// they carry.
template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::postPlanExchange(
  flbcomm::CommPlan<buffer_type>& plan, flbcomm::CommMode Mode,
  std::vector<MPI_Request>& sendRequests,
  std::vector<MPI_Request>& recvRequests) {
  flbcomm::postPlan(
    plan, sendRequests, recvRequests,
    [&](const flbcomm::CommRef& r, buffer_type& buffer, std::size_t offset) {
      if (Mode == flbcomm::CommMode::Partial) {
        BlockLats[r.Block].packNormalSend(r.Comm, buffer, offset);
      } else {
        BlockLats[r.Block].packFullSend(r.Comm, buffer, offset,
                                        Mode == flbcomm::CommMode::All);
      }
    });
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::finishPlanExchange(
  flbcomm::CommPlan<buffer_type>& plan, flbcomm::CommMode Mode,
  std::vector<MPI_Request>& recvRequests) {
  flbcomm::finishRecvs(
    plan, recvRequests,
    [&](const flbcomm::CommRef& r, const buffer_type& buffer, std::size_t offset) {
      if (Mode == flbcomm::CommMode::Partial) {
        BlockLats[r.Block].unpackNormalRecv(r.Comm, buffer, offset);
      } else {
        BlockLats[r.Block].unpackFullRecv(r.Comm, buffer, offset,
                                          Mode == flbcomm::CommMode::All);
      }
    });
}

#endif

template <typename T, typename LatSet, typename TypePack>
template <typename CELLDYNAMICS, flbcomm::CommMode Mode, typename FieldType>
void BlockLatticeManager<T, LatSet, TypePack>::ApplyCellDynamicsWithHiddenComm(
  const BlockFieldManager<FieldType, T, LatSet::d>& BFM) {

  #ifdef MPI_ENABLED
  // Partial + refined blocks is not a supported combination: the direction
  // filter and the cross-level halos are only validated for same-level 2D
  // blocks.  The plain NormalCommunicate(count) path stays available for the
  // amr_comm diagnostic; this guard only closes the hidden path, which has no
  // count overload and therefore cannot sub-cycle.
  if constexpr (Mode == flbcomm::CommMode::Partial) {
    if (this->getMaxLevel() > 0) {
      std::cerr << "[hidden comm] Partial is not supported on refined blocks; "
                   "use CommMode::Full or CommMode::All.\n";
      std::exit(1);
    }
  }
  // Pick the plan that matches the requested exchange.  The three plans differ
  // only in the cell lists and message sizes; the schedule below is the same.
  // The hidden step is never level-filtered, so it always uses the deepest
  // plan, the one the unfiltered overloads use too.
  flbcomm::CommPlan<buffer_type>& planRef =
    planFor(Mode, static_cast<std::size_t>(this->getMaxLevel()));
  constexpr bool All = (Mode == flbcomm::CommMode::All);

  // --- 1. collide the cells adjacent to the overlap ---
  // Every block does this before any halo is touched, so the values packed
  // below are always post-collision values of the neighbour's inner cells.
#if defined(_OPENMP) && !defined(_HIDDEN_PROFILE) && !defined(SingleBlock_OMP)
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
  for (int i = 0; i < static_cast<int>(BlockLats.size()); ++i) {
    BlockLats[i].template ApplyCellDynamicsWithHiddenCommFirst<CELLDYNAMICS,
      typename FieldType::array_type>(BFM.getBlockField(i).getField(0));
  }

  // --- 1b. post receives first, then pack/post this rank's sends ---
  // postSends() packs in parallel and posts the sends serially; the plan owns
  // the buffers until the wait below.
  std::vector<MPI_Request> SendRequests;
  std::vector<MPI_Request> RecvRequests;
#ifdef _HIDDEN_PROFILE
  const double _r0 = MPI_Wtime();
#endif
  flbcomm::postRecvs(planRef, RecvRequests);
#ifdef _HIDDEN_PROFILE
  g_hiddenCommProfile.recv += MPI_Wtime() - _r0;
  const double _p0 = MPI_Wtime();
#endif
  flbcomm::postSends(
    planRef, SendRequests,
    [&](const flbcomm::CommRef& r, buffer_type& buffer, std::size_t offset) {
      if constexpr (Mode == flbcomm::CommMode::Partial) {
        BlockLats[r.Block].packNormalSend(r.Comm, buffer, offset);
      } else {
        BlockLats[r.Block].packFullSend(r.Comm, buffer, offset, All);
      }
    });
#ifdef _HIDDEN_PROFILE
  g_hiddenCommProfile.send += MPI_Wtime() - _p0;
#endif

  // --- 2. collide the remaining interior cells ---
  // The overlap is still stale at this point, exactly as with a normal
  // communicate (all cells are collided before the halos are refreshed).
#if defined(_OPENMP) && !defined(_HIDDEN_PROFILE) && !defined(SingleBlock_OMP)
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
  for (int i = 0; i < static_cast<int>(BlockLats.size()); ++i) {
    BlockLats[i].template ApplyCellDynamicsWithHiddenCommLast<CELLDYNAMICS,
      typename FieldType::array_type>(BFM.getBlockField(i).getField(0));
  }

#ifdef _HIDDEN_PROFILE
  const double _w0 = MPI_Wtime();
#endif
  flbcomm::waitSends(SendRequests);
#ifdef _HIDDEN_PROFILE
  g_hiddenCommProfile.waitall += MPI_Wtime() - _w0;
  g_hiddenCommProfile.sendRequests += static_cast<long>(SendRequests.size());
  ++g_hiddenCommProfile.calls;
  const double _h0 = MPI_Wtime();
#endif

  // --- 3. refresh the halos ---
  // 3a. neighbours on THIS rank: shared-memory copy.  This was missing before,
  //     which silently dropped every same-rank halo: a block whose neighbour
  //     lives on the same rank never had its overlap refreshed, and the error
  //     spread inward from there.  It went unnoticed because the configurations
  //     exercised so far had exactly one block per rank.
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (int i = 0; i < static_cast<int>(BlockLats.size()); ++i) {
    if constexpr (Mode == flbcomm::CommMode::Partial) {
      normalCommunicate(BlockLats[i]);
    } else if constexpr (Mode == flbcomm::CommMode::Full) {
      normalFullCommunicate(BlockLats[i]);
    } else {
      normalAllCommunicate(BlockLats[i]);
    }
  }
#ifdef _HIDDEN_PROFILE
  g_hiddenCommProfile.shared += MPI_Wtime() - _h0;
  const double _s0 = MPI_Wtime();
#endif

  // 3b. neighbours on other ranks: wait for the coalesced messages and unpack
  flbcomm::finishRecvs(
    planRef, RecvRequests,
    [&](const flbcomm::CommRef& r, const buffer_type& buffer, std::size_t offset) {
      if constexpr (Mode == flbcomm::CommMode::Partial) {
        BlockLats[r.Block].unpackNormalRecv(r.Comm, buffer, offset);
      } else {
        BlockLats[r.Block].unpackFullRecv(r.Comm, buffer, offset, All);
      }
    });
#ifdef _HIDDEN_PROFILE
  g_hiddenCommProfile.set += MPI_Wtime() - _s0;
#endif

  #else

  ApplyCellDynamics<CELLDYNAMICS, FieldType>(BFM);

  #endif

}

template <typename T, typename LatSet, typename TypePack>
template <flbcomm::CommMode Mode>
void BlockLatticeManager<T, LatSet, TypePack>::normalCommunicateOverlapped(
  std::size_t depth) {
#ifdef MPI_ENABLED
  flbcomm::CommPlan<buffer_type>& plan = planFor(Mode, depth);
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> recvRequests;
  postPlanExchange(plan, Mode, sendRequests, recvRequests);
#endif
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
  for (int i = 0; i < static_cast<int>(BlockLats.size()); ++i) {
    BLOCKLATTICE& BLat = BlockLats[i];
    const std::size_t deLevel =
      static_cast<std::size_t>(getMaxLevel() - BLat.getLevel());
    if (deLevel > depth) continue;
    if constexpr (Mode == flbcomm::CommMode::Partial) {
      normalCommunicate(BLat);
    } else if constexpr (Mode == flbcomm::CommMode::Full) {
      normalFullCommunicate(BLat);
    } else {
      normalAllCommunicate(BLat);
    }
  }
#ifdef MPI_ENABLED
  flbcomm::waitSends(sendRequests);
  finishPlanExchange(plan, Mode, recvRequests);
#endif
}

// Although the new communication scheme of populations is efficient, it did not communicate all directions,
// when you use other post-stream processors(after stream and before communication) in cells to be communicated,
// if they might use invalid pop in one direction to get pop in anohter direction which will not be communicated
// e.g. the half-way bounce-back boundary condition,
// then you should use Lattice.getField<POP<T, LatSet::q>>().Communicate() instead
template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalCommunicate() {
  normalCommunicateOverlapped<flbcomm::CommMode::Partial>(
    static_cast<std::size_t>(getMaxLevel()));
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalCommunicate(std::int64_t count) {
  normalCommunicateOverlapped<flbcomm::CommMode::Partial>(
    flbcomm::activeDepth(count, static_cast<std::size_t>(getMaxLevel())));
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalFullCommunicate() {
  normalCommunicateOverlapped<flbcomm::CommMode::Full>(
    static_cast<std::size_t>(getMaxLevel()));
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalFullCommunicate(std::int64_t count) {
  normalCommunicateOverlapped<flbcomm::CommMode::Full>(
    flbcomm::activeDepth(count, static_cast<std::size_t>(getMaxLevel())));
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalAllCommunicate() {
  normalCommunicateOverlapped<flbcomm::CommMode::All>(
    static_cast<std::size_t>(getMaxLevel()));
}

template <typename T, typename LatSet, typename TypePack>
void BlockLatticeManager<T, LatSet, TypePack>::NormalAllCommunicate(std::int64_t count) {
  normalCommunicateOverlapped<flbcomm::CommMode::All>(
    flbcomm::activeDepth(count, static_cast<std::size_t>(getMaxLevel())));
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
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(LatSet::d);
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

      // The recv slot addresses the cell this block writes; the geometry built
      // it and the write below cannot be bounds-checked after the fact.
      flbcomm::checkCellIndex(comm.SendRecvCells[i + RecvOffset], BLat.getN(),
                              BLat.getBlock().getBlockId(), comm.SendBlockId,
                              "averCommunicate");
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
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(LatSet::d);
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
  static constexpr std::size_t RecvOffset = flbcomm::AverSrcNum(LatSet::d);
  // static constexpr std::size_t SendRecvOffset = RecvOffset + 1;

  T averRho = getInterpolation<D, T, LatSet::d>(nRhoF, i, comm.SendRecvCells);
  Vector<T, LatSet::d> averU = getInterpolation<D, T, LatSet::d>(nUF, i, comm.SendRecvCells);

  std::array<T, LatSet::q> feq{};
  Equilibrium<T, LatSet>::SecondOrder(feq, averU, averRho);

  flbcomm::checkCellIndex(comm.SendRecvCells[i + RecvOffset], BLat.getN(),
                          BLat.getBlock().getBlockId(), comm.SendBlockId,
                          "popIntpCommunicate");
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
