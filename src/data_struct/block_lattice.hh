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

// block_lattice.hh

#pragma once

#include "block_lattice.h"
#include "data_struct/block_lattice.h"

template <typename T, typename LatSet>
BlockLattice<T, LatSet>::BlockLattice(Block<T, LatSet::d>& blockgeo, AbstractConverter<T>& conv,
                                      VectorFieldAOS<T, LatSet::d>& velocity)
    : BlockGeometry(blockgeo), Velocity(velocity),
      Omega(RefineConverter<T>::getOmegaF(conv.GetOMEGA(), blockgeo.getLevel())),
      Pops(blockgeo.getN()), RhoLattice<T>(conv, blockgeo.getN()) {
  _Omega = T(1) - Omega;
  fOmega = T(1) - T(0.5) * Omega;
  Delta_Index = make_Array<int, LatSet::q>([&](int i) { return LatSet::c[i] * getProjection(); });
  // init populations
  for (int i = 0; i < LatSet::q; ++i) {
    Pops.getField(i).Init(this->Lattice_Rho_Init * LatSet::w[i]);
  }
}

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::Refine(bool refinevelo, std::uint8_t deltalevel) {
  BlockGeometry.Refine(deltalevel);
  Delta_Index = make_Array<int, LatSet::q>([&](int i) { return LatSet::c[i] * getProjection(); });
  this->Rho.Resize(getN());
  if (refinevelo) {
    Velocity.Resize(getN());
  }
  // get refined pop field from interpolation
}

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::communicate() {
  // same level communication
  for (BlockLatCommStru<T, LatSet>& comm : Communicators) {
    BlockLattice<T, LatSet>* nBlockLat = comm.SendBlock;
    int size = comm.getRecvs().size();
    for (int k = 1; k < LatSet::q; ++k) {
      const CyclicArray<T>& nPopsk = nBlockLat->getPopField().getField(k);
      CyclicArray<T>& Popsk = Pops.getField(k);
      for (int i = 0; i < size; ++i) {
        std::size_t idrecv = comm.getRecvs()[i];
        std::size_t idsend = comm.getSends()[i];
        Popsk.set(idrecv, nPopsk[idsend]);
      }
    }
  }
}

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::avercommunicate() {
  // average communication, low level get from high level
  for (InterpBlockLatCommStru<T, LatSet>& comm : AverageComm) {
    BlockLattice<T, LatSet>* nBlockLat = comm.SendBlock;
    int size = comm.getRecvs().size();
    T weight = comm.Comm->UniformWeight;
    T OmegaF = nBlockLat->getOmega();
    for (int i = 0; i < size; ++i) {
      std::size_t idrecv = comm.getRecvs()[i];
      InterpSource<LatSet::d>& idsends = comm.getSends()[i];
      // get averaged(interpolated) rho and u
      T averRho;
      Vector<T, LatSet::d> averU;
      GenericArray<T>& nRhoF = nBlockLat->getRhoField().getField(0);
      GenericArray<Vector<T, LatSet::d>>& nUF = nBlockLat->getVelocityField().getField(0);
      if constexpr (LatSet::d == 2) {
        averRho =
          (nRhoF[idsends[0]] + nRhoF[idsends[1]] + nRhoF[idsends[2]] + nRhoF[idsends[3]]) * weight;
        averU = (nUF[idsends[0]] + nUF[idsends[1]] + nUF[idsends[2]] + nUF[idsends[3]]) * weight;
      } else if constexpr (LatSet::d == 3) {
        averRho = (nRhoF[idsends[0]] + nRhoF[idsends[1]] + nRhoF[idsends[2]] + nRhoF[idsends[3]] +
                   nRhoF[idsends[4]] + nRhoF[idsends[5]] + nRhoF[idsends[6]] + nRhoF[idsends[7]]) *
                  weight;
        averU = (nUF[idsends[0]] + nUF[idsends[1]] + nUF[idsends[2]] + nUF[idsends[3]] +
                 nUF[idsends[4]] + nUF[idsends[5]] + nUF[idsends[6]] + nUF[idsends[7]]) *
                weight;
      }
      // get feq, peparing for pop conversion
      std::array<T, LatSet::q> feq{};
      Equilibrium<T, LatSet>::SecondOrder(feq, averU, averRho);

      for (int k = 0; k < LatSet::q; ++k) {
        const CyclicArray<T>& nPopsk = nBlockLat->getPopField().getField(k);
        CyclicArray<T>& Popsk = Pops.getField(k);
        T averpop;
        // get averaged(interpolated) pop, rho and u
        if constexpr (LatSet::d == 2) {
          averpop =
            (nPopsk[idsends[0]] + nPopsk[idsends[1]] + nPopsk[idsends[2]] + nPopsk[idsends[3]]) *
            weight;
        } else if constexpr (LatSet::d == 3) {
          averpop =
            (nPopsk[idsends[0]] + nPopsk[idsends[1]] + nPopsk[idsends[2]] + nPopsk[idsends[3]] +
             nPopsk[idsends[4]] + nPopsk[idsends[5]] + nPopsk[idsends[6]] + nPopsk[idsends[7]]) *
            weight;
        }
        // convert from fine to coarse
        T popC = RefineConverter<T>::getPopC(averpop, feq[k], OmegaF);
        Popsk.set(idrecv, popC);
      }
    }
  }
}

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::interpcommunicate() {
  // average communication, low level get from high level
  for (InterpBlockLatCommStru<T, LatSet>& comm : InterpComm) {
    BlockLattice<T, LatSet>* nBlockLat = comm.SendBlock;
    int size = comm.getRecvs().size();
    T OmegaC = nBlockLat->getOmega();
    for (int i = 0; i < size; ++i) {
      std::size_t idrecv = comm.getRecvs()[i];
      InterpSource<LatSet::d>& idsends = comm.getSends()[i];
      InterpWeight<T, LatSet::d>& weights = comm.getWeights()[i];
      // get averaged(interpolated) rho and u
      T averRho;
      Vector<T, LatSet::d> averU;
      GenericArray<T>& nRhoF = nBlockLat->getRhoField().getField(0);
      GenericArray<Vector<T, LatSet::d>>& nUF = nBlockLat->getVelocityField().getField(0);
      if constexpr (LatSet::d == 2) {
        averRho = (nRhoF[idsends[0]] * weights[0] + nRhoF[idsends[1]] * weights[1] +
                   nRhoF[idsends[2]] * weights[2] + nRhoF[idsends[3]] * weights[3]);
        averU = (nUF[idsends[0]] * weights[0] + nUF[idsends[1]] * weights[1] +
                 nUF[idsends[2]] * weights[2] + nUF[idsends[3]] * weights[3]);
      } else if constexpr (LatSet::d == 3) {
        averRho = (nRhoF[idsends[0]] * weights[0] + nRhoF[idsends[1]] * weights[1] +
                   nRhoF[idsends[2]] * weights[2] + nRhoF[idsends[3]] * weights[3] +
                   nRhoF[idsends[4]] * weights[4] + nRhoF[idsends[5]] * weights[5] +
                   nRhoF[idsends[6]] * weights[6] + nRhoF[idsends[7]] * weights[7]);
        averU = (nUF[idsends[0]] * weights[0] + nUF[idsends[1]] * weights[1] +
                 nUF[idsends[2]] * weights[2] + nUF[idsends[3]] * weights[3] +
                 nUF[idsends[4]] * weights[4] + nUF[idsends[5]] * weights[5] +
                 nUF[idsends[6]] * weights[6] + nUF[idsends[7]] * weights[7]);
      }
      // get feq, peparing for pop conversion
      std::array<T, LatSet::q> feq{};
      Equilibrium<T, LatSet>::SecondOrder(feq, averU, averRho);

      for (int k = 0; k < LatSet::q; ++k) {
        const CyclicArray<T>& nPopsk = nBlockLat->getPopField().getField(k);
        CyclicArray<T>& Popsk = Pops.getField(k);
        T averpop;
        if constexpr (LatSet::d == 2) {
          averpop = (nPopsk[idsends[0]] * weights[0] + nPopsk[idsends[1]] * weights[1] +
                     nPopsk[idsends[2]] * weights[2] + nPopsk[idsends[3]] * weights[3]);
        } else if constexpr (LatSet::d == 3) {
          averpop = (nPopsk[idsends[0]] * weights[0] + nPopsk[idsends[1]] * weights[1] +
                     nPopsk[idsends[2]] * weights[2] + nPopsk[idsends[3]] * weights[3] +
                     nPopsk[idsends[4]] * weights[4] + nPopsk[idsends[5]] * weights[5] +
                     nPopsk[idsends[6]] * weights[6] + nPopsk[idsends[7]] * weights[7]);
        }
        // convert from coarse to fine
        T popF = RefineConverter<T>::getPopF(averpop, feq[k], OmegaC);
        Popsk.set(idrecv, popF);
      }
    }
  }
}

#ifdef MPI_ENABLED

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::MPIcommunicate(MPIBlockCommStru& MPIComm) {
  Mpi().barrier();
  // add to send buffer
  for (int i = 0; i < MPIComm.Senders.size(); ++i) {
    std::vector<T>& buffer = MPIBuffer.SendBuffers[i];
    const std::vector<std::size_t>& sendcells = MPIComm.Senders[i].SendCells;
    int id_buf = 0;
    for (int k = 0; k < LatSet::q; ++k) {
      const CyclicArray<T>& arrk = Pops.getField(k);
      for (std::size_t id : sendcells) {
        buffer[id_buf] = arrk[id];
        ++id_buf;
      }
    }
  }
  // non-blocking send
  std::vector<MPI_Request> requests;
  for (int i = 0; i < MPIComm.Senders.size(); ++i) {
    std::vector<T>& buffer = MPIBuffer.SendBuffers[i];
    MPI_Request request;
    Mpi().iSend(buffer.data(), buffer.size(), MPIComm.Senders[i].RecvRank, &request, 0,
                MPI_COMM_WORLD);
    requests.push_back(request);
  }
  // non-blocking recv
  for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
    std::vector<T>& buffer = MPIBuffer.RecvBuffers[i];
    MPI_Request request;
    Mpi().iRecv(buffer.data(), buffer.size(), MPIComm.Recvers[i].SendRank, &request, 0,
                MPI_COMM_WORLD);
    requests.push_back(request);
  }
  // wait
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  // set field from recv buffer
  for (int i = 0; i < MPIComm.Recvers.size(); ++i) {
    const std::vector<T>& buffer = MPIBuffer.RecvBuffers[i];
    const std::vector<std::size_t>& recvcells = MPIComm.Recvers[i].RecvCells;
    int id_buf = 0;
    for (int k = 0; k < LatSet::q; ++k) {
      CyclicArray<T>& arrk = Pops.getField(k);
      for (std::size_t id : recvcells) {
        arrk.set(id, buffer[id_buf]);
        ++id_buf;
      }
    }
  }
  Mpi().barrier();
}

#endif

template <typename T, typename LatSet>
template <typename flagtype>
void BlockLattice<T, LatSet>::UpdateRho(const GenericArray<flagtype>& flagarr, std::uint8_t flag) {
  for (std::size_t id = 0; id < getN(); ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BasicCell<T, LatSet> cell(id, *this);
      moment::Rho<T, LatSet>::apply(cell, this->Rho.get(id));
    }
  }
}

template <typename T, typename LatSet>
template <typename flagtype>
void BlockLattice<T, LatSet>::UpdateRho_Source(const GenericArray<flagtype>& flagarr,
                                               std::uint8_t flag, const GenericArray<T>& source) {
  for (std::size_t id = 0; id < getN(); ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BasicCell<T, LatSet> cell(id, *this);
      moment::Rho<T, LatSet>::apply(cell, this->Rho.get(id), source[id]);
    }
  }
}

template <typename T, typename LatSet>
template <typename flagtype>
void BlockLattice<T, LatSet>::UpdateU(const GenericArray<flagtype>& flagarr, std::uint8_t flag) {
  for (std::size_t id = 0; id < getN(); ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BasicCell<T, LatSet> cell(id, *this);
      moment::Velocity<T, LatSet>::apply(cell, Velocity.get(id));
    }
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename flagtype>
void BlockLattice<T, LatSet>::BGK(const GenericArray<flagtype>& flagarr, std::uint8_t flag) {
  for (std::size_t id = 0; id < getN(); ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BCell<T, LatSet> cell(id, *this);
      collision::BGK<T, LatSet>::template apply<GetFeq>(cell);
    }
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename flagtype>
void BlockLattice<T, LatSet>::BGK_Source(const GenericArray<flagtype>& flagarr, std::uint8_t flag,
                                         const GenericArray<T>& source) {
  T fOmega = T(1) - Omega * T(0.5);
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (std::size_t id = 0; id < getN(); ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BCell<T, LatSet> cell(id, *this);
      collision::BGK<T, LatSet>::template applySource<GetFeq>(cell, source[id]);
    }
  }
}
template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::Stream() {
  for (int i = 1; i < LatSet::q; ++i) {
    Pops.getField(i).rotate(Delta_Index[i]);
  }
}

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::EnableToleranceRho(T rhores) {
  RhoRes = rhores;
  RhoOld.reserve(getN());
  for (int i = 0; i < getN(); ++i) RhoOld.push_back(this->Rho.get(i));
}

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::EnableToleranceU(T ures) {
  URes = ures;
  UOld.reserve(getN());
  for (int i = 0; i < getN(); ++i) UOld.push_back(Velocity.get(i));
}

template <typename T, typename LatSet>
T BlockLattice<T, LatSet>::getToleranceRho() {
  T res;
  T maxres = T(0);
  for (int i = 0; i < getN(); ++i) {
    res = std::abs(this->Rho.get(i) - RhoOld[i]);
    maxres = std::max(res, maxres);
    RhoOld[i] = this->Rho.get(i);
  }
  return maxres;
}

template <typename T, typename LatSet>
T BlockLattice<T, LatSet>::getToleranceU() {
  T res0, res1, res;
  T maxres = T(0);
  for (int i = 0; i < getN(); ++i) {
    res0 = std::abs(Velocity.get(i)[0] - UOld[i][0]);
    res1 = std::abs(Velocity.get(i)[1] - UOld[i][1]);
    res = std::max(res0, res1);
    maxres = std::max(res, maxres);
    // set UOld
    UOld[i][0] = Velocity.get(i)[0];
    UOld[i][1] = Velocity.get(i)[1];
  }
  return maxres;
}

template <typename T, typename LatSet>
T BlockLattice<T, LatSet>::getTolRho(int shift) {
  T res;
  T maxres = T(0);

  if constexpr (LatSet::d == 2) {
    for (int j = shift; j < getNy() - shift; ++j) {
      for (int i = shift; i < getNx() - shift; ++i) {
        std::size_t id = j * getNx() + i;
        res = std::abs(this->Rho.get(id) - RhoOld[id]);
        maxres = std::max(res, maxres);
        RhoOld[id] = this->Rho.get(id);
      }
    }
  } else if constexpr (LatSet::d == 3) {
    int NxNy = getNx() * getNy();
    for (int k = shift; k < getNz() - shift; ++k) {
      for (int j = shift; j < getNy() - shift; ++j) {
        for (int i = shift; i < getNx() - shift; ++i) {
          std::size_t id = k * NxNy + j * getNx() + i;
          res = std::abs(this->Rho.get(id) - RhoOld[id]);
          maxres = std::max(res, maxres);
          RhoOld[id] = this->Rho.get(id);
        }
      }
    }
  }
  return maxres;
}

template <typename T, typename LatSet>
T BlockLattice<T, LatSet>::getTolU(int shift) {
  T res0, res1, res;
  T maxres = T(0);
  if constexpr (LatSet::d == 2) {
    for (int j = shift; j < getNy() - shift; ++j) {
      for (int i = shift; i < getNx() - shift; ++i) {
        std::size_t id = j * getNx() + i;
        res0 = std::abs(Velocity.get(id)[0] - UOld[id][0]);
        res1 = std::abs(Velocity.get(id)[1] - UOld[id][1]);
        res = std::max(res0, res1);
        maxres = std::max(res, maxres);
        // set UOld
        UOld[id][0] = Velocity.get(id)[0];
        UOld[id][1] = Velocity.get(id)[1];
      }
    }
  } else if constexpr (LatSet::d == 3) {
    int NxNy = getNx() * getNy();
    for (int k = shift; k < getNz() - shift; ++k) {
      for (int j = shift; j < getNy() - shift; ++j) {
        for (int i = shift; i < getNx() - shift; ++i) {
          std::size_t id = k * NxNy + j * getNx() + i;
          res0 = std::abs(Velocity.get(id)[0] - UOld[id][0]);
          res1 = std::abs(Velocity.get(id)[1] - UOld[id][1]);
          res = std::max(res0, res1);
          maxres = std::max(res, maxres);
          // set UOld
          UOld[id][0] = Velocity.get(id)[0];
          UOld[id][1] = Velocity.get(id)[1];
        }
      }
    }
  }
  return maxres;
}

// BlockLatticeManager

template <typename T, typename LatSet>
BlockLatticeManager<T, LatSet>::BlockLatticeManager(BlockGeometry<T, LatSet::d>& blockgeo,
                                                    AbstractConverter<T>& conv,
                                                    BlockVectFieldAOS<T, LatSet::d>& blockvelocity)
    : BlockGeo(blockgeo), BlockVelocity(blockvelocity) {
  // init block lattices
  for (int i = 0; i < BlockGeo.getBlocks().size(); ++i) {
    BlockLats.emplace_back(BlockGeo.getBlock(i), conv, BlockVelocity.getBlockField(i));
  }
  InitCommunicators();
  InitAverComm();
  InitIntpComm();
  UpdateMaxLevel();
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::InitCommunicators() {
  // init based on block geometry
  for (BlockLattice<T, LatSet>& BlockLat : BlockLats) {
    Block<T, LatSet::d>& block = BlockLat.getGeo();
    std::vector<BlockLatCommStru<T, LatSet>>& latcomm = BlockLat.getCommunicators();
    latcomm.clear();
    for (BlockCommStru<T, LatSet::d>& comm : block.getCommunicators()) {
      latcomm.emplace_back(&BlockLats[comm.SendBlock->getBlockId()], &comm);
    }
  }
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::InitAverComm() {
  for (BlockLattice<T, LatSet>& BlockLat : BlockLats) {
    Block<T, LatSet::d>& block = BlockLat.getGeo();
    std::vector<InterpBlockLatCommStru<T, LatSet>>& latcomm = BlockLat.getAverageComm();
    latcomm.clear();
    for (InterpBlockCommStru<T, LatSet::d>& comm : block.getAverageBlockComm()) {
      latcomm.emplace_back(&BlockLats[comm.SendBlock->getBlockId()], &comm);
    }
  }
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::InitIntpComm() {
  for (BlockLattice<T, LatSet>& BlockLat : BlockLats) {
    Block<T, LatSet::d>& block = BlockLat.getGeo();
    std::vector<InterpBlockLatCommStru<T, LatSet>>& latcomm = BlockLat.getInterpComm();
    latcomm.clear();
    for (InterpBlockCommStru<T, LatSet::d>& comm : block.getInterpBlockComm()) {
      latcomm.emplace_back(&BlockLats[comm.SendBlock->getBlockId()], &comm);
    }
  }
}

template <typename T, typename LatSet>
std::vector<RhoLattice<T>*> BlockLatticeManager<T, LatSet>::getRhoLattices() {
  std::vector<RhoLattice<T>*> RhoLatticeVec;
  RhoLatticeVec.reserve(BlockLats.size());
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    RhoLatticeVec.push_back(&BLat);
  }
  return RhoLatticeVec;
}
template <typename T, typename LatSet>
std::vector<ScalerField<T>*> BlockLatticeManager<T, LatSet>::getRhoField() {
  std::vector<ScalerField<T>*> RhoFieldVec;
  RhoFieldVec.reserve(BlockLats.size());
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    RhoFieldVec.push_back(&(BLat.getRhoField()));
  }
  return RhoFieldVec;
}

template <typename T, typename LatSet>
std::vector<PopulationField<T, LatSet::q>*> BlockLatticeManager<T, LatSet>::getPopField() {
  std::vector<PopulationField<T, LatSet::q>*> PopFieldVec;
  PopFieldVec.reserve(BlockLats.size());
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    PopFieldVec.push_back(&(BLat.getPopField()));
  }
  return PopFieldVec;
}

template <typename T, typename LatSet>
template <typename flagtype>
void BlockLatticeManager<T, LatSet>::UpdateRho(std::int64_t count,
                                               std::vector<ScalerField<flagtype>*>& field,
                                               std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < BlockLats.size(); ++i) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BlockLats[i].getLevel())))) == 0)
      BlockLats[i].UpdateRho(field[i]->getField(), flag);
  }
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::UpdateRho(std::int64_t count, std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.UpdateRho(BLat.getGeo().getGeoFlagField().getField(), flag);
  }
}

template <typename T, typename LatSet>
template <typename flagtype>
void BlockLatticeManager<T, LatSet>::UpdateRho_Source(std::int64_t count,
                                                      std::vector<ScalerField<flagtype>*>& field,
                                                      std::uint8_t flag,
                                                      std::vector<ScalerField<T>*>& source) {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < BlockLats.size(); ++i) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BlockLats[i].getLevel())))) == 0)
      BlockLats[i].UpdateRho_Source(field[i]->getField(), flag, source[i]->getField());
  }
}

template <typename T, typename LatSet>
template <typename flagtype>
void BlockLatticeManager<T, LatSet>::UpdateU(std::int64_t count,
                                             std::vector<ScalerField<flagtype>*>& field,
                                             std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < BlockLats.size(); ++i) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BlockLats[i].getLevel())))) == 0)
      BlockLats[i].UpdateU(field[i]->getField(), flag);
  }
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::UpdateU(std::int64_t count, std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.UpdateU(BLat.getGeo().getGeoFlagField().getField(), flag);
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename flagtype>
void BlockLatticeManager<T, LatSet>::BGK(std::int64_t count,
                                         std::vector<ScalerField<flagtype>*>& field,
                                         std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < BlockLats.size(); ++i) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BlockLats[i].getLevel())))) == 0)
      BlockLats[i].template BGK<GetFeq, flagtype>(field[i]->getField(), flag);
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename flagtype>
void BlockLatticeManager<T, LatSet>::BGK(std::int64_t count, std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.template BGK<GetFeq, flagtype>(BLat.getGeo().getGeoFlagField().getField(), flag);
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename flagtype>
void BlockLatticeManager<T, LatSet>::BGK_Source(std::int64_t count,
                                                std::vector<ScalerField<flagtype>*>& field,
                                                std::uint8_t flag,
                                                std::vector<ScalerField<T>*>& source) {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < BlockLats.size(); ++i) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BlockLats[i].getLevel())))) == 0)
      BlockLats[i].template BGK_Source<GetFeq, flagtype>(field[i]->getField(), flag,
                                                         source[i]->getField());
  }
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::Stream(std::int64_t count) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.Stream();
  }
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::Communicate(std::int64_t count) {
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.communicate();
  }
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.avercommunicate();
  }
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.interpcommunicate();
  }
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::EnableToleranceRho(T rhores) {
  for (BlockLattice<T, LatSet>& BLat : BlockLats) BLat.EnableToleranceRho(rhores);
}
template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::EnableToleranceU(T ures) {
  for (BlockLattice<T, LatSet>& BLat : BlockLats) BLat.EnableToleranceU(ures);
}

template <typename T, typename LatSet>
T BlockLatticeManager<T, LatSet>::getToleranceRho(int shift) {
  int size = BlockLats.size();
  std::vector<T> res(size);
  int i = 0;
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (shift == 0) {
      res[i] = BLat.getToleranceRho();
    } else {
      res[i] = BLat.getTolRho(shift);
    }
    ++i;
  }
  T maxres = T(0);
  for (T r : res) maxres = std::max(r, maxres);
  return maxres;
}

template <typename T, typename LatSet>
T BlockLatticeManager<T, LatSet>::getToleranceU(int shift) {
  int size = BlockLats.size();
  std::vector<T> res(size);
  int i = 0;
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (shift == 0) {
      res[i] = BLat.getToleranceU();
    } else {
      res[i] = BLat.getTolU(shift);
    }
    ++i;
  }
  T maxres = T(0);
  for (T r : res) maxres = std::max(r, maxres);
  return maxres;
}

// DynamicBlockLatticeHelper2D

template <typename T, typename LatSet>
void DynamicBlockLatticeHelper2D<T, LatSet>::ComputeGradNorm2() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int i = 0; i < BlockGeoHelper.getBlockCells().size(); ++i) {
    // rho field
    GenericArray<T>& GradNorm2 = _GradNorm2s[i].getField(0);
    // basic block
    BasicBlock<T, 2>& block = BlockGeoHelper.getBlockCell(i);
    Vector<T, 2> minCentre = block.getMinCenter();
    T voxsize = block.getVoxelSize();
    // block centre
    Vector<T, 2> centre = block.getCenter();
    // find corresponding block
    int j = 0;
    for (BlockLattice<T, LatSet>& blocklat : BlockLatMan.getBlockLats()) {
      if (blocklat.getGeo().getBaseBlock().isInside(centre)) {
        break;
      } else {
        ++j;
      }
    }
    // get rho field and block2d
    ScalerField<T>& RhoF = BlockLatMan.getBlockLat(j).getRhoField();
    Block2D<T>& block2d = BlockLatMan.getBlockLat(j).getGeo();
    FDM2D<T> FDM(block2d.getNx(), block2d.getNy(), RhoF.getField(0));
    // start index of block in block2d
    Vector<T, 2> Ext = minCentre - block2d.getMinCenter();
    int x = static_cast<int>(std::round(Ext[0] / voxsize));
    int y = static_cast<int>(std::round(Ext[1] / voxsize));
    // get gradnorm2
    for (int iy = 0; iy < block.getNy(); ++iy) {
      for (int ix = 0; ix < block.getNx(); ++ix) {
        std::size_t id = iy * block.getNx() + ix;
        std::size_t idF = (iy + y) * block2d.getNx() + ix + x;
        GradNorm2[id] = FDM.gradnorm2(idF);
      }
    }
  }
  // experimental
  // UpdateGradNorm2F();
}

template <typename T, typename LatSet>
void DynamicBlockLatticeHelper2D<T, LatSet>::UpdateMaxGradNorm2() {
  for (int i = 0; i < _GradNorm2s.size(); ++i) {
    GenericArray<T>& GradNorm2A = _GradNorm2s[i].getField(0);
    // get max of gradnorm2
    T maxGradNorm2 = T(0);
    for (std::size_t id = 0; id < GradNorm2A.size(); ++id) {
      if (GradNorm2A[id] > maxGradNorm2) maxGradNorm2 = GradNorm2A[id];
    }
    _MaxGradNorm2s[i] = maxGradNorm2;
  }
}

template <typename T, typename LatSet>
void DynamicBlockLatticeHelper2D<T, LatSet>::GeoRefineAndCoarsen(int OptProcNum, int MaxProcNum,
                                                                 bool enforce) {
  ComputeGradNorm2();
  UpdateMaxGradNorm2();

  for (int i = 0; i < _GradNorm2s.size(); ++i) {
    BasicBlock<T, 2>& block = BlockGeoHelper.getBlockCell(i);
    int level = static_cast<int>(block.getLevel());

    if (_MaxGradNorm2s[i] > _RefineThresholds[level] && level < _MaxRefineLevel) {
      block.refine();
    } else if (_MaxGradNorm2s[i] < _CoarsenThresholds[level] && level > 0) {
      block.coarsen();
    }
  }
  // post refine
  BlockGeoHelper.PostRefine();
  // update _BasicBlocks in GeoHelper
  BlockGeoHelper.CreateBlocks();
  BlockGeoHelper.AdaptiveOptimization(OptProcNum, MaxProcNum, enforce);
}

template <typename T, typename LatSet>
void DynamicBlockLatticeHelper2D<T, LatSet>::LatticeRefineAndCoarsen() {
  // lattice data operation
}

// template <typename T, typename LatSet>
// void DynamicBlockLatticeHelper2D<T, LatSet>::UpdateGradNorm2F() {
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
