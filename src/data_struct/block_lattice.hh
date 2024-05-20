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

#include "block_lattice.h"
#include "data_struct/block_lattice.h"

template <typename T, typename LatSet>
BlockLattice<T, LatSet>::BlockLattice(Block<T, LatSet::d>& block, ScalerField<T>& rho,
                                      VectorFieldAOS<T, LatSet::d>& velocity,
                                      PopulationField<T, LatSet::q>& pops,
                                      AbstractConverter<T>& conv, bool initpop)
    : BlockGeo(block), BlockRhoLattice<T>(conv, rho), Velocity(velocity), Pops(pops),
      Omega(RefineConverter<T>::getOmegaF(conv.getOMEGA(), block.getLevel())) {
  _Omega = T(1) - Omega;
  fOmega = T(1) - T(0.5) * Omega;
  Delta_Index =
    make_Array<int, LatSet::q>([&](int i) { return LatSet::c[i] * getProjection(); });
  // init populations
  if (initpop) {
    for (int i = 0; i < LatSet::q; ++i) {
      Pops.getField(i).Init(this->getLatRhoInit() * LatSet::w[i]);
    }
  }
}

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::PopConvFineToCoarse() {
  const auto& VeloArr = Velocity.getField(0);
  const auto& RhoArr = this->Rho.getField(0);
  // post average comm
  for (InterpBlockLatComm<T, LatSet>& comm : AverageComm) {
    const T OmegaF = comm.SendBlock->getOmega();
    for (std::size_t idrecv : comm.getRecvs()) {
      std::array<T, LatSet::q> feq{};
      Equilibrium<T, LatSet>::SecondOrder(feq, VeloArr[idrecv], RhoArr[idrecv]);
      // convert
      RefineConverter<T>::template computePopC<LatSet::q>(
        Pops.getArray(idrecv), feq, OmegaF);
    }
  }
}

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::PopConvCoarseToFine() {
  const auto& VeloArr = Velocity.getField(0);
  const auto& RhoArr = this->Rho.getField(0);
  // post interp comm
  for (InterpBlockLatComm<T, LatSet>& comm : InterpComm) {
    const T OmegaC = comm.SendBlock->getOmega();
    for (std::size_t idrecv : comm.getRecvs()) {
      std::array<T, LatSet::q> feq{};
      Equilibrium<T, LatSet>::SecondOrder(feq, VeloArr[idrecv], RhoArr[idrecv]);
      // convert
      RefineConverter<T>::template computePopF<LatSet::q>(
        Pops.getArray(idrecv), feq, OmegaC);
    }
  }
}

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::communicate() {
  // same level communication
  for (BlockLatComm<T, LatSet>& comm : Communicators) {
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
  for (InterpBlockLatComm<T, LatSet>& comm : AverageComm) {
    BlockLattice<T, LatSet>* nBlockLat = comm.SendBlock;
    int size = comm.getRecvs().size();
    constexpr T weight = comm.Comm->getUniformWeight();
    const T OmegaF = nBlockLat->getOmega();
    GenericArray<T>& nRhoF = nBlockLat->getRhoField().getField(0);
    GenericArray<Vector<T, LatSet::d>>& nUF = nBlockLat->getVelocityField().getField(0);
    for (std::size_t i = 0; i < size; ++i) {
      const InterpSource<LatSet::d>& sends = comm.getSends()[i];
      // get averaged(interpolated) rho and u
      T averRho = getAverage<T, LatSet::d>(nRhoF, sends);
      Vector<T, LatSet::d> averU = getAverage<T, LatSet::d>(nUF, sends);
      // get feq, peparing for pop conversion
      std::array<T, LatSet::q> feq{};
      Equilibrium<T, LatSet>::SecondOrder(feq, averU, averRho);

      std::array<T*, LatSet::q> CellPop = Pops.getArray(comm.getRecvs()[i]);
      for (unsigned int k = 0; k < LatSet::q; ++k) {
        T averpop = getAverage<T, LatSet::d>(nBlockLat->getPopField().getField(k), sends);
        // convert from fine to coarse
        *(CellPop[k]) = RefineConverter<T>::getPopC(averpop, feq[k], OmegaF);
      }
    }
  }
}

template <typename T, typename LatSet>
void BlockLattice<T, LatSet>::interpcommunicate() {
  // interp communication, high level get from low level
  for (InterpBlockLatComm<T, LatSet>& comm : InterpComm) {
    BlockLattice<T, LatSet>* nBlockLat = comm.SendBlock;
    int size = comm.getRecvs().size();
    const T OmegaC = nBlockLat->getOmega();
    GenericArray<T>& nRhoF = nBlockLat->getRhoField().getField(0);
    GenericArray<Vector<T, LatSet::d>>& nUF = nBlockLat->getVelocityField().getField(0);
    if constexpr (LatSet::d == 2) {
      for (std::size_t i = 0; i < size; i += 4) {
        {
          const InterpSource<LatSet::d>& sends = comm.getSends()[i];
          // get averaged(interpolated) rho and u
          T averRho = getInterpolation<0, T, LatSet::d>(nRhoF, sends);
          Vector<T, LatSet::d> averU = getInterpolation<0, T, LatSet::d>(nUF, sends);
          // get feq, peparing for pop conversion
          std::array<T, LatSet::q> feq{};
          Equilibrium<T, LatSet>::SecondOrder(feq, averU, averRho);

          std::array<T*, LatSet::q> CellPop =
            Pops.getArray(comm.getRecvs()[i]);
          for (unsigned int k = 0; k < LatSet::q; ++k) {
            T averpop = getInterpolation<0, T, LatSet::d>(
              nBlockLat->getPopField().getField(k), sends);
            // convert from coarse to fine
            *(CellPop[k]) = RefineConverter<T>::getPopF(averpop, feq[k], OmegaC);
          }
        }
        {
          const InterpSource<LatSet::d>& sends = comm.getSends()[i + 1];
          // get averaged(interpolated) rho and u
          T averRho = getInterpolation<1, T, LatSet::d>(nRhoF, sends);
          Vector<T, LatSet::d> averU = getInterpolation<1, T, LatSet::d>(nUF, sends);
          // get feq, peparing for pop conversion
          std::array<T, LatSet::q> feq{};
          Equilibrium<T, LatSet>::SecondOrder(feq, averU, averRho);

          std::array<T*, LatSet::q> CellPop =
            Pops.getArray(comm.getRecvs()[i + 1]);
          for (unsigned int k = 0; k < LatSet::q; ++k) {
            T averpop = getInterpolation<1, T, LatSet::d>(
              nBlockLat->getPopField().getField(k), sends);
            // convert from coarse to fine
            *(CellPop[k]) = RefineConverter<T>::getPopF(averpop, feq[k], OmegaC);
          }
        }
        {
          const InterpSource<LatSet::d>& sends = comm.getSends()[i + 2];
          // get averaged(interpolated) rho and u
          T averRho = getInterpolation<2, T, LatSet::d>(nRhoF, sends);
          Vector<T, LatSet::d> averU = getInterpolation<2, T, LatSet::d>(nUF, sends);
          // get feq, peparing for pop conversion
          std::array<T, LatSet::q> feq{};
          Equilibrium<T, LatSet>::SecondOrder(feq, averU, averRho);

          std::array<T*, LatSet::q> CellPop =
            Pops.getArray(comm.getRecvs()[i + 2]);
          for (unsigned int k = 0; k < LatSet::q; ++k) {
            T averpop = getInterpolation<2, T, LatSet::d>(
              nBlockLat->getPopField().getField(k), sends);
            // convert from coarse to fine
            *(CellPop[k]) = RefineConverter<T>::getPopF(averpop, feq[k], OmegaC);
          }
        }
        {
          const InterpSource<LatSet::d>& sends = comm.getSends()[i + 3];
          // get averaged(interpolated) rho and u
          T averRho = getInterpolation<3, T, LatSet::d>(nRhoF, sends);
          Vector<T, LatSet::d> averU = getInterpolation<3, T, LatSet::d>(nUF, sends);
          // get feq, peparing for pop conversion
          std::array<T, LatSet::q> feq{};
          Equilibrium<T, LatSet>::SecondOrder(feq, averU, averRho);

          std::array<T*, LatSet::q> CellPop =
            Pops.getArray(comm.getRecvs()[i + 3]);
          for (unsigned int k = 0; k < LatSet::q; ++k) {
            T averpop = getInterpolation<3, T, LatSet::d>(
              nBlockLat->getPopField().getField(k), sends);
            // convert from coarse to fine
            *(CellPop[k]) = RefineConverter<T>::getPopF(averpop, feq[k], OmegaC);
          }
        }
      }
    }
  }
}


template <typename T, typename LatSet>
template <typename ArrayType>
void BlockLattice<T, LatSet>::UpdateRho(const ArrayType& flagarr,
                                        std::uint8_t flag) {
  for (std::size_t id = 0; id < getN(); ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BasicCell<T, LatSet> cell(id, *this);
      moment::Rho<T, LatSet>::apply(cell, this->Rho.get(id));
    }
  }
}

template <typename T, typename LatSet>
template <typename ArrayType>
void BlockLattice<T, LatSet>::UpdateRho_Source(const ArrayType& flagarr,
                                               std::uint8_t flag,
                                               const GenericArray<T>& source) {
  for (std::size_t id = 0; id < getN(); ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BasicCell<T, LatSet> cell(id, *this);
      moment::Rho<T, LatSet>::apply(cell, this->Rho.get(id), source[id]);
    }
  }
}

template <typename T, typename LatSet>
template <typename ArrayType>
void BlockLattice<T, LatSet>::UpdateU(const ArrayType& flagarr,
                                      std::uint8_t flag) {
  for (std::size_t id = 0; id < getN(); ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BasicCell<T, LatSet> cell(id, *this);
      moment::Velocity<T, LatSet>::apply(cell, Velocity.get(id));
    }
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename ArrayType>
void BlockLattice<T, LatSet>::BGK(const ArrayType& flagarr,
                                  std::uint8_t flag) {
  for (std::size_t id = 0; id < getN(); ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BCell<T, LatSet> cell(id, *this);
      collision::BGK<T, LatSet>::template apply<GetFeq>(cell);
    }
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename ArrayType>
void BlockLattice<T, LatSet>::BGK_Source(const ArrayType& flagarr,
                                         std::uint8_t flag,
                                         const GenericArray<T>& source) {
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
  T res0, res1, res2, res;
  T maxres = T(0);
  for (int i = 0; i < getN(); ++i) {
    res0 = std::abs(Velocity.get(i)[0] - UOld[i][0]);
    res1 = std::abs(Velocity.get(i)[1] - UOld[i][1]);
    if constexpr (LatSet::d == 3) {
      res2 = std::abs(Velocity.get(i)[2] - UOld[i][2]);
      res1 = std::max(res1, res2);
    }
    res = std::max(res0, res1);
    maxres = std::max(res, maxres);
    // set UOld
    UOld[i][0] = Velocity.get(i)[0];
    UOld[i][1] = Velocity.get(i)[1];
    if constexpr (LatSet::d == 3) UOld[i][2] = Velocity.get(i)[2];
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
  T res0, res1, res2, res;
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
          res2 = std::abs(Velocity.get(id)[2] - UOld[id][2]);
          res1 = std::max(res1, res2);
          res = std::max(res0, res1);
          maxres = std::max(res, maxres);
          // set UOld
          UOld[id][0] = Velocity.get(id)[0];
          UOld[id][1] = Velocity.get(id)[1];
          UOld[id][2] = Velocity.get(id)[2];
        }
      }
    }
  }
  return maxres;
}

// BlockLatticeManager

template <typename T, typename LatSet>
BlockLatticeManager<T, LatSet>::BlockLatticeManager(
  BlockGeometry<T, LatSet::d>& blockgeo, AbstractConverter<T>& conv,
  BlockFieldManager<VectorFieldAOS<T, LatSet::d>, T, LatSet::d>& blockvelocity)
    : BlockGeo(blockgeo), VelocityFM(blockvelocity), Conv(conv),
      RhoFM(blockgeo, conv.getLatRhoInit()), PopsFM(blockgeo) {
  // init block lattices
  for (int i = 0; i < BlockGeo.getBlocks().size(); ++i) {
    BlockLats.emplace_back(BlockGeo.getBlock(i), RhoFM.getBlockField(i).getField(),
                           VelocityFM.getBlockField(i).getField(),
                           PopsFM.getBlockField(i).getField(), Conv);
  }
  InitComm();
  InitAverComm();
  InitIntpComm();
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::Init() {
  BlockLats.clear();
  for (int i = 0; i < BlockGeo.getBlocks().size(); ++i) {
    BlockLats.emplace_back(BlockGeo.getBlock(i), RhoFM.getBlockField(i).getField(),
                           VelocityFM.getBlockField(i).getField(),
                           PopsFM.getBlockField(i).getField(), Conv, false);
  }
  InitComm();
  InitAverComm();
  InitIntpComm();
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::InitComm() {
  // init based on block geometry
  for (BlockLattice<T, LatSet>& BlockLat : BlockLats) {
    Block<T, LatSet::d>& block = BlockLat.getGeo();
    std::vector<BlockLatComm<T, LatSet>>& latcomm = BlockLat.getCommunicators();
    latcomm.clear();
    for (BlockComm<T, LatSet::d>& comm : block.getCommunicators()) {
      latcomm.emplace_back(&findBlockLat(comm.getSendId()), &comm);
    }
  }
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::InitAverComm() {
  for (BlockLattice<T, LatSet>& BlockLat : BlockLats) {
    Block<T, LatSet::d>& block = BlockLat.getGeo();
    std::vector<InterpBlockLatComm<T, LatSet>>& latcomm = BlockLat.getAverageComm();
    latcomm.clear();
    for (InterpBlockComm<T, LatSet::d>& comm : block.getAverageBlockComm()) {
      latcomm.emplace_back(&findBlockLat(comm.getSendId()), &comm);
    }
  }
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::InitIntpComm() {
  for (BlockLattice<T, LatSet>& BlockLat : BlockLats) {
    Block<T, LatSet::d>& block = BlockLat.getGeo();
    std::vector<InterpBlockLatComm<T, LatSet>>& latcomm = BlockLat.getInterpComm();
    latcomm.clear();
    for (InterpBlockComm<T, LatSet::d>& comm : block.getInterpBlockComm()) {
      latcomm.emplace_back(&findBlockLat(comm.getSendId()), &comm);
    }
  }
}

template <typename T, typename LatSet>
template <typename FieldType>
void BlockLatticeManager<T, LatSet>::UpdateRho(
  std::int64_t count, std::uint8_t flag,
  const BlockFieldManager<FieldType, T, LatSet::d>& BFM) {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(pow(2, deLevel))) == 0)
      BlockLats[i].UpdateRho(BFM.getBlockField(i).getField().getField(0), flag);
  }
}


template <typename T, typename LatSet>
template <typename FieldType>
void BlockLatticeManager<T, LatSet>::UpdateRho_Source(
  std::int64_t count, std::uint8_t flag,
  const BlockFieldManager<FieldType, T, LatSet::d>& BFM,
  const BlockFieldManager<ScalerField<T>, T, LatSet::d>& source) {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(pow(2, deLevel))) == 0)
      BlockLats[i].UpdateRho_Source(BFM.getBlockField(i).getField().getField(0), flag,
                                    source.getBlockField(i).getField().getField(0));
  }
}

template <typename T, typename LatSet>
template <typename FieldType>
void BlockLatticeManager<T, LatSet>::UpdateU(
  std::int64_t count, std::uint8_t flag,
  const BlockFieldManager<FieldType, T, LatSet::d>& BFM) {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(pow(2, deLevel))) == 0)
      BlockLats[i].UpdateU(BFM.getBlockField(i).getField().getField(0), flag);
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename FieldType>
void BlockLatticeManager<T, LatSet>::BGK(
  std::int64_t count, std::uint8_t flag,
  const BlockFieldManager<FieldType, T, LatSet::d>& BFM) {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(pow(2, deLevel))) == 0)
      BlockLats[i].template BGK<GetFeq, typename FieldType::array_type>(
        BFM.getBlockField(i).getField().getField(0), flag);
  }
}


template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename FieldType>
void BlockLatticeManager<T, LatSet>::BGK_Source(
  std::int64_t count, std::uint8_t flag,
  const BlockFieldManager<FieldType, T, LatSet::d>& BFM,
  const BlockFieldManager<ScalerField<T>, T, LatSet::d>& source) {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < BlockLats.size(); ++i) {
    const int deLevel = static_cast<int>(getMaxLevel() - BlockLats[i].getLevel());
    if (count % (static_cast<int>(pow(2, deLevel))) == 0)
      BlockLats[i].template BGK_Source<GetFeq, typename FieldType::array_type>(
        BFM.getBlockField(i).getField().getField(0), flag,
        source.getBlockField(i).getField().getField(0));
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

#ifdef MPI_ENABLED
template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::MPIAverComm(std::int64_t count) {
  mpi().barrier();
  // pop conversion before communication
  std::vector<MPI_Request> SendRequestsPop;
  int ifield = 0;
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    // send data to lower level, deLevel+1
    const int deLevel = static_cast<int>(getMaxLevel() - BLat.getLevel()) + 1;
    if (BLat.getLevel() != std::uint8_t(0)) {
      if ((count % (static_cast<int>(pow(2, deLevel))) == 0) &&
          BLat.getGeo()._NeedMPIComm) {
        const MPIInterpBlockComm<T, LatSet::d>& MPIComm =
          BLat.getGeo().getMPIAverBlockComm();

        const T OmegaF = BLat.getOmega();

        for (int i = 0; i < MPIComm.Senders.size(); ++i) {
          const std::vector<InterpSource<LatSet::d>>& sendcells =
            MPIComm.Senders[i].SendCells;

          std::vector<T>& rhobuffer =
            RhoFM.getBlockField(ifield).getMPIAverBuffer().SendBuffers[i];
          std::vector<Vector<T, LatSet::d>>& ubuffer =
            VelocityFM.getBlockField(ifield).getMPIAverBuffer().SendBuffers[i];
          std::vector<T>& popbuffer =
            PopsFM.getBlockField(ifield).getMPIAverBuffer().SendBuffers[i];

          const auto& RhoArray = RhoFM.getBlockField(ifield).getField().getField(0);
          const auto& UArray = VelocityFM.getBlockField(ifield).getField().getField(0);

          std::size_t bufidx = 0;
          for (const InterpSource<LatSet::d>& sends : sendcells) {
            rhobuffer[bufidx] = getAverage<T, LatSet::d>(RhoArray, sends);
            ubuffer[bufidx] = getAverage<T, LatSet::d>(UArray, sends);
            ++bufidx;
          }
          bufidx = 0;
          for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
            const auto& PopArray = PopsFM.getBlockField(ifield).getField().getField(iArr);
            for (const InterpSource<LatSet::d>& sends : sendcells) {
              popbuffer[bufidx] = getAverage<T, LatSet::d>(PopArray, sends);
              ++bufidx;
            }
          }

          // pop conversion
          const std::size_t shift = sendcells.size();
          for (std::size_t id = 0; id < shift; ++id) {
            std::array<T, LatSet::q> feq{};
            Equilibrium<T, LatSet>::SecondOrder(feq, ubuffer[id], rhobuffer[id]);
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              popbuffer[id + shift * k] =
                RefineConverter<T>::getPopC(popbuffer[id + shift * k], feq[k], OmegaF);
            }
          }
        }
        // non-blocking send pop buffer
        for (int i = 0; i < MPIComm.Senders.size(); ++i) {
          MPI_Request request;
          std::vector<T>& popbuffer =
            PopsFM.getBlockField(ifield).getMPIAverBuffer().SendBuffers[i];
          mpi().iSend(popbuffer.data(), popbuffer.size(), MPIComm.Senders[i].RecvRank,
                      &request, MPIComm.Senders[i].RecvBlockid);
          SendRequestsPop.push_back(request);
        }
      }
    }
    ++ifield;
  }
  // recv pop buffer
  std::vector<MPI_Request> RecvRequestsPop;
  PopsFM.MPIAverRecv(count, RecvRequestsPop);
  // wait for all requests
  MPI_Waitall(SendRequestsPop.size(), SendRequestsPop.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequestsPop.size(), RecvRequestsPop.data(), MPI_STATUSES_IGNORE);
  PopsFM.MPIAverSet(count, RecvRequestsPop);
}

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::MPIInterpComm(std::int64_t count) {
  mpi().barrier();
  // pop conversion before communication
  std::vector<MPI_Request> SendRequestsPop;
  int ifield = 0;
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    // send data to higher level, deLevel-1
    const int deLevel = static_cast<int>(getMaxLevel() - BLat.getLevel()) - 1;
    if (deLevel != -1) {
      if ((count % (static_cast<int>(pow(2, deLevel))) == 0) &&
          BLat.getGeo()._NeedMPIComm) {
        const MPIInterpBlockComm<T, LatSet::d>& MPIComm =
          BLat.getGeo().getMPIInterpBlockComm();

        const T OmegaC = BLat.getOmega();

        for (int i = 0; i < MPIComm.Senders.size(); ++i) {
          const std::vector<InterpSource<LatSet::d>>& sendcells =
            MPIComm.Senders[i].SendCells;

          std::vector<T>& rhobuffer =
            RhoFM.getBlockField(ifield).getMPIInterpBuffer().SendBuffers[i];
          std::vector<Vector<T, LatSet::d>>& ubuffer =
            VelocityFM.getBlockField(ifield).getMPIInterpBuffer().SendBuffers[i];
          std::vector<T>& popbuffer =
            PopsFM.getBlockField(ifield).getMPIInterpBuffer().SendBuffers[i];

          const auto& RhoArray = RhoFM.getBlockField(ifield).getField().getField(0);
          const auto& UArray = VelocityFM.getBlockField(ifield).getField().getField(0);

          const std::size_t size = sendcells.size();
          std::size_t bufidx = 0;
          for (std::size_t i = 0; i < size;) {
            getInterpolation<T, LatSet::d>(RhoArray, sendcells, i, rhobuffer, bufidx);
          }
          bufidx = 0;
          for (std::size_t i = 0; i < size;) {
            getInterpolation<T, LatSet::d>(UArray, sendcells, i, ubuffer, bufidx);
          }
          bufidx = 0;
          for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
            const auto& PopArray = PopsFM.getBlockField(ifield).getField().getField(iArr);
            for (std::size_t i = 0; i < size;) {
              getInterpolation<T, LatSet::d>(PopArray, sendcells, i, popbuffer, bufidx);
            }
          }
          // pop conversion
          const std::size_t shift = sendcells.size();
          for (std::size_t id = 0; id < shift; ++id) {
            std::array<T, LatSet::q> feq{};
            Equilibrium<T, LatSet>::SecondOrder(feq, ubuffer[id], rhobuffer[id]);
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              popbuffer[id + shift * k] =
                RefineConverter<T>::getPopF(popbuffer[id + shift * k], feq[k], OmegaC);
            }
          }
        }
        // non-blocking send pop buffer
        for (int i = 0; i < MPIComm.Senders.size(); ++i) {
          MPI_Request request;
          std::vector<T>& popbuffer =
            PopsFM.getBlockField(ifield).getMPIInterpBuffer().SendBuffers[i];
          mpi().iSend(popbuffer.data(), popbuffer.size(), MPIComm.Senders[i].RecvRank,
                      &request, MPIComm.Senders[i].RecvBlockid);
          SendRequestsPop.push_back(request);
        }
      }
    }
    ++ifield;
  }
  // recv pop buffer
  std::vector<MPI_Request> RecvRequestsPop;
  PopsFM.MPIInterpRecv(count, RecvRequestsPop);
  // wait for all requests
  MPI_Waitall(SendRequestsPop.size(), SendRequestsPop.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitall(RecvRequestsPop.size(), RecvRequestsPop.data(), MPI_STATUSES_IGNORE);
  // set
  PopsFM.MPIInterpSet(count, RecvRequestsPop);
}

#endif

template <typename T, typename LatSet>
void BlockLatticeManager<T, LatSet>::Communicate(std::int64_t count) {
  // --- noraml communication ---
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.communicate();
  }

#ifdef MPI_ENABLED
  PopsFM.MPINormalCommunicate(count);
#endif

  // --- average communication ---
  // RhoFM.AverCommunicate(count);
  // VelocityFM.AverCommunicate(count);
  // PopsFM.AverCommunicate(count);
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.avercommunicate();
    // BLat.PopConvFineToCoarse();
  }

#ifdef MPI_ENABLED
  MPIAverComm(count);
#endif

  // --- interpolation communication ---
  // RhoFM.InterpCommunicate(count);
  // VelocityFM.InterpCommunicate(count);
  // PopsFM.InterpCommunicate(count);
#pragma omp parallel for num_threads(Thread_Num)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
    if (count % (static_cast<int>(pow(2, int(getMaxLevel() - BLat.getLevel())))) == 0)
      BLat.interpcommunicate();
    // BLat.PopConvCoarseToFine();
  }

#ifdef MPI_ENABLED
  MPIInterpComm(count);
#endif

mpi().barrier();
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
  T maxres = T(0);
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : maxres)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
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

template <typename T, typename LatSet>
T BlockLatticeManager<T, LatSet>::getToleranceU(int shift) {
  T maxres = T(0);
#pragma omp parallel for num_threads(Thread_Num) schedule(static) reduction(max : maxres)
  for (BlockLattice<T, LatSet>& BLat : BlockLats) {
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

// DynamicBlockLatticeHelper2D

template <typename T, typename LatSet>
void DynamicBlockLatticeHelper2D<T, LatSet>::ComputeGradNorm2() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int icell = 0; icell < BlockGeoHelper.getBlockCells().size(); ++icell) {
    // rho field
    GenericArray<T>& GradNorm2 = _GradNorm2s[icell].getField(0);

    const BasicBlock<T, 2>& cellblock = BlockGeoHelper.getBlockCell(icell);
    T voxsize = cellblock.getVoxelSize();
    // find corresponding block
    for (int iblock = 0; iblock < BlockGeo.getBlockNum(); ++iblock) {
      const BasicBlock<T, 2>& block = BlockGeo.getBlock(iblock);
      const BasicBlock<T, 2>& baseblock = BlockGeo.getBlock(iblock).getBaseBlock();
      if (isOverlapped(cellblock, baseblock)) {
        const ScalerField<T>& RhoF = BlockLatMan.getBlockLat(iblock).getRhoField();
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

template <typename T, typename LatSet>
void DynamicBlockLatticeHelper2D<T, LatSet>::UpdateMaxGradNorm2() {
  for (int i = 0; i < _GradNorm2s.size(); ++i) {
    GenericArray<T>& GradNorm2A = _GradNorm2s[i].getField(0);
    // get max of gradnorm2
    T maxGradNorm2 = T(0);
    for (int id = 0; id < GradNorm2A.size(); ++id) {
      if (GradNorm2A[id] > maxGradNorm2) maxGradNorm2 = GradNorm2A[id];
    }
    _MaxGradNorm2s[i] = maxGradNorm2;
  }
}

template <typename T, typename LatSet>
bool DynamicBlockLatticeHelper2D<T, LatSet>::WillRefineOrCoarsen() {
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

template <typename T, typename LatSet>
void DynamicBlockLatticeHelper2D<T, LatSet>::GeoRefine(int OptProcNum, int MaxProcNum,
                                                       bool enforce) {
  // post refine
  BlockGeoHelper.CheckRefine();
  // update BasicBlocks in GeoHelper
  BlockGeoHelper.CreateBlocks();
  BlockGeoHelper.AdaptiveOptimization(OptProcNum, MaxProcNum, enforce);
  BlockGeoHelper.LoadBalancing();
}


template <typename T, typename LatSet>
void DynamicBlockLatticeHelper2D<T, LatSet>::PopFieldInit() {
  // rho field data transfer could be done here
  BlockLatMan.getRhoFM().InitAndComm(BlockGeoHelper,
                                     BlockLatMan.getConverter().getLatRhoInit());
  // pop field data, this assumes that other fields like rho and velocity have been
  // transfered
  BlockLatMan.getPopsFM().Init(BlockGeoHelper);
  // now pops field data is transferred, but conversion of distribution function
  // between blocks of different refinement levels is not done yet
#pragma omp parallel for num_threads(Thread_Num)
  for (int inewblock = 0; inewblock < BlockLatMan.getPopsFM().getBlockFields().size();
       ++inewblock) {
    BlockField<PopulationField<T, LatSet::q>, T, LatSet::d>& PopsField =
      BlockLatMan.getPopsFM().getBlockField(inewblock);
    const BlockField<ScalerField<T>, T, LatSet::d>& RhoField =
      BlockLatMan.getRhoFM().getBlockField(inewblock);
    const BlockField<VectorFieldAOS<T, LatSet::d>, T, LatSet::d>& VelocityField =
      BlockLatMan.getVelocityFM().getBlockField(inewblock);

    const BasicBlock<T, 2>& newblock = BlockGeo.getBlock(inewblock);
    const BasicBlock<T, 2>& newbaseblock = BlockGeo.getBlock(inewblock).getBaseBlock();
    std::uint8_t Level = newblock.getLevel();
    // find overlapped old block field
    for (int iblock = 0; iblock < BlockGeo.getBlockNum(); ++iblock) {
      const BasicBlock<T, 2>& baseblock = BlockGeoHelper.getAllOldBasicBlock(iblock);
      int overlap = (baseblock.getLevel() != std::uint8_t(0)) ? 2 : 1;
      const BasicBlock<T, 2>& block = baseblock.getExtBlock(overlap);
      if (isOverlapped(newbaseblock, baseblock)) {
        // get omega
        T omega = BlockLatMan.getBlockLat(iblock).getOmega();
        if (Level > block.getLevel()) {
          PopConversionCoarseToFine(RhoField.getField(), VelocityField.getField(),
                                    PopsField.getField(), baseblock, newblock,
                                    newbaseblock, omega);
        } else if (Level < block.getLevel()) {
          PopConversionFineToCoarse(RhoField.getField(), VelocityField.getField(),
                                    PopsField.getField(), baseblock, newblock,
                                    newbaseblock, omega);
        }
      }
    }
  }
  BlockLatMan.getPopsFM().CommunicateAll();
}

template <typename T, typename LatSet>
void DynamicBlockLatticeHelper2D<T, LatSet>::PopConversionFineToCoarse(
  const ScalerField<T>& RhoF, const VectorFieldAOS<T, LatSet::d>& VelocityF,
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
      std::array<T*, LatSet::q> Pops = PopsF.getArray(id);
      // get feq, peparing for pop conversion
      std::array<T, LatSet::q> feq{};
      Equilibrium<T, LatSet>::SecondOrder(feq, VelocityArr[id], RhoArr[id]);
      // convert from fine to coarse
      for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
        T popC = RefineConverter<T>::getPopC(*(Pops[iArr]), feq[iArr], OmegaF);
        *(Pops[iArr]) = popC;
      }
    }
  }
}


template <typename T, typename LatSet>
void DynamicBlockLatticeHelper2D<T, LatSet>::PopConversionCoarseToFine(
  const ScalerField<T>& RhoF, const VectorFieldAOS<T, LatSet::d>& VelocityF,
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
      std::array<T*, LatSet::q> Pops = PopsF.getArray(id);
      // get feq, peparing for pop conversion
      std::array<T, LatSet::q> feq{};
      Equilibrium<T, LatSet>::SecondOrder(feq, VelocityArr[id], RhoArr[id]);
      // convert from fine to coarse
      for (unsigned int iArr = 0; iArr < LatSet::q; ++iArr) {
        T popF = RefineConverter<T>::getPopF(*(Pops[iArr]), feq[iArr], OmegaC);
        *(Pops[iArr]) = popF;
      }
    }
  }
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
