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

// block_lattice.h

#pragma once

#include "data_struct/block_lattice_base.h"
// AbstractConverter
#include "lbm/unit_converter.h"

#include "data_struct/voxel_map.h"

// block structure for refined lattice
template <typename T, typename LatSet, typename TypePack>
class BlockLattice : public BlockLatticeBase<T, LatSet, TypePack> {
 public:
  using CellType = Cell<T, LatSet, TypePack>;
  using LatticeSet = LatSet;
  using FloatType = T;
  using FieldTypePack = TypePack;

  using GenericRho = typename FindGenericRhoType<T, TypePack>::type;

#ifdef __CUDACC__
  using cudev_TypePack = typename ExtractCudevFieldPack<TypePack>::cudev_pack;
  using cudev_BlockLattice = cudev::BlockLattice<T, LatSet, cudev_TypePack>;
#endif

 protected:
  // --- lattice communication structure ---
  std::vector<std::vector<unsigned int>> CommDirection;
#ifdef MPI_ENABLED
  std::vector<std::vector<unsigned int>> RecvDirection;
  std::vector<std::vector<unsigned int>> SendDirection;
#endif

  // omega = 1 / tau
  T Omega;
  // 1 - omega
  T _Omega;
  // 1 - omega/2
  T fOmega;

  // tolerance
  T RhoRes;
  GenericArray<T> RhoOld;
  T URes;
  GenericArray<Vector<T, LatSet::d>> UOld;

#ifdef __CUDACC__
  T* devOmega;
  T* dev_Omega;
  T* dev_fOmega;
  cudev_BlockLattice* dev_BlockLat;
#endif

 public:
  template <typename... FIELDPTRS>
  BlockLattice(Block<T, LatSet::d>& block, AbstractConverter<T>& conv,
               std::tuple<FIELDPTRS...> fieldptrs);

#ifdef __CUDACC__
  void InitDeviceData() {
    devOmega = cuda_malloc<T>(1);
    dev_Omega = cuda_malloc<T>(1);
    dev_fOmega = cuda_malloc<T>(1);
    host_to_device(devOmega, &Omega, 1);
    host_to_device(dev_Omega, &_Omega, 1);
    host_to_device(dev_fOmega, &fOmega, 1);
    constructInDevice();
  }
  void constructInDevice() {
    dev_BlockLat = cuda_malloc<cudev_BlockLattice>(1);
    cudev_BlockLattice temp(this->dev_Delta_Index, this->dev_Fields, devOmega, dev_Omega,
                            dev_fOmega);
    host_to_device(dev_BlockLat, &temp, 1);
  }
  cudev_BlockLattice* get_devObj() { return dev_BlockLat; }
#endif

  std::array<T*, LatSet::q> getPop(std::size_t id) {
    return this->template getField<POP<T, LatSet::q>>().getArray(id);
  }

  inline T getOmega() const { return Omega; }
  inline T get_Omega() const { return _Omega; }
  inline T getfOmega() const { return fOmega; }

  const std::vector<std::vector<unsigned int>>& getCommDirection () const { return CommDirection; }
#ifdef MPI_ENABLED
  const std::vector<std::vector<unsigned int>>& getMPIRecvDirection() const { return RecvDirection; }
  const std::vector<std::vector<unsigned int>>& getMPISendDirection() const { return SendDirection; }
#endif

  void Stream();

  template <typename CELLDYNAMICS, typename ArrayType>
  void ApplyCellDynamics(const ArrayType& flagarr);

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics();

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics(const Genericvector<std::size_t>& Idx);

  template <typename DYNAMICS, typename elementType>
  void ApplyDynamics(const Genericvector<elementType>& Idx);

  template <typename CELLDYNAMICS, typename ArrayType>
  void ApplyInnerCellDynamics(const ArrayType& flagarr);

  template <typename CELLDYNAMICS>
  void ApplyInnerCellDynamics();

#ifdef __CUDACC__

  void CuDevStream();

  template <typename CELLDYNAMICS, typename ArrayType>
  void CuDevApplyCellDynamics(ArrayType& flagarr);

  template <typename CELLDYNAMICS>
  void CuDevApplyCellDynamics();

#endif

  // tolerance
  void EnableToleranceRho(T rhores = T(1e-5));
  void EnableToleranceU(T ures = T(1e-5));
  T getToleranceRho();
  T getToleranceU();
  // get inner block tolerance
  T getTolRho(int shift = 1);
  T getTolU(int shift = 1);

#ifdef MPI_ENABLED

  void mpiNormalSend(std::vector<MPI_Request>& SendRequests, 
  std::vector<std::vector<T>>& SendBuffers, const std::vector<DistributedComm>& MPISends);

  void mpiNormalFullSend(std::vector<MPI_Request>& SendRequests, 
  std::vector<std::vector<T>>& SendBuffers, const std::vector<DistributedComm>& MPISends);

  void mpiAverSend(std::vector<MPI_Request>& SendRequests,
  std::vector<std::vector<T>>& SendBuffers, const std::vector<DistributedComm>& MPISends);

  void mpiIntpSend(std::vector<MPI_Request>& SendRequests, 
  std::vector<std::vector<T>>& SendBuffers, const std::vector<DistributedComm>& MPISends);

  void mpiNormalRecv(std::vector<MPI_Request>& RecvRequests, 
  std::vector<std::vector<T>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs);

  void mpiFullRecv(std::vector<MPI_Request>& RecvRequests, 
  std::vector<std::vector<T>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs);

  void mpiNormalSet(int& reqidx, std::vector<MPI_Request>& RecvRequests,
  const std::vector<std::vector<T>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs);

  void mpiNormalFullSet(int& reqidx, std::vector<MPI_Request>& RecvRequests,
  const std::vector<std::vector<T>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs);

  void mpiAverIntpSet(int& reqidx, std::vector<MPI_Request>& RecvRequests,
  const std::vector<std::vector<T>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs);

  template <unsigned int D>
  void popIntp(T OmegaC, const std::vector<std::size_t>& sends, 
  const GenericArray<T>& RhoArr, const GenericArray<Vector<T, LatSet::d>>& UArr, std::size_t i,
  std::vector<T>& buffer, std::size_t& bufidx);
  
#endif
};

#ifdef __CUDACC__

template <typename T, typename LatSet, typename TypePack>
__global__ void CuDevStreamKernel(cudev::BlockLattice<T, LatSet, TypePack>* blocklat) {
  unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
  blocklat->Stream(id);
}

template <typename T, typename LatSet, typename TypePack, typename CELLDYNAMICS, typename ArrayType>
__global__ void CuDevApplyCellDynamicsKernel(cudev::BlockLattice<T, LatSet, TypePack>* blocklat, ArrayType* flagarr) {
  std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  cudev::Cell<T, LatSet, TypePack> cell(idx, blocklat);
  CELLDYNAMICS::Execute(flagarr->operator[](idx), cell);
}

template <typename T, typename LatSet, typename TypePack, typename CELLDYNAMICS>
__global__ void CuDevApplyCellDynamicsKernel(cudev::BlockLattice<T, LatSet, TypePack>* blocklat) {
  std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  cudev::Cell<T, LatSet, TypePack> cell(idx, blocklat);
  CELLDYNAMICS::apply(cell);
}

// old version of unaligned memory
template <typename T, typename LatSet, typename TypePack, typename CELLDYNAMICS, typename ArrayType>
__global__ void CuDevApplyCellDynamicsKernel(cudev::BlockLattice<T, LatSet, TypePack>* blocklat, ArrayType* flagarr, std::size_t N) {
  std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < N) {
    cudev::Cell<T, LatSet, TypePack> cell(idx, blocklat);
    CELLDYNAMICS::Execute(flagarr->operator[](idx), cell);
  }
}

template <typename T, typename LatSet, typename TypePack, typename CELLDYNAMICS>
__global__ void CuDevApplyCellDynamicsKernel(cudev::BlockLattice<T, LatSet, TypePack>* blocklat, std::size_t N) {
  std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < N) {
    cudev::Cell<T, LatSet, TypePack> cell(idx, blocklat);
    CELLDYNAMICS::apply(cell);
  }
}

#endif

// block lattice manager
template <typename T, typename LatSet, typename TypePack>
class BlockLatticeManager : public BlockLatticeManagerBase<T, LatSet, TypePack> {
 public:
  using FIELDS = typename ExtractFieldPack<TypePack>::pack1;
  using FIELDPTRS = typename ExtractFieldPack<TypePack>::pack2;
  using ALLFIELDS = typename ExtractFieldPack<TypePack>::mergedpack;

  using BLOCKLATTICE = BlockLattice<T, LatSet, ALLFIELDS>;
  using CellType = Cell<T, LatSet, ALLFIELDS>;
  using FloatType = T;
  using LatticeSet = LatSet;

  using GenericRho = typename GetGenericRhoType<T, FIELDS>::type;

 private:
  std::vector<BlockLattice<T, LatSet, ALLFIELDS>> BlockLats;
  AbstractConverter<T>& Conv;

 public:
  template <typename... FIELDPTRTYPES>
  BlockLatticeManager(BlockGeometry<T, LatSet::d>& blockgeo, AbstractConverter<T>& conv,
                      FIELDPTRTYPES*... fieldptrs);
  template <typename INITVALUEPACK, typename... FIELDPTRTYPES>
  BlockLatticeManager(BlockGeometry<T, LatSet::d>& blockgeo, INITVALUEPACK& initvalues,
                      AbstractConverter<T>& conv, FIELDPTRTYPES*... fieldptrs);

  void Init();
  template <typename... InitValues>
  void Init(BlockGeometryHelper<T, LatSet::d>& GeoHelper,
            const std::tuple<InitValues...>& initvalues);

  template <typename Func>
  void ForEachField(Func&& func) {
    this->Fields.ForEachField(std::forward<Func>(func));
  }

  template <typename FieldType>
  void addField(BlockFieldManager<FieldType, T, LatSet::d>& field) {
    this->FieldPtrs.template addField<FieldType>(field);
    // constexpr std::size_t ith = this->FieldPtrs.is_at_index<FieldType>() +
    // FIELDS::size;
    int i = 0;
    for (auto& blocklat : BlockLats) {
      // blocklat.template addField<FieldType, ith>(field);
      blocklat.template addField<FieldType>(field.getBlockField(i));
      ++i;
    }
  }

  AbstractConverter<T>& getConverter() { return Conv; }

  // get block lattice with block id
  BlockLattice<T, LatSet, ALLFIELDS>& findBlockLat(int blockid) {
    for (BlockLattice<T, LatSet, ALLFIELDS>& blocklat : BlockLats) {
      if (blocklat.getBlock().getBlockId() == blockid) return blocklat;
    }
    std::cerr << "[BlockLatticeManager]: can't find blocklat with blockid " << blockid
              << std::endl;
    exit(1);
  }

  inline std::uint8_t getMaxLevel() const { return this->BlockGeo.getMaxLevel(); }

  BlockLattice<T, LatSet, ALLFIELDS>& getBlockLat(int i) { return BlockLats[i]; }
  const BlockLattice<T, LatSet, ALLFIELDS>& getBlockLat(int i) const {
    return BlockLats[i];
  }
  std::vector<BlockLattice<T, LatSet, ALLFIELDS>>& getBlockLats() { return BlockLats; }
  const std::vector<BlockLattice<T, LatSet, ALLFIELDS>>& getBlockLats() const {
    return BlockLats;
  }

  void Stream(std::int64_t count);
  void Stream();

  template <typename CELLDYNAMICS, typename FieldType>
  void ApplyCellDynamics(std::int64_t count,
                         const BlockFieldManager<FieldType, T, LatSet::d>& BFM);
  template <typename CELLDYNAMICS, typename FieldType>
  void ApplyCellDynamics(const BlockFieldManager<FieldType, T, LatSet::d>& BFM);

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics(std::int64_t count);
  template <typename CELLDYNAMICS>
  void ApplyCellDynamics();

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics(std::int64_t count,
                         const GenericvectorManager<std::size_t>& blockids);
  template <typename CELLDYNAMICS>
  void ApplyCellDynamics(const GenericvectorManager<std::size_t>& blockids);

  template <typename DYNAMICS, typename elementType>
  void ApplyDynamics(std::int64_t count,
                     const GenericvectorManager<elementType>& blockids);
  template <typename DYNAMICS, typename elementType>
  void ApplyDynamics(const GenericvectorManager<elementType>& blockids);


  template <typename CELLDYNAMICS, typename FieldType>
  void ApplyInnerCellDynamics(std::int64_t count,
                              const BlockFieldManager<FieldType, T, LatSet::d>& BFM);
  template <typename CELLDYNAMICS, typename FieldType>
  void ApplyInnerCellDynamics(const BlockFieldManager<FieldType, T, LatSet::d>& BFM);

  template <typename CELLDYNAMICS>
  void ApplyInnerCellDynamics(std::int64_t count);
  template <typename CELLDYNAMICS>
  void ApplyInnerCellDynamics();

#ifdef __CUDACC__

  void CuDevStream();

  template <typename CELLDYNAMICS, typename FieldType>
  void CuDevApplyCellDynamics(BlockFieldManager<FieldType, T, LatSet::d>& BFM);

  template <typename CELLDYNAMICS>
  void CuDevApplyCellDynamics();

#endif

  template <typename Func>
  void ForEachBlockLattice(Func&& func) {
    for (auto& blocklat : BlockLats) {
      func(blocklat);
    }
  }

  template <typename Func>
  void ForEachBlockLattice(std::int64_t count, Func&& func) {
    for (auto& blocklat : BlockLats) {
      const int deLevel = static_cast<int>(getMaxLevel() - blocklat.getLevel());
      if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
        func(blocklat);
      }
    }
  }


#ifdef MPI_ENABLED

  void MPINormalCommunicate();
  void MPINormalCommunicate(std::int64_t count);

  void MPINormalFullCommunicate();
  void MPINormalFullCommunicate(std::int64_t count);

  void MPINormalAllCommunicate();
  void MPINormalAllCommunicate(std::int64_t count);

  void MPIAverageCommunicate();
  void MPIAverageCommunicate(std::int64_t count);

  void MPIInterpolateCommunicate();
  void MPIInterpolateCommunicate(std::int64_t count);

#endif

  void NormalCommunicate();
  void NormalCommunicate(std::int64_t count);

  void NormalFullCommunicate();
  void NormalFullCommunicate(std::int64_t count);

  void NormalAllCommunicate();
  void NormalAllCommunicate(std::int64_t count);

  void AverageCommunicate();
  void AverageCommunicate(std::int64_t count);

  void InterpolateCommunicate();
  void InterpolateCommunicate(std::int64_t count);

  void Communicate();
  void Communicate(std::int64_t count);

  void FullCommunicate();
  void FullCommunicate(std::int64_t count);

  void AllCommunicate();
  void AllCommunicate(std::int64_t count);

  // tolerance
  void EnableToleranceRho(T rhores = T(1e-5));
  void EnableToleranceU(T ures = T(1e-5));
  T getToleranceRho(int shift = 0);
  T getToleranceU(int shift = 0);

 private:
  void normalCommunicate(BLOCKLATTICE& BLat);
  void normalFullCommunicate(BLOCKLATTICE& BLat);
  void normalAllCommunicate(BLOCKLATTICE& BLat);
  void averCommunicate(BLOCKLATTICE& BLat, const std::vector<SharedComm>& comms);
  void intpCommunicate(BLOCKLATTICE& BLat, const std::vector<SharedComm>& comms);
  template <unsigned int D>
  void popIntpCommunicate(BLOCKLATTICE& BLat, const BLOCKLATTICE& nBLat, T OmegaC, const SharedComm& comm, 
  const GenericArray<T>& nRhoF, const GenericArray<Vector<T, LatSet::d>>& nUF, std::size_t i);
};

// coupling block lattice manager
template <typename BlockLatManager0, typename BlockLatManager1>
class BlockLatManagerCoupling {
 public:
  using T = typename BlockLatManager0::FloatType;
  using CELL0 = typename BlockLatManager0::CellType;
  using CELL1 = typename BlockLatManager1::CellType;
  using LatSet0 = typename BlockLatManager0::LatticeSet;
  using LatSet1 = typename BlockLatManager1::LatticeSet;

  BlockLatManagerCoupling(BlockLatManager0& blocklatman0, BlockLatManager1& blocklatman1)
      : BlockLatMan0(blocklatman0), BlockLatMan1(blocklatman1) {}

  BlockLatManager0& getLat0() {return BlockLatMan0;}
  BlockLatManager1& getLat1() {return BlockLatMan1;}

  template <typename CELLDYNAMICS, typename BLOCKFIELDMANAGER>
  void ApplyCellDynamics(std::int64_t count, const BLOCKFIELDMANAGER& BFM) {
    int size = static_cast<int>(BlockLatMan0.getBlockLats().size());
#pragma omp parallel for num_threads(Thread_Num)
    for (int iblock = 0; iblock < size; ++iblock) {
      auto& blocklat0 = BlockLatMan0.getBlockLat(iblock);
      auto& blocklat1 = BlockLatMan1.getBlockLat(iblock);
      const int deLevel = static_cast<int>(BlockLatMan0.getMaxLevel() - blocklat0.getLevel());
      if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
        const auto& flagArray = BFM.getBlockField(iblock).getField(0);
        CELL0 cell0(0, blocklat0);
        CELL1 cell1(0, blocklat1);
        const std::size_t voxNum = blocklat0.getVoxNum();
        #ifdef _VOX_ENABLED
        const VoxelMap& map = blocklat0.getBlock().getVoxelMap();
        #endif
        for (std::size_t id = 0; id < voxNum; ++id) {
          #ifdef _VOX_ENABLED
          cell0.setId(map[id]);
          cell1.setId(map[id]);
          CELLDYNAMICS::Execute(flagArray[map[id]], cell0, cell1);
          #else
          cell0.setId(id);
          cell1.setId(id);
          CELLDYNAMICS::Execute(flagArray[id], cell0, cell1);
          #endif
        }
      }
    }
  }

  template <typename CELLDYNAMICS, typename BLOCKFIELDMANAGER>
  void ApplyCellDynamics(const BLOCKFIELDMANAGER& BFM) {
    int size = static_cast<int>(BlockLatMan0.getBlockLats().size());
#pragma omp parallel for num_threads(Thread_Num)
    for (int iblock = 0; iblock < size; ++iblock) {
      auto& blocklat0 = BlockLatMan0.getBlockLat(iblock);
      auto& blocklat1 = BlockLatMan1.getBlockLat(iblock);
      const auto& flagArray = BFM.getBlockField(iblock).getField(0);
      CELL0 cell0(0, blocklat0);
      CELL1 cell1(0, blocklat1);
      const std::size_t voxNum = blocklat0.getVoxNum();
      #ifdef _VOX_ENABLED
      const VoxelMap& map = blocklat0.getBlock().getVoxelMap();
      #endif
      for (std::size_t id = 0; id < voxNum; ++id) {
        #ifdef _VOX_ENABLED
        cell0.setId(map[id]);
        cell1.setId(map[id]);
        CELLDYNAMICS::Execute(flagArray[map[id]], cell0, cell1);
        #else
        cell0.setId(id);
        cell1.setId(id);
        CELLDYNAMICS::Execute(flagArray[id], cell0, cell1);
        #endif
      }
    }
  }

  template <typename DYNAMICS>
  void ApplyCellDynamics() {
    int size = static_cast<int>(BlockLatMan0.getBlockLats().size());
#pragma omp parallel for num_threads(Thread_Num)
    for (int iblock = 0; iblock < size; ++iblock) {
      auto& blocklat0 = BlockLatMan0.getBlockLat(iblock);
      auto& blocklat1 = BlockLatMan1.getBlockLat(iblock);
      CELL0 cell0(0, blocklat0);
      CELL1 cell1(0, blocklat1);
      const std::size_t voxNum = blocklat0.getVoxNum();
      #ifdef _VOX_ENABLED
      const VoxelMap& map = blocklat0.getBlock().getVoxelMap();
      #endif
      for (std::size_t id = 0; id < voxNum; ++id) {
        #ifdef _VOX_ENABLED
        cell0.setId(map[id]);
        cell1.setId(map[id]);
        #else
        cell0.setId(id);
        cell1.setId(id);
        #endif
        DYNAMICS::apply(cell0, cell1);
      }
    }
  }

  template <typename CELLDYNAMICS, typename BLOCKFIELDMANAGER>
  void ApplyInnerCellDynamics(std::int64_t count, const BLOCKFIELDMANAGER& BFM) {
    int size = static_cast<int>(BlockLatMan0.getBlockLats().size());
#pragma omp parallel for num_threads(Thread_Num)
    for (int iblock = 0; iblock < size; ++iblock) {
      auto& blocklat0 = BlockLatMan0.getBlockLat(iblock);
      auto& blocklat1 = BlockLatMan1.getBlockLat(iblock);
      const int deLevel = static_cast<int>(BlockLatMan0.getMaxLevel() - blocklat0.getLevel());
      if (count % (static_cast<int>(std::pow(2, deLevel))) == 0) {
        const auto& flagArray = BFM.getBlockField(iblock).getField(0);
        CELL0 cell0(0, blocklat0);
        CELL1 cell1(0, blocklat1);
        if constexpr (LatSet0::d == 2) {
          for (int j = blocklat0.getOverlap(); j < blocklat0.getNy() - blocklat0.getOverlap(); ++j) {
            std::size_t id = j * blocklat0.getNx() + blocklat0.getOverlap();
            for (int i = blocklat0.getOverlap(); i < blocklat0.getNx() - blocklat0.getOverlap(); ++i) {
              cell0.setId(id);
              cell1.setId(id);
              CELLDYNAMICS::Execute(flagArray[id], cell0, cell1);
              ++id;
            }
          }
        } else if constexpr (LatSet0::d == 3) {
          for (int k = blocklat0.getOverlap(); k < blocklat0.getNz() - blocklat0.getOverlap(); ++k) {
            for (int j = blocklat0.getOverlap(); j < blocklat0.getNy() - blocklat0.getOverlap(); ++j) {
              std::size_t id = k * blocklat0.getProjection()[2] + j * blocklat0.getProjection()[1] + blocklat0.getOverlap();
              for (int i = blocklat0.getOverlap(); i < blocklat0.getNx() - blocklat0.getOverlap(); ++i) {
                cell0.setId(id);
                cell1.setId(id);
                CELLDYNAMICS::Execute(flagArray[id], cell0, cell1);
                ++id;
              }
            }
          }
        }
      }
    }
  }

  template <typename CELLDYNAMICS, typename BLOCKFIELDMANAGER>
  void ApplyInnerCellDynamics(const BLOCKFIELDMANAGER& BFM) {
    int size = static_cast<int>(BlockLatMan0.getBlockLats().size());
#pragma omp parallel for num_threads(Thread_Num)
    for (int iblock = 0; iblock < size; ++iblock) {
      auto& blocklat0 = BlockLatMan0.getBlockLat(iblock);
      auto& blocklat1 = BlockLatMan1.getBlockLat(iblock);
      const auto& flagArray = BFM.getBlockField(iblock).getField(0);
      CELL0 cell0(0, blocklat0);
      CELL1 cell1(0, blocklat1);
      if constexpr (LatSet0::d == 2) {
        for (int j = blocklat0.getOverlap(); j < blocklat0.getNy() - blocklat0.getOverlap(); ++j) {
          std::size_t id = j * blocklat0.getNx() + blocklat0.getOverlap();
          for (int i = blocklat0.getOverlap(); i < blocklat0.getNx() - blocklat0.getOverlap(); ++i) {
            cell0.setId(id);
            cell1.setId(id);
            CELLDYNAMICS::Execute(flagArray[id], cell0, cell1);
            ++id;
          }
        }
      } else if constexpr (LatSet0::d == 3) {
        for (int k = blocklat0.getOverlap(); k < blocklat0.getNz() - blocklat0.getOverlap(); ++k) {
          for (int j = blocklat0.getOverlap(); j < blocklat0.getNy() - blocklat0.getOverlap(); ++j) {
            std::size_t id = k * blocklat0.getProjection()[2] + j * blocklat0.getProjection()[1] + blocklat0.getOverlap();
            for (int i = blocklat0.getOverlap(); i < blocklat0.getNx() - blocklat0.getOverlap(); ++i) {
              cell0.setId(id);
              cell1.setId(id);
              CELLDYNAMICS::Execute(flagArray[id], cell0, cell1);
              ++id;
            }
          }
        }
      }
    }
  }

  template <typename DYNAMICS>
  void ApplyInnerCellDynamics() {
    int size = static_cast<int>(BlockLatMan0.getBlockLats().size());
#pragma omp parallel for num_threads(Thread_Num)
    for (int iblock = 0; iblock < size; ++iblock) {
      auto& blocklat0 = BlockLatMan0.getBlockLat(iblock);
      auto& blocklat1 = BlockLatMan1.getBlockLat(iblock);
      CELL0 cell0(0, blocklat0);
      CELL1 cell1(0, blocklat1);
      if constexpr (LatSet0::d == 2) {
        for (int j = blocklat0.getOverlap(); j < blocklat0.getNy() - blocklat0.getOverlap(); ++j) {
          std::size_t id = j * blocklat0.getNx() + blocklat0.getOverlap();
          for (int i = blocklat0.getOverlap(); i < blocklat0.getNx() - blocklat0.getOverlap(); ++i) {
            cell0.setId(id);
            cell1.setId(id);
            DYNAMICS::apply(cell0, cell1);
            ++id;
          }
        }
      } else if constexpr (LatSet0::d == 3) {
        for (int k = blocklat0.getOverlap(); k < blocklat0.getNz() - blocklat0.getOverlap(); ++k) {
          for (int j = blocklat0.getOverlap(); j < blocklat0.getNy() - blocklat0.getOverlap(); ++j) {
            std::size_t id = k * blocklat0.getProjection()[2] + j * blocklat0.getProjection()[1] + blocklat0.getOverlap();
            for (int i = blocklat0.getOverlap(); i < blocklat0.getNx() - blocklat0.getOverlap(); ++i) {
              cell0.setId(id);
              cell1.setId(id);
              DYNAMICS::apply(cell0, cell1);
              ++id;
            }
          }
        }
      }
    }
  }

 private:
  BlockLatManager0& BlockLatMan0;
  BlockLatManager1& BlockLatMan1;
};


// dynamic block lattice, refine and coarsen based on gradient of rho
template <typename T, typename LatSet, typename TypePack>
class DynamicBlockLatticeHelper2D {
 public:
  using FIELDS = typename ExtractFieldPack<TypePack>::pack1;
  using GenericRho = typename GetGenericRhoType<T, FIELDS>::type;

  using FIELDPTRS = typename ExtractFieldPack<TypePack>::pack2;
  using ALLFIELDS = typename ExtractFieldPack<TypePack>::mergedpack;

 private:
  BlockLatticeManager<T, LatSet, TypePack>& BlockLatMan;
  BlockGeometry2D<T>& BlockGeo;
  BlockGeometryHelper2D<T>& BlockGeoHelper;
  // square norm of gradient of rho field, each BlockCell corresponds to a scalar field
  std::vector<ScalarField<T>> _GradNorm2s;
  // lattice refine threshold
  std::vector<T> _RefineTholds;
  // lattice coarsen threshold
  std::vector<T> _CoarsenTholds;
  // uer defined max refine level
  int _MaxRefineLevel;

  std::vector<T> _MaxGradNorm2s;
  // ScalarField<T> _GradNorm2F;


 public:
  DynamicBlockLatticeHelper2D(BlockLatticeManager<T, LatSet, TypePack>& blocklatman,
                              BlockGeometryHelper2D<T>& geohelper,
                              const std::vector<T>& refineth,
                              const std::vector<T>& coarsenth, int MaxRefineLevel = 2)
      : BlockLatMan(blocklatman), BlockGeo(blocklatman.getGeo()),
        BlockGeoHelper(geohelper), _RefineTholds(refineth), _CoarsenTholds(coarsenth),
        _MaxRefineLevel(MaxRefineLevel)
  // ,_GradNorm2F(BlockGeo.getBaseBlock().getN(), T(0)) {
  {
    // int minsize = std::min(_RefineTholds.size(), _CoarsenTholds.size());
    // _MaxRefineLevel = std::min(_MaxRefineLevel, minsize);
    // init gradnorm2
    for (BasicBlock<T, 2>& block : BlockGeoHelper.getBlockCells()) {
      _GradNorm2s.emplace_back(block.getN(), T(0));
      _MaxGradNorm2s.push_back(T(0));
    }
  }

  void ComputeGradNorm2();
  void UpdateMaxGradNorm2();

  bool WillRefineOrCoarsen();

  void GeoRefine(int OptProcNum, int MaxProcNum = -1, bool enforce = true);

  // field data transfer based on GeoHelper, before re-init geometry
  // pop date transfer between blocks with different level should be treated separately
  // other field data like rho and velo transfer should use Init() in BlockFieldManager
  // you should transfer pop data after other field data hs been transferred
  void PopFieldInit();


  void PopConversionFineToCoarse(const ScalarField<T>& RhoF,
                                 const VectorFieldAOS<T, LatSet::d>& VelocityF,
                                 PopulationField<T, LatSet::q>& PopsF,
                                 const BasicBlock<T, 2>& FBaseBlock,
                                 const BasicBlock<T, 2>& CBlock,
                                 const BasicBlock<T, 2>& CBaseBlock, T OmegaF);

  void PopConversionCoarseToFine(const ScalarField<T>& RhoF,
                                 const VectorFieldAOS<T, LatSet::d>& VelocityF,
                                 PopulationField<T, LatSet::q>& PopsF,
                                 const BasicBlock<T, 2>& CBaseBlock,
                                 const BasicBlock<T, 2>& FBlock,
                                 const BasicBlock<T, 2>& FBaseBlock, T OmegaC);

  // experimental
  std::vector<T>& getMaxGradNorm2s() { return _MaxGradNorm2s; }

  // gather all ScalarField<T> in _GradNorm2s to _GradNorm2F, only used for testing
  // uinform grid void UpdateGradNorm2F(); ScalarField<T>& getGradNorm2F() { return
  // _GradNorm2F; }
};
