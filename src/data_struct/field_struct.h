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

#ifdef __CUDACC__
#include "data_struct/cuda_field_struct.h"
#endif


template <typename T>
class Octree;

// BlockField with communication structure
template <typename FieldType, typename FloatType, unsigned int Dim>
class BlockField : public FieldType {
 public:
  using datatype = typename FieldType::value_type;
  static constexpr unsigned int array_dim = FieldType::array_dim;
  static constexpr bool isField = FieldType::isField;

 private:
  // block(geometry) structure of the field
  Block<FloatType, Dim>& _Block;

#ifdef __CUDACC__
  using cudev_FieldType = typename FieldType::cudev_FieldType;
  using cudev_BlockFieldType = cudev::BlockField<cudev_FieldType, FloatType, Dim>;

  cudev_BlockFieldType* dev_BlockField;
#endif

 public:
  BlockField() = default;
  BlockField(Block<FloatType, Dim>& block);
  BlockField(Block<FloatType, Dim>& block, datatype initvalue);
  // copy constructor
  BlockField(const BlockField& blockF);
  // move constructor
  BlockField(BlockField&& blockF) noexcept;

  ~BlockField();

  FieldType& getFieldType() { return *this; }
  const FieldType& getFieldType() const { return *this; }

  Block<FloatType, Dim>& getBlock() { return _Block; }
  const Block<FloatType, Dim>& getBlock() const { return _Block; }

  void constructInDevice();
#ifdef __CUDACC__
  cudev_BlockFieldType* get_devObj() { return dev_BlockField; }
#endif

#ifdef MPI_ENABLED
  // send data for mpi normal communication
  void mpiNormalSend(std::vector<MPI_Request>& SendRequests, 
    std::vector<std::vector<datatype>>& SendBuffers, const std::vector<DistributedComm>& MPISends);

  // send data for mpi average communication
  void mpiAverSend(std::vector<MPI_Request>& SendRequests,
    std::vector<std::vector<datatype>>& SendBuffers, const std::vector<DistributedComm>& MPISends);

  // send data for mpi interp communication
  void mpiIntpSend(std::vector<MPI_Request>& SendRequests, 
    std::vector<std::vector<datatype>>& SendBuffers, const std::vector<DistributedComm>& MPISends);
  
  // recv data for mpi communication
  void mpiRecv(std::vector<MPI_Request>& RecvRequests, 
    std::vector<std::vector<datatype>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs);

  // set data for mpi communication
  void mpiSet(int& reqidx, std::vector<MPI_Request>& RecvRequests,
    const std::vector<std::vector<datatype>>& RecvBuffers, const std::vector<DistributedComm>& MPIRecvs);
#endif

  // copy assignment operator
  BlockField& operator=(const BlockField& blockF) {
    if (this != &blockF) {
      FieldType::operator=(blockF);
      _Block = blockF._Block;
    }
    return *this;
  }
  // move assignment operator
  BlockField& operator=(BlockField&& blockF) noexcept {
    if (this != &blockF) {
      FieldType::operator=(std::move(blockF));
      _Block = blockF._Block;
    }
    return *this;
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

#ifdef __CUDACC__
  using cudev_FieldType = typename FieldType::cudev_FieldType;
  using cudev_BlockFieldType = cudev::BlockField<cudev_FieldType, FloatType, Dim>;
  using cudev_BlockFieldManagerType =
    cudev::BlockFieldManager<cudev_FieldType, FloatType, Dim>;

  cudev_BlockFieldType** dev_Fields;
  cudev_BlockFieldManagerType* dev_BlockFieldManager;
#endif

 public:
  using datatype = typename FieldType::value_type;
  using array_type = typename FieldType::array_type;
  using field_type = FieldType;
  using float_type = FloatType;
  static constexpr unsigned int dim = Dim;
  static constexpr bool isField = FieldType::isField;

  BlockFieldManager(BlockGeometry<FloatType, Dim>& blockgeometry);
  BlockFieldManager(BlockGeometry<FloatType, Dim>& blockgeometry, datatype initvalue);
  // copy constructor
  BlockFieldManager(const BlockFieldManager& blockFManager);
  // move constructor
  BlockFieldManager(BlockFieldManager&& blockFManager) noexcept;

  ~BlockFieldManager() = default;

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

  std::size_t size() const { return _Fields.size(); }

  // get block field with block id
  BlockField<FieldType, FloatType, Dim>* findBlockField(int blockid);

  void InitDeviceData();

#ifdef __CUDACC__
  void copyToDevice() {
    for (auto& field : _Fields) field.copyToDevice();
  }
  void copyToHost() {
    for (auto& field : _Fields) field.copyToHost();
  }
  cudev_BlockFieldType** get_devptr() { return dev_Fields; }
  cudev_BlockFieldManagerType* get_devObj() { return dev_BlockFieldManager; }
  void constructInDevice();
#endif


  void InitAndComm();
  void InitAndComm(datatype initvalue);
  void InitAndComm(BlockGeometryHelper<FloatType, Dim>& GeoHelper);
  void InitAndComm(BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue);
  // pure init
  void Init();
  // pure init with initvalue
  void Init(datatype initvalue);
  void NonFieldInit(datatype initvalue);
  // this assumes that the BlockGeo is already initialized
  void Init(BlockGeometryHelper<FloatType, Dim>& GeoHelper);
  // init with initvalue and data transfer
  void Init(BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue);
  // init with initvalue and data transfer
  template <typename FlagFieldType, typename Func>
  void InitCopy(BlockGeometryHelper<FloatType, Dim>& GeoHelper, datatype initvalue,
                const BlockFieldManager<FlagFieldType, FloatType, Dim>& FlagFManager,
                std::uint8_t flag, const Func& func);

  // construct new field with Geohelper and copy data from old field, then swap
  void FieldDataTransfer(BlockGeometryHelper<FloatType, Dim>& GeoHelper,
                         std::vector<BlockField<FieldType, FloatType, Dim>>& NewFields);

  // construct new field with Geohelper and transfer with average, then swap
  void FieldDataCopyTransfer(
    BlockGeometryHelper<FloatType, Dim>& GeoHelper,
    std::vector<BlockField<FieldType, FloatType, Dim>>& NewFields);

  template <typename LatSet>
  void SetupBoundary(const AABB<FloatType, Dim>& block, datatype bdvalue);

  template <typename LatSet>
  void SetupBoundary(datatype fromvalue, datatype voidvalue, datatype bdvalue);

  void ReadOctree(Octree<FloatType>* tree, datatype stlflag);

  // for flag field, clean lonely flags after reading from octree/stl
  // this should be called after ReadOctree and before SetupBoundary
  // so at this time it is supposed that there are only 2 kinds of flags: void(1) and stl(2)
  // lonelyth is the threshold when cell has less than lonelyth neighbors of the same flag
  // each lonely cell will be set to voidflag after all cells have been checked
  // if recursive is true, the process will be repeated until no lonely cell is found
  template <typename LatSet>
  void CleanLonelyFlags(std::uint8_t flag = std::uint8_t(2), std::uint8_t voidflag = std::uint8_t(1), 
    unsigned int lonelyth = 2, bool recursive = true);


  template <typename Func>
  void forEachField(const Func& func);

  void InitValue(datatype initvalue);

  template <typename Func>
  void forEach(const Func& func);

  template <typename Func>
  void forEachInner(const Func& func);

  template <typename Func, typename Func1>
  void forEach_TransFlag(const Func& func, const Func1& func1);

  // call forEach(AABBs, [&](FieldType& field, std::size_t id){});
  template <typename Func>
  void forEach(const AABB<FloatType, Dim>& AABBs, const Func& func);

  // call forEach(AABBs, FlagFManager, flag, [&](FieldType& field, std::size_t id){});
  template <typename FlagFieldType, typename Func>
  void forEach(const AABB<FloatType, Dim>& AABBs,
               const BlockFieldManager<FlagFieldType, FloatType, Dim>& FlagFManager,
               std::uint8_t flag, const Func& func);

  // call forEach(FlagFManager, flag, [&](FieldType& field, std::size_t id){});
  template <typename FlagFieldType, typename Func>
  void forEach(const BlockFieldManager<FlagFieldType, FloatType, Dim>& FlagFManager,
               std::uint8_t flag, const Func& func);

  // call forEach(FlagFManager, [&](FieldType& field, std::size_t id){});
  template <typename FieldTypeX, typename Func>
  void forEach(const BlockFieldManager<FieldTypeX, FloatType, Dim>& BlockFM, const Func& func);

  // ------ communication ------

#ifdef MPI_ENABLED

  void MPINormalCommunicate();
  void MPIAverCommunicate();
  void MPIIntpCommunicate();

  void MPINormalCommunicate(std::int64_t count);
  void MPIAverCommunicate(std::int64_t count);
  void MPIIntpCommunicate(std::int64_t count);

  void AllMPINormalCommunicate();
  void AllMPIAverCommunicate();
  void AllMPIIntpCommunicate();

  void AllMPINormalCommunicate(std::int64_t count);
  void AllMPIAverCommunicate(std::int64_t count);
  void AllMPIIntpCommunicate(std::int64_t count);
  
#endif

  void NormalCommunicate();
  void AverCommunicate();
  void IntpCommunicate();
  void Communicate();

  void NormalCommunicate(std::int64_t count);
  void AverCommunicate(std::int64_t count);
  void IntpCommunicate(std::int64_t count);
  void Communicate(std::int64_t count);

  void AllNormalCommunicate();
  void AllAverCommunicate();
  void AllIntpCommunicate();
  void AllCommunicate();

  void AllNormalCommunicate(std::int64_t count);
  void AllAverCommunicate(std::int64_t count);
  void AllIntpCommunicate(std::int64_t count);
  void AllCommunicate(std::int64_t count);

// copy assignment operator
  BlockFieldManager& operator=(const BlockFieldManager& blockFManager) {
    if (this != &blockFManager) {
      _Fields = blockFManager._Fields;
      _BlockGeo = blockFManager._BlockGeo;
    }
#ifdef __CUDACC__
    device_to_device(dev_Fields, blockFManager.dev_Fields, _Fields.size());
#endif
    return *this;
  }
  // move assignment operator
  BlockFieldManager& operator=(BlockFieldManager&& blockFManager) noexcept {
    if (this != &blockFManager) {
      _Fields = std::move(blockFManager._Fields);
      _BlockGeo = blockFManager._BlockGeo;
    }
#ifdef __CUDACC__
    dev_Fields = blockFManager.dev_Fields;
#endif
    return *this;
  }
};


// vector of Genericvector
template <typename T>
class GenericvectorManager {
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
  GenericvectorManager(std::size_t size, FLAGFM& FlagFM, FLAGTYPE Flag,
                       std::string name = "GenericvectorManager")
      : _vectors(size) {
    Init(FlagFM, Flag);
    // get vector info
    std::size_t sumsize{};
    for (auto& vec : _vectors) {
      sumsize += vec.size();
    }
    std::cout << "[" << name << "Num of indices]: " << sumsize << std::endl;
  }
  // copy constructor
  GenericvectorManager(const GenericvectorManager& vecManager)
      : _vectors(vecManager._vectors) {}
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
  void Init(FLAGFM& FlagFM, FLAGTYPE Flag) {
    FlagFM.forEachField([&](auto& blockfield) {
      auto& block = blockfield.getBlock();
      auto& VecIds = getvector(block.getBlockId()).getvector();
      block.forEach([&](std::size_t id) {
        if (util::isFlag(blockfield.get(id), Flag)) {
          VecIds.push_back(id);
        }
      });
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