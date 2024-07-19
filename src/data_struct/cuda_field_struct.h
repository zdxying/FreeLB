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

#pragma once

#ifdef __CUDACC__

#include "parallel/communicator.h"
#include "utils/alias.h"
#include "utils/cuda_device.h"
#include "utils/util.h"

namespace cudev {

template <typename T, unsigned int D>
struct BlockComm;

template <typename T, unsigned int D>
struct IntpBlockComm;

template <typename FieldType, typename FloatType, unsigned int Dim>
class BlockField;

template <typename FieldType, typename FloatType, unsigned int Dim>
struct BlockFieldComm {
  BlockField<FieldType, FloatType, Dim>* BlockF;
  BlockComm<FloatType, Dim>* Comm;

  __any__ BlockFieldComm(BlockField<FieldType, FloatType, Dim>* blockF,
                         BlockComm<FloatType, Dim>* comm)
      : BlockF(blockF), Comm(comm) {}

  __device__ const std::size_t* getSends() const { return Comm->SendCells; }
  __device__ std::size_t getSend(std::size_t i) const { return (Comm->SendCells)[i]; }

  __device__ const std::size_t* getRecvs() const { return Comm->RecvCells; }
  __device__ std::size_t getRecv(std::size_t i) const { return (Comm->RecvCells)[i]; }
};

template <typename FieldType, typename FloatType, unsigned int Dim>
struct IntpBlockFieldComm {
  BlockField<FieldType, FloatType, Dim>* BlockF;
  IntpBlockComm<FloatType, Dim>* Comm;

  __any__ IntpBlockFieldComm(BlockField<FieldType, FloatType, Dim>* blockF,
                               IntpBlockComm<FloatType, Dim>* comm)
      : BlockF(blockF), Comm(comm) {}

  __device__ const std::size_t* getSends() const { return Comm->SendCells; }
  __device__ std::size_t getSend(std::size_t i) const { return (Comm->SendCells)[i]; }

  __device__ const std::size_t* getRecvs() const { return Comm->RecvCells; }
  __device__ std::size_t getRecv(std::size_t i) const { return (Comm->RecvCells)[i]; }
};


template <typename FieldType, typename FloatType, unsigned int Dim>
class BlockField : public FieldType {
 public:
  using datatype = typename FieldType::value_type;
  static constexpr unsigned int array_dim = FieldType::array_dim;
  static constexpr bool isField = FieldType::isField;
  using ArrayType = typename FieldType::array_type;

 private:
  // BlockFieldComm<FieldType, FloatType, Dim>* Comms;
  // IntpBlockFieldComm<FieldType, FloatType, Dim>* AverComms;
  // IntpBlockFieldComm<FieldType, FloatType, Dim>* IntpComms;

 public:
  __any__ BlockField(ArrayType** data) : FieldType(data) {}

  FieldType& getFieldType() { return *this; }
  const FieldType& getFieldType() const { return *this; }

  // BlockFieldComm<FieldType, FloatType, Dim>* getComms() { return Comms; }
  // IntpBlockFieldComm<FieldType, FloatType, Dim>* getAverComms() { return AverComms; }
  // IntpBlockFieldComm<FieldType, FloatType, Dim>* getIntpComms() { return IntpComms; }
};


template <typename FieldType, typename FloatType, unsigned int Dim>
class BlockFieldManager {
 private:
  // block fields
  BlockField<FieldType, FloatType, Dim>** _Fields;

 public:
  using datatype = typename FieldType::value_type;
  using array_type = typename FieldType::array_type;
  using field_type = FieldType;
  using float_type = FloatType;
  static constexpr unsigned int dim = Dim;
  static constexpr bool isField = FieldType::isField;

  __any__ BlockFieldManager(BlockField<FieldType, FloatType, Dim>** Fields)
      : _Fields(Fields) {}
  // copy constructor
  __device__ BlockFieldManager(const BlockFieldManager& blockFManager)
      : _Fields(blockFManager._Fields) {}
  // move constructor
  __device__ BlockFieldManager(BlockFieldManager&& blockFManager) noexcept
      : _Fields(blockFManager._Fields){}
  // copy assignment operator
  __device__ BlockFieldManager& operator=(const BlockFieldManager& blockFManager) {
    if (this != &blockFManager) {
      _Fields = blockFManager._Fields;
    }
    return *this;
  }
  // move assignment operator
  __device__ BlockFieldManager& operator=(BlockFieldManager&& blockFManager) noexcept {
    if (this != &blockFManager) {
      _Fields = blockFManager._Fields;
    }
    return *this;
  }

  BlockField<FieldType, FloatType, Dim>& getBlockField(int i) { return _Fields[i]; }
  const BlockField<FieldType, FloatType, Dim>& getBlockField(int i) const {
    return _Fields[i];
  }

};

}  // namespace cudev

#endif