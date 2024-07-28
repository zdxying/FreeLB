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

// cell.h
// an interface to access distribution functions
#pragma once

#include <array>
#include <cstdint>

#include "data_struct/field_struct.h"
#include "lbm/equilibrium.h"


template <typename T, typename LatSet>
class PopLattice;

template <typename T, typename LatSet, typename TypePack>
class BlockLattice;

template <typename T, typename LatSet, typename TypePack>
class BlockLatticeBase;

template <typename T, typename LatSet>
class BasicPopCell {
 protected:
  // populations of distribution functions
  std::array<T*, LatSet::q> Pop;

 public:
  BasicPopCell(std::size_t id, PopLattice<T, LatSet>& lat) : Pop(lat.getPop(id)) {}
  template <typename BLOCKLATTICE>
  BasicPopCell(std::size_t id, BLOCKLATTICE& lat) : Pop(lat.getPop(id)) {}

  // access to pop[i]
  const T& operator[](int i) const { return *Pop[i]; }
  T& operator[](int i) { return *Pop[i]; }
};

template <typename T, typename LatSet>
class PopCell final : public BasicPopCell<T, LatSet> {
 protected:
  // global cell index to access field data and distribution functions
  std::size_t Id;
  // reference to lattice
  PopLattice<T, LatSet>& Lat;

 public:
  using FloatType = T;
  using LatticeSet = LatSet;
  PopCell(std::size_t id, PopLattice<T, LatSet>& lat)
      : BasicPopCell<T, LatSet>(id, lat), Id(id), Lat(lat) {}

  PopCell<T, LatSet> getNeighbor(int i) const { return Lat.getNeighbor(*this, i); }
  PopCell<T, LatSet> getNeighbor(const Vector<int, LatSet::d>& direction) const {
    return Lat.getNeighbor(*this, direction);
  }
  PopLattice<T, LatSet>& getLattice() { return Lat; }
  int getNeighborId(int i) const { return Id + Lat.getDelta_Index(i); }

  // get cell index
  std::size_t getId() const { return Id; }
  // get population before streaming
  T& getPrevious(int i) const { return Lat.getPopField().getField(i).getPrevious(Id); }

  // get field
  const T& getRho() const { return Lat.getRhoField().get(Id); }
  T& getRho() { return Lat.getRhoField().get(Id); }
  // Lat.getOmega()
  inline T getOmega() const { return Lat.getOmega(); }
  // Lat.get_Omega()
  inline T get_Omega() const { return Lat.get_Omega(); }
  // Lat.getfOmega()
  inline T getfOmega() const { return Lat.getfOmega(); }

  const Vector<T, LatSet::d>& getVelocity() const { return Lat.getVelocity(Id); }
  Vector<T, LatSet::d>& getVelocity() { return Lat.getVelocity(Id); }
};

// cell interface for block lattice
template <typename T, typename LatSet, typename TypePack>
class Cell {
 protected:
  // global cell index to access field data and distribution functions
  std::size_t Id;
  // reference to lattice
  BlockLattice<T, LatSet, TypePack>& Lat;

 public:
  using FloatType = T;
  using LatticeSet = LatSet;
  using BLOCKLATTICE = BlockLattice<T, LatSet, TypePack>;
  using GenericRho = typename BLOCKLATTICE::GenericRho;

  Cell(std::size_t id, BlockLattice<T, LatSet, TypePack>& lat)
      : Id(id), Lat(lat) {}
  
  // get population
  const T& operator[](int i) const { return Lat.template getField<POP<T, LatSet::q>>().getField(i)[Id]; }
  T& operator[](int i) { return Lat.template getField<POP<T, LatSet::q>>().getField(i)[Id]; }

  template <typename FieldType, unsigned int i = 0>
  auto& get() {
    if constexpr (FieldType::isField) {
      return Lat.template getField<FieldType>().template get<i>(Id);
    } else {
      return Lat.template getField<FieldType>().template get<i>();
    }
  }
  template <typename FieldType, unsigned int i = 0>
  const auto& get() const {
    if constexpr (FieldType::isField) {
      return Lat.template getField<FieldType>().template get<i>(Id);
    } else {
      return Lat.template getField<FieldType>().template get<i>();
    }
  }
  template <typename FieldType>
  auto& get(unsigned int i) {
    if constexpr (FieldType::isField) {
      return Lat.template getField<FieldType>().get(Id, i);
    } else {
      return Lat.template getField<FieldType>().get(i);
    }
  }
  template <typename FieldType>
  const auto& get(unsigned int i) const {
    if constexpr (FieldType::isField) {
      return Lat.template getField<FieldType>().get(Id, i);
    } else {
      return Lat.template getField<FieldType>().get(i);
    }
  }

  template <typename FieldType>
  auto& getField() {
    return Lat.template getField<FieldType>();
  }
  template <typename FieldType>
  const auto& getField() const {
    return Lat.template getField<FieldType>();
  }

  template <typename FieldType>
  static constexpr bool hasField() {
    return BLOCKLATTICE::template hasField<FieldType>();
  }

  Cell<T, LatSet, TypePack> getNeighbor(int i) const {
    return Cell<T, LatSet, TypePack>(Id + Lat.getDelta_Index()[i], Lat);
  }
  Cell<T, LatSet, TypePack> getNeighbor(const Vector<int, LatSet::d>& direction) const {
    return Cell<T, LatSet, TypePack>(Id + direction * Lat.getProjection());
  }

  void setId(std::size_t id) { Id = id; }
  // ++id
  void operator++() { ++Id; }

  std::size_t getId() const { return Id; }
  std::size_t getNeighborId(int i) const { return Id + Lat.getDelta_Index()[i]; }

  // get population before streaming
  T& getPrevious(int i) const {
    return Lat.template getField<POP<T, LatSet::q>>().getField(i).getPrevious(Id);
  }
  // Lat.getOmega()
  inline T getOmega() const { return Lat.getOmega(); }
  // Lat.get_Omega()
  inline T get_Omega() const { return Lat.get_Omega(); }
  // Lat.getfOmega()
  inline T getfOmega() const { return Lat.getfOmega(); }
};

// a generic cell interface for block lattice structure, can't access pops through []
// operator
template <typename T, typename LatSet, typename TypePack>
class GenericCell {
 protected:
  // global cell index to access field data and distribution functions
  std::size_t Id;
  // reference to lattice
  BlockLatticeBase<T, LatSet, TypePack>& Lat;

 public:
  using FloatType = T;
  using LatticeSet = LatSet;
  using BLOCKLATTICE = BlockLatticeBase<T, LatSet, TypePack>;

  using GenericRho = typename BLOCKLATTICE::GenericRho;

  GenericCell(std::size_t id, BlockLatticeBase<T, LatSet, TypePack>& lat)
      : Id(id), Lat(lat) {}

  template <typename FieldType, unsigned int i = 0>
  auto& get() {
    if constexpr (FieldType::isField) {
      return Lat.template getField<FieldType>().template get<i>(Id);
    } else {
      return Lat.template getField<FieldType>().template get<i>();
    }
  }
  template <typename FieldType, unsigned int i = 0>
  const auto& get() const {
    if constexpr (FieldType::isField) {
      return Lat.template getField<FieldType>().template get<i>(Id);
    } else {
      return Lat.template getField<FieldType>().template get<i>();
    }
  }
  template <typename FieldType>
  auto& get(unsigned int i) {
    if constexpr (FieldType::isField) {
      return Lat.template getField<FieldType>().get(Id, i);
    } else {
      return Lat.template getField<FieldType>().get(i);
    }
  }
  template <typename FieldType>
  const auto& get(unsigned int i) const {
    if constexpr (FieldType::isField) {
      return Lat.template getField<FieldType>().get(Id, i);
    } else {
      return Lat.template getField<FieldType>().get(i);
    }
  }

  GenericCell<T, LatSet, TypePack> getNeighbor(int i) const {
    return GenericCell<T, LatSet, TypePack>(Id + Lat.getDelta_Index()[i], Lat);
  }
  GenericCell<T, LatSet, TypePack> getNeighbor(
    const Vector<int, LatSet::d>& direction) const {
    return GenericCell<T, LatSet, TypePack>(Id + direction * Lat.getProjection());
  }

  // get cell index
  std::size_t getId() const { return Id; }
  std::size_t getNeighborId(int i) const { return Id + Lat.getDelta_Index()[i]; }
};


// gpu representation of cell interface
namespace cudev {

#ifdef __CUDACC__

template <typename T, typename LatSet, typename TypePack>
class BlockLattice;

template <typename T, typename LatSet, typename TypePack>
class BlockLatticeBase;

template <typename T, typename LatSet, typename TypePack>
class Cell {
 protected:
  // global cell index to access field data and distribution functions
  std::size_t Id;
  // reference to lattice
  BlockLattice<T, LatSet, TypePack>* Lat;

 public:
  using FloatType = T;
  using LatticeSet = LatSet;
  using BLOCKLATTICE = BlockLattice<T, LatSet, TypePack>;
  using GenericRho = typename BLOCKLATTICE::GenericRho;

  __device__ Cell(std::size_t id, BlockLattice<T, LatSet, TypePack>* lat)
      : Id(id), Lat(lat) {}

  // get population
  __device__ const T& operator[](int i) const { return Lat->template getField<POP<T, LatSet::q>>().getField(i)[Id]; }
  __device__ T& operator[](int i) { return Lat->template getField<POP<T, LatSet::q>>().getField(i)[Id]; }

  template <typename FieldType, unsigned int i = 0>
  __device__ auto& get() {
    using cudev_FieldType = typename GetCuDevFieldType<FieldType>::type;
    if constexpr (cudev_FieldType::isField) {
      return Lat->template getField<cudev_FieldType>().template get<i>(Id);
    } else {
      return Lat->template getField<cudev_FieldType>().template get<i>();
    }
  }
  template <typename FieldType, unsigned int i = 0>
  __device__ const auto& get() const {
    using cudev_FieldType = typename GetCuDevFieldType<FieldType>::type;
    if constexpr (cudev_FieldType::isField) {
      return Lat->template getField<cudev_FieldType>().template get<i>(Id);
    } else {
      return Lat->template getField<cudev_FieldType>().template get<i>();
    }
  }
  template <typename FieldType>
  __device__ auto& get(unsigned int i) {
    using cudev_FieldType = typename GetCuDevFieldType<FieldType>::type;
    if constexpr (cudev_FieldType::isField) {
      return Lat->template getField<cudev_FieldType>().get(Id, i);
    } else {
      return Lat->template getField<cudev_FieldType>().get(i);
    }
  }
  template <typename FieldType>
  __device__ const auto& get(unsigned int i) const {
    using cudev_FieldType = typename GetCuDevFieldType<FieldType>::type;
    if constexpr (cudev_FieldType::isField) {
      return Lat->template getField<cudev_FieldType>().get(Id, i);
    } else {
      return Lat->template getField<cudev_FieldType>().get(i);
    }
  }

  // template <typename FieldType>
  // static constexpr bool hasField() {
  //   return BLOCKLATTICE::template hasField<FieldType>();
  // }

  __device__ Cell<T, LatSet, TypePack> getNeighbor(int i) const {
    return Cell<T, LatSet, TypePack>(Id + Lat->getDelta_Index()[i], Lat);
  }
  __device__ Cell<T, LatSet, TypePack> getNeighbor(const Vector<int, LatSet::d>& direction) const {
    return Cell<T, LatSet, TypePack>(Id + direction * Lat->getProjection());
  }

  __device__ inline void setId(std::size_t id) { Id = id; }

  __device__ std::size_t getId() const { return Id; }
  __device__ std::size_t getNeighborId(int i) const { return Id + Lat->getDelta_Index()[i]; }

  // get population before streaming
  __device__ T& getPrevious(int i) const {
    return Lat->template getField<POP<T, LatSet::q>>().getField(i).getPrevious(Id);
  }
  __device__ inline T getOmega() const { return Lat->getOmega(); }
  __device__ inline T get_Omega() const { return Lat->get_Omega(); }
  __device__ inline T getfOmega() const { return Lat->getfOmega(); }
};

#endif
} // namespace cudev
