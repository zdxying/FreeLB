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

// #include "data_struct/field_struct.h"
#include "data_struct/field_statistics.h"
#include "lbm/equilibrium.ur.h"


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

  Cell<T, LatSet, TypePack> getNeighbor(unsigned int i) const {
    return Cell<T, LatSet, TypePack>(Id + Lat.getDelta_Index()[i], Lat);
  }
  Cell<T, LatSet, TypePack> getNeighbor(const Vector<int, LatSet::d>& direction) const {
    return Cell<T, LatSet, TypePack>(Id + direction * Lat.getProjection(), Lat);
  }

  void setId(std::size_t id) { Id = id; }
  // ++id
  void operator++() { ++Id; }

  std::size_t getId() const { return Id; }
  std::size_t getNeighborId(unsigned int i) const { return Id + Lat.getDelta_Index()[i]; }

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

// proxy reference for reduced-precision POP storage (FP16): converts to the
// compute type on read and back to the storage type on write, so dynamics
// code (cell[i] = omega * feq[i] + _omega * cell[i]) compiles unchanged.
// Only instantiated when the storage type differs from the compute type.
template <typename ST, typename CT>
struct PopRef {
  ST* ptr;
  __device__ PopRef(ST* p) : ptr(p) {}
  __device__ operator CT() const { return static_cast<CT>(*ptr); }
  __device__ PopRef& operator=(CT v) {
    *ptr = static_cast<ST>(v);
    return *this;
  }
};

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

  // get population; with reduced-precision storage (FP16) this returns a
  // converting proxy instead of a bare reference
  __device__ decltype(auto) operator[](int i) const {
    if constexpr (std::is_same_v<PopStorage<T>, T>) {
      return Lat->template getField<POP<T, LatSet::q>>().getField(i)[Id];
    } else {
      return PopRef<PopStorage<T>, T>{
          &Lat->template getField<POP<T, LatSet::q>>().getField(i)[Id]};
    }
  }
  __device__ decltype(auto) operator[](int i) {
    if constexpr (std::is_same_v<PopStorage<T>, T>) {
      return Lat->template getField<POP<T, LatSet::q>>().getField(i)[Id];
    } else {
      return PopRef<PopStorage<T>, T>{
          &Lat->template getField<POP<T, LatSet::q>>().getField(i)[Id]};
    }
  }

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

// ---------------------------------------------------------------------------
// A cudev::Cell whose distribution functions live in registers.
//
// collision::BGK reads every population twice: once for the moment pass
// (MomentaScheme::apply) and again inside the collision loop
// (cell[i] = omega * feq[i] + _omega * cell[i]).  cudev::Cell::operator[] only
// ever hands back a reference into global memory, so nothing survives between
// the two passes.  On the 100^3 lid-driven cavity (sm_86) the collision kernel
// therefore issued 438 LDG against only 55 STG for a payload of 19 loads + 19
// stores, which left it at ~70% of DRAM peak.
//
// RegCell resolves the q element addresses once and keeps both the addresses
// and the q values in registers, so the two passes share one set of loads.
// It derives from Cell and overrides only operator[]/flush(), so the collision,
// moment and equilibrium templates are reused verbatim and the physics cannot
// diverge from the baseline.  Result: 78 LDG / 19 STG, ~91% of DRAM peak, and
// 1.33x faster end to end.
//
// Generic over the POP container: the addresses come from
// BlockLatticeBase::getPopArray(), so StreamMapArray and CyclicArray
// both work unchanged.
//
// IMPORTANT: the surrounding task collection must be rebuilt for this cell type
// (tmp::RebindSelector).  Handing a RegCell to a collection built around
// cudev::Cell still compiles and still produces correct results, but slices the
// cell in TaskSelector::Execute and silently falls back to the unoptimised
// baseline -- inspect the SASS load count to tell the two apart.
// ---------------------------------------------------------------------------
template <typename T, typename LatSet, typename TypePack>
class RegCell : public Cell<T, LatSet, TypePack> {
 public:
  using CELL = Cell<T, LatSet, TypePack>;
  using FloatType = T;
  using LatticeSet = LatSet;
  using BLOCKLATTICE = BlockLattice<T, LatSet, TypePack>;
  using GenericRho = typename CELL::GenericRho;
  static constexpr unsigned int Q = LatSet::q;

 private:
  // addresses of the q distribution functions of this cell (storage type;
  // loads/stores convert through __half's implicit float conversions)
  PopStorage<T>* pop_[Q];
  // the q values themselves -- these are what stay in registers
  T cache_[Q];

 public:
  __device__ RegCell(std::size_t id, BLOCKLATTICE* lat) : CELL(id, lat) {
    lat->getPopArray(id, pop_);
#pragma unroll
    for (unsigned int d = 0; d < Q; ++d) cache_[d] = pop_[d][0];
  }

  __device__ T& operator[](int i) { return cache_[i]; }
  __device__ const T& operator[](int i) const { return cache_[i]; }

  // write the q values back in one coalesced pass
  __device__ void flush() {
#pragma unroll
    for (unsigned int d = 0; d < Q; ++d) pop_[d][0] = cache_[d];
  }
};

#endif
} // namespace cudev
