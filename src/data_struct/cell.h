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
#include "data_struct/pop_policy.h"
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


// ---------------------------------------------------------------------------
// POP storage strategy.
//
// The strategy is a compile-time template parameter of Cell rather than a base
// class, because the CSE-generated .ur.h specialisations are partial
// specialisations on the exact cell type: `rhoUImpl<CELL<T, D3Q19<T>, TypePack,
// POPPOLICY>, W>`.  A derived class such as a former `RegCell` would not match
// that pattern (C++ partial-specialisation matching does not consider base
// classes), silently falling back to the loop-based primary template.  With the
// strategy as a parameter, one emitted specialisation covers every strategy and
// csegen needs no change to accommodate a new one.
//
// The two strategies attack different bottlenecks and are complementary --
// measured on D3Q19/sm_86, see docs/CELL_POLICY_PLAN.md:
//   DirectPop  direct global loads/stores through the POP container; the
//              CSE-generated arithmetic still applies.
//   RegPop     the q populations are pulled into registers for the whole
//              duration of the cell dynamics, cutting POP traffic to exactly
//              q reads + q writes per cell, at the cost of ~q registers.
// ---------------------------------------------------------------------------

// gpu representation of cell interface
namespace cudev {

#ifdef __CUDACC__

template <typename T, typename LatSet, typename TypePack>
class BlockLattice;

template <typename T, typename LatSet, typename TypePack>
class BlockLatticeBase;

// The register cache lives in its own class so that the DirectPop
// instantiation is an empty base and costs zero bytes and zero registers.
// A plain empty *member* would still occupy space: the project is C++17, so
// [[no_unique_address]] is not available.  Cell inherits PopCache privately to
// get the empty-base optimisation.
template <typename T, unsigned int Q, typename POPPOLICY>
struct PopCache {};

template <typename T, unsigned int Q>
struct PopCache<T, Q, RegPop> {
  // The q element addresses, resolved ONCE by the cell constructor.  Re-deriving
  // them per direction inside the load/flush loops looks tempting (40 fewer
  // registers) but measured a 1.75x slowdown on FP16/D3Q19: the data traffic is
  // unchanged while the 64-bit address loads go from 90 to 217, and with FP16 the
  // payload is only half the size so that address traffic dominates.  Resolve
  // once, read/write the payloads from the registers.
  PopStorage<T>* p[Q];
  // the q values themselves -- these are what stay in registers
  T v[Q];

  __device__ void load() {
#pragma unroll
    for (unsigned int d = 0; d < Q; ++d) v[d] = p[d][0];
  }
  __device__ void store() {
#pragma unroll
    for (unsigned int d = 0; d < Q; ++d) p[d][0] = v[d];
  }
};

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
  // A reference proxy has to be assign-transparent, and this is the only way to
  // get it: `cell[i] = cell[j]` binds a PopRef prvalue to the implicitly
  // declared copy-assignment (an exact match) in preference to the converting
  // operator= above, which would need a user-defined conversion.  Without this
  // overload the store lands in a discarded temporary and never reaches memory
  // -- silently, with no diagnostic.  collision::BounceBack is written exactly
  // that way (`cell[i] = cell[iopp]`), so the bounce-back walls were never
  // written at all under FP16 and the field went NaN within 100 steps.
  __device__ PopRef& operator=(const PopRef& rhs) {
    *ptr = *rhs.ptr;
    return *this;
  }
};

// No default argument for POPPOLICY here: it is given on the forward
// declaration in lbm/equilibrium.h, and C++ allows it in only one place.
template <typename T, typename LatSet, typename TypePack, typename POPPOLICY>
class Cell : private PopCache<T, LatSet::q, POPPOLICY> {
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
  static constexpr unsigned int Q = LatSet::q;
  // The launch geometry follows from the storage strategy: the register
  // footprint of RegPop wants a different block size than the streaming
  // baseline.  Re-exported here because that is the one thing the launcher in
  // block_lattice.hh needs, and it keeps the policy tag itself free of any
  // dependency on the CUDA launch configuration.
  static constexpr unsigned int block_size = POPPOLICY::block_size;

  // Under RegPop the constructor resolves the q element addresses once and
  // pulls the q values into PopCache; flush() writes them back.  Under DirectPop
  // there is nothing to do here and flush() is a no-op, which is what lets a
  // single kernel serve both strategies.
  __device__ Cell(std::size_t id, BlockLattice<T, LatSet, TypePack>* lat)
      : Id(id), Lat(lat) {
    if constexpr (std::is_same_v<POPPOLICY, RegPop>) {
      lat->getPopArray(id, this->p);
      this->load();
    }
  }

  // Write the cached populations back.  A no-op under DirectPop, which is what
  // lets a single kernel serve both strategies.
  __device__ void flush() {
    if constexpr (std::is_same_v<POPPOLICY, RegPop>) {
      this->store();
    }
  }

  // get population; under RegPop this reads the register cache instead of
  // global memory.  With reduced-precision storage (FP16) the DirectPop
  // instantiation returns a converting proxy instead of a bare reference.
  __device__ decltype(auto) operator[](int i) const {
    if constexpr (std::is_same_v<POPPOLICY, RegPop>) {
      return this->v[i];
    } else if constexpr (std::is_same_v<PopStorage<T>, T>) {
      return Lat->template getField<POP<T, LatSet::q>>().getField(i)[Id];
    } else {
      return PopRef<PopStorage<T>, T>{
          &Lat->template getField<POP<T, LatSet::q>>().getField(i)[Id]};
    }
  }
  __device__ decltype(auto) operator[](int i) {
    if constexpr (std::is_same_v<POPPOLICY, RegPop>) {
      return this->v[i];
    } else if constexpr (std::is_same_v<PopStorage<T>, T>) {
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

  // A neighbour NEVER inherits the strategy: it is returned by value and the
  // caller normally reads a single direction from it, so a register cache would
  // cost q extra registers for no reuse.  Schemes that pull several neighbours
  // inside one kernel (BounceBackMovingWall, the free-surface schemes) would
  // multiply that.  Read-once neighbours therefore use DirectPop even when the
  // cell itself is RegPop.
  //
  // The asymmetry is deliberate but invisible, so it is spelled out here rather
  // than at getNeighbor(): this alias is what a reader meets in an error message
  // or a decltype, and "NeighborCell is not the same kind of cell as `this`" is
  // exactly the thing that would otherwise look like a bug.
  //
  // Note the host Cell::getNeighbor overloads in freeSurface.h and
  // bounce_back_boundary.h (26 call sites) are the CPU counterparts of these; they
  // return the host cell, which has no such distinction.
  using NeighborCell = Cell<T, LatSet, TypePack, DirectPop>;

  __device__ NeighborCell getNeighbor(int i) const {
    return NeighborCell(Id + Lat->getDelta_Index()[i], Lat);
  }
  __device__ NeighborCell getNeighbor(const Vector<int, LatSet::d>& direction) const {
    return NeighborCell(Id + direction * Lat->getProjection(), Lat);
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
// Register-resident cell.
//
// collision::BGK reads every population twice: once for the moment pass
// (MomentaScheme::apply) and again inside the collision loop
// (cell[i] = omega * feq[i] + _omega * cell[i]).  Under DirectPop
// Cell::operator[] only ever hands back a reference into global memory, so
// nothing survives between the two passes.  On the 100^3 lid-driven cavity
// (sm_86) the collision kernel therefore issued 438 LDG against only 55 STG for
// a payload of 19 loads + 19 stores, which left it at ~70% of DRAM peak.
//
// RegPop pulls the q values into registers for the duration of the cell
// dynamics, so the two passes share one set of loads: 78 LDG / 19 STG, ~91% of
// DRAM peak, 1.33x faster end to end.  The addresses are re-derived per
// direction inside the unrolled load/flush loops rather than being pinned in a
// q-sized pointer array, which keeps them from staying live across the
// dynamics.
//
// RegCell is only a spelling of Cell<..., RegPop>; it is an alias, not a
// derived class.  That distinction is load-bearing: CSE emits its
// specialisations as partial specialisations on the exact cell type
// (`rhoUImpl<CELL<T, D3Q19<T>, TypePack, POPPOLICY>, W>`), and a derived type
// would not match -- it would silently fall back to the loop-based primary
// template.  See docs/CELL_POLICY_PLAN.md.
//
// The task collection must still be rebuilt for this cell type
// (tmp::RebindSelector), because every task bakes its own cell type in.
// ---------------------------------------------------------------------------
template <typename T, typename LatSet, typename TypePack>
using RegCell = Cell<T, LatSet, TypePack, RegPop>;

#endif
} // namespace cudev
