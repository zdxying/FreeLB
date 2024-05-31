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
#include <cstddef>
#include <type_traits>

#include "data_struct/field_struct.h"
#include "lbm/equilibrium.h"
#include "lbm/unit_converter.h"

template <typename T, typename LatSet>
class BasicLattice;

template <typename T, typename LatSet, typename TypePack>
class BlockLattice;

template <typename T, typename LatSet>
class BasicCell {
 protected:
  // populations of distribution functions
  std::array<T*, LatSet::q> Pop;

 public:
  BasicCell(std::size_t id, BasicLattice<T, LatSet>& lat) : Pop(lat.getPop(id)) {}
  template <typename BLOCKLATTICE>
  BasicCell(std::size_t id, BLOCKLATTICE& lat) : Pop(lat.getPop(id)) {}

  // access to pop[i]
  const T& operator[](int i) const { return *Pop[i]; }
  T& operator[](int i) { return *Pop[i]; }
  // get pointer to pop[i]
  T* getPopPtr(int i) { return Pop[i]; }
  // get pop array: std::array<T*, LatSet::q>
  std::array<T*, LatSet::q>& getPops() { return Pop; }
};

template <typename T, typename LatSet>
class Cell final : public BasicCell<T, LatSet> {
 protected:
  // global cell index to access field data and distribution functions
  std::size_t Id;
  // reference to lattice
  BasicLattice<T, LatSet>& Lat;

 public:
  using FloatType = T;
  using LatticeSet = LatSet;
  Cell(std::size_t id, BasicLattice<T, LatSet>& lat)
      : Id(id), Lat(lat), BasicCell<T, LatSet>(id, lat) {}

  Cell<T, LatSet> getNeighbor(int i) const { return Lat.getNeighbor(*this, i); }
  Cell<T, LatSet> getNeighbor(const Vector<int, LatSet::d>& direction) const {
    return Lat.getNeighbor(*this, direction);
  }
  BasicLattice<T, LatSet>& getLattice() { return Lat; }
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
class BCell final : public BasicCell<T, LatSet> {
 protected:
  // global cell index to access field data and distribution functions
  std::size_t Id;
  // reference to lattice
  BlockLattice<T, LatSet, TypePack>& Lat;

 public:
  using FloatType = T;
  using LatticeSet = LatSet;
  using BLOCKLATTICE = BlockLattice<T, LatSet, TypePack>;

  BCell(std::size_t id, BlockLattice<T, LatSet, TypePack>& lat)
      : Id(id), Lat(lat), BasicCell<T, LatSet>(id, lat) {}

  template <typename FieldType, unsigned int i = 0>
  auto& get() {
    return Lat.template getField<FieldType>().template get<i>(Id);
  }
  template <typename FieldType, unsigned int i = 0>
  const auto& get() const {
    return Lat.template getField<FieldType>().template get<i>(Id);
  }
  template <typename FieldType>
  auto& get(unsigned int i) {
    return Lat.template getField<FieldType>().get(Id, i);
  }
  template <typename FieldType>
  const auto& get(unsigned int i) const {
    return Lat.template getField<FieldType>().get(Id, i);
  }

  BCell<T, LatSet, TypePack> getNeighbor(int i) const {
    return BCell<T, LatSet, TypePack>(Id + Lat.getDelta_Index()[i], Lat);
  }
  BCell<T, LatSet, TypePack> getNeighbor(const Vector<int, LatSet::d>& direction) const {
    return BCell<T, LatSet, TypePack>(Id + direction * Lat.getProjection());
  }

  // get cell index
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

  const Vector<T, LatSet::d>& getVelocity() const {
    return Lat.template getField<VELOCITY<T, LatSet::d>>().get(Id);
  }
  Vector<T, LatSet::d>& getVelocity() {
    return Lat.template getField<VELOCITY<T, LatSet::d>>().get(Id);
  }
};
