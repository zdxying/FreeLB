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

// block_lattice_base.h

#pragma once

// Cell
#include "data_struct/cell.h"
// Block and BlockGeometry
#include "utils/alias.h"


// base class for block lattice
// without special treatment for population field
// without lattice communication structure
template <typename T, typename LatSet, typename TypePack>
class BlockLatticeBase {
 public:
  using CellType = GenericCell<T, LatSet, TypePack>;
  using LatticeSet = LatSet;
  using FloatType = T;

  using GenericRho = typename FindGenericRhoType<T, TypePack>::type;

 protected:
  // nbr index
  std::array<int, LatSet::q> Delta_Index;
  // geometry
  Block<T, LatSet::d>& BlockGeo;
  // field
  FieldPtrCollection<TypePack> Fields;

 public:
  template <typename... FIELDPTRS>
  BlockLatticeBase(Block<T, LatSet::d>& block, std::tuple<FIELDPTRS...> fieldptrs)
      : BlockGeo(block), Fields(fieldptrs) {
    Delta_Index =
      make_Array<int, LatSet::q>([&](int i) { return LatSet::c[i] * getProjection(); });
  }

  template <typename FieldType>
  auto& getField() {
    return Fields.template getField<FieldType>();
  }
  template <typename FieldType>
  const auto& getField() const {
    return Fields.template getField<FieldType>();
  }
  // check if field type is in Fields
  template <typename FieldType>
  static constexpr bool hasField() {
    return isTypeInTuple<FieldType, TypePack>::value;
  }

  // template <typename FieldType, unsigned int i>
  // void addField(FieldType& field) {
  //   Fields.template addField<FieldType, i>(field);
  // }
  template <typename FieldType>
  void addField(FieldType& field) {
    Fields.template addField<FieldType>(field);
  }

  std::size_t getNbrId(std::size_t id, int dir) const { return id + Delta_Index[dir]; }

  Block<T, LatSet::d>& getGeo() { return BlockGeo; }
  const Block<T, LatSet::d>& getGeo() const { return BlockGeo; }
  int getNx() const { return BlockGeo.getNx(); }
  int getNy() const { return BlockGeo.getNy(); }
  int getNz() const { return BlockGeo.getNz(); }
  std::size_t getN() const { return BlockGeo.getN(); }
  int getOverlap() const { return BlockGeo.getOverlap(); }
  std::uint8_t getLevel() const { return BlockGeo.getLevel(); }
  const Vector<int, LatSet::d>& getProjection() const { return BlockGeo.getProjection(); }
  const std::array<int, LatSet::q>& getDelta_Index() const { return Delta_Index; }

  template <typename CELLDYNAMICS, typename ArrayType>
  void ApplyCellDynamics(const ArrayType& flagarr) {
    for (std::size_t id = 0; id < getN(); ++id) {
      CellType cell(id, *this);
      CELLDYNAMICS::Execute(flagarr[id], cell);
    }
  }

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics() {
    for (std::size_t id = 0; id < getN(); ++id) {
      CellType cell(id, *this);
      CELLDYNAMICS::apply(cell);
    }
  }

  template <typename CELLDYNAMICS, typename ArrayType>
  void ApplyInnerCellDynamics(const ArrayType& flagarr) {
    if constexpr(LatSet::d == 2){
    for (int j = getOverlap(); j < getNy() - getOverlap(); ++j) {
      std::size_t id = j * getNx() + getOverlap();
      for (int i = getOverlap(); i < getNx() - getOverlap(); ++i) {
        CellType cell(id, *this);
        CELLDYNAMICS::Execute(flagarr[id], cell);
        ++id;
      }
    }
    } else if constexpr(LatSet::d == 3){
      for (int k = getOverlap(); k < getNz() - getOverlap(); ++k) {
        for (int j = getOverlap(); j < getNy() - getOverlap(); ++j) {
          std::size_t id = k * getNx() * getNy() + j * getNx() + getOverlap();
          for (int i = getOverlap(); i < getNx() - getOverlap(); ++i) {
            CellType cell(id, *this);
            CELLDYNAMICS::Execute(flagarr[id], cell);
            ++id;
          }
        }
      }
    }
  }

  template <typename CELLDYNAMICS>
  void ApplyInnerCellDynamics() {
    if constexpr(LatSet::d == 2){
    for (int j = getOverlap(); j < getNy() - getOverlap(); ++j) {
      std::size_t id = j * getNx() + getOverlap();
      for (int i = getOverlap(); i < getNx() - getOverlap(); ++i) {
        CellType cell(id, *this);
        CELLDYNAMICS::apply(cell);
        ++id;
      }
    }
    } else if constexpr(LatSet::d == 3){
      for (int k = getOverlap(); k < getNz() - getOverlap(); ++k) {
        for (int j = getOverlap(); j < getNy() - getOverlap(); ++j) {
          std::size_t id = k * getNx() * getNy() + j * getNx() + getOverlap();
          for (int i = getOverlap(); i < getNx() - getOverlap(); ++i) {
            CellType cell(id, *this);
            CELLDYNAMICS::apply(cell);
            ++id;
          }
        }
      }
    }
  }

};

// block lattice manager base
template <typename T, typename LatSet, typename TypePack>
class BlockLatticeManagerBase {
 public:
  using FIELDS = typename ExtractFieldPack<TypePack>::pack1;
  using FIELDPTRS = typename ExtractFieldPack<TypePack>::pack2;
  using ALLFIELDS = typename ExtractFieldPack<TypePack>::mergedpack;

  using CellType = GenericCell<T, LatSet, ALLFIELDS>;

  using FloatType = T;

 protected:
  BlockGeometry<T, LatSet::d>& BlockGeo;

  // field
  BlockFieldManagerCollection<T, LatSet, FIELDS> Fields;
  BlockFieldManagerPtrCollection<T, LatSet, FIELDPTRS> FieldPtrs;

 public:
  template <typename... FIELDPTRTYPES>
  BlockLatticeManagerBase(BlockGeometry<T, LatSet::d>& blockgeo,
                          FIELDPTRTYPES*... fieldptrs)
      : BlockGeo(blockgeo), Fields(blockgeo), FieldPtrs(fieldptrs...) {}

  template <typename INITVALUEPACK, typename... FIELDPTRTYPES>
  BlockLatticeManagerBase(BlockGeometry<T, LatSet::d>& blockgeo,
                          INITVALUEPACK& initvalues, FIELDPTRTYPES*... fieldptrs)
      : BlockGeo(blockgeo), Fields(blockgeo, initvalues.values), FieldPtrs(fieldptrs...) {
  }

  template <typename FieldType>
  auto& getField() {
    if constexpr (isTypeInTuple<FieldType, FIELDS>::value) {
      return Fields.template getField<FieldType>();
    } else if constexpr (isTypeInTuple<FieldType, FIELDPTRS>::value) {
      return FieldPtrs.template getField<FieldType>();
    }
  }
  template <typename FieldType>
  const auto& getField() const {
    if constexpr (isTypeInTuple<FieldType, FIELDS>::value) {
      return Fields.template getField<FieldType>();
    } else if constexpr (isTypeInTuple<FieldType, FIELDPTRS>::value) {
      return FieldPtrs.template getField<FieldType>();
    }
  }
  // check if the field type is in Fields or FieldPtrs
  template <typename FieldType>
  static constexpr bool hasField() {
    return isTypeInTuple<FieldType, FIELDS>::value ||
           isTypeInTuple<FieldType, FIELDPTRS>::value;
  }

  inline std::uint8_t getMaxLevel() const { return BlockGeo.getMaxLevel(); }
  BlockGeometry<T, LatSet::d>& getGeo() { return BlockGeo; }
  const BlockGeometry<T, LatSet::d>& getGeo() const { return BlockGeo; }
};