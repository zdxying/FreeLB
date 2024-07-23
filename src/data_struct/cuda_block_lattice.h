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

// Cell
#include "data_struct/cell.h"
// Block and BlockGeometry
#include "utils/alias.h"


namespace cudev {

template <typename Pack>
class FieldPtrCollection;

template <typename... Fields>
class FieldPtrCollection<TypePack<Fields...>> {
 public:
  __any__ FieldPtrCollection(std::tuple<Fields*...> fieldrefs) : fields(fieldrefs) {}

  void Init(Fields&... fieldrefs) {
    fields = std::make_tuple(&fieldrefs...);
    // fields = make_Tuple(&fieldrefs...);
  }

  template <typename FieldType>
  __device__ FieldType& getField() {
    return *(std::get<FieldType*>(fields));
    // return *(fields.template get<FieldType*>());
  }
  template <typename FieldType>
  __device__ const FieldType& getField() const {
    return *(std::get<FieldType*>(fields));
    // return *(fields.template get<FieldType*>());
  }

  // assign a field pointer to the ith position
  // template <typename FieldType, unsigned int i>
  // void addField(FieldType& field) {
  //   std::get<i>(fields) = &field;
  // }
  template <typename FieldType>
  __device__ void addField(FieldType& field) {
    std::get<FieldType*>(fields) = &field;
    // fields.template get<FieldType*>() = &field;
  }

 private:
  std::tuple<Fields*...> fields;
  // Tuple<Fields*...> fields;
};


template <typename T, typename LatSet, typename TypePack>
class BlockLatticeBase {
 public:
  using CellType = GenericCell<T, LatSet, TypePack>;
  using LatticeSet = LatSet;
  using FloatType = T;

  using GenericRho = typename FindGenericRhoType<T, TypePack>::type;

 protected:
  // nbr index
  int* Delta_Index;
  // field
  FieldPtrCollection<TypePack>* Fields;

 public:
  __any__ BlockLatticeBase(int* delta_index, FieldPtrCollection<TypePack>* fields)
      : Delta_Index(delta_index), Fields(fields) {}

  template <typename FieldType>
  __device__ auto& getField() {
    return Fields->template getField<FieldType>();
  }
  template <typename FieldType>
  __device__ const auto& getField() const {
    return Fields->template getField<FieldType>();
  }
  // check if field type is in Fields
  template <typename FieldType>
  __device__ static constexpr bool hasField() {
    return isTypeInTuple<FieldType, TypePack>::value;
  }

  __device__ std::size_t getNbrId(std::size_t id, int dir) const {
    return id + Delta_Index[dir];
  }
  __device__ const int* getDelta_Index() const { return Delta_Index; }
};


template <typename T, typename LatSet, typename TypePack>
class BlockLattice : public BlockLatticeBase<T, LatSet, TypePack> {
 public:
  using CellType = Cell<T, LatSet, TypePack>;
  using LatticeSet = LatSet;
  using FloatType = T;

  using GenericRho = typename FindGenericRhoType<T, TypePack>::type;

 protected:
  // --- lattice communication structure ---

  // omega = 1 / tau
  T* Omega;
  // 1 - omega
  T* _Omega;
  // 1 - omega/2
  T* fOmega;

  // tolerance
  // T* RhoRes;
  // GenericArray<T> RhoOld;
  // T URes;
  // GenericArray<Vector<T, LatSet::d>> UOld;

 public:
  template <typename... FIELDPTRS>
  __any__ BlockLattice(int* delta_index, FieldPtrCollection<TypePack>* fields, T* omega,
                       T* _omega, T* fomega)
      : BlockLatticeBase<T, LatSet, TypePack>(delta_index, fields), Omega(omega),
        _Omega(_omega), fOmega(fomega) {}

  __device__ inline T getOmega() const { return *Omega; }
  __device__ inline T get_Omega() const { return *_Omega; }
  __device__ inline T getfOmega() const { return *fOmega; }

  __device__ void Stream() {
    for (unsigned int i = 1; i < LatSet::q; ++i) {
      this->template getField<POP<T, LatSet::q>>().getField(i).rotate(this->Delta_Index[i]);
    }
  }

};

}  // namespace cudev

#endif
