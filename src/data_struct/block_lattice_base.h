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

#ifdef __CUDACC__
#include "data_struct/cuda_block_lattice.h"
#endif

template <typename Pack>
class FieldPtrCollection;

template <typename... Fields>
class FieldPtrCollection<TypePack<Fields...>> {
 public:
  FieldPtrCollection(Fields&... fieldrefs) : fields(&fieldrefs...) {}
  FieldPtrCollection(std::tuple<Fields*...> fieldrefs) : fields(fieldrefs) {}

  void Init(Fields&... fieldrefs) {
    fields = std::make_tuple(&fieldrefs...);
    // fields = make_Tuple(&fieldrefs...);
  }

  template <typename FieldType>
  FieldType& getField() {
    return *(std::get<FieldType*>(fields));
    // return *(fields.template get<FieldType*>());
  }
  template <typename FieldType>
  const FieldType& getField() const {
    return *(std::get<FieldType*>(fields));
    // return *(fields.template get<FieldType*>());
  }

  // assign a field pointer to the ith position
  // template <typename FieldType, unsigned int i>
  // void addField(FieldType& field) {
  //   std::get<i>(fields) = &field;
  // }
  template <typename FieldType>
  void addField(FieldType& field) {
    std::get<FieldType*>(fields) = &field;
    // fields.template get<FieldType*>() = &field;
  }

#ifdef __CUDACC__
  using cudev_TypePack = typename ExtractCudevFieldPack<TypePack<Fields...>>::cudev_pack;

  template <std::size_t... Is>
  auto get_dev_FPtrCol(std::index_sequence<Is...>) {
    return std::make_tuple((std::get<Is>(fields)->get_devObj())...);
  }
  // get cuda device version of FieldPtrCollection
  auto get_dev_FPtrCol() {
    // construct a tuple of device pointers
    auto tup = get_dev_FPtrCol(std::make_index_sequence<sizeof...(Fields)>());
    return cudev::FieldPtrCollection<cudev_TypePack>(tup);
  }
#endif

 private:
  std::tuple<Fields*...> fields;
  // Tuple<Fields*...> fields;
};


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

#ifdef __CUDACC__
  using cudev_TypePack = typename ExtractCudevFieldPack<TypePack>::cudev_pack;
  using cudev_BlockLatticeBase = cudev::BlockLatticeBase<T, LatSet, cudev_TypePack>;
#endif

 protected:
  // nbr index
  std::array<int, LatSet::q> Delta_Index;
  // geometry
  Block<T, LatSet::d>& BlockGeo;
  // field
  FieldPtrCollection<TypePack> Fields;

#ifdef __CUDACC__
  int* dev_Delta_Index;
  cudev::FieldPtrCollection<cudev_TypePack>* dev_Fields;
  cudev_BlockLatticeBase* dev_this;
#endif

 public:
  template <typename... FIELDPTRS>
  BlockLatticeBase(Block<T, LatSet::d>& block, std::tuple<FIELDPTRS...> fieldptrs)
      : BlockGeo(block), Fields(fieldptrs) {
    Delta_Index =
      make_Array<int, LatSet::q>([&](int i) { return latset::c<LatSet>(i) * getProjection(); });
#ifdef __CUDACC__
    InitDeviceData();
#endif
  }

  ~BlockLatticeBase() {
#ifdef __CUDACC__
    if (dev_Delta_Index) cuda_free(dev_Delta_Index);
    if (dev_Fields) cuda_free(dev_Fields);
    if (dev_this) cuda_free(dev_this);
#endif
  }

#ifdef __CUDACC__
  void InitDeviceData() {
    dev_Delta_Index = cuda_malloc<int>(LatSet::q);
    host_to_device(dev_Delta_Index, Delta_Index.data(), LatSet::q);

    dev_Fields = cuda_malloc<cudev::FieldPtrCollection<cudev_TypePack>>(1);
    cudev::FieldPtrCollection<cudev_TypePack> host_Data = Fields.get_dev_FPtrCol();
    host_to_device(dev_Fields, &host_Data, 1);
    constructInDevice();
  }
  void constructInDevice() {
    dev_this = cuda_malloc<cudev_BlockLatticeBase>(1);
    cudev_BlockLatticeBase temp(dev_Delta_Index, dev_Fields);
    host_to_device(dev_this, &temp, 1);
  }
  cudev_BlockLatticeBase* get_devObj() { return dev_this; }

#endif

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

  Block<T, LatSet::d>& getBlock() { return BlockGeo; }
  const Block<T, LatSet::d>& getBlock() const { return BlockGeo; }
  int getNx() const { return BlockGeo.getNx(); }
  int getNy() const { return BlockGeo.getNy(); }
  int getNz() const { return BlockGeo.getNz(); }
  std::size_t getN() const { return BlockGeo.getN(); }
  std::size_t getVoxNum() const { return BlockGeo.getVoxNum(); }
  int getOverlap() const { return BlockGeo.getOverlap(); }
  std::uint8_t getLevel() const { return BlockGeo.getLevel(); }
  const Vector<int, LatSet::d>& getProjection() const { return BlockGeo.getProjection(); }
  const std::array<int, LatSet::q>& getDelta_Index() const { return Delta_Index; }

};


template <typename T, typename LatSet, typename Pack>
class BlockFieldManagerCollection;

template <typename T, typename LatSet, typename... Fields>
class BlockFieldManagerCollection<T, LatSet, TypePack<Fields...>> {
 public:
  static constexpr std::size_t FieldNum = sizeof...(Fields);

  BlockFieldManagerCollection() : fields(std::make_tuple(Fields()...)) {}
  // : fields(make_Tuple(Fields()...)) {}
  template <typename BLOCKGEOMETRY>
  BlockFieldManagerCollection(BLOCKGEOMETRY& blockgeo)
      : fields(std::make_tuple(BlockFieldManager<Fields, T, LatSet::d>(blockgeo)...)) {}
  // : fields(make_Tuple(BlockFieldManager<Fields, T, LatSet::d>(blockgeo)...)) {}

  template <typename BLOCKGEOMETRY, typename... InitValues>
  BlockFieldManagerCollection(BLOCKGEOMETRY& blockgeo,
                              const std::tuple<InitValues...>& initvalues)
      : fields(make_fields(blockgeo, initvalues, std::make_index_sequence<FieldNum>())) {}

  ~BlockFieldManagerCollection() = default;

  template <typename FieldType>
  auto& getField() {
    return std::get<BlockFieldManager<FieldType, T, LatSet::d>>(fields);
    // return fields.template get<BlockFieldManager<FieldType, T, LatSet::d>>();
  }
  template <typename FieldType>
  const auto& getField() const {
    return std::get<BlockFieldManager<FieldType, T, LatSet::d>>(fields);
    // return fields.template get<BlockFieldManager<FieldType, T, LatSet::d>>();
  }

  auto& getFields() { return fields; }
  template <std::size_t... Is>
  auto get_ith(std::index_sequence<Is...>, int i) {
    return std::make_tuple(&(std::get<Is>(fields).getBlockField(i))...);
    // return make_Tuple(&(fields.template get<Is>().getBlockField(i).getFieldType())...);
  }
  auto get_ith(int i) { return get_ith(std::make_index_sequence<FieldNum>(), i); }

  // for each field, call func
  template <typename Func>
  void forEachField(Func&& func) {
    forEachFieldImpl(std::forward<Func>(func), std::make_index_sequence<FieldNum>());
  }
  template <typename Func, std::size_t... Is>
  void forEachFieldImpl(Func&& func, std::index_sequence<Is...>) {
    (func(std::get<Is>(fields), std::integral_constant<std::size_t, Is>{}), ...);
    // (func(fields.template get<Is>(), std::integral_constant<std::size_t, Is>{}), ...);
  }

 private:
  std::tuple<BlockFieldManager<Fields, T, LatSet::d>...> fields;
  // Tuple<BlockFieldManager<Fields, T, LatSet::d>...> fields;

  template <typename BLOCKGEOMETRY, typename... InitValues, std::size_t... Is>
  auto make_fields(BLOCKGEOMETRY& blockgeo, const std::tuple<InitValues...>& initvalues,
                   std::index_sequence<Is...>) {
    return std::make_tuple(
      BlockFieldManager<Fields, T, LatSet::d>(blockgeo, std::get<Is>(initvalues))...);
    // return make_Tuple(BlockFieldManager<Fields, T, LatSet::d>(blockgeo,
    // initvalues.template get<Is>())...);
  }
};

template <typename T, typename LatSet, typename Pack>
class BlockFieldManagerPtrCollection;

template <typename T, typename LatSet, typename... Fields>
class BlockFieldManagerPtrCollection<T, LatSet, TypePack<Fields...>> {
 public:
  static constexpr std::size_t FieldNum = sizeof...(Fields);

  BlockFieldManagerPtrCollection() : fields() {}
  // fields(make_Tuple((BlockFieldManager<Fields, T, LatSet::d>*)nullptr...)) {}
  template <typename First, typename... Rest>
  BlockFieldManagerPtrCollection(First* first, Rest*... rest) : fields(first, rest...) {}
  ~BlockFieldManagerPtrCollection() = default;

  void Init(Fields&... fieldrefs) {
    fields = std::make_tuple(&fieldrefs...);
    // fields = make_Tuple(&fieldrefs...);
  }

  template <typename FieldType>
  static constexpr unsigned int is_at_index() {
    return is_same_at_index<FieldType, Fields...>::value;
  }
  // not real 'add', just assign the valid field pointer to an existing null pointer
  template <typename FieldType>
  void addField(BlockFieldManager<FieldType, T, LatSet::d>& field) {
    // std::get<is_same_at_index<FieldType, Fields...>::value>(fields) = &field;
    std::get<BlockFieldManager<FieldType, T, LatSet::d>*>(fields) = &field;
    // fields.template get<BlockFieldManager<FieldType, T, LatSet::d>*>() = &field;
  }

  template <typename FieldType>
  auto& getField() {
    return *(std::get<BlockFieldManager<FieldType, T, LatSet::d>*>(fields));
    // return *(fields.template get<BlockFieldManager<FieldType, T, LatSet::d>*>());
  }
  template <typename FieldType>
  const auto& getField() const {
    return *(std::get<BlockFieldManager<FieldType, T, LatSet::d>*>(fields));
    // return *(fields.template get<BlockFieldManager<FieldType, T, LatSet::d>*>());
  }

  auto& getFields() { return fields; }

  template <std::size_t... Is>
  auto get_ith(std::index_sequence<Is...>, int i) {
    // return std::make_tuple(&(std::get<Is>(fields)->getBlockField(i))...);
    return std::make_tuple(
      (std::get<Is>(fields) ? &(std::get<Is>(fields)->getBlockField(i)) : nullptr)...);
    // return make_Tuple((fields.template get<Is>() ? &(fields.template
    // get<Is>()->getBlockField(i).getFieldType()) : nullptr)...);
  }
  auto get_ith(int i) { return get_ith(std::make_index_sequence<FieldNum>(), i); }

 private:
  std::tuple<BlockFieldManager<Fields, T, LatSet::d>*...> fields;
  // Tuple<BlockFieldManager<Fields, T, LatSet::d>*...> fields;
};

template <typename T, typename LatSet, typename Pack>
struct ExtractFieldPtrs {};

template <typename T, typename LatSet, typename... Fields, typename... FieldPtrs>
struct ExtractFieldPtrs<T, LatSet,
                        TypePack<TypePack<Fields...>, TypePack<FieldPtrs...>>> {
  static auto getFieldPtrTuple(
    int i, BlockFieldManagerCollection<T, LatSet, TypePack<Fields...>>& BFMCol,
    BlockFieldManagerPtrCollection<T, LatSet, TypePack<FieldPtrs...>>& BFMPtrCol) {
    return std::tuple_cat(BFMCol.get_ith(i), BFMPtrCol.get_ith(i));
    // return Tuple_cat(BFMCol.get_ith(i), BFMPtrCol.get_ith(i));
  }
};

template <typename T, typename LatSet, typename... Fields>
struct ExtractFieldPtrs<T, LatSet, TypePack<Fields...>> {
  static auto getFieldPtrTuple(
    int i, BlockFieldManagerCollection<T, LatSet, TypePack<Fields...>>& BFMCol) {
    return BFMCol.get_ith(i);
  }
  static auto getFieldPtrTuple(
    int i, BlockFieldManagerCollection<T, LatSet, TypePack<Fields...>>& BFMCol,
    BlockFieldManagerPtrCollection<T, LatSet, TypePack<>>& BFMPtrCol) {
    return BFMCol.get_ith(i);
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