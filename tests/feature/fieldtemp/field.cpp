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

#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;

using LatSet = D2Q9<T>;

struct TESTFIELDBase : public FieldBase<1> {};
template <unsigned int D>
struct TESTFIELD1Base : public FieldBase<D> {};
template <unsigned int Q>
struct TESTFIELD2Base : public FieldBase<Q> {};

// field identifier class
// normal scalar field
struct TESTFIELD{
template <typename LatSet>
using type = GenericField<GenericArray<typename LatSet::FloatType>, TESTFIELDBase>;
};

// vector field in SOA format
struct TESTFIELD1{
template <typename LatSet>
using type = GenericField<GenericArray<typename LatSet::FloatType>, TESTFIELD1Base<LatSet::d>>;
};

// population field
struct TESTFIELD2{
template <typename LatSet>
using type = GenericField<CyclicArray<typename LatSet::FloatType>, TESTFIELD1Base<LatSet::q>>;
};

// block field manager collection
template <typename T, typename LatSet, typename Pack>
class BFMCollection;

template <typename T, typename LatSet, typename... Fields>
class BFMCollection<T, LatSet, TypePack<Fields...>> {
 public:
  static constexpr std::size_t FieldNum = sizeof...(Fields);

  BFMCollection(BlockGeometry<T, LatSet::d>& blockgeo)
      : fields(std::make_tuple(BlockFieldManager<typename Fields::type<LatSet>, T, LatSet::d>(blockgeo)...)) {}

  template <typename... InitValues>
  BFMCollection(BlockGeometry<T, LatSet::d>& blockgeo, const std::tuple<InitValues...>& initvalues)
      : fields(make_fields(blockgeo, initvalues, std::make_index_sequence<FieldNum>())) {}

  ~BFMCollection() = default;

  template <typename FieldType>
  auto& getField() {
    return std::get<BlockFieldManager<typename FieldType::type<LatSet>, T, LatSet::d>>(fields);
  }
  template <typename FieldType>
  const auto& getField() const {
    return std::get<BlockFieldManager<typename FieldType::type<LatSet>, T, LatSet::d>>(fields);
  }

  auto& getFields() { return fields; }
  template <std::size_t... Is>
  auto get_ith(std::index_sequence<Is...>, int i) {
    return std::make_tuple(&(std::get<Is>(fields).getBlockField(i))...);
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
  }

 private:
  std::tuple<BlockFieldManager<typename Fields::type<LatSet>, T, LatSet::d>...> fields;

  template <typename... InitValues, std::size_t... Is>
  auto make_fields(BlockGeometry<T, LatSet::d>& blockgeo, const std::tuple<InitValues...>& initvalues,
                   std::index_sequence<Is...>) {
    return std::make_tuple(BlockFieldManager<typename Fields::type<LatSet>, T, LatSet::d>(blockgeo, std::get<Is>(initvalues))...);
  }
};

int Ni = 100;
int Nj = 100;
int Cell_Len = 1;
int Thread_Num = 4;

int main() {
	AABB<T, 2> cavity(Vector<T, 2>(T(0), T(0)),
										Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
	BlockGeometry2D<T> Geo(Ni, Nj, Thread_Num, cavity, Cell_Len);

	BaseConverter<T> BaseConv(LatSet::cs2);
	BaseConv.SimplifiedConverterFromRT(Ni, T{1}, T{1});

	using FIELDS = TypePack<TESTFIELD1, TESTFIELD2>;
	using CELL = Cell<T, LatSet, FIELDS>;
	ValuePack InitValues(T{2}, T{3});

	BFMCollection<T, LatSet, FIELDS> NSLattice(Geo, InitValues.values);

	// get compiler error
	auto fvar0 = NSLattice.template getField<TESTFIELD>().getBlockField(0).get(0);
	// correct
	auto fvar1 = NSLattice.template getField<TESTFIELD1>().getBlockField(0).template get<1>(0);
	auto fvar2 = NSLattice.template getField<TESTFIELD2>().getBlockField(0).template get<3>(0);

	std::cout << fvar0 << std::endl;
	std::cout << fvar1 << std::endl;
	std::cout << fvar2 << std::endl;


	// BlockLatticeManager<T, LatSet, FIELDS> NSLattice(Geo, InitValues, BaseConv);

}