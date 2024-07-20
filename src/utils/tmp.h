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

#include <type_traits>
#include <tuple>
#include <utility>

// #include "head.h"

// template meta-programming library
namespace util {
bool isFlag(std::uint8_t flag1, std::uint8_t flag2);
bool isFlag(std::uint16_t flag1, std::uint16_t flag2);
}  // namespace util

template <typename Type, typename... Types>
struct isOneOf {
  static constexpr bool value = (std::is_same_v<Type, Types> || ...);
};


// fixed-size collection of heterogeneous values, simple implementation of tuple
// Primary template for Tuple
template <typename... Types>
struct Tuple;

// Specialization for Tuple with at least one type
template <typename T, typename... Types>
struct Tuple<T, Types...> {

  using Head = T;
  using Tail = Tuple<Types...>;
  static constexpr std::size_t size = 1 + sizeof...(Types);

  constexpr Tuple() = default;

  __any__ constexpr Tuple(const T& head, const Types&... tail) : _head(head), _tail(tail...) {}

  __any__ constexpr Tuple(T&& head, Types&&... tail)
    : _head(std::forward<T>(head)), _tail(std::forward<Types>(tail)...) {}

  template <typename... Ts, typename... Us>
  __any__ constexpr Tuple(const Tuple<Ts...>& t1, const Tuple<Us...>& t2) : _head(t1._head), _tail(t1._tail, t2) {}
  template <typename... Ts>
  __any__ constexpr Tuple(const Tuple<> &t1, const Tuple<Ts...>& t2) : _head(t2._head), _tail(t2._tail) {}

  __any__ constexpr Tuple(const Tuple &other) : _head(other._head), _tail(other._tail) {}
  __any__ constexpr Tuple(Tuple &&other) : _head(std::move(other._head)), _tail(std::move(other._tail)) {}
  __any__ constexpr Tuple &operator=(const Tuple &other) {
    _head = other._head;
    _tail = other._tail;
    return *this;
  }
  __any__ constexpr Tuple &operator=(Tuple &&other) {
    _head = std::move(other._head);
    _tail = std::move(other._tail);
    return *this;
  }

  __any__ constexpr Tail& getTail() { return _tail; }
  __any__ constexpr const Tail& getTail() const { return _tail; }

  // get ith element
  template <std::size_t N>
  __any__ constexpr auto &get() {
    if constexpr (N == 0) {
      return _head;
    } else {
      return _tail.template get<N - 1>();
    }
  }
  template <std::size_t N>
  __any__ constexpr const auto &get() const {
    if constexpr (N == 0) {
      return _head;
    } else {
      return _tail.template get<N - 1>();
    }
  }

  // get by type
  template <typename U>
  __any__ constexpr auto &get() {
    static_assert(!isOneOf<T, Types...>::value, "[Tuple]: Duplicate type in Tuple");
    if constexpr (std::is_same_v<T, U>) {
      return _head;
    } else {
      return _tail.template get<U>();
    }
  }
  template <typename U>
  __any__ constexpr const auto &get() const {
    static_assert(!isOneOf<T, Types...>::value, "[Tuple]: Duplicate type in Tuple");
    if constexpr (std::is_same_v<T, U>) {
      return _head;
    } else {
      return _tail.template get<U>();
    }
  }

  Head _head;
  Tail _tail;
};

// Specialization for an empty Tuple
template <>
struct Tuple<> {
  static constexpr std::size_t size = 0;

  template <std::size_t N>
  __any__ void get() const {
    static_assert(N != N, "[Tuple]: Attempt to access element in an empty Tuple");
  }

  template <typename U>
  __any__ void get() const {
    static_assert(sizeof(U) == 0, "[Tuple]: Attempt to access element in an empty Tuple");
  }
};

// helper function to create Tuple
template <typename... Types>
__any__ constexpr auto make_Tuple(Types... args) {
  return Tuple<Types...>(args...);
}
// handle empty Tuple
__any__ constexpr auto make_Tuple() {
  return Tuple<>();
}


// helper function to merge Tuple
// Base case
template<typename... Types>
__any__ constexpr Tuple<Types...> Tuple_cat(Tuple<Types...> t) {
  return t;
}
template <typename... Types1, typename... Types2>
__any__ constexpr Tuple<Types1..., Types2...> Tuple_cat(Tuple<Types1...> t1, Tuple<Types2...> t2) {
  return Tuple<Types1..., Types2...>(t1, t2);
}
template <typename... Types1, typename... Types2, typename... Types3, typename... Rest>
__any__ constexpr Tuple<Types1..., Types2..., Types3...> Tuple_cat(Tuple<Types1...> t1, Tuple<Types2...> t2,
                                                 Tuple<Types3...> t3, Rest... rest) {
  return Tuple_cat(Tuple_cat(t1, t2), t3, rest...);
}


namespace tmp {
// a sequence of integers
template <auto I, auto... Is>
struct int_sequence {
  static constexpr auto first = I;
  using rest = int_sequence<Is...>;
  static constexpr std::size_t rest_size = sizeof...(Is);
};
template <auto I>
struct int_sequence<I> {
  static constexpr auto first = I;
  static constexpr std::size_t rest_size = 0;
};

// need c++20
// lambda-expression in template-argument only available with ‘-std=c++20’ or
// ‘-std=gnu++20’
template <typename Func, std::size_t... Indices>
auto make_int_sequence_impl(Func &&f, std::index_sequence<Indices...>) {
  return int_sequence<f(Indices)...>{};
}
template <std::size_t size, typename Func>
auto make_int_sequence(Func &&f) {
  return make_int_sequence_impl(std::forward<Func>(f), std::make_index_sequence<size>{});
}

// key-value pair
template <auto KeyValue, typename Type>
struct Key_TypePair {
  // key is typically an int value (uint8_t, uint16_t, etc.)
  static constexpr auto key = KeyValue;
  // type is the type of the struct
  using key_type = decltype(KeyValue);
  using type = Type;
};
// a wrapper class for tuple
template <typename... Types>
struct TupleWrapper {
 private:
  template <auto KeyValue, std::size_t I, std::size_t... Is>
  __any__ static auto get_by_key(std::index_sequence<I, Is...>) {
    if constexpr (get_by_index<I>::key == KeyValue)
      return get_by_index<I>();
    else if constexpr (sizeof...(Is) > 0)
      return get_by_key<KeyValue>(std::index_sequence<Is...>{});
  }

 public:
  using tuple = std::tuple<Types...>;
  // access ith element in tuple
  template <std::size_t idx>
  using get_by_index = std::tuple_element_t<idx, tuple>;
  // access element by key
  template <auto KeyValue>
  using get =
    decltype(get_by_key<KeyValue>(std::make_index_sequence<sizeof...(Types)>{}));
  // size
  static constexpr std::size_t size = sizeof...(Types);
  // key sequence
  using make_key_sequence = int_sequence<Types::key...>;
};

// helper struct to achieve if-else structure
template <typename TUPLE, typename FlagType, typename CELL, typename KeySequence>
struct SelectTask {
  __any__ static void execute(FlagType flag, CELL &cell) {
    if (util::isFlag(flag, KeySequence::first)) {
      TUPLE::template get<KeySequence::first>::type::apply(cell);
    } else if constexpr (KeySequence::rest_size > 0) {
      SelectTask<TUPLE, FlagType, CELL, typename KeySequence::rest>::execute(flag, cell);
    }
  }
};

template <typename TUPLE, typename FlagType, typename CELL>
struct TaskSelector {
  __any__ static void Execute(FlagType flag, CELL &cell) {
    SelectTask<TUPLE, FlagType, CELL, typename TUPLE::make_key_sequence>::execute(flag,
                                                                                  cell);
  }
};

// helper struct to achieve if-else structure for coupled tasks
template <typename TUPLE, typename FlagType, typename CELL0, typename CELL1,
          typename KeySequence>
struct SelectCoupledTask {
  __any__ static void execute(FlagType flag, CELL0 &cell0, CELL1 &cell1) {
    if (util::isFlag(flag, KeySequence::first)) {
      TUPLE::template get<KeySequence::first>::type::apply(cell0, cell1);
    } else if constexpr (KeySequence::rest_size > 0) {
      SelectCoupledTask<TUPLE, FlagType, CELL0, CELL1,
                        typename KeySequence::rest>::execute(flag, cell0, cell1);
    }
  }
};

template <typename TUPLE, typename FlagType, typename CELL0, typename CELL1>
struct CoupledTaskSelector {
  __any__ static void Execute(FlagType flag, CELL0 &cell0, CELL1 &cell1) {
    SelectCoupledTask<TUPLE, FlagType, CELL0, CELL1,
                      typename TUPLE::make_key_sequence>::execute(flag, cell0, cell1);
  }
};


}  // namespace tmp

// TASK SELECTOR
template <typename FlagType, typename CELL, typename FirstTask, typename... RestTasks>
struct SelectTask {
  __any__ static void execute(FlagType flag, CELL &cell) {
    if (util::isFlag(flag, FirstTask::key)) {
      FirstTask::type::apply(cell);
    } else if constexpr (sizeof...(RestTasks) > 0) {
      SelectTask<FlagType, CELL, RestTasks...>::execute(flag, cell);
    }
  }
};

template <typename FlagType, typename CELL, typename FirstTask>
struct SelectTask<FlagType, CELL, FirstTask>{
  __any__ static void execute(FlagType flag, CELL &cell) {
    if (util::isFlag(flag, FirstTask::key)) {
      FirstTask::type::apply(cell);
    }
  }
};

template <typename FlagType, typename CELL, typename... Tasks>
struct TaskSelector {
  __any__ static void Execute(FlagType flag, CELL &cell) {
    SelectTask<FlagType, CELL, Tasks...>::execute(flag, cell);
  }
};

// COUPLED TASK SELECTOR
template <typename FlagType, typename CELL0, typename CELL1, typename FirstTask, typename... RestTasks>
struct SelectCoupledTask {
  __any__ static void execute(FlagType flag, CELL0 &cell0, CELL1 &cell1) {
    if (util::isFlag(flag, FirstTask::key)) {
      FirstTask::type::apply(cell0, cell1);
    } else if constexpr (sizeof...(RestTasks) > 0) {
      SelectCoupledTask<FlagType, CELL0, CELL1, RestTasks...>::execute(flag, cell0, cell1);
    }
  }
};

template <typename FlagType, typename CELL0, typename CELL1, typename FirstTask>
struct SelectCoupledTask<FlagType, CELL0, CELL1, FirstTask>{
  __any__ static void execute(FlagType flag, CELL0 &cell0, CELL1 &cell1) {
    if (util::isFlag(flag, FirstTask::key)) {
      FirstTask::type::apply(cell0, cell1);
    }
  }
};

template <typename FlagType, typename CELL0, typename CELL1, typename... Tasks>
struct CoupledTaskSelector {
  __any__ static void Execute(FlagType flag, CELL0 &cell0, CELL1 &cell1) {
    SelectCoupledTask<FlagType, CELL0, CELL1, Tasks...>::execute(flag, cell0, cell1);
  }
};

template <typename... Parameter>
struct TypePack {
  using types = std::tuple<Parameter...>;
  // using types = Tuple<Parameter...>;
  static constexpr std::size_t size = sizeof...(Parameter);
};

template <typename... Params1, typename... Params2>
struct TypePack<TypePack<Params1...>, TypePack<Params2...>> {
  using types = std::tuple<Params1..., Params2...>;
  // using types = Tuple<Params1..., Params2...>;
  static constexpr std::size_t size = sizeof...(Params1) + sizeof...(Params2);
};


template <typename Type, typename Pack>
struct isTypeInTuple {
  static constexpr bool value = false;
};

// If Type is the same as any type in Parameter..., value is true
template <typename Type, typename... Parameter>
struct isTypeInTuple<Type, TypePack<Parameter...>> {
  static constexpr bool value = (std::is_same_v<Type, Parameter> || ...);
};

template <typename Type, typename First, typename... Rest>
struct is_same_at_index {
  static constexpr unsigned int value =
    std::is_same_v<Type, First> ? 0 : 1 + is_same_at_index<Type, Rest...>::value;
};

template <typename Type, typename First>
struct is_same_at_index<Type, First>{
  static constexpr unsigned int value = 0;
};


template <typename Pack>
struct ExtractFieldPack;

template <typename... Params1, typename... Params2>
struct ExtractFieldPack<TypePack<TypePack<Params1...>, TypePack<Params2...>>> {
  using pack1 = TypePack<Params1...>;
  using pack2 = TypePack<Params2...>;
  using mergedpack = TypePack<Params1..., Params2...>;
};

template <typename... Params1>
struct ExtractFieldPack<TypePack<Params1...>> {
  using pack1 = TypePack<Params1...>;
  using pack2 = TypePack<>;
  using mergedpack = TypePack<Params1...>;
};

// extract cudev field pack from TypePack
template <typename Pack>
struct ExtractCudevFieldPack;

template <typename... Params1>
struct ExtractCudevFieldPack<TypePack<Params1...>> {
  using cudev_pack = TypePack<typename Params1::cudev_FieldType...>;
};

template <typename... Packs>
struct MergeFieldPack;

template <typename... Params1, typename... Params2>
struct MergeFieldPack<TypePack<Params1...>, TypePack<Params2...>> {
  using mergedpack = TypePack<Params1..., Params2...>;
};

template <typename... Params1, typename... Params2, typename... Params3>
struct MergeFieldPack<TypePack<Params1...>, TypePack<Params2...>, TypePack<Params3...>> {
  using mergedpack = TypePack<Params1..., Params2..., Params3...>;
};


template <typename... Parameter>
struct ValuePack {
  std::tuple<Parameter...> values;
  // Tuple<Parameter...> values;

  ValuePack(Parameter... value) : values(value...) {}
  ValuePack(const std::tuple<Parameter...>& value) : values(value) {}
  // ValuePack(const Tuple<Parameter...>& value) : values(value) {}
};

// merge two ValuePack
template <typename... Params1, typename... Params2>
ValuePack<Params1..., Params2...> mergeValuePack(ValuePack<Params1...> &pack1,
                                                 ValuePack<Params2...> &pack2) {
  return ValuePack<Params1..., Params2...>(std::tuple_cat(pack1.values, pack2.values));
  // return ValuePack<Params1..., Params2...>(Tuple_cat(pack1.values, pack2.values));
}
template <typename... Params1, typename... Params2, typename... Params3>
ValuePack<Params1..., Params2..., Params3...> mergeValuePack(ValuePack<Params1...> &pack1,
                                                             ValuePack<Params2...> &pack2,
                                                             ValuePack<Params3...> &pack3) {
  return ValuePack<Params1..., Params2..., Params3...>(std::tuple_cat(pack1.values, pack2.values, pack3.values));
  // return ValuePack<Params1..., Params2..., Params3...>(Tuple_cat(pack1.values, pack2.values, pack3.values));
}


template <typename T, typename Pack>
struct GetGenericRhoType {
  using type = std::conditional_t<
    isTypeInTuple<RHO<T>, Pack>::value, RHO<T>,
    std::conditional_t<
      isTypeInTuple<TEMP<T>, Pack>::value, TEMP<T>,
      std::conditional_t<isTypeInTuple<CONC<T>, Pack>::value, CONC<T>, void>>>;
};

template <typename T, typename Pack>
struct FindGenericRhoType;

template <typename T, typename First, typename... Rest>
struct FindGenericRhoType<T, TypePack<First, Rest...>> {
  using type =
    std::conditional_t<isOneOf<First, RHO<T>, TEMP<T>, CONC<T>>::value, First,
                       typename FindGenericRhoType<T, TypePack<Rest...>>::type>;
};

template <typename T, typename First>
struct FindGenericRhoType<T, TypePack<First>> {
  using type =
    std::conditional_t<isOneOf<First, RHO<T>, TEMP<T>, CONC<T>>::value, First, void>;
};

namespace cudev{

#ifdef __CUDACC__

template <typename T, typename Pack>
struct FindGenericRhoType;

template <typename T, typename First, typename... Rest>
struct FindGenericRhoType<T, TypePack<First, Rest...>> {
  using type =
    std::conditional_t<isOneOf<First, RHO<T>, TEMP<T>, CONC<T>>::value, First,
                       typename FindGenericRhoType<T, TypePack<Rest...>>::type>;
};

template <typename T, typename First>
struct FindGenericRhoType<T, TypePack<First>> {
  using type =
    std::conditional_t<isOneOf<First, RHO<T>, TEMP<T>, CONC<T>>::value, First, void>;
};

#endif

}  // namespace cudev


// unroll a for loop, use: unroll_for<x, y>([&](unsigned int i) { ... });
template <unsigned int start, unsigned int end, typename Func>
void unroll_for(Func &&f) {
  if constexpr (start < end) {
    f(start);
    unroll_for<start + 1, end>(std::forward<Func>(f));
  }
}
