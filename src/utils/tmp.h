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
#include <utility>

#include "head.h"

// template meta-programming library
namespace util {
bool isFlag(std::uint8_t flag1, std::uint8_t flag2);
bool isFlag(std::uint16_t flag1, std::uint16_t flag2);
}  // namespace util

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
  static auto get_by_key(std::index_sequence<I, Is...>) {
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
  static void execute(FlagType flag, CELL &cell) {
    if (util::isFlag(flag, KeySequence::first)) {
      TUPLE::template get<KeySequence::first>::type::apply(cell);
    } else if constexpr (KeySequence::rest_size > 0) {
      SelectTask<TUPLE, FlagType, CELL, typename KeySequence::rest>::execute(flag, cell);
    }
  }
};

template <typename TUPLE, typename FlagType, typename CELL>
struct TaskSelector {
  static void Execute(FlagType flag, CELL &cell) {
    SelectTask<TUPLE, FlagType, CELL, typename TUPLE::make_key_sequence>::execute(flag,
                                                                                  cell);
  }
};

// helper struct to achieve if-else structure for coupled tasks
template <typename TUPLE, typename FlagType, typename CELL0, typename CELL1,
          typename KeySequence>
struct SelectCoupledTask {
  static void execute(FlagType flag, CELL0 &cell0, CELL1 &cell1) {
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
  static void Execute(FlagType flag, CELL0 &cell0, CELL1 &cell1) {
    SelectCoupledTask<TUPLE, FlagType, CELL0, CELL1,
                      typename TUPLE::make_key_sequence>::execute(flag, cell0, cell1);
  }
};


}  // namespace tmp

// TASK SELECTOR
template <typename FlagType, typename CELL, typename FirstTask, typename... RestTasks>
struct SelectTask {
  static void execute(FlagType flag, CELL &cell) {
    if (util::isFlag(flag, FirstTask::key)) {
      FirstTask::type::apply(cell);
    } else if constexpr (sizeof...(RestTasks) > 0) {
      SelectTask<FlagType, CELL, RestTasks...>::execute(flag, cell);
    }
  }
};

template <typename FlagType, typename CELL, typename FirstTask>
struct SelectTask<FlagType, CELL, FirstTask>{
  static void execute(FlagType flag, CELL &cell) {
    if (util::isFlag(flag, FirstTask::key)) {
      FirstTask::type::apply(cell);
    }
  }
};

template <typename FlagType, typename CELL, typename... Tasks>
struct TaskSelector {
  static void Execute(FlagType flag, CELL &cell) {
    SelectTask<FlagType, CELL, Tasks...>::execute(flag, cell);
  }
};

// COUPLED TASK SELECTOR
template <typename FlagType, typename CELL0, typename CELL1, typename FirstTask, typename... RestTasks>
struct SelectCoupledTask {
  static void execute(FlagType flag, CELL0 &cell0, CELL1 &cell1) {
    if (util::isFlag(flag, FirstTask::key)) {
      FirstTask::type::apply(cell0, cell1);
    } else if constexpr (sizeof...(RestTasks) > 0) {
      SelectCoupledTask<FlagType, CELL0, CELL1, RestTasks...>::execute(flag, cell0, cell1);
    }
  }
};

template <typename FlagType, typename CELL0, typename CELL1, typename FirstTask>
struct SelectCoupledTask<FlagType, CELL0, CELL1, FirstTask>{
  static void execute(FlagType flag, CELL0 &cell0, CELL1 &cell1) {
    if (util::isFlag(flag, FirstTask::key)) {
      FirstTask::type::apply(cell0, cell1);
    }
  }
};

template <typename FlagType, typename CELL0, typename CELL1, typename... Tasks>
struct CoupledTaskSelector {
  static void Execute(FlagType flag, CELL0 &cell0, CELL1 &cell1) {
    SelectCoupledTask<FlagType, CELL0, CELL1, Tasks...>::execute(flag, cell0, cell1);
  }
};

template <typename... Parameter>
struct TypePack {
  using types = std::tuple<Parameter...>;
  static constexpr std::size_t size = sizeof...(Parameter);
};

template <typename... Params1, typename... Params2>
struct TypePack<TypePack<Params1...>, TypePack<Params2...>> {
  using types = std::tuple<Params1..., Params2...>;
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

template <typename Type, typename... Types>
struct isOneOf {
  static constexpr bool value = (std::is_same_v<Type, Types> || ...);
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

template <typename... Parameter>
struct ValuePack {
  std::tuple<Parameter...> values;

  ValuePack(Parameter... value) : values(value...) {}
};


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
