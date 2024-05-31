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
  using CELL = typename Type::CELL;
  // key is typically an int value (uint8_t, uint16_t, etc.)
  static constexpr auto key = KeyValue;
  // type is the type of the struct
  static auto apply(CELL &cell) { return Type::apply(cell); }
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
      TUPLE::template get<KeySequence::first>::apply(cell);
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

// a simple task executor
template <typename Type>
struct TaskExecutor {
  using CELL = typename Type::CELL;
  static void Execute(CELL &cell) { Type::apply(cell); }
};


}  // namespace tmp


template <typename... Parameter>
struct TypePack {
  using types = std::tuple<Parameter...>;
};

template <typename... Params1, typename... Params2>
struct TypePack<TypePack<Params1...>, TypePack<Params2...>> {
  using types = std::tuple<Params1..., Params2...>;
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


template <typename Pack>
struct ExtractFieldPack;

template <typename... Params1, typename... Params2>
struct ExtractFieldPack<TypePack<TypePack<Params1...>, TypePack<Params2...>>> {
  using pack1 = TypePack<Params1...>;
  using pack2 = TypePack<Params2...>;
  using mergedpack = TypePack<Params1..., Params2...>;
};

template <typename... Parameter>
struct ValuePack {
  std::tuple<Parameter...> values;

  ValuePack(Parameter... value) : values(value...) {}
};