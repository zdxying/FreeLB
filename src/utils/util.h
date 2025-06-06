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

// utility functions
#pragma once

#include <cmath>
#include <omp.h>

// copy
#include <algorithm>
#include <array>
#include <iostream>
// numeric_limits
#include <limits>
#include <vector>
// store any callable object and used in unit converter
#include <functional>

#include "head.h"
#include "utils/tmp.h"


template <typename T, unsigned D, typename F, std::size_t... INDICES>
std::array<T, D> make_array_impl(F &&f, std::index_sequence<INDICES...>) {
  return {(static_cast<void>(INDICES), f())...};
}

template <typename T, unsigned D, typename F>
std::array<T, D> make_array(F &&f) {
  return make_array_impl<T, D>(std::forward<F>(f), std::make_index_sequence<D>{});
}

template <typename T, unsigned D, typename F, std::size_t... INDICES>
std::array<T, D> make_Array_impl(F &&f, std::index_sequence<INDICES...>) {
  return {f(INDICES)...};
}
template <typename T, unsigned D, typename F>
std::array<T, D> make_Array(F &&f) {
  return make_Array_impl<T, D>(std::forward<F>(f), std::make_index_sequence<D>{});
}

template <typename T, unsigned D, typename F, std::size_t... INDICES>
__any__ void make_PlainArray_impl(F &&f, std::index_sequence<INDICES...>) {
  (f(INDICES), ...);
}
template <typename T, unsigned D, typename F>
__any__ void make_PlainArray(F &&f) {
  make_PlainArray_impl<T, D>(std::forward<F>(f), std::make_index_sequence<D>{});
}

// /// construct a std::array<T,D> where T is an object of class
// template <typename T, unsigned D, typename F, std::size_t... INDICES>
// std::array<T,D> make_array(F&& f, std::index_sequence<INDICES...>)
// {
//   return std::array<T,D> {(static_cast<void>(INDICES), f())...};
// }

// /// construct a std::array<T,D> where T is an object of class
// template <typename T, unsigned D, typename F>
// std::array<T,D> make_array(F&& f)
// {
//   return make_array<T,D,F>(std::forward<F&&>(f), std::make_index_sequence<D> {});
// }


template <typename T = int>
struct Fraction {
  T numerator;
  T denominator;

  constexpr Fraction(T num, T denum) : numerator(num), denominator(denum) {}
  constexpr Fraction(T num) : Fraction(num, 1) {}
  constexpr Fraction() : Fraction(0) {}

  // return result of numerator/denominator
  template <typename U = T>
  constexpr U operator()() const {
    if (U(denominator) == U(0)) {
      std::cout << "denominator is zero !" << std::endl;
      exit(-1);
    }
    return U(numerator) / U(denominator);
  }
  
  void set(T num, T denum) {
    numerator = num;
    denominator = denum;
  }
  // seperately add to numerator and denominator
  void add(T num, T denum) {
    numerator += num;
    denominator += denum;
  }
};

struct OMP_helper {
  static void check_Thread_Num() {
    int max_threads = omp_get_max_threads();
    if (Thread_Num > max_threads) {
      std::cout << "Thread_Num is too large, max_threads = " << max_threads
                << ", but Thread_Num = " << Thread_Num << std::endl;
      std::cout << "press any key to continue..." << std::endl;
      std::cin.get();
    }
  }
  static int get_revised_Thread_Num() {
    int max_threads = omp_get_max_threads();
    if (Thread_Num > omp_get_max_threads())
      return max_threads;
    else
      return Thread_Num;
  }

  // performing reserve for all vectors in vec_th, and push_back for all
  template <typename T>
  static void Setup_Thread_Vec(std::vector<std::vector<T>> &vec_th, int TotalReserveSize,
                               int SubVecNum = Thread_Num) {
    int size_thread = TotalReserveSize / SubVecNum + 1;
    for (int i = 0; i < SubVecNum; i++) {
      std::vector<T> vec_thread_i;
      vec_thread_i.reserve(size_thread);
      vec_th.emplace_back(vec_thread_i);
    }
  }

  // create vector of vectors for performing parallel computation
  // return a reference to vec_th
  // need to delete vec_th after use
  template <typename T>
  static std::vector<std::vector<T>> &Create_Thread_Vec(int N) {
    std::vector<std::vector<T>> *vec_th = new std::vector<std::vector<T>>;
    Setup_Thread_Vec(*vec_th, N);
    return *vec_th;
  }

  // merge all vectors in vec_th to vec, do not clear vec_th
  template <typename T>
  static void Merge(std::vector<T> &vec, const std::vector<std::vector<T>> &vec_th) {
    for (int i = 0; i < vec_th.size(); i++) {
      vec.insert(vec.end(), vec_th[i].begin(), vec_th[i].end());
    }
  }

  // clear all sub-vectors in vec_th
  template <typename T>
  static void Clear(std::vector<std::vector<T>> &vec_th) {
    for (int i = 0; i < vec_th.size(); i++) {
      vec_th[i].clear();
    }
  }

  template <typename T>
  static void Clear_and_Merge(std::vector<T> &vec, std::vector<std::vector<T>> &vec_th) {
    vec.clear();
    for (int i = 0; i < vec_th.size(); i++) {
      vec.insert(vec.end(), vec_th[i].begin(), vec_th[i].end());
    }
  }

  // merge all vectors in vec_th to vec, and clear vec_th
  template <typename T>
  static void Merge_and_Clear(std::vector<T> &vec, std::vector<std::vector<T>> &vec_th) {
    for (int i = 0; i < vec_th.size(); i++) {
      vec.insert(vec.end(), vec_th[i].begin(), vec_th[i].end());
      vec_th[i].clear();
    }
  }

  // a simple load balancer:
  // get max and min size of sub-vectors
  // if max/min > x, then re-divide vector into sub-vectors
  template <typename T>
  static void LoadBalancer(std::vector<T> &vec, std::vector<std::vector<T>> &vec_th,
                           int load_ratio = 2, int Num = Thread_Num) {
    // get max and min size of sub-vectors
    int max = 0, min = vec.size();
    for (int i = 0; i < vec_th.size(); i++) {
      if (vec_th[i].size() > max) max = vec_th[i].size();
      if (vec_th[i].size() < min) min = vec_th[i].size();
    }
    // check if re-divide is needed
    if (max / min > load_ratio) Divide(vec, vec_th, Num);
  }

  // a simple load balancer:
  // get max and min size of sub-vectors
  // if max/min > x, then re-divide vector into sub-vectors
  // just need to pass sub-vectors
  template <typename T>
  static void LoadBalancer(std::vector<std::vector<T>> &vec_th, int load_ratio = 2,
                           int Num = Thread_Num) {
    // get max and min size of sub-vectors
    int max = 0, min = vec_th[0].size(), total_size = 0;
    for (int i = 0; i < vec_th.size(); i++) {
      int thread_size = vec_th[i].size();
      total_size += vec_th[i].size();
      max = thread_size > max ? thread_size : max;
      min = thread_size < min ? thread_size : min;
    }
    min = min == 0 ? 1 : min;
    // check if re-divide is needed
    if ((double(max) / double(min)) > load_ratio) {
      std::vector<T> vec;
      vec.reserve(total_size);
      Merge(vec, vec_th);
      Divide(vec, vec_th, Num);
    }
  }

  // set Index to access a container in parallel
  // start and end will be used like:
  // Func(start[i], end[i] ...){
  // for (int i = start; i < end; i++) {...}}
  static void Divide_Index(std::vector<int> &start, std::vector<int> &end, int num,
                           int Npart = Thread_Num) {
    start.clear();
    end.clear();
    int size_th = num / Npart;
    if (size_th == 0) {
      for (int i = 0; i < num; i++) {
        // 0 - num-1: 1
        // num - Npart-1: 0
        start.emplace_back(i);
        end.emplace_back(i + 1);
      }
      for (int i = num; i < Npart; i++) {
        start.emplace_back(0);
        end.emplace_back(0);
      }
    } else {
      int _start = 0;
      for (int i = 0; i < Npart; i++) {
        int _end = (i == Npart - 1) ? num : _start + size_th;
        start.emplace_back(_start);
        end.emplace_back(_end);
        _start = _end;
      }
    }
  }

  // divide a certain num to N(default: Thread_Num) parts
  // return a vector of N parts
  static std::vector<int> Divide_Num(int num, int Npart = Thread_Num) {
    std::vector<int> vec;
    Divide_Num(vec, num, Npart);
    return vec;
  }

  // divide a certain num to N(default: Thread_Num) parts
  // and store the result in vec
  static void Divide_Num(std::vector<int> &vec, int num, int Npart = Thread_Num) {
    vec.clear();
    vec.reserve(Npart);
    int size_th = num / Npart;
    if (size_th == 0) {
      // 0 - num-1: 1
      // num - Npart-1: 0
      for (int i = 0; i < num; i++) vec.emplace_back(1);
      for (int i = num; i < Npart; i++) vec.emplace_back(0);
    } else {
      int start = 0;
      for (int i = 0; i < Npart; i++) {
        int end = (i == Npart - 1) ? num : start + size_th;
        vec.emplace_back(end - start);
        start = end;
      }
    }
  }

  // Divide vector into N(default: Thread_Num) vectors
  // need to define std::vector<std::vector<T>> vec_th first
  // and setup vec_th using Setup_Thread_Vec<>(vec_th, N)
  template <typename T>
  static void Divide(const std::vector<T> &vec, std::vector<std::vector<T>> &vec_th,
                     int Npart = Thread_Num) {
    std::vector<int> vec_th_size = Divide_Num(vec.size(), Npart);
    auto it = vec.begin();
    for (int i = 0; i < Npart; i++) {
      if (vec_th_size[i] == 0) continue;
      auto _end = it + vec_th_size[i];
      vec_th[i].insert(vec_th[i].end(), it, _end);
      it = _end;
    }
    // Legacy
    //   int size_th = vec.size() / Npart;
    //   if (size_th == 0) {
    //     std::cout << "subvector size < 1" << std::endl;
    //     exit(-1);
    //   }
    //   auto it = vec.begin();
    //   for (int i = 0; i < Npart; i++) {
    //     auto _end = (i == Npart - 1) ? vec.end() : it + size_th;
    //     vec_th[i].insert(vec_th[i].end(), it, _end);
    //     it = _end;
    //   }
  }

  // Divide vector into N(default: Thread_Num) vectors
  // need to define std::vector<std::vector<T>> vec_th first
  // setup is integrated into this function
  template <typename T>
  static void Divide_and_Setup(const std::vector<T> &vec,
                               std::vector<std::vector<T>> &vec_th,
                               int Npart = Thread_Num) {
    Setup_Thread_Vec(vec_th, vec.size(), Npart);
    Divide(vec, vec_th, Npart);
  }
};


template <typename T>
class GenericArray;
template <typename T, unsigned int D>
class Vector;

namespace util {
// return true if a is close to zero
template <typename T>
inline bool nearZero(T a) {
  if (a == T()) {
    return true;
  }
  constexpr T EPSILON = std::numeric_limits<T>::epsilon();
  if (a > -EPSILON && a < EPSILON) {
    return true;
  } else {
    return false;
  }
}
// return true if a is close to epsilon
template <typename T>
inline bool nearZero(T a, T epsilon) {
  if (a > -epsilon && a < epsilon) {
    return true;
  } else {
    return false;
  }
}

// equality of two floating point numbers
template <typename T>
inline bool fpequal(T x, T y, T eps = std::numeric_limits<T>::epsilon()) {
  // This is an attempt to treat underflow.
  if (std::fabs(x) <= eps) x = T{};
  // This is an attempt to treat underflow.
  if (std::fabs(y) <= eps) y = T{};

  if (x == T{} || y == T{}) {
    return (std::fabs(x - y) <= eps);
  } else {
    return (std::fabs(x - y) <= eps * std::fabs(x) && std::fabs(x - y) <= eps * std::fabs(y));
  }
}

template<typename T, unsigned int D>
inline bool Vectorfpequal(const Vector<T, D>& x, const Vector<T, D>& y, T eps = std::numeric_limits<T>::epsilon()) {
	for (unsigned int i = 0; i < D; ++i) {
		if (!fpequal(x[i], y[i], eps)) return false;
	}
	return true;
}

// get underlying value of enum 
template <typename T>
constexpr auto underlyingVal(T e) noexcept {
  return static_cast<std::underlying_type_t<T>>(e);
}
// get underlying reference of enum
template <typename T>
constexpr auto& underlyingRef(T& e) noexcept {
  return reinterpret_cast<std::underlying_type_t<T>&>(e);
}
template <typename T>
constexpr const auto& underlyingRef(const T& e) noexcept {
  return reinterpret_cast<const std::underlying_type_t<T>&>(e);
}

// compare 2 uint8_t flags by bitwise AND operation
// return static_cast<bool>(uint8_t flag1 & uint8_t flag2);
__any__ inline bool isFlag(std::uint8_t flag1, std::uint8_t flag2) {
  return static_cast<bool>(flag1 & flag2);
}
// compare 2 uint8_t flags by bitwise AND operation
// return static_cast<bool>(uint16_t flag1 & uint16_t flag2);
// __any__ inline bool isFlag(std::uint16_t flag1, std::uint16_t flag2) {
//   return static_cast<bool>(flag1 & flag2);
// }
// add flagx to existing flag
inline void addFlag(std::uint8_t src, std::uint8_t &dst) { dst |= src; }
// remove flagx from existing flag
inline void removeFlag(std::uint8_t src, std::uint8_t &dst) { dst &= ~src; }

// get effective number of elements in a symmetric matrix
template <unsigned int N>
constexpr unsigned int SymmetricMatrixSize() {
  return N * (N + 1) / 2;
}

// get norm square of a tensor
template <typename T, unsigned int N>
T NormSquare(const std::array<T, SymmetricMatrixSize<N>()> &tensor) {
  T normsquare{};
  if constexpr (N == 2) {
    normsquare = tensor[0] * tensor[0] 
               + tensor[1] * tensor[1] * 2
               + tensor[2] * tensor[2];
  } else if constexpr (N == 3) {
    normsquare = tensor[0] * tensor[0] 
               + tensor[1] * tensor[1] * 2
               + tensor[2] * tensor[2] * 2
               + tensor[3] * tensor[3] 
               + tensor[4] * tensor[4] * 2
               + tensor[5] * tensor[5];
  }
  return normsquare;
}



template <typename ArrayType, unsigned int Dim>
void CopyFromFieldArrayImp(const Vector<int, Dim> &Mesh, int Overlap, const ArrayType &Array,
                        typename ArrayType::value_type *dst) {
	if (Overlap == 0) {
		const std::size_t Size = Array.size();
		std::copy(Array.getdataPtr(), Array.getdataPtr(Size), dst);
		return;
	}

  int delta_x = Mesh[0] - 2 * Overlap;
  std::size_t id = 0;
  if constexpr (Dim == 2) {
    for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
      std::size_t start = Overlap + j * Mesh[0];
      std::copy(Array.getdataPtr(start), Array.getdataPtr(start + delta_x), dst + id);
      id += delta_x;
    }
  } else if constexpr (Dim == 3) {
    std::size_t XY = Mesh[0] * Mesh[1];
    for (int k = Overlap; k < Mesh[2] - Overlap; ++k) {
      for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
        std::size_t start = Overlap + j * Mesh[0] + k * XY;
        std::copy(Array.getdataPtr(start), Array.getdataPtr(start + delta_x), dst + id);
        id += delta_x;
      }
    }
  }
}

template <typename ArrayType, unsigned int Dim, typename datatype>
void CastFromFieldArrayImp(const Vector<int, Dim> &Mesh, int Overlap, const ArrayType &Array, datatype *dst) {
  if (Overlap == 0) {
    const std::size_t Size = Array.size();
    for (std::size_t i = 0; i < Size; ++i) {
      dst[i] = static_cast<datatype>(Array[i]);
    }
    return;
  }

  int delta_x = Mesh[0] - 2 * Overlap;
  std::size_t id{};
  if constexpr (Dim == 2) {
    for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
      std::size_t start = Overlap + j * Mesh[0];
      for (int i = 0; i < delta_x; ++i) {
        dst[id] = static_cast<datatype>(Array[start + i]);
        ++id;
      }
    }
  } else if constexpr (Dim == 3) {
    std::size_t XY = Mesh[0] * Mesh[1];
    for (int k = Overlap; k < Mesh[2] - Overlap; ++k) {
      for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
        std::size_t start = Overlap + j * Mesh[0] + k * XY;
        for (int i = 0; i < delta_x; ++i) {
          dst[id] = static_cast<datatype>(Array[start + i]);
          ++id;
        }
      }
    }
  }
}

// copy data from field
template <typename ArrayType, unsigned int Dim, typename datatype>
void CopyFromFieldArray(const Vector<int, Dim> &Mesh, int Overlap, const ArrayType &Array, datatype *dst) {
  // check if datatype is the same as ArrayType::value_type and decide which scheme to use
  if constexpr (std::is_same_v<datatype, typename ArrayType::value_type>) {
    // just copy data
    CopyFromFieldArrayImp(Mesh, Overlap, Array, dst);
  } else {
    // copy and convert data
    CastFromFieldArrayImp(Mesh, Overlap, Array, dst);
  }
}


template <typename ArrayType, unsigned int Dim>
void CopyFromFieldArrayImp(const Vector<int, Dim> &Mesh, int Overlap, const ArrayType &Array, typename ArrayType::value_type *dst, 
	std::function<typename ArrayType::value_type(typename ArrayType::value_type)> f) {
  if (Overlap == 0) {
    const std::size_t Size = Array.size();
    for (std::size_t i = 0; i < Size; ++i) {
      dst[i] = f(Array[i]);
    }
    return;
  }
  
  int delta_x = Mesh[0] - 2 * Overlap;
  std::size_t id{};
  if constexpr (Dim == 2) {
    for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
      std::size_t start = Overlap + j * Mesh[0];
      for (int i = 0; i < delta_x; ++i) {
        dst[id] = f(Array[start + i]);
        ++id;
      }
    }
  } else if constexpr (Dim == 3) {
    std::size_t XY = Mesh[0] * Mesh[1];
    for (int k = Overlap; k < Mesh[2] - Overlap; ++k) {
      for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
        std::size_t start = Overlap + j * Mesh[0] + k * XY;
        for (int i = 0; i < delta_x; ++i) {
          dst[id] = f(Array[start + i]);
          ++id;
        }
      }
    }
  }
}

template <typename ArrayType, unsigned int Dim, typename datatype>
void CastFromFieldArrayImp(const Vector<int, Dim> &Mesh, int Overlap, const ArrayType &Array, datatype *dst, 
	std::function<typename ArrayType::value_type(typename ArrayType::value_type)> f) {
  if (Overlap == 0) {
    const std::size_t Size = Array.size();
    for (std::size_t i = 0; i < Size; ++i) {
      dst[i] = static_cast<datatype>(f(Array[i]));
    }
    return;
  }
  
  int delta_x = Mesh[0] - 2 * Overlap;
  std::size_t id{};
  if constexpr (Dim == 2) {
    for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
      std::size_t start = Overlap + j * Mesh[0];
      for (int i = 0; i < delta_x; ++i) {
        dst[id] = static_cast<datatype>(f(Array[start + i]));
        ++id;
      }
    }
  } else if constexpr (Dim == 3) {
    std::size_t XY = Mesh[0] * Mesh[1];
    for (int k = Overlap; k < Mesh[2] - Overlap; ++k) {
      for (int j = Overlap; j < Mesh[1] - Overlap; ++j) {
        std::size_t start = Overlap + j * Mesh[0] + k * XY;
        for (int i = 0; i < delta_x; ++i) {
          dst[id] = static_cast<datatype>(f(Array[start + i]));
          ++id;
        }
      }
    }
  }
}

// copy data from field with unit conversion
template <typename ArrayType, unsigned int Dim, typename datatype>
void CopyFromFieldArray(const Vector<int, Dim> &Mesh, int Overlap, const ArrayType &Array, datatype *dst, 
	std::function<typename ArrayType::value_type(typename ArrayType::value_type)> f) {
  // check if datatype is the same as ArrayType::value_type and decide which scheme to use
  if constexpr (std::is_same_v<datatype, typename ArrayType::value_type>) {
    // copy data after unit conversion
    CopyFromFieldArrayImp(Mesh, Overlap, Array, dst, f);
  } else {
    // copy and convert data after unit conversion
    CastFromFieldArrayImp(Mesh, Overlap, Array, dst, f);
  }
}



#include <string>

  // trim spaces from both ends of a string
std::string trim(const std::string& str) {
  std::size_t first = str.find_first_not_of(' ');
  if (std::string::npos == first) {
    return str;
  }
  std::size_t last = str.find_last_not_of(' ');
  return str.substr(first, (last - first + 1));
}

template<typename T>
__any__ inline T clamp(T val, T min, T max) {
  return val < min ? min : (val > max ? max : val);
}


}  // namespace util
