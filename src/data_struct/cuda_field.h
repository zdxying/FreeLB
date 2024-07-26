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

// field for cuda

#pragma once

#ifdef __CUDACC__

#include <cstddef>

#include "data_struct/Vector.h"
#include "utils/alias.h"
#include "utils/cuda_device.h"
#include "utils/util.h"

namespace cudev {

// single data
template <typename T, typename Base>
class Data {
 public:
  using value_type = T;
  static constexpr unsigned int array_dim = 1;
  using array_type = Data<T, Base>;
  static constexpr bool isField = false;
  static constexpr bool isCuDevField = true;
  using cudev_FieldType = Data<T, Base>;

 public:
  T* _data;

  __device__ Data() : _data(nullptr) {}
  // construct in host, then copy to device
  __any__ Data(T* data) : _data(data) {}
  // Copy constructor
  __device__ Data(const Data& data) : _data(data._data) {}
  // Move constructor
  __device__ Data(Data&& data) noexcept : _data(data._data) {}
  // Copy assignment operator
  __device__ Data& operator=(const Data& data) {
    if (&data == this) return *this;
    _data = data._data;
    return *this;
  }
  // Move assignment operator
  __device__ Data& operator=(Data&& data) noexcept {
    if (&data == this) return *this;
    _data = data._data;
    return *this;
  }

  template <unsigned int i = 0>
  __device__ auto& get() {
    return *_data;
  }
  template <unsigned int i = 0>
  __device__ const auto& get() const {
    return *_data;
  }
  __device__ auto& get(unsigned int i) { return *_data; }
  __device__ const auto& get(unsigned int i) const { return *_data; }

  // argument i will not be used
  __device__ auto& getField(std::size_t i = 0) { return *this; }
  __device__ const auto& getField(std::size_t i = 0) const { return *this; }

  template <unsigned int i = 0>
  __device__ void SetField(T value) {
    *_data = value;
  }
  __device__ void SetField(int i, T value) { *_data = value; }

  __device__ static constexpr unsigned int Size() { return 1; }
};


// a single std::array, type T couldn't dynamically allocate memory
template <typename T, typename Base>
class Array {
 public:
  using value_type = T;
  static constexpr unsigned int array_dim = Base::array_dim;
  using array_type = Array<T, Base>;
  static constexpr bool isField = false;
  static constexpr bool isCuDevField = true;
  using cudev_FieldType = Array<T, Base>;

 private:
  T* _data;

 public:
  __device__ Array() : _data(nullptr) {}
  __any__ Array(T* data) : _data(data) {}
  // Copy constructor
  __device__ Array(const Array& arr) : _data(arr._data) {}
  // Move constructor
  __device__ Array(Array&& arr) noexcept : _data(arr._data) {}
  // Copy assignment operator
  __device__ Array& operator=(const Array& arr) {
    if (&arr == this) return *this;
    _data = arr._data;
    return *this;
  }
  // Move assignment operator
  __device__ Array& operator=(Array&& arr) noexcept {
    if (&arr == this) return *this;
    _data = arr._data;
    return *this;
  }

  template <unsigned int i = 0>
  __device__ auto& get() {
    return _data[i];
  }
  template <unsigned int i = 0>
  __device__ const auto& get() const {
    return _data[i];
  }
  __device__ auto& get(unsigned int i) { return _data[i]; }
  __device__ const auto& get(unsigned int i) const { return _data[i]; }

  __device__ T* getArray() { return _data; }
  __device__ const T* getArray() const { return _data; }

  // argument i will not be used
  __device__ auto& getField(std::size_t i = 0) { return *this; }
  __device__ const auto& getField(std::size_t i = 0) const { return *this; }

  template <unsigned int i = 0>
  __device__ void SetField(T value) {
    _data[i] = value;
  }
  __device__ void SetField(int i, T value) { _data[i] = value; }

  __device__ static constexpr unsigned int Size() { return array_dim; }
};


// a class containing std::vector with some methods
template <typename T>
class Genericvector {
 private:
  T* _data;
  std::size_t* _count;
  std::size_t _capacity;

 public:
  using value_type = T;

  __device__ Genericvector() : _data(nullptr), _count(nullptr), _capacity(0) {}
  __any__ Genericvector(T* data, std::size_t* size)
      : _data(data), _count(size), _capacity(0) {
#ifdef __CUDA_ARCH__
    _capacity = *size;
#endif
  }
  // Copy constructor
  __device__ Genericvector(const Genericvector& arr)
      : _data(arr._data), _count(arr._count), _capacity(arr._capacity) {}
  // Move constructor
  __device__ Genericvector(Genericvector&& arr) noexcept
      : _data(arr._data), _count(arr._count), _capacity(arr._capacity) {}
  // Copy assignment operator
  __device__ Genericvector& operator=(const Genericvector& arr) {
    if (&arr == this) return *this;
    _data = arr._data;
    _count = arr._count;
    _capacity = arr._capacity;
    return *this;
  }
  // Move assignment operator
  __device__ Genericvector& operator=(Genericvector&& arr) noexcept {
    if (&arr == this) return *this;
    _data = arr._data;
    _count = arr._count;
    _capacity = arr._capacity;
    return *this;
  }

  __device__ const T& operator[](std::size_t i) const { return _data[i]; }
  __device__ T& operator[](std::size_t i) { return _data[i]; }

  // push back without checking the size
  __device__ void push_back(T value) {
    if (*_count < _capacity) {
      _data[*_count] = value;
      (*_count)++;
    }
  }
  __device__ void clear() { *_count = 0; }

  // return the pointer of ith element
  __device__ T* getdataPtr(std::size_t i = 0) { return _data + i; }
  __device__ const T* getdataPtr(std::size_t i = 0) const { return _data + i; }

  __device__ inline void set(std::size_t i, T value) { _data[i] = value; }

  __device__ std::size_t size() const { return *_count; }
  __device__ std::size_t capacity() const { return _capacity; }
};


template <typename T>
class GenericArray {
 private:
  // number of elements
  std::size_t* _count;
  // pointer to the data
  T* _data;

 public:
  using value_type = T;

  __device__ GenericArray() : _count(nullptr), _data(nullptr) {}
  __any__ GenericArray(std::size_t* size, T* data) : _count(size), _data(data) {}
  // Copy constructor
  __device__ GenericArray(const GenericArray& arr)
      : _count(arr._count), _data(arr._data) {}
  // Move constructor
  __device__ GenericArray(GenericArray&& arr) noexcept
      : _count(arr._count), _data(arr._data) {}
  // Copy assignment operator
  __device__ GenericArray& operator=(const GenericArray& arr) {
    if (&arr == this) return *this;
    _count = arr._count;
    _data = arr._data;
    return *this;
  }
  // Move assignment operator
  __device__ GenericArray& operator=(GenericArray&& arr) noexcept {
    if (&arr == this) return *this;
    _count = arr._count;
    _data = arr._data;
    return *this;
  }

  __device__ void Init(T InitValue) {
    for (std::size_t i = 0; i < *_count; ++i) {
      _data[i] = InitValue;
    }
  }

  __device__ const T& operator[](std::size_t i) const { return _data[i]; }
  __device__ T& operator[](std::size_t i) { return _data[i]; }

  // return pointer to ith element
  __device__ T* getdataPtr(std::size_t i = 0) { return _data + i; }
  __device__ const T* getdataPtr(std::size_t i = 0) const { return _data + i; }

  __device__ inline void set(std::size_t i, T value) { _data[i] = value; }

  __device__ std::size_t size() const { return *_count; }
};


// a modified version of CyclicArray
// rotate is handled at host side

template <typename T>
class StreamArray {
 private:
  // number of elements
  std::size_t* count;
  // base pointer to the data
  T* _data;
  // shift
  std::ptrdiff_t* _shift;
  T* _start;
  // facilitate the access of data before the last shift(rotate)
  std::ptrdiff_t* Offset;

 public:
  using value_type = T;

  __device__ StreamArray()
      : count(0), _data(nullptr), _shift(0), _start(nullptr), Offset(0) {}
  __any__ StreamArray(std::size_t* size, T* data, std::ptrdiff_t* shift, T* start,
                      std::ptrdiff_t* offset)
      : count(size), _data(data), _shift(shift), _start(start), Offset(offset) {
    _start = _data;
  }
  // Copy constructor
  __device__ StreamArray(const StreamArray& arr)
      : count(arr.count), _data(arr._data), _shift(arr._shift), _start(arr._start),
        Offset(arr.Offset) {}
  // Move constructor
  __device__ StreamArray(StreamArray&& arr) noexcept
      : count(arr.count), _data(arr._data), _shift(arr._shift), _start(arr._start),
        Offset(arr.Offset) {}
  // Copy assignment operator
  __device__ StreamArray& operator=(const StreamArray& arr) {
    if (&arr == this) return *this;
    count = arr.count;
    _data = arr._data;
    _shift = arr._shift;
    _start = arr._start;
    Offset = arr.Offset;
    return *this;
  }
  // Move assignment operator
  __device__ StreamArray& operator=(StreamArray&& arr) noexcept {
    if (&arr == this) return *this;
    count = arr.count;
    _data = arr._data;
    _shift = arr._shift;
    _start = arr._start;
    Offset = arr.Offset;
    return *this;
  }

  __device__ void setOffset(int offset) { *Offset = offset; }

  __device__ const T& operator[](std::size_t i) const { return _start[i]; }
  __device__ T& operator[](std::size_t i) { return _start[i]; }

  __device__ inline void set(std::size_t i, T value) { _start[i] = value; }
  __device__ std::size_t size() const { return *count; }
  // return the pointer of ith element
  __device__ T* getdataPtr(std::size_t i = 0) { return _start + i; }
  __device__ const T* getdataPtr(std::size_t i = 0) const { return _start + i; }

  // get data before the last shift(rotate), used in bcs
  __device__ T& getPrevious(std::size_t i) {
    std::ptrdiff_t prevIndex = i + *Offset;
    if (prevIndex < 0) {
      prevIndex += *count;
    } else if (prevIndex >= static_cast<std::ptrdiff_t>(*count)) {
      prevIndex -= *count;
    }
    return _start[static_cast<std::size_t>(prevIndex)];
  }

  // experimental
  __device__ void rotate() {
    const std::ptrdiff_t n = *count;
    std::ptrdiff_t shift = *_shift;
    shift -= *Offset;
    if (shift >= n) {
      shift -= n;
      copyToFront(shift);
    } else if (shift < 0) {
      shift += n;
      copyToBack(shift);
    }
    *_shift = shift;
    _start = _data + shift;
  }
  __device__ void rotate(std::ptrdiff_t offset) {
    const std::ptrdiff_t n = *count;
    std::ptrdiff_t shift = *_shift;
    shift -= offset;
    if (shift >= n) {
      shift -= n;
      copyToFront(shift);
    } else if (shift < 0) {
      shift += n;
      copyToBack(shift);
    }
    *_shift = shift;
    _start = _data + shift;
  }
  __device__ void copyToBack(std::ptrdiff_t endoffset = 0) {
    T* const base = _data;
    endoffset = endoffset == 0 ? *count : endoffset;
    dev_copy(base + *count, base, endoffset);
  }
  __device__ void copyToFront(std::ptrdiff_t startoffset = 0) {
    T* const base = _data;
    dev_copy(base + startoffset, base + *count + startoffset, *count - startoffset);
  }
};

template <typename T>
__global__ void Stream_kernel(cudev::StreamArray<T>* arr) {
  arr->rotate();
}

template <typename T>
class StreamMapArray {
 private:
  // number of elements
  std::size_t* count;
  // base pointer to the data
  T* _data;
  // shift
  std::ptrdiff_t* _shift;
  T* _start;
  // facilitate the access of data before the last shift(rotate)
  std::ptrdiff_t* Offset;

 public:
  using value_type = T;

  __device__ StreamMapArray()
      : count(0), _data(nullptr), _shift(0), _start(nullptr), Offset(0) {}
  __any__ StreamMapArray(std::size_t* size, T* data, std::ptrdiff_t* shift, T* start,
                      std::ptrdiff_t* offset)
      : count(size), _data(data), _shift(shift), _start(start), Offset(offset) {
    _start = _data;
  }
  // Copy constructor
  __device__ StreamMapArray(const StreamMapArray& arr)
      : count(arr.count), _data(arr._data), _shift(arr._shift), _start(arr._start),
        Offset(arr.Offset) {}
  // Move constructor
  __device__ StreamMapArray(StreamMapArray&& arr) noexcept
      : count(arr.count), _data(arr._data), _shift(arr._shift), _start(arr._start),
        Offset(arr.Offset) {}
  // Copy assignment operator
  __device__ StreamMapArray& operator=(const StreamMapArray& arr) {
    if (&arr == this) return *this;
    count = arr.count;
    _data = arr._data;
    _shift = arr._shift;
    _start = arr._start;
    Offset = arr.Offset;
    return *this;
  }
  // Move assignment operator
  __device__ StreamMapArray& operator=(StreamMapArray&& arr) noexcept {
    if (&arr == this) return *this;
    count = arr.count;
    _data = arr._data;
    _shift = arr._shift;
    _start = arr._start;
    Offset = arr.Offset;
    return *this;
  }

  __device__ void setOffset(int offset) { *Offset = offset; }

  __device__ const T& operator[](std::size_t i) const { return _start[i]; }
  __device__ T& operator[](std::size_t i) { return _start[i]; }

  __device__ inline void set(std::size_t i, T value) { _start[i] = value; }
  __device__ std::size_t size() const { return *count; }
  // return the pointer of ith element
  __device__ T* getdataPtr(std::size_t i = 0) { return _start + i; }
  __device__ const T* getdataPtr(std::size_t i = 0) const { return _start + i; }

  // get data before the last shift(rotate), used in bcs
  __device__ T& getPrevious(std::size_t i) {
    std::ptrdiff_t prevIndex = i + *Offset;
    if (prevIndex < 0) {
      prevIndex += *count;
    } else if (prevIndex >= static_cast<std::ptrdiff_t>(*count)) {
      prevIndex -= *count;
    }
    return _start[static_cast<std::size_t>(prevIndex)];
  }

  // experimental
  __device__ void rotate() {
    const std::ptrdiff_t n = *count;
    std::ptrdiff_t shift = *_shift;
    shift -= *Offset;
    if (shift >= n) {
      shift -= n;
    } else if (shift < 0) {
      shift += n;
    }
    *_shift = shift;
    _start = _data + shift;
  }
};

template <typename T>
__global__ void Stream_kernel(cudev::StreamMapArray<T>* arr) {
  arr->rotate();
}

template <typename ArrayType, unsigned int D>
class GenericArrayField {
 public:
  using array_type = ArrayType;
  using value_type = typename ArrayType::value_type;
  static constexpr unsigned int array_dim = D;
  static constexpr bool isField = true;
  static constexpr bool isCuDevField = true;
  using cudev_FieldType = GenericArrayField<ArrayType, D>;

 private:
  ArrayType** _data;

 public:
  __device__ GenericArrayField() : _data{} {}
  __any__ GenericArrayField(ArrayType** data) : _data(data) {}
  // Copy constructor
  __device__ GenericArrayField(const GenericArrayField& genF) : _data(genF._data) {}
  // Move constructor
  __device__ GenericArrayField(GenericArrayField&& genF) noexcept : _data(genF._data) {}
  // Copy assignment operator
  __device__ GenericArrayField& operator=(const GenericArrayField& genF) {
    if (&genF == this) return *this;
    _data = genF._data;
    return *this;
  }
  // Move assignment operator
  __device__ GenericArrayField& operator=(GenericArrayField&& genF) noexcept {
    if (&genF == this) return *this;
    _data = genF._data;
    return *this;
  }

  __device__ ArrayType& getField(std::size_t i = 0) { return *_data[i]; }
  __device__ const ArrayType& getField(std::size_t i = 0) const { return *_data[i]; }

  // get<i>(id): return _data[i][id];
  template <unsigned int i = 0>
  __device__ auto& get(std::size_t id) {
    return (*_data[i])[id];
  }
  template <unsigned int i = 0>
  __device__ const auto& get(std::size_t id) const {
    return (*_data[i])[id];
  }
  __device__ auto& get(std::size_t id, unsigned int dir) { return (*_data[dir])[id]; }
  __device__ const auto& get(std::size_t id, unsigned int dir) const {
    return (*_data[dir])[id];
  }

  // get pointer to ith data in all arrays
  __device__ void getArray(std::size_t id, value_type** ptr_arr) {
    for (unsigned int i = 0; i < D; ++i) ptr_arr[i] = _data[i]->getdataPtr(id);
  }

  template <unsigned int i = 0>
  __device__ void SetField(std::size_t id, value_type value) {
    _data[i]->set(id, value);
  }
  __device__ void SetField(int i, std::size_t id, value_type value) {
    _data[i]->set(id, value);
  }

  __device__ static constexpr unsigned int Size() { return D; }
};

template <typename ArrayType, typename Base>
class GenericField : public GenericArrayField<ArrayType, Base::array_dim> {
 public:
  static constexpr unsigned int array_dim = Base::array_dim;
  using array_type = ArrayType;
  using value_type = typename ArrayType::value_type;
  using cudev_FieldType = GenericField<ArrayType, Base>;

  __any__ GenericField(ArrayType** data)
      : GenericArrayField<ArrayType, array_dim>(data) {}
};

}  // namespace cudev


#endif