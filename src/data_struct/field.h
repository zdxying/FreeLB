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

// physical field

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "data_struct/Vector.h"
// getIntpWeight / getInterpolation / getAverage (shared with parallel/communicator.h)
#include "data_struct/interpolation.h"
#include "utils/alias.h"
#include "utils/util.h"

#ifdef __CUDACC__
// #include "utils/cuda_device.h"
#include "data_struct/cuda_field.h"
#endif


// single data, T could be Vector but should not dynamically allocate memory
template <typename T, typename Base>
class Data {
 public:
  using value_type = T;
  static constexpr unsigned int array_dim = 1;
  using array_type = Data<T, Base>;
  static constexpr bool isField = false;
  static constexpr bool isCuDevField = false;

#ifdef __CUDACC__
  using cudev_array_type = cudev::Data<T, Base>;
#endif

 private:
  T _data;
#ifdef __CUDACC__
  T* dev_data;
  cudev::Data<T, Base>* dev_Data;
#endif

 public:
  Data() : _data{} { InitDeviceData(); }
  // argument size will not be used
  Data(std::size_t size) : Data() {}
  Data(std::size_t size, T initialValue) : _data(initialValue) { InitDeviceData(); }
  // Copy constructor
  Data(const Data& data) : _data(data._data) {
#ifdef __CUDACC__
    dev_data = cuda_malloc<T>(1);
    device_to_device(dev_data, data.dev_data, 1);
    constructInDevice();
#endif
  }
  // Move constructor
  Data(Data&& data) noexcept : _data(std::move(data._data)) {
    data._data = T{};
#ifdef __CUDACC__
    dev_data = data.dev_data;
    data.dev_data = nullptr;
    constructInDevice();
#endif
  }
  // Copy assignment operator
  Data& operator=(const Data& data) {
    if (&data == this) return *this;
    _data = data._data;
#ifdef __CUDACC__
    device_to_device(dev_data, data.dev_data, 1);
#endif
    return *this;
  }
  // Move assignment operator
  Data& operator=(Data&& data) noexcept {
    if (&data == this) return *this;
    _data = std::move(data._data);
    data._data = T{};
#ifdef __CUDACC__
    dev_data = data.dev_data;
    data.dev_data = nullptr;
#endif
    return *this;
  }

  ~Data() {
#ifdef __CUDACC__
    if (dev_data) cuda_free(dev_data);
    if (dev_Data) cuda_free(dev_Data);
#endif
  }

  void InitDeviceData() {
#ifdef __CUDACC__
    dev_data = cuda_malloc<T>(1);
    copyToDevice();
    constructInDevice();
#endif
  }

#ifdef __CUDACC__
  void copyToDevice() { host_to_device(dev_data, &_data, 1); }
  void copyToHost() { device_to_host(&_data, dev_data, 1); }
  T* get_devptr() { return dev_data; }
  cudev::Data<T, Base>* get_devObj() { return dev_Data; }
  void constructInDevice() {
    dev_Data = cuda_malloc<cudev::Data<T, Base>>(1);
    // temp host object
    cudev::Data<T, Base> temp(dev_data);
    // copy to device
    host_to_device(dev_Data, &temp, 1);
  }
#endif

  template <unsigned int i = 0>
  auto& get() {
    return _data;
  }
  template <unsigned int i = 0>
  const auto& get() const {
    return _data;
  }
  auto& get(unsigned int i) { return _data; }
  const auto& get(unsigned int i) const { return _data; }

  // argument i will not be used
  auto& getField(std::size_t i = 0) { return *this; }
  const auto& getField(std::size_t i = 0) const { return *this; }

  // init
  void Init(T value = T{}) {
    _data = value;
#ifdef __CUDACC__
    copyToDevice();
#endif
  }
  void Reset(T value = T{}) { Init(value); }
  template <unsigned int i = 0>
  void SetField(T value) {
    _data = value;
  }
  void SetField(int i, T value) { _data = value; }

  static constexpr unsigned int Size() { return 1; }
};


// a class containing std::vector with some methods
template <typename T>
class Genericvector {
 private:
  std::vector<T> data;
#ifdef __CUDACC__
  // remember to malloc enought memory on device,
  // dynamic memory allocation on device is not suggested
  T* dev_data;
  std::size_t* dev_count;
  cudev::Genericvector<T>* dev_Genericvector;
#endif

 public:
  using value_type = T;
  using array_type = Genericvector<T>;
#ifdef __CUDACC__
  using cudev_array_type = cudev::Genericvector<T>;
#endif

  Genericvector() {
#ifdef __CUDACC__
    dev_data = nullptr;
    dev_count = nullptr;
    dev_Genericvector = nullptr;
#endif
  }
  Genericvector(std::size_t size) {
    data.resize(size);
#ifdef __CUDACC__
    InitDeviceData(size);
#endif
  }
  Genericvector(std::size_t size, T InitValue) {
    data.resize(size, InitValue);
#ifdef __CUDACC__
    InitDeviceData(size);
#endif
  }
  // Copy constructor
  Genericvector(const Genericvector& arr) : data(arr.data) {
#ifdef __CUDACC__
    dev_data = cuda_malloc<T>(data.size());
    dev_count = cuda_malloc<std::size_t>(1);
    std::size_t dev_count_temp = arr.get_devcount();
    device_to_device(dev_data, arr.dev_data, dev_count_temp);
    device_to_device(dev_count, &dev_count_temp, 1);
    constructInDevice();
#endif
  }
  // Move constructor
  Genericvector(Genericvector&& arr) noexcept : data(std::move(arr.data)) {
#ifdef __CUDACC__
    dev_data = arr.dev_data;
    dev_count = arr.dev_count;
    arr.dev_data = nullptr;
    arr.dev_count = nullptr;
    constructInDevice();
#endif
  }
  // Copy assignment operator
  Genericvector& operator=(const Genericvector& arr) {
    if (&arr == this) return *this;
#ifdef __CUDACC__
    if (data.size() != arr.get_devcount()) {
      cuda_free(dev_data);
      dev_data = cuda_malloc<T>(arr.get_devcount());
    }
    std::size_t dev_count_temp = arr.get_devcount();
    device_to_device(dev_data, arr.dev_data, dev_count_temp);
    device_to_device(dev_count, &dev_count_temp, 1);
#endif
    data = arr.data;
    return *this;
  }
  // Move assignment operator
  Genericvector& operator=(Genericvector&& arr) noexcept {
    if (&arr == this) return *this;
    data = std::move(arr.data);
#ifdef __CUDACC__
    dev_data = arr.dev_data;
    dev_count = arr.dev_count;
    arr.dev_data = nullptr;
    arr.dev_count = nullptr;
#endif
    return *this;
  }

  ~Genericvector() {
#ifdef __CUDACC__
    if (dev_data) cuda_free(dev_data);
    if (dev_count) cuda_free(dev_count);
    if (dev_Genericvector) cuda_free(dev_Genericvector);
#endif
  }

  void InitDeviceData(std::size_t size) {
#ifdef __CUDACC__
    dev_data = cuda_malloc<T>(size);
    dev_count = cuda_malloc<std::size_t>(1);
    copyToDevice();
    constructInDevice();
#endif
  }

#ifdef __CUDACC__
  void copyToDevice() {
    host_to_device(dev_data, data.data(), data.size());
    std::size_t dev_count_temp = data.size();
    host_to_device(dev_count, &dev_count_temp, 1);
  }
  void copyToHost() {
    std::size_t dev_count_temp = get_devcount();
    data.resize(dev_count_temp);
    device_to_host(data.data(), dev_data, dev_count_temp);
  }
  T* get_devptr() { return dev_data; }
  std::size_t get_devcount() const {
    std::size_t temp;
    device_to_host(&temp, dev_count, 1);
    return temp;
  }
  cudev::Genericvector<T>* get_devObj() { return dev_Genericvector; }
  void constructInDevice() {
    dev_Genericvector = cuda_malloc<cudev::Genericvector<T>>(1);
    // temp host object
    cudev::Genericvector<T> temp(dev_data, dev_count);
    // copy to device
    host_to_device(dev_Genericvector, &temp, 1);
  }
#endif

  void Init(T InitValue) {
    std::fill(data.begin(), data.end(), InitValue);
#ifdef __CUDACC__
    copyToDevice();
#endif
  }

  void Resize(std::size_t size) {
    if (size == data.size()) return;
    data.resize(size);
#ifdef __CUDACC__
    cuda_free(dev_data);
    dev_data = cuda_malloc<T>(size);
    copyToDevice();
#endif
  }

  const T& operator[](std::size_t i) const { return data[i]; }
  T& operator[](std::size_t i) { return data[i]; }

  // return the pointer of ith element
  T* getdataPtr(std::size_t i = 0) { return data.data() + i; }
  const T* getdataPtr(std::size_t i = 0) const { return data.data() + i; }
  // return reference to the data
  std::vector<T>& getvector() { return data; }
  const std::vector<T>& getvector() const { return data; }

  inline void set(std::size_t i, T value) { data[i] = value; }

  std::size_t size() const { return data.size(); }
};


template <typename T>
class GenericArray {
 private:
  // number of elements
  std::size_t count;
  // base pointer to the data
  T* data;
#ifdef __CUDACC__
  std::size_t* dev_count;
  // device pointer to the data
  T* dev_data;
  cudev::GenericArray<T>* dev_GenericArray;
#endif

 public:
  using value_type = T;
  using array_type = GenericArray<T>;
#ifdef __CUDACC__
  using cudev_array_type = cudev::GenericArray<T>;
#endif

  GenericArray() : count(0), data(nullptr) {
#ifdef __CUDACC__
    dev_data = nullptr;
    dev_count = nullptr;
    dev_GenericArray = nullptr;
#endif
  }
  GenericArray(std::size_t size) : count(size), data(new T[size]{}) {
    InitDeviceData();
  }
  GenericArray(std::size_t size, T InitValue) : count(size), data(new T[size]) {
    std::fill(data, data + size, InitValue);
    InitDeviceData();
  }
  // Copy constructor
  GenericArray(const GenericArray& arr) : count(arr.count), data(new T[arr.count]) {
    std::copy(arr.data, arr.data + count, data);
#ifdef __CUDACC__
    dev_data = cuda_malloc<T>(count);
    dev_count = cuda_malloc<std::size_t>(1);
    device_to_device(dev_data, arr.dev_data, arr.count);
    device_to_device(dev_count, arr.dev_count, 1);
    constructInDevice();
#endif
  }
  // Move constructor
  GenericArray(GenericArray&& arr) noexcept : count(arr.count), data(arr.data) {
    // Reset 'arr'
    arr.count = 0;
    arr.data = nullptr;
#ifdef __CUDACC__
    dev_data = arr.dev_data;
    dev_count = arr.dev_count;
    arr.dev_data = nullptr;
    arr.dev_count = nullptr;
    constructInDevice();
#endif
  }
  // Copy assignment operator
  GenericArray& operator=(const GenericArray& arr) {
    if (&arr == this) return *this;
    const bool reallocated = (count != arr.count);
    if (reallocated) {
      delete[] data;
      data = new T[arr.count];
      count = arr.count;
#ifdef __CUDACC__
      if (dev_data) cuda_free(dev_data);
      if (dev_count) cuda_free(dev_count);
      if (dev_GenericArray) cuda_free(dev_GenericArray);
      dev_data = (count > 0) ? cuda_malloc<T>(count) : nullptr;
      dev_count = cuda_malloc<std::size_t>(1);
      dev_GenericArray = nullptr;
#endif
    }
    std::copy(arr.data, arr.data + arr.count, data);
#ifdef __CUDACC__
    if (count > 0) {
      device_to_device(dev_data, arr.dev_data, count);
      std::size_t n = count;
      host_to_device(dev_count, &n, 1);
    }
    if (reallocated || dev_GenericArray == nullptr) constructInDevice();
#endif
    return *this;
  }
  // Move assignment operator
  GenericArray& operator=(GenericArray&& arr) noexcept {
    if (&arr == this) return *this;
    delete[] data;
    // Steal the data from 'arr'
    count = arr.count;
    data = arr.data;
    // Reset 'arr'
    arr.count = 0;
    arr.data = nullptr;
#ifdef __CUDACC__
    if (dev_data) cuda_free(dev_data);
    if (dev_count) cuda_free(dev_count);
    if (dev_GenericArray) cuda_free(dev_GenericArray);
    dev_data = arr.dev_data;
    dev_count = arr.dev_count;
    dev_GenericArray = arr.dev_GenericArray;
    arr.dev_data = nullptr;
    arr.dev_count = nullptr;
    arr.dev_GenericArray = nullptr;
#endif
    return *this;
  }

  ~GenericArray() {
    delete[] data;
#ifdef __CUDACC__
    if (dev_data) cuda_free(dev_data);
    if (dev_count) cuda_free(dev_count);
    if (dev_GenericArray) cuda_free(dev_GenericArray);
#endif
  }

  void InitDeviceData() {
#ifdef __CUDACC__
    dev_data = cuda_malloc<T>(count);
    dev_count = cuda_malloc<std::size_t>(1);
    copyToDevice();
    constructInDevice();
#endif
  }

#ifdef __CUDACC__

  void copyToDevice() {
    host_to_device(dev_data, data, count);
    host_to_device(dev_count, &count, 1);
  }
  void copyToHost() {
    device_to_host(data, dev_data, count);
    device_to_host(&count, dev_count, 1);
  }
  T* get_devptr() { return dev_data; }
  std::size_t get_devcount() const {
    std::size_t temp;
    device_to_host(&temp, dev_count, 1);
    return temp;
  }
  cudev::GenericArray<T>* get_devObj() { return dev_GenericArray; }
  void constructInDevice() {
    dev_GenericArray = cuda_malloc<cudev::GenericArray<T>>(1);
    // temp host object
    cudev::GenericArray<T> temp(dev_count, dev_data);
    // copy to device
    host_to_device(dev_GenericArray, &temp, 1);
  }

#endif

  void Init(T InitValue) {
    std::fill(data, data + count, InitValue);
#ifdef __CUDACC__
    copyToDevice();
#endif
  }

  // Re-initialise in place; this is Init() for arrays without rotation state.
  void Reset(T InitValue) { Init(InitValue); }

  void Resize(std::size_t size) {
    if (size == count) return;
    delete[] data;
    data = new T[size]{};
    count = size;
#ifdef __CUDACC__
    if (dev_data) cuda_free(dev_data);
    if (dev_count) cuda_free(dev_count);
    if (dev_GenericArray) cuda_free(dev_GenericArray);
    InitDeviceData();
#endif
  }

  const T& operator[](std::size_t i) const { return data[i]; }
  T& operator[](std::size_t i) { return data[i]; }

  // return pointer to ith element
  T* getdataPtr(std::size_t i = 0) { return data + i; }
  const T* getdataPtr(std::size_t i = 0) const { return data + i; }

  inline void set(std::size_t i, T value) { data[i] = value; }

  std::size_t size() const { return count; }
};


// Kummerländer A, Dorn M, Frank M, Krause MJ. Implicit propagation of directly
// addressed grids in lattice Boltzmann methods. Concurrency Computat Pract
// Exper. 2023;35(8):e7509. doi: 10.1002/cpe.7509 a cyclic array to present periodic
// shift(PS) pattern

template <typename T>
class CyclicArray {
 private:
  // number of elements
  std::size_t count;
  // base pointer to the data
  T* data;
  // shift
  std::ptrdiff_t shift;
  std::size_t remainder;
  std::array<T*, 2> start;
  // facilitate the access of data before the last shift(rotate)
  std::ptrdiff_t Offset;
#ifdef __CUDACC__
  // device mirror handles; see the device block at the end of this class
  cudev::CyclicArray<T>* dev_CyclicArray;
  std::size_t* dev_count;
  T* dev_data;
  std::ptrdiff_t* dev_shift;
  std::size_t* dev_remainder;
  std::ptrdiff_t* dev_Offset;
#endif

 public:
  using value_type = T;
  using array_type = CyclicArray<T>;
#ifdef __CUDACC__
  using cudev_array_type = cudev::CyclicArray<T>;
#endif

  CyclicArray()
      : count(0), data(nullptr), shift(0), remainder(0), start{}, Offset(0) {
#ifdef __CUDACC__
    nullDeviceData();
    InitDeviceData();
#endif
  }
  CyclicArray(std::size_t size)
      : count(size), data(new T[size]{}), shift(0), remainder(size), Offset(0) {
    refresh();
#ifdef __CUDACC__
    nullDeviceData();
    InitDeviceData();
#endif
  }
  CyclicArray(std::size_t size, T InitValue)
      : count(size), data(new T[size]), shift(0), remainder(size), Offset(0) {
    std::fill(data, data + size, InitValue);
    refresh();
#ifdef __CUDACC__
    nullDeviceData();
    InitDeviceData();
#endif
  }
  // Copy constructor
  CyclicArray(const CyclicArray& arr)
      : count(arr.count), data(new T[arr.count]), shift(arr.shift),
        remainder(arr.remainder), Offset(arr.Offset) {
    std::copy(arr.data, arr.data + count, data);
    refresh();
#ifdef __CUDACC__
    nullDeviceData();
    InitDeviceData();
#endif
  }
  // Move constructor
  CyclicArray(CyclicArray&& arr) noexcept
      : count(0), data(nullptr), shift(0), remainder(0), start{}, Offset(0) {
    // Steal the data from 'arr'
    count = arr.count;
    data = arr.data;
    shift = arr.shift;
    remainder = arr.remainder;
    start = arr.start;
    Offset = arr.Offset;
    refresh();
    // Reset 'arr'
    arr.count = 0;
    arr.data = nullptr;
    arr.shift = 0;
    arr.remainder = 0;
    arr.start = {};
    arr.Offset = 0;
    refresh();
#ifdef __CUDACC__
    nullDeviceData();
    InitDeviceData();
#endif
  }
  // Copy assignment operator
  CyclicArray& operator=(const CyclicArray& arr) {
    if (&arr == this) return *this;
#ifdef __CUDACC__
    freeDeviceData();
    nullDeviceData();
#endif
    delete[] data;
    count = arr.count;
    data = new T[arr.count];
    std::copy(arr.data, arr.data + count, data);
    shift = arr.shift;
    remainder = arr.remainder;
    Offset = arr.Offset;
    refresh();
#ifdef __CUDACC__
    InitDeviceData();
#endif
    return *this;
  }
  // Move assignment operator
  CyclicArray& operator=(CyclicArray&& arr) noexcept {
    if (&arr == this) return *this;
#ifdef __CUDACC__
    freeDeviceData();
    nullDeviceData();
#endif
    delete[] data;
    // Steal the data from 'arr'
    count = arr.count;
    data = arr.data;
    shift = arr.shift;
    remainder = arr.remainder;
    start = arr.start;
    Offset = arr.Offset;
    // Reset 'arr'
    arr.count = 0;
    arr.data = nullptr;
    arr.shift = 0;
    arr.remainder = 0;
    arr.start = {};
    arr.Offset = 0;
    refresh();
#ifdef __CUDACC__
    InitDeviceData();
#endif
    return *this;
  }

  void Init(T InitValue, int offset = 0) { 
    std::fill(data, data + count, InitValue); 
    Offset = offset;
  }
  // Re-initialise in place and drop the rotation state, so the result matches
  // a freshly constructed array (Init() alone intentionally keeps shift).
  void Reset(T InitValue) {
    std::fill(data, data + count, InitValue);
    shift = 0;
    Offset = 0;
    refresh();
#ifdef __CUDACC__
    copyToDevice();
    refreshDeviceObject();
#endif
  }
  void setOffset(int offset) {
    Offset = offset;
#ifdef __CUDACC__
    copyToDevice();
#endif
  }

  // Reset the rotation state to the unrotated view (shift = 0), keeping the
  // data and Offset.  Used by the checkpoint loader (io/checkpoint.hh) so that
  // a restored logical view maps 1:1 onto physical indices: the container's
  // future evolution depends only on the logical view and Offset (rotate()
  // applies Offset to shift, getPrevious() indexes by Offset), never on the
  // absolute shift, so this normalization is behaviour preserving.
  void normalize() {
    shift = 0;
    refresh();
  }

  void Resize(std::size_t size) {
    if (size == count) return;
    delete[] data;
    data = new T[size]{};
    count = size;
    shift = 0;
    remainder = size;
    Offset = 0;
    refresh();
#ifdef __CUDACC__
    freeDeviceData();
    nullDeviceData();
    InitDeviceData();
#endif
  }

  ~CyclicArray() {
    delete[] data;
#ifdef __CUDACC__
    freeDeviceData();
#endif
  }

  // get more info to achieve higher performance
  std::size_t getRemainder() const { return remainder; }
  const T* getStart(int i) const { return start[i]; }
  // end get more info

  const T& operator[](std::size_t i) const {
    return (i > remainder ? start[1] : start[0])[i];
  }
  T& operator[](std::size_t i) { return (i > remainder ? start[1] : start[0])[i]; }

  inline void set(std::size_t i, T value) {
    (i > remainder ? start[1] : start[0])[i] = value;
  }
  std::size_t size() const { return count; }
  // return the pointer of ith element
  T* getdataPtr(std::size_t i = 0) {
    return i > remainder ? start[1] + i : start[0] + i;
    // return &((i > remainder ? start[1] : start[0])[i]);
  }
  const T* getdataPtr(std::size_t i = 0) const {
    return i > remainder ? start[1] + i : start[0] + i;
    // return &((i > remainder ? start[1] : start[0])[i]);
  }

  // get data before the last shift(rotate), used in bcs
  T& getPrevious(std::size_t i) {
    std::ptrdiff_t prevIndex = i + Offset;
    if (prevIndex < 0) {
      prevIndex += count;
    } else if (prevIndex >= static_cast<std::ptrdiff_t>(count)) {
      prevIndex -= count;
    }
    return (static_cast<std::size_t>(prevIndex) > remainder ? start[1]
                                                           : start[0])
      [static_cast<std::size_t>(prevIndex)];
  }
  const T& getPrevious(std::size_t i) const {
    std::ptrdiff_t prevIndex = i + Offset;
    if (prevIndex < 0) {
      prevIndex += count;
    } else if (prevIndex >= static_cast<std::ptrdiff_t>(count)) {
      prevIndex -= count;
    }
    return (static_cast<std::size_t>(prevIndex) > remainder ? start[1]
                                                           : start[0])
      [static_cast<std::size_t>(prevIndex)];
  }

  void refresh() {
    const std::ptrdiff_t n = count;
    T* const base = data;
    if (shift >= 0) {
      remainder = n - shift - 1;
      // base - remainder - 1 + n
      start[0] = base + shift;
      // base - remainder - 1
      start[1] = base - (n - shift);
    } else {
      remainder = -shift - 1;
      // base - remainder - 1 + n
      start[0] = base + (n + shift);
      // base - remainder - 1
      start[1] = base + shift;
    }
  }

  void rotate() {
    const std::ptrdiff_t n = count;
    shift -= Offset;
    if (shift >= n) {
      shift -= n;
    } else if (shift <= -n) {
      shift += n;
    }
    refresh();
  }

  void rotate(std::ptrdiff_t offset) {
    Offset = offset;
    const std::ptrdiff_t n = count;
    shift -= offset;
    if (shift >= n) {
      shift -= n;
    } else if (shift <= -n) {
      shift += n;
    }
    refresh();
  }

#ifdef __CUDACC__
  // -------------------------------------------------------------------------
  // device mirror (host-driven), same pattern as the other containers: the device
  // object aliases device copies of the scalars it needs, so the caller must
  // call copyToDevice() after any host mutation of count/shift/remainder/start
  // (Resize, Init, setOffset, rotate).  Additive only -- no existing member or
  // behaviour is changed, and none of this is compiled on the CPU path.
  // -------------------------------------------------------------------------

  T* get_devptr() { return dev_data; }
  cudev::CyclicArray<T>* get_devObj() { return dev_CyclicArray; }
  const cudev::CyclicArray<T>* get_devObj() const { return dev_CyclicArray; }

  // Re-upload the device mirror object.  The cudev mirror caches start[] inside
  // the object itself (it is NOT one of the scalars copyToDevice() pushes), so
  // a host-side normalize() -- which changes what start[] must be -- has to be
  // re-pushed this way.  See io/checkpoint.hh.
  void refreshDeviceObject() {
#ifdef __CUDACC__
    if (dev_CyclicArray) {
      cuda_free(dev_CyclicArray);
      dev_CyclicArray = nullptr;
    }
    constructInDevice();
#endif
  }

  void nullDeviceData() {
    dev_count = nullptr;
    dev_data = nullptr;
    dev_shift = nullptr;
    dev_remainder = nullptr;
    dev_Offset = nullptr;
    dev_CyclicArray = nullptr;
  }
  void freeDeviceData() {
    if (dev_count) cuda_free(dev_count);
    if (dev_data) cuda_free(dev_data);
    if (dev_shift) cuda_free(dev_shift);
    if (dev_remainder) cuda_free(dev_remainder);
    if (dev_Offset) cuda_free(dev_Offset);
    if (dev_CyclicArray) cuda_free(dev_CyclicArray);
  }

  void InitDeviceData() {
    dev_count = cuda_malloc<std::size_t>(1);
    dev_data = cuda_malloc<T>(count);
    dev_shift = cuda_malloc<std::ptrdiff_t>(1);
    dev_remainder = cuda_malloc<std::size_t>(1);
    dev_Offset = cuda_malloc<std::ptrdiff_t>(1);
    copyToDevice();
    constructInDevice();
  }

  // push the host scalars and the data buffer to the device
  void copyToDevice() {
    if (!dev_data) return;
    host_to_device(dev_count, &count, 1);
    host_to_device(dev_data, data, count);
    host_to_device(dev_shift, &shift, 1);
    host_to_device(dev_remainder, &remainder, 1);
    host_to_device(dev_Offset, &Offset, 1);
  }
  // pull the data buffer and the scalars back; start[] is re-derived on the host
  void copyToHost() {
    if (!dev_data) return;
    device_to_host(&count, dev_count, 1);
    device_to_host(data, dev_data, count);
    device_to_host(&shift, dev_shift, 1);
    device_to_host(&remainder, dev_remainder, 1);
    device_to_host(&Offset, dev_Offset, 1);
    refresh();
  }
  // build the device object; start[] is NOT copied -- the host pointers are
  // meaningless on the device, so the device object starts unrotated from
  // dev_data and is rotated on the device from there, exactly like
  // constructInDevice(): assemble the mirror object from device buffers
  void constructInDevice() {
    dev_CyclicArray = cuda_malloc<cudev::CyclicArray<T>>(1);
    T* dev_base = dev_data;
    cudev::CyclicArray<T> temp(count, dev_base, dev_shift, dev_remainder,
                               dev_Offset);
    host_to_device(dev_CyclicArray, &temp, 1);
  }

#endif
};



template <typename ArrayType, unsigned int D>
class GenericFieldBase {
 public:
  using array_type = ArrayType;
  using value_type = typename ArrayType::value_type;
  static constexpr unsigned int array_dim = D;
  static constexpr bool isField = true;
  static constexpr bool isCuDevField = false;

#ifdef __CUDACC__
  using cudev_ArrayType = typename ArrayType::cudev_array_type;
  using cudev_GenericArrayFieldType = cudev::GenericFieldBase<cudev_ArrayType, D>;
#endif

 private:
  // field data
  std::array<ArrayType, D> _data;
#ifdef __CUDACC__
  cudev_GenericArrayFieldType* dev_GenericArrayField = nullptr;
  cudev_ArrayType** dev_data = nullptr;
#endif

 public:
  GenericFieldBase() : _data{} {}
  GenericFieldBase(std::size_t size)
      : _data(make_array<ArrayType, D>([&]() { return ArrayType(size); })) {
    InitDeviceData();
  }
  GenericFieldBase(std::size_t size, value_type initialValue)
      : _data(make_array<ArrayType, D>([&]() { return ArrayType(size, initialValue); })) {
    InitDeviceData();
  }
  // Copy constructor
  GenericFieldBase(const GenericFieldBase& genF) : _data(genF._data) {
    InitDeviceData();
  }
  // Move constructor
  GenericFieldBase(GenericFieldBase&& genF) noexcept
      : _data(std::move(genF._data)) {
    InitDeviceData();
  }
  // Copy assignment operator
  GenericFieldBase& operator=(const GenericFieldBase& genF) {
    if (&genF == this) return *this;
    _data = genF._data;
#ifdef __CUDACC__
    ensureDeviceDataTable();
    uploadDataTable();
    refreshDeviceObject();
#endif
    return *this;
  }
  // Move assignment operator
  GenericFieldBase& operator=(GenericFieldBase&& genF) noexcept {
    if (&genF == this) return *this;
    _data = std::move(genF._data);
#ifdef __CUDACC__
    ensureDeviceDataTable();
    uploadDataTable();
    refreshDeviceObject();
#endif
    return *this;
  }

  ~GenericFieldBase() {
#ifdef __CUDACC__
    if (dev_data) cuda_free(dev_data);
    if (dev_GenericArrayField) cuda_free(dev_GenericArrayField);
#endif
  }

  void InitDeviceData() {
#ifdef __CUDACC__
    ensureDeviceDataTable();
    uploadDataTable();
    constructInDevice();
#endif
  }

#ifdef __CUDACC__
  // The device field reads _data[i] mirrors through this table; it must be
  // rebuilt from THIS object's arrays, never copied from the source object.
  void ensureDeviceDataTable() {
    if (dev_data == nullptr) dev_data = cuda_malloc<cudev_ArrayType*>(D);
  }
  void uploadDataTable() {
    cudev_ArrayType* host_Data[D];
    for (unsigned int i = 0; i < D; ++i) host_Data[i] = _data[i].get_devObj();
    host_to_device(dev_data, host_Data, D);
  }
  void copyToDevice() {
    for (unsigned int i = 0; i < D; ++i) _data[i].copyToDevice();
  }
  void copyToHost() {
    for (unsigned int i = 0; i < D; ++i) _data[i].copyToHost();
  }
  cudev_ArrayType** get_devptr() { return dev_data; }
  cudev_GenericArrayFieldType* get_devObj() { return dev_GenericArrayField; }
  void constructInDevice() {
    dev_GenericArrayField = cuda_malloc<cudev_GenericArrayFieldType>(1);
    // temp host object
    cudev_GenericArrayFieldType temp(dev_data);
    // copy to device
    host_to_device(dev_GenericArrayField, &temp, 1);
  }
  void refreshDeviceObject() {
    if (dev_GenericArrayField) {
      cuda_free(dev_GenericArrayField);
      dev_GenericArrayField = nullptr;
    }
    constructInDevice();
  }
#endif

  ArrayType& getField(std::size_t i = 0) { return _data[i]; }
  const ArrayType& getField(std::size_t i = 0) const { return _data[i]; }

  // get<i>(id): return _data[i][id];
  template <unsigned int i = 0>
  auto& get(std::size_t id) {
    return _data[i][id];
  }
  template <unsigned int i = 0>
  const auto& get(std::size_t id) const {
    return _data[i][id];
  }
  auto& get(std::size_t id, unsigned int dir) { return _data[dir][id]; }
  const auto& get(std::size_t id, unsigned int dir) const { return _data[dir][id]; }

  // get pointer to ith data in all arrays
  std::array<value_type*, D> getArray(std::size_t id = 0) {
    std::array<value_type*, D> data{};
    for (unsigned int i = 0; i < D; ++i) data[i] = _data[i].getdataPtr(id);
    return data;
  }

  template <unsigned int i = 0>
  void SetField(std::size_t id, value_type value) {
    _data[i].set(id, value);
  }
  void SetField(int i, std::size_t id, value_type value) { _data[i].set(id, value); }
  // resize each array/field
  void Resize(std::size_t size) {
    for (unsigned int i = 0; i < D; ++i) _data[i].Resize(size);
  }
  // init
  void Init(value_type value = value_type{}) {
    for (unsigned int i = 0; i < D; ++i) _data[i].Init(value);
  }
  // re-initialise in place, including each array's rotation state
  void Reset(value_type value = value_type{}) {
    for (unsigned int i = 0; i < D; ++i) _data[i].Reset(value);
  }

  static constexpr unsigned int Size() { return D; }
};

template <typename ArrayType, typename Base>
class GenericField : public GenericFieldBase<ArrayType, Base::array_dim> {
 public:
  static constexpr unsigned int array_dim = Base::array_dim;
  using array_type = ArrayType;
  using value_type = typename ArrayType::value_type;


#ifdef __CUDACC__
  using cudev_ArrayType = typename ArrayType::cudev_array_type;
  using cudev_FieldType = cudev::GenericField<cudev_ArrayType, Base>;

  cudev_FieldType* dev_GenericField = nullptr;
#endif

  GenericField() = default;
  GenericField(std::size_t size) : GenericFieldBase<ArrayType, array_dim>(size) {
    constructInDevice();
  }
  GenericField(std::size_t size, value_type initialValue)
      : GenericFieldBase<ArrayType, array_dim>(size, initialValue) {
    constructInDevice();
  }
  GenericField(const GenericField& genF)
      : GenericFieldBase<ArrayType, array_dim>(genF) {
    constructInDevice();
  }
  GenericField(GenericField&& genF) noexcept
      : GenericFieldBase<ArrayType, array_dim>(std::move(genF)) {
    constructInDevice();
  }
  GenericField& operator=(const GenericField& genF) {
    if (&genF == this) return *this;
    GenericFieldBase<ArrayType, array_dim>::operator=(genF);
    refreshDeviceObject();
    return *this;
  }
  GenericField& operator=(GenericField&& genF) noexcept {
    if (&genF == this) return *this;
    GenericFieldBase<ArrayType, array_dim>::operator=(std::move(genF));
    refreshDeviceObject();
    return *this;
  }

  ~GenericField() {
#ifdef __CUDACC__
    if (dev_GenericField) cuda_free(dev_GenericField);
#endif
  }


  void constructInDevice() {
#ifdef __CUDACC__
    dev_GenericField = cuda_malloc<cudev_FieldType>(1);
    // temp host object
    cudev_FieldType temp(this->get_devptr());
    // copy to device
    host_to_device(dev_GenericField, &temp, 1);
#endif
  }
  void refreshDeviceObject() {
#ifdef __CUDACC__
    if (dev_GenericField) {
      cuda_free(dev_GenericField);
      dev_GenericField = nullptr;
    }
    constructInDevice();
#endif
  }
#ifdef __CUDACC__
  cudev_FieldType* get_devObj() { return dev_GenericField; }
#endif
};


template <typename T, unsigned int D>
class BasicBlock;

// field data communication: interpolate from coarse to fine
template <typename FIELD, typename FloatType>
void FieldInterpolation2D(const FIELD& CField, FIELD& FField,
                          const BasicBlock<FloatType, 2>& CBlock,
                          const BasicBlock<FloatType, 2>& CBaseBlock,
                          const BasicBlock<FloatType, 2>& FBlock,
                          const BasicBlock<FloatType, 2>& FBaseBlock) {
  const FloatType Cvoxsize = CBlock.getVoxelSize();
  const FloatType Fvoxsize = FBlock.getVoxelSize();
  // get intersection
  const AABB<FloatType, 2> intsec = getIntersection(CBaseBlock, FBaseBlock);
  int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
  int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
  // get start index of intsec in CBlock
  Vector<FloatType, 2> startC = intsec.getMin() - CBlock.getMin();
  int startCx_ = static_cast<int>(std::round(startC[0] / Cvoxsize));
  int startCy_ = static_cast<int>(std::round(startC[1] / Cvoxsize));
  // shift 1 voxel to left bottom for interpolation
  int startCx = startCx_ - 1;
  int startCy = startCy_ - 1;
  // start index of intsec in FBlock
  Vector<FloatType, 2> startF = intsec.getMin() - FBlock.getMin();
  int startFx = static_cast<int>(std::round(startF[0] / Fvoxsize));
  int startFy = static_cast<int>(std::round(startF[1] / Fvoxsize));
  // shift 1 voxel to right top for interpolation
  int startFx_ = startFx + 1;
  int startFy_ = startFy + 1;

  for (unsigned int iArr = 0; iArr < FIELD::array_dim; ++iArr) {
    const auto& CArray = CField.getField(iArr);
    auto& FArray = FField.getField(iArr);
    // interpolation
    // 0
    for (int iy = 0; iy < CNy; ++iy) {
      for (int ix = 0; ix < CNx; ++ix) {
        std::size_t Cid0 = (iy + startCy) * CBlock.getNx() + ix + startCx;
        std::size_t Cid1 = Cid0 + 1;
        std::size_t Cid2 = Cid0 + CBlock.getNx();
        std::size_t Cid3 = Cid2 + 1;
        std::size_t Fid = (iy * 2 + startFy) * FBlock.getNx() + ix * 2 + startFx;
        FArray[Fid] = CArray[Cid0] * FloatType(0.0625) +
                      CArray[Cid1] * FloatType(0.1875) +
                      CArray[Cid2] * FloatType(0.1875) + CArray[Cid3] * FloatType(0.5625);
      }
    }
    // 1
    for (int iy = 0; iy < CNy; ++iy) {
      for (int ix = 0; ix < CNx; ++ix) {
        std::size_t Cid0 = (iy + startCy) * CBlock.getNx() + ix + startCx_;
        std::size_t Cid1 = Cid0 + 1;
        std::size_t Cid2 = Cid0 + CBlock.getNx();
        std::size_t Cid3 = Cid2 + 1;
        std::size_t Fid = (iy * 2 + startFy) * FBlock.getNx() + ix * 2 + startFx_;
        FArray[Fid] = CArray[Cid0] * FloatType(0.1875) +
                      CArray[Cid1] * FloatType(0.0625) +
                      CArray[Cid2] * FloatType(0.5625) + CArray[Cid3] * FloatType(0.1875);
      }
    }
    // 2
    for (int iy = 0; iy < CNy; ++iy) {
      for (int ix = 0; ix < CNx; ++ix) {
        std::size_t Cid0 = (iy + startCy_) * CBlock.getNx() + ix + startCx;
        std::size_t Cid1 = Cid0 + 1;
        std::size_t Cid2 = Cid0 + CBlock.getNx();
        std::size_t Cid3 = Cid2 + 1;
        std::size_t Fid = (iy * 2 + startFy_) * FBlock.getNx() + ix * 2 + startFx;
        FArray[Fid] = CArray[Cid0] * FloatType(0.1875) +
                      CArray[Cid1] * FloatType(0.5625) +
                      CArray[Cid2] * FloatType(0.0625) + CArray[Cid3] * FloatType(0.1875);
      }
    }
    // 3
    for (int iy = 0; iy < CNy; ++iy) {
      for (int ix = 0; ix < CNx; ++ix) {
        std::size_t Cid0 = (iy + startCy_) * CBlock.getNx() + ix + startCx_;
        std::size_t Cid1 = Cid0 + 1;
        std::size_t Cid2 = Cid0 + CBlock.getNx();
        std::size_t Cid3 = Cid2 + 1;
        std::size_t Fid = (iy * 2 + startFy_) * FBlock.getNx() + ix * 2 + startFx_;
        FArray[Fid] = CArray[Cid0] * FloatType(0.5625) +
                      CArray[Cid1] * FloatType(0.1875) +
                      CArray[Cid2] * FloatType(0.1875) + CArray[Cid3] * FloatType(0.0625);
      }
    }
  }
}

// field data communication: average from fine to coarse
template <typename FIELD, typename FloatType>
void FieldAverage2D(const FIELD& FField, FIELD& CField,
                    const BasicBlock<FloatType, 2>& FBlock,
                    const BasicBlock<FloatType, 2>& FBaseBlock,
                    const BasicBlock<FloatType, 2>& CBlock,
                    const BasicBlock<FloatType, 2>& CBaseBlock) {
  const FloatType Cvoxsize = CBlock.getVoxelSize();
  const FloatType Fvoxsize = FBlock.getVoxelSize();
  // get intersection
  const AABB<FloatType, 2> intsec = getIntersection(CBaseBlock, FBaseBlock);
  int CNx = static_cast<int>(std::round(intsec.getExtension()[0] / Cvoxsize));
  int CNy = static_cast<int>(std::round(intsec.getExtension()[1] / Cvoxsize));
  // get start index of intsec in CBlock
  Vector<FloatType, 2> startC = intsec.getMin() - CBlock.getMin();
  int startCx = static_cast<int>(std::round(startC[0] / Cvoxsize));
  int startCy = static_cast<int>(std::round(startC[1] / Cvoxsize));
  // start index of intsec in FBlock
  Vector<FloatType, 2> startF = intsec.getMin() - FBlock.getMin();
  int startFx = static_cast<int>(std::round(startF[0] / Fvoxsize));
  int startFy = static_cast<int>(std::round(startF[1] / Fvoxsize));

  for (unsigned int iArr = 0; iArr < FIELD::array_dim; ++iArr) {
    auto& CArray = CField.getField(iArr);
    const auto& FArray = FField.getField(iArr);

    for (int iy = 0; iy < CNy; ++iy) {
      for (int ix = 0; ix < CNx; ++ix) {
        std::size_t Cid = (iy + startCy) * CBlock.getNx() + ix + startCx;
        std::size_t Fid0 = (iy * 2 + startFy) * FBlock.getNx() + ix * 2 + startFx;
        std::size_t Fid1 = Fid0 + 1;
        std::size_t Fid2 = Fid0 + FBlock.getNx();
        std::size_t Fid3 = Fid2 + 1;
        CArray[Cid] =
          (FArray[Fid0] + FArray[Fid1] + FArray[Fid2] + FArray[Fid3]) * FloatType(0.25);
      }
    }
  }
}

// field data communication: simple copy
template <typename FIELD, typename FloatType>
void FieldCopy2D(const FIELD& FromField, FIELD& ToField,
                 const BasicBlock<FloatType, 2>& FromBlock,
                 const BasicBlock<FloatType, 2>& FromBaseBlock,
                 const BasicBlock<FloatType, 2>& ToBlock,
                 const BasicBlock<FloatType, 2>& ToBaseBlock) {
  const FloatType voxsize = FromBlock.getVoxelSize();
  // get intersection
  const AABB<FloatType, 2> intsec = getIntersection(FromBaseBlock, ToBaseBlock);
  int Nx = static_cast<int>(std::round(intsec.getExtension()[0] / voxsize));
  int Ny = static_cast<int>(std::round(intsec.getExtension()[1] / voxsize));
  // get start index of intsec in FromBlock
  Vector<FloatType, 2> startFrom = intsec.getMin() - FromBlock.getMin();
  int startFromx = static_cast<int>(std::round(startFrom[0] / voxsize));
  int startFromy = static_cast<int>(std::round(startFrom[1] / voxsize));
  // start index of intsec in ToBlock
  Vector<FloatType, 2> startTo = intsec.getMin() - ToBlock.getMin();
  int startTox = static_cast<int>(std::round(startTo[0] / voxsize));
  int startToy = static_cast<int>(std::round(startTo[1] / voxsize));

  for (unsigned int iArr = 0; iArr < FIELD::array_dim; ++iArr) {
    auto& ToArray = ToField.getField(iArr);
    const auto& FromArray = FromField.getField(iArr);

    for (int iy = 0; iy < Ny; ++iy) {
      for (int ix = 0; ix < Nx; ++ix) {
        std::size_t Fromid = (iy + startFromy) * FromBlock.getNx() + ix + startFromx;
        std::size_t Toid = (iy + startToy) * ToBlock.getNx() + ix + startTox;
        ToArray[Toid] = FromArray[Fromid];
      }
    }
  }
}

// average / interpolation kernels + weight tables live in
// data_struct/interpolation.h (shared with parallel/communicator.h).

template <typename FloatType, unsigned int Dim, typename ArrayType>
void Interpolation(const ArrayType& Arr, const std::vector<IntpSource<Dim>>& srcs,
                   std::size_t& srcidx,
                   std::vector<typename ArrayType::value_type>& Buffer,
                   std::size_t& Bufferidx) {
  if constexpr (Dim == 2) {
    Buffer[Bufferidx++] = getInterpolation<0, FloatType, Dim>(Arr, srcs[srcidx++]);
    Buffer[Bufferidx++] = getInterpolation<1, FloatType, Dim>(Arr, srcs[srcidx++]);
    Buffer[Bufferidx++] = getInterpolation<2, FloatType, Dim>(Arr, srcs[srcidx++]);
    Buffer[Bufferidx++] = getInterpolation<3, FloatType, Dim>(Arr, srcs[srcidx++]);
  } else if constexpr (Dim == 3) {
    Buffer[Bufferidx++] = getInterpolation<0, FloatType, Dim>(Arr, srcs[srcidx++]);
    Buffer[Bufferidx++] = getInterpolation<1, FloatType, Dim>(Arr, srcs[srcidx++]);
    Buffer[Bufferidx++] = getInterpolation<2, FloatType, Dim>(Arr, srcs[srcidx++]);
    Buffer[Bufferidx++] = getInterpolation<3, FloatType, Dim>(Arr, srcs[srcidx++]);
    Buffer[Bufferidx++] = getInterpolation<4, FloatType, Dim>(Arr, srcs[srcidx++]);
    Buffer[Bufferidx++] = getInterpolation<5, FloatType, Dim>(Arr, srcs[srcidx++]);
    Buffer[Bufferidx++] = getInterpolation<6, FloatType, Dim>(Arr, srcs[srcidx++]);
    Buffer[Bufferidx++] = getInterpolation<7, FloatType, Dim>(Arr, srcs[srcidx++]);
  }
}

template <typename FloatType, unsigned int Dim, typename ArrayType>
void Interpolation(ArrayType& Arr, const ArrayType& nArr,
                   const std::vector<IntpSource<Dim>>& sends,
                   const std::vector<std::size_t>& recvs, std::size_t& sendidx,
                   std::size_t& recvidx) {
  if constexpr (Dim == 2) {
    Arr.set(recvs[recvidx++],
            getInterpolation<0, FloatType, Dim>(nArr, sends[sendidx++]));
    Arr.set(recvs[recvidx++],
            getInterpolation<1, FloatType, Dim>(nArr, sends[sendidx++]));
    Arr.set(recvs[recvidx++],
            getInterpolation<2, FloatType, Dim>(nArr, sends[sendidx++]));
    Arr.set(recvs[recvidx++],
            getInterpolation<3, FloatType, Dim>(nArr, sends[sendidx++]));
  } else if constexpr (Dim == 3) {
    Arr.set(recvs[recvidx++],
            getInterpolation<0, FloatType, Dim>(nArr, sends[sendidx++]));
    Arr.set(recvs[recvidx++],
            getInterpolation<1, FloatType, Dim>(nArr, sends[sendidx++]));
    Arr.set(recvs[recvidx++],
            getInterpolation<2, FloatType, Dim>(nArr, sends[sendidx++]));
    Arr.set(recvs[recvidx++],
            getInterpolation<3, FloatType, Dim>(nArr, sends[sendidx++]));
    Arr.set(recvs[recvidx++],
            getInterpolation<4, FloatType, Dim>(nArr, sends[sendidx++]));
    Arr.set(recvs[recvidx++],
            getInterpolation<5, FloatType, Dim>(nArr, sends[sendidx++]));
    Arr.set(recvs[recvidx++],
            getInterpolation<6, FloatType, Dim>(nArr, sends[sendidx++]));
    Arr.set(recvs[recvidx++],
            getInterpolation<7, FloatType, Dim>(nArr, sends[sendidx++]));
  }
}


