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

// device-side check for cudev::CyclicArray and its host device-mirror hooks.
//
// For each (count, offset) inside CyclicArray's contract (|offset| < count):
// fill a host CyclicArray, push it with copyToDevice(), rotate ON the device
// through the mirror object, read back, and compare against the exact-modulo
// reference model.  Also exercises getPrevious() and set() on the device, and
// the copyToHost() round-trip.

#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "freelb.h"
#include "freelb.hh"

using T = float;

static int g_checks = 0;
static int g_failures = 0;

#define CUDA_CHECK(call)                                                        \
  do {                                                                          \
    const cudaError_t err = (call);                                             \
    if (err != cudaSuccess) {                                                   \
      std::printf("CUDA error %s at %s:%d\n", cudaGetErrorString(err), __FILE__, \
                  __LINE__);                                                    \
      std::exit(1);                                                             \
    }                                                                           \
  } while (0)

// Reference model of CyclicArray addressing: after a rotate sequence the view
// is content[(shift + i) mod count], with the window start below.
__host__ __device__ inline std::size_t refSlot(std::size_t count, std::size_t index,
                                               std::size_t i) {
  return (index + i) % count;
}

// rotate on the device through the mirror object
__global__ void rotateKernel(cudev::CyclicArray<T>* arr, std::ptrdiff_t offset) {
  arr->rotate(offset);
}

// every thread reads one element of the current view
__global__ void readKernel(const cudev::CyclicArray<T>* arr, T* out) {
  const std::size_t n = arr->size();
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    out[i] = (*arr)[i];
  }
}

// read through getPrevious (the view at position i + Offset)
__global__ void readPreviousKernel(cudev::CyclicArray<T>* arr, T* out) {
  const std::size_t n = arr->size();
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    out[i] = arr->getPrevious(i);
  }
}

// write every element of the current view through set()
__global__ void writeKernel(cudev::CyclicArray<T>* arr) {
  const std::size_t n = arr->size();
  const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    arr->set(i, (*arr)[i] + T(1000));
  }
}

// upload a freshly filled host array and return the mirror object pointer
cudev::CyclicArray<T>* upload(const std::vector<T>& host, CyclicArray<T>& arr) {
  for (std::size_t k = 0; k < host.size(); ++k) arr[k] = host[k];
  arr.copyToDevice();
  return arr.get_devObj();
}

std::vector<T> download(const cudev::CyclicArray<T>* dev, T* dev_out,
                        std::size_t count) {
  const int threads = 128;
  const int blocks = static_cast<int>((count + threads - 1) / threads);
  readKernel<<<blocks, threads>>>(dev, dev_out);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
  std::vector<T> got(count);
  CUDA_CHECK(cudaMemcpy(got.data(), dev_out, count * sizeof(T),
                        cudaMemcpyDeviceToHost));
  return got;
}

void checkView(const std::vector<T>& got, const std::vector<T>& host,
               std::size_t index, const char* what, std::size_t count,
               std::ptrdiff_t off) {
  for (std::size_t i = 0; i < count; ++i) {
    ++g_checks;
    const T want = host[refSlot(count, index, i)];
    if (got[i] != want) {
      std::printf("  %s MISMATCH count=%zu off=%td i=%zu: got=%g want=%g\n", what,
                  count, off, i, static_cast<double>(got[i]),
                  static_cast<double>(want));
      ++g_failures;
      return;
    }
  }
}

int main() {
  int dev = 0;
  CUDA_CHECK(cudaGetDevice(&dev));
  cudaDeviceProp prop{};
  CUDA_CHECK(cudaGetDeviceProperties(&prop, dev));
  std::printf("device: %s (sm_%d%d)\n", prop.name, prop.major, prop.minor);
  std::printf("cudev::CyclicArray<T> size = %zu\n", sizeof(cudev::CyclicArray<T>));
  std::printf("\n");

  // counts deliberately include non-power-of-two values; offsets are filtered
  // to CyclicArray's contract (|offset| < count)
  const std::size_t counts[] = {1, 2, 3, 5, 8, 9, 16, 17, 100, 1000};
  const std::ptrdiff_t offsets[] = {0, 1, 2, 7, -1, -2, -7, 3, -3};

  T* dev_out = nullptr;
  CUDA_CHECK(cudaMalloc(&dev_out, 1000 * sizeof(T)));

  // -------------------------------------------------------------------------
  // rotate() executed ON the device, then read back
  // -------------------------------------------------------------------------
  for (std::size_t count : counts) {
    for (std::ptrdiff_t off : offsets) {
      if (static_cast<std::size_t>(off < 0 ? -off : off) >= count) continue;

      // content: slot k holds k + 1; window start after rotate(off) is
      // (0 - off) mod count
      std::vector<T> host(count, T{});
      for (std::size_t k = 0; k < count; ++k) host[k] = T(k) + T(1);
      const std::size_t start_index =
          static_cast<std::size_t>((count - static_cast<std::size_t>(off)) % count);

      CyclicArray<T> arr(count);
      cudev::CyclicArray<T>* dev_arr = upload(host, arr);

      rotateKernel<<<1, 1>>>(dev_arr, off);
      CUDA_CHECK(cudaGetLastError());
      CUDA_CHECK(cudaDeviceSynchronize());

      const std::vector<T> got = download(dev_arr, dev_out, count);
      checkView(got, host, start_index, "rotate", count, off);
    }
  }

  // -------------------------------------------------------------------------
  // getPrevious on device: after rotate(off) the pre-rotate view is content[i]
  // -------------------------------------------------------------------------
  {
    const std::size_t count = 64;
    std::vector<T> host(count, T{});
    for (std::size_t k = 0; k < count; ++k) host[k] = T(k) + T(1);

    CyclicArray<T> arr(count);
    cudev::CyclicArray<T>* dev_arr = upload(host, arr);

    rotateKernel<<<1, 1>>>(dev_arr, 5);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    readPreviousKernel<<<1, 64>>>(dev_arr, dev_out);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<T> got(count);
    CUDA_CHECK(cudaMemcpy(got.data(), dev_out, count * sizeof(T),
                          cudaMemcpyDeviceToHost));
    checkView(got, host, 0, "getPrevious", count, 5);
  }

  // -------------------------------------------------------------------------
  // set() on device: write through set(), read back through operator[]
  // -------------------------------------------------------------------------
  {
    const std::size_t count = 100;
    std::vector<T> host(count, T{});
    for (std::size_t k = 0; k < count; ++k) host[k] = T(k) + T(1);

    CyclicArray<T> arr(count);
    cudev::CyclicArray<T>* dev_arr = upload(host, arr);

    writeKernel<<<1, 128>>>(dev_arr);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    const std::vector<T> got = download(dev_arr, dev_out, count);
    for (std::size_t i = 0; i < count; ++i) {
      ++g_checks;
      const T want = T(i) + T(1) + T(1000);
      if (got[i] != want) {
        std::printf("  SET MISMATCH i=%zu got=%g want=%g\n", i,
                    static_cast<double>(got[i]), static_cast<double>(want));
        ++g_failures;
        break;
      }
    }
  }

  // -------------------------------------------------------------------------
  // copyToHost round-trip: the device-rotated state must come back intact
  // -------------------------------------------------------------------------
  for (std::size_t count : {std::size_t(5), std::size_t(17), std::size_t(100)}) {
    for (std::ptrdiff_t off : {std::ptrdiff_t(1), std::ptrdiff_t(-3)}) {
      std::vector<T> host(count, T{});
      for (std::size_t k = 0; k < count; ++k) host[k] = T(k) + T(1);
      const std::size_t start_index =
          static_cast<std::size_t>((count - static_cast<std::size_t>(off)) % count);

      CyclicArray<T> arr(count);
      cudev::CyclicArray<T>* dev_arr = upload(host, arr);
      rotateKernel<<<1, 1>>>(dev_arr, off);
      CUDA_CHECK(cudaGetLastError());
      CUDA_CHECK(cudaDeviceSynchronize());

      arr.copyToHost();
      const std::vector<T> got = download(dev_arr, dev_out, count);
      checkView(got, host, start_index, "roundtrip-device", count, off);
      ++g_checks;
      for (std::size_t i = 0; i < count; ++i) {
        if (arr[i] != got[i]) {
          std::printf("  ROUNDTRIP MISMATCH count=%zu off=%td i=%zu: host=%g "
                      "device=%g\n", count, off, i, static_cast<double>(arr[i]),
                      static_cast<double>(got[i]));
          ++g_failures;
          break;
        }
      }
    }
  }

  CUDA_CHECK(cudaFree(dev_out));

  std::printf("\nchecks: %d, failures: %d\n", g_checks, g_failures);
  std::printf("%s\n", g_failures ? "FAILED" : "ALL PASSED");
  return g_failures ? 1 : 0;
}
