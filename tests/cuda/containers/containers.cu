/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 */

// cuda/containers: device-side checks for the cudev container wrappers.
//
// Merged from the former cuda/cyclicarray and cuda/genericarray targets, which
// both exercised the same host<->device plumbing (copyToDevice / get_devObj /
// copyToHost) -- the part that silently rots when a container layout changes.
//
//   Section A  cudev::CyclicArray  (rotate / getPrevious / set / copyToHost)
//   Section B  GenericArray / Genericvector / Data / GenericField
//
// Section A is the higher-value one: for each (count, offset) inside the
// container's contract it rotates ON the device through the mirror object and
// compares the read-back view against an exact-modulo reference model.

#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "freelb.h"
#include "freelb.hh"

#include "../../test_util.h"

using T = FLOAT;

#define CUDA_CHECK(call)                                                        \
  do {                                                                          \
    const cudaError_t err = (call);                                             \
    if (err != cudaSuccess) {                                                   \
      std::printf("CUDA error %s at %s:%d\n", cudaGetErrorString(err), __FILE__, \
                  __LINE__);                                                    \
      std::exit(1);                                                             \
    }                                                                           \
  } while (0)

// ======================== Section A: cudev::CyclicArray =====================
namespace cyclic {

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
               std::size_t index, std::size_t count) {
  for (std::size_t i = 0; i < count; ++i) {
    CHECK_EQ(got[i], host[refSlot(count, index, i)]);
  }
}

// the cudev object caches count by value; read it on the device
__global__ void mirrorSizeKernel(const cudev::CyclicArray<T>* arr, std::size_t* out) {
  *out = arr->size();
}

void run() {
  std::printf("\n--- cudev::CyclicArray ---\n");

  // Resize must rebuild the device mirror, not only the host buffer.
  {
    const std::size_t n0 = 32;
    const std::size_t n1 = 128;
    CyclicArray<T> arr(n0, T{3});
    std::size_t* dev_n = nullptr;
    CUDA_CHECK(cudaMalloc(&dev_n, sizeof(std::size_t)));

    mirrorSizeKernel<<<1, 1>>>(arr.get_devObj(), dev_n);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    std::size_t before = 0;
    CUDA_CHECK(cudaMemcpy(&before, dev_n, sizeof(before), cudaMemcpyDeviceToHost));

    arr.Resize(n1);
    mirrorSizeKernel<<<1, 1>>>(arr.get_devObj(), dev_n);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    std::size_t after = 0;
    CUDA_CHECK(cudaMemcpy(&after, dev_n, sizeof(after), cudaMemcpyDeviceToHost));

    std::printf("[resize] host=%zu mirror before=%zu after=%zu\n", arr.size(),
                before, after);
    CHECK_EQ(before, n0);
    CHECK_EQ(arr.size(), n1);
    CHECK_EQ(after, n1);
    CUDA_CHECK(cudaFree(dev_n));
  }

  // counts deliberately include non-power-of-two values; offsets are filtered
  // to CyclicArray's contract (|offset| < count)
  const std::size_t counts[] = {1, 2, 3, 5, 8, 9, 16, 17, 100, 1000};
  const std::ptrdiff_t offsets[] = {0, 1, 2, 7, -1, -2, -7, 3, -3};

  T* dev_out = nullptr;
  CUDA_CHECK(cudaMalloc(&dev_out, 1000 * sizeof(T)));

  // ------------------------------------------------------------------------
  // rotate() executed ON the device, then read back
  // ------------------------------------------------------------------------
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
      checkView(got, host, start_index, count);
    }
  }

  // ------------------------------------------------------------------------
  // getPrevious on device: after rotate(off) the pre-rotate view is content[i]
  // ------------------------------------------------------------------------
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
    checkView(got, host, 0, count);
  }

  // ------------------------------------------------------------------------
  // set() on device: write through set(), read back through operator[]
  // ------------------------------------------------------------------------
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
      CHECK_EQ(got[i], T(i) + T(1) + T(1000));
    }
  }

  // ------------------------------------------------------------------------
  // copyToHost round-trip: the device-rotated state must come back intact
  // ------------------------------------------------------------------------
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
      checkView(got, host, start_index, count);

      bool same = true;
      for (std::size_t i = 0; i < count; ++i) {
        if (arr[i] != got[i]) {
          same = false;
          break;
        }
      }
      CHECK(same);
    }
  }

  CUDA_CHECK(cudaFree(dev_out));
}

}  // namespace cyclic

// ===================== Section B: cudev generic containers ==================
namespace gen {

// ---------------------------------------------------------------- GenericArray
template <typename U>
__any__ void addGenericArrayImp(cudev::GenericArray<U>& a, std::size_t id, U value) {
  a[id] += value;
}
template <typename U>
__global__ void addGenericArray_kernel(cudev::GenericArray<U>* a, U value) {
  const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < a->size()) addGenericArrayImp(*a, idx, value);
}
template <typename U>
void addGenericArray(GenericArray<U>& a, U value) {
  const unsigned int bs = THREADS_PER_BLOCK;
  const unsigned int bn = static_cast<unsigned int>((a.size() + bs - 1) / bs);
  addGenericArray_kernel<<<bn, bs>>>(a.get_devObj(), value);
}

// ------------------------------------------------------------- Genericvector
template <typename U>
__any__ void addvectorImp(cudev::Genericvector<U>& a, std::size_t id, U value) {
  a[id] += value;
}
template <typename U>
__global__ void addvector_kernel(cudev::Genericvector<U>* a, U value) {
  const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < a->size()) addvectorImp(*a, idx, value);
}
template <typename U>
void addvector(Genericvector<U>& a, U value) {
  const unsigned int bs = THREADS_PER_BLOCK;
  const unsigned int bn = static_cast<unsigned int>((a.size() + bs - 1) / bs);
  addvector_kernel<<<bn, bs>>>(a.get_devObj(), value);
}

// ------------------------------------------------------------------ Data
template <typename U, typename Base>
__global__ void addData_kernel(cudev::Data<U, Base>* a, U value) {
  a->get() += value;
}
template <typename U, typename Base>
void addData(Data<U, Base>& a, U value) {
  addData_kernel<<<1, 1>>>(a.get_devObj(), value);
}

// ------------------------------------------------------------- GenericField
template <typename ArrayType, typename U, typename Base>
__any__ void addGenericFieldImp(cudev::GenericField<ArrayType, Base>& a, std::size_t id,
                                U value) {
  a.get(id) += value;
}
template <typename ArrayType, typename U, typename Base>
__global__ void addGenericField_kernel(cudev::GenericField<ArrayType, Base>* a, U value) {
  const std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < a->getField().size()) addGenericFieldImp(*a, idx, value);
}
template <typename ArrayType, typename U, typename Base>
void addGenericField(GenericField<ArrayType, Base>& a, U value) {
  const unsigned int bs = THREADS_PER_BLOCK;
  const unsigned int bn =
      static_cast<unsigned int>((a.getField().size() + bs - 1) / bs);
  addGenericField_kernel<<<bn, bs>>>(a.get_devObj(), value);
}

// A GenericField keeps a device table of cudev array mirrors.  Each object must
// point at its own mirrors; copying the source table sends writes to the source.
template <typename Field>
cudev::GenericArray<T>* firstMirror(Field& field) {
  cudev::GenericArray<T>* entry = nullptr;
  CUDA_CHECK(cudaMemcpy(&entry, field.get_devptr(), sizeof(entry),
                        cudaMemcpyDeviceToHost));
  return entry;
}

void testFieldMirrorOwnership() {
  std::printf("\n--- GenericField mirror ownership ---\n");
  RHO<T> src(64, T{1});
  CHECK(firstMirror(src) == src.getField(0).get_devObj());

  // Pre-initialise the targets: a default-constructed target also exposes a
  // separate null dev_count in GenericArray's copy assignment.
  RHO<T> copy(64, T{0});
  copy = src;
  std::printf("[copy] own=%d source=%d field_mirror_shared=%d\n",
              firstMirror(copy) == copy.getField(0).get_devObj(),
              firstMirror(copy) == src.getField(0).get_devObj(),
              copy.get_devObj() == src.get_devObj());
  CHECK(firstMirror(copy) == copy.getField(0).get_devObj());
  CHECK(copy.get_devObj() != src.get_devObj());

  addGenericField(copy, T{10});
  copy.copyToHost();
  src.copyToHost();
  CHECK_NEAR(copy.get(0), T{11}, 0.0);
  CHECK_NEAR(src.get(0), T{1}, 0.0);

  RHO<T> moved(64, T{0});
  moved = std::move(src);
  std::printf("[move] own=%d source=%d field_mirror_shared=%d\n",
              firstMirror(moved) == moved.getField(0).get_devObj(),
              firstMirror(moved) == src.getField(0).get_devObj(),
              moved.get_devObj() == src.get_devObj());
  CHECK(firstMirror(moved) == moved.getField(0).get_devObj());
  CHECK(moved.get_devObj() != src.get_devObj());
}

// Construct and destroy a manager: its destructor owns dev_Fields and the
// mirror object.  compute-sanitizer --leak-check reports the table/mirror if
// either is dropped.
void testManagerLifetime() {
  std::printf("\n--- BlockFieldManager lifetime ---\n");
  AABB<T, 2> domain(Vector<T, 2>{T(0), T(0)}, Vector<T, 2>{T(32), T(32)});
  {
    BlockGeometry2D<T> geo(32, 32, 4, domain, T(1), 1);
    BlockFieldManager<FLAG, T, 2> fm(geo, std::uint8_t(1));
    fm.copyToDevice();
  }
  CHECK(true);
}

void run() {
  std::printf("\n--- cudev generic containers ---\n");

  testFieldMirrorOwnership();
  testManagerLifetime();

  const std::size_t N = 10;

  // GenericArray: element-wise device add
  {
    GenericArray<T> a(N, T{1});
    addGenericArray(a, T{1});
    a.copyToHost();
    for (std::size_t i = 0; i < N; ++i) CHECK_NEAR(a[i], T{2}, 0.0);
  }

  // Genericvector: element-wise device add
  {
    Genericvector<T> v(N, T{3});
    addvector(v, T{2});
    v.copyToHost();
    for (std::size_t i = 0; i < N; ++i) CHECK_NEAR(v[i], T{5}, 0.0);
  }

  // Data (CONSTRHO): single-value device add
  {
    CONSTRHO<T> d(1, T{4});
    addData(d, T{1});
    d.copyToHost();
    CHECK_NEAR(d.get(), T{5}, 0.0);
  }

  // GenericField (VELOCITY): per-cell vector device add
  {
    VELOCITY<T, 3> vf(N, Vector<T, 3>{T{1}, T{2}, T{3}});
    addGenericField(vf, Vector<T, 3>{T{1}, T{1}, T{1}});
    vf.copyToHost();
    for (std::size_t i = 0; i < N; ++i) {
      CHECK_NEAR(vf.get(i)[0], T{2}, 0.0);
      CHECK_NEAR(vf.get(i)[1], T{3}, 0.0);
      CHECK_NEAR(vf.get(i)[2], T{4}, 0.0);
    }
  }
}

}  // namespace gen

int main() {
  int dev = 0;
  CUDA_CHECK(cudaGetDevice(&dev));
  cudaDeviceProp prop{};
  CUDA_CHECK(cudaGetDeviceProperties(&prop, dev));
  std::printf("device: %s (sm_%d%d)\n", prop.name, prop.major, prop.minor);
  std::printf("cudev::CyclicArray<T> size = %zu\n", sizeof(cudev::CyclicArray<T>));

  cyclic::run();
  gen::run();

  return flbtest::summary("cuda/containers");
}
