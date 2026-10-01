/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 *
 * Array container evaluation harness: CyclicArray vs StreamMapArray on GPU.
 *
 * Scenarios (Q = 19 populations, sizes 1M / 4M cells):
 *   mem      : device memory footprint of the 19-pop field (cudaMemGetInfo)
 *   write    : linear store through set()                 (collide output)
 *   read     : linear load through operator[]             (momenta input)
 *   gather   : random-index load                          (boundary cell lists)
 *   readprev : getPrevious() sweep                        (bounce-back pattern)
 *   rotate   : pure rotate() arithmetic, 1 thread, 200k steps (stream step)
 *   transfer : copyToDevice / copyToHost (init & VTK output steps)
 *   correct  : bit-identical cross-check of the two containers after
 *              device-side rotation, against an exact-modulo reference
 *
 * Output lines are prefixed with RESULT for easy aggregation.
 */

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <random>
#include <vector>

#include "freelb.h"
#include "freelb.hh"

using T = float;

static constexpr unsigned int Q = 19;      // D3Q19 populations
static constexpr unsigned int BLOCK = 256;
static constexpr unsigned int MAX_BLOCKS = 512;  // for per-thread-accumulator kernels

#define CUDA_CHECK(call)                                                        \
  do {                                                                          \
    const cudaError_t err = (call);                                             \
    if (err != cudaSuccess) {                                                   \
      std::printf("CUDA error %s at %s:%d\n", cudaGetErrorString(err), __FILE__, \
                  __LINE__);                                                    \
      std::exit(1);                                                             \
    }                                                                           \
  } while (0)

struct EvTimer {
  cudaEvent_t beg, end;
  EvTimer() { CUDA_CHECK(cudaEventCreate(&beg)); CUDA_CHECK(cudaEventCreate(&end)); }
  void start() { CUDA_CHECK(cudaEventRecord(beg)); }
  float stop_ms() {
    CUDA_CHECK(cudaEventRecord(end));
    CUDA_CHECK(cudaEventSynchronize(end));
    float ms = 0.f;
    CUDA_CHECK(cudaEventElapsedTime(&ms, beg, end));
    return ms;
  }
};

// ---------------------------------------------------------------------------
// kernels (templated over the container)
// ---------------------------------------------------------------------------
template <typename Arr>
__global__ void writeLinearKernel(Arr* a, std::size_t n) {
  const std::size_t stride = (std::size_t)gridDim.x * blockDim.x;
  for (std::size_t i = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += stride) {
    a->set(i, T(1) + T(i % 13) * T(0.125));
  }
}

// block-level reduction: one atomicAdd per block instead of per thread
struct BlockSink {
  __device__ static void add(T* sink, T acc) {
    __shared__ T sacc;
    if (threadIdx.x == 0) sacc = T(0);
    __syncthreads();
    atomicAdd(&sacc, acc);
    __syncthreads();
    if (threadIdx.x == 0) atomicAdd(sink, sacc);
  }
};

template <typename Arr>
__global__ void readLinearKernel(const Arr* a, std::size_t n, T* sink) {
  T acc = T(0);
  const std::size_t stride = (std::size_t)gridDim.x * blockDim.x;
  for (std::size_t i = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += stride) {
    acc += (*a)[i];
  }
  BlockSink::add(sink, acc);
}

template <typename Arr>
__global__ void gatherRandomKernel(const Arr* a, const std::size_t* idx,
                                   std::size_t n, T* sink) {
  T acc = T(0);
  const std::size_t stride = (std::size_t)gridDim.x * blockDim.x;
  for (std::size_t j = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x; j < n;
       j += stride) {
    acc += (*a)[idx[j]];
  }
  BlockSink::add(sink, acc);
}

template <typename Arr>
__global__ void readPreviousKernel(Arr* a, std::size_t n, T* sink) {
  T acc = T(0);
  const std::size_t stride = (std::size_t)gridDim.x * blockDim.x;
  for (std::size_t i = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += stride) {
    acc += a->getPrevious(i);
  }
  BlockSink::add(sink, acc);
}

// pure rotate arithmetic: one thread, many steps.  The offset comes from the
// shared *Offset buffer (pushed by the host before the launch); both
// containers only expose the no-arg rotate() on the device with identical
// semantics (shift -= Offset, single-pass reduction).
template <typename Arr>
__global__ void rotateLoopKernel(Arr* a, int steps) {
  if (blockIdx.x != 0 || threadIdx.x != 0) return;
  for (int s = 0; s < steps; ++s) a->rotate();
}

// full-view download for the correctness cross-check
template <typename Arr>
__global__ void dumpViewKernel(const Arr* a, std::size_t n, T* out) {
  const std::size_t i = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) out[i] = (*a)[i];
}

// ---------------------------------------------------------------------------
// per-container test fixture: host containers own the buffers and push state,
// kernels operate on the cudev mirror objects (via cudev_array_type)
// ---------------------------------------------------------------------------
template <typename HostArr>
struct Fixture {
  using DevArr = typename HostArr::cudev_array_type;
  std::size_t n;
  std::vector<std::unique_ptr<HostArr>> pops;
  std::vector<DevArr*> dev;
  T* sink = nullptr;
  std::size_t* idx = nullptr;  // device gather indices

  explicit Fixture(std::size_t cells) : n(cells) {
    for (unsigned int q = 0; q < Q; ++q) {
      pops.emplace_back(new HostArr(n, T{}));
      for (std::size_t i = 0; i < n; ++i) {
        (*pops[q])[i] = T(q) + T(i % 17) * T(0.25) + T(1);
      }
      pops[q]->copyToDevice();
      dev.push_back(pops[q]->get_devObj());
    }
    CUDA_CHECK(cudaMalloc(&sink, sizeof(T)));
    // shuffled indices, identical across containers
    std::vector<std::size_t> h(n);
    for (std::size_t i = 0; i < n; ++i) h[i] = i;
    std::mt19937 rng(12345u);
    std::shuffle(h.begin(), h.end(), rng);
    CUDA_CHECK(cudaMalloc(&idx, n * sizeof(std::size_t)));
    CUDA_CHECK(cudaMemcpy(idx, h.data(), n * sizeof(std::size_t),
                          cudaMemcpyHostToDevice));
  }
  ~Fixture() {
    CUDA_CHECK(cudaFree(sink));
    CUDA_CHECK(cudaFree(idx));
  }

  void resetSink() { CUDA_CHECK(cudaMemset(sink, 0, sizeof(T))); }

  std::size_t blocks(std::size_t elems) const {
    return std::min<std::size_t>(MAX_BLOCKS, (elems + BLOCK - 1) / BLOCK);
  }

  // one sweep = Q kernel launches
  void runWrite() {
    const std::size_t b = blocks(n);
    for (unsigned int q = 0; q < Q; ++q)
      writeLinearKernel<<<b, BLOCK>>>(dev[q], n);
  }
  void runRead() {
    const std::size_t b = blocks(n);
    for (unsigned int q = 0; q < Q; ++q)
      readLinearKernel<<<b, BLOCK>>>(dev[q], n, sink);
  }
  void runGather() {
    const std::size_t b = blocks(n);
    for (unsigned int q = 0; q < Q; ++q)
      gatherRandomKernel<<<b, BLOCK>>>(dev[q], idx, n, sink);
  }
  void runReadPrev() {
    const std::size_t b = blocks(n);
    for (unsigned int q = 0; q < Q; ++q)
      readPreviousKernel<<<b, BLOCK>>>(dev[q], n, sink);
  }
  void runRotate(int steps, std::ptrdiff_t d) {
    // push the rotation offset into the shared device buffer; this also
    // resets the device rotation state so every rep starts identically
    for (unsigned int q = 0; q < Q; ++q) {
      pops[q]->setOffset((int)d);
      pops[q]->copyToDevice();
    }
    for (unsigned int q = 0; q < Q; ++q)
      rotateLoopKernel<<<1, 1>>>(dev[q], steps);
  }
};

// median of a small sample
float median(std::vector<float> v) {
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}

// ---------------------------------------------------------------------------
// scenarios
// ---------------------------------------------------------------------------
template <typename Arr>
void runPerf(const char* name, std::size_t n, int reps) {
  EvTimer ev;
  const double bytes = (double)n * sizeof(T);

  // warmup + timing sweeps; alternating containers is done by the caller
  {
    Fixture<Arr> fx(n);
    fx.runWrite();  // warmup
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<float> ts;
    for (int r = 0; r < reps; ++r) {
      ev.start();
      fx.runWrite();
      ts.push_back(ev.stop_ms());
    }
    std::printf("RESULT scenario=write container=%s n=%zu ms=%.3f gbps=%.1f\n",
                name, n, median(ts), Q * bytes / median(ts) / 1e6);
  }
  {
    Fixture<Arr> fx(n);
    fx.runRead();
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<float> ts;
    for (int r = 0; r < reps; ++r) {
      fx.resetSink();
      ev.start();
      fx.runRead();
      ts.push_back(ev.stop_ms());
    }
    std::printf("RESULT scenario=read container=%s n=%zu ms=%.3f gbps=%.1f\n",
                name, n, median(ts), Q * bytes / median(ts) / 1e6);
  }
  {
    Fixture<Arr> fx(n);
    fx.runGather();
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<float> ts;
    for (int r = 0; r < reps; ++r) {
      fx.resetSink();
      ev.start();
      fx.runGather();
      ts.push_back(ev.stop_ms());
    }
    std::printf("RESULT scenario=gather container=%s n=%zu ms=%.3f gbps=%.1f\n",
                name, n, median(ts), Q * bytes / median(ts) / 1e6);
  }
  {
    Fixture<Arr> fx(n);
    fx.runReadPrev();
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<float> ts;
    for (int r = 0; r < reps; ++r) {
      fx.resetSink();
      ev.start();
      fx.runReadPrev();
      ts.push_back(ev.stop_ms());
    }
    std::printf("RESULT scenario=readprev container=%s n=%zu ms=%.3f gbps=%.1f\n",
                name, n, median(ts), Q * bytes / median(ts) / 1e6);
  }
  // rotate: single launch, in-kernel loop; ms per million stream steps
  {
    Fixture<Arr> fx(n);
    const int steps = 200000;
    fx.runRotate(100, 3);  // warmup
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<float> ts;
    for (int r = 0; r < reps; ++r) {
      ev.start();
      fx.runRotate(steps, 3);
      ts.push_back(ev.stop_ms());
    }
    std::printf(
        "RESULT scenario=rotate container=%s n=%zu ms=%.3f ns_per_rotate=%.1f\n",
        name, n, median(ts), median(ts) * 1e6 / ((double)steps * Q));
  }
  // transfer: H2D and D2H, host-chrono (synchronous copies)
  {
    Fixture<Arr> fx(n);
    auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < 10; ++r)
      for (unsigned int q = 0; q < Q; ++q) fx.pops[q]->copyToDevice();
    auto t1 = std::chrono::steady_clock::now();
    for (int r = 0; r < 10; ++r)
      for (unsigned int q = 0; q < Q; ++q) fx.pops[q]->copyToHost();
    auto t2 = std::chrono::steady_clock::now();
    const double h2d = std::chrono::duration<double, std::milli>(t1 - t0).count() / 10;
    const double d2h = std::chrono::duration<double, std::milli>(t2 - t1).count() / 10;
    std::printf("RESULT scenario=transfer container=%s n=%zu h2d_ms=%.2f d2h_ms=%.2f\n",
                name, n, h2d, d2h);
  }
}

// device memory footprint of the Q-pop field
template <typename Arr>
void runMem(const char* name, std::size_t n) {
  std::size_t free0 = 0, total = 0;
  CUDA_CHECK(cudaMemGetInfo(&free0, &total));
  {
    Fixture<Arr> fx(n);
    std::size_t free1 = 0;
    CUDA_CHECK(cudaMemGetInfo(&free1, &total));
    const double mb = (free0 - free1) / 1024.0 / 1024.0;
    std::printf("RESULT scenario=mem container=%s n=%zu mb=%.1f kb_per_cell=%.1f\n",
                name, n, mb, mb * 1024.0 / n);
  }
}

// correctness: identical fill + identical device-side rotation must give
// bit-identical views in both containers, and match the exact-modulo reference
void runCorrect(std::size_t n) {
  const std::ptrdiff_t d = 3;
  const int steps = 101;
  // every rotate subtracts d from the shift: net shift = -steps*d, so the
  // view base is (n - steps*d) mod n for BOTH containers (verified identities)
  const std::ptrdiff_t total = (std::ptrdiff_t)steps * d;
  const std::size_t start_index =
      (std::size_t)(((std::ptrdiff_t)n - total % (std::ptrdiff_t)n) %
                    (std::ptrdiff_t)n);

  std::vector<T> host(n);
  for (std::size_t i = 0; i < n; ++i) host[i] = T(i % 17) * T(0.25) + T(1);

  T* dev_out = nullptr;
  CUDA_CHECK(cudaMalloc(&dev_out, n * sizeof(T)));
  std::vector<T> got(n);

  // reference: view[i] = host[(start_index + i) % n]
  auto check = [&](const std::vector<T>& v, const char* name) {
    int bad = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const T want = host[(start_index + i) % n];
      if (v[i] != want) {
        if (++bad <= 3)
          std::printf("  CORRECT MISMATCH %s i=%zu got=%g want=%g\n", name, i,
                      (double)v[i], (double)want);
      }
    }
    std::printf("RESULT scenario=correct container=%s n=%zu mismatches=%d\n", name,
                n, bad);
  };

  {
    CyclicArray<T> a(n, T{});
    for (std::size_t i = 0; i < n; ++i) a[i] = host[i];
    a.setOffset((int)d);
    a.copyToDevice();
    rotateLoopKernel<<<1, 1>>>(a.get_devObj(), steps);
    CUDA_CHECK(cudaDeviceSynchronize());
    dumpViewKernel<<<(unsigned)((n + 255) / 256), 256>>>(a.get_devObj(), n, dev_out);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(got.data(), dev_out, n * sizeof(T),
                          cudaMemcpyDeviceToHost));
    check(got, "cyclic");
    // host readback round-trip through copyToHost must agree with the device
    a.copyToHost();
    int bad = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (a[i] != got[i]) ++bad;
    std::printf("RESULT scenario=roundtrip container=cyclic n=%zu mismatches=%d\n",
                n, bad);
  }
  {
    StreamMapArray<T> a(n, T{});
    for (std::size_t i = 0; i < n; ++i) a[i] = host[i];
    a.setOffset((int)d);
    a.copyToDevice();
    rotateLoopKernel<<<1, 1>>>(a.get_devObj(), steps);
    CUDA_CHECK(cudaDeviceSynchronize());
    dumpViewKernel<<<(unsigned)((n + 255) / 256), 256>>>(a.get_devObj(), n, dev_out);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(got.data(), dev_out, n * sizeof(T),
                          cudaMemcpyDeviceToHost));
    check(got, "streammap");
    a.copyToHost();
    int bad = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (a[i] != got[i]) ++bad;
    std::printf("RESULT scenario=roundtrip container=streammap n=%zu mismatches=%d\n",
                n, bad);
  }
  CUDA_CHECK(cudaFree(dev_out));
}

int main() {
  int dev = 0;
  CUDA_CHECK(cudaGetDevice(&dev));
  cudaDeviceProp prop{};
  CUDA_CHECK(cudaGetDeviceProperties(&prop, dev));
  std::size_t free_b = 0, total_b = 0;
  CUDA_CHECK(cudaMemGetInfo(&free_b, &total_b));
  std::printf("device: %s (sm_%d%d), %.1f GB free / %.1f GB total\n", prop.name,
              prop.major, prop.minor, free_b / 1073741824.0,
              total_b / 1073741824.0);

  runCorrect(1000);  // non-page-aligned: map_count (1024) != count
  runCorrect(1024);  // page-aligned control: map_count == count

  for (std::size_t n : {std::size_t(1048576), std::size_t(4194304)}) {
    std::printf("\n===== N = %zu =====\n", n);
    // memory first (fixtures are constructed inside)
    runMem<CyclicArray<T>>("cyclic", n);
    runMem<StreamMapArray<T>>("streammap", n);
    // perf: alternate containers to cancel thermal drift
    for (int r = 0; r < 1; ++r) {
      runPerf<CyclicArray<T>>("cyclic", n, 5);
      runPerf<StreamMapArray<T>>("streammap", n, 5);
    }
  }
  return 0;
}
