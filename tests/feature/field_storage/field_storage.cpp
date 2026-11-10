/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * FreeLB. If not, see <https://www.gnu.org/licenses/>.
 */

// field_storage.cpp
//
// Storage-level regressions for data_struct/field.h and
// BlockFieldManager's rebuild paths:
//   - BlockFieldManager::Init(): reuses the existing arrays when the block set
//     did not change
//   - CyclicArray::getPrevious() const-correctness
//
// --bench  : GenericArray/CyclicArray construction cost (the start-up write
//            sweep); run before/after a change to compare, not an assertion.
// --peak   : VmHWM delta across Init(GeoHelper); the old and new field sets
//            coexist during the transfer, so the delta is the field footprint.

#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <new>
#include <string>
#include <type_traits>
#include <vector>

#include "freelb.h"
#include "freelb.hh"

#include "../../test_util.h"

using T = FLOAT;

int Thread_Num = 1;

// Count every array allocation so a rebuild path can be checked without
// relying on the allocator handing back a different address.
static std::atomic<std::size_t> g_array_bytes{0};

void* operator new[](std::size_t n) {
  g_array_bytes.fetch_add(n, std::memory_order_relaxed);
  void* p = std::malloc(n == 0 ? 1 : n);
  if (p == nullptr) throw std::bad_alloc();
  return p;
}
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

static std::size_t vmHwmKb() {
  std::ifstream in("/proc/self/status");
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("VmHWM:", 0) == 0) {
      return static_cast<std::size_t>(std::stoul(line.substr(6)));
    }
  }
  return 0;
}

// ---------------------------------------------------------------------------
// BlockFieldManager::Init() reuse
// ---------------------------------------------------------------------------

static void testFieldManagerReuse() {
  AABB<T, 2> domain(Vector<T, 2>{T(0), T(0)}, Vector<T, 2>{T(40), T(40)});
  BlockGeometry2D<T> geo(40, 40, 4, domain, T(1), 1);
  BlockFieldManager<RHO<T>, T, 2> fm(geo, T(1));

  // mark every field so a rebuild that silently leaves the arrays unreset is
  // visible even if the allocator reuses the same addresses
  for (std::size_t b = 0; b < fm.size(); ++b) {
    fm.getBlockField(b).getField(0).Init(T(1));
  }
  const std::size_t allocated_before = g_array_bytes.load();

  fm.Init(T(7));

  const std::size_t allocated_after = g_array_bytes.load();
  std::printf("[fieldmanager] Init() array bytes allocated: %zu\n",
              allocated_after - allocated_before);
  // A changed block set must rebuild; an unchanged one must not.
  CHECK_EQ(allocated_after, allocated_before);
  for (std::size_t b = 0; b < fm.size(); ++b) {
    const auto& field = fm.getBlockField(b).getField(0);
    CHECK_NEAR(field[0], T(7), 0.0);
  }

  // A geometry rebuild must invalidate the reuse path and allocate again.
  AABB<T, 2> domain2(Vector<T, 2>{T(0), T(0)}, Vector<T, 2>{T(40), T(40)});
  BlockGeometryHelper2D<T> helper(40, 40, domain2, T(1), 10);
  helper.CreateBlocks();
  helper.LoadBalancing(1);
  geo.Init(helper);
  const std::size_t allocated_rebuild = g_array_bytes.load();
  fm.Init(T(3));
  std::printf("[fieldmanager] Init() after geometry rebuild allocated: %zu\n",
              g_array_bytes.load() - allocated_rebuild);
  CHECK(g_array_bytes.load() > allocated_rebuild);
}

// ---------------------------------------------------------------------------
// CyclicArray const-correctness
// ---------------------------------------------------------------------------

static void testConstCorrectness() {
  CyclicArray<T> arr(8);
  const CyclicArray<T>& carr = arr;
  CHECK(std::is_const_v<std::remove_reference_t<decltype(carr.getPrevious(0))>>);
  CHECK(!std::is_const_v<std::remove_reference_t<decltype(arr.getPrevious(0))>>);
}

// ---------------------------------------------------------------------------
// --bench: start-up write sweep of GenericArray / CyclicArray
// ---------------------------------------------------------------------------

template <typename Arr>
static void benchArray(const char* name, std::size_t n, int iters) {
  volatile T sink = T(0);
  const auto t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < iters; ++i) {
    Arr a(n, T(1));
    sink = sink + a[0] + a[n - 1];
  }
  const auto t1 = std::chrono::steady_clock::now();
  const double sec =
    std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
  const double bytes = static_cast<double>(n) * sizeof(T) * iters;
  std::printf("[bench] %-14s n=%zu iters=%d  %.3f s  %.2f GB/s (1 write/elem)\n",
              name, n, iters, sec, bytes / sec / 1e9);
  (void)sink;
}

static int runBench() {
  const std::size_t n = std::size_t(8) << 20;  // 64 MiB at T=double
  std::printf("[bench] sizeof(T)=%zu\n", sizeof(T));
  benchArray<GenericArray<T>>("GenericArray", n, 8);
  benchArray<CyclicArray<T>>("CyclicArray", n, 8);
  return 0;
}

// ---------------------------------------------------------------------------
// --peak: VmHWM across Init(GeoHelper), i.e. old + new field sets
// ---------------------------------------------------------------------------

static int runPeak() {
  constexpr int N = 2000;
  constexpr int BlockLen = 500;
  AABB<T, 2> domain(Vector<T, 2>{T(0), T(0)}, Vector<T, 2>{T(N), T(N)});

  BlockGeometryHelper2D<T> helper(N, N, domain, T(1), BlockLen);
  helper.CreateBlocks();
  helper.LoadBalancing(1);
  BlockGeometry2D<T> geo(helper);
  BlockFieldManager<RHO<T>, T, 2> fm(geo, T(1));

  std::size_t field_bytes = 0;
  for (const auto& block : geo.getBlocks()) field_bytes += block.getN() * sizeof(T);

  // Real rebuild order: the geometry replaces its blocks first, then the field
  // manager transfers into a freshly allocated set.
  geo.Init(helper);

  const std::size_t before = vmHwmKb();
  fm.Init(helper);
  const std::size_t after = vmHwmKb();
  const std::size_t delta = after > before ? after - before : 0;

  std::printf("[peak] field=%.1f MiB, VmHWM before=%.1f MiB after=%.1f MiB "
              "delta=%.1f MiB (%.2fx)\n",
              field_bytes / 1048576.0, before / 1024.0, after / 1024.0,
              delta / 1024.0, field_bytes ? static_cast<double>(delta) * 1024.0 /
                                              static_cast<double>(field_bytes)
                                          : 0.0);
  return 0;
}

int main(int argc, char** argv) {
  if (argc > 1 && std::strcmp(argv[1], "--bench") == 0) return runBench();
  if (argc > 1 && std::strcmp(argv[1], "--peak") == 0) return runPeak();

  std::printf("--- BlockFieldManager ---\n");
  testFieldManagerReuse();

  std::printf("--- containers ---\n");
  testConstCorrectness();

  return flbtest::summary("feature/field_storage");
}
