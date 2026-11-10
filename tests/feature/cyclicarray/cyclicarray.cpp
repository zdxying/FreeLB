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

// behavioral test for CyclicArray, the CPU-side POP container.
//
// The container is checked against RefRing, a naive exact-modulo reference
// model, over a wide sweep of sizes, offset magnitudes and rotate counts.
// All deltas stay inside CyclicArray's documented contract (|offset| < count,
// accumulated |shift| < count); the solver only ever drives that regime.  The
// former RingBufferArray differential suite (which this replaces) proved
// CyclicArray == RingBufferArray in-contract and lives in git history.

#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;

static int g_failures = 0;
static int g_checks = 0;

// ---------------------------------------------------------------------------
// reference model: exact-modulo ring
// ---------------------------------------------------------------------------

// Semantics matched to CyclicArray's refresh():
//   view[i]        = data[(shift + i) mod n]
//   getPrevious(i) = data[(shift + (i + Offset) mod n) mod n]
//   remainder      = shift >= 0 ? n - shift - 1 : -shift - 1
// rotate(offset) sets Offset = offset and shifts the window; rotate() consumes
// the current Offset without resetting it; Resize restarts at shift 0.
template <typename T>
class RefRing {
 public:
  explicit RefRing(std::size_t n) { reset(n); }
  RefRing(std::size_t n, T v) { reset(n); Init(v); }
  RefRing(const RefRing&) = default;
  RefRing(RefRing&&) = default;
  RefRing& operator=(const RefRing&) = default;
  RefRing& operator=(RefRing&&) = default;

  std::size_t size() const { return data_.size(); }
  void Init(T v) {
    std::fill(data_.begin(), data_.end(), v);
    offset_ = 0;  // CyclicArray::Init(InitValue) defaults offset to 0
  }
  void Init(T v, int off) {
    Init(v);
    offset_ = off;
  }
  void setOffset(int off) { offset_ = off; }

  std::size_t getRemainder() const {
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(data_.size());
    return static_cast<std::size_t>(shift_ >= 0 ? n - shift_ - 1 : -shift_ - 1);
  }

  T& operator[](std::size_t i) { return data_[mod(base() + i)]; }
  const T& operator[](std::size_t i) const { return data_[mod(base() + i)]; }
  void set(std::size_t i, T value) { (*this)[i] = value; }

  T& getPrevious(std::size_t i) {
    return data_[mod(base() + mod(static_cast<std::ptrdiff_t>(i) + offset_))];
  }
  T getPrevious(std::size_t i) const {
    return data_[mod(base() + mod(static_cast<std::ptrdiff_t>(i) + offset_))];
  }

  void rotate(std::ptrdiff_t offset) {
    offset_ = static_cast<int>(offset);
    shift_ -= offset;
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(data_.size());
    shift_ %= n;
  }
  void rotate() { rotate(offset_); }

  void Resize(std::size_t n) {
    // match CyclicArray's early return: Resize to the same size keeps the
    // current rotation state instead of resetting it
    if (n == data_.size()) return;
    reset(n);
  }

 private:
  void reset(std::size_t n) {
    data_.assign(n, T{});
    shift_ = 0;
    offset_ = 0;
  }
  // window start into data_
  std::size_t base() const {
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(data_.size());
    std::ptrdiff_t b = shift_ % n;
    if (b < 0) b += n;
    return static_cast<std::size_t>(b);
  }
  std::size_t mod(std::ptrdiff_t x) const {
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(data_.size());
    x %= n;
    if (x < 0) x += n;
    return static_cast<std::size_t>(x);
  }

  std::vector<T> data_;
  std::ptrdiff_t shift_ = 0;
  int offset_ = 0;
};

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

// deterministic content
template <typename Arr>
void fill(Arr& arr, T base) {
  for (std::size_t i = 0; i < arr.size(); ++i) {
    arr[i] = base + T(i) * T(0.5) + T(i % 7) * T(0.25);
  }
}

template <typename A, typename B>
bool compareElements(const A& got, const B& want, const char* what,
                     std::size_t size, std::ptrdiff_t delta, int iter) {
  for (std::size_t i = 0; i < want.size(); ++i) {
    ++g_checks;
    if (got[i] != want[i]) {
      std::cout << "  MISMATCH [" << what << "] size=" << size
                << " delta=" << delta << " iter=" << iter << " at i=" << i
                << ": got=" << got[i] << " want=" << want[i] << std::endl;
      ++g_failures;
      return false;
    }
  }
  return true;
}

template <typename A, typename B>
bool comparePrevious(const A& got, const B& want, const char* what,
                     std::size_t size, std::ptrdiff_t delta, int iter) {
  for (std::size_t i = 0; i < want.size(); ++i) {
    ++g_checks;
    if (got.getPrevious(i) != want.getPrevious(i)) {
      std::cout << "  MISMATCH [" << what << "] size=" << size
                << " delta=" << delta << " iter=" << iter << " at i=" << i
                << ": got=" << got.getPrevious(i)
                << " want=" << want.getPrevious(i) << std::endl;
      ++g_failures;
      return false;
    }
  }
  return true;
}

// ---------------------------------------------------------------------------
// 1. construction / Init / size
// ---------------------------------------------------------------------------
void testConstruction() {
  std::cout << "[1] construction, Init, size" << std::endl;
  const std::size_t sizes[] = {0, 1, 2, 3, 5, 8, 9, 16, 17, 100};

  for (std::size_t n : sizes) {
    CyclicArray<T> c(n, T(7));
    RefRing<T> m(n, T(7));
    ++g_checks;
    if (c.size() != m.size()) {
      std::cout << "  MISMATCH size for n=" << n << ": cyclic=" << c.size()
                << " model=" << m.size() << std::endl;
      ++g_failures;
    }
    // default-constructed Init(InitValue) path, as used by GenericFieldBase
    CyclicArray<T> c2(n);
    RefRing<T> m2(n);
    c2.Init(T(3));
    m2.Init(T(3));
    compareElements(c2, m2, "Init", n, 0, 0);
  }
}

// ---------------------------------------------------------------------------
// 2. element access after construction
// ---------------------------------------------------------------------------
void testAccess() {
  std::cout << "[2] element access (operator[], set, getdataPtr)" << std::endl;
  const std::size_t sizes[] = {1, 2, 7, 8, 64, 100};

  for (std::size_t n : sizes) {
    CyclicArray<T> c(n);
    RefRing<T> m(n);
    fill(c, T(1));
    fill(m, T(1));
    compareElements(c, m, "fill-read", n, 0, 0);

    // set() must write the logical slot
    for (std::size_t i = 0; i < n; ++i) {
      c.set(i, T(i) * T(2));
      m.set(i, T(i) * T(2));
    }
    compareElements(c, m, "set", n, 0, 0);

    // getdataPtr must address the logical element of the current view
    for (std::size_t i = 0; i < n; ++i) {
      ++g_checks;
      if (*c.getdataPtr(i) != m[i]) {
        std::cout << "  MISMATCH getdataPtr n=" << n << " i=" << i << std::endl;
        ++g_failures;
      }
    }
    ++g_checks;
    if (n > 0 && *c.getdataPtr() != m[0]) {
      std::cout << "  MISMATCH getdataPtr() n=" << n << std::endl;
      ++g_failures;
    }
  }
}

// ---------------------------------------------------------------------------
// 3. rotate with a single offset, then read back
//
// Dense sweep of in-contract deltas (|d| < count), the regime the class
// comment defines as the compatibility contract and the one the solver
// actually drives.
// ---------------------------------------------------------------------------
void testRotateSingle() {
  std::cout << "[3] rotate(offset) single-shot (in-range deltas)" << std::endl;
  const std::size_t sizes[] = {5, 8, 16, 17, 64, 100, 1000};

  for (std::size_t n : sizes) {
    std::vector<std::ptrdiff_t> deltas;
    for (std::ptrdiff_t d = -std::ptrdiff_t(n) + 1; d <= std::ptrdiff_t(n) - 1; ++d) {
      deltas.push_back(d);
    }
    for (std::ptrdiff_t d : deltas) {
      CyclicArray<T> c(n);
      RefRing<T> m(n);
      for (std::size_t i = 0; i < n; ++i) {
        c[i] = T(i + 1);
        m[i] = T(i + 1);
      }
      c.rotate(d);
      m.rotate(d);
      if (!compareElements(c, m, "after-rotate", n, d, 0)) return;
      if (!comparePrevious(c, m, "after-rotate", n, d, 0)) return;
      ++g_checks;
      if (m.getRemainder() != c.getRemainder()) {
        std::cout << "  MISMATCH getRemainder n=" << n << " d=" << d
                  << ": model=" << m.getRemainder()
                  << " cyclic=" << c.getRemainder() << std::endl;
        ++g_failures;
        return;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// 4. repeated rotations, both the no-arg and offset forms
//
// Deltas and iteration counts chosen so the accumulated shift stays inside
// CyclicArray's well-behaved range (|shift| < count).
// ---------------------------------------------------------------------------
void testRotateRepeated() {
  std::cout << "[4] rotate repeated / interleaved forms (in-range)" << std::endl;
  const std::size_t sizes[] = {4, 8, 16, 64, 1000};

  for (std::size_t n : sizes) {
    std::vector<std::ptrdiff_t> deltas;
    for (std::ptrdiff_t d : {std::ptrdiff_t(1), std::ptrdiff_t(-1),
                             std::ptrdiff_t(3), std::ptrdiff_t(-5)}) {
      if (std::size_t(d < 0 ? -d : d) < n) deltas.push_back(d);
    }

    // (a) repeated rotate(offset), sign alternating so the running shift stays
    //     bounded
    for (std::ptrdiff_t d : deltas) {
      CyclicArray<T> c(n);
      RefRing<T> m(n);
      for (std::size_t i = 0; i < n; ++i) {
        c[i] = T(i + 1);
        m[i] = T(i + 1);
      }
      for (int it = 0; it < 200; ++it) {
        const std::ptrdiff_t step = (it % 2 == 0) ? d : -d;
        c.rotate(step);
        m.rotate(step);
        if (!compareElements(c, m, "repeated", n, step, it)) return;
        if (!comparePrevious(c, m, "repeated", n, step, it)) return;
        ++g_checks;
        if (m.getRemainder() != c.getRemainder()) {
          std::cout << "  MISMATCH getRemainder repeated n=" << n
                    << " step=" << step << " iter=" << it << ": model="
                    << m.getRemainder() << " cyclic=" << c.getRemainder()
                    << std::endl;
          ++g_failures;
          return;
        }
      }
    }

    // (b) no-arg rotate() consuming the offset set by a previous rotate(offset)
    for (std::ptrdiff_t d : deltas) {
      CyclicArray<T> c(n);
      RefRing<T> m(n);
      for (std::size_t i = 0; i < n; ++i) {
        c[i] = T(i + 1);
        m[i] = T(i + 1);
      }
      for (int it = 0; it < 100; ++it) {
        const std::ptrdiff_t step = (it % 2 == 0) ? d : -d;
        c.rotate(step);  // sets Offset and rotates
        m.rotate(step);
        c.rotate();      // consumes Offset
        m.rotate();
        if (!compareElements(c, m, "noarg", n, step, it)) return;
        if (!comparePrevious(c, m, "noarg", n, step, it)) return;
      }
    }

    // (c) mixed rotation sizes, mimicking the per-direction displacements of a
    //     D3Q19 stream
    CyclicArray<T> c(n);
    RefRing<T> m(n);
    for (std::size_t i = 0; i < n; ++i) {
      c[i] = T(i + 1);
      m[i] = T(i + 1);
    }
    const std::ptrdiff_t lim = std::ptrdiff_t(n) / 8 + 1;
    std::uint32_t seed = 12345u;
    for (int it = 0; it < 500; ++it) {
      seed = seed * 1664525u + 1013904223u;
      const std::ptrdiff_t d =
        static_cast<std::ptrdiff_t>(seed % (2 * lim + 1)) - lim;
      c.rotate(d);
      m.rotate(d);
      if (!compareElements(c, m, "mixed", n, d, it)) return;
      if (!comparePrevious(c, m, "mixed", n, d, it)) return;
    }
  }
}

// ---------------------------------------------------------------------------
// 4a. write-all-then-forward-rotate (regression guard)
//
// The exact sequence a solver performs:
//   1. write every logical slot once, unrotated
//   2. rotate FORWARD (positive delta)
//   3. read the view back
// ---------------------------------------------------------------------------
void testWriteThenForwardRotate() {
  std::cout << "[4a] write-all then forward rotate" << std::endl;
  const std::size_t sizes[] = {5, 8, 16, 17, 64, 100, 1000};
  for (std::size_t n : sizes) {
    const std::ptrdiff_t nc = static_cast<std::ptrdiff_t>(n);
    const std::ptrdiff_t deltas[] = {1, 2, 5, nc / 2, nc - 1};
    for (std::ptrdiff_t d : deltas) {
      CyclicArray<T> c(n);
      RefRing<T> m(n);
      for (std::size_t i = 0; i < n; ++i) {
        c[i] = T(i + 1);
        m[i] = T(i + 1);
      }
      c.rotate(d);
      m.rotate(d);
      if (!compareElements(c, m, "write-then-fwd-rotate", n, d, 0)) return;
      if (!comparePrevious(c, m, "write-then-fwd-rotate", n, d, 0)) return;
    }
  }
}

// ---------------------------------------------------------------------------
// 5. setOffset / Init with offset, then rotate()
// ---------------------------------------------------------------------------
void testOffset() {
  std::cout << "[5] setOffset and Init(value, offset)" << std::endl;
  const std::size_t sizes[] = {8, 16, 64, 100};
  // offsets kept small so the runs stay inside the well-behaved range
  const int offsets[] = {0, 1, -1, 5, -5};

  for (std::size_t n : sizes) {
    for (int off : offsets) {
      CyclicArray<T> c(n);
      RefRing<T> m(n);
      c.Init(T(2), off);
      m.Init(T(2), off);
      for (int it = 0; it < 50; ++it) {
        c.rotate();
        m.rotate();
        if (!compareElements(c, m, "setOffset", n, off, it)) return;
        if (!comparePrevious(c, m, "setOffset", n, off, it)) return;
      }

      CyclicArray<T> c2(n);
      RefRing<T> m2(n);
      c2.setOffset(off);
      m2.setOffset(off);
      for (std::size_t i = 0; i < n; ++i) {
        c2[i] = T(i + 1);
        m2[i] = T(i + 1);
      }
      for (int it = 0; it < 50; ++it) {
        c2.rotate();
        m2.rotate();
        if (!compareElements(c2, m2, "setOffset2", n, off, it)) return;
        if (!comparePrevious(c2, m2, "setOffset2", n, off, it)) return;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// 6. Resize
// ---------------------------------------------------------------------------
void testResize() {
  std::cout << "[6] Resize" << std::endl;
  const std::size_t seq[] = {10, 10, 20, 16, 5, 100, 3, 64};

  CyclicArray<T> c(seq[0]);
  RefRing<T> m(seq[0]);
  for (std::size_t n : seq) {
    c.Resize(n);
    m.Resize(n);
    ++g_checks;
    if (c.size() != m.size()) {
      std::cout << "  MISMATCH Resize size n=" << n << ": cyclic=" << c.size()
                << " model=" << m.size() << std::endl;
      ++g_failures;
    }
    ++g_checks;
    if (c.getRemainder() != m.getRemainder()) {
      std::cout << "  MISMATCH Resize remainder n=" << n << ": cyclic="
                << c.getRemainder() << " model=" << m.getRemainder()
                << std::endl;
      ++g_failures;
    }
    c.Init(T(1));
    m.Init(T(1));
    for (std::size_t i = 0; i < n; ++i) {
      c[i] = T(i + 1);
      m[i] = T(i + 1);
    }
    // after Resize both restart at shift 0, so a rotate must agree
    c.rotate(3);
    m.rotate(3);
    compareElements(c, m, "resize-rotate", n, 3, 0);
    comparePrevious(c, m, "resize-rotate", n, 3, 0);
  }
}

// ---------------------------------------------------------------------------
// 7. copy / move semantics
// ---------------------------------------------------------------------------
void testCopyMove() {
  std::cout << "[7] copy and move" << std::endl;
  const std::size_t n = 64;

  for (std::ptrdiff_t d : {std::ptrdiff_t(0), std::ptrdiff_t(5), std::ptrdiff_t(-5)}) {
    CyclicArray<T> c(n);
    RefRing<T> m(n);
    for (std::size_t i = 0; i < n; ++i) {
      c[i] = T(i + 1);
      m[i] = T(i + 1);
    }
    c.rotate(d);
    m.rotate(d);

    // copy construct: must preserve the rotated view
    CyclicArray<T> ccopy(c);
    RefRing<T> mcopy(m);
    compareElements(ccopy, mcopy, "copy-ctor", n, d, 0);
    comparePrevious(ccopy, mcopy, "copy-ctor", n, d, 0);

    // copy assign
    CyclicArray<T> cassign(n);
    RefRing<T> massign(n);
    cassign = c;
    massign = m;
    compareElements(cassign, massign, "copy-assign", n, d, 0);
    comparePrevious(cassign, massign, "copy-assign", n, d, 0);

    // move construct
    CyclicArray<T> cmoved(std::move(ccopy));
    RefRing<T> mmoved(std::move(mcopy));
    compareElements(cmoved, mmoved, "move-ctor", n, d, 0);
    comparePrevious(cmoved, mmoved, "move-ctor", n, d, 0);

    // move assign
    CyclicArray<T> cmoved2(n);
    RefRing<T> mmoved2(n);
    cmoved2 = std::move(cassign);
    mmoved2 = std::move(massign);
    compareElements(cmoved2, mmoved2, "move-assign", n, d, 0);
    comparePrevious(cmoved2, mmoved2, "move-assign", n, d, 0);
  }

  // self-assignment must be a no-op
  CyclicArray<T> c(n);
  RefRing<T> m(n);
  for (std::size_t i = 0; i < n; ++i) {
    c[i] = T(i + 1);
    m[i] = T(i + 1);
  }
  c.rotate(7);
  m.rotate(7);
  c = c;
  m = m;
  compareElements(c, m, "self-assign", n, 7, 0);
}

// ---------------------------------------------------------------------------
// 8. the real access pattern: stream then read through a cell at a given id
//    A real solver alternates positive and negative displacements across
//    directions, so the running shift oscillates rather than accumulating.
// ---------------------------------------------------------------------------
void testStreamPattern() {
  std::cout << "[8] D2Q9/D3Q19 style alternating stream pattern" << std::endl;

  using T2 = T;
  const std::size_t n = 1000;
  const std::ptrdiff_t offsets[] = {0, 1, -1, 4, -4, 10, -10, 99, -99};

  for (std::ptrdiff_t off : offsets) {
    CyclicArray<T2> c(n);
    RefRing<T2> m(n);
    for (std::size_t i = 0; i < n; ++i) {
      c[i] = T2(i + 1);
      m[i] = T2(i + 1);
    }
    // several stream steps, reading through getPrevious each time, which is how
    // bounce-back boundaries consume the pre-rotate view
    for (int step = 0; step < 60; ++step) {
      for (std::size_t i = 0; i < 200; ++i) {
        ++g_checks;
        if (c.getPrevious(i) != m.getPrevious(i)) {
          std::cout << "  MISMATCH stream getPrevious off=" << off
                    << " step=" << step << " i=" << i << std::endl;
          ++g_failures;
          return;
        }
      }
      // alternate the direction so the running shift stays bounded
      const std::ptrdiff_t d = (step % 2 == 0) ? off : -off;
      c.rotate(d);
      m.rotate(d);
      if (!compareElements(c, m, "stream", n, d, step)) return;
    }
  }
}

// ---------------------------------------------------------------------------
int main() {
  std::cout << "CyclicArray behavioral test (vs exact-modulo reference model)"
            << std::endl;
  std::cout << "sizeof(CyclicArray<T>) = " << sizeof(CyclicArray<T>) << std::endl;
  std::cout << std::endl;

  testConstruction();
  testAccess();
  testRotateSingle();
  testRotateRepeated();
  testWriteThenForwardRotate();
  testOffset();
  testResize();
  testCopyMove();
  testStreamPattern();

  std::cout << std::endl;
  std::cout << "checks: " << g_checks << ", failures: " << g_failures << std::endl;
  if (g_failures == 0) {
    std::cout << "ALL PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
  }
  return g_failures == 0 ? 0 : 1;
}
