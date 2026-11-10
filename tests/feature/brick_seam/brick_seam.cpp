/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation; either
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

// brick_seam.cpp
// Two cheap regressions for the brick/storage work:
//   1. an RCB decomposition must produce a geometry that BlockGeometry3D can
//      consume directly (the basic-block pointer table is derived state);
//   2. a two-brick 1D cyclic array plus explicit seam exchange must reproduce
//      the dense cyclic view for both lattice directions.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "freelb.h"
#include "freelb.hh"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

int Thread_Num = 1;

namespace {

// RCBOptimization/TRCBOptimization replace getAllBasicBlocks(); the separate
// LoadBalancing step assigns blocks to ranks and builds the rank-local pointer
// table that BlockGeometry3D(helper) consumes.  Keep that step at the call
// site, like the private-cc examples: the block count and the process count
// may differ.
template <typename Decompose>
void CheckDecomposition(Decompose decompose) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t InsideFlag = std::uint8_t(2);

  StlReader<T> reader("../geo_stl/stl/helical_duct.stl", T(1), T(1), 2);
  BlockGeometryHelper3D<T> helper(reader, 8, 1, 1, true);
  helper.getFlagField().Init(reader, 1, InsideFlag, VoidFlag);
  helper.TagBlockCells(VoidFlag);
  decompose(helper);
  CHECK(helper.getAllBlockIndexPtrs().empty());
  helper.LoadBalancing(8);

  BlockGeometry3D<T> geo(helper);
  CHECK_EQ(geo.getBlocks().size(), std::size_t(8));
  CHECK(geo.getTotalCellNum() > 0);
  CHECK(geo.getBaseCellNum() > 0);

  std::vector<bool> seen(geo.getBlocks().size(), false);
  for (const auto& block : geo.getBlocks()) {
    CHECK(block.getN() > 0);
    const int id = block.getBlockId();
    CHECK(id >= 0 && static_cast<std::size_t>(id) < seen.size());
    if (id >= 0 && static_cast<std::size_t>(id) < seen.size()) {
      CHECK(!seen[static_cast<std::size_t>(id)]);
      seen[static_cast<std::size_t>(id)] = true;
    }
  }
}

void CheckRCB() {
  CheckDecomposition([](BlockGeometryHelper3D<T>& helper) { helper.RCBOptimization(8); });
}

void CheckTRCB() {
  CheckDecomposition([](BlockGeometryHelper3D<T>& helper) { helper.TRCBOptimization(8); });
}

// The pointer table is derived state; every block replacement clears it.
template <typename Decompose>
void CheckInvalidated(Decompose decompose) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t InsideFlag = std::uint8_t(2);

  StlReader<T> reader("../geo_stl/stl/helical_duct.stl", T(1), T(1), 2);
  BlockGeometryHelper3D<T> helper(reader, 8, 1, 1, true);
  helper.getFlagField().Init(reader, 1, InsideFlag, VoidFlag);
  helper.TagBlockCells(VoidFlag);
  decompose(helper);
  CHECK(helper.getAllBlockIndexPtrs().empty());
  helper.LoadBalancing(8);

  BlockGeometry3D<T> geo(helper);
  CHECK(geo.getTotalCellNum() > 0);
  CHECK(geo.getBaseCellNum() > 0);
}

void CheckCreateBlocks() {
  CheckInvalidated([](BlockGeometryHelper3D<T>& helper) {
    helper.CreateBlocks(true, false);
  });
}

void CheckOptimize() {
  CheckInvalidated([](BlockGeometryHelper3D<T>& helper) {
    helper.CreateBlocks(true, false);
    helper.Optimize(8, false);
  });
}

// Two bricks, dense periodic ring as the reference.  For direction +1 the
// receiver's wrap slot is the other brick's last pre-rotate value; for -1 it
// is the other brick's first pre-rotate value.  The source is read before the
// rotate and the destination written after it.
void CheckSeam1D(int offset) {
  constexpr int kBrick = 16;
  constexpr int kCount = 2 * kBrick;

  CyclicArray<T> dense(kCount);
  CyclicArray<T> a(kBrick);
  CyclicArray<T> b(kBrick);
  for (int i = 0; i < kCount; ++i) dense[i] = T(i) * T(0.5) + T(1);
  for (int i = 0; i < kBrick; ++i) {
    a[i] = dense[i];
    b[i] = dense[kBrick + i];
  }
  dense.setOffset(offset);
  a.setOffset(offset);
  b.setOffset(offset);

  for (int step = 0; step < 7; ++step) {
    const T aSource = (offset > 0) ? a[kBrick - 1] : a[0];
    const T bSource = (offset > 0) ? b[kBrick - 1] : b[0];

    dense.rotate();
    a.rotate();
    b.rotate();

    if (offset > 0) {
      a[0] = bSource;
      b[0] = aSource;
    } else {
      a[kBrick - 1] = bSource;
      b[kBrick - 1] = aSource;
    }

    for (int i = 0; i < kBrick; ++i) {
      CHECK_NEAR(a[i], dense[i], 1e-12);
      CHECK_NEAR(b[i], dense[kBrick + i], 1e-12);
    }
  }
}

}  // namespace

int main() {
  CheckRCB();
  CheckTRCB();
  CheckCreateBlocks();
  CheckOptimize();
  CheckSeam1D(+1);
  CheckSeam1D(-1);
  return flbtest::summary("feature/brick_seam");
}
