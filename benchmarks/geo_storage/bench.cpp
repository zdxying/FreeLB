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
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// bench.cpp
// B5 storage comparison scaffold.  It measures the dense AABB table and the
// compact VoxelMap iteration on the same fill-aware decomposition, and prints
// first-order byte estimates for an active-slot field and for brick-sparse
// storage.  The latter two are models, not solvers: B5 chooses whether B6a
// should implement one of them.

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "freelb.h"
#include "freelb.hh"
#include "utils/memory_usage.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

int Thread_Num = 1;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t BulkFlag = std::uint8_t(2);
constexpr int Steps = 50;

using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
using CELL = Cell<T, LatSet, FIELDS>;
using TaskBulk =
  tmp::Key_TypePair<BulkFlag,
                    collision::BGK<moment::rhoU<CELL>,
                                   equilibrium::SecondOrder<CELL>>>;
using TaskCollection = tmp::TupleWrapper<TaskBulk>;
using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;

struct Fixture {
  const char* name;
  int n;
};

flbgeo::VoxelGeometry<T> makeRecord(const Fixture& fixture) {
  const T n = T(fixture.n);
  const AABB<T, 3> domain(Vector<T, 3>{}, Vector<T, 3>{n, n, n});
  flbgeo::VoxelGeometryBuilder<T> builder(fixture.n, fixture.n, fixture.n,
                                          domain, T(1), T(1.5));
  builder.addBox(domain, 1, 1, flbgeo::VoxelRole::Domain);
  if (std::string(fixture.name) == "bulk") {
    builder.addBox(AABB<T, 3>(Vector<T, 3>{T(4), T(4), T(4)},
                              Vector<T, 3>{T(28), T(28), T(28)}),
                   2, 2, flbgeo::VoxelRole::Solid);
  } else {
    const T thickness = std::string(fixture.name) == "thin" ? T(1) : T(4);
    // Match tests/feature/geo_thick (16-cell plates) and geo_thin (24-cell
    // plates) so the B5 numbers are the same fixtures as the B3 tests.
    const T lo = std::string(fixture.name) == "thick" ? T(8) : T(4);
    const T hi = std::string(fixture.name) == "thick" ? T(24) : T(28);
    builder.addBox(AABB<T, 3>(Vector<T, 3>{T(8), lo, lo},
                              Vector<T, 3>{T(8) + thickness, hi, hi}),
                   2, 2, flbgeo::VoxelRole::Solid);
    builder.addBox(AABB<T, 3>(Vector<T, 3>{T(20), lo, lo},
                              Vector<T, 3>{T(20) + thickness, hi, hi}),
                   3, 3, flbgeo::VoxelRole::Solid);
  }
  return builder.build();
}

std::size_t seamCells(const BlockGeometry3D<T>& geo) {
  std::size_t cells = 0;
  for (const Block3D<T>& block : geo.getBlocks()) {
    const Communicator& comm = block.getCommunicator();
    for (const SharedComm& link : comm.Comm.Comms) {
      cells += link.SendRecvCells.size() / 2;
    }
#ifdef MPI_ENABLED
    for (const DistributedComm& link : comm.MPIComm.Recvs) {
      cells += link.Cells.size();
    }
    for (const DistributedComm& link : comm.MPIComm.Sends) {
      cells += link.Cells.size();
    }
#endif
  }
  return cells;
}

struct BrickEstimate {
  std::size_t bricks{};
  std::size_t cells{};
  std::size_t bytes{};
  std::size_t seamBytes{};
  double fill{};
};

BrickEstimate estimateBricks(const flbgeo::VoxelGeometry<T>& record,
                             std::size_t active, int b) {
  const int nx = (record.getNx() + b - 1) / b;
  const int ny = (record.getNy() + b - 1) / b;
  const int nz = (record.getNz() + b - 1) / b;
  std::vector<std::uint8_t> present(static_cast<std::size_t>(nx) * ny * nz, 0);
  for (int z = 0; z < record.getNz(); ++z) {
    for (int y = 0; y < record.getNy(); ++y) {
      for (int x = 0; x < record.getNx(); ++x) {
        const std::size_t id = record.getIndex(Vector<int, 3>{x, y, z});
        if (!record.isDomain(id)) continue;
        const std::size_t bid =
          (static_cast<std::size_t>(z / b) * ny +
           static_cast<std::size_t>(y / b)) *
            nx +
          static_cast<std::size_t>(x / b);
        present[bid] = 1;
      }
    }
  }
  BrickEstimate est;
  auto at = [&](int bx, int by, int bz) -> bool {
    if (bx < 0 || bx >= nx || by < 0 || by >= ny || bz < 0 || bz >= nz) {
      return false;
    }
    return present[(static_cast<std::size_t>(bz) * ny + by) * nx + bx] != 0;
  };
  std::size_t faces = 0;
  for (int z = 0; z < nz; ++z) {
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nx; ++x) {
        if (!at(x, y, z)) continue;
        ++est.bricks;
        if (at(x + 1, y, z)) ++faces;
        if (at(x, y + 1, z)) ++faces;
        if (at(x, y, z + 1)) ++faces;
      }
    }
  }
  const std::size_t fieldBytes = LatSet::q * sizeof(T) + sizeof(T) +
                                 LatSet::d * sizeof(T);
  est.cells = est.bricks * static_cast<std::size_t>(b) * b * b;
  est.bytes = est.cells * fieldBytes;
  est.seamBytes = faces * static_cast<std::size_t>(b) * b * LatSet::q * sizeof(T);
  est.fill = est.cells ? static_cast<double>(active) / est.cells : 0.;
  return est;
}

struct Metrics {
  std::size_t blocks{};
  std::size_t allocated{};
  std::size_t owned{};
  std::size_t seam{};
  std::size_t denseBytes{};
  std::size_t slotBytes{};
  double mlups{};
  std::size_t rss{};
  std::size_t peakRss{};
};

Metrics runDense(const flbgeo::VoxelGeometry<T>& record, int maxBlocks) {
  BlockGeometryHelper3D<T> helper(
    record.getNx(), record.getNy(), record.getNz(), record.getDomain(),
    record.getVoxelSize(), 8, 1, 0, false);
  helper.getFlagField().Init(record, BulkFlag, VoidFlag);
  helper.CreateBlocks(2, 2, 2);
  helper.FillAwareDecomposition(maxBlocks, T(0), false);
  helper.LoadBalancing(1);

  BlockGeometry3D<T> geo(helper);
  BlockFieldManager<FLAG, T, 3> flagFM(geo, VoidFlag);
  flbgeo::writeVoxelFlags(record, flagFM, BulkFlag, VoidFlag);

  BaseConverter<T> conv(LatSet::cs2);
  conv.SimplifiedConverterFromRT(record.getNx(), T(0.1), T(0.56));
  ValuePack initValues(conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> lattice(geo, initValues, conv);

  for (int step = 0; step < 2; ++step) {
    lattice.Stream();
    lattice.ApplyCellDynamicsWithHiddenComm<NSTask>(flagFM);
  }
  const auto t0 = std::chrono::steady_clock::now();
  for (int step = 0; step < Steps; ++step) {
    lattice.Stream();
    lattice.ApplyCellDynamicsWithHiddenComm<NSTask>(flagFM);
  }
  const auto t1 = std::chrono::steady_clock::now();

  Metrics m;
  m.blocks = geo.getBlockNum();
  m.allocated = geo.getTotalCellNum();
  m.owned = geo.getTotalOwnedCellNum(flagFM, VoidFlag);
  m.seam = seamCells(geo);
  const std::size_t fieldBytes =
    LatSet::q * sizeof(T) + sizeof(T) + LatSet::d * sizeof(T);
  m.denseBytes = m.allocated * fieldBytes;
  m.slotBytes = 0;
  m.rss = flbmem::currentRssBytes();
  m.peakRss = flbmem::peakRssBytes();
  const double seconds =
    std::chrono::duration<double>(t1 - t0).count();
  m.mlups = seconds > 0. ? m.owned * Steps / seconds / 1e6 : 0.;
  return m;
}

Metrics runCompact(const flbgeo::VoxelGeometry<T>& record, int maxBlocks) {
  BlockGeometryHelper3D<T> helper(
    record.getNx(), record.getNy(), record.getNz(), record.getDomain(),
    record.getVoxelSize(), 8, 1, 0, false);
  helper.getFlagField().Init(record, BulkFlag, VoidFlag);
  helper.CreateBlocks(2, 2, 2);
  helper.FillAwareDecomposition(maxBlocks, T(0), false);
  helper.LoadBalancing(1);

  BlockGeometry3D<T> geo(helper);
  BlockFieldManager<FLAG, T, 3> flagFM(geo, VoidFlag);
  flbgeo::writeVoxelFlags(record, flagFM, BulkFlag, VoidFlag);
  geo.InitVoxelMap(flagFM, VoidFlag);

  BaseConverter<T> conv(LatSet::cs2);
  conv.SimplifiedConverterFromRT(record.getNx(), T(0.1), T(0.56));
  ValuePack initValues(conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> lattice(geo, initValues, conv);

  for (int step = 0; step < 2; ++step) {
    lattice.Stream();
    lattice.ApplyCellDynamicsWithHiddenComm<NSTask>(flagFM);
  }
  const auto t0 = std::chrono::steady_clock::now();
  for (int step = 0; step < Steps; ++step) {
    lattice.Stream();
    lattice.ApplyCellDynamicsWithHiddenComm<NSTask>(flagFM);
  }
  const auto t1 = std::chrono::steady_clock::now();

  Metrics m;
  m.blocks = geo.getBlockNum();
  m.allocated = geo.getTotalCellNum();
  m.owned = geo.getTotalOwnedCellNum(flagFM, VoidFlag);
  m.seam = seamCells(geo);
  const std::size_t fieldBytes =
    LatSet::q * sizeof(T) + sizeof(T) + LatSet::d * sizeof(T);
  m.denseBytes = m.allocated * fieldBytes;
  std::size_t aabb = 0;
  std::size_t vox = 0;
  for (const Block3D<T>& block : geo.getBlocks()) {
    const VoxelMap& map = block.getVoxelMap();
    if (!map.isInitialized()) continue;
    aabb += map.getAabbNum();
    vox += map.getVoxNum();
  }
  // First-order active-slot footprint: active slots for the fields, plus the
  // forward int32 map and the reverse size_t/owned arrays.
  m.slotBytes = vox * fieldBytes + aabb * sizeof(std::int32_t) +
                vox * (sizeof(std::size_t) + sizeof(std::uint8_t));
  m.rss = flbmem::currentRssBytes();
  m.peakRss = flbmem::peakRssBytes();
  const double seconds =
    std::chrono::duration<double>(t1 - t0).count();
  m.mlups = seconds > 0. ? m.owned * Steps / seconds / 1e6 : 0.;
  return m;
}

void report(const Fixture& fixture, const char* mode, const Metrics& m,
            const BrickEstimate& b4, const BrickEstimate& b8) {
  std::cout << "[geo_storage] fixture=" << fixture.name << " mode=" << mode
            << " blocks=" << m.blocks << " allocated=" << m.allocated
            << " owned=" << m.owned
            << " fill=" << (m.allocated ? static_cast<double>(m.owned) / m.allocated : 0.)
            << " seam_cells=" << m.seam
            << " seam_bytes=" << m.seam * LatSet::q * sizeof(T)
            << " dense_bytes=" << m.denseBytes
            << " slot_bytes=" << m.slotBytes << " mlups=" << m.mlups
            << " brick4_blocks=" << b4.bricks
            << " brick4_bytes=" << b4.bytes << " brick4_fill=" << b4.fill
            << " brick4_seam_bytes=" << b4.seamBytes
            << " brick8_blocks=" << b8.bricks
            << " brick8_bytes=" << b8.bytes << " brick8_fill=" << b8.fill
            << " brick8_seam_bytes=" << b8.seamBytes
            << " rss=" << m.rss << " peak_rss=" << m.peakRss << std::endl;
}

}  // namespace

int main() {
  for (const Fixture& fixture : {Fixture{"bulk", 32}, Fixture{"thick", 32},
                                 Fixture{"thin", 32}}) {
    const flbgeo::VoxelGeometry<T> record = makeRecord(fixture);
    const std::size_t active = record.countRole(flbgeo::VoxelRole::Domain);
    const BrickEstimate b4 = estimateBricks(record, active, 4);
    const BrickEstimate b8 = estimateBricks(record, active, 8);

    const Metrics dense = runDense(record, 24);
    report(fixture, "dense", dense, b4, b8);
    const Metrics compact = runCompact(record, 24);
    report(fixture, "compact", compact, b4, b8);
  }
  return 0;
}
