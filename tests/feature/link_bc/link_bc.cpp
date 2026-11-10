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

// link_bc.cpp
// B4 first slice: build the narrow-band link table on the helical duct and run
// static interpolated bounce-back.  The test checks that the table exposes
// fractional links, then measures the total fluid mass over a short run.  The
// curved-wall accuracy baseline is a later B4 step.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

#include "freelb.h"
#include "freelb.hh"
#include "boundary/link_bc.h"

#include "../../test_util.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

int Thread_Num = 1;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t BulkFlag = std::uint8_t(2);
constexpr int Steps = 20;

using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
using CELL = Cell<T, LatSet, FIELDS>;
using BulkTask =
  tmp::Key_TypePair<BulkFlag,
                    collision::BGK<moment::rhoU<CELL>,
                                   equilibrium::SecondOrder<CELL>>>;
using TaskCollection = tmp::TupleWrapper<BulkTask>;
using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;

}  // namespace

int main() {
  StlReader<T> reader("../geo_stl/stl/helical_duct.stl", T(1), T(1), 2);
  CHECK(reader.getMesh().triangleSize() > 0);

  const T h = reader.getVoxelSize();
  const Vector<T, 3> meshMin = reader.getMin();
  const Vector<T, 3> meshMax = reader.getMax();
  const Vector<T, 3> lo = meshMin - Vector<T, 3>{h};
  const int nx = static_cast<int>(std::ceil((meshMax[0] - meshMin[0]) / h)) + 2;
  const int ny = static_cast<int>(std::ceil((meshMax[1] - meshMin[1]) / h)) + 2;
  const int nz = static_cast<int>(std::ceil((meshMax[2] - meshMin[2]) / h)) + 2;
  const AABB<T, 3> domain(
    lo, lo + Vector<T, 3>{T(nx) * h, T(ny) * h, T(nz) * h});

  flbgeo::VoxelGeometryBuilder<T> builder(nx, ny, nz, domain, h, T(1.5));
  builder.addStl(reader, 1, 1, flbgeo::VoxelRole::Domain, true);
  const flbgeo::VoxelGeometry<T> record = builder.build();
  CHECK(record.countRole(flbgeo::VoxelRole::Domain) > 0);
  CHECK(record.getBandSize() > 0);

  BlockGeometry3D<T> geo(nx, ny, nz, 1, record.getDomain(),
                         record.getVoxelSize(), 1);
  CHECK_EQ(geo.getBlockNum(), 1);
  BlockFieldManager<FLAG, T, 3> flagFM(geo, VoidFlag);
  flbgeo::writeVoxelFlags(record, flagFM, BulkFlag, VoidFlag);

  BaseConverter<T> conv(LatSet::cs2);
  conv.SimplifiedConverterFromRT(nx, T(0.1), T(0.56));
  ValuePack initValues(conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(geo, initValues, conv);

  flbbc::InterpolatedBBManager<decltype(NSLattice), decltype(flagFM)> ibb(
    NSLattice, flagFM, record, VoidFlag);
  CHECK(ibb.getLinkCount() > 0);

  std::size_t fractional{};
  std::size_t fallback{};
  std::size_t thin{};
  std::size_t withEntity{};
  for (std::size_t ib = 0; ib < geo.getBlockNum(); ++ib) {
    for (const flbbc::SurfaceLink<T>& link : ibb.getLinks(static_cast<int>(ib))) {
      CHECK(link.Fraction > T(0));
      CHECK(link.Fraction <= T(1));
      if (link.Fraction != T(0.5)) ++fractional;
      if (link.Fallback) ++fallback;
      if (link.Thin) ++thin;
      if (link.Entity != flbbc::kNoEntity) ++withEntity;
    }
  }
  CHECK(fractional > 0);
  CHECK(withEntity > 0);

  auto totalMass = [&]() {
    T sum{};
    for (std::size_t ib = 0; ib < geo.getBlockNum(); ++ib) {
      const auto& flag =
        flagFM.getBlockField(static_cast<int>(ib)).getField(0);
      const auto& rho =
        NSLattice.getBlockLat(static_cast<int>(ib))
          .template getField<RHO<T>>();
      for (std::size_t id = 0; id < flag.size(); ++id) {
        if (flag[id] != VoidFlag) sum += rho.get(id);
      }
    }
    return sum;
  };

  const T before = totalMass();
  for (int step = 0; step < Steps; ++step) {
    NSLattice.Stream();
    NSLattice.ApplyCellDynamics<NSTask>(flagFM);
    ibb.Apply();
  }
  const T after = totalMass();
  const T drift = std::abs(after - before);
  T maxU{};
  for (std::size_t ib = 0; ib < geo.getBlockNum(); ++ib) {
    const auto& flag = flagFM.getBlockField(static_cast<int>(ib)).getField(0);
    const auto& vel = NSLattice.getBlockLat(static_cast<int>(ib))
                        .template getField<VELOCITY<T, LatSet::d>>();
    for (std::size_t id = 0; id < flag.size(); ++id) {
      if (flag[id] != VoidFlag) maxU = std::max(maxU, vel.get(id).getnorm());
    }
  }

  std::cout << "[feature/link_bc] links=" << ibb.getLinkCount()
            << " fractional=" << fractional << " fallback=" << fallback
            << " thin=" << thin << " entity=" << withEntity
            << " fb_nosample=" << ibb.getFallbackNoSample()
            << " fb_nocross=" << ibb.getFallbackNoCross()
            << " fb_range=" << ibb.getFallbackOutOfRange()
            << " mass_before=" << before << " mass_after=" << after
            << " drift=" << drift << " max_u=" << maxU << std::endl;

  CHECK(drift <= T(1e-10) * std::abs(before));
  CHECK(maxU < T(1e-12));

  // Seam fixture: the same record in 2x2x2 blocks.  The link table must see
  // the same owned links through the halo, and the multi-block run must keep
  // the same total mass.
  BlockGeometryHelper3D<T> helper(nx, ny, nz, record.getDomain(),
                                  record.getVoxelSize(), 8, 1, 0, false);
  helper.CreateBlocks(2, 2, 2);
  helper.LoadBalancing(1);
  BlockGeometry3D<T> geo2(helper);
  CHECK(geo2.getBlockNum() > 1);
  BlockFieldManager<FLAG, T, 3> flagFM2(geo2, VoidFlag);
  flbgeo::writeVoxelFlags(record, flagFM2, BulkFlag, VoidFlag);
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice2(geo2, initValues, conv);
  flbbc::InterpolatedBBManager<decltype(NSLattice2), decltype(flagFM2)> ibb2(
    NSLattice2, flagFM2, record, VoidFlag);
  CHECK(ibb2.getLinkCount() > 0);
  CHECK_EQ(ibb2.getLinkCount(), ibb.getLinkCount());

  auto totalMass2 = [&]() {
    T sum{};
    for (int ib = 0; ib < geo2.getBlockNum(); ++ib) {
      const auto& block = geo2.getBlock(ib);
      const int overlap = block.getOverlap();
      const Vector<int, 3> mesh = block.getBaseBlock().getMesh();
      const auto& flag = flagFM2.getBlockField(ib).getField(0);
      const auto& rho = NSLattice2.getBlockLat(ib).template getField<RHO<T>>();
      for (int z = overlap; z < overlap + mesh[2]; ++z) {
        for (int y = overlap; y < overlap + mesh[1]; ++y) {
          for (int x = overlap; x < overlap + mesh[0]; ++x) {
            const std::size_t id = block.getIndex(Vector<int, 3>{x, y, z});
            if (flag[id] != VoidFlag) sum += rho.get(id);
          }
        }
      }
    }
    return sum;
  };

  const T before2 = totalMass2();
  CHECK_NEAR(before2, before, T(1e-12));
  for (int step = 0; step < Steps; ++step) {
    NSLattice2.Stream();
    NSLattice2.NormalFullCommunicate();
    NSLattice2.ApplyCellDynamics<NSTask>(flagFM2);
    ibb2.Apply();
  }
  const T after2 = totalMass2();
  const T drift2 = std::abs(after2 - before2);
  std::cout << "[feature/link_bc] seam blocks=" << geo2.getBlockNum()
            << " links=" << ibb2.getLinkCount()
            << " mass_before=" << before2 << " mass_after=" << after2
            << " drift=" << drift2 << std::endl;
  CHECK(drift2 <= T(1e-10) * std::abs(before2));

  return flbtest::summary("feature/link_bc");
}
