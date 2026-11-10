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

// patch_bc.cpp
// Boundary patches now drive the task/BC flags: the recipe carries two Wall
// patches (a static obstacle and a moving plate) and geometry_build writes them
// into the .flbgeo.  The test covers:
//   - the file round-trip keeps both patches (entity, kind, cells)
//   - the kind mapping paints both, the per-entity override splits them
//   - task dispatch sees each flag (static and moving bounce-back)
//   - an unmapped kind is counted, not silently ignored and not painted

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include "freelb.h"
#include "freelb.hh"
#include "../../test_util.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

namespace {

constexpr std::uint8_t VoidFlag = std::uint8_t(1);
constexpr std::uint8_t BulkFlag = std::uint8_t(2);
constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
constexpr std::uint8_t MovingWallFlag = std::uint8_t(8);

using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;

// task dispatch probe: counts cells whose flag carries each wall bit
struct CountWalls {
  static inline std::size_t bounceback = 0;
  static inline std::size_t moving = 0;
  template <typename CELL>
  static void Execute(std::uint8_t flag, CELL&) {
    if (util::isFlag(flag, BouncebackFlag)) ++bounceback;
    if (util::isFlag(flag, MovingWallFlag)) ++moving;
  }
};

}  // namespace

int main() {
  flbgeo::GeometryFile<T> file = flbgeo::readGeometryFile<T>("patch_bc.flbgeo");
  const flbgeo::VoxelGeometry<T>& record = file.Geometry;

  CHECK_EQ(record.getPatchCount(), std::size_t{2});
  CHECK_EQ(record.getPatch(0).Entity, std::uint16_t{2});
  CHECK_EQ(record.getPatch(1).Entity, std::uint16_t{3});
  for (std::size_t p = 0; p < 2; ++p) {
    CHECK_EQ(static_cast<int>(record.getPatch(p).Kind),
             static_cast<int>(flbgeo::PatchKind::Wall));
  }
  CHECK(!record.getPatchCells(0).empty());
  CHECK(!record.getPatchCells(1).empty());
  std::unordered_set<std::size_t> obstacleCells(record.getPatchCells(0).begin(),
                                                record.getPatchCells(0).end());
  std::unordered_set<std::size_t> plateCells(record.getPatchCells(1).begin(),
                                             record.getPatchCells(1).end());

  BlockGeometry3D<T> geo(record.getNx(), record.getNy(), record.getNz(), 1,
                         record.getDomain(), record.getVoxelSize(), 1);
  BlockFieldManager<FLAG, T, 3> flagFM(geo, VoidFlag);
  flbgeo::writeVoxelFlags(record, flagFM, BulkFlag, VoidFlag);

  // ---- kind mapping paints both; the entity override splits them ----
  flbbc::PatchFlagMap map;
  map.wall = BouncebackFlag;
  std::unordered_map<std::uint16_t, std::uint8_t> byEntity;
  byEntity[3] = MovingWallFlag;
  const flbbc::PatchBindStats stats =
    flbbc::bindPatches(record, flagFM, map, BulkFlag, &byEntity);
  CHECK(stats.wall > 0);
  CHECK_EQ(stats.unmappedPatches, std::size_t{0});
  CHECK_EQ(stats.unmappedCells, std::size_t{0});
  CHECK(stats.interiorCells > 0);

  std::size_t obstaclePainted = 0;
  std::size_t platePainted = 0;
  std::size_t wrong = 0;
  flagFM.forEach([&](auto& bf, std::size_t id) {
    const std::size_t rec = record.findCellId(bf.getBlock().getLoc_t(id));
    if (rec == flbgeo::kInvalidCell || !record.isDomain(rec)) return;
    const std::uint8_t flag = bf.getField(0)[id];
    const bool bb = util::isFlag(flag, BouncebackFlag);
    const bool moving = util::isFlag(flag, MovingWallFlag);
    if (!bb && !moving) return;  // band cells not adjacent to the wall stay bulk
    const bool isObstacle = obstacleCells.count(rec) != 0;
    const bool isPlate = plateCells.count(rec) != 0;
    if (isObstacle && bb && !moving) {
      ++obstaclePainted;
    } else if (isPlate && moving && !bb) {
      ++platePainted;
    } else {
      ++wrong;
    }
  });
  CHECK_EQ(wrong, std::size_t{0});
  CHECK_EQ(stats.wall, obstaclePainted + platePainted);
  CHECK(obstaclePainted > 0);
  CHECK(platePainted > 0);

  // ---- task dispatch sees each flag ----
  BaseConverter<T> conv(LatSet::cs2);
  conv.SimplifiedConverterFromRT(record.getNx(), T(0.1), T(0.56));
  ValuePack initValues(conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> lattice(geo, initValues, conv);
  std::size_t flagBounceback = 0;
  std::size_t flagMoving = 0;
  flagFM.forEach([&](auto& bf, std::size_t id) {
    const std::uint8_t flag = bf.getField(0)[id];
    if (util::isFlag(flag, BouncebackFlag)) ++flagBounceback;
    if (util::isFlag(flag, MovingWallFlag)) ++flagMoving;
  });
  CountWalls::bounceback = 0;
  CountWalls::moving = 0;
  lattice.ApplyCellDynamics<CountWalls>(flagFM);
  CHECK_EQ(CountWalls::bounceback, flagBounceback);
  CHECK_EQ(CountWalls::moving, flagMoving);

  // ---- an unmapped kind is reported, not silently dropped ----
  BlockFieldManager<FLAG, T, 3> flagFM2(geo, VoidFlag);
  flbgeo::writeVoxelFlags(record, flagFM2, BulkFlag, VoidFlag);
  flbbc::PatchFlagMap empty;
  const flbbc::PatchBindStats skipped =
    flbbc::bindPatches(record, flagFM2, empty, BulkFlag);
  CHECK_EQ(skipped.unmappedPatches, std::size_t{2});
  CHECK_EQ(skipped.unmappedCells, obstacleCells.size() + plateCells.size());
  std::size_t leaked = 0;
  flagFM2.forEach([&](auto& bf, std::size_t id) {
    const std::uint8_t flag = bf.getField(0)[id];
    if (util::isFlag(flag, BouncebackFlag) || util::isFlag(flag, MovingWallFlag)) {
      ++leaked;
    }
  });
  CHECK_EQ(leaked, std::size_t{0});

  std::cout << "[feature/patch_bc] obstacle_cells=" << obstacleCells.size()
            << " plate_cells=" << plateCells.size()
            << " bc_cells=" << flagBounceback << " moving_cells=" << flagMoving
            << " unmapped_cells=" << skipped.unmappedCells << std::endl;

  return flbtest::summary("feature/patch_bc");
}
