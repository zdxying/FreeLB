/* This file is part of FreeLB
 * 
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 * 
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 * 
 */

#pragma once

#include "boundary/basic_boundary.h"

// virtual node scheme(cell with void flag)
// impl periodic bcs is to handle the pops in the virtual nodes
// this bcs should execute after the collision step and before the streaming
// step
// template <typename T, typename LatSet>
// class FixedPeriodicBoundaryBlock final : public AbstractBoundary {
//  private:
//   // virtual cells
//   std::vector<std::size_t> VCells;
//   // virtual cells get info from Real cells
//   std::vector<std::size_t> RCells;
//   // boundary cells for macroscopic update
//   std::vector<std::size_t> BdCells;
//   // virtual cells' lattice
//   BlockLattice<T, LatSet> &VLat;
//   // real cells' lattice
//   BlockLattice<T, LatSet> &RLat;

//   // geometry
//   Block<T, LatSet::d> &Geo;
//   // geometry flag
//   GenericArray<std::uint8_t> &Field;
//   // boundary cell flag
//   std::uint8_t BdCellFlag;
//   // boundary flag
//   std::uint8_t voidFlag;
//   // AABB
//   AABB<T, LatSet::d> &Box;
//   // from AABB
//   AABB<T, LatSet::d> &FromBox;
//   std::string _name = "FixedPeriodic";

//  public:
//   FixedPeriodicBoundaryBlock(BlockLattice<T, LatSet> &vlat, BlockLattice<T, LatSet> &rlat,
//   AABB<T, LatSet::d> &box,
//                         AABB<T, LatSet::d> &frombox, std::uint8_t cellflag,
//                         std::uint8_t voidflag = std::uint8_t(1))
//       : VLat(vlat),
//         RLat(rlat),
//         Geo(vlat.getGeo()),
//         Field(vlat.getGeo().getGeoFlagField().getField()),
//         Box(box),
//         FromBox(frombox),
//         BdCellFlag(cellflag),
//         voidFlag(voidflag) {
//     Vector<T, LatSet::d> boxext = Box.getExtension();
//     Vector<int, 3> SIZE(1);
//     if constexpr (LatSet::d == 2) {
//       SIZE[0] = int(boxext[0] / Geo.getVoxelSize()) + 2;
//       SIZE[1] = int(boxext[1] / Geo.getVoxelSize()) + 2;
//     } else if constexpr (LatSet::d == 3) {
//       SIZE[0] = int(boxext[0] / Geo.getVoxelSize()) + 2;
//       SIZE[1] = int(boxext[1] / Geo.getVoxelSize()) + 2;
//       SIZE[2] = int(boxext[2] / Geo.getVoxelSize()) + 2;
//     }
//     int size = SIZE[0] * SIZE[1] * SIZE[2];
//     VCells.reserve(size);
//     RCells.reserve(size);
//     BdCells.reserve(size);
//     Geo.forEachVoxel(Box, BdCellFlag,
//                      [this](int id) { BdCells.push_back(id); });
//     // dist between two AABB
//     Vector<T, LatSet::d> dist = FromBox.getCenter() - Box.getCenter();
//     // get corrected dist
//     for (int i = 0; i < LatSet::d; ++i) {
//       if (dist[i] < T(-(1e-6))) {
//         dist[i] -= Geo.getVoxelSize();
//       } else if (dist[i] > T(1e-6)) {
//         dist[i] += Geo.getVoxelSize();
//       }
//     }
//     // extended box
//     AABB<T, LatSet::d> extBox = Box.getExtended(Geo.getVoxelSize());
//     // get virtual cells and corresponding real cells
//     Geo.forEachVoxel(extBox, voidFlag, [this, &dist](int id) {
//       if (Geo.template hasNeighborFlag<LatSet>(id, BdCellFlag)) {
//         VCells.push_back(id);
//         int idr = Geo.findCellId(Geo.getVoxel(id) + dist);
//         RCells.push_back(idr);
//       }
//     });
//   }

//   void Apply() override {
//     int size = VCells.size();
// #pragma omp parallel for num_threads(Thread_Num) schedule(static)
//     for (int i = 0; i < size; ++i) {
//       BasicPopCell<T, LatSet> vcell(VCells[i], Lat);
//       BasicPopCell<T, LatSet> rcell(RCells[i], Lat);
//       for (int k = 1; k < LatSet::q; ++k) vcell[k] = rcell[k];
//     }
//   }
//   void getinfo() override {
//     std::cout << std::setw(18) << std::left << _name << std::setw(10)
//               << std::left << BdCells.size() << std::endl;
//   }
//   void UpdateRho() override {
// #pragma omp parallel for num_threads(Thread_Num) schedule(static)
//     for (std::size_t id : BdCells) {
// BasicPopCell<T, LatSet> cell(id, Lat);
// moment::Rho<T, LatSet>::apply(cell, Lat.getRho(id));
//     }
//   }
//   void UpdateU() override {
// #pragma omp parallel for num_threads(Thread_Num) schedule(static)
//     for (std::size_t id : BdCells) {
// BasicPopCell<T, LatSet> cell(id, Lat);
      // moment::Velocity<T, LatSet>::apply(cell, Lat.getVelocity(id));
//     }
//   }
// };

template <typename T, typename LatSet>
class FixedPeriodicBoundary final : public AbstractBoundary {
 private:
  // virtual cells
  std::vector<std::size_t> VCells;
  // virtual cells get info from Real cells
  std::vector<std::size_t> RCells;
  // boundary cells for macroscopic update
  std::vector<std::size_t> BdCells;
  // lattice
  PopLattice<T, LatSet> &Lat;
  // geometry
  Geometry<T, LatSet::d> &Geo;
  // geometry flag
  GenericArray<std::uint8_t> &Field;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;
  // AABB
  AABB<T, LatSet::d> &Box;
  // from AABB
  AABB<T, LatSet::d> &FromBox;
  std::string _name = "FixedPeriodic";

 public:
  FixedPeriodicBoundary(PopLattice<T, LatSet> &lat, AABB<T, LatSet::d> &box,
                        AABB<T, LatSet::d> &frombox, std::uint8_t cellflag,
                        std::uint8_t voidflag = std::uint8_t(1))
      : Lat(lat), Geo(lat.getGeo()), Field(lat.getGeo().getGeoFlagField().getField()),
        Box(box), FromBox(frombox), BdCellFlag(cellflag), voidFlag(voidflag) {
    Vector<T, LatSet::d> boxext = Box.getExtension();
    Vector<int, 3> SIZE(1);
    if constexpr (LatSet::d == 2) {
      SIZE[0] = int(boxext[0] / Geo.getVoxelSize()) + 2;
      SIZE[1] = int(boxext[1] / Geo.getVoxelSize()) + 2;
    } else if constexpr (LatSet::d == 3) {
      SIZE[0] = int(boxext[0] / Geo.getVoxelSize()) + 2;
      SIZE[1] = int(boxext[1] / Geo.getVoxelSize()) + 2;
      SIZE[2] = int(boxext[2] / Geo.getVoxelSize()) + 2;
    }
    int size = SIZE[0] * SIZE[1] * SIZE[2];
    VCells.reserve(size);
    RCells.reserve(size);
    BdCells.reserve(size);
    Geo.forEachVoxel(Box, BdCellFlag, [this](int id) { BdCells.push_back(id); });
    // dist between two AABB
    Vector<T, LatSet::d> dist = FromBox.getCenter() - Box.getCenter();
    // get corrected dist
    for (int i = 0; i < LatSet::d; ++i) {
      if (dist[i] < T(-(1e-6))) {
        dist[i] -= Geo.getVoxelSize();
      } else if (dist[i] > T(1e-6)) {
        dist[i] += Geo.getVoxelSize();
      }
    }
    // extended box
    AABB<T, LatSet::d> extBox = Box.getExtended(Geo.getVoxelSize());
    // get virtual cells and corresponding real cells
    Geo.forEachVoxel(extBox, voidFlag, [this, &dist](int id) {
      if (Geo.template hasNeighborFlag<LatSet>(id, BdCellFlag)) {
        VCells.push_back(id);
        int idr = Geo.findCellId(Geo.getVoxel(id) + dist);
        RCells.push_back(idr);
      }
    });
  }

  void Apply() override {
    int size = VCells.size();
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int i = 0; i < size; ++i) {
      BasicPopCell<T, LatSet> vcell(VCells[i], Lat);
      BasicPopCell<T, LatSet> rcell(RCells[i], Lat);
      for (int k = 1; k < LatSet::q; ++k) vcell[k] = rcell[k];
    }
  }
  void getinfo() override {
    std::cout << std::setw(18) << std::left << _name << std::setw(10) << std::left << BdCells.size()
              << std::endl;
  }
  void UpdateRho() override {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (std::size_t id : BdCells) {
      BasicPopCell<T, LatSet> cell(id, Lat);
      moment::Rho<T, LatSet>::apply(cell, Lat.getRho(id));
    }
  }
  void UpdateU() override {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (std::size_t id : BdCells) {
      BasicPopCell<T, LatSet> cell(id, Lat);
      moment::Velocity<T, LatSet>::apply(cell, Lat.getVelocity(id));
    }
  }
};
