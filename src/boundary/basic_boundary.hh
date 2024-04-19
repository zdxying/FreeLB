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

// basic_boundary.hh

#pragma once

#include "boundary/basic_boundary.h"

// --------------------------------------------------------------------------------
// fixed boundary
// --------------------------------------------------------------------------------

template <typename T, typename LatSet, typename flagType>
FixedBoundary<T, LatSet, flagType>::FixedBoundary(BasicLattice<T, LatSet>& lat,
                                                  std::uint8_t cellflag,
                                                  std::uint8_t voidflag)
    : Lat(lat), Geo(lat.getGeo()), Field(lat.getGeo().getGeoFlagField().getField()),
      BdCellFlag(cellflag), voidFlag(voidflag) {
  Setup();
}

template <typename T, typename LatSet, typename flagType>
FixedBoundary<T, LatSet, flagType>::FixedBoundary(BasicLattice<T, LatSet>& lat,
                                                  GenericArray<flagType>& f,
                                                  std::uint8_t cellflag,
                                                  std::uint8_t voidflag)
    : Lat(lat), Geo(lat.getGeo()), Field(f), BdCellFlag(cellflag), voidFlag(voidflag) {
  Setup();
}

template <typename T, typename LatSet, typename flagType>
void FixedBoundary<T, LatSet, flagType>::addtoBd(std::size_t id) {
  BdCells.emplace_back(id, LatSet::q);
  // get reference to the last element
  FixedBdCell& fixedbdcell = BdCells.back();
  // get neighbor
  // Attention: if voidFlag is 0, DO NOT use util::isFlag
  for (int k = 1; k < LatSet::q; ++k) {
    // if (Field[Lat.getNbrId(id, k)] == voidFlag &&
    //     Field[Lat.getNbrId(id, LatSet::opp[k])] != voidFlag) {
    //   fixedbdcell.outflows.push_back(LatSet::opp[k]);
    // }

    // using util::isFlag requires voidFlag to be non-zero
    if (util::isFlag(Field[Lat.getNbrId(id, k)], voidFlag)) {
      fixedbdcell.outflows.push_back(LatSet::opp[k]);
    }
  }
}

template <typename T, typename LatSet, typename flagType>
void FixedBoundary<T, LatSet, flagType>::Setup() {
  for (std::size_t id = 0; id < Geo.getVoxelsNum(); ++id) {
    if (util::isFlag(Field[id], BdCellFlag)) addtoBd(id);
  }
}

// --------------------------------------------------------------------------------
// moving boundary
// --------------------------------------------------------------------------------

template <typename T, typename LatSet, typename flagType>
MovingBoundary<T, LatSet, flagType>::MovingBoundary(BasicLattice<T, LatSet>& lat,
                                                    std::vector<std::size_t>& ids,
                                                    std::uint8_t voidflag,
                                                    std::uint8_t cellflag)
    : Lat(lat), Geo(lat.getGeo()), Ids(ids),
      Field(lat.getGeo().getGeoFlagField().getField()), BdCellFlag(cellflag),
      voidFlag(voidflag) {}

template <typename T, typename LatSet, typename flagType>
MovingBoundary<T, LatSet, flagType>::MovingBoundary(BasicLattice<T, LatSet>& lat,
                                                    std::vector<std::size_t>& ids,
                                                    GenericArray<flagType>& f,
                                                    std::uint8_t voidflag,
                                                    std::uint8_t cellflag)
    : Lat(lat), Geo(lat.getGeo()), Ids(ids), Field(f), BdCellFlag(cellflag),
      voidFlag(voidflag) {}

template <typename T, typename LatSet, typename flagType>
void MovingBoundary<T, LatSet, flagType>::UpdateBdCells() {
  Ids.clear();
  for (std::size_t id = 0; id < Geo.getVoxelsNum(); ++id) {
    if (util::isFlag(Field[id], BdCellFlag)) Ids.push_back(id);
  }
}

// --------------------------------------------------------------------------------
// FixedBoundary for block structure
// --------------------------------------------------------------------------------

template <typename T, typename LatSet, typename flagType>
BlockFixedBoundary<T, LatSet, flagType>::BlockFixedBoundary(
  BlockLattice<T, LatSet>& lat, const GenericArray<flagType>& f, std::uint8_t cellflag,
  std::uint8_t voidflag)
    : Lat(lat), Field(f), BdCellFlag(cellflag), voidFlag(voidflag) {
  Setup();
}

template <typename T, typename LatSet, typename flagType>
void BlockFixedBoundary<T, LatSet, flagType>::addtoBd(std::size_t id) {
  BdCells.emplace_back(id, LatSet::q);
  // get reference to the last element
  FixedBdCell& fixedbdcell = BdCells.back();
  // get neighbor
  // Attention: if voidFlag is 0, DO NOT use util::isFlag
  for (int k = 1; k < LatSet::q; ++k) {
    // if (util::isFlag(Field[Lat.getNbrId(id, k)], voidFlag) &&
    //     !util::isFlag(Field[Lat.getNbrId(id, LatSet::opp[k])], voidFlag)) {

    // using util::isFlag requires voidFlag to be non-zero
    if (util::isFlag(Field[Lat.getNbrId(id, k)], voidFlag)) {
      fixedbdcell.outflows.push_back(LatSet::opp[k]);
    }
  }
}

template <typename T, typename LatSet, typename flagType>
void BlockFixedBoundary<T, LatSet, flagType>::Setup() {
  std::size_t reserveSize;
  if constexpr (LatSet::d == 2) {
    reserveSize = (Lat.getGeo().getNx() + Lat.getGeo().getNy()) * 2;
  } else if constexpr (LatSet::d == 3) {
    reserveSize = (Lat.getGeo().getNx() * Lat.getGeo().getNy() +
                   Lat.getGeo().getNx() * Lat.getGeo().getNz() +
                   Lat.getGeo().getNy() * Lat.getGeo().getNz()) *
                  2;
  }
  BdCells.reserve(reserveSize);
  // add inner cells
  // in overlapped region, access neighbor cells may s
  if constexpr (LatSet::d == 2) {
    for (int iy = Lat.getGeo().getOverlap();
         iy < Lat.getGeo().getNy() - Lat.getGeo().getOverlap(); ++iy) {
      for (int ix = Lat.getGeo().getOverlap();
           ix < Lat.getGeo().getNx() - Lat.getGeo().getOverlap(); ++ix) {
        std::size_t id = ix + iy * Lat.getGeo().getNx();
        if (util::isFlag(Field[id], BdCellFlag)) addtoBd(id);
      }
    }
  } else if constexpr (LatSet::d == 3) {
    for (int iz = Lat.getGeo().getOverlap();
         iz < Lat.getGeo().getNz() - Lat.getGeo().getOverlap(); ++iz) {
      for (int iy = Lat.getGeo().getOverlap();
           iy < Lat.getGeo().getNy() - Lat.getGeo().getOverlap(); ++iy) {
        for (int ix = Lat.getGeo().getOverlap();
             ix < Lat.getGeo().getNx() - Lat.getGeo().getOverlap(); ++ix) {
          std::size_t id = ix + iy * Lat.getGeo().getProjection()[1] +
                           iz * Lat.getGeo().getProjection()[2];
          if (util::isFlag(Field[id], BdCellFlag)) addtoBd(id);
        }
      }
    }
  }
  // for (std::size_t id = 0; id < Lat.getGeo().getN(); ++id) {
  //   if (util::isFlag(Field[id], BdCellFlag)) addtoBd(id);
  // }
  // shrink capacity to actual size
  BdCells.shrink_to_fit();
}

// --------------------------------------------------------------------------------
// MovingBoundary for block structure
// --------------------------------------------------------------------------------

template <typename T, typename LatSet, typename flagType>
BlockMovingBoundary<T, LatSet, flagType>::BlockMovingBoundary(
  BlockLattice<T, LatSet>& lat, std::vector<std::size_t>& ids, GenericArray<flagType>& f,
  std::uint8_t voidflag, std::uint8_t cellflag)
    : Lat(lat), Ids(ids), Field(f), BdCellFlag(cellflag), voidFlag(voidflag) {}

template <typename T, typename LatSet, typename flagType>
void BlockMovingBoundary<T, LatSet, flagType>::UpdateBdCells() {
  Ids.clear();
  for (std::size_t id = 0; id < Lat.getGeo().getN(); ++id) {
    if (util::isFlag(Field[id], BdCellFlag)) Ids.push_back(id);
  }
}
