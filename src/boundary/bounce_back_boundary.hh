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

// bounce_back_boundary.hh

#pragma once

#include "boundary/bounce_back_boundary.h"


// --------------------------------------------------------------------------
// Fixed BounceBackLike Boundary
// --------------------------------------------------------------------------

template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType>
void BBLikeFixedBoundary<T, LatSet, BBLikemethod, flagType>::Apply() {
#pragma omp parallel for num_threads(Thread_Num)
  for (const auto &bdcell : this->BdCells) {
    Cell<T, LatSet> cell(bdcell.Id, this->Lat);
    for (int k : bdcell.outflows) {
      BBLikemethod(cell, k);
    }
  }
}
template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType>
void BBLikeFixedBoundary<T, LatSet, BBLikemethod, flagType>::getinfo() {
  std::cout << std::setw(18) << std::left << _name << std::setw(10) << std::left
            << this->BdCells.size() << std::endl;
}

template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType>
void BBLikeFixedBoundary<T, LatSet, BBLikemethod, flagType>::UpdateRho() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (const auto &bdcell : this->BdCells) {
    BasicCell<T, LatSet> cell(bdcell.Id, this->Lat);
    moment::Rho<T, LatSet>::apply(cell, this->Lat.getRho(bdcell.Id));
  }
}

template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType>
void BBLikeFixedBoundary<T, LatSet, BBLikemethod, flagType>::UpdateU() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (const auto &bdcell : this->BdCells) {
    BasicCell<T, LatSet> cell(bdcell.Id, this->Lat);
    moment::Velocity<T, LatSet>::apply(cell, this->Lat.getVelocity(bdcell.Id));
  }
}


// --------------------------------------------------------------------------
// BBLikeMovingBoundary
// --------------------------------------------------------------------------

template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType>
BBLikeMovingBoundary<T, LatSet, BBLikemethod, flagType>::BBLikeMovingBoundary(
  std::string name, BasicLattice<T, LatSet> &lat, std::vector<std::size_t> &ids,
  std::uint8_t voidflag, std::uint8_t cellflag)
    : MovingBoundary<T, LatSet, flagType>(lat, ids, voidflag, cellflag), _name(name) {}

template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType>
void BBLikeMovingBoundary<T, LatSet, BBLikemethod, flagType>::Apply() {
  for (std::size_t id : this->Ids) {
    Cell<T, LatSet> cell(id, this->Lat);
    for (int k = 1; k < LatSet::q; ++k) {
      if (util::isFlag(this->Field[this->Lat.getNbrId(id, k)], this->voidFlag)) {
        BBLikemethod(cell, LatSet::opp[k]);
      }
    }
  }
}

template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType>
void BBLikeMovingBoundary<T, LatSet, BBLikemethod, flagType>::getinfo() {
  std::cout << std::setw(18) << std::left << _name << std::setw(10) << std::left
            << this->Ids.size() << std::endl;
}

template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType>
void BBLikeMovingBoundary<T, LatSet, BBLikemethod, flagType>::UpdateRho() {
  for (std::size_t id : this->Ids) {
    BasicCell<T, LatSet> cell(id, this->Lat);
    moment::Rho<T, LatSet>::apply(cell, this->Lat.getRho(id));
  }
}

template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType>
void BBLikeMovingBoundary<T, LatSet, BBLikemethod, flagType>::UpdateU() {
  for (std::size_t id : this->Ids) {
    BasicCell<T, LatSet> cell(id, this->Lat);
    moment::Velocity<T, LatSet>::apply(cell, this->Lat.getVelocity(id));
  }
}


// --------------------------------------------------------------------------
// BBLikeFixedBlockBoundary
// --------------------------------------------------------------------------

template <typename T, typename LatSet, typename flagType>
BBLikeFixedBlockBoundary<T, LatSet, flagType>::BBLikeFixedBlockBoundary(
  BlockLattice<T, LatSet> &lat, const GenericArray<flagType> &f, std::uint8_t cellflag,
  std::uint8_t voidflag)
    : BlockFixedBoundary<T, LatSet, flagType>(lat, f, cellflag, voidflag) {}

template <typename T, typename LatSet, typename flagType>
template <void (*BBLikemethod)(BCell<T, LatSet> &, int)>
void BBLikeFixedBlockBoundary<T, LatSet, flagType>::Apply() {
  for (const auto &bdcell : this->BdCells) {
    BCell<T, LatSet> cell(bdcell.Id, this->Lat);
    for (int k : bdcell.outflows) {
      BBLikemethod(cell, k);
    }
  }
}

template <typename T, typename LatSet, typename flagType>
void BBLikeFixedBlockBoundary<T, LatSet, flagType>::UpdateRho() {
  for (const auto &bdcell : this->BdCells) {
    BasicCell<T, LatSet> cell(bdcell.Id, this->Lat);
    moment::Rho<T, LatSet>::apply(cell, this->Lat.getRho(bdcell.Id));
  }
}

template <typename T, typename LatSet, typename flagType>
void BBLikeFixedBlockBoundary<T, LatSet, flagType>::UpdateU() {
  for (const auto &bdcell : this->BdCells) {
    BasicCell<T, LatSet> cell(bdcell.Id, this->Lat);
    moment::Velocity<T, LatSet>::apply(cell, this->Lat.getVelocity(bdcell.Id));
  }
}


// --------------------------------------------------------------------------
// BBLikeFixedBlockBdManager
// --------------------------------------------------------------------------

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType>
BBLikeFixedBlockBdManager<T, LatSet, BBLikemethod, flagType>::BBLikeFixedBlockBdManager(
  std::string name, BlockLatticeManager<T, LatSet> &lat,
  BlockFieldManager<GenericField<GenericArray<flagType>, 1>, T, LatSet::d> &BlockFM,
  std::uint8_t cellflag, std::uint8_t voidflag)
    : _name(name), BdCellFlag(cellflag), voidFlag(voidflag), LatMan(lat),
      BlockFManager(BlockFM) {
  Init();
}

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType>
void BBLikeFixedBlockBdManager<T, LatSet, BBLikemethod, flagType>::Init() {
  BdBlocks.clear();
  // for each blocks in blocklat
  for (int i = 0; i < LatMan.getGeo().getBlockNum(); ++i) {
    BdBlocks.emplace_back(LatMan.getBlockLat(i),
                          BlockFManager.getBlockField(i).getField().getField(0),
                          BdCellFlag, voidFlag);
  }
}

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType>
void BBLikeFixedBlockBdManager<T, LatSet, BBLikemethod, flagType>::Apply(
  std::int64_t count) {
  std::uint8_t MaxLevel = LatMan.getMaxLevel();
#pragma omp parallel for num_threads(Thread_Num)
  for (auto &bdBlock : BdBlocks) {
    if (count % (static_cast<int>(pow(2, int(MaxLevel - bdBlock.getLat().getLevel())))) ==
        0)
      bdBlock.template Apply<BBLikemethod>();
  }
}

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType>
void BBLikeFixedBlockBdManager<T, LatSet, BBLikemethod, flagType>::UpdateRho(
  std::int64_t count) {
  std::uint8_t MaxLevel = LatMan.getMaxLevel();
#pragma omp parallel for num_threads(Thread_Num)
  for (auto &bdBlock : BdBlocks) {
    if (count % (static_cast<int>(pow(2, int(MaxLevel - bdBlock.getLat().getLevel())))) ==
        0)
      bdBlock.UpdateRho();
  }
}

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType>
void BBLikeFixedBlockBdManager<T, LatSet, BBLikemethod, flagType>::UpdateU(
  std::int64_t count) {
  std::uint8_t MaxLevel = LatMan.getMaxLevel();
#pragma omp parallel for num_threads(Thread_Num)
  for (auto &bdBlock : BdBlocks) {
    if (count % (static_cast<int>(pow(2, int(MaxLevel - bdBlock.getLat().getLevel())))) ==
        0)
      bdBlock.UpdateU();
  }
}


// --------------------------------------------------------------------------
// BBLikeMovingBlockBoundary
// --------------------------------------------------------------------------

template <typename T, typename LatSet, typename flagType>
BBLikeMovingBlockBoundary<T, LatSet, flagType>::BBLikeMovingBlockBoundary(
  BlockLattice<T, LatSet> &lat, std::vector<std::size_t> *ids, GenericArray<flagType> &f,
  std::uint8_t voidflag, std::uint8_t cellflag)
    : BlockMovingBoundary<T, LatSet, flagType>(lat, *ids, f, voidflag, cellflag) {}

template <typename T, typename LatSet, typename flagType>
template <void (*BBLikemethod)(BCell<T, LatSet> &, int)>
void BBLikeMovingBlockBoundary<T, LatSet, flagType>::Apply() {
  for (std::size_t id : this->Ids) {
    BCell<T, LatSet> cell(id, this->Lat);
    for (int k = 1; k < LatSet::q; ++k) {
      if (util::isFlag(this->Field[this->Lat.getNbrId(id, k)], this->voidFlag)) {
        BBLikemethod(cell, LatSet::opp[k]);
      }
    }
  }
}

template <typename T, typename LatSet, typename flagType>
void BBLikeMovingBlockBoundary<T, LatSet, flagType>::UpdateRho() {
  for (std::size_t id : this->Ids) {
    BasicCell<T, LatSet> cell(id, this->Lat);
    moment::Rho<T, LatSet>::apply(cell, this->Lat.getRho(id));
  }
}

template <typename T, typename LatSet, typename flagType>
void BBLikeMovingBlockBoundary<T, LatSet, flagType>::UpdateU() {
  for (std::size_t id : this->Ids) {
    BasicCell<T, LatSet> cell(id, this->Lat);
    moment::Velocity<T, LatSet>::apply(cell, this->Lat.getVelocity(id));
  }
}

// --------------------------------------------------------------------------
// BBLikeBlockMovingBd Manager
// --------------------------------------------------------------------------

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType>
BBLikeMovingBlockBdManager<T, LatSet, BBLikemethod, flagType>::BBLikeMovingBlockBdManager(
  std::string name, BlockLatticeManager<T, LatSet> &lat,
  std::vector<std::vector<std::size_t> *> &idss,
  BlockFieldManager<ScalerField<flagType>, T, LatSet::d> &BlockFM, std::uint8_t voidflag,
  std::uint8_t cellflag)
    : _name(name), BdCellFlag(cellflag), voidFlag(voidflag), LatMan(lat), IDss(idss),
      BlockFManager(BlockFM) {
  Init();
}

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType>
void BBLikeMovingBlockBdManager<T, LatSet, BBLikemethod, flagType>::Init() {
  BdBlocks.clear();
  for (int i = 0; i < LatMan.getBlockLats().size(); ++i) {
    BdBlocks.emplace_back(LatMan.getBlockLat(i), IDss[i],
                          BlockFManager.getBlockField(i).getField().getField(0), voidFlag,
                          BdCellFlag);
  }
}

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType>
void BBLikeMovingBlockBdManager<T, LatSet, BBLikemethod, flagType>::Apply(
  std::int64_t count) {
  std::uint8_t MaxLevel = LatMan.getMaxLevel();
#pragma omp parallel for num_threads(Thread_Num)
  for (auto &bdBlock : BdBlocks) {
    if (count % (static_cast<int>(pow(2, int(MaxLevel - bdBlock.getLat().getLevel())))) ==
        0)
      bdBlock.template Apply<BBLikemethod>();
  }
}

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType>
void BBLikeMovingBlockBdManager<T, LatSet, BBLikemethod, flagType>::UpdateRho(
  std::int64_t count) {
  std::uint8_t MaxLevel = LatMan.getMaxLevel();
#pragma omp parallel for num_threads(Thread_Num)
  for (auto &bdBlock : BdBlocks) {
    if (count % (static_cast<int>(pow(2, int(MaxLevel - bdBlock.getLat().getLevel())))) ==
        0)
      bdBlock.UpdateRho();
  }
}

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType>
void BBLikeMovingBlockBdManager<T, LatSet, BBLikemethod, flagType>::UpdateU(
  std::int64_t count) {
  std::uint8_t MaxLevel = LatMan.getMaxLevel();
#pragma omp parallel for num_threads(Thread_Num)
  for (auto &bdBlock : BdBlocks) {
    if (count % (static_cast<int>(pow(2, int(MaxLevel - bdBlock.getLat().getLevel())))) ==
        0)
      bdBlock.UpdateU();
  }
}