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

// bounceback like boundary condition

#pragma once

#include "boundary/basic_boundary.h"

template <typename T, typename LatSet>
struct BounceBackLikeMethod {
  static inline void normal_bounceback(Cell<T, LatSet> &cell, int k) {
    cell[k] = cell.getPrevious(LatSet::opp[k]);
  }
  static inline void anti_bounceback_simplified(Cell<T, LatSet> &cell, int k) {
    cell[k] = 2 * cell.getRho() * LatSet::w[k] - cell.getPrevious(LatSet::opp[k]);
  }
  // however the velocity here is NOT Wall velocity but fluid velocity relative to wall
  // so you don't need to specify the wall velocity when calling this function
  static inline void movingwall_bounceback(Cell<T, LatSet> &cell, int k) {
    cell[k] = cell.getPrevious(LatSet::opp[k]) + 2 * LatSet::InvCs2 * LatSet::w[k] *
                                                   cell.getRho() *
                                                   (cell.getVelocity() * LatSet::c[k]);
  }
  static inline void anti_bounceback_O1(Cell<T, LatSet> &cell, int k) {
    cell[k] = 2 * Equilibrium<T, LatSet>::Order1(k, cell.getVelocity(), cell.getRho()) -
              cell.getPrevious(LatSet::opp[k]);
  }
  static inline void anti_bounceback_O2(Cell<T, LatSet> &cell, int k) {
    cell[k] = 2 * Equilibrium<T, LatSet>::Order2(k, cell.getVelocity(), cell.getRho(),
                                                 cell.getVelocity().getnorm2()) -
              cell.getPrevious(LatSet::opp[k]);
  }
  static inline void anti_bounceback_pressure(Cell<T, LatSet> &cell, int k) {
    // get the interpolated velocity
    const Vector<T, LatSet::d> uwall =
      cell.getVelocity() +
      T(0.5) * (cell.getVelocity() - cell.getNeighbor(LatSet::opp[k]).getVelocity());
    cell[k] = 2 * cell.getRho() * LatSet::w[k] *
                (T(1) + pow((uwall * LatSet::c[k]), 2) * T(0.5) * LatSet::InvCs4 -
                 uwall.getnorm2() * T(0.5) * LatSet::InvCs2) -
              cell.getPrevious(LatSet::opp[k]);
  }
};

// Fixed BounceBackLike Boundary
template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType = std::uint8_t>
class BBLikeFixedBoundary final : public FixedBoundary<T, LatSet, flagType> {
 private:
  std::string _name;

 public:
  BBLikeFixedBoundary(std::string name, BasicLattice<T, LatSet> &lat,
                      std::uint8_t cellflag, std::uint8_t voidflag = std::uint8_t(1))
      : FixedBoundary<T, LatSet, flagType>(lat, cellflag, voidflag), _name(name) {}

  void Apply() override;
  void getinfo() override;
  void UpdateRho() override;
  void UpdateU() override;
};

// Moving BounceBackLike Boundary
template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType = std::uint8_t>
class BBLikeMovingBoundary final : public MovingBoundary<T, LatSet, flagType> {
 private:
  std::string _name;

 public:
  BBLikeMovingBoundary(std::string name, BasicLattice<T, LatSet> &lat,
                       std::vector<std::size_t> &ids, std::uint8_t voidflag,
                       std::uint8_t cellflag = std::uint8_t(0));

  void Apply() override;
  void getinfo() override;
  void UpdateRho() override;
  void UpdateU() override;
};

// --------------------------------------------------------------------------------------
// -----------------------------------BlockBoundary--------------------------------------
// --------------------------------------------------------------------------------------

namespace bounceback {
template <typename CELL>
struct normal {
  using LatSet = typename CELL::LatticeSet;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = cell.getPrevious(LatSet::opp[k]);
  }
};

template <typename CELL>
struct anti_simple {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] =
      2 * cell.template get<GenericRho>() * LatSet::w[k] - cell.getPrevious(LatSet::opp[k]);
  }
};
template <typename CELL>
struct anti_O1 {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] =
      2 * equilibrium::FirstOrder<CELL>::get(
            k, cell.template get<VELOCITY<T, LatSet::d>>(), cell.template get<GenericRho>()) -
      cell.getPrevious(LatSet::opp[k]);
  }
};
template <typename CELL>
struct anti_O2 {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] =
      2 * equilibrium::SecondOrder<CELL>::get(
            k, cell.template get<VELOCITY<T, LatSet::d>>(), cell.template get<GenericRho>(),
            cell.template get<VELOCITY<T, LatSet::d>>().getnorm2()) -
      cell.getPrevious(LatSet::opp[k]);
  }
};
template <typename CELL>
struct anti_pressure {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    // get the interpolated velocity
    const Vector<T, LatSet::d> uwall =
      cell.template get<VELOCITY<T, LatSet::d>>() +
      T{0.5} * (cell.template get<VELOCITY<T, LatSet::d>>() -
                cell.getNeighbor(LatSet::opp[k]).template get<VELOCITY<T, LatSet::d>>());
    cell[k] = 2 * cell.template get<GenericRho>() * LatSet::w[k] *
                (T{1} + std::pow((uwall * LatSet::c[k]), 2) * T{0.5} * LatSet::InvCs4 -
                 uwall.getnorm2() * T(0.5) * LatSet::InvCs2) -
              cell.getPrevious(LatSet::opp[k]);
  }
};

template <typename CELL>
struct movingwall {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = cell.getPrevious(LatSet::opp[k]) +
              2 * LatSet::InvCs2 * LatSet::w[k] * cell.template get<GenericRho>() *
                (cell.template get<VELOCITY<T, LatSet::d>>() * LatSet::c[k]);
  }
  static inline void apply(CELL &cell, unsigned int k,
                           const Vector<T, LatSet::d> &wall_velocity) {
    cell[k] = cell.getPrevious(LatSet::opp[k]) - 2 * LatSet::InvCs2 * LatSet::w[k] *
                                                   cell.template get<GenericRho>() *
                                                   (wall_velocity * LatSet::c[k]);
  }
};

}  // namespace bounceback

// --------------------------------------------------------------------------
// BBLike Fixed BlockBoundary
// --------------------------------------------------------------------------
template <typename BLOCKLATTICE, typename ArrayType>
class BBLikeFixedBlockBoundary : public BlockFixedBoundary<BLOCKLATTICE, ArrayType> {
 public:
  using CELL = typename BLOCKLATTICE::CellType;

  BBLikeFixedBlockBoundary(BLOCKLATTICE &lat, const ArrayType &f, std::uint8_t cellflag,
                           std::uint8_t voidflag);
  template <typename CELLDYNAMICS>
  void ApplyCellDynamics();
};

template <typename CELLDYNAMICS, typename BLOCKLATTICEMANAGER, typename BLOCKFIELDMANAGER>
class BBLikeFixedBlockBdManager final : public AbstractBlockBoundary {
 public:
  using BLOCKLATTICE = typename BLOCKLATTICEMANAGER::BLOCKLATTICE;
  using ArrayType = typename BLOCKFIELDMANAGER::array_type;

 private:
  std::string _name;
  std::vector<BBLikeFixedBlockBoundary<BLOCKLATTICE, ArrayType>> BdBlocks;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;

  BLOCKLATTICEMANAGER &LatMan;

  BLOCKFIELDMANAGER &BlockFManager;

 public:
  BBLikeFixedBlockBdManager(std::string name, BLOCKLATTICEMANAGER &lat,
                            BLOCKFIELDMANAGER &BlockFM, std::uint8_t cellflag,
                            std::uint8_t voidflag = std::uint8_t(1));

  void Init();
  void Apply(std::int64_t count) override;
};


// --------------------------------------------------------------------------
// BBLike Moving BlockBoundary
// --------------------------------------------------------------------------
template <typename BLOCKLATTICE, typename ArrayType>
class BBLikeMovingBlockBoundary : public BlockMovingBoundary<BLOCKLATTICE, ArrayType> {
 public:
  using CELL = typename BLOCKLATTICE::CellType;
  using LatSet = typename CELL::LatticeSet;

  BBLikeMovingBlockBoundary(BLOCKLATTICE &lat, std::vector<std::size_t> *ids,
                            ArrayType &f, std::uint8_t voidflag, std::uint8_t cellflag);

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics();
};

template <typename CELLDYNAMICS, typename BLOCKLATTICEMANAGER, typename BLOCKFIELDMANAGER>
class BBLikeMovingBlockBdManager final : public AbstractBlockBoundary {
 public:
  using BLOCKLATTICE = typename BLOCKLATTICEMANAGER::BLOCKLATTICE;
  using ArrayType = typename BLOCKFIELDMANAGER::array_type;

  BBLikeMovingBlockBdManager(std::string name, BLOCKLATTICEMANAGER &lat,
                             std::vector<std::vector<std::size_t> *> &idss,
                             BLOCKFIELDMANAGER &BlockFM, std::uint8_t voidflag,
                             std::uint8_t cellflag = std::uint8_t(0));

  void Init();
  void Apply(std::int64_t count) override;

 private:
  std::string _name;
  std::vector<BBLikeMovingBlockBoundary<BLOCKLATTICE, ArrayType>> BdBlocks;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;
  BLOCKLATTICEMANAGER &LatMan;
  // ids
  std::vector<std::vector<std::size_t> *> &IDss;

  BLOCKFIELDMANAGER &BlockFManager;
};
