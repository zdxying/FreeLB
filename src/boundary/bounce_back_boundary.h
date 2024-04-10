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

// bounceback like boundary condition

#pragma once

#include "boundary/basic_boundary.h"

template <typename T, typename LatSet>
struct BounceBackLikeMethod {
  static inline void normal_bounceback(Cell<T, LatSet> &cell, int k) {
    cell[k] = cell.getPrevious(LatSet::opp[k]);
  }
  static inline void normal_bounceback(BCell<T, LatSet> &cell, int k) {
    cell[k] = cell.getPrevious(LatSet::opp[k]);
  }

  static inline void anti_bounceback_simplified(Cell<T, LatSet> &cell, int k) {
    cell[k] = 2 * cell.getRho() * LatSet::w[k] - cell.getPrevious(LatSet::opp[k]);
  }
  static inline void anti_bounceback_simplified(BCell<T, LatSet> &cell, int k) {
    cell[k] = 2 * cell.getRho() * LatSet::w[k] - cell.getPrevious(LatSet::opp[k]);
  }

  static inline void movingwall_bounceback(Cell<T, LatSet> &cell, int k) {
    cell[k] = cell.getPrevious(LatSet::opp[k]) + 2 * LatSet::InvCs2 * LatSet::w[k] * cell.getRho() *
                                                   (cell.getVelocity() * LatSet::c[k]);
  }
  static inline void movingwall_bounceback(BCell<T, LatSet> &cell, int k) {
    cell[k] = cell.getPrevious(LatSet::opp[k]) + 2 * LatSet::InvCs2 * LatSet::w[k] * cell.getRho() *
                                                   (cell.getVelocity() * LatSet::c[k]);
  }

  static inline void anti_bounceback_O1(Cell<T, LatSet> &cell, int k) {
    cell[k] = 2 * Equilibrium<T, LatSet>::Order1(k, cell.getVelocity(), cell.getRho()) -
              cell.getPrevious(LatSet::opp[k]);
  }
  static inline void anti_bounceback_O1(BCell<T, LatSet> &cell, int k) {
    cell[k] = 2 * Equilibrium<T, LatSet>::Order1(k, cell.getVelocity(), cell.getRho()) -
              cell.getPrevious(LatSet::opp[k]);
  }

  static inline void anti_bounceback_O2(Cell<T, LatSet> &cell, int k) {
    cell[k] = 2 * Equilibrium<T, LatSet>::Order2(k, cell.getVelocity(), cell.getRho(),
                                                 cell.getVelocity().getnorm2()) -
              cell.getPrevious(LatSet::opp[k]);
  }
  static inline void anti_bounceback_O2(BCell<T, LatSet> &cell, int k) {
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
  static inline void anti_bounceback_pressure(BCell<T, LatSet> &cell, int k) {
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
  BBLikeFixedBoundary(std::string name, BasicLattice<T, LatSet> &lat, std::uint8_t cellflag,
                      std::uint8_t voidflag = std::uint8_t(0))
      : FixedBoundary<T, LatSet, flagType>(lat, cellflag, voidflag), _name(name) {}

  void Apply() override;
  void getinfo() override;
  void UpdateRho() override;
  void UpdateU() override;
};

// Block Fixed BounceBackLike Boundary
template <typename T, typename LatSet, typename flagType>
class BBLikeBlockFixedBoundary : public BlockFixedBoundary<T, LatSet, flagType> {
 public:
  BBLikeBlockFixedBoundary(BlockLattice<T, LatSet> &lat, std::uint8_t cellflag,
                           std::uint8_t voidflag);
  BBLikeBlockFixedBoundary(BlockLattice<T, LatSet> &lat, const GenericArray<flagType> &f,
                           std::uint8_t cellflag, std::uint8_t voidflag);

  template <void (*BBLikemethod)(BCell<T, LatSet> &, int)>
  void Apply();
  void UpdateRho();
  void UpdateU();
};

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType = std::uint8_t>
class BBLikeBlockFixedBdManager final : public AbstractBlockBoundary {
 private:
  std::string _name;
  std::vector<BBLikeBlockFixedBoundary<T, LatSet, flagType>> BdBlocks;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;
  BlockLatticeManager<T, LatSet> &LatMan;

 public:
  BBLikeBlockFixedBdManager(std::string name, BlockLatticeManager<T, LatSet> &lat,
                            std::uint8_t cellflag, std::uint8_t voidflag = std::uint8_t(1));
  BBLikeBlockFixedBdManager(std::string name, BlockLatticeManager<T, LatSet> &lat,
                            std::vector<GenericArray<flagType> *> fs, std::uint8_t cellflag,
                            std::uint8_t voidflag = std::uint8_t(1));
  // this version is for mpi, may be removed
  BBLikeBlockFixedBdManager(std::string name, std::vector<BlockLattice<T, LatSet> *> lats,
                            std::uint8_t cellflag, std::uint8_t voidflag = std::uint8_t(1));

  void Setup(BlockLatticeManager<T, LatSet> &lat);
  void Setup(BlockLatticeManager<T, LatSet> &lat, std::vector<GenericArray<flagType> *> fs);
  void Setup(std::vector<BlockLattice<T, LatSet> *> lats);
  void Apply(std::int64_t count) override;
  void UpdateRho(std::int64_t count) override;
  void UpdateU(std::int64_t count) override;
};

// Moving BounceBackLike Boundary
template <typename T, typename LatSet, void (*BBLikemethod)(Cell<T, LatSet> &, int),
          typename flagType = std::uint8_t>
class MovingBounceBackLikeBoundary final : public MovingBoundary<T, LatSet, flagType> {
 private:
  std::string _name;

 public:
  MovingBounceBackLikeBoundary(std::string name, BasicLattice<T, LatSet> &lat,
                               std::vector<std::size_t> &ids, std::uint8_t voidflag,
                               std::uint8_t cellflag = std::uint8_t(0));

  void Apply() override;
  void getinfo() override;
  void UpdateRho() override;
  void UpdateU() override;
};

// MovingBounceBackLikeBoundaryBlock
template <typename T, typename LatSet, typename flagType>
class BBLikeBlockMovingBoundary : public BlockMovingBoundary<T, LatSet, flagType> {
 public:
  BBLikeBlockMovingBoundary(BlockLattice<T, LatSet> &lat, std::vector<std::size_t> *ids,
                            std::uint8_t voidflag, std::uint8_t cellflag);
  BBLikeBlockMovingBoundary(BlockLattice<T, LatSet> &lat, std::vector<std::size_t> *ids,
                            const GenericArray<flagType> *f, std::uint8_t voidflag,
                            std::uint8_t cellflag);

  template <void (*BBLikemethod)(BCell<T, LatSet> &, int)>
  void Apply();
  void UpdateRho();
  void UpdateU();
};

template <typename T, typename LatSet, void (*BBLikemethod)(BCell<T, LatSet> &, int),
          typename flagType = std::uint8_t>
class BBLikeBlockMovingBdManager final : public AbstractBlockBoundary {
 private:
  std::string _name;
  std::vector<BBLikeBlockMovingBoundary<T, LatSet, flagType>> BdBlocks;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;
  BlockLatticeManager<T, LatSet> &LatMan;

 public:
  BBLikeBlockMovingBdManager(std::string name, BlockLatticeManager<T, LatSet> &lat,
                             std::vector<std::vector<std::size_t> *> idss, std::uint8_t voidflag,
                             std::uint8_t cellflag = std::uint8_t(0));
  BBLikeBlockMovingBdManager(std::string name, BlockLatticeManager<T, LatSet> &lat,
                             std::vector<std::vector<std::size_t> *> idss,
                             std::vector<GenericArray<flagType> *> fs, std::uint8_t voidflag,
                             std::uint8_t cellflag = std::uint8_t(0));
  // this version is for mpi, may be removed
  // BBLikeBlockMovingBdManager(std::string name, std::vector<BlockLattice<T, LatSet> *> lats,
  // std::vector<std::vector<std::size_t>*> idss,
  //                           std::uint8_t cellflag, std::uint8_t voidflag = std::uint8_t(0));

  void Setup(std::vector<std::vector<std::size_t> *> idss);
  void Setup(std::vector<std::vector<std::size_t> *> idss,
             std::vector<GenericArray<flagType> *> fs);
  // void Setup(std::vector<BlockLattice<T, LatSet> *> lats);
  void Apply(std::int64_t count) override;
  void UpdateRho(std::int64_t count) override;
  void UpdateU(std::int64_t count) override;
};
