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
#include "lbm/equilibrium.ur.h"
#include "lbm/moment.ur.h"

template <typename T, typename LatSet>
struct BounceBackLikeMethod {
  static inline void normal_bounceback(PopCell<T, LatSet> &cell, int k) {
    cell[k] = cell.getPrevious(latset::opp<LatSet>(k));
  }
  static inline void anti_bounceback_simplified(PopCell<T, LatSet> &cell, int k) {
    cell[k] = 2 * cell.getRho() * latset::w<LatSet>(k) - cell.getPrevious(latset::opp<LatSet>(k));
  }
  // however the velocity here is NOT Wall velocity but fluid velocity relative to wall
  // so you don't need to specify the wall velocity when calling this function
  static inline void movingwall_bounceback(PopCell<T, LatSet> &cell, int k) {
    cell[k] = cell.getPrevious(latset::opp<LatSet>(k)) + 2 * LatSet::InvCs2 * latset::w<LatSet>(k) *
                                                   cell.getRho() *
                                                   (cell.getVelocity() * latset::c<LatSet>(k));
  }
  static inline void anti_bounceback_O1(PopCell<T, LatSet> &cell, int k) {
    cell[k] = 2 * Equilibrium<T, LatSet>::Order1(k, cell.getVelocity(), cell.getRho()) -
              cell.getPrevious(latset::opp<LatSet>(k));
  }
  static inline void anti_bounceback_O2(PopCell<T, LatSet> &cell, int k) {
    cell[k] = 2 * Equilibrium<T, LatSet>::Order2(k, cell.getVelocity(), cell.getRho(),
                                                 cell.getVelocity().getnorm2()) -
              cell.getPrevious(latset::opp<LatSet>(k));
  }
  static inline void anti_bounceback_pressure(PopCell<T, LatSet> &cell, int k) {
    // get the interpolated velocity
    const Vector<T, LatSet::d> uwall =
      cell.getVelocity() +
      T(0.5) * (cell.getVelocity() - cell.getNeighbor(latset::opp<LatSet>(k)).getVelocity());
    cell[k] = 2 * cell.getRho() * latset::w<LatSet>(k) *
                (T(1) + std::pow((uwall * latset::c<LatSet>(k)), 2) * T(0.5) * LatSet::InvCs4 -
                 uwall.getnorm2() * T(0.5) * LatSet::InvCs2) -
              cell.getPrevious(latset::opp<LatSet>(k));
  }
};

// Fixed BounceBackLike Boundary
template <typename T, typename LatSet, void (*BBLikemethod)(PopCell<T, LatSet> &, int),
          typename flagType = std::uint8_t>
class BBLikeFixedBoundary final : public FixedBoundary<T, LatSet, flagType> {
 private:
  std::string _name;

 public:
  BBLikeFixedBoundary(std::string name, PopLattice<T, LatSet> &lat, std::uint8_t cellflag,
                      std::uint8_t voidflag = std::uint8_t(1))
      : FixedBoundary<T, LatSet, flagType>(lat, cellflag, voidflag), _name(name) {}

  void Apply() override {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
    for (const auto &bdcell : this->BdCells) {
      PopCell<T, LatSet> cell(bdcell.Id, this->Lat);
      for (int k : bdcell.outflows) {
        BBLikemethod(cell, k);
      }
    }
  }
  void getinfo() override {
    std::cout << std::setw(18) << std::left << _name << std::setw(10) << std::left
              << this->BdCells.size() << std::endl;
  }
  void UpdateRho() override {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
    for (const auto &bdcell : this->BdCells) {
      BasicPopCell<T, LatSet> cell(bdcell.Id, this->Lat);
      moment::Rho<T, LatSet>::apply(cell, this->Lat.getRho(bdcell.Id));
    }
  }
  void UpdateU() override {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
    for (const auto &bdcell : this->BdCells) {
      BasicPopCell<T, LatSet> cell(bdcell.Id, this->Lat);
      moment::Velocity<T, LatSet>::apply(cell, this->Lat.getVelocity(bdcell.Id));
    }
  }
};

// Moving BounceBackLike Boundary
template <typename T, typename LatSet, void (*BBLikemethod)(PopCell<T, LatSet> &, int),
          typename flagType = std::uint8_t>
class BBLikeMovingBoundary final : public MovingBoundary<T, LatSet, flagType> {
 private:
  std::string _name;

 public:
  BBLikeMovingBoundary(std::string name, PopLattice<T, LatSet> &lat,
                       std::vector<std::size_t> &ids, std::uint8_t voidflag,
                       std::uint8_t cellflag = std::uint8_t(0))
      : MovingBoundary<T, LatSet, flagType>(lat, ids, voidflag, cellflag), _name(name) {}

  void Apply() override {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
    for (std::size_t id : this->Ids) {
      PopCell<T, LatSet> cell(id, this->Lat);
      for (unsigned int k = 1; k < LatSet::q; ++k) {
        if (util::isFlag(this->Field[this->Lat.getNbrId(id, k)], this->voidFlag)) {
          BBLikemethod(cell, latset::opp<LatSet>(k));
        }
      }
    }
  }
  void getinfo() override {
    std::cout << std::setw(18) << std::left << _name << std::setw(10) << std::left
              << this->Ids.size() << std::endl;
  }
  void UpdateRho() override {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
    for (std::size_t id : this->Ids) {
      BasicPopCell<T, LatSet> cell(id, this->Lat);
      moment::Rho<T, LatSet>::apply(cell, this->Lat.getRho(id));
    }
  }
  void UpdateU() override {
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
#endif
    for (std::size_t id : this->Ids) {
      BasicPopCell<T, LatSet> cell(id, this->Lat);
      moment::Velocity<T, LatSet>::apply(cell, this->Lat.getVelocity(id));
    }
  }
};


// --------------------------------------------------------------------------------------
// -----------------------------------BlockBoundary--------------------------------------
// --------------------------------------------------------------------------------------

namespace bounceback {

template <typename CELL>
struct normal {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = cell.getPrevious(latset::opp<LatSet>(k));
  }

  static inline void apply(CyclicArray<T>& arr, const CyclicArray<T>& arrk, std::size_t id) {
    arr[id] = arrk.getPrevious(id);
  }
};

template <typename CELL>
struct anti_simple {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = 2 * cell.template get<GenericRho>() * latset::w<LatSet>(k) -
              cell.getPrevious(latset::opp<LatSet>(k));
  }
};
template <typename CELL>
struct anti_O1 {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = 2 * equilibrium::FirstOrder<CELL>::get(
                    k, cell.template get<VELOCITY<T, LatSet::d>>(),
                    cell.template get<GenericRho>()) -
              cell.getPrevious(latset::opp<LatSet>(k));
  }
};
template <typename CELL>
struct anti_O2 {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = 2 * equilibrium::SecondOrder<CELL>::get(
                    k, cell.template get<VELOCITY<T, LatSet::d>>(),
                    cell.template get<GenericRho>(),
                    cell.template get<VELOCITY<T, LatSet::d>>().getnorm2()) -
              cell.getPrevious(latset::opp<LatSet>(k));
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
                cell.getNeighbor(latset::opp<LatSet>(k)).template get<VELOCITY<T, LatSet::d>>());
    cell[k] = 2 * cell.template get<GenericRho>() * latset::w<LatSet>(k) *
                (T{1} + std::pow((uwall * latset::c<LatSet>(k)), 2) * T{0.5} * LatSet::InvCs4 -
                 uwall.getnorm2() * T(0.5) * LatSet::InvCs2) -
              cell.getPrevious(latset::opp<LatSet>(k));
  }
};

template <typename CELL>
struct movingwall {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = cell.getPrevious(latset::opp<LatSet>(k)) +
              2 * LatSet::InvCs2 * latset::w<LatSet>(k) * cell.template get<GenericRho>() *
                (cell.template get<VELOCITY<T, LatSet::d>>() * latset::c<LatSet>(k));
  }
  static inline void apply(CELL &cell, unsigned int k, const Vector<T, LatSet::d> &wall_velocity) {
    cell[k] = cell.getPrevious(latset::opp<LatSet>(k)) - 2 * LatSet::InvCs2 * latset::w<LatSet>(k) *
                                                   cell.template get<GenericRho>() *
                                                   (wall_velocity * latset::c<LatSet>(k));
  }
};

}  // namespace bounceback
