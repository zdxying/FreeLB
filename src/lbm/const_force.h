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

// const_force.h
#pragma once

#include "data_struct/block_lattice.h"

// const force for blocklattice

template <typename T, typename LatSet>
class BlockConstForce {
 private:
  BlockLattice<T, LatSet> &Lat;

  VectorFieldAOS<T, LatSet::d> &Velocity;

  Vector<T, LatSet::d> Force;
  // omega
  T Omega;
  // 1 - omega
  T _Omega;
  // 1 - omega/2
  T fOmega;

 public:
  BlockConstForce(BlockLattice<T, LatSet> &lat, VectorFieldAOS<T, LatSet::d> &velocity,
                  Vector<T, LatSet::d> force)
      : Lat(lat), Velocity(velocity), Force(force), Omega(lat.getOmega()),
        _Omega(lat.get_Omega()), fOmega(T(1) - lat.getOmega() * T(0.5)) {}
  std::uint8_t getLevel() const { return Lat.getLevel(); }
  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T),
            typename flagtype>
  void BGK_U(const GenericArray<flagtype> &flagarr, std::uint8_t flag) {
    for (std::size_t id = 0; id < Lat.getN(); ++id) {
      BCell<T, LatSet> cell(id, Lat);
      std::array<T, LatSet::q> Fi{};
      force::Force<T, LatSet>::ComputeForcePop(Fi, Velocity.get(id), Force);
      moment::Velocity<T, LatSet>::apply(cell, Velocity.get(id), Fi);
      collision::BGK<T, LatSet>::template applySource<GetFeq>(cell, Fi);
    }
  }
};

template <typename T, typename LatSet>
class ConstForceManager {
  std::vector<BlockConstForce<T, LatSet> > ConstForces;

  BlockLatticeManager<T, LatSet> &LatMan;

  BlockFieldManager<VectorFieldAOS<T, LatSet::d>, T, LatSet::d> &Velocity;

  Vector<T, LatSet::d> Force;

 public:
  ConstForceManager(
    BlockLatticeManager<T, LatSet> &latman,
    BlockFieldManager<VectorFieldAOS<T, LatSet::d>, T, LatSet::d> &velocity,
    Vector<T, LatSet::d> force)
      : LatMan(latman), Velocity(velocity), Force(force) {
    for (int i = 0; i < LatMan.getBlockLats().size(); ++i) {
      ConstForces.emplace_back(LatMan.getBlockLat(i),
                               Velocity.getBlockField(i).getField(), Force);
    }
  }

  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T),
            typename flagtype>
  void BGK_U(std::int64_t count, std::uint8_t flag,
             const BlockFieldManager<ScalerField<flagtype>, T, LatSet::d> &BFM) {
    std::uint8_t MaxLevel = LatMan.getMaxLevel();
#pragma omp parallel for num_threads(Thread_Num)
    for (int i = 0; i < ConstForces.size(); ++i) {
      if (count % (static_cast<int>(pow(2, int(MaxLevel - ConstForces[i].getLevel())))) ==
          0)
        ConstForces[i].template BGK_U<GetFeq>(BFM.getBlockField(i).getField().getField(0),
                                              flag);
    }
  }
};
