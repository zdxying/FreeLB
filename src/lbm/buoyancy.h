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

// buoyancy.h
#pragma once

#include "data_struct/block_lattice.h"

template <typename T, typename LatSet>
class BlockBuoyancy {
 private:
  // NS lattice
  BlockLattice<T, LatSet> &NSLat;
  // source lattice, e.g., thermal and solute lattice
  std::vector<BlockRhoLattice<T> *> Source;

  VectorFieldAOS<T, LatSet::d> &Velocity;

  ScalerField<T> Force;
  // omega
  T Omega;
  // 1 - omega
  T _Omega;
  // 1 - omega/2
  T fOmega;

 public:
  BlockBuoyancy(BlockLattice<T, LatSet> &lat, VectorFieldAOS<T, LatSet::d> &velocity)
      : NSLat(lat), Velocity(velocity), Omega(lat.getOmega()), _Omega(lat.get_Omega()),
        fOmega(T(1) - lat.getOmega() * T(0.5)), Force(lat.getN(), T(0)) {}
  template <typename... Args>
  void AddSource(BlockRhoLattice<T> *lat, Args... args) {
    Source.push_back(lat);
    AddSource(args...);
  }
  void AddSource(BlockRhoLattice<T> *lat) { Source.push_back(lat); }

  std::uint8_t getLevel() const { return NSLat.getLevel(); }

  template <typename flagtype>
  void GetBuoyancy(const GenericArray<flagtype> &flagarr, std::uint8_t flag) {
    // reset force
    Force.getField(0).Init(T(0));
    // add to buoyancy
    for (BlockRhoLattice<T> *lat : Source) {
      T latRhoInit = lat->getLatRhoInit();
      T latgbeta = lat->getLatgBeta(NSLat.getLevel());
      for (std::size_t id = 0; id < NSLat.getN(); ++id) {
        if (util::isFlag(flagarr[id], flag))
          Force.get(id) += (lat->getRho(id) - latRhoInit) * latgbeta;
      }
    }
  }
  Vector<T, LatSet::d> getForce(std::size_t id) const {
    if constexpr (LatSet::d == 3) {
      return Vector<T, LatSet::d>{T(0), T(0), Force.get(id)};
    } else if constexpr (LatSet::d == 2) {
      return Vector<T, LatSet::d>{T(0), Force.get(id)};
    }
  }
  // BGK with FORCE term
  // update force and u at the same time
  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T),
            typename flagtype>
  void BGK_U(const GenericArray<flagtype> &flagarr, std::uint8_t flag) {
    for (std::size_t id = 0; id < NSLat.getN(); ++id) {
      if (util::isFlag(flagarr[id], flag)) {
        BCell<T, LatSet> cell(id, NSLat);
        std::array<T, LatSet::q> Fi{};
        // force::Force<T, LatSet>::ComputeForcePop(Fi, Velocity.get(id), getForce(id));
        force::Force<T, LatSet>::ComputeForcePop(Fi, Velocity.get(id), Force.get(id));
        moment::Velocity<T, LatSet>::apply(cell, Velocity.get(id), Fi);
        collision::BGK<T, LatSet>::template applySource<GetFeq>(cell, Fi);
      }
    }
  }
  // BGK with FORCE term
  // do not update u
  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T),
            typename flagtype>
  void BGK(const GenericArray<flagtype> &flagarr, std::uint8_t flag) {
    for (std::size_t id = 0; id < NSLat.getN(); ++id) {
      if (util::isFlag(flagarr[id], flag)) {
        BCell<T, LatSet> cell(id, NSLat);
        std::array<T, LatSet::q> Fi{};
        // force::Force<T, LatSet>::ComputeForcePop(Fi, Velocity.get(id), getForce(id));
        force::Force<T, LatSet>::ComputeForcePop(Fi, Velocity.get(id), Force.get(id));
        collision::BGK<T, LatSet>::template applySource<GetFeq>(cell, Fi);
      }
    }
  }
};

template <typename T, typename LatSet>
class BlockBuoyancyManager {
 private:
  std::vector<BlockBuoyancy<T, LatSet> > _Buoyancys;

  BlockLatticeManager<T, LatSet> &NSLatMan;

  BlockFieldManager<VectorFieldAOS<T, LatSet::d>, T, LatSet::d> &Velocity;

 public:
  BlockBuoyancyManager(
    BlockLatticeManager<T, LatSet> &nsLatMan,
    BlockFieldManager<VectorFieldAOS<T, LatSet::d>, T, LatSet::d> &velocity)
      : NSLatMan(nsLatMan), Velocity(velocity) {
    for (int i = 0; i < NSLatMan.getBlockLats().size(); ++i) {
      _Buoyancys.emplace_back(NSLatMan.getBlockLat(i),
                              velocity.getBlockField(i).getField());
    }
  }

  void Init() {
    _Buoyancys.clear();
    for (int i = 0; i < NSLatMan.getBlockLats().size(); ++i) {
      _Buoyancys.emplace_back(NSLatMan.getBlockLat(i),
                              Velocity.getBlockField(i).getField());
    }
  }

  BlockBuoyancy<T, LatSet> &getBuoyancy(int i) { return _Buoyancys[i]; }
  const BlockBuoyancy<T, LatSet> &getBuoyancy(int i) const { return _Buoyancys[i]; }

  template <typename LatSet1>
  void AddSource(BlockLatticeManager<T, LatSet1> &LatMan) {
    for (int i = 0; i < _Buoyancys.size(); ++i) {
      _Buoyancys[i].AddSource(&(LatMan.getBlockLat(i)));
    }
  }
  template <typename LatSet1, typename... Args>
  void AddSource(BlockLatticeManager<T, LatSet1> &LatMan, Args... args) {
    AddSource(LatMan);
    AddSource(args...);
  }

  template <typename flagtype>
  void GetBuoyancy(std::int64_t count, std::uint8_t flag,
                   const BlockFieldManager<ScalerField<flagtype>, T, LatSet::d> &BFM) {
    std::uint8_t MaxLevel = NSLatMan.getMaxLevel();
#pragma omp parallel for num_threads(Thread_Num)
    for (int i = 0; i < _Buoyancys.size(); ++i) {
      if (count % (static_cast<int>(pow(2, int(MaxLevel - _Buoyancys[i].getLevel())))) ==
          0)
        _Buoyancys[i].GetBuoyancy(BFM.getBlockField(i).getField().getField(0), flag);
    }
  }

  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T),
            typename flagtype>
  void BGK_U(std::int64_t count, std::uint8_t flag,
             const BlockFieldManager<ScalerField<flagtype>, T, LatSet::d> &BFM) {
    std::uint8_t MaxLevel = NSLatMan.getMaxLevel();
#pragma omp parallel for num_threads(Thread_Num)
    for (int i = 0; i < _Buoyancys.size(); ++i) {
      if (count % (static_cast<int>(pow(2, int(MaxLevel - _Buoyancys[i].getLevel())))) ==
          0)
        _Buoyancys[i].template BGK_U<GetFeq>(BFM.getBlockField(i).getField().getField(0),
                                             flag);
    }
  }

  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T),
            typename flagtype>
  void BGK(std::int64_t count, std::uint8_t flag,
           const BlockFieldManager<ScalerField<flagtype>, T, LatSet::d> &BFM) {
    std::uint8_t MaxLevel = NSLatMan.getMaxLevel();
#pragma omp parallel for num_threads(Thread_Num)
    for (int i = 0; i < _Buoyancys.size(); ++i) {
      if (count % (static_cast<int>(pow(2, int(MaxLevel - _Buoyancys[i].getLevel())))) ==
          0)
        _Buoyancys[i].template BGK<GetFeq>(BFM.getBlockField(i).getField().getField(0),
                                           flag);
    }
  }
};