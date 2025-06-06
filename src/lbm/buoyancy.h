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

#include "data_struct/lattice.h"

// buoyancy for basiclattice
template <typename T, typename LatSet>
class Buoyancy {
 private:
  // NS lattice
  PopLattice<T, LatSet> &NSLat;
  // source lattice, e.g., thermal and solute lattice
  std::vector<RhoLattice<T> *> Source;
  VectorFieldAOS<T, LatSet::d> &Velocity;
  
  // omega
  T Omega;
  // 1 - omega
  T _Omega;
  // 1 - omega/2
  T fOmega;

  ScalarField<T> Force;

 public:
  Buoyancy(PopLattice<T, LatSet> &lat, VectorFieldAOS<T, LatSet::d> &velocity)
      : NSLat(lat), Velocity(velocity), Omega(lat.getOmega()), _Omega(lat.get_Omega()),
        fOmega(T(1) - lat.getOmega() * T(0.5)), Force(lat.getN(), T(0)) {
  }
  template <typename... Args>
  void AddSource(RhoLattice<T> *lat, Args*... args) {
    Source.push_back(lat);
    AddSource(args...);
  }
  void AddSource(RhoLattice<T> *lat) { Source.push_back(lat); }
  void GetBuoyancy() {
    // reset force
    Force.getField(0).Init(T(0));
    // add to buoyancy
    for (RhoLattice<T> *lat : Source) {
      T latRhoInit = lat->getLatRhoInit();
      T latgbeta = lat->getLatgBeta();
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
      for (int id : NSLat.getIndex()) {
        Force.get(id) += (lat->getRho(id) - latRhoInit) * latgbeta;
      }
    }
  }
  template <typename ArrayType>
  void GetBuoyancy(const ArrayType &flagarr, std::uint8_t flag) {
    // reset force
    Force.getField(0).Init(T(0));
    // add to buoyancy
    for (RhoLattice<T> *lat : Source) {
      T latRhoInit = lat->getLatRhoInit();
      T latgbeta = lat->getLatgBeta();
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
      for (int id = 0; id < NSLat.getN(); ++id) {
        if (util::isFlag(flagarr[id], flag))
          Force.get(id) += (lat->getRho(id) - latRhoInit) * latgbeta;
      }
    }
  }
  Vector<T, LatSet::d> getForce(int id) const {
    if constexpr (LatSet::d == 3) {
      return Vector<T, LatSet::d>{T(0), T(0), Force.get(id)};
    } else if constexpr (LatSet::d == 2) {
      return Vector<T, LatSet::d>{T(0), Force.get(id)};
    }
  }
  // BGK with FORCE term
  // update force and u at the same time
  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T)>
  void BGK_U(const std::vector<int> &index) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id : index) {
      PopCell<T, LatSet> cell(id, NSLat);
      std::array<T, LatSet::q> Fi{};
      force::ForcePop<T, LatSet>::compute(Fi, Velocity.get(id), getForce(id));
      moment::Velocity<T, LatSet>::apply(cell, Velocity.get(id), Fi);
      legacy::BGK<T, LatSet>::template applySource<GetFeq>(cell, Fi);
    }
  }
  // BGK with FORCE term
  // do not update u
  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T)>
  void BGK(const std::vector<int> &index) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id : index) {
      PopCell<T, LatSet> cell(id, NSLat);
        std::array<T, LatSet::q> Fi{};
        force::ForcePop<T, LatSet>::compute(Fi, Velocity.get(id), getForce(id));
        legacy::BGK<T, LatSet>::template applySource<GetFeq>(cell, Fi);
    }
  }
  // BGK with FORCE term
  // update force and u at the same time
  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T),
            typename ArrayType>
  void BGK_U(const ArrayType &flagarr, std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id = 0; id < NSLat.getN(); ++id) {
      if (util::isFlag(flagarr[id], flag)) {
        PopCell<T, LatSet> cell(id, NSLat);
        std::array<T, LatSet::q> Fi{};
        force::ForcePop<T, LatSet>::compute(Fi, Velocity.get(id), getForce(id));
        moment::Velocity<T, LatSet>::apply(cell, Velocity.get(id), Fi);
        legacy::BGK<T, LatSet>::template applySource<GetFeq>(cell, Fi);
      }
    }
  }
  // BGK with FORCE term
  // do not update u
  template <void (*GetFeq)(std::array<T, LatSet::q> &, const Vector<T, LatSet::d> &, T), typename ArrayType>
  void BGK(const ArrayType &flagarr, std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id = 0; id < NSLat.getN(); ++id) {
      if (util::isFlag(flagarr[id], flag)) {
        PopCell<T, LatSet> cell(id, NSLat);
        std::array<T, LatSet::q> Fi{};
        force::ForcePop<T, LatSet>::compute(Fi, Velocity.get(id), getForce(id));
        legacy::BGK<T, LatSet>::template applySource<GetFeq>(cell, Fi);
      }
    }
  }
  void UpdateU(const std::vector<int> &index) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id : index) {
      BasicPopCell<T, LatSet> cell(id, NSLat);
      std::array<T, LatSet::q> Fi{};
      force::ForcePop<T, LatSet>::compute(Fi, Velocity.get(id), getForce(id));
      moment::Velocity<T, LatSet>::apply(cell, Velocity.get(id), Fi);
    }
  }
};