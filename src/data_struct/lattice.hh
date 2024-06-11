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

// lattice.hh

#pragma once

#include "data_struct/lattice.h"
#include "lattice.h"

template <typename T, typename LatSet>
PopLattice<T, LatSet>::PopLattice(Geometry<T, LatSet::d>& geo, AbstractConverter<T>& conv,
                                      VectorFieldAOS<T, LatSet::d>& velocity, bool InitIdx)
    : Geo(geo), Nx(geo.getNx()), Ny(geo.getNy()), Nz(geo.getNz()), Pops(geo.getVoxelsNum()),
      Omega(conv.getOMEGA()), Velocity(velocity), RhoLattice<T>(conv, geo.getVoxelsNum()) {
  _Omega = T(1) - Omega;
  fOmega = T(1) - T(0.5) * Omega;
  if constexpr (LatSet::d == 3) {
    Projection = Vector<int, 3>(1, Nx, Nx * Ny);
    InnerIndex.reserve((Nx - 2) * (Ny - 2) * (Nz - 2));
  } else {
    Projection = Vector<int, 2>(1, Nx);
    InnerIndex.reserve((Nx - 2) * (Ny - 2));
  }
  Delta_Index = make_Array<int, LatSet::q>([&](int i) { return LatSet::c[i] * Projection; });
  N = Nx * Ny * Nz;
  // init populations
  for (int i = 0; i < LatSet::q; ++i) {
    Pops.getField(i).Init(this->Lattice_Rho_Init * LatSet::w[i]);
  }
  if (InitIdx) {
    Index.reserve(N);
    // init index
    Geo.getGeoFlagField().getField().for_isNotflag(Geo.getVoidflag(),
                                                   [this](int id) { Index.push_back(id); });
    // init inner index
    Geo.getGeoFlagField().getField().for_isflag(Geo.getAABBflag(),
                                                [this](int id) { InnerIndex.push_back(id); });
  }
}
template <typename T, typename LatSet>
void PopLattice<T, LatSet>::UpdateRho(const std::vector<int>& index) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : index) {
    BasicPopCell<T, LatSet> cell(id, *this);
    moment::Rho<T, LatSet>::apply(cell, this->Rho.get(id));
  }
}


template <typename T, typename LatSet>
template <typename ArrayType>
void PopLattice<T, LatSet>::UpdateRho(const ArrayType& flagarr, std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BasicPopCell<T, LatSet> cell(id, *this);
      moment::Rho<T, LatSet>::apply(cell, this->Rho.get(id));
    }
  }
}

template <typename T, typename LatSet>
template <typename ArrayType>
void PopLattice<T, LatSet>::UpdateRho_Source(const ArrayType& flagarr,
                                               std::uint8_t flag, const GenericArray<T>& source) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BasicPopCell<T, LatSet> cell(id, *this);
      moment::Rho<T, LatSet>::apply(cell, this->Rho.get(id), source[id]);
    }
  }
}

template <typename T, typename LatSet>
template <typename ArrayType>
void PopLattice<T, LatSet>::UpdateU(const ArrayType& flagarr, std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      BasicPopCell<T, LatSet> cell(id, *this);
      moment::Velocity<T, LatSet>::apply(cell, Velocity.get(id));
    }
  }
}

template <typename T, typename LatSet>
void PopLattice<T, LatSet>::UpdateU(const std::vector<int>& index) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : index) {
    BasicPopCell<T, LatSet> cell(id, *this);
    moment::Velocity<T, LatSet>::apply(cell, Velocity.get(id));
  }
}

template <typename T, typename LatSet>
void PopLattice<T, LatSet>::UpdateRhoU(const std::vector<int>& index) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : index) {
    BasicPopCell<T, LatSet> cell(id, *this);
    moment::RhoVelocity<T, LatSet>::apply(cell, this->Rho.get(id), Velocity.get(id));
  }
}


template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
void PopLattice<T, LatSet>::BGK() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    PopCell<T, LatSet> cell(id, *this);
    collision::BGK<T, LatSet>::template apply<GetFeq>(cell);
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
void PopLattice<T, LatSet>::BGK(const std::vector<int>& index) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : index) {
    PopCell<T, LatSet> cell(id, *this);
    collision::BGK<T, LatSet>::template apply<GetFeq>(cell);
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename ArrayType>
void PopLattice<T, LatSet>::BGK(const ArrayType& flagarr, std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      PopCell<T, LatSet> cell(id, *this);
      collision::BGK<T, LatSet>::template apply<GetFeq>(cell);
    }
  }
}

template <typename T, typename LatSet>
template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
          typename ArrayType>
void PopLattice<T, LatSet>::BGK_Source(const ArrayType& flagarr, std::uint8_t flag,
                                         const GenericArray<T>& source) {
  T fOmega = T(1) - Omega * T(0.5);
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    if (util::isFlag(flagarr[id], flag)) {
      PopCell<T, LatSet> cell(id, *this);
      collision::BGK<T, LatSet>::template applySource<GetFeq>(cell,
                                                              source[id]);
    }
  }
}

template <typename T, typename LatSet>
void PopLattice<T, LatSet>::Stream() {
  for (int i = 1; i < LatSet::q; ++i) {
    Pops.getField(i).rotate(LatSet::c[i] * Projection);
  }
}

template <typename T, typename LatSet>
inline void PopLattice<T, LatSet>::EnableToleranceRho(T rhores) {
  RhoRes = rhores;
  RhoOld.reserve(N);
  for (int i = 0; i < N; ++i) RhoOld.push_back(this->Rho.get(i));
}

template <typename T, typename LatSet>
inline void PopLattice<T, LatSet>::EnableToleranceU(T ures) {
  URes = ures;
  UOld.reserve(N);
  for (int i = 0; i < N; ++i) UOld.push_back(Velocity.get(i));
}

template <typename T, typename LatSet>
inline T PopLattice<T, LatSet>::getToleranceRho() {
  T res;
  T maxres = T(0);
  // rho and u with void flag should never be updated
  for (int i = 0; i < N; ++i) {
    res = std::abs(this->Rho.get(i) - RhoOld[i]);
    maxres = std::max(res, maxres);
    RhoOld[i] = this->Rho.get(i);
  }
  return maxres;
}

template <typename T, typename LatSet>
inline T PopLattice<T, LatSet>::getToleranceU() {
  T res0, res1, res;
  T maxres = T(0);
  for (int i = 0; i < N; ++i) {
    res0 = std::abs(Velocity.get(i)[0] - UOld[i][0]);
    res1 = std::abs(Velocity.get(i)[1] - UOld[i][1]);
    res = std::max(res0, res1);
    maxres = std::max(res, maxres);
    // set UOld
    UOld[i][0] = Velocity.get(i)[0];
    UOld[i][1] = Velocity.get(i)[1];
  }
  return maxres;
}

template <typename T, typename LatSet>
T PopLattice<T, LatSet>::getTolRho(int shift) {
  T res;
  T maxres = T(0);

  if constexpr (LatSet::d == 2) {
    for (int j = shift; j < Ny - shift; ++j) {
      for (int i = shift; i < Nx - shift; ++i) {
        int id = j * Nx + i;
        res = std::abs(this->Rho.get(id) - RhoOld[id]);
        maxres = std::max(res, maxres);
        RhoOld[id] = this->Rho.get(id);
      }
    }
  } else if constexpr (LatSet::d == 3) {
    int NxNy = Nx * Ny;
    for (int k = shift; k < Nz - shift; ++k) {
      for (int j = shift; j < Ny - shift; ++j) {
        for (int i = shift; i < Nx - shift; ++i) {
          int id = k * NxNy + j * Nx + i;
          res = std::abs(this->Rho.get(id) - RhoOld[id]);
          maxres = std::max(res, maxres);
          RhoOld[id] = this->Rho.get(id);
        }
      }
    }
  }
  return maxres;
}

template <typename T, typename LatSet>
T PopLattice<T, LatSet>::getTolU(int shift) {
  T res0, res1, res;
  T maxres = T(0);
  if constexpr (LatSet::d == 2) {
    for (int j = shift; j < Ny - shift; ++j) {
      for (int i = shift; i < Nx - shift; ++i) {
        int id = j * Nx + i;
        res0 = std::abs(Velocity.get(id)[0] - UOld[id][0]);
        res1 = std::abs(Velocity.get(id)[1] - UOld[id][1]);
        res = std::max(res0, res1);
        maxres = std::max(res, maxres);
        // set UOld
        UOld[id][0] = Velocity.get(id)[0];
        UOld[id][1] = Velocity.get(id)[1];
      }
    }
  } else if constexpr (LatSet::d == 3) {
    int NxNy = Nx * Ny;
    for (int k = shift; k < Nz - shift; ++k) {
      for (int j = shift; j < Ny - shift; ++j) {
        for (int i = shift; i < Nx - shift; ++i) {
          int id = k * NxNy + j * Nx + i;
          res0 = std::abs(Velocity.get(id)[0] - UOld[id][0]);
          res1 = std::abs(Velocity.get(id)[1] - UOld[id][1]);
          res = std::max(res0, res1);
          maxres = std::max(res, maxres);
          // set UOld
          UOld[id][0] = Velocity.get(id)[0];
          UOld[id][1] = Velocity.get(id)[1];
        }
      }
    }
  }
  return maxres;
}
