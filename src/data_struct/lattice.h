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

// lattice.h

#pragma once

#include "data_struct/cell.h"
#include "geometry/geometry2d.h"
#include "geometry/geometry3d.h"
#include "lbm/unit_converter.h"
#include "lbm/collision.h"
#include "lbm/moment.h"
#include "lbm/force.h"
#include "utils/alias.h"

// TODO: use a unified geometry reference like BasicBlock

// a base lattice for PopLattice
template <typename T>
class RhoLattice {
 protected:
  // converter
  AbstractConverter<T>& Conv;
  // rho init
  T Lattice_Rho_Init;
  // buoyancy
  T Lattice_gbeta;

  ScalarField<T> Rho;

 public:
  RhoLattice(AbstractConverter<T>& conv, std::size_t size)
      : Conv(conv), Lattice_Rho_Init(conv.getLatRhoInit()), Lattice_gbeta(conv.getLattice_gbeta()),
        Rho(size, conv.getLatRhoInit()) {}

  ScalarField<T>& getRhoField() { return Rho; }
  const T& getRho(int i) const { return Rho.get(i); }
  T& getRho(int i) { return Rho.get(i); }
  void SetRhoField(std::size_t id, T value) { Rho.SetField(id, value); }

  T getLatRhoInit() const { return Lattice_Rho_Init; }
  T getLatgBeta() const { return Lattice_gbeta; }
};

template <typename T, typename LatSet>
class PopLattice : public RhoLattice<T> {
 protected:
  // domain size
  int Nx;
  int Ny;
  int Nz;
  int N;
  // geometry
  Geometry<T, LatSet::d>& Geo;
  // velocity field
  VectorFieldAOS<T, LatSet::d>& Velocity;
  // population field
  PopulationField<T, LatSet::q> Pops;

  // omega
  T Omega;
  // 1 - omega
  T _Omega;
  // 1 - omega/2
  T fOmega;

  // projection to 1d array
  Vector<int, LatSet::d> Projection;
  // nbr index
  std::array<int, LatSet::q> Delta_Index;

  // tolerance
  T RhoRes;
  std::vector<T> RhoOld;
  T URes;
  std::vector<Vector<T, LatSet::d>> UOld;

  // index
  std::vector<int> Index;
  // index of inner cells
  std::vector<int> InnerIndex;

 public:
  PopLattice(Geometry<T, LatSet::d>& geo, AbstractConverter<T>& conv,
               VectorFieldAOS<T, LatSet::d>& velocity, bool InitIdx = true);

  void InitPop(int Id, T rho) {
    for (int i = 0; i < LatSet::q; ++i) Pops.getField(i)[Id] = rho * LatSet::w[i];
  }
  std::array<T*, LatSet::q> getPop(std::size_t id) { return Pops.getArray(id); }

  T& getPopdir(std::size_t id, int dir) { return Pops.getField(dir)[id]; }
  const T& getPopdir(std::size_t id, int dir) const { return Pops.getField(dir)[id]; }

  PopCell<T, LatSet> getNeighbor(const PopCell<T, LatSet>& cell, int i) const {
    return PopCell<T, LatSet>(cell.getId() + Delta_Index[i], *this);
  }
  PopCell<T, LatSet> getNeighbor(const PopCell<T, LatSet>& cell,
                              const Vector<int, LatSet::d>& direction) const {
    return PopCell<T, LatSet>(cell.getId() + direction * Projection, *this);
  }

  const std::array<int, LatSet::q>& getDelta_Index() const { return Delta_Index; }
  int getDelta_Index(int dir) const { return Delta_Index[dir]; }
  std::size_t getNbrId(std::size_t id, int dir) const { return id + Delta_Index[dir]; }
  PopulationField<T, LatSet::q>& getPopField() { return Pops; }
  VectorFieldAOS<T, LatSet::d>& getVelocityField() { return Velocity; }
  Geometry<T, LatSet::d>& getGeo() { return Geo; }
  std::vector<int>& getIndex() { return Index; }
  std::vector<int>& getInnerIndex() { return InnerIndex; }

  int getNx() const { return Nx; }
  int getNy() const { return Ny; }
  int getNz() const { return Nz; }
  int getN() const { return N; }
  inline T getOmega() const { return Omega; }
  inline T get_Omega() const { return _Omega; }
  inline T getfOmega() const { return fOmega; }
  const Vector<int, LatSet::d>& getProjection() const { return Projection; }
  Vector<T, LatSet::d>& getVelocity(int i) { return Velocity.get(i); }

  void UpdateRho(const std::vector<int>& index);
  // update rho based on a flag field of flagtype(usually std::uint8_t or enum
  // of std::uint8_t)
  template <typename ArrayType>
  void UpdateRho(const ArrayType& flagarr, std::uint8_t flag);
  // update rho based on a flag field of flagtype(usually std::uint8_t or enum
  // of std::uint8_t) from a source
  template <typename ArrayType>
  void UpdateRho_Source(const ArrayType& flagarr, std::uint8_t flag,
                        const GenericArray<T>& source);

  void UpdateU(const std::vector<int>& index);
  template <typename ArrayType>
  void UpdateU(const ArrayType& flagarr, std::uint8_t flag);

  void UpdateRhoU(const std::vector<int>& index);

  // dynamics
  template <void (*get_feq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  void BGK();
  template <void (*get_feq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  void BGK(const std::vector<int>& index);
  template <void (*get_feq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T), typename ArrayType>
  void BGK(const ArrayType& flagarr, std::uint8_t flag);
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T), typename ArrayType>
  void BGK_Source(const ArrayType& flagarr, std::uint8_t flag,
                  const GenericArray<T>& source);

  void Stream();
  // tolerance
  void EnableToleranceRho(T rhores = T(1e-5));
  void EnableToleranceU(T ures = T(1e-5));
  T getToleranceRho();
  T getToleranceU();
  // experimental
  T getTolRho(int shift = 1);
  T getTolU(int shift = 1);
};
