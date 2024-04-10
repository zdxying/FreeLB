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

// free surface model
#pragma once

#include "legacy/legacy_lbm/lbm2D.h"
#include "data_struct/lattice.h"

struct NbrInfo {
  bool hasFluidNbr = false;
  bool hasGasNbr = false;
  int InterfaceNbrNum = 0;
};

enum Type : std::uint8_t { Gas = 1, Interface = 2, Fluid = 4, Solid = 8 };

enum Flag : std::uint8_t { None = 0, ToGas = 1, ToFluid = 2, NewInterface = 4 };

template <typename T, typename LatSet>
class FreeSurface2D {
 private:
  BasicLattice<T, LatSet>& NS;
  Geometry2D<T>& Geo;
  int N;

  // index
  std::vector<int> InterfaceIdx;
  std::vector<int> FluidIdx;

  // fluid fraction
  ScalerField<T> Fv;
  // cell type
  ScalerField<Type> Types;
  // cell flag
  ScalerField<Flag> Flags;
  // mass
  ScalerField<T> Mass;
  VectorFieldAOS<T, LatSet::q> TEMP_MASS_EXCHANGE;
  VectorFieldAOS<T, LatSet::d> PREVIOUS_VELOCITY;

  // FREE SURFACE PARAMETERS
  T surface_tension_parameter;  // surface_tension_coefficient_factor *
                                // surface_tension_coefficient
  T lonely_threshold = 1.0;
  T transitionThreshold = 1e-3;
  bool drop_isolated_cells = true;
  bool has_surface_tension = false;

 public:
  FreeSurface2D(BasicLattice<T, LatSet>& ns)
      : NS(ns),
        Geo(ns.getGeo()),
        N(ns.getGeo().getVoxelsNum()),
        Fv(ns.getGeo().getVoxelsNum(), T(0)),
        Types(ns.getGeo().getVoxelsNum(), Type::Gas),
        Flags(ns.getGeo().getVoxelsNum(), Flag::None),
        Mass(ns.getGeo().getVoxelsNum(), T(0)),
        TEMP_MASS_EXCHANGE(ns.getGeo().getVoxelsNum(), Vector<T, LatSet::q>{}),
        PREVIOUS_VELOCITY(ns.getGeo().getVoxelsNum(), Vector<T, LatSet::d>{}) {
    InterfaceIdx.reserve(N);
    FluidIdx.reserve(N);
  }
  void Initialize();
  void setSrufaceTension(T value) {
    if (value > 1e-3) {
      std::cout << "Surface tension: "<<value<<" may be too large" << std::endl;
    }
    surface_tension_parameter = value;
    has_surface_tension = true;
  }
  // DO NOT call this on computation boundary, i.e. solid cells
  NbrInfo getNbrInfo(int id);
  bool hasNeighborType(int id, std::uint8_t type);
  bool hasNeighborFlag(int id, std::uint8_t flag);
  bool isType(int id, std::uint8_t type);
  bool isFlag(int id, std::uint8_t flag);
  void setType(int id, std::uint8_t type);
  void setFlag(int id, std::uint8_t flag);

  inline T getF(int id) { return std::max(std::min(Fv.get(id), T(1)), T(0)); }
  T calculateSurfaceTensionCurvature2D(int id);
  Vector<T, LatSet::d> computeParkerYoungInterfaceNormal(int id);
  T calculateCubeOffset(T volume, const Vector<T, LatSet::d>& normal);
  template <int S>
  std::array<T, S> solvePivotedLU(std::array<std::array<T, S>, S>& matrix,
                                  const std::array<T, S>& b, int N);

  void Prepare();
  void ToFluidCellConversion();
  void ToGasCellConversion();
  void MassExcess();
  void FinalizeConversion();

  // free surface lbm, this should be called after the stream and bcs step
  void Post_Stream();
  
  void Apply(){   
    Post_Stream();
    ToFluidCellConversion();
    ToGasCellConversion();
    MassExcess();
    FinalizeConversion();
    // Prepare();
  }

  // get
  ScalerField<Type>& getType() { return Types; }
  ScalerField<Flag>& getFlag() { return Flags; }
  ScalerField<T>& getMass() { return Mass; }
  ScalerField<T>& getF() { return Fv; }
  const std::vector<int>& getInterfaceIdx() { return InterfaceIdx; }
  const std::vector<int>& getFluidIdx() { return FluidIdx; }

  // set
};

namespace FreeSurface {

enum Type : std::uint16_t {
  Fluid = 1,
  Interface = 2,
  Gas = 4,
  Solid = 8,
  NO_FLUID_NEIGH = 16,
  NO_EMPTY_NEIGH = 32,
  NO_IFACE_NEIGH = 64,
  TO_FLUID = 128,
  TO_GAS = 256
};

template <typename T, typename LatSet>
class FreeSurface2D {
 private:
  Genericlbm2D<T, LatSet>& NS;
  VoxelGeometry2D<T>& Geo;
  std::vector<Voxel<T, LatSet::d>>& Voxels;

  // index
  // std::vector<int> InterfaceIdx;
  // std::vector<int> FluidIdx;

  // fluid fraction
  std::vector<FreeSurface::Type> _Type;
  // std::vector<Flag> _Flag;
  std::vector<T> _Mass;
  std::vector<T> mass_prev;
  std::vector<T> rho_prev;
  std::vector<Vector<T, LatSet::d>> u_prev;
  std::vector<FreeSurface::Type> cell_type_prev;
  std::vector<Vector<T, LatSet::q>> TEMP_MASS_EXCHANGE;
  std::vector<Vector<T, LatSet::d>> PREVIOUS_VELOCITY;

  // FREE SURFACE PARAMETERS
  T surface_tension_parameter;  // surface_tension_coefficient_factor *
                                // surface_tension_coefficient
  T lonely_threshold = 0.1;
  T transitionThreshold = 0.03;
  bool drop_isolated_cells = true;
  bool has_surface_tension = false;

 public:
  FreeSurface2D(Genericlbm2D<T, LatSet>& ns)
      : NS(ns), Geo(ns.getGeo()), Voxels(ns.getGeo().getVoxels()) {
    _Type.resize(Voxels.size(), FreeSurface::Type::Gas);
    _Mass.resize(Voxels.size(), T(0));
    TEMP_MASS_EXCHANGE.resize(Voxels.size(), Vector<T, LatSet::q>{});
    PREVIOUS_VELOCITY.resize(Voxels.size(), Vector<T, LatSet::d>{});
    mass_prev.resize(Voxels.size(), T(0));
    rho_prev.resize(Voxels.size(), T(0));
    u_prev.resize(Voxels.size(), Vector<T, LatSet::d>{});
    cell_type_prev.resize(Voxels.size(), FreeSurface::Type::Gas);
  }

  void Initialize();
  void setSrufaceTension(T value) {
    if (value > 1e-3) {
      std::cout << "Surface tension may be too large" << std::endl;
    }
    surface_tension_parameter = value;
    has_surface_tension = true;
  }

  // NbrInfo getNbrInfo(int id);
  bool hasNeighborType(int id, const FreeSurface::Type& type);
  // bool hasNeighborFlag(int id, const Flag& flag);
  bool isType(int id, const FreeSurface::Type& type);
  // bool isFlag(int id, const Flag& flag);
  void setType(int id, const FreeSurface::Type& type);
  // void setFlag(int id, const Flag& flag);

  T getEpsilon(const FreeSurface::Type& tp, const T rho_prev,
               const T mass_prev) const {
    // cell type:
    //  const FreeSurface::Type& type = _Type[id];
    if (static_cast<bool>(
            tp & (FreeSurface::Type::Fluid | FreeSurface::Type::Solid))) {
      return T(1);
    } else if (static_cast<bool>(tp & FreeSurface::Type::Gas)) {
      return T(0);
    } else {
      if (rho_prev > 0) {
        T eps = mass_prev / rho_prev;
        return std::max(std::min(eps, T(1)), T(0));
      } else {
        return T(0.5);
      }
    }
  }

  void Collide(T omega, T _omega, const Vector<T, LatSet::d>& force);
  void Stream();
  void Conversion();
  Vector<T, LatSet::d> getNormal(int id);
  Vector<T, LatSet::d> get_normal(int id,
                                  const std::vector<FreeSurface::Type>& types,
                                  const std::vector<T>& rho_prev,
                                  const std::vector<T>& mass_prev);

  void Apply() {
    // Prepare() will be called at the beginning of each time step
    // Collide();
    Stream();
    Conversion();
  }

  // get
  std::vector<FreeSurface::Type>& getType() { return _Type; }
  std::vector<T>& getMass() { return _Mass; }

  // set
};
}  // namespace FreeSurface
