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

#include "data_struct/lattice.h"

// namespace FreeSurface

namespace FS {

enum FSType : std::uint16_t {
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
  Block2D<T>& Geo;
  std::vector<Voxel<T, LatSet::d>>& Voxels;

  // index
  // std::vector<int> InterfaceIdx;
  // std::vector<int> FluidIdx;

  // fluid fraction
  std::vector<FSType> _Type;
  // std::vector<Flag> _Flag;
  std::vector<T> _Mass;
  std::vector<T> mass_prev;
  std::vector<T> rho_prev;
  std::vector<Vector<T, LatSet::d>> u_prev;
  std::vector<FSType> cell_type_prev;
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
    _Type.resize(Voxels.size(), FSType::Gas);
    _Mass.resize(Voxels.size(), T(0));
    TEMP_MASS_EXCHANGE.resize(Voxels.size(), Vector<T, LatSet::q>{});
    PREVIOUS_VELOCITY.resize(Voxels.size(), Vector<T, LatSet::d>{});
    mass_prev.resize(Voxels.size(), T(0));
    rho_prev.resize(Voxels.size(), T(0));
    u_prev.resize(Voxels.size(), Vector<T, LatSet::d>{});
    cell_type_prev.resize(Voxels.size(), FSType::Gas);
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
  bool hasNeighborType(int id, const FSType& type);
  // bool hasNeighborFlag(int id, const Flag& flag);
  bool isType(int id, const FSType& type);
  // bool isFlag(int id, const Flag& flag);
  void setType(int id, const FSType& type);
  // void setFlag(int id, const Flag& flag);

  T getEpsilon(const FSType& tp, const T rho_prev,
               const T mass_prev) const {
    // cell type:
    //  const FSType& type = _Type[id];
    if (static_cast<bool>(
            tp & (FSType::Fluid | FSType::Solid))) {
      return T(1);
    } else if (static_cast<bool>(tp & FSType::Gas)) {
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
                                  const std::vector<FSType>& types,
                                  const std::vector<T>& rho_prev,
                                  const std::vector<T>& mass_prev);

  void Apply() {
    // Prepare() will be called at the beginning of each time step
    // Collide();
    Stream();
    Conversion();
  }

  // get
  std::vector<FSType>& getType() { return _Type; }
  std::vector<T>& getMass() { return _Mass; }

  // set
};

template <typename T, typename LatSet>
class FreeSurface2DManager {
};


}  // namespace FreeSurface
