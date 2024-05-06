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

#include "lbm/freeSurface.h"

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::Initialize() {
  for (int Id = 0; Id < N; ++Id) {
    Mass.SetField(Id, T(0));
    if (isType(Id, FSType::Interface)) {
      Fv.SetField(Id, T(0.5));
      Mass.SetField(Id, T(0.5));
    } else if (isType(Id, FSType::Gas)) {
      Fv.SetField(Id, T(0));
      Mass.SetField(Id, T(0));
    } else if (isType(Id, FSType::Fluid)) {
      Fv.SetField(Id, T(1));
      Mass.SetField(Id, T(1));
    }
  }
  Prepare();
}

template <typename T, typename LatSet>
NbrInfo FreeSurface2D<T, LatSet>::getNbrInfo(int id) {
  NbrInfo nbrInfo;
  for (int i = 1; i < LatSet::q; ++i) {
    int nbrId = Geo.template getNeighborId<LatSet>(id, i);
    if (isType(nbrId, FSType::Gas)) {
      nbrInfo.hasGasNbr = true;
    } else if (isType(nbrId, FSType::Fluid)) {
      nbrInfo.hasFluidNbr = true;
    } else if (isType(nbrId, FSType::Interface)) {
      ++nbrInfo.InterfaceNbrNum;
    }
  }
  return nbrInfo;
}
template <typename T, typename LatSet>
inline bool FreeSurface2D<T, LatSet>::hasNeighborType(int id,
                                                      std::uint8_t type) {
  for (int i = 1; i < LatSet::q; ++i) {
    if (isType(Geo.template getNeighborId<LatSet>(id, i), type)) return true;
  }
  return false;
}
template <typename T, typename LatSet>
inline bool FreeSurface2D<T, LatSet>::hasNeighborFlag(int id,
                                                      std::uint8_t flag) {
  for (int i = 1; i < LatSet::q; ++i) {
    if (isFlag(NS.getNbrId(id, i), flag)) return true;
  }
  return false;
}

template <typename T, typename LatSet>
inline bool FreeSurface2D<T, LatSet>::isType(int id, std::uint8_t type) {
  return static_cast<bool>(Types.get(id) & type);
}
template <typename T, typename LatSet>
inline bool FreeSurface2D<T, LatSet>::isFlag(int id, std::uint8_t flag) {
  return static_cast<bool>(Flags.get(id) & flag);
}
template <typename T, typename LatSet>
inline void FreeSurface2D<T, LatSet>::setType(int id, std::uint8_t type) {
  Types.SetField(id, type);
}
template <typename T, typename LatSet>
inline void FreeSurface2D<T, LatSet>::setFlag(int id, std::uint8_t flag) {
  Flags.SetField(id, flag);
}

template <typename T, typename LatSet>
T FreeSurface2D<T, LatSet>::calculateSurfaceTensionCurvature2D(int id) {
  Vector<T, LatSet::d> normal = computeParkerYoungInterfaceNormal(id);

  T norm = normal.getnorm();
  if (norm < 1e-6) {
    return T(0);
  }
  normal = normal / norm;
  // Rotation matrix:
  // ( n1 | -n0 )
  // ( n0 |  n1 )

  constexpr int S = 2;
  std::array<std::array<T, S>, S> lq_matrix;
  std::array<T, S> b_rhs;
  for (int i = 0; i < S; ++i) {
    for (int j = 0; j < S; ++j) {
      lq_matrix[i][j] = T(0);
    }
    b_rhs[i] = T(0);
  }
  // Offset for the plic correction
  T origin_offset = T(0);
  T fill_level = getF(id);
  origin_offset = calculateCubeOffset(fill_level, normal);
  // end offset

  int healthy_interfaces = 0;

  for (int i = 1; i < LatSet::q; ++i) {
    int nbrId = Geo.template getNeighborId<LatSet>(id, i);

    if (!isType(nbrId, FSType::Interface) || !hasNeighborType(nbrId, FSType::Gas)) {
      continue;
    }

    ++healthy_interfaces;

    T fill_level = getF(nbrId);

    T cube_offset = calculateCubeOffset(fill_level, normal);

    T x_pos = LatSet::c[i][0];
    T y_pos = LatSet::c[i][1];

    // Rotation
    T rot_x_pos = x_pos * normal[1] - y_pos * normal[0];
    T rot_y_pos =
        x_pos * normal[0] + y_pos * normal[1] + (cube_offset - origin_offset);

    T rot_x_pos_2 = rot_x_pos * rot_x_pos;
    T rot_x_pos_3 = rot_x_pos_2 * rot_x_pos;
    T rot_x_pos_4 = rot_x_pos_3 * rot_x_pos;

    lq_matrix[1][1] += rot_x_pos_2;
    lq_matrix[1][0] += rot_x_pos_3;
    lq_matrix[0][0] += rot_x_pos_4;

    b_rhs[0] += rot_x_pos_2 * (rot_y_pos);
    b_rhs[1] += rot_x_pos * (rot_y_pos);
  }

  lq_matrix[0][1] = lq_matrix[1][0];

  // Thikonov regularization parameter
  T alpha = 0.0;
  for (int i = 0; i < LatSet::d; ++i) {
    lq_matrix[i][i] += alpha;
  }

  // It is 2 because of the fitting parameters. Not dependent on the dimension
  std::array<T, S> solved_fit =
      solvePivotedLU<S>(lq_matrix, b_rhs, healthy_interfaces);

  // signed curvature -> kappa = y'' / ( (1 + y'²)^(3/2) )
  T denom = std::sqrt(1. + solved_fit[1] * solved_fit[1]);
  denom = denom * denom * denom;
  T curvature = 2. * solved_fit[0] / denom;
  return std::max(-1., std::min(1., curvature));
}

template <typename T, typename LatSet>
Vector<T, LatSet::d>
FreeSurface2D<T, LatSet>::computeParkerYoungInterfaceNormal(int id) {
  Vector<T, LatSet::d> normal{};
  for (int i = 1; i < LatSet::q; ++i) {
    int omega_weight = 1;
    if (LatSet::c[i][0] != 0) {
      omega_weight *= 2;
    }
    if (LatSet::c[i][1] != 0) {
      omega_weight *= 2;
    }
    if constexpr (LatSet::d == 3) {
      if (LatSet::c[i][2] != 0) {
        omega_weight *= 2;
      }
    }
    omega_weight /= 2;
    T f = getF(Geo.template getNeighborId<LatSet>(id, i));
    normal[0] -= omega_weight * LatSet::c[i][0] * f;
    normal[1] -= omega_weight * LatSet::c[i][1] * f;
    if constexpr (LatSet::d == 3) {
      normal[2] -= omega_weight * LatSet::c[i][2] * f;
    }
  }
  return normal;
}

template <typename T, typename LatSet>
T FreeSurface2D<T, LatSet>::calculateCubeOffset(
    T volume, const Vector<T, LatSet::d> &normal) {
  std::vector<T> abs_normal(LatSet::d, T{0});

  for (int i = 0; i < LatSet::d; i++) {
    abs_normal[i] = std::abs(normal[i]);
  }

  T volume_symmetry = 0.5 - std::abs(volume - 0.5);

  std::sort(abs_normal.begin(), abs_normal.end());

  if constexpr (LatSet::d == 2) {
    abs_normal[0] = std::max(normal[0], 1e-5);
  } else if (LatSet::d == 3) {
    abs_normal[0] = std::max(normal[0], 1e-12);
    abs_normal[1] = std::max(normal[1], 1e-12);
  }
  T d = T(0);
  // offsetHelper
  T d2 = volume_symmetry * abs_normal[1] + 0.5 * abs_normal[0];
  if (d2 >= abs_normal[0]) {
    d = d2;
  } else {
    d = std::sqrt(2. * abs_normal[0] * abs_normal[1] * volume_symmetry);
  }
  // end offsetHelper

  T sorted_normal_acc = 0;
  for (int i = 0; i < LatSet::d; i++) {
    sorted_normal_acc += abs_normal[i];
  }

  return std::copysign(d - 0.5 * sorted_normal_acc, volume - 0.5);
}

template <typename T, typename LatSet>
template <int S>
inline std::array<T, S> FreeSurface2D<T, LatSet>::solvePivotedLU(
    std::array<std::array<T, S>, S> &matrix, const std::array<T, S> &b, int N) {
  std::array<T, S> x;
  std::array<T, S> pivots;
  for (int i = 0; i < S; ++i) {
    pivots[i] = i;
    x[i] = 0.;
  }

  N = std::min(N, S);

  for (int i = 0; i < N; ++i) {
    T max = 0.;
    int max_index = i;

    for (int j = i; j < N; ++j) {
      T abs = std::abs(matrix[pivots[j]][i]);
      if (abs > max) {
        max_index = j;
        max = abs;
      }
    }

    if (max_index != i) {
      int tmp_index = pivots[i];
      pivots[i] = pivots[max_index];
      pivots[max_index] = tmp_index;
    }

    for (int j = i + 1; j < N; ++j) {
      matrix[pivots[j]][i] /= matrix[pivots[i]][i];

      for (int k = i + 1; k < N; ++k) {
        matrix[pivots[j]][k] -= matrix[pivots[j]][i] * matrix[pivots[i]][k];
      }
    }
  }

  for (int i = 0; i < N; ++i) {
    x[i] = b[pivots[i]];

    for (int j = 0; j < i; ++j) {
      x[i] -= matrix[pivots[i]][j] * x[j];
    }
  }

  for (int i = N; i > 0; --i) {
    for (int j = i; j < N; ++j) {
      x[i - 1] -= matrix[pivots[i - 1]][j] * x[j];
    }

    x[i - 1] /= matrix[pivots[i - 1]][i - 1];
  }

  return x;
}

template <typename T, typename LatSet>
inline void FreeSurface2D<T, LatSet>::Prepare() {
  InterfaceIdx.clear();
  FluidIdx.clear();
  for (int Id = 0; Id < N; ++Id) {
    // Remove all cell flags
    setFlag(Id, Flag::None);
    // add to index
    if (isType(Id, FSType::Fluid))
      FluidIdx.push_back(Id);
    else if (isType(Id, FSType::Interface))
      InterfaceIdx.push_back(Id);
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::Post_Stream() {
#pragma omp parallel for num_threads(Thread_Num)
  for (int Id : FluidIdx) {
    const BasicCell<T, LatSet> cell(Id, NS);
    for (int i = 1; i < LatSet::q; ++i) {
      int nbrId = Geo.template getNeighborId<LatSet>(Id, i);
      if (isType(nbrId, static_cast<FSType>(FSType::Fluid | FSType::Interface))) {
        const BasicCell<T, LatSet> celli(nbrId, NS);
        Mass.get(Id) += cell[LatSet::opp[i]] - celli[i];
      }
    }
  }
#pragma omp parallel for num_threads(Thread_Num)
  for (int Id : InterfaceIdx) {
    BasicCell<T, LatSet> cell(Id, NS);
    NbrInfo nbrInfo = getNbrInfo(Id);
    T mass = Mass.get(Id);
    // prepare for curvature calculation
    T curvature = T(0);
    if (has_surface_tension) {
      if (nbrInfo.hasGasNbr) curvature = calculateSurfaceTensionCurvature2D(Id);
    }
    T gas_pressure = 1. - 6. * surface_tension_parameter * curvature;
    for (int i = 1; i < LatSet::q; ++i) {
      int nbrId = Geo.template getNeighborId<LatSet>(Id, i);
      const BasicCell<T, LatSet> celli(nbrId, NS);
      int iopp = LatSet::opp[i];
      if (isType(nbrId, FSType::Fluid)) {
        mass += cell[iopp] - celli[i];
      } else if (isType(nbrId, FSType::Interface)) {
        T massflow = cell[iopp] - celli[i];
        mass += massflow * T(0.5) * (getF(Id) + getF(nbrId));
      } else if (isType(nbrId, FSType::Gas)) {
        // reconstruct interface pops streamed from gas cells
        Vector<T, LatSet::d> u = PREVIOUS_VELOCITY.get(Id);
        T u2 = u.getnorm2();
        cell[iopp] = Equilibrium<T, LatSet>::Order2(i, u, gas_pressure, u2) +
                     Equilibrium<T, LatSet>::Order2(iopp, u, gas_pressure, u2) -
                     celli[i];
      }
    }
    Mass.SetField(Id, mass);
    // set conversion flag
    T rho = moment::Rho<T, LatSet>::get(cell);
    if (mass < -transitionThreshold * rho ||
        (mass < lonely_threshold * rho && !nbrInfo.hasFluidNbr)) {
      setFlag(Id, Flag::ToGas);
    } else if (mass > (1. + transitionThreshold) * rho ||
               (mass > (1 - lonely_threshold) * rho && !nbrInfo.hasGasNbr)) {
      setFlag(Id, Flag::ToFluid);
    } else if (drop_isolated_cells && (nbrInfo.InterfaceNbrNum == 0)) {
      if (!nbrInfo.hasGasNbr) {
        setFlag(Id, Flag::ToFluid);
      }
    }
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::ToFluidCellConversion() {
/*
 * Initializing new interface cells with DFs from surrounding fluid and
 * interface cells Just takes the arithmetic average.
 */
#pragma omp parallel for num_threads(Thread_Num)
  for (int Id = 0; Id < N; ++Id) {
    if (isType(Id, FSType::Gas)) {
      if (hasNeighborFlag(Id, Flag::ToFluid)) {
        BasicCell<T, LatSet> cell(Id, NS);
        setFlag(Id, Flag::NewInterface);
        T rho_avg = 0.;
        Vector<T, LatSet::d> u_avg{};
        int ctr = 0;
        for (int i = 1; i < LatSet::q; i++) {
          int nbrId = Geo.template getNeighborId<LatSet>(Id, i);
          if (isType(nbrId, static_cast<FSType>(FSType::Fluid | FSType::Interface))) {
            T rho_tmp = 0.;
            Vector<T, LatSet::d> u_tmp{};
            ++ctr;
            BasicCell<T, LatSet> celli(nbrId, NS);
            moment::RhoVelocity<T, LatSet>::apply(celli, rho_tmp, u_tmp);
            rho_avg += rho_tmp;
            u_avg = u_tmp + u_avg;
          }
        }

        if (ctr > 0) {
          rho_avg /= static_cast<T>(ctr);
          u_avg = u_avg / static_cast<T>(ctr);
        }

        // initialize DFs with feq
        cell.template Compute_Equilibrium<
            Equilibrium<T, LatSet>::Feq_secondOrder>(rho_avg, u_avg);
        // end initialize DFs with feq
      }
    } else if (isFlag(Id, Flag::ToGas)) {
      // If a toGas cell has a neighbouring toFluid cell, unset the toGas flag
      if (hasNeighborFlag(Id, Flag::ToFluid)) setFlag(Id, Flag::None);
    }
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::ToGasCellConversion() {
  /*
   * For the to be converted toGas cells, set the neighbours to interface cells
   */
#pragma omp parallel for num_threads(Thread_Num)
  for (int Id = 0; Id < N; ++Id) {
    if (isType(Id, FSType::Fluid) && hasNeighborFlag(Id, Flag::ToGas)) {
      BasicCell<T, LatSet> cell(Id, NS);
      setFlag(Id, Flag::NewInterface);
      Mass.SetField(Id, moment::Rho<T, LatSet>::get(cell));
    }
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::MassExcess() {
#pragma omp parallel for num_threads(Thread_Num)
  for (int Id : InterfaceIdx) {
    BasicCell<T, LatSet> cell(Id, NS);
    T rho = moment::Rho<T, LatSet>::get(cell);
    T mass = Mass.get(Id);
    T mass_excess = 0.;

    Vector<T, LatSet::d> normal = computeParkerYoungInterfaceNormal(Id);
    // redistribute excess mass

    /// @hint EPSILON of neighbours used here
    /// @hint Mass can be set in this processor, but not epsilon since it is
    /// needed for the normal computation. epsilon is set in the next
    /// processor Became untrue due to code section removal, but epsilon is
    /// still set in the next part because of pushed mass excess
    if (isFlag(Id, Flag::ToGas)) {
      mass_excess = mass;
      Mass.SetField(Id, T(0));
      normal = {-normal[0], -normal[1]};
    } else if (isFlag(Id, Flag::ToFluid)) {
      mass_excess = mass - rho;
      Mass.SetField(Id, rho);
    } else {
      continue;
    }

    int product_total = 0;

    for (int i = 1; i < LatSet::q; ++i) {
      int nbrId = Geo.template getNeighborId<LatSet>(Id, i);
      if (isType(nbrId, FSType::Interface) && (isFlag(nbrId, Flag::None)))
        ++product_total;
    }

    /* Prepare Mass excess push */
    Vector<T, LatSet::q> mass_excess_vector{};
    mass_excess_vector[0] = 0.;
    if (product_total > 0) {
      T product_fraction = 1. / product_total;
      for (int i = 1; i < LatSet::q; ++i) {
        int nbrId = Geo.template getNeighborId<LatSet>(Id, i);
        if (isType(nbrId, FSType::Interface) && (isFlag(nbrId, Flag::None))) {
          mass_excess_vector[i] = mass_excess * product_fraction;
        } else {
          mass_excess_vector[i] = 0.;
        }
      }  // end for i
      TEMP_MASS_EXCHANGE.SetField(Id, mass_excess_vector);
    } else {
      mass_excess_vector[0] = mass_excess;
      for (int i = 1; i < LatSet::q; ++i) {
        mass_excess_vector[i] = 0.;
      }
      TEMP_MASS_EXCHANGE.SetField(Id, mass_excess_vector);
    }
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::FinalizeConversion() {
#pragma omp parallel for num_threads(Thread_Num)
  for (int Id = 0; Id < N; ++Id) {
    Flag flags = Flags.get(Id);

    switch (flags) {
      case Flag::ToFluid: {
        setType(Id, FSType::Fluid);
        Fv.SetField(Id, T(1));
        Mass.get(Id) += TEMP_MASS_EXCHANGE.get(Id)[0];
      } break;
      case Flag::ToGas: {
        setType(Id, FSType::Gas);
        Fv.SetField(Id, T(0));
        Mass.get(Id) += TEMP_MASS_EXCHANGE.get(Id)[0];
      } break;
      case Flag::NewInterface: {
        setType(Id, FSType::Interface);
      } break;
      default:
        break;
    }  // end switch flags

    FSType type = Types.get(Id);

    /* Collection of mass excess in a pulling step */
    switch (type) {
      case FSType::Interface: {
        T collected_excess = 0.;
        for (int i = 1; i < LatSet::q; ++i) {
          int nbrId = Geo.template getNeighborId<LatSet>(Id, i);
          const auto &tempMassExchange = TEMP_MASS_EXCHANGE.get(nbrId);

          if (isFlag(nbrId, static_cast<Flag>(Flag::ToFluid | Flag::ToGas))) {
            int iopp = LatSet::opp[i];
            collected_excess += tempMassExchange[iopp];
          }
        }

        T mass_tmp = Mass.get(Id);

        mass_tmp += collected_excess;
        T rho;
        Vector<T, LatSet::d> u;
        BasicCell<T, LatSet> cell(Id, NS);
        moment::RhoVelocity<T, LatSet>::apply(cell, rho, u);

        Fv.SetField(Id, mass_tmp / rho);
        Mass.SetField(Id, mass_tmp);
        PREVIOUS_VELOCITY.SetField(Id, u);
      } break;
      case FSType::Fluid: {
        T collected_excess = 0.;

        for (int i = 1; i < LatSet::q; ++i) {
          int nbrId = Geo.template getNeighborId<LatSet>(Id, i);
          const auto &tempMassExchange = TEMP_MASS_EXCHANGE.get(nbrId);

          if (isFlag(nbrId, static_cast<Flag>(Flag::ToFluid | Flag::ToGas))) {
            int iopp = LatSet::opp[i];
            collected_excess += tempMassExchange[iopp];
          }
        }
        Mass.get(Id) += collected_excess;
      } break;
      default:
        break;
    }  // end switch type
    // Flags.SetField(Id, Flag::None);
  }
  // end for each voxel
}

///////////////////////////////////
// namespace FreeSurface
namespace FS {
template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::Collide(T omega, T _omega,
                                       const Vector<T, LatSet::d> &force) {
#pragma omp parallel for num_threads(Thread_Num)
  for (const Voxel<T, LatSet::d> &vox : Voxels) {
    int Id = vox.getId();
    if (isType(Id, static_cast<FSType>(FSType::Solid |
                                                  FSType::Gas)))
      continue;
    NS.getPop(Id).Compute_rhoU();
    NS.getPop(Id).BGK_O2(omega, _omega);
    NS.getPop(Id).addForceO2(NS.getPop(Id).rho * force, omega);
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::Stream() {
  // copy _Mass to mass_prev and rho_prev
  mass_prev = _Mass;
  rho_prev = _Mass;

#pragma omp parallel for num_threads(Thread_Num)
  for (const Voxel<T, LatSet::d> &vox : Voxels) {
    int Id = vox.getId();

    if (isType(Id, static_cast<FSType>(FSType::Solid |
                                                  FSType::Gas)))
      continue;

    // fluid
    if (isType(Id, FSType::Fluid)) {
      // local streaming
      NS.getPop(Id).f[0] = NS.getPop(Id).fpostcol[0];

      for (int i = 1; i < LatSet::q; ++i) {
        int nbrId = vox.getNeighborId(i);
        int iopp = LatSet::opp[i];
        // streaming
        if (isType(nbrId, static_cast<FSType>(
                              FSType::Fluid |
                              FSType::Interface))) {
          // streaming
          NS.getPop(Id).f[iopp] = NS.getPop(nbrId).fpostcol[iopp];
          // mass exchange
          _Mass[Id] +=
              NS.getPop(nbrId).fpostcol[iopp] - NS.getPop(Id).fpostcol[i];
        } else {
          // bounce back
          NS.getPop(Id).f[iopp] = NS.getPop(Id).fpostcol[i];
        }
      }
    } else if (isType(Id, FSType::Interface)) {
      // Interface
      // get fraction filled
      T epsilon = getEpsilon(_Type[Id], rho_prev[Id], mass_prev[Id]);
      // get equilibrium
      T gas_pressure = 1.;
      T feq[LatSet::q] = {T(0)};
      Equilibrium<T, LatSet>::Feq_secondOrder(feq, NS.getPop(Id).getVelocity(),
                                              gas_pressure);

      // local streaming
      NS.getPop(Id).f[0] = NS.getPop(Id).fpostcol[0];

      for (int i = 1; i < LatSet::q; ++i) {
        int nbrId = vox.getNeighborId(i);
        int iopp = LatSet::opp[i];
        // streaming
        // fluid cell
        if (isType(nbrId, FSType::Fluid)) {
          // streaming
          NS.getPop(Id).f[iopp] = NS.getPop(nbrId).fpostcol[iopp];
          // mass exchange
          _Mass[Id] +=
              NS.getPop(nbrId).fpostcol[iopp] - NS.getPop(Id).fpostcol[i];

          // interface cell
        } else if (isType(nbrId, FSType::Interface)) {
          // streaming
          NS.getPop(Id).f[iopp] = NS.getPop(nbrId).fpostcol[iopp];
          // mass exchange
          T epsilon_nei =
              getEpsilon(_Type[nbrId], rho_prev[nbrId], mass_prev[nbrId]);

          T massflow = T(0);
          if (isType(Id, FSType::NO_FLUID_NEIGH)) {
            if (isType(nbrId, FSType::NO_FLUID_NEIGH)) {
              massflow =
                  NS.getPop(nbrId).fpostcol[iopp] - NS.getPop(Id).fpostcol[i];
            } else {
              massflow = -NS.getPop(Id).fpostcol[i];
            }
          } else if (isType(Id, FSType::NO_EMPTY_NEIGH)) {
            if (isType(nbrId, FSType::NO_EMPTY_NEIGH)) {
              massflow =
                  NS.getPop(nbrId).fpostcol[iopp] - NS.getPop(Id).fpostcol[i];
            } else {
              massflow = NS.getPop(nbrId).fpostcol[iopp];
            }
          } else {
            if (isType(nbrId, FSType::NO_FLUID_NEIGH)) {
              massflow = NS.getPop(nbrId).fpostcol[iopp];
            } else if (isType(nbrId, FSType::NO_EMPTY_NEIGH)) {
              massflow = NS.getPop(Id).fpostcol[i];
            } else {
              massflow =
                  NS.getPop(nbrId).fpostcol[iopp] - NS.getPop(Id).fpostcol[i];
            }
          }  // end of massflow calculation
          _Mass[Id] += massflow * T(0.5) * (epsilon + epsilon_nei);

          // gas cell
        } else if (isType(nbrId, FSType::Gas)) {
          NS.getPop(Id).f[iopp] =
              feq[iopp] + feq[i] - NS.getPop(Id).fpostcol[i];
          // solid cell
        } else {
          // bounce back
          NS.getPop(Id).f[iopp] = NS.getPop(Id).fpostcol[i];
        }
      }  // end for i in LatSet::q

      // correct for surface normal
      Vector<T, LatSet::d> normal = get_normal(Id, _Type, rho_prev, mass_prev);
      for (int i = 1; i < LatSet::q; ++i) {
        int iopp = LatSet::opp[i];
        if (LatSet::c[iopp] * normal > 0) {
          NS.getPop(Id).f[i] =
              feq[i] + feq[iopp] - NS.getPop(Id).fpostcol[iopp];
        }
      }
      // end correct for surface normal
    }  // end if interface

    // calculate density and velocity
    NS.getPop(Id).Compute_rhoU();
    if (isType(Id, FSType::Fluid)) NS.getPop(Id).rho = _Mass[Id];
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::Conversion() {
  // copy
  for (const Voxel<T, LatSet::d> &vox : Voxels) {
    int Id = vox.getId();
    mass_prev[Id] = _Mass[Id];
    rho_prev[Id] = NS.getPop(Id).rho;
    u_prev[Id] = NS.getPop(Id).getVelocity();
    cell_type_prev[Id] = _Type[Id];
  }

  T fill_offset = 0.001;
  T lonely_tresh = 0.1;
  // set the to_fluid/gas flags
#pragma omp parallel for num_threads(Thread_Num)
  for (const Voxel<T, LatSet::d> &vox : Voxels) {
    int Id = vox.getId();
    if (isType(Id, FSType::Interface)) {
      T mass = _Mass[Id];
      T rho = NS.getPop(Id).rho;
      // T rho = 1;
      if (mass > (1 + fill_offset) * rho ||
          (mass >= (1 - lonely_tresh) * rho &&
           isType(Id, FSType::NO_EMPTY_NEIGH))) {
        setType(Id, FSType::TO_FLUID);
      } else if (mass < -fill_offset * rho ||
                 (mass <= lonely_tresh * rho &&
                  isType(Id, FSType::NO_FLUID_NEIGH)) ||
                 isType(Id, static_cast<FSType>(
                                FSType::NO_IFACE_NEIGH |
                                FSType::NO_FLUID_NEIGH))) {
        setType(Id, FSType::TO_GAS);
      }
    }
    // remove neighbor flags

    _Type[Id] = static_cast<FSType>(
        _Type[Id] & (~(FSType::NO_FLUID_NEIGH |
                       FSType::NO_EMPTY_NEIGH |
                       FSType::NO_IFACE_NEIGH)));
  }

// interface -> fluid
#pragma omp parallel for num_threads(Thread_Num)
  for (const Voxel<T, LatSet::d> &vox : Voxels) {
    int Id = vox.getId();
    if (isType(Id, FSType::TO_FLUID)) {
      for (int i = 1; i < LatSet::q; ++i) {
        int nbrId = vox.getNeighborId(i);
        if (isType(nbrId, FSType::Gas)) {
          setType(nbrId, FSType::Interface);

          T rho_avg = 0.;
          Vector<T, LatSet::d> u_avg{};

          const Voxel<T, LatSet::d> &nbrvox = Voxels[nbrId];

          int ctr = 0;
          for (int iPop = 1; iPop < LatSet::q; iPop++) {
            int nbrnbrId = nbrvox.getNeighborId(iPop);
            if (isType(nbrnbrId, static_cast<FSType>(
                                     FSType::Fluid |
                                     FSType::Interface |
                                     FSType::TO_FLUID))) {
              ++ctr;
              rho_avg += rho_prev[nbrnbrId];
              u_avg = u_avg + u_prev[nbrnbrId];
            }
          }
          if (ctr > 0) {
            rho_avg /= static_cast<T>(ctr);
            u_avg = u_avg / static_cast<T>(ctr);
          }

          NS.getPop(nbrId).rho = rho_avg;
          NS.getPop(nbrId).getVelocity() = u_avg;
          // if (u_avg[0] > 0.1) u_avg[0] = 0.1;
          // if (u_avg[1] > 0.1) u_avg[1] = 0.1;

          // initialize DFs with feq
          T u2 = u_avg.getnorm2();
          for (int i = 0; i < LatSet::q; ++i) {
            NS.getPop(nbrId).f[i] =
                Equilibrium<T, LatSet>::Order2(i, u_avg, rho_avg, u2);
          }
          // end initialize DFs with feq
        }
      }
    }
  }

// interface-> gas
#pragma omp parallel for num_threads(Thread_Num)
  for (const Voxel<T, LatSet::d> &vox : Voxels) {
    int Id = vox.getId();
    if (isType(Id, FSType::TO_GAS)) {
      for (int i = 1; i < LatSet::q; ++i) {
        int nbrId = vox.getNeighborId(i);
        if (isType(nbrId, FSType::Fluid)) {
          setType(nbrId, FSType::Interface);
        }
      }
    }
  }
// distribute excess mass
#pragma omp parallel for num_threads(Thread_Num)
  for (const Voxel<T, LatSet::d> &vox : Voxels) {
    int Id = vox.getId();
    if (isType(Id, FSType::Solid)) continue;
    Vector<T, LatSet::d> normal =
        get_normal(Id, cell_type_prev, rho_prev, mass_prev);
    T mex = 0.;
    if (isType(Id, FSType::TO_FLUID)) {
      mex = _Mass[Id] - NS.getPop(Id).rho;
      _Mass[Id] = NS.getPop(Id).rho;
    } else if (isType(Id, FSType::TO_GAS)) {
      mex = _Mass[Id];
      normal = {-normal[0], -normal[1]};
      _Mass[Id] = 0.;
    } else {
      continue;
    }

    std::array<T, LatSet::q> eta{0};
    std::array<bool, LatSet::q> isIF{false};
    T eta_total = 0.;
    T IF_total = 0.;

    for (int iPop = 1; iPop < LatSet::q; iPop++) {
      int nbrId = vox.getNeighborId(iPop);

      if (isType(nbrId, FSType::Interface)) {
        eta[iPop] = normal * LatSet::c[iPop];
        if (eta[iPop] <= 0) {
          eta[iPop] = 0.;
        }
        eta_total += eta[iPop];
        isIF[iPop] = true;
        IF_total += 1;
      }
    }  // end for iPop

    if (eta_total > 0) {
      T eta_frac = T(1) / eta_total;
      for (int iPop = 1; iPop < LatSet::q; iPop++) {
        TEMP_MASS_EXCHANGE[Id][iPop] = mex * eta[iPop] * eta_frac;
      }
    } else if (IF_total > 0) {
      T mex_rel = mex / IF_total;
      for (int iPop = 1; iPop < LatSet::q; iPop++) {
        if (isIF[iPop]) {
          TEMP_MASS_EXCHANGE[Id][iPop] = mex_rel;
        } else {
          TEMP_MASS_EXCHANGE[Id][iPop] = 0.;
        }
      }
    }
  }

// collect distributed mass and finalize cell flags
#pragma omp parallel for num_threads(Thread_Num)
  for (const Voxel<T, LatSet::d> &vox : Voxels) {
    int Id = vox.getId();
    if (isType(Id, FSType::Interface)) {
      for (int iPop = 1; iPop < LatSet::q; iPop++) {
        int nbrId = vox.getNeighborId(iPop);
        _Mass[Id] += TEMP_MASS_EXCHANGE[nbrId][iPop];
      }
    } else if (isType(Id, FSType::TO_FLUID)) {
      setType(Id, FSType::Fluid);
    } else if (isType(Id, FSType::TO_GAS)) {
      setType(Id, FSType::Gas);
    }
  }

// set neighborhood flags
#pragma omp parallel for num_threads(Thread_Num)
  for (const Voxel<T, LatSet::d> &vox : Voxels) {
    int Id = vox.getId();
    if (isType(Id, FSType::Solid)) continue;

    _Type[Id] = static_cast<FSType>(
        _Type[Id] |
        (FSType::NO_FLUID_NEIGH | FSType::NO_EMPTY_NEIGH |
         FSType::NO_IFACE_NEIGH));
    for (int iPop = 1; iPop < LatSet::q; iPop++) {
      int nbrId = vox.getNeighborId(iPop);
      if (isType(nbrId, FSType::Fluid)) {
        _Type[Id] = static_cast<FSType>(
            _Type[Id] & (~FSType::NO_FLUID_NEIGH));
      } else if (isType(nbrId, FSType::Interface)) {
        _Type[Id] = static_cast<FSType>(
            _Type[Id] & (~FSType::NO_IFACE_NEIGH));
      } else if (isType(nbrId, FSType::Gas)) {
        _Type[Id] = static_cast<FSType>(
            _Type[Id] & (~FSType::NO_EMPTY_NEIGH));
      }
    }
  }

  // // ToFluidCell
  // #pragma omp parallel for num_threads(Thread_Num)
  //   for (const Voxel<T, LatSet::d> &vox : Voxels) {
  //     int Id = vox.getId();
  //     if (isType(Id, FSType::Gas)) {
  //       if (hasNeighborFlag(Id, Flag::ToFluid)) {
  //         setFlag(Id, Flag::NewInterface);
  //         T rho_avg = 0.;
  //         Vector<T, LatSet::d> u_avg{};

  //         int ctr = 0;
  //         for (int iPop = 1; iPop < LatSet::q; iPop++) {
  //           int nbrId = vox.getNeighborId(iPop);
  //           if (isType(nbrId, FSType::Fluid) || isType(nbrId,
  //           FSType::Interface))
  //           {
  //             T rho_tmp = 0.;
  //             Vector<T, LatSet::d> u_tmp{};
  //             ++ctr;
  //             NS.getPop(nbrId).Compute_rhoU(rho_tmp, u_tmp);
  //             rho_avg += rho_tmp;
  //             u_avg = u_tmp + u_avg;
  //           }
  //         }
  //         if (ctr > 0) {
  //           rho_avg /= static_cast<T>(ctr);
  //           u_avg = u_avg / static_cast<T>(ctr);
  //         }

  //         // initialize DFs with feq
  //         T u2 = u_avg.getnorm2();
  //         for (int i = 0; i < LatSet::q; ++i) {
  //           NS.getPop(Id).f[i] =
  //               Equilibrium<T, LatSet>::Order2(i, u_avg, rho_avg, u2);
  //         }
  //         // end initialize DFs with feq
  //       }
  //     } else if (isFlag(Id, Flag::ToGas)) {
  //       // If a toGas cell has a neighbouring toFluid cell, unset the toGas
  //       flag if (hasNeighborFlag(Id, Flag::ToFluid)) setFlag(Id,
  //       Flag::None);
  //     }
  //   }
  // // end ToFluidCell

  // // ToGasCell
  // #pragma omp parallel for num_threads(Thread_Num)
  //   for (const Voxel<T, LatSet::d> &vox : Voxels) {
  //     int Id = vox.getId();
  //     if (isType(Id, FSType::Fluid) && hasNeighborFlag(Id, Flag::ToGas)) {
  //       setFlag(Id, Flag::NewInterface);
  //       T rho = T(0);
  //       NS.getPop(Id).Compute_rho(rho);
  //       _Mass[Id] = rho;
  //     }
  //   }
  // // end ToGasCell

  // // MassExcess
  // #pragma omp parallel for num_threads(Thread_Num)
  //   for (int Id : InterfaceIdx) {
  //     const Voxel<T, LatSet::d> &vox = Voxels[Id];
  //     T rho = T(0);
  //     NS.getPop(Id).Compute_rho(rho);
  //     T mass = _Mass[Id];
  //     T mass_excess = 0.;
  //     Vector<T, LatSet::d> normal = getNormal(Id);

  //     if (isFlag(Id, Flag::ToGas)) {
  //       mass_excess = mass;
  //       _Mass[Id] = 0.;
  //       normal = {-normal[0], -normal[1]};
  //     } else if (isFlag(Id, Flag::ToFluid)) {
  //       mass_excess = mass - rho;
  //       _Mass[Id] = rho;
  //     } else {
  //       continue;
  //     }

  //     std::array<T, LatSet::q> eta{0};
  //     std::array<bool, LatSet::q> isIF{false};
  //     T eta_total = 0.;
  //     T IF_total = 0.;

  //     for (int iPop = 1; iPop < LatSet::q; iPop++) {
  //       int nbrId = vox.getNeighborId(iPop);

  //       if (isType(nbrId, FSType::Interface)) {
  //         eta[iPop] = normal * LatSet::c[iPop];
  //         if (eta[iPop] <= 0) {
  //           eta[iPop] = 0.;
  //         }
  //         eta_total += eta[iPop];
  //         isIF[iPop] = true;
  //         IF_total += 1;
  //       }
  //     }

  //     if (eta_total > 0) {
  //       T eta_frac = T(1) / eta_total;
  //       for (int iPop = 1; iPop < LatSet::q; iPop++) {
  //         TEMP_MASS_EXCHANGE[Id][iPop] = mass_excess * eta[iPop] *
  //         eta_frac;
  //       }
  //     } else if (IF_total > 0) {
  //       T mex_rel = mass_excess / IF_total;
  //       for (int iPop = 1; iPop < LatSet::q; iPop++) {
  //         if (isIF[iPop]) {
  //           TEMP_MASS_EXCHANGE[Id][iPop] = mex_rel;
  //         } else {
  //           TEMP_MASS_EXCHANGE[Id][iPop] = 0.;
  //         }
  //       }
  //     }
  //   }
  // // end MassExcess
  // // collect distributed mass and finalize cell flags
  // #pragma omp parallel for num_threads(Thread_Num) schedule(static)
  //   for (int Id : InterfaceIdx) {
  //     const Voxel<T, LatSet::d> &vox = Voxels[Id];
  //     for (int iPop = 1; iPop < LatSet::q; iPop++) {
  //       int nbrId = vox.getNeighborId(iPop);
  //       const auto &tempMassExchange = TEMP_MASS_EXCHANGE[nbrId];
  //       _Mass[Id] += tempMassExchange[iPop];
  //     }

  //     T rho = T(0);
  //     Vector<T, LatSet::d> u{};
  //     NS.getPop(Id).Compute_rhoU(rho, u);

  //     _F[Id] = _Mass[Id] / rho;
  //     PREVIOUS_VELOCITY[Id] = u;

  //     if (isFlag(Id, Flag::ToFluid)) {
  //       setType(Id, FSType::Fluid);
  //       _F[Id] = T(1);
  //     } else if (isFlag(Id, Flag::ToGas)) {
  //       setType(Id, FSType::Gas);
  //       _F[Id] = T(0);
  //     } else if (isFlag(Id, Flag::NewInterface)) {
  //       setType(Id, FSType::Interface);
  //     }
  //   }

  // end collect distributed mass and finalize cell flags

  // FinalizeConversion
  // #pragma omp parallel for num_threads(Thread_Num) schedule(static)
  //   for (const Voxel<T, LatSet::d> &vox : Voxels) {
  //     int Id = vox.getId();
  //     Flag flags = _Flag[Id];

  //     switch (flags) {
  //       case Flag::ToFluid: {
  //         setType(Id, FSType::Fluid);
  //         _F[Id] = T(1);
  //         _Mass[Id] += TEMP_MASS_EXCHANGE[Id][0];
  //       } break;
  //       case Flag::ToGas: {
  //         setType(Id, FSType::Gas);
  //         _F[Id] = T(0);
  //         _Mass[Id] += TEMP_MASS_EXCHANGE[Id][0];
  //       } break;
  //       case Flag::NewInterface: {
  //         setType(Id, FSType::Interface);
  //       } break;
  //       default:
  //         break;
  //     }  // end switch flags

  //     FSType type = _Type[Id];
  //     /* Collection of mass excess in a pulling step */
  //     switch (type) {
  //       case FSType::Interface: {
  //         T collected_excess = 0.;
  //         for (int iPop = 1; iPop < LatSet::q; ++iPop) {
  //           int nbrId = vox.getNeighborId(iPop);
  //           const auto &tempMassExchange = TEMP_MASS_EXCHANGE[nbrId];

  //           if (isFlag(nbrId, Flag::ToFluid) || isFlag(nbrId, Flag::ToGas))
  //           {
  //             int iPop_op = LatSet::opp[iPop];
  //             collected_excess += tempMassExchange[iPop_op];
  //           }
  //         }

  //         T mass_tmp = _Mass[Id];
  //         mass_tmp += collected_excess;

  //         T rho;
  //         Vector<T, LatSet::d> u;
  //         NS.getPop(Id).Compute_rhoU(rho, u);

  //         _F[Id] = mass_tmp / rho;
  //         _Mass[Id] = mass_tmp;
  //         PREVIOUS_VELOCITY[Id] = u;

  //       } break;
  //       case FSType::Fluid: {
  //         T collected_excess = 0.;

  //         for (int iPop = 1; iPop < LatSet::q; ++iPop) {
  //           int nbrId = vox.getNeighborId(iPop);
  //           const auto &tempMassExchange = TEMP_MASS_EXCHANGE[nbrId];
  //           if (isFlag(nbrId, Flag::ToFluid) || isFlag(nbrId, Flag::ToGas))
  //           {
  //             int iPop_op = LatSet::opp[iPop];
  //             collected_excess += tempMassExchange[iPop_op];
  //           }
  //         }
  //         _Mass[Id] += collected_excess;
  //       } break;
  //       default:
  //         break;
  //     }  // end switch type
  //   }
  // end FinalizeConversion
}

template <typename T, typename LatSet>
Vector<T, LatSet::d> FreeSurface2D<T, LatSet>::get_normal(
    int id, const std::vector<FSType> &types,
    const std::vector<T> &rho_prev, const std::vector<T> &mass_prev) {
  int id1 = Voxels[id].getNeighborId(1);
  int id2 = Voxels[id].getNeighborId(2);
  int id3 = Voxels[id].getNeighborId(3);
  int id4 = Voxels[id].getNeighborId(4);
  return T(0.5) *
         Vector<T, LatSet::d>(
             getEpsilon(types[id3], rho_prev[id3], mass_prev[id3]) -
                 getEpsilon(types[id1], rho_prev[id1], mass_prev[id1]),
             getEpsilon(types[id4], rho_prev[id4], mass_prev[id4]) -
                 getEpsilon(types[id2], rho_prev[id2], mass_prev[id2]));
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::Initialize() {
  for (const Voxel<T, LatSet::d> &vox : Voxels) {
    int Id = vox.getId();
    _Mass[Id] = NS.getPoprho(Id);
    // if (isType(Id, FSType::Interface)) {
    //   _F[Id] = T(0.5);
    // } else if (isType(Id, FSType::Gas)) {
    //   _F[Id] = T(0);
    // } else {
    //   _F[Id] = T(1);
    // }
  }
}

template <typename T, typename LatSet>
bool FreeSurface2D<T, LatSet>::hasNeighborType(int id,
                                               const FSType &type) {
  const Voxel<T, 2> &voxel = Voxels[id];
  for (int i = 1; i < LatSet::q; ++i) {
    if (voxel.getNeighbor(i) != nullptr) {
      int nbrId = voxel.getNeighborId(i);
      if (isType(nbrId, type)) return true;
    }
  }
  return false;
}

template <typename T, typename LatSet>
bool FreeSurface2D<T, LatSet>::isType(int id, const FSType &type) {
  return static_cast<bool>(_Type[id] & type);
}
template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::setType(int id, const FSType &type) {
  _Type[id] = type;
}

}  // namespace FreeSurface
// namespace FreeSurface
