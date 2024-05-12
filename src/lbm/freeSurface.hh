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

#include "freeSurface.h"
#include "lbm/freeSurface.h"

///////////////////////////////////
// namespace FreeSurface
namespace FS {

template <typename T, typename LatSet>
inline bool FreeSurface2D<T, LatSet>::hasNeighborType(std::size_t id,
                                                      FSType fstype) const {
  for (int i = 1; i < LatSet::q; ++i) {
    if (util::isFlag(State.get(id + NS.getDelta_Index()[i]), fstype)) return true;
  }
  return false;
}

template <typename T, typename LatSet>
T ComputeCurvature(BCell<T, LatSet>& cell) {
  return T{};
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::MassTransfer() {
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      std::size_t id = i + j * NS.getNx();
      // for interface cells
      if (util::isFlag(State.get(id), FSType::Interface)) {
        T deltamass = T{};
        BCell<T, LatSet> cell(id, NS);
        // find neighbor cells
        for (int k = 1; k < LatSet::q; ++k) {
          std::size_t idn = id + NS.getDelta_Index()[k];
          const BCell<T, LatSet> celln(idn, NS);
          if (util::isFlag(State.get(idn), FSType::Fluid)) {
            deltamass += cell[LatSet::opp[k]] - celln[k];
          } else if (util::isFlag(State.get(idn), FSType::Interface)) {
            deltamass += (cell[LatSet::opp[k]] - celln[k]) * T(0.5) *
                         (getClampedVOF(id) + getClampedVOF(idn));
          }
        }
        Mass.get(id) += deltamass;

        // reconstruct pop streamed in from a gas cell
        T curvature = T{};
        if (Surface_Tension_Enabled) {
          if (hasNeighborType(id, FSType::Gas)) curvature = ComputeCurvature(cell);
        }
        T rho_gas = T(1) - T(6) * surface_tension_parameter * curvature;
        const Vector<T, LatSet::d>& u = cell.getVelocity();
        T u2 = u.getnorm2();
        for (int k = 1; k < LatSet::q; ++k) {
          std::size_t idn = id + NS.getDelta_Index()[k];
          if (util::isFlag(State.get(idn), FSType::Gas)) {
            const BCell<T, LatSet> celln(idn, NS);
            // fiopp = feqiopp(rho_gas) + feqi(rho_gas) - fi(x+ei)
            cell[LatSet::opp[k]] =
              Equilibrium<T, LatSet>::Order2(k, u, rho_gas, u2) +
              Equilibrium<T, LatSet>::Order2(LatSet::opp[k], u, rho_gas, u2) - celln[k];
          }
        }

        // transition flag for interface cell
        T rho = moment::Rho<T, LatSet>::get(cell);

        // transition by mass criterion
        if (Mass.get(id) > (T(1) + VOF_Trans_Threshold) * rho) {
          util::addFlag(FSType::To_Fluid, State.getField(0).getUnderlying(id));
          continue;
        } else if (Mass.get(id) < -VOF_Trans_Threshold * rho) {
          util::addFlag(FSType::To_Gas, State.getField(0).getUnderlying(id));
          continue;
        }
        // transition by lonely criterion
        if (Mass.get(id) > (T(1) - Lonely_Threshold) * rho) {
          if (!hasNeighborType(id, FSType::Gas)) {
            util::addFlag(FSType::To_Fluid, State.getField(0).getUnderlying(id));
            continue;
          }
        } else if (Mass.get(id) < Lonely_Threshold * rho) {
          if (!hasNeighborType(id, FSType::Fluid)) {
            util::addFlag(FSType::To_Gas, State.getField(0).getUnderlying(id));
            continue;
          }
        }
        // deal with isolated interface cells
        if (!hasNeighborType(id, FSType::Interface)) {
          if (!hasNeighborType(id, FSType::Fluid)) {
            util::addFlag(FSType::To_Gas, State.getField(0).getUnderlying(id));
          } else if (!hasNeighborType(id, FSType::Gas)) {
            util::addFlag(FSType::To_Fluid, State.getField(0).getUnderlying(id));
          }
        }
      }
    }
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::ToFluidNbrConversion() {
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();
      if (util::isFlag(State.get(id), FSType::To_Fluid)) {
        // check neighbors
        for (int k = 1; k < LatSet::q; ++k) {
          const std::size_t idn = id + NS.getDelta_Index()[k];
          if (util::isFlag(State.get(idn), FSType::To_Gas)) {
            // remove to_gas flag
            util::removeFlag(FSType::To_Gas, State.getField(0).getUnderlying(idn));
          } else if (util::isFlag(State.get(idn), FSType::Gas)) {
            // set to_interface for gas neighbor cells
            util::addFlag(FSType::To_Interface, State.getField(0).getUnderlying(idn));
          }
        }
      }
    }
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::GasToInterfacePopInit() {
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();
      if (util::isFlag(State.get(id), FSType::To_Interface)) {
        BCell<T, LatSet> cell(id, NS);
        // init fi using [fluid|interface] neighbor cells
        T averho = T{};
        Vector<T, LatSet::d> aveu = Vector<T, LatSet::d>{};
        int count = 0;
        for (int k = 1; k < LatSet::q; ++k) {
          const std::size_t idn = id + NS.getDelta_Index()[k];
          const BCell<T, LatSet> celln(idn, NS);
          if (util::isFlag(State.get(idn), (FSType::Fluid | FSType::Interface))) {
            averho += moment::Rho<T, LatSet>::get(celln);
            aveu += moment::Velocity<T, LatSet>::get(celln);
            ++count;
          }
        }
        averho /= count;
        aveu /= count;
        // set fi
        T aveu2 = aveu.getnorm2();
        for (int k = 0; k < LatSet::q; ++k) {
          cell[k] = Equilibrium<T, LatSet>::Order2(k, aveu, averho, aveu2);
        }
      }
    }
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::ToGasNbrConversion() {
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();
      if (util::isFlag(State.get(id), FSType::To_Gas)) {
        // check neighbors
        for (int k = 1; k < LatSet::q; ++k) {
          const std::size_t idn = id + NS.getDelta_Index()[k];
          if (util::isFlag(State.get(idn), FSType::Fluid)) {
            // to interface
            util::addFlag(FSType::To_Interface, State.getField(0).getUnderlying(idn));
          }
        }
      }
    }
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::InterfaceExcessMass() {
  util::Parker_YoungsNormal2D<T, LatSet> PYNormal(NS.getNx(), NS.getNy(),
                                                  VolumeFrac.getField(0));
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();

      // for interface cells to be converted to fluid or gas
      // excess mass is distributed to interface neighbors
      T excessmass{};
      if (util::isFlag(State.get(id), FSType::To_Fluid)) {
        BCell<T, LatSet> cell(id, NS);
        T rho = moment::Rho<T, LatSet>::get(cell);
        excessmass = Mass.get(id) - rho;
        Mass.get(id) = rho;
      } else if (util::isFlag(State.get(id), FSType::To_Gas)) {
        excessmass = Mass.get(id);
        Mass.get(id) = T{};
      } else {
        continue;
      }

      // find neighbors
      int count{};
      for (int k = 1; k < LatSet::q; ++k) {
        const std::size_t idn = id + NS.getDelta_Index()[k];
        if (util::isFlag(State.get(idn), FSType::Interface) &&
            !util::isFlag(State.get(idn), FSType::To_Gas | FSType::To_Fluid)) {
          ++count;
        }
      }

      std::array<T*, LatSet::q> exmasscell = ExcessMass.getArray(id);
      *(exmasscell[0]) = T{};
      if (count > 0) {
        T excessmassk = excessmass / count;
        for (int k = 1; k < LatSet::q; ++k) {
          const std::size_t idn = id + NS.getDelta_Index()[k];
          if (util::isFlag(State.get(idn), FSType::Interface) &&
              !util::isFlag(State.get(idn), FSType::To_Gas | FSType::To_Fluid)) {
            *(exmasscell[k]) = excessmassk;
          }
        }
      } else {
        *(exmasscell[0]) = excessmass;
      }
    }
  }
}


template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::FinalizeConversion() {
  // update state
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();
      BCell<T, LatSet> cell(id, NS);
      if (util::isFlag(State.get(id), FSType::To_Fluid)) {
        State.SetField(id, FSType::Fluid);
        VolumeFrac.SetField(id, T(1));
        Mass.get(id) += ExcessMass.template get<0>(id);
      } else if (util::isFlag(State.get(id), FSType::To_Gas)) {
        State.SetField(id, FSType::Gas);
        VolumeFrac.SetField(id, T{});
        Mass.get(id) += ExcessMass.template get<0>(id);
        cell.getVelocity().clear();
      } else if (util::isFlag(State.get(id), FSType::To_Interface)) {
        State.SetField(id, FSType::Interface);
      }
    }
  }
}

template <typename T, typename LatSet>
void FreeSurface2D<T, LatSet>::CollectExcessMass() {
  // stream
  for (int i = 1; i < LatSet::q; ++i) {
    ExcessMass.getField(i).rotate(NS.getDelta_Index()[i]);
  }
  // collect
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();
      if (util::isFlag(State.get(id), FSType::Interface | FSType::Fluid)) {
                T exmass_sum = ExcessMass.template get<0>(id);
        for (int k = 1; k < LatSet::q; ++k) {
          exmass_sum += ExcessMass.get(id, k);
        }
        Mass.get(id) += exmass_sum;

        if (util::isFlag(State.get(id), FSType::Interface)) {
          BCell<T, LatSet> cell(id, NS);
          T rho = moment::Rho<T, LatSet>::get(cell);
          VolumeFrac.SetField(id, Mass.get(id) / rho);
        }
      }
    }
  }
  // clear
  ExcessMass.Init();
}

}  // namespace FS
// namespace FreeSurface
