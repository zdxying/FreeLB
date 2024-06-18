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

template <typename T, typename LatSet, typename TypePack>
void FreeSurface2D<T, LatSet, TypePack>::MassTransfer() {
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      std::size_t id = i + j * NS.getNx();
      // for interface cells
      if (util::isFlag(this->template getField<STATE>().get(id), FSType::Interface)) {
        T deltamass{};
        T massflow{};
        // cell's nbr info
        Cell<T, LatSet, TypePack> cell(id, NS);
        NbrInfo cellNbrInfo;
        getNbrInfo(cell.getId(), cellNbrInfo);
        // find neighbor cells
        for (int k = 1; k < LatSet::q; ++k) {
          std::size_t idn = id + this->Delta_Index[k];
          Cell<T, LatSet, TypePack> celln(idn, NS);
          int kopp = LatSet::opp[k];
          if (util::isFlag(this->template getField<STATE>().get(idn), FSType::Fluid)) {
            deltamass += cell[kopp] - celln[k];
          } else if (util::isFlag(this->template getField<STATE>().get(idn), FSType::Interface)) {
          // celln's nbr info
          NbrInfo cellnNbrInfo;
          getNbrInfo(celln.getId(), cellnNbrInfo);
          if (!cellNbrInfo.fluid_nbr){
            if (!cellnNbrInfo.fluid_nbr){
              if (cellNbrInfo.interface_nbrs < cellnNbrInfo.interface_nbrs){
                massflow = -celln[k];
              } else if (cellNbrInfo.interface_nbrs > cellnNbrInfo.interface_nbrs){
                massflow = cell[kopp];
              } else {
                massflow = cell[kopp] - celln[k];
              }
            } else {
              massflow = -celln[k];
            }
          } else if (!cellNbrInfo.gas_nbr){
            if (!cellnNbrInfo.gas_nbr){
              if (cellNbrInfo.interface_nbrs < cellnNbrInfo.interface_nbrs){
                massflow = cell[kopp];
              } else if (cellNbrInfo.interface_nbrs > cellnNbrInfo.interface_nbrs){
                massflow = -celln[k];
              } else {
                massflow = cell[kopp] - celln[k];
              }
            } else {
              massflow = cell[kopp];
            }
          } else {
            if (!cellnNbrInfo.fluid_nbr){
              massflow = cell[kopp];
            } else if (!cellnNbrInfo.gas_nbr){
              massflow = -celln[k];
            } else {
              massflow = cell[kopp] - celln[k];
            } 
          }

          deltamass += massflow * T(0.5) * (getClampedVOF(cell.getId()) + getClampedVOF(celln.getId()));
          }
        }
        this->template getField<MASS<T>>().get(id) += deltamass;

        // reconstruct pop streamed in from a gas cell
        T curvature{};
        // if (Surface_Tension_Enabled) {
        //   if (hasNeighborType(id, FSType::Gas)) curvature = ComputeCurvature(cell);
        // }
        T rho_gas = T(1) - T(6) * surface_tension_parameter * curvature;
        const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T,LatSet::d>>();
        T u2 = u.getnorm2();
        for (int k = 1; k < LatSet::q; ++k) {
          std::size_t idn = id + this->Delta_Index[k];
          if (util::isFlag(this->template getField<STATE>().get(idn), FSType::Gas)) {
            Cell<T, LatSet, TypePack> celln(idn, NS);
            // fiopp = feqiopp(rho_gas) + feqi(rho_gas) - fi(x+ei)
            cell[LatSet::opp[k]] =
              Equilibrium<T, LatSet>::Order2(k, u, rho_gas, u2) +
              Equilibrium<T, LatSet>::Order2(LatSet::opp[k], u, rho_gas, u2) - celln[k];
          }
        }

        // transition flag for interface cell
        T rho = moment::rho<Cell<T, LatSet, TypePack>>::get(cell);

        // transition by mass criterion
        if (this->template getField<MASS<T>>().get(id) > (T(1) + VOF_Trans_Threshold) * rho) {
          util::addFlag(FSType::To_Fluid, util::underlyingRef(this->template getField<STATE>().get(id)));
          continue;
        } else if (this->template getField<MASS<T>>().get(id) < -VOF_Trans_Threshold * rho) {
          util::addFlag(FSType::To_Gas, util::underlyingRef(this->template getField<STATE>().get(id)));
          continue;
        }
        // transition by lonely criterion
        if (this->template getField<MASS<T>>().get(id) > (T(1) - Lonely_Threshold) * rho) {
          if (!hasNeighborType(id, FSType::Gas)) {
            util::addFlag(FSType::To_Fluid, util::underlyingRef(this->template getField<STATE>().get(id)));
            continue;
          }
        } else if (this->template getField<MASS<T>>().get(id) < Lonely_Threshold * rho) {
          if (!hasNeighborType(id, FSType::Fluid)) {
            util::addFlag(FSType::To_Gas, util::underlyingRef(this->template getField<STATE>().get(id)));
            continue;
          }
        }
        // deal with isolated interface cells
        if (!hasNeighborType(id, FSType::Interface)) {
          if (!hasNeighborType(id, FSType::Fluid)) {
            util::addFlag(FSType::To_Gas, util::underlyingRef(this->template getField<STATE>().get(id)));
          } else if (!hasNeighborType(id, FSType::Gas)) {
            util::addFlag(FSType::To_Fluid, util::underlyingRef(this->template getField<STATE>().get(id)));
          }
        }
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
void FreeSurface2D<T, LatSet, TypePack>::ToFluidNbrConversion() {
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();
      if (util::isFlag(this->template getField<STATE>().get(id), FSType::To_Fluid)) {
        // check neighbors
        for (int k = 1; k < LatSet::q; ++k) {
          const std::size_t idn = id + this->Delta_Index[k];
          if (util::isFlag(this->template getField<STATE>().get(idn), FSType::To_Gas)) {
            // remove to_gas flag
            util::removeFlag(FSType::To_Gas, util::underlyingRef(this->template getField<STATE>().get(idn)));
          } else if (util::isFlag(this->template getField<STATE>().get(idn), FSType::Gas)) {
            // set to_interface for gas neighbor cells
            util::addFlag(FSType::To_Interface, util::underlyingRef(this->template getField<STATE>().get(idn)));
          }
        }
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
void FreeSurface2D<T, LatSet, TypePack>::GasToInterfacePopInit() {
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();
      if (util::isFlag(this->template getField<STATE>().get(id), FSType::To_Interface)) {
        Cell<T, LatSet, TypePack> cell(id, NS);
        // init fi using [fluid|interface] neighbor cells
        T averho = T{};
        Vector<T, LatSet::d> aveu = Vector<T, LatSet::d>{};
        int count = 0;
        for (int k = 1; k < LatSet::q; ++k) {
          const std::size_t idn = id + this->Delta_Index[k];
          Cell<T, LatSet, TypePack> celln(idn, NS);
          if (util::isFlag(this->template getField<STATE>().get(idn), (FSType::Fluid | FSType::Interface))) {
            averho += moment::rho<Cell<T, LatSet, TypePack>>::get(celln);
            aveu += moment::u<Cell<T, LatSet, TypePack>>::get(celln);
            ++count;
          }
        }
        averho /= count;
        aveu /= count;
        // set fi
        T aveu2 = aveu.getnorm2();
        for (int k = 0; k < LatSet::q; ++k) {
          cell[k] = equilibrium::SecondOrder<Cell<T, LatSet, TypePack>>::get(k, aveu, averho, aveu2);
        }
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
void FreeSurface2D<T, LatSet, TypePack>::ToGasNbrConversion() {
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();
      if (util::isFlag(this->template getField<STATE>().get(id), FSType::To_Gas)) {
        // check neighbors
        for (int k = 1; k < LatSet::q; ++k) {
          const std::size_t idn = id + this->Delta_Index[k];
          if (util::isFlag(this->template getField<STATE>().get(idn), FSType::Fluid)) {
            // to interface
            util::addFlag(FSType::To_Interface, util::underlyingRef(this->template getField<STATE>().get(idn)));
          }
        }
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
void FreeSurface2D<T, LatSet, TypePack>::InterfaceExcessMass() {
  // util::Parker_YoungsNormal2D<T, LatSet> PYNormal(NS.getNx(), NS.getNy(),
  //                                                 this->template getField<VOLUMEFRAC<T>>().getField(0));
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();

      // for interface cells to be converted to fluid or gas
      // excess mass is distributed to interface neighbors
      T excessmass{};
      if (util::isFlag(this->template getField<STATE>().get(id), FSType::To_Fluid)) {
        Cell<T, LatSet, TypePack> cell(id, NS);
        T rho = moment::rho<Cell<T, LatSet, TypePack>>::get(cell);
        excessmass = this->template getField<MASS<T>>().get(id) - rho;
        this->template getField<MASS<T>>().get(id) = rho;
      } else if (util::isFlag(this->template getField<STATE>().get(id), FSType::To_Gas)) {
        excessmass = this->template getField<MASS<T>>().get(id);
        this->template getField<MASS<T>>().get(id) = T{};
      } else {
        continue;
      }

      // find neighbors
      int count{};
      for (int k = 1; k < LatSet::q; ++k) {
        const std::size_t idn = id + this->Delta_Index[k];
        if (util::isFlag(this->template getField<STATE>().get(idn), FSType::Interface) &&
            !util::isFlag(this->template getField<STATE>().get(idn), FSType::To_Gas | FSType::To_Fluid)) {
          ++count;
        }
      }

      std::array<T*, LatSet::q> exmasscell = this->template getField<EXCESSMASS<T,LatSet::q>>().getArray(id);
      *(exmasscell[0]) = T{};
      if (count > 0) {
        T excessmassk = excessmass / count;
        for (int k = 1; k < LatSet::q; ++k) {
          const std::size_t idn = id + this->Delta_Index[k];
          if (util::isFlag(this->template getField<STATE>().get(idn), FSType::Interface) &&
              !util::isFlag(this->template getField<STATE>().get(idn), FSType::To_Gas | FSType::To_Fluid)) {
            *(exmasscell[k]) = excessmassk;
          }
        }
      } else {
        *(exmasscell[0]) = excessmass;
      }
    }
  }
}


template <typename T, typename LatSet, typename TypePack>
void FreeSurface2D<T, LatSet, TypePack>::FinalizeConversion() {
  // update state
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();
      Cell<T, LatSet, TypePack> cell(id, NS);
      if (util::isFlag(this->template getField<STATE>().get(id), FSType::To_Fluid)) {
        this->template getField<STATE>().SetField(id, FSType::Fluid);
        this->template getField<VOLUMEFRAC<T>>().SetField(id, T(1));
        this->template getField<MASS<T>>().get(id) += this->template getField<EXCESSMASS<T,LatSet::q>>().template get<0>(id);
      } else if (util::isFlag(this->template getField<STATE>().get(id), FSType::To_Gas)) {
        this->template getField<STATE>().SetField(id, FSType::Gas);
        this->template getField<VOLUMEFRAC<T>>().SetField(id, T{});
        this->template getField<MASS<T>>().get(id) += this->template getField<EXCESSMASS<T,LatSet::q>>().template get<0>(id);
        cell.template get<VELOCITY<T,LatSet::d>>().clear();
      } else if (util::isFlag(this->template getField<STATE>().get(id), FSType::To_Interface)) {
        this->template getField<STATE>().SetField(id, FSType::Interface);
      }
    }
  }
}

template <typename T, typename LatSet, typename TypePack>
void FreeSurface2D<T, LatSet, TypePack>::CollectExcessMass() {
  // stream
  for (int i = 1; i < LatSet::q; ++i) {
    this->template getField<EXCESSMASS<T,LatSet::q>>().getField(i).rotate(this->Delta_Index[i]);
  }
  // collect
  for (int j = NS.getOverlap(); j < NS.getNy() - NS.getOverlap(); ++j) {
    for (int i = NS.getOverlap(); i < NS.getNx() - NS.getOverlap(); ++i) {
      const std::size_t id = i + j * NS.getNx();
      if (util::isFlag(this->template getField<STATE>().get(id), FSType::Interface | FSType::Fluid)) {
          T exmass_sum = this->template getField<EXCESSMASS<T,LatSet::q>>().template get<0>(id);
        for (int k = 1; k < LatSet::q; ++k) {
          exmass_sum += this->template getField<EXCESSMASS<T,LatSet::q>>().get(id, k);
        }
        this->template getField<MASS<T>>().get(id) += exmass_sum;

        if (util::isFlag(this->template getField<STATE>().get(id), FSType::Interface)) {
          Cell<T, LatSet, TypePack> cell(id, NS);
          T rho = moment::rho<Cell<T, LatSet, TypePack>>::get(cell);
          this->template getField<VOLUMEFRAC<T>>().SetField(id, this->template getField<MASS<T>>().get(id) / rho);
        }
      }
    }
  }
  // clear
  this->template getField<EXCESSMASS<T,LatSet::q>>().Init();
}

}
// namespace FreeSurface
