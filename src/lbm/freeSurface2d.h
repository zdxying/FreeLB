/* This file is part of FreeLB, modified from openLB with the following copyright notice:
 *
 * // start of the original OpenLB's copyright notice
 * 
 * This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 * 
 * // end of the original OpenLB's copyright notice
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

// free surface model
#pragma once

#include "data_struct/block_lattice.h"
#include "lbm/freeSurface.h"

// namespace olbfs: openlb's implementation of free surface model

namespace olbfs {

template <typename LATTICEMANTYPE>
struct FreeSurfaceHelper2D {
  using LATTICEMAN = LATTICEMANTYPE;
  using T = typename LATTICEMAN::FloatType;
  using LatSet = typename LATTICEMAN::LatticeSet;

  static void Init(LATTICEMAN& latman) {
    // set interface cells
    latman.template getField<STATE>().forEach([&](auto& field, std::size_t id) {
      const auto& block = field.getBlock();
      if (util::isFlag(field.get(id), FSType::Fluid)) {
        for (int i = 1; i < LatSet::q; ++i) {
          std::size_t idn = id + latset::c<LatSet>(i) * block.getProjection();
          Vector<T, LatSet::d> loc_t = block.getLoc_t(idn);
          if (block.IsInside(loc_t)) {
            if (util::isFlag(field.get(idn), FSType::Gas)) {
              util::removeFlag(FSType::Gas, util::underlyingRef(field.get(idn)));
              util::addFlag(FSType::Interface, util::underlyingRef(field.get(idn)));
            }
          }
        }
      }
    });

    latman.template getField<STATE>().NormalCommunicate();
#ifdef MPI_ENABLED
    latman.template getField<STATE>().MPINormalCommunicate(count);
#endif

    // set mass and volume fraction
    latman.template getField<MASS<T>>().forEach(
      latman.template getField<STATE>(), FSType::Fluid,
      [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
    latman.template getField<MASS<T>>().forEach(
      latman.template getField<STATE>(), FSType::Interface,
      [&](auto& field, std::size_t id) { field.SetField(id, T{0.5}); });

    latman.template getField<MASS<T>>().NormalCommunicate();
#ifdef MPI_ENABLED
    latman.template getField<MASS<T>>().MPINormalCommunicate(count);
#endif

    latman.template getField<VOLUMEFRAC<T>>().forEach(
      latman.template getField<STATE>(), FSType::Fluid,
      [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
    latman.template getField<VOLUMEFRAC<T>>().forEach(
      latman.template getField<STATE>(), FSType::Interface,
      [&](auto& field, std::size_t id) { field.SetField(id, T{0.5}); });

    latman.template getField<VOLUMEFRAC<T>>().NormalCommunicate();
#ifdef MPI_ENABLED
    latman.template getField<VOLUMEFRAC<T>>().MPINormalCommunicate(count);
#endif
  }
};


// mass transfer
template <typename CELLTYPE>
struct MassTransfer2D {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  // this is from openLB's struct NeighbourInfo in FreeSurfaceHelpers.h
  struct NbrInfo {
    bool fluid_nbr = false;
    bool gas_nbr = false;
    int interface_nbrs = 0;

    NbrInfo(CELL& cell) {
      for (int k = 1; k < LatSet::q; ++k) {
        auto iflag = cell.template getField<STATE>().get(cell.getNeighborId(k));
        if (util::isFlag(iflag, FSType::Fluid)) {
          fluid_nbr = true;
        } else if (util::isFlag(iflag, FSType::Gas)) {
          gas_nbr = true;
        } else if (util::isFlag(iflag, FSType::Interface)) {
          ++interface_nbrs;
        }
      }
    }
  };

  static void apply(CELL& cell) {
    // openlb removes all cells' transition flags here
    cell.template get<FLAG>() = FSFlag::None;

    if (util::isFlag(cell.template get<STATE>(), FSType::Interface)) {
      // mass transfer
      T mass_tmp = cell.template get<MASS<T>>();
      // cell's nbr info
      NbrInfo cellNbrInfo(cell);

      for (int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        int kopp = latset::opp<LatSet>(k);

        if (util::isFlag(celln.template get<STATE>(), FSType::Fluid)) {
          mass_tmp += cell[kopp] - celln[k];
        } else if (util::isFlag(celln.template get<STATE>(), FSType::Interface)) {
          // celln's nbr info
          NbrInfo cellnNbrInfo(celln);
          T massflow{};
          // openlb deletes latset::w<LatSet>(k) term cause it is already contained in fi
          // not that latset::w<LatSet>(k) = latset::w<LatSet>(kopp)
          if (!cellNbrInfo.fluid_nbr) {
            if (!cellnNbrInfo.fluid_nbr) {
              if (cellNbrInfo.interface_nbrs < cellnNbrInfo.interface_nbrs) {
                massflow = -celln[k] + latset::w<LatSet>(k);
              } else if (cellNbrInfo.interface_nbrs > cellnNbrInfo.interface_nbrs) {
                massflow = cell[kopp] - latset::w<LatSet>(k);
              } else {
                massflow = cell[kopp] - celln[k];
              }
            } else {
              massflow = -celln[k] + latset::w<LatSet>(k);
            }
          } else if (!cellNbrInfo.gas_nbr) {
            if (!cellnNbrInfo.gas_nbr) {
              if (cellNbrInfo.interface_nbrs < cellnNbrInfo.interface_nbrs) {
                massflow = cell[kopp] - latset::w<LatSet>(k);
              } else if (cellNbrInfo.interface_nbrs > cellnNbrInfo.interface_nbrs) {
                massflow = -celln[k] + latset::w<LatSet>(k);
              } else {
                massflow = cell[kopp] - celln[k];
              }
            } else {
              massflow = cell[kopp] - latset::w<LatSet>(k);
            }
          } else {
            if (!cellnNbrInfo.fluid_nbr) {
              massflow = cell[kopp] - latset::w<LatSet>(k);
            } else if (!cellnNbrInfo.gas_nbr) {
              massflow = -celln[k] + latset::w<LatSet>(k);
            } else {
              massflow = cell[kopp] - celln[k];
            }
          }

          mass_tmp += massflow * T(0.5) * (getClampedVOF(cell) + getClampedVOF(celln));
        }
      }

      cell.template get<MASS<T>>() = mass_tmp;

      // reconstruct pop streamed in from a gas cell
      T curvature{};
      if (cell.template get<Surface_Tension_Enabled>()) {
        if (cellNbrInfo.gas_nbr) curvature = ComputeCurvature(cell);
      }
      T rho_gas = T(1) - T(6) * cell.template get<Surface_Tension_Parameter<T>>() * curvature;
      const Vector<T, LatSet::d>& u = cell.template get<PREVIOUS_VELOCITY<T, LatSet::d>>();
      T u2 = u.getnorm2();
      for (int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        if (util::isFlag(celln.template get<STATE>(), FSType::Gas)) {
          // fiopp = feqiopp(rho_gas) + feqi(rho_gas) - fi(x+ei)
          cell[latset::opp<LatSet>(k)] =
            equilibrium::SecondOrder<CELL>::get(k, u, rho_gas, u2) +
            equilibrium::SecondOrder<CELL>::get(latset::opp<LatSet>(k), u, rho_gas, u2) -
            celln[k];
        }
      }

      // openlb uses cell.computeRho();
      // however FreeLB can't get dynamics from a single cell
      // however this should be ok since if there's no mass source, rho computation should be the same
      T rho = moment::rho<CELL>::get(cell);

      // transition by mass criterion
      if (mass_tmp < -cell.template get<VOF_Trans_Th<T>>() * rho ||
      (mass_tmp < cell.template get<Lonely_Th<T>>() * rho && !cellNbrInfo.fluid_nbr)) {
        cell.template get<FLAG>() = FSFlag::To_Gas;
        return;
      } else if (mass_tmp > (T(1) + cell.template get<VOF_Trans_Th<T>>()) * rho ||
      (mass_tmp > (T(1) - cell.template get<Lonely_Th<T>>()) * rho && !cellNbrInfo.gas_nbr)) {
        cell.template get<FLAG>() = FSFlag::To_Fluid;
        return;
      } else if (cellNbrInfo.interface_nbrs == 0) {
        if (!cellNbrInfo.gas_nbr) {
          cell.template get<FLAG>() = FSFlag::To_Fluid;
        }
      }
    }
  }
};

// to fluid neighbor conversion
template <typename CELLTYPE>
struct ToFluidNbrConversion2D {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {

    if (util::isFlag(cell.template get<STATE>(), FSType::Gas)) {
      if (hasNeighborFlag(cell, FSFlag::To_Fluid)) {
        // set to_interface
        cell.template get<FLAG>() = FSFlag::To_Interface;
				// init fi using [fluid|interface] neighbor cells
				T averho{};
				Vector<T, LatSet::d> aveu{};
				int count{};
				for (int k = 1; k < LatSet::q; ++k) {
					CELL celln = cell.getNeighbor(k);
					if (util::isFlag(celln.template get<STATE>(), (FSType::Fluid | FSType::Interface))) {
						// openlb uses: cellC.computeRhoU(rho_tmp, u_tmp);
            // however FreeLB can't get dynamics from a single cell
            // we have to use forceRhou here
            T rho{};
            Vector<T, LatSet::d> u{};
            moment::forceRhou<CELL, force::ConstForce<CELL>>::apply(celln, rho, u);
            averho += rho;
            aveu += u;
						++count;
					}
				}
				if (count > 0) {
					averho /= count;
					aveu /= count;
				}
				// set fi, openlb uses: cell.iniEquilibrium(rho_avg, u_avg);
				T aveu2 = aveu.getnorm2();
				for (int k = 0; k < LatSet::q; ++k) {
					cell[k] = equilibrium::SecondOrder<CELL>::get(k, aveu, averho, aveu2);
				}
      }
    } else if (util::isFlag(cell.template get<FLAG>(), FSFlag::To_Gas)) {
      if (hasNeighborFlag(cell, FSFlag::To_Fluid)) {
        // remove to_gas flag, in openlb: clear transition flag
        cell.template get<FLAG>() = FSFlag::None;
      }
    }
  }
};

// to gas neighbor conversion
template <typename CELLTYPE>
struct ToGasNbrConversion2D {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {

    if (util::isFlag(cell.template get<STATE>(), FSType::Fluid)) {
      if (hasNeighborFlag(cell, FSFlag::To_Gas)) {
        // set to_interface
        cell.template get<FLAG>() = FSFlag::To_Interface;
				cell.template get<MASS<T>>() = moment::rho<CELL>::get(cell);
      }
    }
  }
};

// interface excess mass
template <typename CELLTYPE>
struct InterfaceExcessMass2D {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    // for interface cells to be converted to fluid or gas
    // excess mass is distributed to interface neighbors
		if (util::isFlag(cell.template get<STATE>(), FSType::Interface)) {

			T rho = moment::rho<CELL>::get(cell);
			T mass = cell.template get<MASS<T>>();
			T mass_excess{};

			if (util::isFlag(cell.template get<FLAG>(), FSFlag::To_Gas)) {
				mass_excess = mass;
				cell.template get<MASS<T>>() = T{};
    	} else if (util::isFlag(cell.template get<FLAG>(), FSFlag::To_Fluid)) {
      	mass_excess = mass - rho;
				cell.template get<MASS<T>>() = rho;
			} else {
				return;
			}

			// find neighbors
			int count{};
			for (int k = 1; k < LatSet::q; ++k) {
				CELL celln = cell.getNeighbor(k);
				if (util::isFlag(celln.template get<STATE>(), FSType::Interface) &&
						!util::isFlag(celln.template get<FLAG>(), std::uint8_t(255))) {
					++count;
				}
			}

      Vector<T, LatSet::q>& mass_ex_vec = cell.template get<MASSEX<T, LatSet::q>>();
      mass_ex_vec[0] = T{};
			if (count > 0) {
        const T mass_excess_frac = mass_excess / count;
        
        for (int k = 1; k < LatSet::q; ++k) {
          CELL celln = cell.getNeighbor(k);
          if (util::isFlag(celln.template get<STATE>(), FSType::Interface) &&
              !util::isFlag(celln.template get<FLAG>(), std::uint8_t(255))) {
            mass_ex_vec[k] = mass_excess_frac;
          } else {
            mass_ex_vec[k] = T{};
          }
        }
      } else {
				mass_ex_vec[0] = mass_excess;
        for (int k = 1; k < LatSet::q; ++k) {
          mass_ex_vec[k] = T{};
        }
			}
		}
	}
};

// finalize conversion
template <typename CELLTYPE>
struct FinalizeConversion2D {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    // update state
    if (util::isFlag(cell.template get<FLAG>(), FSFlag::To_Fluid)) {
      cell.template get<STATE>() = FSType::Fluid;
      cell.template get<VOLUMEFRAC<T>>() = T(1);
      cell.template get<MASS<T>>() += cell.template get<MASSEX<T, LatSet::q>>()[0];
    } else if (util::isFlag(cell.template get<FLAG>(), FSFlag::To_Gas)) {
      cell.template get<STATE>() = FSType::Gas;
      cell.template get<VOLUMEFRAC<T>>() = T{};
      cell.template get<MASS<T>>() += cell.template get<MASSEX<T, LatSet::q>>()[0];
    } else if (util::isFlag(cell.template get<FLAG>(), FSFlag::To_Interface)) {
      cell.template get<STATE>() = FSType::Interface;
    }
		// collect excess mass
		if (util::isFlag(cell.template get<STATE>(), FSType::Interface)) {
      T collected_excess{};
			for (int k = 1; k < LatSet::q; ++k) {
				CELL celln = cell.getNeighbor(k);
        if (util::isFlag(celln.template get<FLAG>(), FSFlag::To_Fluid | FSFlag::To_Gas))
				collected_excess += celln.template get<MASSEX<T, LatSet::q>>()[latset::opp<LatSet>(k)];
			}
			T mass_tmp = cell.template get<MASS<T>>();
			mass_tmp += collected_excess;

      // openlb uses: cell.computeRhoU(rho, u_tmp);
      // however FreeLB can't get dynamics from a single cell
      // we have to use forceRhou here
			T rho;
			Vector<T, LatSet::d> u{};
			moment::forceRhou<CELL, force::ConstForce<CELL>>::apply(cell, rho, u);

			cell.template get<MASS<T>>() = mass_tmp;
			cell.template get<VOLUMEFRAC<T>>() = mass_tmp / rho;
			cell.template get<PREVIOUS_VELOCITY<T, LatSet::d>>() = u;

    } else if (util::isFlag(cell.template get<STATE>(), FSType::Fluid)) {
      T collected_excess{};
			for (int k = 1; k < LatSet::q; ++k) {
				CELL celln = cell.getNeighbor(k);
				if (util::isFlag(celln.template get<FLAG>(), FSFlag::To_Fluid | FSFlag::To_Gas))
				collected_excess += celln.template get<MASSEX<T, LatSet::q>>()[latset::opp<LatSet>(k)];
			}
			cell.template get<MASS<T>>() += collected_excess;
    }
  }
};

// apply all the free surface dynamics
template <typename LATTICEMANAGERTYPE>
struct FreeSurfaceApply2D {
  using CELL = typename LATTICEMANAGERTYPE::CellType;
  using BLOCKLAT = typename LATTICEMANAGERTYPE::BLOCKLATTICE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void Apply(LATTICEMANAGERTYPE& latManager, std::int64_t count) {
    // mass transfer
    latManager.template ApplyInnerCellDynamics<MassTransfer2D<CELL>>(count);

    latManager.template getField<FLAG>().NormalCommunicate(count);
#ifdef MPI_ENABLED
    latManager.template getField<FLAG>().MPINormalCommunicate(count);
#endif
    latManager.template getField<MASS<T>>().CommunicateAll(count);

    // communicate reconstructed pops streamed in from a gas cell
    // this is NOT a post-stream process, so we must communicate fi in each direction
		latManager.FullCommunicate(count);


    // to fluid neighbor conversion
    latManager.template ApplyInnerCellDynamics<ToFluidNbrConversion2D<CELL>>(count);

    latManager.template getField<FLAG>().NormalCommunicate(count);
#ifdef MPI_ENABLED
    latManager.template getField<FLAG>().MPINormalCommunicate(count);
#endif

    // communicate equilibrium fi from nbr Fluid/Interface cells' rho and u for a Gas->Interface cell
    // this is NOT a post-stream process, so we must communicate fi in each direction
    latManager.FullCommunicate(count);


    // to gas neighbor conversion
    latManager.template ApplyInnerCellDynamics<ToGasNbrConversion2D<CELL>>(count);

    latManager.template getField<FLAG>().NormalCommunicate(count);
#ifdef MPI_ENABLED
    latManager.template getField<FLAG>().MPINormalCommunicate(count);
#endif
		latManager.template getField<MASS<T>>().CommunicateAll(count);


    // interface excess mass
    latManager.template ApplyInnerCellDynamics<InterfaceExcessMass2D<CELL>>(count);

    latManager.template getField<MASS<T>>().CommunicateAll(count);
    latManager.template getField<MASSEX<T, LatSet::q>>().CommunicateAll(count);


    // finalize conversion
    latManager.template ApplyInnerCellDynamics<FinalizeConversion2D<CELL>>(count);

    latManager.template getField<STATE>().NormalCommunicate(count);
#ifdef MPI_ENABLED
    latManager.template getField<STATE>().MPINormalCommunicate(count);
#endif
    latManager.template getField<MASS<T>>().CommunicateAll(count);
    latManager.template getField<VOLUMEFRAC<T>>().CommunicateAll(count);
    latManager.template getField<PREVIOUS_VELOCITY<T,LatSet::d>>().CommunicateAll(count);


    // clear EXCESSMASS<T,LatSet::q>
    latManager.ForEachBlockLattice(
      count, [&](auto& blocklat) { blocklat.template getField<MASSEX<T, LatSet::q>>().Init(Vector<T,LatSet::q>{}); });
  }
};

}  // namespace olbfs