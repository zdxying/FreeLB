/* This file is part of FreeLB, modified from openLB and FluidX3D with the following copyright notice:
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
 * FluidX3D: https://github.com/ProjectPhysX/FluidX3Ds
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

#include "lbm/freeSurfaceUtils.h"


// namespace olbfs: openlb's implementation of free surface model

namespace olbfs {

template <typename LATTICEMANTYPE>
struct FreeSurfaceHelper {
  using LATTICEMAN = LATTICEMANTYPE;
  using T = typename LATTICEMAN::FloatType;
  using LatSet = typename LATTICEMAN::LatticeSet;

  static void Init(LATTICEMAN& latman) {
    // set interface cells
    latman.template getField<STATE>().forEach([&](auto& field, std::size_t id) {
      const auto& block = field.getBlock();
      if (util::isFlag(field.get(id), FSType::Fluid)) {
        for (unsigned int i = 1; i < LatSet::q; ++i) {
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

    latman.template getField<STATE>().AllNormalCommunicate();

    // set mass and volume fraction
    latman.template getField<MASS<T>>().forEach(
      latman.template getField<STATE>(), FSType::Fluid,
      [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
    latman.template getField<MASS<T>>().forEach(
      latman.template getField<STATE>(), FSType::Interface,
      [&](auto& field, std::size_t id) { field.SetField(id, T{0.5}); });

    latman.template getField<MASS<T>>().AllNormalCommunicate();

    latman.template getField<VOLUMEFRAC<T>>().forEach(
      latman.template getField<STATE>(), FSType::Fluid | FSType::Wall,
      [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
    latman.template getField<VOLUMEFRAC<T>>().forEach(
      latman.template getField<STATE>(), FSType::Interface,
      [&](auto& field, std::size_t id) { field.SetField(id, T{0.5}); });

    latman.template getField<VOLUMEFRAC<T>>().AllNormalCommunicate();
  }
};


// mass transfer
template <typename CELLTYPE>
struct MassTransfer {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  // this is from openLB's struct NeighbourInfo in FreeSurfaceHelpers.h
  struct NbrInfo {
    bool fluid_nbr = false;
    bool gas_nbr = false;
    int interface_nbrs = 0;

    NbrInfo(CELL& cell) {
      for (unsigned int k = 1; k < LatSet::q; ++k) {
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

      for (unsigned int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        int kopp = latset::opp<LatSet>(k);

        if (util::isFlag(celln.template get<STATE>(), FSType::Fluid)) {
          mass_tmp += cell[kopp] - celln[k];
        } else if (util::isFlag(celln.template get<STATE>(), FSType::Interface)) {
          // celln's nbr info
          // 2 layers of overlaped cells is needed to avoid accessing non-existing cell nbr' nbr 
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
      for (unsigned int k = 1; k < LatSet::q; ++k) {
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
template <typename CELLTYPE, unsigned int scalardir>
struct ToFluidNbrConversion {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static constexpr bool hasForce = (
    CELL::template hasField<FORCE<T, LatSet::d>>() || 
    CELL::template hasField<SCALARFORCE<T>>() || 
    CELL::template hasField<CONSTFORCE<T, LatSet::d>>() || 
    CELL::template hasField<SCALARCONSTFORCE<T>>());
  // here only 2 force schemes are considered
  using ForceScheme = std::conditional_t<CELL::template hasField<FORCE<T, LatSet::d>>(), force::Force<CELL>, 
                        std::conditional_t<CELL::template hasField<SCALARFORCE<T>>(), force::ScalarForce<CELL, scalardir>, 
                          std::conditional_t<CELL::template hasField<CONSTFORCE<T, LatSet::d>>(), force::ConstForce<CELL>, 
                            std::conditional_t<CELL::template hasField<SCALARCONSTFORCE<T>>(), force::ScalarConstForce<CELL, scalardir>, void>>>>;

  static void apply(CELL& cell) {

    if (util::isFlag(cell.template get<STATE>(), FSType::Gas)) {
      if (hasNeighborFlag(cell, FSFlag::To_Fluid)) {
        // set to_interface
        cell.template get<FLAG>() = FSFlag::To_Interface;
				// init fi using [fluid|interface] neighbor cells
				T averho{};
				Vector<T, LatSet::d> aveu{};
				int count{};
				for (unsigned int k = 1; k < LatSet::q; ++k) {
					CELL celln = cell.getNeighbor(k);
					if (util::isFlag(celln.template get<STATE>(), (FSType::Fluid | FSType::Interface))) {
						// openlb uses: cellC.computeRhoU(rho_tmp, u_tmp);
            // however FreeLB can't get dynamics from a single cell
            // we have to use forcerhoU here
            T rho{};
            Vector<T, LatSet::d> u{};
            if constexpr (hasForce) moment::forcerhoU<CELL, ForceScheme>::apply(celln, rho, u);
            else moment::rhoU<CELL>::apply(celln, rho, u);
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
				for (unsigned int k = 0; k < LatSet::q; ++k) {
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
struct ToGasNbrConversion {
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
struct InterfaceExcessMass {
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
			for (unsigned int k = 1; k < LatSet::q; ++k) {
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
        
        for (unsigned int k = 1; k < LatSet::q; ++k) {
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
        for (unsigned int k = 1; k < LatSet::q; ++k) {
          mass_ex_vec[k] = T{};
        }
			}
		}
	}
};

// finalize conversion
template <typename CELLTYPE>
struct FinalizeConversion {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  // here only 2 force schemes are considered
  using ForceScheme = std::conditional_t<CELL::template hasField<CONSTFORCE<T, LatSet::d>>(), 
  force::ConstForce<CELL>, force::ScalarConstForce<CELL>>;

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
			for (unsigned int k = 1; k < LatSet::q; ++k) {
				CELL celln = cell.getNeighbor(k);
        if (util::isFlag(celln.template get<FLAG>(), FSFlag::To_Fluid | FSFlag::To_Gas))
				collected_excess += celln.template get<MASSEX<T, LatSet::q>>()[latset::opp<LatSet>(k)];
			}
			T mass_tmp = cell.template get<MASS<T>>();
			mass_tmp += collected_excess;

      // openlb uses: cell.computeRhoU(rho, u_tmp);
      // however FreeLB can't get dynamics from a single cell
      // we have to use forcerhoU here
			T rho;
			Vector<T, LatSet::d> u{};
			moment::forcerhoU<CELL, ForceScheme>::apply(cell, rho, u);

			cell.template get<MASS<T>>() = mass_tmp;
			cell.template get<VOLUMEFRAC<T>>() = mass_tmp / rho;
			cell.template get<PREVIOUS_VELOCITY<T, LatSet::d>>() = u;

    } else if (util::isFlag(cell.template get<STATE>(), FSType::Fluid)) {
      T collected_excess{};
			for (unsigned int k = 1; k < LatSet::q; ++k) {
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
struct FreeSurfaceApply {
  using CELL = typename LATTICEMANAGERTYPE::CellType;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  template <unsigned int scalardir = 2>
  static void Apply(LATTICEMANAGERTYPE& latManager) {
    // mass transfer
    latManager.template ApplyInnerCellDynamics<MassTransfer<CELL>>();
    latManager.template getField<FLAG>().AllNormalCommunicate();
    latManager.template getField<MASS<T>>().AllNormalCommunicate();
    // communicate reconstructed pops streamed in from a gas cell
    // this is NOT a post-stream process, so we must communicate fi in each direction
    latManager.NormalAllCommunicate();


    // to fluid neighbor conversion
    latManager.template ApplyInnerCellDynamics<ToFluidNbrConversion<CELL, scalardir>>();
    latManager.template getField<FLAG>().AllNormalCommunicate();
    // communicate equilibrium fi from nbr Fluid/Interface cells' rho and u for a Gas->Interface cell
    // this is NOT a post-stream process, so we must communicate fi in each direction
    latManager.NormalAllCommunicate();


    // to gas neighbor conversion
    latManager.template ApplyInnerCellDynamics<ToGasNbrConversion<CELL>>();
    latManager.template getField<FLAG>().AllNormalCommunicate();
		latManager.template getField<MASS<T>>().AllNormalCommunicate();


    // interface excess mass
    latManager.template ApplyInnerCellDynamics<InterfaceExcessMass<CELL>>();
    latManager.template getField<MASS<T>>().AllNormalCommunicate();
    latManager.template getField<MASSEX<T, LatSet::q>>().AllNormalCommunicate();


    // finalize conversion
    latManager.template ApplyInnerCellDynamics<FinalizeConversion<CELL>>();
    latManager.template getField<STATE>().AllNormalCommunicate();
    latManager.template getField<MASS<T>>().AllNormalCommunicate();
    latManager.template getField<VOLUMEFRAC<T>>().AllNormalCommunicate();
    latManager.template getField<PREVIOUS_VELOCITY<T,LatSet::d>>().AllNormalCommunicate();


    // clear EXCESSMASS<T,LatSet::q>
    latManager.ForEachBlockLattice([&](auto& blocklat) { blocklat.template getField<MASSEX<T, LatSet::q>>().Init(Vector<T,LatSet::q>{}); });
  }

  template <unsigned int scalardir = 2>
  static void Apply(LATTICEMANAGERTYPE& latManager, std::int64_t count) {
    // mass transfer
    latManager.template ApplyInnerCellDynamics<MassTransfer<CELL>>(count);
    latManager.template getField<FLAG>().AllNormalCommunicate(count);
    latManager.template getField<MASS<T>>().AllNormalCommunicate(count);
    // communicate reconstructed pops streamed in from a gas cell
    // this is NOT a post-stream process, so we must communicate fi in each direction
		latManager.NormalAllCommunicate(count);


    // to fluid neighbor conversion
    latManager.template ApplyInnerCellDynamics<ToFluidNbrConversion<CELL, scalardir>>(count);
    latManager.template getField<FLAG>().AllNormalCommunicate(count);
    // communicate equilibrium fi from nbr Fluid/Interface cells' rho and u for a Gas->Interface cell
    // this is NOT a post-stream process, so we must communicate fi in each direction
    latManager.NormalAllCommunicate(count);


    // to gas neighbor conversion
    latManager.template ApplyInnerCellDynamics<ToGasNbrConversion<CELL>>(count);
    latManager.template getField<FLAG>().AllNormalCommunicate(count);
		latManager.template getField<MASS<T>>().AllNormalCommunicate(count);


    // interface excess mass
    latManager.template ApplyInnerCellDynamics<InterfaceExcessMass<CELL>>(count);
    latManager.template getField<MASS<T>>().AllNormalCommunicate(count);
    latManager.template getField<MASSEX<T, LatSet::q>>().AllNormalCommunicate(count);


    // finalize conversion
    latManager.template ApplyInnerCellDynamics<FinalizeConversion<CELL>>(count);
    latManager.template getField<STATE>().AllNormalCommunicate(count);
    latManager.template getField<MASS<T>>().AllNormalCommunicate(count);
    latManager.template getField<VOLUMEFRAC<T>>().AllNormalCommunicate(count);
    latManager.template getField<PREVIOUS_VELOCITY<T,LatSet::d>>().AllNormalCommunicate(count);


    // clear EXCESSMASS<T,LatSet::q>
    latManager.ForEachBlockLattice(
      count, [&](auto& blocklat) { blocklat.template getField<MASSEX<T, LatSet::q>>().Init(Vector<T,LatSet::q>{}); });
  }
};

// to fluid neighbor conversion
template <typename NSCELLTYPE, typename XCELLTYPE, unsigned int scalardir>
struct CoupledToFluidNbrConversion {
  using CELL = NSCELLTYPE;
  using XCELL = XCELLTYPE;
  using XGenericRho = typename XCELL::GenericRho;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static constexpr bool hasForce = (
    CELL::template hasField<FORCE<T, LatSet::d>>() || 
    CELL::template hasField<SCALARFORCE<T>>() || 
    CELL::template hasField<CONSTFORCE<T, LatSet::d>>() || 
    CELL::template hasField<SCALARCONSTFORCE<T>>());
  static_assert(CELL::template hasField<XGenericRho>(), "CoupledToFluidNbrConversion: CELL must has ref of XCELL's GenericRho");
  // here only 2 force schemes are considered
  using ForceScheme = 
    std::conditional_t<CELL::template hasField<FORCE<T, LatSet::d>>(), force::Force<CELL>, 
      std::conditional_t<CELL::template hasField<SCALARFORCE<T>>(), force::ScalarForce<CELL, scalardir>, 
        std::conditional_t<CELL::template hasField<CONSTFORCE<T, LatSet::d>>(), force::ConstForce<CELL>, 
          std::conditional_t<CELL::template hasField<SCALARCONSTFORCE<T>>(), force::ScalarConstForce<CELL, scalardir>, void>>>>;

  static void apply(CELL& cell, XCELL& xcell) {

    if (util::isFlag(cell.template get<STATE>(), FSType::Gas)) {
      if (hasNeighborFlag(cell, FSFlag::To_Fluid)) {
        // set to_interface
        cell.template get<FLAG>() = FSFlag::To_Interface;
				// init fi using [fluid|interface] neighbor cells
				T averho{};
        T avexrho{};
				Vector<T, LatSet::d> aveu{};
				int count{};
				for (unsigned int k = 1; k < LatSet::q; ++k) {
					CELL celln = cell.getNeighbor(k);
					if (util::isFlag(celln.template get<STATE>(), (FSType::Fluid | FSType::Interface))) {
						// openlb uses: cellC.computeRhoU(rho_tmp, u_tmp);
            // however FreeLB can't get dynamics from a single cell
            // we have to use forcerhoU here
            T rho{};
            Vector<T, LatSet::d> u{};
            if constexpr (hasForce) moment::forcerhoU<CELL, ForceScheme>::apply(celln, rho, u);
            else moment::rhoU<CELL>::apply(celln, rho, u);
            averho += rho;
            avexrho += celln.template get<XGenericRho>();
            aveu += u;
						++count;
					}
				}
				if (count > 0) {
					averho /= count;
					aveu /= count;
          avexrho /= count;
				}
        cell.template get<XGenericRho>() = avexrho;
				// set fi, openlb uses: cell.iniEquilibrium(rho_avg, u_avg);
				T aveu2 = aveu.getnorm2();
				for (unsigned int k = 0; k < LatSet::q; ++k) {
					cell[k] = equilibrium::SecondOrder<CELL>::get(k, aveu, averho, aveu2);
          xcell[k] = equilibrium::SecondOrder<XCELL>::get(k, aveu, avexrho, aveu2);
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

template <typename LatManagerCouplingTYPE>
struct CoupledFreeSurfaceApply {
  using CELL = typename LatManagerCouplingTYPE::CELL0;
  using XCELL = typename LatManagerCouplingTYPE::CELL1;
  using GenericRho = typename CELL::GenericRho;
  using XGenericRho = typename XCELL::GenericRho;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  template <unsigned int scalardir = 2>
  static void Apply(LatManagerCouplingTYPE& latManagerCoupling) {
    auto& latManager = latManagerCoupling.getLat0();
    auto& xlatManager = latManagerCoupling.getLat1();
    // mass transfer
    latManager.template ApplyInnerCellDynamics<MassTransfer<CELL>>();
    latManager.template getField<FLAG>().AllNormalCommunicate();
    latManager.template getField<MASS<T>>().AllNormalCommunicate();
    // communicate reconstructed pops streamed in from a gas cell
    // this is NOT a post-stream process, so we must communicate fi in each direction
    latManager.NormalAllCommunicate();


    // to fluid neighbor conversion
    latManagerCoupling.template ApplyInnerCellDynamics<CoupledToFluidNbrConversion<CELL, XCELL, scalardir>>();
    latManager.template getField<FLAG>().AllNormalCommunicate();
    xlatManager.template getField<XGenericRho>().AllNormalCommunicate();
    // communicate equilibrium fi from nbr Fluid/Interface cells' rho and u for a Gas->Interface cell
    // this is NOT a post-stream process, so we must communicate fi in each direction
    latManager.NormalAllCommunicate();
    xlatManager.NormalAllCommunicate();


    // to gas neighbor conversion
    latManager.template ApplyInnerCellDynamics<ToGasNbrConversion<CELL>>();
    latManager.template getField<FLAG>().AllNormalCommunicate();
		latManager.template getField<MASS<T>>().AllNormalCommunicate();


    // interface excess mass
    latManager.template ApplyInnerCellDynamics<InterfaceExcessMass<CELL>>();
    latManager.template getField<MASS<T>>().AllNormalCommunicate();
    latManager.template getField<MASSEX<T, LatSet::q>>().AllNormalCommunicate();


    // finalize conversion
    latManager.template ApplyInnerCellDynamics<FinalizeConversion<CELL>>();
    latManager.template getField<STATE>().AllNormalCommunicate();
    latManager.template getField<MASS<T>>().AllNormalCommunicate();
    latManager.template getField<VOLUMEFRAC<T>>().AllNormalCommunicate();
    latManager.template getField<PREVIOUS_VELOCITY<T,LatSet::d>>().AllNormalCommunicate();


    // clear EXCESSMASS<T,LatSet::q>
    latManager.ForEachBlockLattice([&](auto& blocklat) { blocklat.template getField<MASSEX<T, LatSet::q>>().Init(Vector<T,LatSet::q>{}); });
  }
};

}  // namespace olbfs


// FluidX3D's implementation of free surface model

namespace fx3dfs {

template <typename LATTICEMANTYPE>
struct FreeSurfaceHelper {
  using LATTICEMAN = LATTICEMANTYPE;
  using T = typename LATTICEMAN::FloatType;
  using LatSet = typename LATTICEMAN::LatticeSet;

  static void Init(LATTICEMAN& latman) {
    // set interface cells
    latman.template getField<STATE>().forEach([&](auto& field, std::size_t id) {
      const auto& block = field.getBlock();
      if (util::isFlag(field.get(id), FSType::Fluid)) {
        for (unsigned int i = 1; i < LatSet::q; ++i) {
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

    // set mass and volume fraction
    latman.template getField<MASS<T>>().forEach(
      latman.template getField<STATE>(), FSType::Fluid,
      [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
    latman.template getField<MASS<T>>().forEach(
      latman.template getField<STATE>(), FSType::Interface,
      [&](auto& field, std::size_t id) { field.SetField(id, T{0.5}); });

    latman.template getField<MASS<T>>().NormalCommunicate();

    latman.template getField<VOLUMEFRAC<T>>().forEach(
      latman.template getField<STATE>(), FSType::Fluid | FSType::Wall,
      [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
    latman.template getField<VOLUMEFRAC<T>>().forEach(
      latman.template getField<STATE>(), FSType::Interface,
      [&](auto& field, std::size_t id) { field.SetField(id, T{0.5}); });

    latman.template getField<VOLUMEFRAC<T>>().NormalCommunicate();
  }
};

// mass transfer post process
template <typename CELLTYPE>
struct surface_post_process0 {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    // for interface and fluid cells
    if (util::isFlag(cell.template get<STATE>(), FSType::Interface | FSType::Fluid)) {
      // mass transfer
      T deltamass{};

      for (unsigned int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        deltamass += celln.template get<MASSEX<T>>();
      }

      if (util::isFlag(cell.template get<STATE>(), FSType::Fluid)) {
        for (unsigned int k = 1; k < LatSet::q; ++k) {
          CELL celln = cell.getNeighbor(k);
          deltamass += cell[latset::opp<LatSet>(k)] - celln[k];
        }
      } else {
        const T cellvof = computeVOF(cell);
        for (unsigned int k = 1; k < LatSet::q; ++k) {
          CELL celln = cell.getNeighbor(k);
          const int kopp = latset::opp<LatSet>(k);

          if (util::isFlag(celln.template get<STATE>(), FSType::Fluid)) {
            deltamass += cell[kopp] - celln[k];
          } else if (util::isFlag(celln.template get<STATE>(), FSType::Interface)) {
            deltamass +=
              (cell[kopp] - celln[k]) * T(0.5) * (cellvof + computeVOF(celln));
          }
        }
      }

      cell.template get<MASS<T>>() += deltamass;

      // reconstruct pop streamed in from a gas cell
      if (util::isFlag(cell.template get<STATE>(), FSType::Interface)) {

        if (hasNeighborType(cell, FSType::Gas)) {
          T curvature{};
          if (cell.template get<Surface_Tension_Enabled>()) {
            curvature = ComputeCurvature(cell);
          }
          T rho_gas =
            T(1) - T(6) * cell.template get<Surface_Tension_Parameter<T>>() * curvature;
          // in Lehmann's thesis: stream-collide-update rhoU - free surface,
          // the current velocity is used to calculate the reconstructed fi
          // FreeLB uses: update rhoU-colide-stream-free surface
          // so before call free surface post process, rho and u should be updated
          const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
          const T u2 = u.getnorm2();
          for (unsigned int k = 1; k < LatSet::q; ++k) {
            CELL celln = cell.getNeighbor(k);
            if (util::isFlag(celln.template get<STATE>(), FSType::Gas)) {
              // fiopp = feqiopp(rho_gas) + feqi(rho_gas) - fi(x+ei) reconstructed f
              cell[latset::opp<LatSet>(k)] =
                equilibrium::SecondOrder<CELL>::get(k, u, rho_gas, u2) +
                equilibrium::SecondOrder<CELL>::get(latset::opp<LatSet>(k), u, rho_gas,
                                                    u2) -
                celln[k];
            }
          }
        }

        // transition flag for interface cell
        const T epsilon = cell.template get<VOF_Trans_Th<T>>();
        const T rho = cell.template get<RHO<T>>();
        const T hi_th = rho * (T{1} + epsilon);
        const T lo_th = -rho * epsilon;
        if (cell.template get<MASS<T>>() > hi_th){
          util::addFlag(FSType::To_Fluid, util::underlyingRef(cell.template get<STATE>()));
          return;
        } else if (cell.template get<MASS<T>>() < lo_th){
          // util::addFlag(FSType::To_Gas, util::underlyingRef(cell.template get<STATE>()));
          // set to [to_gas] flag and cancel [interface] flag
          cell.template get<STATE>() = FSType::To_Gas;
          return;
        }
        // isolated interface cells
        if (!hasNeighborType(cell, FSType::Gas)) {
          util::addFlag(FSType::To_Fluid, util::underlyingRef(cell.template get<STATE>()));
          return;
        } else if (!hasNeighborType(cell, FSType::Fluid)) {
          // util::addFlag(FSType::To_Gas, util::underlyingRef(cell.template get<STATE>()));
          // set to [to_gas] flag and cancel [interface] flag
          cell.template get<STATE>() = FSType::To_Gas;
          return;
        }
      }
    }
  }
};
// possible cell type at the end of surface_post_process0
// [Fluid] [Interface] [Gas] [To_Fluid | Interface] [To_Gas]

// prevent neighbors from interface->fluid nodes to become/be gas nodes
template <typename CELLTYPE>
struct surface_post_process1 {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    // for correct communication, cell should be [To_Gas] or [Gas] 
    // while celln should be [To_Fluid]

    // [To_Fluid] neighbor conversion
    if (util::isFlag(cell.template get<STATE>(), FSType::To_Gas)) {
      if (hasNeighborType(cell, FSType::To_Fluid)) {
        // remove to_gas flag
        // util::removeFlag(FSType::To_Gas, util::underlyingRef(cell.template get<STATE>()));
        // set to [interface]
        cell.template get<STATE>() = FSType::Interface;
      }
    } else if (util::isFlag(cell.template get<STATE>(), FSType::Gas)) {
      if (hasNeighborType(cell, FSType::To_Fluid)) {
        // set to_interface
        util::addFlag(FSType::To_Interface, util::underlyingRef(cell.template get<STATE>()));
      }
    }
  }
};
// possible cell type at the end of surface_post_process1
// [Fluid] [Interface] [Gas] [To_Fluid | Interface] [To_Gas] [To_Interface | Gas]

// gas->interface and interface->gas
template <typename CELLTYPE>
struct surface_post_process2 {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    // initialize the fi of gas nodes that should become interface
    // [to_interface] flag is only for gas cells, fluid cells directly set to [interface]
    if (util::isFlag(cell.template get<STATE>(), FSType::To_Interface)) {
      // init fi using [fluid|interface] neighbor cells
      T averho{};
      Vector<T, LatSet::d> aveu{};
      int count{};
      for (unsigned int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        if (util::isFlag(celln.template get<STATE>(),
                         (FSType::Fluid | FSType::Interface))) {
          averho += celln.template get<RHO<T>>();
          aveu += celln.template get<VELOCITY<T, LatSet::d>>();
          ++count;
        }
      }
      averho /= count;
      aveu /= count;
      // set fi
      T aveu2 = aveu.getnorm2();
      for (unsigned int k = 0; k < LatSet::q; ++k) {
        cell[k] = equilibrium::SecondOrder<CELL>::get(k, aveu, averho, aveu2);
      }
    } else if (util::isFlag(cell.template get<STATE>(), FSType::Fluid)) {
      if (hasNeighborType(cell, FSType::To_Gas)) {
        // it said that fluid cell should be set to interface at once
        cell.template get<STATE>() = FSType::Interface;
      }
    } else if (util::isFlag(cell.template get<STATE>(), FSType::To_Fluid)){
      if (hasNeighborType(cell, FSType::To_Gas)) {
        // remove to_fluid flag
        util::removeFlag(FSType::To_Fluid, util::underlyingRef(cell.template get<STATE>()));
      }
    }

    // for correct communication, cell should be [To_Fluid] or [Fluid] 
    // while celln should be [To_Gas]
  }
};

// possible cell type at the end of surface_post_process1
// [Fluid] [Interface] [Gas] [To_Fluid | Interface] [To_Gas] [To_Interface | Gas]

// apply flag changes and calculate excess mass
template <typename CELLTYPE>
struct surface_post_process3 {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    T excessmass{};
    T vof{};
    T mass = cell.template get<MASS<T>>();
    const T rho = cell.template get<RHO<T>>();
    // note that transition flag like [to_interface] [to_gas] [to_fluid]
    // should be checked before state flag [interface] and [gas]
    if (util::isFlag(cell.template get<STATE>(), FSType::Fluid)) {
      excessmass = mass - rho;
      mass = rho;
      vof = T{1};
    } else if (util::isFlag(cell.template get<STATE>(), FSType::To_Fluid)) {
      cell.template get<STATE>() = FSType::Fluid;
      excessmass = mass - rho;
      mass = rho;
      vof = T{1};
    } else if (util::isFlag(cell.template get<STATE>(), FSType::To_Gas)) {
      cell.template get<STATE>() = FSType::Gas;
      excessmass = mass;
      mass = T{};
      vof = T{};
    } else if (util::isFlag(cell.template get<STATE>(), FSType::To_Interface)) {
      cell.template get<STATE>() = FSType::Interface;
      excessmass = (mass > rho) ? mass - rho : (mass < T{}) ? mass : T{};
      mass = std::clamp(mass, T{}, rho);
      vof = computeVOF(mass, rho);
    } else if (util::isFlag(cell.template get<STATE>(), FSType::Interface)) {
      excessmass = (mass > rho) ? mass - rho : (mass < T{}) ? mass : T{};
      mass = std::clamp(mass, T{}, rho);
      vof = computeVOF(mass, rho);
    } else if (util::isFlag(cell.template get<STATE>(), FSType::Gas)) {
      excessmass = mass;
      mass = T{};
      vof = T{};
    } 
    // distribute excess mass equally to all interface and fluid neighbors
    int count{};
    for (unsigned int k = 1; k < LatSet::q; ++k) {
      CELL celln = cell.getNeighbor(k);
      if (util::isFlag(celln.template get<STATE>(), 
      FSType::Interface | FSType::Fluid | FSType::To_Fluid | FSType::To_Interface)) {
        count++;
      }
    }
    // if excess mass can't be distributed to neighboring interface or fluid nodes, add
    // it to local mass (ensure mass conservation)
    mass += count > 0 ? T{} : excessmass;
    // divide excess mass up for all interface or fluid neighbors
    excessmass = count > 0 ? excessmass / count : T{};

    // update field values
    cell.template get<MASS<T>>() = mass;
    cell.template get<VOLUMEFRAC<T>>() = vof;
    cell.template get<MASSEX<T>>() = excessmass;
  }
};


// apply all the free surface dynamics
template <typename LATTICEMANAGERTYPE>
struct FreeSurfaceApply {
  using CELL = typename LATTICEMANAGERTYPE::CellType;
  using BLOCKLAT = typename LATTICEMANAGERTYPE::BLOCKLATTICE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(LATTICEMANAGERTYPE& latManager, std::int64_t count) {
    // mass transfer
    latManager.template ApplyInnerCellDynamics<surface_post_process0<CELL>>(count);

    latManager.template getField<STATE>().NormalCommunicate(count);

    latManager.template getField<MASS<T>>().AllCommunicate(count);

    // to fluid neighbor conversion
    latManager.template ApplyInnerCellDynamics<surface_post_process1<CELL>>(count);

    latManager.template getField<STATE>().NormalCommunicate(count);



    // gas to interface pop init and to gas neighbor conversion
    latManager.template ApplyInnerCellDynamics<surface_post_process2<CELL>>(count);

    latManager.template getField<STATE>().NormalCommunicate(count);


    latManager.Communicate(count);


    // finalize
    latManager.template ApplyInnerCellDynamics<surface_post_process3<CELL>>(count);

    latManager.template getField<STATE>().NormalCommunicate(count);


    latManager.template getField<MASS<T>>().AllCommunicate(count);
    latManager.template getField<MASSEX<T>>().AllCommunicate(count);
    latManager.template getField<VOLUMEFRAC<T>>().AllCommunicate(count);
  }
};

}  // namespace fx3dfs
