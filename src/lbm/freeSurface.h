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

// free surface model
#pragma once

#include "data_struct/block_lattice.h"

// namespace FreeSurface

namespace FS {

enum FSType : std::uint8_t {
  Void = 0,
  Gas = 1,
  Interface = 2,
  Fluid = 4,
  To_Fluid = 8,
  To_Gas = 16,
  To_Interface = 32
};

template <typename T, typename LatSet>
T ComputeCurvature(BCell<T, LatSet>& cell);

template <typename T, typename LatSet>
class FreeSurface2D {
 private:
  // interface cells
  std::vector<std::size_t> Interface;

  BlockLattice<T, LatSet>& NS;
  ScalerField<FSType>& State;
  ScalerField<T>& Mass;
  PopulationField<T, LatSet::q>& ExcessMass;
  ScalerField<T>& VolumeFrac;

  T Lonely_Threshold;
  T VOF_Trans_Threshold;

  bool Surface_Tension_Enabled;
  T surface_tension_parameter;  // coefficient_factor * coefficient

 public:
  FreeSurface2D(BlockLattice<T, LatSet>& ns, ScalerField<FSType>& type,
                ScalerField<T>& mass, PopulationField<T, LatSet::q>& exmass,
                ScalerField<T>& vf, T lth, T vtth)
      : NS(ns), State(type), Mass(mass), ExcessMass(exmass), VolumeFrac(vf),
        Lonely_Threshold(lth), VOF_Trans_Threshold(vtth) {}


  inline bool hasNeighborType(std::size_t id, FSType fstype) const;

  inline T getClampedVOF(std::size_t id) const {
    return std::clamp(VolumeFrac.get(id), T(0), T(1));
  }

  // mass transfer/ advection
  void MassTransfer();
  // prepare for To_Fluid conversion
  void ToFluidNbrConversion();
  // prepare for To_Interface conversion
  void GasToInterfacePopInit();
  // prepare for To_Gas conversion
  void ToGasNbrConversion();
  // excess mass
  void InterfaceExcessMass();
  // finalize conversion
  void FinalizeConversion();
  // collect excess mass
  void CollectExcessMass();
};

template <typename T, typename LatSet>
class FreeSurface2DManager {
 private:
  std::vector<FreeSurface2D<T, LatSet>> BlockFS;
  // free surface state
  BlockFieldManager<ScalerField<FSType>, T, 2> StateFM;
  // mass = rho * volumefraction
  BlockFieldManager<ScalerField<T>, T, 2> MassFM;
  // Excess mass
  BlockFieldManager<PopulationField<T, LatSet::q>, T, 2> ExcessMassFM;
  // fill level/ volume fraction in VOF
  BlockFieldManager<ScalerField<T>, T, 2> VolumeFracFM;

  BlockLatticeManager<T, LatSet>& LatMan;

  T Lonely_Threshold = 0.1;
  // vof transition threshold
  T VOF_Trans_Threshold = 0.003;

  bool Surface_Tension_Enabled = false;
  T surface_tension_parameter;  // coefficient_factor * coefficient

 public:
  FreeSurface2DManager(BlockLatticeManager<T, LatSet>& lm)
      : LatMan(lm), StateFM(lm.getGeo(), FSType::Void), MassFM(lm.getGeo(), T{}),
        ExcessMassFM(lm.getGeo(), T{}), VolumeFracFM(lm.getGeo(), T{}) {
    // init FreeSurface2D
    for (int i = 0; i < LatMan.getGeo().getBlockNum(); ++i) {
      BlockFS.emplace_back(
        lm.getBlockLat(i), StateFM.getBlockField(i).getField(),
        MassFM.getBlockField(i).getField(), ExcessMassFM.getBlockField(i).getField(),
        VolumeFracFM.getBlockField(i).getField(), Lonely_Threshold, VOF_Trans_Threshold);
    }
  }

  // get
  BlockFieldManager<ScalerField<FSType>, T, 2>& getStateFM() { return StateFM; }
  BlockFieldManager<ScalerField<T>, T, 2>& getMassFM() { return MassFM; }
  BlockFieldManager<PopulationField<T, LatSet::q>, T, 2>& getExcessMassFM() {
    return ExcessMassFM;
  }
  BlockFieldManager<ScalerField<T>, T, 2>& getVolumeFracFM() { return VolumeFracFM; }

  void setSrufaceTension(T value) {
    if (value > 1e-3)
      std::cout << "Surface tension: " << value << " may be too large" << std::endl;
    surface_tension_parameter = value;
    Surface_Tension_Enabled = true;
  }

  void Init() {
    // set interface cells
    StateFM.forEach([&](auto& blockfield, std::size_t id) {
      auto& field = blockfield.getField();
      const auto& block = blockfield.getBlock();
      if (util::isFlag(field.get(id), FSType::Fluid)) {
        for (int i = 1; i < LatSet::q; ++i) {
          std::size_t idn = id + LatSet::c[i] * block.getProjection();
          if (util::isFlag(field.get(idn), FSType::Gas)) {
            util::removeFlag(FSType::Gas, field.getField(0).getUnderlying(idn));
            util::addFlag(FSType::Interface, field.getField(0).getUnderlying(idn));
          }
        }
      }
    });
    // set mass and volume fraction
    MassFM.forEach(StateFM, FSType::Fluid,
                   [&](auto& field, std::size_t id) { field.SetField(id, T(1)); });
    MassFM.forEach(StateFM, FSType::Interface,
                   [&](auto& field, std::size_t id) { field.SetField(id, T(0.5)); });
    VolumeFracFM.forEach(StateFM, FSType::Fluid,
                         [&](auto& field, std::size_t id) { field.SetField(id, T(1)); });
    VolumeFracFM.forEach(StateFM, FSType::Interface, [&](auto& field, std::size_t id) {
      field.SetField(id, T(0.5));
    });
  }
  void Apply() {
    for (FreeSurface2D<T, LatSet>& fs : BlockFS) {
      fs.MassTransfer();
      // for cells with to_fluid flag, check neighbors and set transition flag
      fs.ToFluidNbrConversion();
      // for (gas) cells with to_interface flag, init pop using nbr fluid/interface cells
      fs.GasToInterfacePopInit();
      // for cells with to_gas flag, check neighbors and set transition flag
      fs.ToGasNbrConversion();
      // excess mass
      fs.InterfaceExcessMass();
      fs.FinalizeConversion();
      fs.CollectExcessMass();
    }
  }
};

}  // namespace FS
