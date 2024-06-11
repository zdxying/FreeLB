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
  Solid = 1,
  Wall = 2,
  Gas = 4,
  Interface = 8,
  Fluid = 16,
  To_Fluid = 32,
  To_Gas = 64,
  To_Interface = 128
};

template <typename CELL>
typename CELL::FloatType ComputeCurvature(CELL& cell){
  using T = typename CELL::FloatType;
  return T{};
}

// define unique FS Field
struct STATEBase : public FieldBase<1> {};
struct MASSBase : public FieldBase<1> {};
struct VOLUMEFRACBase : public FieldBase<1> {};
template <unsigned int q>
struct EXCESSMASSBase : public FieldBase<q> {};

// free surface state, init with Solid
using STATE = GenericField<GenericArray<FSType>, STATEBase>;
// mass = rho * volumefraction
template <typename T>
using MASS = GenericField<GenericArray<T>, MASSBase>;
// fill level/ volume fraction in VOF
template <typename T>
using VOLUMEFRAC = GenericField<GenericArray<T>, VOLUMEFRACBase>;
// Excess mass
template <typename T, unsigned int q>
using EXCESSMASS = GenericField<CyclicArray<T>, EXCESSMASSBase<q>>;

// define FS parameters as single data stored in Array
struct Lonely_ThBase : public FieldBase<1> {};
struct VOF_Trans_ThBase : public FieldBase<1> {};
struct Surface_Tension_EnabledBase : public FieldBase<1> {};
struct Surface_Tension_ParameterBase : public FieldBase<1> {};

// lonely threshold in mass transfer
template <typename T>
using Lonely_Th = Array<T, Lonely_ThBase>;
// vof transition threshold
template <typename T>
using VOF_Trans_Th = Array<T, VOF_Trans_ThBase>;
// surface tension enabled
using Surface_Tension_Enabled = Array<bool, Surface_Tension_EnabledBase>;
// surface tension parameter
template <typename T>
using Surface_Tension_Parameter = Array<T, Surface_Tension_ParameterBase>;


template <typename T, typename LatSet>
using FSFIELDS = TypePack<STATE, MASS<T>, VOLUMEFRAC<T>, EXCESSMASS<T, LatSet::q>>;

template <typename T>
using FSPARAMS = TypePack<Lonely_Th<T>, VOF_Trans_Th<T>, Surface_Tension_Enabled, Surface_Tension_Parameter<T>>;

// incorporate free surface fields into BlockLattice is possible
// (considered as a post-process of NS Lattice)
// which may be done in the future
template <typename T, typename LatSet, typename TypePack>
class FreeSurface2D : public BlockLatticeBase<T, LatSet, FSFIELDS<T, LatSet>>{
 private:
  // interface cells
  std::vector<std::size_t> Interface;

  BlockLattice<T, LatSet, TypePack>& NS;

  T Lonely_Threshold;
  T VOF_Trans_Threshold;

  bool Surface_Tension_Enabled;
  T surface_tension_parameter;  // coefficient_factor * coefficient

 public:
 template <typename... FIELDPTRS>
  FreeSurface2D(BlockLattice<T, LatSet, TypePack>& ns, std::tuple<FIELDPTRS...> fieldptrs, 
  T lth, T vtth)
  : BlockLatticeBase<T, LatSet, FSFIELDS<T, LatSet>>(ns.getGeo(), fieldptrs),
  NS(ns), Lonely_Threshold(lth), VOF_Trans_Threshold(vtth) {}

  inline bool hasNeighborType(std::size_t id, FSType fstype) const{
  for (int i = 1; i < LatSet::q; ++i) {
    if (util::isFlag(this->template getField<STATE>().get(id + this->Delta_Index[i]), fstype)) return true;
  }
  return false;
  }

  inline T getClampedVOF(std::size_t id) const {
    return std::clamp(this->template getField<VOLUMEFRAC<T>>().get(id), T(0), T{1});
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

template <typename T, typename LatSet, typename TypePack>
class FreeSurface2DManager : public BlockLatticeManagerBase<T, LatSet, FSFIELDS<T, LatSet>>{
 private:
  std::vector<FreeSurface2D<T, LatSet, TypePack>> BlockFS;

  BlockLatticeManager<T, LatSet, TypePack>& NSLatMan;

  T Lonely_Threshold = 0.1;
  // vof transition threshold
  T VOF_Trans_Threshold = 0.003;

  bool Surface_Tension_Enabled = false;
  T surface_tension_parameter;  // coefficient_factor * coefficient

 public:
  //  FSType::Solid, T{}, T{}, T{}
  template <typename INITVALUEPACK>
  FreeSurface2DManager(BlockLatticeManager<T, LatSet, TypePack>& lm, INITVALUEPACK& initvalues)
      : BlockLatticeManagerBase<T, LatSet, FSFIELDS<T, LatSet>>(lm.getGeo(), initvalues),
      NSLatMan(lm){
    // init FreeSurface2D
    for (int i = 0; i < this->BlockGeo.getBlockNum(); ++i){
      BlockFS.emplace_back(
        NSLatMan.getBlockLat(i), 
        ExtractFieldPtrs<T, LatSet, FSFIELDS<T, LatSet>>::getFieldPtrTuple(i, this->Fields), 
        Lonely_Threshold, VOF_Trans_Threshold);
    }
  }

  void setSrufaceTension(T value) {
    if (value > 1e-3)
      std::cout << "Surface tension: " << value << " may be too large" << std::endl;
    surface_tension_parameter = value;
    Surface_Tension_Enabled = true;
  }

  void Init() {
    // set interface cells
    this->template getField<STATE>().forEach([&](auto& field, std::size_t id) {
      const auto& block = field.getBlock();
      if (util::isFlag(field.get(id), FSType::Fluid)) {
        for (int i = 1; i < LatSet::q; ++i) {
          std::size_t idn = id + LatSet::c[i] * block.getProjection();
          if (util::isFlag(field.get(idn), FSType::Gas)) {
            util::removeFlag(FSType::Gas, util::underlyingRef(field.get(idn)));
            util::addFlag(FSType::Interface, util::underlyingRef(field.get(idn)));
          }
        }
      }
    });
    // set mass and volume fraction
    this->template getField<MASS<T>>().forEach(this->template getField<STATE>(), FSType::Fluid,
                   [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
    this->template getField<MASS<T>>().forEach(this->template getField<STATE>(), FSType::Interface,
                   [&](auto& field, std::size_t id) { field.SetField(id, T{0.5}); });
    
    this->template getField<VOLUMEFRAC<T>>().forEach(this->template getField<STATE>(), FSType::Fluid,
                         [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
    this->template getField<VOLUMEFRAC<T>>().forEach(this->template getField<STATE>(), FSType::Interface, [&](auto& field, std::size_t id) {
      field.SetField(id, T{0.5});
    });
  }

  void Apply() {
    for (auto& fs : BlockFS) {
      // int deLevel = static_cast<int>(this->getMaxLevel() - fs.getLevel());
      // if (count % (static_cast<int>(pow(2, deLevel))) == 0)
      
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


template <typename LATTICEMANTYPE>
struct FreeSurfaceHelper{
  using LATTICEMAN = LATTICEMANTYPE;
  using T = typename LATTICEMAN::FloatType;
  using LatSet = typename LATTICEMAN::LatticeSet;

  static void Init(LATTICEMAN& latman){
  // set interface cells
  latman.template getField<STATE>().forEach([&](auto& field, std::size_t id) {
    const auto& block = field.getBlock();
    if (util::isFlag(field.get(id), FSType::Fluid)) {
      for (int i = 1; i < LatSet::q; ++i) {
        std::size_t idn = id + LatSet::c[i] * block.getProjection();
        if (util::isFlag(field.get(idn), FSType::Gas)) {
          util::removeFlag(FSType::Gas, util::underlyingRef(field.get(idn)));
          util::addFlag(FSType::Interface, util::underlyingRef(field.get(idn)));
        }
      }
    }
  });
  // set mass and volume fraction
  latman.template getField<MASS<T>>().forEach(latman.template getField<STATE>(), FSType::Fluid,
                  [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
  latman.template getField<MASS<T>>().forEach(latman.template getField<STATE>(), FSType::Interface,
                  [&](auto& field, std::size_t id) { field.SetField(id, T{0.5}); });
  
  latman.template getField<VOLUMEFRAC<T>>().forEach(latman.template getField<STATE>(), FSType::Fluid,
                        [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
  latman.template getField<VOLUMEFRAC<T>>().forEach(latman.template getField<STATE>(), FSType::Interface, [&](auto& field, std::size_t id) {
    field.SetField(id, T{0.5});
  });
  }
};

}  // namespace FS
