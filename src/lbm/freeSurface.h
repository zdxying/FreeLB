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

// define FS parameters as single data stored in Data
struct Lonely_ThBase : public FieldBase<1> {};
struct VOF_Trans_ThBase : public FieldBase<1> {};
struct Surface_Tension_EnabledBase : public FieldBase<1> {};
struct Surface_Tension_ParameterBase : public FieldBase<1> {};

// lonely threshold in mass transfer
template <typename T>
using Lonely_Th = Data<T, Lonely_ThBase>;
// vof transition threshold
template <typename T>
using VOF_Trans_Th = Data<T, VOF_Trans_ThBase>;
// surface tension enabled
using Surface_Tension_Enabled = Data<bool, Surface_Tension_EnabledBase>;
// surface tension parameter
template <typename T>
using Surface_Tension_Parameter = Data<T, Surface_Tension_ParameterBase>;


template <typename T, typename LatSet>
using FSFIELDS = TypePack<STATE, MASS<T>, VOLUMEFRAC<T>, EXCESSMASS<T, LatSet::q>>;

template <typename T>
using FSPARAMS = TypePack<Lonely_Th<T>, VOF_Trans_Th<T>, Surface_Tension_Enabled, Surface_Tension_Parameter<T>>;

struct NbrInfo {
  bool fluid_nbr = false;
  bool gas_nbr = false;
  int interface_nbrs = 0;
};
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

  void getNbrInfo(std::size_t id, NbrInfo& nbrinfo) const {
    for (int i = 1; i < LatSet::q; ++i) {
      auto iflag = this->template getField<STATE>().get(id + this->Delta_Index[i]);
      if (util::isFlag(iflag, FSType::Fluid)) {
        nbrinfo.fluid_nbr = true;
      } else if (util::isFlag(iflag, FSType::Gas)) {
        nbrinfo.gas_nbr = true;
      } else if (util::isFlag(iflag, FSType::Interface)) {
        ++nbrinfo.interface_nbrs;
      }
    }
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

  T Lonely_Threshold;
  // vof transition threshold
  T VOF_Trans_Threshold;

  bool Surface_Tension_Enabled = false;
  T surface_tension_parameter;  // coefficient_factor * coefficient

 public:
  //  FSType::Solid, T{}, T{}, T{}
  template <typename INITVALUEPACK>
  FreeSurface2DManager(BlockLatticeManager<T, LatSet, TypePack>& lm, INITVALUEPACK& initvalues, 
  T lonely_th = 0.3, T vof_trans_th = 0.01, T surface_tension = false, T surface_tension_param = 0.0)
      : BlockLatticeManagerBase<T, LatSet, FSFIELDS<T, LatSet>>(lm.getGeo(), initvalues),
      NSLatMan(lm), Lonely_Threshold(lonely_th), VOF_Trans_Threshold(vof_trans_th), 
      Surface_Tension_Enabled(surface_tension), surface_tension_parameter(surface_tension_param) {
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
          std::size_t idn = id + latset::c<LatSet>(i) * block.getProjection();
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
        std::size_t idn = id + latset::c<LatSet>(i) * block.getProjection();
        Vector<T, LatSet::d> loc_t = block.getLoc_t(idn);
        if (block.IsInside(loc_t)){
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
  latman.template getField<MASS<T>>().forEach(latman.template getField<STATE>(), FSType::Fluid,
                  [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
  latman.template getField<MASS<T>>().forEach(latman.template getField<STATE>(), FSType::Interface,
                  [&](auto& field, std::size_t id) { field.SetField(id, T{0.5}); });

  latman.template getField<MASS<T>>().NormalCommunicate();
  #ifdef MPI_ENABLED
  latman.template getField<MASS<T>>().MPINormalCommunicate(count);
  #endif
  
  latman.template getField<VOLUMEFRAC<T>>().forEach(latman.template getField<STATE>(), FSType::Fluid,
                        [&](auto& field, std::size_t id) { field.SetField(id, T{1}); });
  latman.template getField<VOLUMEFRAC<T>>().forEach(latman.template getField<STATE>(), FSType::Interface, 
                  [&](auto& field, std::size_t id) {field.SetField(id, T{0.5});});

  latman.template getField<VOLUMEFRAC<T>>().NormalCommunicate();
  #ifdef MPI_ENABLED
  latman.template getField<VOLUMEFRAC<T>>().MPINormalCommunicate(count);
  #endif

  }

};


// free surface as a post process of NS Lattice

// functions for free surface

template <typename CELL>
typename CELL::FloatType getClampedVOF(CELL& cell) {
  using T = typename CELL::FloatType;
  return std::clamp(cell.template get<VOLUMEFRAC<T>>(), T{}, T{1});
}

template <typename CELL>
static bool hasNeighborType(CELL& cell, FSType fstype) {
  using LatSet = typename CELL::LatticeSet;
  for (int i = 1; i < LatSet::q; ++i) {
    if (util::isFlag(cell.template getField<STATE>().get(cell.getNeighborId(i)), fstype)) return true;
  }
  return false;
}

// Parker-Youngs normal 
// Parker-Youngs weights, a corrected version:
// |c|^2  1  2  3
// weight 4  2  1 
template <typename LatSet>
constexpr std::array<int, LatSet::q> Parker_YoungsWeights() {
  return make_Array<int, LatSet::q>([&](unsigned int i) {
    int weight = 8;
    if (latset::c<LatSet>(i)[0] != 0) weight /= 2;
    if (latset::c<LatSet>(i)[1] != 0) weight /= 2;
    if constexpr (LatSet::d == 3) {
      if (latset::c<LatSet>(i)[2] != 0) weight /= 2;
    }
    return weight;
  });
}

template <typename CELL>
void computeParker_YoungsNormal(CELL& cell, Vector<typename CELL::FloatType, CELL::LatticeSet::d>& normal){
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  for (int i = 1; i < LatSet::q; ++i) {
    CELL celln = cell.getNeighbor(i);
    T clampedvof = getClampedVOF(celln);
    normal -= Parker_YoungsWeights<LatSet>()[i] * latset::c<LatSet>(i) * clampedvof;
  }
}

// openlb offset calculation
template<typename T>
T offsetHelper2D(T volume, const std::vector<T>& sorted_normal) {
  T d2 = volume * sorted_normal[1] + 0.5 * sorted_normal[0];
  if(d2 >= sorted_normal[0]){
    return d2;
  }
  T d1 = std::sqrt(2. * sorted_normal[0] * sorted_normal[1] * volume);
  return d1;
}

// openlb
// A lot of magic numbers are happening here. Optimized algorithm taken from Moritz Lehmann
template<typename T>
T offsetHelper3D(T volume, const std::vector<T>& sorted_normal) {
  T sn0_plus_sn1 = sorted_normal[0] + sorted_normal[1];
  T sn0_times_sn1 = sorted_normal[0] * sorted_normal[1];
  T sn2_volume = sorted_normal[2] * volume;

  T min_sn0_plus_sn1_and_sn2 = std::min(sn0_plus_sn1, sorted_normal[2]);

  T d5 = sn2_volume + 0.5 * sn0_plus_sn1;
  if(d5 > min_sn0_plus_sn1_and_sn2 && d5 <= sorted_normal[2]){
    return d5;
  }

  T d2 = 0.5 * sorted_normal[0] + 0.28867513 * std::sqrt(std::max(0., 24. * sorted_normal[1] * sn2_volume - sorted_normal[0]*sorted_normal[0]) );

  if(d2 > sorted_normal[0] && d2 <= sorted_normal[1]){
    return d2;
  }
  // cubic root
  T d1  = std::cbrt(6.0 * sn0_times_sn1 * sn2_volume);
  if(d1 <= sorted_normal[0]){
    return d1;
  }

  T x3 = 81.0  * sn0_times_sn1 * (sn0_plus_sn1 - 2. * sn2_volume);
  T y3 = std::sqrt(std::max(0., 23328. * sn0_times_sn1*sn0_times_sn1*sn0_times_sn1 - x3*x3 ));
  T u3 = std::cbrt(x3*x3 + y3*y3);
  T d3 = sn0_plus_sn1 - (7.5595264 * sn0_times_sn1 + 0.26456684 * u3) * (1./std::sqrt(u3)) * std::sin(0.5235988 - 0.3333334 * std::atan(y3 / x3));
  if(d3 > sorted_normal[1] && d3 <= min_sn0_plus_sn1_and_sn2){
    return d3;
  }

  T t4 = 9. * std::pow(sn0_plus_sn1 + sorted_normal[2], 2) - 18.;
  T x4 = std::max(sn0_times_sn1 * sorted_normal[2] * (324. - 648. * volume), 1.1754944e-38);
  T y4 = std::sqrt(std::max(4. * t4*t4*t4 - x4*x4, 0.));
  T u4 = std::cbrt(x4*x4 + y4*y4);
  T d4  = 0.5 * (sn0_plus_sn1 + sorted_normal[2]) - (0.20998684 * t4 + 0.13228342 * u4) * (1./std::sqrt(u4)) * std::sin(0.5235988- 0.3333334 * std::atan(y4/x4));

  return d4;
}

// openlb
// A lot of magic numbers are happening here. Optimized algorithm taken from Moritz Lehmann
template<typename T, typename LatSet>
T offsetHelperOpt(T vol, const Vector<T,LatSet::d>& sn) {
  const T sn0_p_sn1 = sn[0] + sn[1];
  const T sn2_t_V = sn[2] * vol;

  if(sn0_p_sn1 <= 2. * sn2_t_V){
    return sn2_t_V + 0.5 * sn0_p_sn1;
  }

  const T sq_sn0 = std::pow(sn[0],2), sn1_6 = 6. * sn[1], v1 = sq_sn0 / sn1_6;

  if(v1 <= sn2_t_V && sn2_t_V < v1 + 0.5 * (sn[1]-sn[0])){
    return 0.5 *(sn[0] + std::sqrt(sq_sn0 + 8.0 * sn[1] * (sn2_t_V - v1)));
  }

  const T v6 = sn[0] * sn1_6 * sn2_t_V;
  if(sn2_t_V < v1){
    return std::cbrt(v6);
  }

  const T v3 = sn[2] < sn0_p_sn1 ? (std::pow(sn[2],2) * (3. * sn0_p_sn1 - sn[2]) + sq_sn0 *(sn[0] - 3.0 * sn[2]) + std::pow(sn[1],2)*(sn[1]-3.0 * sn[2])) / (sn[0] * sn1_6) : 0.5 * sn0_p_sn1;

  const T sq_sn0_sq_sn1 = sq_sn0 + std::pow(sn[1],2), v6_cb_sn0_sn1 = v6 - std::pow(sn[0],3) - std::pow(sn[1],3);

  const bool case34 = sn2_t_V < v3;
  const T a = case34 ? v6_cb_sn0_sn1 : 0.5 * (v6_cb_sn0_sn1 - std::pow(sn[2], 3));
  const T b = case34 ? sq_sn0_sq_sn1 : 0.5 * (sq_sn0_sq_sn1 + std::pow(sn[2], 2));
  const T c = case34 ? sn0_p_sn1 : 0.5;
  const T t = std::sqrt(std::pow(c,2) - b);
  return c - 2.0 * t * std::sin(0.33333334 * std::asin((std::pow(c,3) - 0.5 * a - 1.5 * b * c) / std::pow(t,3)));
}

// openlb 
template<typename T, typename LatSet>
T calculateCubeOffset(T volume, const Vector<T,LatSet::d>& normal) {
  std::vector<T> abs_normal(LatSet::d, T{});
  for(unsigned int i = 0; i < LatSet::d; i++){
    abs_normal[i] = std::abs(normal[i]);
  }

  T volume_symmetry = 0.5 - std::abs(volume - 0.5);

  std::sort(abs_normal.begin(), abs_normal.end());

  if constexpr (LatSet::d == 2) {
    abs_normal[0] = std::max(normal[0], 1e-5);
  } else if (LatSet::d == 3){
    abs_normal[0] = std::max(normal[0], 1e-12);
    abs_normal[1] = std::max(normal[1], 1e-12);
  }

  T d{};
  if constexpr (LatSet::d == 2) {
    d = offsetHelper2D<T>(volume_symmetry, abs_normal);
  } else if (LatSet::d == 3){
    d = offsetHelper3D<T>(volume_symmetry, abs_normal);
  }

  T sorted_normal_acc = 0;
  for(unsigned int i = 0; i < LatSet::d; i++){
    sorted_normal_acc += abs_normal[i];
  }

  return std::copysign(d - 0.5 * sorted_normal_acc, volume - 0.5);
}

// openlb Optimized version of calculateCubeOffset
template<typename T, typename LatSet>
T calculateCubeOffsetOpt(T volume, const Vector<T,LatSet::d>& normal) {
  Vector<T, LatSet::d> abs_normal = normal.getabs();

  T a_l1 = abs_normal[0] + abs_normal[1] + abs_normal[2];

  T volume_symmetry = 0.5 - std::abs(volume - 0.5);

  Vector<T,LatSet::d> sorted_normal{};
  sorted_normal[0] = std::min(std::min(abs_normal[0], abs_normal[1]), abs_normal[2]) / a_l1;
  sorted_normal[1] = 0.;
  sorted_normal[2] = std::max(std::max(abs_normal[0], abs_normal[1]), abs_normal[2]) / a_l1;

  sorted_normal[1] = std::max(1. - sorted_normal[0] - sorted_normal[2], 0.);

  T d = offsetHelperOpt<T,LatSet>(volume_symmetry, sorted_normal);

  return a_l1 * std::copysign(0.5 - d, volume - 0.5);
}

// openlb pivoted LU solver
template<typename T, size_t S>
std::array<T,S> solvePivotedLU(std::array<std::array<T,S>,S>& matrix, const std::array<T,S>& b, size_t N) {
  std::array<T,S> x;
  std::array<T,S> pivots;
  for(size_t i = 0; i < S; ++i){
    pivots[i] = i;
    x[i] = 0.;
  }

  N = std::min(N,S);

  for(size_t i = 0; i < N; ++i){

    T max = 0.;
    size_t max_index = i;

    for(size_t j = i; j < N; ++j){
      T abs = std::abs(matrix[pivots[j]][i]);
      if(abs > max){
        max_index = j;
        max = abs;
      }
    }

    if(max_index != i){
      size_t tmp_index = pivots[i];
      pivots[i] = pivots[max_index];
      pivots[max_index] = tmp_index;
    }

    for(size_t j = i + 1; j < N; ++j){
      matrix[pivots[j]][i] /= matrix[pivots[i]][i];

      for(size_t k = i + 1; k < N; ++k){

        matrix[pivots[j]][k] -= matrix[pivots[j]][i] * matrix[pivots[i]][k];
      }
    }
  }

  for(size_t i = 0; i  < N; ++i){
    x[i] = b[pivots[i]];

    for(size_t j = 0; j < i; ++j){
      x[i] -= matrix[pivots[i]][j] * x[j];
    }
  }

  for(size_t i = N; i > 0; --i){
    for(size_t j = i; j < N; ++j){
      x[i-1] -= matrix[pivots[i-1]][j] * x[j];
    }

    x[i-1] /= matrix[pivots[i-1]][i-1];
  }

  return x;
}


// openlb 
template <typename CELL>
typename CELL::FloatType ComputeCurvature2D(CELL& cell) {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  // Parker-Youngs Normal
  Vector<T, 2> normal{};
  computeParker_YoungsNormal(cell, normal);
  T norm = normal.getnorm();
  if (norm < 1e-6) return T{};
  normal /= norm;

  // Rotation matrix is
  // ( n1 | -n0 )
  // ( n0 |  n1 )

  // It is 2 because of the amount of fitting parameters. Not because of the dimension
  constexpr size_t S = 2;
  std::array<std::array<T, S>, S> lq_matrix;
  std::array<T, S> b_rhs;
  for(size_t i = 0; i < S; ++i){
    for(size_t j = 0; j < S; ++j){
      lq_matrix[i][j] = 0.;
    }
    b_rhs[i] = 0.;
  }

  // Offset for the plic correction
  T origin_offset{};
  T VOF = getClampedVOF(cell);
  origin_offset = calculateCubeOffset<T, LatSet>(VOF, normal);

  std::size_t healthy_interfaces = 0;

  for (int iPop=1; iPop < LatSet::q; ++iPop) {
    CELL cellC = cell.getNeighbor(iPop);
    if(!util::isFlag(cellC.template get<STATE>(), FSType::Interface) || !hasNeighborType(cellC, FSType::Gas)) {
      continue;
    }

    ++healthy_interfaces;

    T fill_level = getClampedVOF(cellC);

    T cube_offset = calculateCubeOffset<T, LatSet>(fill_level, normal);

    T x_pos = latset::c<LatSet>(iPop)[0];
    T y_pos = latset::c<LatSet>(iPop)[1];

    // Rotation
    T rot_x_pos = x_pos * normal[1] - y_pos * normal[0];
    T rot_y_pos = x_pos * normal[0] + y_pos * normal[1] + (cube_offset - origin_offset);

    T rot_x_pos_2 = rot_x_pos * rot_x_pos;
    T rot_x_pos_3 = rot_x_pos_2 * rot_x_pos;
    T rot_x_pos_4 = rot_x_pos_3 * rot_x_pos;

    lq_matrix[1][1] += rot_x_pos_2;
    lq_matrix[1][0] += rot_x_pos_3;
    lq_matrix[0][0] += rot_x_pos_4;

    b_rhs[0] += rot_x_pos_2*(rot_y_pos);
    b_rhs[1] += rot_x_pos*(rot_y_pos);
  }

  lq_matrix[0][1] = lq_matrix[1][0];

  // Thikonov regularization parameter
  T alpha{};
  for(size_t i = 0; i < LatSet::d; ++i){
    lq_matrix[i][i] += alpha;
  }

  // It is 2 because of the fitting parameters. Not dependent on the dimension
  std::array<T,S> solved_fit = solvePivotedLU<T,S>(lq_matrix, b_rhs, healthy_interfaces);

  // signed curvature -> kappa = y'' / ( (1 + y'²)^(3/2) )
  T denom = std::sqrt(1. + solved_fit[1]*solved_fit[1]);
  denom = denom * denom * denom;
  T curvature = 2.*solved_fit[0] / denom;
  return std::max(-1., std::min(1., curvature));
}


// openlb 
template <typename CELL>
typename CELL::FloatType ComputeCurvature3D(CELL& cell){
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  // This is b_z
  // Parker-Youngs Normal
  Vector<T, 3> normal{};
  computeParker_YoungsNormal(cell, normal);
  T norm = normal.getnorm();
  if (norm < 1e-6) return T{};
  normal /= norm;

  std::array<T,3> r_vec{0.56270900, 0.32704452, 0.75921047};
  /*
  std::array<T,DESCRIPTOR::d> r_vec{
    0.,0.,1.
  };
  */
  std::array<std::array<T,3>,3> rotation{{
    {{0., 0., 0.}},
    //{{normal[1], -normal[0], 0.}},
    {{normal[1] * r_vec[2] - normal[2] * r_vec[1], normal[2] * r_vec[0] - normal[0] * r_vec[2], normal[0] * r_vec[1] - normal[1] * r_vec[0]}},
    {{normal[0], normal[1], normal[2]}}
  }};

  // Cross product with (0,0,1) x normal
  // This is b_y

  // (normal[0], normal[1], normal[2])

  T cross_norm = 0.;
  for(size_t i = 0; i < LatSet::d; ++i){
    cross_norm += rotation[1][i] * rotation[1][i];
  }

  // If too close too each other use the identity matrix
  if(cross_norm > 1e-6){
    cross_norm = std::sqrt(cross_norm);
    for(size_t i = 0; i <LatSet::d; ++i){
      rotation[1][i] /= cross_norm;
    }
  }else {
    rotation[1] = {{
      -normal[2],
      0.,
      normal[0]
    }};

    cross_norm = 0.;
    for(size_t i = 0; i < LatSet::d; ++i){
      cross_norm += rotation[1][i] * rotation[1][i];
    }

    cross_norm = std::sqrt(cross_norm);

    for(size_t i = 0; i <LatSet::d; ++i){
      rotation[1][i] /= cross_norm;
    }
  }

  // Cross product of ((0,0,1) x normal / | (0,0,1) x normal |) x normal
  // This is b_x
  rotation[0] = {{
    rotation[1][1] * normal[2] - rotation[1][2] * normal[1],
    rotation[1][2] * normal[0] - rotation[1][0] * normal[2],
    rotation[1][0] * normal[1] - rotation[1][1] * normal[0]
  }};

  // These three form a matrix and are entered into each row
  // ( b_x )
  // ( b_y )
  // ( b_z )

  constexpr size_t S = 5;
  std::array<std::array<T,S>, S> lq_matrix;
  std::array<T,S> b_rhs;
  for(size_t i = 0; i < S; ++i){
    for(size_t j = 0; j < S; ++j){
      lq_matrix[i][j] = 0.;
    }
    b_rhs[i] = 0.;
  }
  T origin_offset{};
  {
    T fill_level = getClampedVOF(cell);
    origin_offset = calculateCubeOffsetOpt<T,LatSet>(fill_level, normal);
  }

  size_t healthy_interfaces = 0;
  for(int iPop = 1; iPop < LatSet::q; iPop++){
    auto cellC = cell.getNeighbor(iPop);

    if(!util::isFlag(cellC.template get<STATE>(), FSType::Interface) || !hasNeighborType(cellC, FSType::Gas)) {
      continue;
    }

    ++healthy_interfaces;

    T fill_level = getClampedVOF(cellC);

    T cube_offset = calculateCubeOffsetOpt<T, LatSet>(fill_level, normal);

    int i = latset::c<LatSet>(iPop)[0];
    int j = latset::c<LatSet>(iPop)[1];
    int k = latset::c<LatSet>(iPop)[2];

    std::array<T,3> pos{static_cast<T>(i),static_cast<T>(j),static_cast<T>(k)};
    std::array<T,3> r_pos{0.,0.,cube_offset - origin_offset};

    for(size_t a = 0; a < LatSet::d; ++a){
      for(size_t b = 0; b < LatSet::d; ++b){
        r_pos[a] += rotation[a][b] * pos[b];
      }
    }

    T r_x_2 = r_pos[0] * r_pos[0];
    T r_x_3 = r_x_2 * r_pos[0];
    T r_x_4 = r_x_3 * r_pos[0];

    T r_y_2 = r_pos[1] * r_pos[1];
    T r_y_3 = r_y_2 * r_pos[1];
    T r_y_4 = r_y_3 * r_pos[1];

    T r_x_2_y_2 = r_x_2 * r_y_2;
    T r_x_3_y = r_x_3 * r_pos[1];
    T r_x_2_y = r_x_2 * r_pos[1];

    T r_x_y_3 = r_pos[0] * r_y_3;
    T r_x_y_2 = r_pos[0] * r_y_2;

    T r_x_y = r_pos[0] * r_pos[1];

    lq_matrix[0][0] += r_x_4;
    lq_matrix[1][1] += r_y_4;
    lq_matrix[2][2] += r_x_2_y_2;
    lq_matrix[3][3] += r_x_2;
    lq_matrix[4][4] += r_y_2;

    // skip [1][0] copy later from [2][2]
    lq_matrix[2][0] += r_x_3_y;
    lq_matrix[3][0] += r_x_3;
    lq_matrix[4][0] += r_x_2_y;

    lq_matrix[2][1] += r_x_y_3;
    lq_matrix[3][1] += r_x_y_2;
    lq_matrix[4][1] += r_y_3;

    // skip [3][2] copy from [4][0]
    // skip [4][2] copy from [3][1]

    lq_matrix[4][3] += r_x_y;

    b_rhs[0] +=  r_x_2 * r_pos[2];
    b_rhs[1] +=  r_y_2 * r_pos[2];
    b_rhs[2] +=  r_x_y * r_pos[2];
    b_rhs[3] +=  r_pos[0] * r_pos[2];
    b_rhs[4] +=  r_pos[1] * r_pos[2];
  }

  lq_matrix[1][0] = lq_matrix[2][2];
  lq_matrix[3][2] = lq_matrix[4][0];
  lq_matrix[4][2] = lq_matrix[3][1];

  for(size_t i = 0; i < S; ++i){
    for(size_t j = i + 1; j < S; ++j){
      lq_matrix[i][j] = lq_matrix[j][i];
    }
  }

  // Consider using Thikonov regularization?
  //T alpha = 1e-8;
  T alpha{};
  for(size_t i = 0; i < S; ++i){
    lq_matrix[i][i] += alpha;
  }

  std::array<T,S> solved_fit = solvePivotedLU<T,S>(lq_matrix, b_rhs, healthy_interfaces);

  T denom = std::sqrt(1. + solved_fit[3]*solved_fit[3] + solved_fit[4]*solved_fit[4]);
  denom = denom * denom * denom;
  T curvature = ( (1.+solved_fit[4]*solved_fit[4]) * solved_fit[0] + (1. + solved_fit[3]*solved_fit[3] ) * solved_fit[1] - solved_fit[3] * solved_fit[4] * solved_fit[2] ) / denom;

  return std::max(-1., std::min(1., curvature));
}


template <typename CELL>
typename CELL::FloatType ComputeCurvature(CELL& cell) {
  using LatSet = typename CELL::LatticeSet;
  if constexpr (LatSet::d == 2) {
    return ComputeCurvature2D(cell);
  } else if constexpr (LatSet::d == 3) {
    return ComputeCurvature3D(cell);
  }
}


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

    NbrInfo(CELL& cell){
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
    if (util::isFlag(cell.template get<STATE>(), FSType::Interface)) {
      // mass transfer
      T deltamass{};
      T massflow{};
      // cell's nbr info
      NbrInfo cellNbrInfo(cell);

      for (int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        int kopp = latset::opp<LatSet>(k);

        if (util::isFlag(celln.template get<STATE>(), FSType::Fluid)) {
          deltamass += cell[kopp] - celln[k];
        } else if (util::isFlag(celln.template get<STATE>(), FSType::Interface)) {
          // celln's nbr info
          NbrInfo cellnNbrInfo(celln);
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

          deltamass += massflow * T(0.5) * (getClampedVOF(cell) + getClampedVOF(celln));
        }

      }
      
      cell.template get<MASS<T>>() += deltamass;

      // reconstruct pop streamed in from a gas cell
      T curvature{};
      if (cell.template get<Surface_Tension_Enabled>()) {
        if (hasNeighborType(cell, FSType::Gas)) curvature = ComputeCurvature(cell);
      }
      T rho_gas = T(1) - T(6) * cell.template get<Surface_Tension_Parameter<T>>() * curvature;
      const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T,LatSet::d>>();
      T u2 = u.getnorm2();
      for (int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        if (util::isFlag(celln.template get<STATE>(), FSType::Gas)) {
          // fiopp = feqiopp(rho_gas) + feqi(rho_gas) - fi(x+ei)
          cell[latset::opp<LatSet>(k)] =
            equilibrium::SecondOrder<CELL>::get(k, u, rho_gas, u2) +
            equilibrium::SecondOrder<CELL>::get(latset::opp<LatSet>(k), u, rho_gas, u2) 
            - celln[k];
        }
      }

      // transition flag for interface cell
      T rho = moment::rho<CELL>::get(cell);

      // transition by mass criterion
      if (cell.template get<MASS<T>>() > (T(1) + cell.template get<VOF_Trans_Th<T>>()) * rho) {
        util::addFlag(FSType::To_Fluid, util::underlyingRef(cell.template get<STATE>()));
        return;
      } else if (cell.template get<MASS<T>>() < -cell.template get<VOF_Trans_Th<T>>() * rho) {
        util::addFlag(FSType::To_Gas, util::underlyingRef(cell.template get<STATE>()));
        return;
      }
      // transition by lonely criterion
      if (cell.template get<MASS<T>>() > (T(1) - cell.template get<Lonely_Th<T>>()) * rho) {
        if (!hasNeighborType(cell, FSType::Gas)) {
          util::addFlag(FSType::To_Fluid, util::underlyingRef(cell.template get<STATE>()));
          return;
        }
      } else if (cell.template get<MASS<T>>() < cell.template get<Lonely_Th<T>>() * rho) {
        if (!hasNeighborType(cell, FSType::Fluid)) {
          util::addFlag(FSType::To_Gas, util::underlyingRef(cell.template get<STATE>()));
          return;
        }
      }
      // deal with isolated interface cells
      if (!hasNeighborType(cell, FSType::Interface)) {
        if (!hasNeighborType(cell, FSType::Fluid)) {
          util::addFlag(FSType::To_Gas, util::underlyingRef(cell.template get<STATE>()));
        } else if (!hasNeighborType(cell, FSType::Gas)) {
          util::addFlag(FSType::To_Fluid, util::underlyingRef(cell.template get<STATE>()));
        }
      }
    }
  }
};

// to fluid neighbor conversion
template <typename CELLTYPE>
struct ToFluidNbrConversion {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    // if (util::isFlag(cell.template get<STATE>(), FSType::To_Fluid)) {
    //   // check neighbors
    //   for (int k = 1; k < LatSet::q; ++k) {
    //     CELL celln = cell.getNeighbor(k);
    //     if (util::isFlag(celln.template get<STATE>(), FSType::To_Gas)) {
    //       // remove to_gas flag
    //       util::removeFlag(FSType::To_Gas, util::underlyingRef(celln.template get<STATE>()));
    //     } else if (util::isFlag(celln.template get<STATE>(), FSType::Gas)) {
    //       // set to_interface for gas neighbor cells
    //       util::addFlag(FSType::To_Interface, util::underlyingRef(celln.template get<STATE>()));
    //     }
    //   }
    // }

    if (util::isFlag(cell.template get<STATE>(), FSType::To_Gas)) {
      if (hasNeighborType(cell, FSType::To_Fluid)){
        // remove to_gas flag
        util::removeFlag(FSType::To_Gas, util::underlyingRef(cell.template get<STATE>()));
      }
    } else if (util::isFlag(cell.template get<STATE>(), FSType::Gas)) {
      if (hasNeighborType(cell, FSType::To_Fluid)){
        // set to_interface
        util::addFlag(FSType::To_Interface, util::underlyingRef(cell.template get<STATE>()));
      }
    }
  }

};

// gas to interface pop init
template <typename CELLTYPE>
struct GasToInterfacePopInit {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    if (util::isFlag(cell.template get<STATE>(), FSType::To_Interface)) {
      // init fi using [fluid|interface] neighbor cells
      T averho{};
      Vector<T, LatSet::d> aveu{};
      int count{};
      for (int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        if (util::isFlag(celln.template get<STATE>(), (FSType::Fluid | FSType::Interface))) {
          averho += moment::rho<CELL>::get(celln);
          aveu += moment::u<CELL>::get(celln);
          ++count;
        }
      }
      averho /= count;
      aveu /= count;
      // set fi
      T aveu2 = aveu.getnorm2();
      for (int k = 0; k < LatSet::q; ++k) {
        cell[k] = equilibrium::SecondOrder<CELL>::get(k, aveu, averho, aveu2);
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
    // if (util::isFlag(cell.template get<STATE>(), FSType::To_Gas)) {
    //   // check neighbors
    //   for (int k = 1; k < LatSet::q; ++k) {
    //     CELL celln = cell.getNeighbor(k);
    //     if (util::isFlag(celln.template get<STATE>(), FSType::Fluid)) {
    //       // to interface
    //       util::addFlag(FSType::To_Interface, util::underlyingRef(celln.template get<STATE>()));
    //     }
    //   }
    // }

    if (util::isFlag(cell.template get<STATE>(), FSType::Fluid)) {
      if (hasNeighborType(cell, FSType::To_Gas)){
        // set to_interface
        util::addFlag(FSType::To_Interface, util::underlyingRef(cell.template get<STATE>()));
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
    T excessmass{};
    if (util::isFlag(cell.template get<STATE>(), FSType::To_Fluid)) {
      T rho = moment::rho<CELL>::get(cell);
      excessmass = cell.template get<MASS<T>>() - rho;
      cell.template get<MASS<T>>() = rho;
    } else if (util::isFlag(cell.template get<STATE>(), FSType::To_Gas)) {
      excessmass = cell.template get<MASS<T>>();
      cell.template get<MASS<T>>() = T{};
    } else {
      return;
    }

    // find neighbors
    int count{};
    for (int k = 1; k < LatSet::q; ++k) {
      CELL celln = cell.getNeighbor(k);
      if (util::isFlag(celln.template get<STATE>(), FSType::Interface) &&
          !util::isFlag(celln.template get<STATE>(), FSType::To_Gas | FSType::To_Fluid)) {
        ++count;
      }
    }

    std::array<T*, LatSet::q> exmasscell = 
    cell.template getField<EXCESSMASS<T,LatSet::q>>().getArray(cell.getId());
    *(exmasscell[0]) = T{};
    if (count > 0) {
      T excessmassk = excessmass / count;
      for (int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        if (util::isFlag(celln.template get<STATE>(), FSType::Interface) &&
            !util::isFlag(celln.template get<STATE>(), FSType::To_Gas | FSType::To_Fluid)) {
          *(exmasscell[k]) = excessmassk;
        }
      }
    } else {
      *(exmasscell[0]) = excessmass;
    }
  }

};

// finalize conversion
template <typename CELLTYPE>
struct FinalizeConversion {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    // update state
    if (util::isFlag(cell.template get<STATE>(), FSType::To_Fluid)) {
      cell.template get<STATE>() = FSType::Fluid;
      cell.template get<VOLUMEFRAC<T>>() = T(1);
      cell.template get<MASS<T>>() += cell.template get<EXCESSMASS<T,LatSet::q>, 0>();
    } else if (util::isFlag(cell.template get<STATE>(), FSType::To_Gas)) {
      cell.template get<STATE>() = FSType::Gas;
      cell.template get<VOLUMEFRAC<T>>() = T{};
      cell.template get<MASS<T>>() += cell.template get<EXCESSMASS<T,LatSet::q>, 0>();
      cell.template get<VELOCITY<T,LatSet::d>>().clear();
    } else if (util::isFlag(cell.template get<STATE>(), FSType::To_Interface)) {
      cell.template get<STATE>() = FSType::Interface;
    }
  }

};
// stream EXCESSMASS<T,LatSet::q>
template <typename LATTICETYPE>
struct StreamExcessMass {
  using LATTICE = LATTICETYPE;
  using T = typename LATTICE::FloatType;
  using LatSet = typename LATTICE::LatticeSet;

  static void apply(LATTICE& lattice){
    for (int i = 1; i < LatSet::q; ++i) {
      lattice.template getField<EXCESSMASS<T,LatSet::q>>().getField(i).rotate(lattice.getDelta_Index()[i]);
    }
  }

};

// collect excess mass
template <typename CELLTYPE>
struct CollectExcessMass {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    if (util::isFlag(cell.template get<STATE>(), FSType::Interface | FSType::Fluid)) {
        T exmass_sum = cell.template get<EXCESSMASS<T,LatSet::q>, 0>();
      for (int k = 1; k < LatSet::q; ++k) {
        exmass_sum += cell.template get<EXCESSMASS<T,LatSet::q>>(k);
      }
      cell.template get<MASS<T>>() += exmass_sum;
      if (util::isFlag(cell.template get<STATE>(), FSType::Interface)) {
        T rho = moment::rho<CELL>::get(cell);
        cell.template get<VOLUMEFRAC<T>>() = cell.template get<MASS<T>>() / rho;
      }
    }
  }

}; 

// clear EXCESSMASS<T,LatSet::q>
template <typename LATTICETYPE>
struct ClearExcessMass {
  using LATTICE = LATTICETYPE;
  using T = typename LATTICE::FloatType;
  using LatSet = typename LATTICE::LatticeSet;

  static void apply(LATTICE& lattice){
    lattice.template getField<EXCESSMASS<T,LatSet::q>>().Init();
  }
  
};

// apply all the free surface dynamics
template <typename LATTICEMANAGERTYPE>
struct FreeSurfaceApply{
  using CELL = typename LATTICEMANAGERTYPE::CellType;
  using BLOCKLAT = typename LATTICEMANAGERTYPE::BLOCKLATTICE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void Apply(LATTICEMANAGERTYPE& latManager, std::int64_t count) {

    // mass transfer
    latManager.template ApplyInnerCellDynamics<MassTransfer<CELL>>(count);

    latManager.template getField<STATE>().NormalCommunicate(count);
    #ifdef MPI_ENABLED
    latManager.template getField<STATE>().MPINormalCommunicate(count);
    #endif
    latManager.template getField<MASS<T>>().CommunicateAll(count);


    // to fluid neighbor conversion
    latManager.template ApplyInnerCellDynamics<ToFluidNbrConversion<CELL>>(count);

    latManager.template getField<STATE>().NormalCommunicate(count);
    #ifdef MPI_ENABLED
    latManager.template getField<STATE>().MPINormalCommunicate(count);
    #endif


    // gas to interface pop init
    latManager.template ApplyInnerCellDynamics<GasToInterfacePopInit<CELL>>(count);

    latManager.Communicate(count);


    // to gas neighbor conversion
    latManager.template ApplyInnerCellDynamics<ToGasNbrConversion<CELL>>(count);

    latManager.template getField<STATE>().NormalCommunicate(count);
    #ifdef MPI_ENABLED
    latManager.template getField<STATE>().MPINormalCommunicate(count);
    #endif


    // interface excess mass
    latManager.template ApplyInnerCellDynamics<InterfaceExcessMass<CELL>>(count);

    latManager.template getField<MASS<T>>().CommunicateAll(count);
    latManager.template getField<EXCESSMASS<T, LatSet::q>>().CommunicateAll(count);


    // finalize conversion
    latManager.template ApplyInnerCellDynamics<FinalizeConversion<CELL>>(count);

    latManager.template getField<STATE>().NormalCommunicate(count);
    #ifdef MPI_ENABLED
    latManager.template getField<STATE>().MPINormalCommunicate(count);
    #endif
    latManager.template getField<MASS<T>>().CommunicateAll(count);
    latManager.template getField<VOLUMEFRAC<T>>().CommunicateAll(count);
    // latManager.template getField<VELOCITY<T,LatSet::d>>().CommunicateAll(count);


    // stream EXCESSMASS<T,LatSet::q>
    latManager.ForEachBlockLattice(count, [&](auto& blocklat){StreamExcessMass<BLOCKLAT>::apply(blocklat);});


    // collect excess mass
    latManager.template ApplyInnerCellDynamics<CollectExcessMass<CELL>>(count);

    latManager.template getField<MASS<T>>().CommunicateAll(count);
    latManager.template getField<VOLUMEFRAC<T>>().CommunicateAll(count);


    // clear EXCESSMASS<T,LatSet::q>
    latManager.ForEachBlockLattice(count, [&](auto& blocklat){ClearExcessMass<BLOCKLAT>::apply(blocklat);});
    

  }
};

}  // namespace FS
