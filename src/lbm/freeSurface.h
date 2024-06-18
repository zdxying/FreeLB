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
    if (LatSet::c[i][0] != 0) weight /= 2;
    if (LatSet::c[i][1] != 0) weight /= 2;
    if constexpr (LatSet::d == 3) {
      if (LatSet::c[i][2] != 0) weight /= 2;
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
    normal -= Parker_YoungsWeights<LatSet>()[i] * LatSet::c[i] * clampedvof;
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

template<typename T>
T offsetHelper3D(T volume, const std::vector<T>& sorted_normal) {
  return T{};
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

    T x_pos = LatSet::c[iPop][0];
    T y_pos = LatSet::c[iPop][1];

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
        int kopp = LatSet::opp[k];

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
          cell[LatSet::opp[k]] =
            equilibrium::SecondOrder<CELL>::get(k, u, rho_gas, u2) +
            equilibrium::SecondOrder<CELL>::get(LatSet::opp[k], u, rho_gas, u2) 
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
    if (util::isFlag(cell.template get<STATE>(), FSType::To_Fluid)) {
      // check neighbors
      for (int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        if (util::isFlag(celln.template get<STATE>(), FSType::To_Gas)) {
          // remove to_gas flag
          util::removeFlag(FSType::To_Gas, util::underlyingRef(celln.template get<STATE>()));
        } else if (util::isFlag(celln.template get<STATE>(), FSType::Gas)) {
          // set to_interface for gas neighbor cells
          util::addFlag(FSType::To_Interface, util::underlyingRef(celln.template get<STATE>()));
        }
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
    if (util::isFlag(cell.template get<STATE>(), FSType::To_Gas)) {
      // check neighbors
      for (int k = 1; k < LatSet::q; ++k) {
        CELL celln = cell.getNeighbor(k);
        if (util::isFlag(celln.template get<STATE>(), FSType::Fluid)) {
          // to interface
          util::addFlag(FSType::To_Interface, util::underlyingRef(celln.template get<STATE>()));
        }
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

}  // namespace FS
