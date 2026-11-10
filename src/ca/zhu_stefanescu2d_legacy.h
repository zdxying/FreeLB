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

// zhu_stefanescu2d_legacy.h
//
// Legacy PopLattice-path Zhu-Stefanescu 2D cellular automaton
// (CA::ZhuStefanescu2D).  Only used by examples/legacy/cazs2d; kept out of
// ca/zhu_stefanescu2d.h so the active BlockLattice path does not pull in
// data_struct/lattice.h (PopLattice) -- the only data_struct -> geometry
// module edge.

#pragma once

#include "ca/zhu_stefanescu2d.h"   // CAType/CAFlag, latset, D2Q sets, ZSConverter
#include "data_struct/lattice.h"   // PopLattice, RhoLattice

namespace CA {

template <typename T, typename LatSet>
class ZhuStefanescu2D {
 private:

  Geometry2D<T>& Geo;
  ZSConverter<T>& ConvCA;
  RhoLattice<T>& lbmSO;
  RhoLattice<T>& lbmTH;

  int Ni;
  int Nj;
  int N;

  // preferred growth angle to x-axis, manually chosen
  T Theta;
  // nucleation sites
  int SiteId;
  // anisotropy coefficient, manually chosen
  T delta;
  // Gibbs-Thomson coefficient
  T GT;
  // initial composition C0
  T C0;
  // equilibrium liquidus temperature at the initial composition C0
  T Tl_eq;
  // slope of liquidus line
  T m_l;
  // partition coefficient
  T Part_Coef;
  // (1 - partition coefficient)
  T _Part_Coef;

  // interface cells
  std::vector<std::size_t> Interface;


  VectorFieldAOS<T, LatSet::d>& Velocity;

  // state field std::uint8_t
  ScalarField<CAType> State;
  // flag field std::uint8_t
  ScalarField<CAFlag> Flag;
  // solid fraction
  ScalarField<T> Fs;
  // delta solid fraction
  ScalarField<T> Delta_Fs;
  // curvature of Solid-Liquid interface
  ScalarField<T> Curvature;
  // Solid phase composition
  ScalarField<T> C_Solids;
  // excess rho
  VectorFieldAOS<T, LatSet::q> ExcessC;
  // colleted excess rho
  ScalarField<T> ExcessC_;
  // solid count
  std::size_t SolidCount;

  // nbr index
  std::array<int, LatSet::q> Delta_Index;

 public:
  ZhuStefanescu2D(ZSConverter<T>& convca, RhoLattice<T>& lbmso, RhoLattice<T>& lbmth,
                  PopLattice<T, LatSet>& lbmns, T delta_, T theta, int siteid,
                  int num = LatSet::q);
  ~ZhuStefanescu2D() {}

  // get field data
  std::vector<std::size_t>& getInterface() { return Interface; }
  ScalarField<CAType>& getState() { return State; }
  ScalarField<CAFlag>& getFlag() { return Flag; }
  ScalarField<T>& getFs() { return Fs; }
  ScalarField<T>& getDeltaFs() { return Delta_Fs; }
  ScalarField<T>& getCurvature() { return Curvature; }
  ScalarField<T>& getCSolids() { return C_Solids; }
  VectorFieldAOS<T, LatSet::q>& getExcessC() { return ExcessC; }
  ScalarField<T>& getExcessC_() { return ExcessC_; }
  // setup
  void Setup(int id, int num);
  void UpdateInterface();
  void UpdateCurvature(T limit = 4);
  // get C_eq, interface equilibrium composition
  T getC_eq(int id);
  // get anisotropy function g(theta)
  T getanisotropy(int id);
  // get interface normal to x-axis
  T getPhi(int id);

  void UpdateDeltaFs();
  // grow, omp parallel is enabled, no erase operation
  // if cells are large, use Grow_erase_s may improve performance
  void Grow();
  void DistributeExcessC(int id, T excessC);
  void CollectExcessC();
  void SimpleCapture();

  void TypeConversion();

  void PreCapture_s(std::vector<int>& precaps);
  void PreCapture_omp(std::vector<int>& precaps, int start, int end);
  // void SLICapture();
  // void slicapture(ZScell2D<T>& cells);
  // void SLICapture_s(const std::vector<int>& precaps,
  //                   std::vector<ZScell2D<T>>& cells);
  T getSolidCountFracton();
  // calc position of the vertices of the polygonal for capture
  // get Lphi
  // inline T getLphi(T phi);
  // get vertice of polygonal
  // return 3 state -1, 0, 1, compared with epsilon and 0
  inline int sign(T x);
  // get position of neighbors with State = -1
  // template <int start = 0>
  // inline void getNbrs_pos(int id, std::vector<int>& nbrs);

  // experimental
  inline T getAveNbrPopRho(int id);
  inline T getStatisticalPopRho();

  void apply_SimpleCapture() {
    UpdateInterface();
    UpdateCurvature();
    UpdateDeltaFs();
    Grow();
    CollectExcessC();
    SimpleCapture();
    TypeConversion();
  }
  // void apply_SLICapture() {
  //   UpdateCurvature();
  //   UpdateDeltaFs();
  //   Grow();
  //   CollectExcessC();
  //   TypeConversion();
  // }

  bool hasNeighborType(int id, std::uint8_t type) const;
  bool hasNeighborFlag(int id, std::uint8_t flag) const;
};


}  // namespace CA
