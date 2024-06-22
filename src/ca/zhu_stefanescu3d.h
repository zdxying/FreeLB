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

// Zhu-Stefanescu (Z-S) model for dendrite growth
// A three-dimensional sharp interface model for the
// quantitative simulation of solutal dendritic growth 2010
#pragma once

#include <atomic>

#include "data_struct/field_struct.h"
#include "ca/cazs.h"
#include "lbm/lattice_set.h"
#include "lbm/unit_converter.h"
#include "utils/fdm_solver.h"
#include "utils/util.h"

//----[driving force]----
// the driving force for dendritic growth is considered to be controlled by
// the difference between local interface equilibrium composition and
// local actual liquid composition
// deltaf is the increased solid fraction during one time step, given by:
// deltaf = (C_eq - Cl)/(C_eq*(1-k)), where:
// C_eq is the interface equilibrium composition, given by:
// C_eq = C0 + [(T_interface - Tl_eq) + GT* wmc ] / m
// C0 is the initial composition,
// Cl is the actual liquid composition
// k is the partition coefficient

// ----[weighted mean curvature]----
// wmc is weighted mean curvature, given by:
// (use eps = \epsilon, Vec_n = grad(f)/|grad(f)| = (nx, ny, nz))
// wmc = (3 eps - 1) (partial_x nx + partial_y ny + partial_z nz)
//      - 48 eps (nx^2 partial_x nx + ny^2 partial_y ny + nz^2 partial_z nz)
//      + 12 eps Q (partial_x nx + partial_y ny + partial_z nz)
//      + 12 eps (nx partial_x Q + ny partial_y Q + nz partial_z Q)
// where Q = nx^4 + ny^4 + nz^4

// However, the coordinate system use above is not the original one
// ----[coordinate transformation]----
// the original coordinate system is: [x,y,z]
// the transformed coordinate system is: [x',y',z']
// with axes parallel to [100] direction of the dendrite growth
// Three Euler angles are used to describe the transformation:
// theta, phi, psi
// x -> x' : psi + phi
// z -> z' : theta
// here we have orthogonal transformation matrix S[3][3]:
// cos(psi)cos(phi)-cos(theta)sin(phi)sin(psi)  |
// cos(psi)sin(phi)+cos(theta)cos(phi)sin(psi)  |  sin(psi)sin(theta)
// -sin(psi)cos(phi)-cos(theta)sin(phi)cos(psi) |
// -sin(psi)sin(phi)+cos(theta)cos(phi)cos(psi) |  cos(psi)sin(theta)
// sin(theta)sin(phi)                           |  -sin(theta)cos(phi) |
// cos(theta)
// i.e.: (cos(Ѱ) cos(Φ) - sin(Ѱ) cos(θ) sin(Φ) | sin(Ѱ) cos(θ) cos(Φ)
// + cos(Ѱ) sin(Φ) | sin(Ѱ) sin(θ)
// -cos(Ѱ) cos(θ) sin(Φ) - sin(Ѱ) cos(Φ) |
// cos(Ѱ) cos(θ) cos(Φ) - sin(Ѱ) sin(Φ) | cos(Ѱ) sin(θ)
// sin(θ) sin(Φ) | sin(θ)(-cos(Φ)) | cos(θ))
// which come from: (z(Ѱ)-x(θ)-z(Φ))
// (cos(Ѱ) | sin(Ѱ) | 0
// -sin(Ѱ) | cos(Ѱ) | 0
// 0 | 0 | 1).
// (1 | 0 | 0
// 0 | cos(θ) | sin(θ)
// 0 | -sin(θ) | cos(θ)).
// (cos(Φ) | sin(Φ) | 0
// -sin(Φ) | cos(Φ) | 0
// 0 | 0 | 1)

// Thus:
// grad' = S(transpose) * grad
// H' = S(transpose) * H * S
// where H is Hessian matrix:
// H(f) = [partial_xx f, partial_xy f, partial_xz f]
//        [partial_yx f, partial_yy f, partial_yz f]
//        [partial_zx f, partial_zy f, partial_zz f]

// rejected solute in an interface cell at each timestep
// deltaC = C_l(1-k)*deltaf

// nucleation: use single nucleation
// capture: too complex to extend 2D virtual interface tracking scheme
// use cellular automaton capture(simple capture):
// When an interface cell has been fully solidified, its neighboring liquid
// cells are captured as interface cells
namespace CA {

// multi-nuclei is NOT supported inside this class, to use multi-nuclei,
// it is recommended to instantiate multiple objects of this class
template <typename T, typename LatSet, typename RhoLatSet = LatSet>
class ZhuStefanescu3D {
 private:
  int Ni;
  int Nj;
  int Nk;
  int N;
  // preferred growth angle around z-axis 1st rotation
  T Psi;
  // preferred growth angle around x-axis 2nd rotation
  T Theta;
  // preferred growth angle around z-axis 3rd rotation
  T Phi;
  // orthogonal transformation matrix
  std::array<Vector<T, 3>, 3> S;
  std::array<Vector<T, 3>, 3> S_transpose;
  // degree of anisotropy of the surface energy
  T Epsilon;
  // nucleation sites
  int SiteId;
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
  // solid count
  std::size_t SolidCount;

  // interface cells
  std::vector<std::size_t> Interface;

  Geometry3D<T>& Geo;
  ZSConverter<T>& ConvCA;
  RhoLattice<T>& lbmSO;
  RhoLattice<T>& lbmTH;
  VectorFieldAOS<T, LatSet::d>& Velocity;

  // state field std::uint8_t
  ScalarField<CAType> State;
  // flag field std::uint8_t
  ScalarField<CAFlag> Flag;
  // solid fraction
  ScalarField<T> Fs;
  // delta solid fraction
  ScalarField<T> Delta_Fs;
  // Solid phase composition
  ScalarField<T> C_Solids;
  // excess rho
  VectorFieldAOS<T, RhoLatSet::q> ExcessC;
  // colleted excess rho
  ScalarField<T> ExcessC_;

  // weighted mean curvature
  ScalarField<T> WMC;
  // Q
  ScalarField<T> Q;
  // normalized gradient of solid fraction: nx, ny, nz
  VectorFieldSoA<T, 3> NGradFs;

   // nbr index
  std::array<int, LatSet::q> Delta_Index;

 public:
  ZhuStefanescu3D(ZSConverter<T>& convca, RhoLattice<T>& lbmso,
                  RhoLattice<T>& lbmth, PopLattice<T, LatSet>& lbmns, T psi,
                  T theta, T phi, T epsilon, int siteid, int num = LatSet::q);
  ~ZhuStefanescu3D() {}

  // get field data
  std::vector<std::size_t>& getInterface() { return Interface; }
  ScalarField<CAType>& getState() { return State; }
  ScalarField<CAFlag>& getFlag() { return Flag; }
  ScalarField<T>& getFs() { return Fs; }
  ScalarField<T>& getDeltaFs() { return Delta_Fs; }
  ScalarField<T>& getCSolids() { return C_Solids; }
  VectorFieldAOS<T, RhoLatSet::q>& getExcessC() { return ExcessC; }
  ScalarField<T>& getExcessC_() { return ExcessC_; }
  ScalarField<T>& getWMC() { return WMC; }
  ScalarField<T>& getQ() { return Q; }
  VectorFieldSoA<T, 3>& getNGradFs() { return NGradFs; }

  void Setup(int id, int num);
  void UpdateInterface();

  inline T getC_eq(int id);
  void UpdateDeltaFs();

  void Grow();

  void DistributeExcessC(int id, T excessC);
  void CollectExcessC();

  void SimpleCapture();
  void TypeConversion();

  // calculate curvature of Solid-Liquid interface in the whole domain
  void UpdateNGradFs();
  void UpdateWMC();
  void UpdateWMC_();

  T getSolidFraction() const { return T(SolidCount) / T(Ni-2)/ T(Nj-2)/ T(Nk-2); }

  // experimental
  inline T getAveNbrPopRho(int id);
  inline T getStatisticalPopRho();

  void apply_SimpleCapture() {
    UpdateInterface();
    UpdateWMC_();
    UpdateDeltaFs();
    Grow();
    CollectExcessC();
    SimpleCapture();
    TypeConversion();
  }
  T getSolidCountFracton() {
    SolidCount = 0;
    for (int id = 0; id < N; ++id) {
      if (isType(id, CAType::Solid)) ++SolidCount;
    }
    return T(SolidCount) / T(N);
  }
  inline bool isType(int id, std::uint8_t type) {
    return static_cast<bool>(State.get(id) & type);
  }
  inline bool isFlag(int id, std::uint8_t flag) {
    return static_cast<bool>(Flag.get(id) & flag);
  }
  void setType(int id, std::uint8_t type) { State.SetField(id, type); }
  void setFlag(int id, std::uint8_t flag) { Flag.SetField(id, flag); }
  bool hasNeighborType(int id, std::uint8_t type) {
    for (int i = 1; i < LatSet::q; ++i) {
      if (isType(id + Delta_Index[i], type)) return true;
    }
    return false;
  }
  bool hasNeighborFlag(int id, std::uint8_t flag) {
    for (int i = 1; i < LatSet::q; ++i) {
      if (isFlag(id + Delta_Index[i], flag)) return true;
    }
    return false;
  }
};

//   inline bool has_SolidNbr_simd(int id) {
//     bool inactive = false;
// #pragma omp simd reduction(|| : inactive)
//     for (int i = 1; i < LatSet::q; i++) {
//       if (CA.State[Id + LatSet::Nbr[i]] == -1) inactive = true;
//     }
//     return inactive;
//   }

}  // namespace CA