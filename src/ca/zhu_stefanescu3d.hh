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

// Zhu-Stefanescu (Z-S) model for dendrite growth, implementations

#include "ca/zhu_stefanescu3d.h"

namespace CA {

template <typename T, typename LatSet, typename RhoLatSet>
ZhuStefanescu3D<T, LatSet, RhoLatSet>::ZhuStefanescu3D(
    ZSConverter<T>& convca, RhoLattice<T>& lbmso, RhoLattice<T>& lbmth,
    BasicLattice<T, LatSet>& lbmns, T psi, T theta, T phi, T epsilon,
    int siteid, int num)
    : Geo(lbmns.getGeo()),
      ConvCA(convca),
      lbmSO(lbmso),
      lbmTH(lbmth),
      Psi(psi),
      Theta(theta),
      Phi(phi),
      Epsilon(epsilon),
      SiteId(siteid),
      GT(convca.Lattice_GT_Coef),
      C0(lbmso.getLatRhoInit()),
      Tl_eq(convca.get_LatTliq(lbmso.getLatRhoInit())),
      m_l(convca.Lattice_m_Liq),
      Part_Coef(convca.Part_Coef),
      _Part_Coef(T(1) - convca.Part_Coef),
      Velocity(lbmns.getVelocityField()),
      Ni(lbmns.getGeo().getNx()),
      Nj(lbmns.getGeo().getNy()),
      Nk(lbmns.getGeo().getNz()),
      N(lbmns.getGeo().getVoxelsNum()),
      State(lbmns.getGeo().getVoxelsNum(), CAType::Fluid),
      Flag(lbmns.getGeo().getVoxelsNum(), CAFlag::None),
      Fs(lbmns.getGeo().getVoxelsNum(), T(0)),
      Delta_Fs(lbmns.getGeo().getVoxelsNum(), T(0)),
      WMC(lbmns.getGeo().getVoxels().size(), T(0)),
      Q(lbmns.getGeo().getVoxels().size(), T(0)),
      NGradFs(lbmns.getGeo().getVoxels().size(), T(0)),
      C_Solids(lbmns.getGeo().getVoxelsNum(), T(0)),
      ExcessC_(lbmns.getGeo().getVoxelsNum(), T(0)),
      ExcessC(lbmns.getGeo().getVoxelsNum(), Vector<T, RhoLatSet::q>{}) {
  Delta_Index = make_Array<int, LatSet::q>(
        [&](int i) { return LatSet::c[i] * Geo.getProjection(); });
        
  Interface.reserve(Ni * Nj * 2 + Ni * Nk * 2 + Nj * Nk * 2);
  // (cos(Ѱ) cos(Φ) - sin(Ѱ) cos(θ) sin(Φ) | sin(Ѱ) cos(θ) cos(Φ)
  // + cos(Ѱ) sin(Φ) | sin(Ѱ) sin(θ)
  // -cos(Ѱ) cos(θ) sin(Φ) - sin(Ѱ) cos(Φ) |
  // cos(Ѱ) cos(θ) cos(Φ) - sin(Ѱ) sin(Φ) | cos(Ѱ) sin(θ)
  // sin(θ) sin(Φ) | sin(θ)(-cos(Φ)) | cos(θ))
  S[0] = Vector<T, 3>(cos(Psi) * cos(Phi) - sin(Psi) * cos(Theta) * sin(Phi),
                      sin(Psi) * cos(Theta) * cos(Phi) + cos(Psi) * sin(Phi),
                      sin(Psi) * sin(Theta));
  S[1] = Vector<T, 3>(-cos(Psi) * cos(Theta) * sin(Phi) - sin(Psi) * cos(Phi),
                      cos(Psi) * cos(Theta) * cos(Phi) - sin(Psi) * sin(Phi),
                      cos(Psi) * sin(Theta));
  S[2] =
      Vector<T, 3>(sin(Theta) * sin(Phi), -sin(Theta) * cos(Phi), cos(Theta));

  S_transpose[0] = Vector<T, 3>(S[0][0], S[1][0], S[2][0]);
  S_transpose[1] = Vector<T, 3>(S[0][1], S[1][1], S[2][1]);
  S_transpose[2] = Vector<T, 3>(S[0][2], S[1][2], S[2][2]);

  // setup
  Setup(SiteId, num);
}

template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::Setup(int id, int num) {
  // set state field from geometry flag
  Geo.forEachVoxel(Geo.getVoidflag(),
                   [&](int id) { State.SetField(id, CAType::Boundary); });
  std::cout << "[Zhu-Stefanescu CA 3D]"
            << "\n";
  // get preferred growth angle to x-axis
  std::cout << "Preferred growth angle: Pi(" << Psi / M_PI << ","
            << Theta / M_PI << "," << Phi / M_PI << ") (Z-X-Z)"
            << "\n"
            << "Part_Coef: " << Part_Coef << "\n";
  // get initial undercooling
  T Tl = lbmTH.getLatRhoInit();
  T deltaT = Tl_eq - Tl;
  std::cout << "Initial undercooling: " << deltaT << " | "
            << ConvCA.TempConv.getPhysDTemp(deltaT) << std::endl;
  const Vector<T, 3>& vox = Geo.getVoxel(id);
  // set to solid phase
  if (isType(id, CAType::Boundary)) {
    std::cout << "Warning: Setup at (" << vox[0] << ", " << vox[1] << ", "
              << vox[2] << "), CAType = Boundary" << std::endl;
  } else {
    State.SetField(id, CAType::Solid);
    Fs.SetField(id, T(1));
    Velocity.SetField(id, Vector<T, LatSet::d>{});
    lbmSO.getRhoField().SetField(id, lbmSO.getRho(id) * Part_Coef);
    std::cout << "Setup at (" << vox[0] << ", " << vox[1] << ", " << vox[2]
              << "), id = " << id << " succeeded" << std::endl;
  }
  // num neighbors
  int count = 0;
  for (int i = 1; i < num; i++) {
    int idn = Geo.template getNeighborId<LatSet>(id, i);
    if (isType(idn, CAType::Fluid)) {
      State.SetField(idn, CAType::Interface);
      Velocity.SetField(idn, Vector<T, LatSet::d>{});
      Interface.push_back(idn);
      count++;
    }
  }
  if (count != num - 1) {
    std::cout << "Warning: Incomplete neighbors setup: " << count << "/" << num
              << std::endl;
  }
}

template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::UpdateInterface() {
  Interface.clear();
  for (int id = 0; id < N; ++id) {
    if (isType(id, CAType::Interface)) Interface.push_back(id);
  }
  // clear excessC_
  ExcessC_.getField(0).Init(T(0));
}

template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::UpdateDeltaFs() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : Interface) {
    // deltaf is the increased solid fraction during one time step, given by:
    // deltaf = (C_eq - Cl)/(C_eq*(1-k))
    T C_eq = getC_eq(id);
    T deltaf = (C_eq - lbmSO.getRho(id)) / (C_eq * _Part_Coef);
    // getAveNbrPopRho(id)
    //  T deltaf = (C_eq - getAveNbrPopRho(id)) / (C_eq * _Part_Coef);
    if (deltaf > 1) {
      const Vector<T, 3>& vox = Geo.getVoxel(id);
      std::cout << "error: at (" << vox[0] << "," << vox[1] << "," << vox[2]
                << "), id = " << id << " ,deltaf = " << deltaf << ",\t"
                << "Ceq = " << C_eq << ",\t"
                << "Cl = " << lbmSO.getRho(id) << std::endl;
      exit(-1);
    }
    deltaf = deltaf < 0 ? 0 : deltaf;
    Delta_Fs.SetField(id, deltaf);
  }
}
template <typename T, typename LatSet, typename RhoLatSet>
inline T ZhuStefanescu3D<T, LatSet, RhoLatSet>::getC_eq(int id) {
  // NOTE that m_l > 0
  T Ceq = C0 + ((Tl_eq - lbmTH.getRho(id)) - GT * WMC.get(id)) / m_l;
  if (Ceq > 1 || Ceq < 0) {
    const Vector<T, 3>& vox = Geo.getVoxel(id);
    std::cout << "error: at (" << vox[0] << "," << vox[1] << "," << vox[2]
              << "), id = " << id << " ,Ceq = " << Ceq << ",\t"
              << "Tl_eq - Tl = " << Tl_eq - lbmTH.getRho(id) << ",\t"
              << "WMC = " << WMC.get(id) << ",\t"
              << "GT * WMC.get(id) = " << GT * WMC.get(id) << std::endl;
    exit(-1);
  }
  return Ceq;
}

template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::Grow() {
#pragma omp parallel for num_threads(Thread_Num)
  for (int id : Interface) {
    T delta_Fs = Delta_Fs.get(id);
    T Fs_temp = Fs.get(id) + delta_Fs;
    if (Fs_temp >= T(1) ||
        !hasNeighborType(
            id, static_cast<std::uint8_t>(CAType::Fluid | CAType::Interface))) {
      // get modified delta_Fs
      delta_Fs = T(1) - Fs.get(id);
      Fs.SetField(id, T(1));
      Delta_Fs.SetField(id, T(0));
      C_Solids.get(id) += Part_Coef * delta_Fs * lbmSO.getRho(id);
      // distribute excess solute to neighbors
      DistributeExcessC(id, _Part_Coef * delta_Fs * lbmSO.getRho(id));
      lbmSO.getRho(id) = C_Solids.get(id);
      Flag.SetField(id, CAFlag::toSolid);
    } else {
      Fs.SetField(id, Fs_temp);
      C_Solids.get(id) += Part_Coef * delta_Fs * lbmSO.getRho(id);
      ExcessC_.SetField(id, _Part_Coef * delta_Fs * lbmSO.getRho(id));
    }
  }
}

template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::DistributeExcessC(int id,
                                                              T excessC) {
  std::vector<int> nbrs;
  std::vector<int> dirs;
  nbrs.reserve(RhoLatSet::q);
  dirs.reserve(RhoLatSet::q);
  T sum = T(0);
  for (int i = 1; i < RhoLatSet::q; i++) {
    int idn = id + Delta_Index[i];
    if (isType(idn, CAType::Fluid)) {
      nbrs.push_back(idn);
      dirs.push_back(i);
      sum += RhoLatSet::w[i];
    }
  }
  T inv_sum = T(1) / sum;
  if (nbrs.size() == 0) return;
  for (int i = 0; i < nbrs.size(); i++) {
    ExcessC.get(nbrs[i])[dirs[i]] = excessC * RhoLatSet::w[dirs[i]] * inv_sum;
  }
}

template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::CollectExcessC() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    for (int i = 1; i < RhoLatSet::q; i++) {
      ExcessC_.get(id) += ExcessC.get(id)[i];
    }
  }
  // clear ExcessC
  ExcessC.getField(0).Init(Vector<T, RhoLatSet::q>{});
}

template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::SimpleCapture() {
#pragma omp parallel for num_threads(Thread_Num)
  for (int id = 0; id < N; ++id) {
    if (isType(id, CAType::Fluid) && hasNeighborFlag(id, CAFlag::toSolid)) {
      Flag.SetField(id, CAFlag::toInterface);
      Velocity.SetField(id, Vector<T, LatSet::d>{});
    }
  }
}

template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::TypeConversion() {
#pragma omp parallel for num_threads(Thread_Num)
  for (int id = 0; id < N; ++id) {
    if (isFlag(id, CAFlag::toInterface)) {
      State.SetField(id, CAType::Interface);
    } else if (isFlag(id, CAFlag::toSolid)) {
      State.SetField(id, CAType::Solid);
    }
    // clear flag
    Flag.SetField(id, CAFlag::None);
  }
}

// get wmc based on neighbor solid fraction gradient
template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::UpdateNGradFs() {
  FDM3D<T> FDM3D_gradFs(Geo, Fs.getField());
  T lim = std::numeric_limits<T>::epsilon() * 1e3;
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    // gradient of solid fraction
    Vector<T, 3> gradFs = FDM3D_gradFs.grad(id);
    T normgradFs = gradFs.getnorm();
    // if (normgradFs < std::numeric_limits<T>::epsilon())
    if (normgradFs < lim) {
      NGradFs.template SetField<0>(id, T(0));
      NGradFs.template SetField<1>(id, T(0));
      NGradFs.template SetField<2>(id, T(0));
      Q.SetField(id, T(0));
    } else {
      NGradFs.template SetField<0>(id, gradFs[0] / normgradFs);
      NGradFs.template SetField<1>(id, gradFs[1] / normgradFs);
      NGradFs.template SetField<2>(id, gradFs[2] / normgradFs);
      Q.SetField(id, pow(NGradFs.getField(0)[id], 4) +
                         pow(NGradFs.getField(1)[id], 4) +
                         pow(NGradFs.getField(2)[id], 4));
    }
  }
}

template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::UpdateWMC() {
  // (if not specified, column vector is used)
  // nabla = (partial_x, partial_y, partial_z)
  // Inabla = (partial_x, 0, 0; 0, partial_y, 0; 0, 0, partial_z)
  // gradf = nablaf = (partial_x f, partial_y f, partial_z f)
  // gradfnorm = |gradf| = sqrt(gradf * gradf)
  // Inablagradf = Inabla * gradf
  //             = (partial_xx f, partial_yy f, partial_zz f)
  // divgradf = partial_xx f + partial_yy f + partial_zz f
  // n = gradf / |gradf| = (nx, ny, nz)
  //   = (partial_x f, partial_y f, partial_z f)/ gradfnorm
  // Inablan = Inabla * n = (partial_x nx, partial_y ny, partial_z nz)
  //         = (partial_xx f, partial_yy f, partial_zz f)/ gradfnorm
  // grad' = S(transpose) * grad
  // H' = S(transpose) * H * S
  // where H is Hessian matrix:
  // H(f) = [partial_xx f, partial_xy f, partial_xz f]
  //        [partial_yx f, partial_yy f, partial_yz f]
  //        [partial_zx f, partial_zy f, partial_zz f]
  // wmc is weighted mean curvature, given by:
  // (use eps = \epsilon, Vec_n = grad(f)/|grad(f)| = (nx, ny, nz))
  // wmc = (3 eps - 1) (partial_x nx + partial_y ny + partial_z nz)
  //      - 48 eps (nx^2 partial_x nx + ny^2 partial_y ny + nz^2 partial_z nz)
  //      + 12 eps Q (partial_x nx + partial_y ny + partial_z nz)
  //      + 12 eps (nx partial_x Q + ny partial_y Q + nz partial_z Q)
  // where Q = nx^4 + ny^4 + nz^4

  // get field of normalized gradient f
  UpdateNGradFs();

  FDM3D<T> FDM3D_nx(Geo, NGradFs.getField(0));
  FDM3D<T> FDM3D_ny(Geo, NGradFs.getField(1));
  FDM3D<T> FDM3D_nz(Geo, NGradFs.getField(2));
  FDM3D<T> FDM3D_Q(Geo, Q.getField());

  // ny_x = nx_y, nz_x = nx_z, nz_y = ny_z;
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : Interface) {
    // get transformed n
    Vector<T, 3> n =
        Vector<T, 3>(NGradFs.getField(0)[id], NGradFs.getField(1)[id],
                     NGradFs.getField(2)[id]);
    Vector<T, 3> ntrans = Vector<T, 3>(S_transpose[0] * n, S_transpose[1] * n,
                                       S_transpose[2] * n);
    // get partial derivatives of n
    T nx_x = FDM3D_nx.p_x(id);
    T ny_y = FDM3D_ny.p_y(id);
    T nz_z = FDM3D_nz.p_z(id);
    T nx_y = FDM3D_nx.p_y(id);
    T nx_z = FDM3D_nx.p_z(id);
    T ny_z = FDM3D_ny.p_z(id);
    // Hessian matrix
    Vector<T, 3> Hcolumn1 = Vector<T, 3>(nx_x, nx_y, nx_z);
    Vector<T, 3> Hcolumn2 = Vector<T, 3>(nx_y, ny_y, ny_z);
    Vector<T, 3> Hcolumn3 = Vector<T, 3>(nx_z, ny_z, nz_z);

    T nx_xtrans = S_transpose[0] * Hcolumn1 * S[0][0] +
                  S_transpose[0] * Hcolumn2 * S[1][0] +
                  S_transpose[0] * Hcolumn3 * S[2][0];
    T ny_ytrans = S_transpose[1] * Hcolumn1 * S[0][1] +
                  S_transpose[1] * Hcolumn2 * S[1][1] +
                  S_transpose[1] * Hcolumn3 * S[2][1];
    T nz_ztrans = S_transpose[2] * Hcolumn1 * S[0][2] +
                  S_transpose[2] * Hcolumn2 * S[1][2] +
                  S_transpose[2] * Hcolumn3 * S[2][2];
    T nx_ytrans = S_transpose[0] * Hcolumn1 * S[0][1] +
                  S_transpose[0] * Hcolumn2 * S[1][1] +
                  S_transpose[0] * Hcolumn3 * S[2][1];
    T nx_ztrans = S_transpose[0] * Hcolumn1 * S[0][2] +
                  S_transpose[0] * Hcolumn2 * S[1][2] +
                  S_transpose[0] * Hcolumn3 * S[2][2];
    T ny_ztrans = S_transpose[1] * Hcolumn1 * S[0][2] +
                  S_transpose[1] * Hcolumn2 * S[1][2] +
                  S_transpose[1] * Hcolumn3 * S[2][2];

    T nx3 = 4 * pow(ntrans[0], 3);
    T ny3 = 4 * pow(ntrans[1], 3);
    T nz3 = 4 * pow(ntrans[2], 3);

    // T wmc_ =
    //     (3 * Epsilon - 1) * (nx_xtrans + ny_ytrans + nz_ztrans) -
    //     48 * Epsilon *
    //         (pow(ntrans[0], 2) * nx_xtrans + pow(ntrans[1], 2) * ny_ytrans +
    //          pow(ntrans[2], 2) * nz_ztrans) +
    //     12 * Epsilon *
    //         (pow(ntrans[0], 4) + pow(ntrans[1], 4) + pow(ntrans[2], 4)) *
    //         (nx_xtrans + ny_ytrans + nz_ztrans) +
    //     12 * Epsilon *
    //         (ntrans[0] * (nx3 * nx_xtrans + ny3 * nx_ytrans + nz3 *
    //         nx_ztrans) +
    //          ntrans[1] * (nx3 * nx_ytrans + ny3 * ny_ytrans + nz3 *
    //          ny_ztrans) + ntrans[2] * (nx3 * nx_ztrans + ny3 * ny_ztrans +
    //          nz3 * nz_ztrans));

    T wmc_ =
        (3 * Epsilon - 1) * (nx_xtrans + ny_ytrans + nz_ztrans) -
        48 * Epsilon *
            (pow(ntrans[0], 2) * nx_xtrans + pow(ntrans[1], 2) * ny_ytrans +
             pow(ntrans[2], 2) * nz_ztrans) +
        12 * Epsilon * Q.get(id) * (nx_xtrans + ny_ytrans + nz_ztrans) +
        12 * Epsilon *
            (ntrans[0] * FDM3D_Q.p_x(id) + ntrans[1] * FDM3D_Q.p_y(id) +
             ntrans[2] * FDM3D_Q.p_z(id));
    WMC.SetField(id, wmc_);
  }
}

template <typename T, typename LatSet, typename RhoLatSet>
void ZhuStefanescu3D<T, LatSet, RhoLatSet>::UpdateWMC_() {
  UpdateNGradFs();

  FDM3D<T> FDM3D_nx(Geo, NGradFs.getField(0));
  FDM3D<T> FDM3D_ny(Geo, NGradFs.getField(1));
  FDM3D<T> FDM3D_nz(Geo, NGradFs.getField(2));
  FDM3D<T> FDM3D_Q(Geo, Q.getField());
  // T lim = T(3);
  // ny_x = nx_y, nz_x = nx_z, nz_y = ny_z;
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : Interface) {
    T nx = NGradFs.getField(0)[id];
    T ny = NGradFs.getField(1)[id];
    T nz = NGradFs.getField(2)[id];
    // get partial derivatives of n
    T nx_x = FDM3D_nx.p_x(id);
    T ny_y = FDM3D_ny.p_y(id);
    T nz_z = FDM3D_nz.p_z(id);

    T wmc_ = (3 * Epsilon - 1) * (nx_x + ny_y + nz_z) -
             48 * Epsilon * (nx * nx * nx_x + ny * ny * ny_y + nz * nz * nz_z) +
             12 * Epsilon * Q.get(id) * (nx_x + ny_y + nz_z) +
             12 * Epsilon *
                 (nx * FDM3D_Q.p_x(id) + ny * FDM3D_Q.p_y(id) +
                  nz * FDM3D_Q.p_z(id));
    // wmc_ = wmc_ > lim ? lim : wmc_;
    // wmc_ = wmc_ < -lim ? -lim : wmc_;
    WMC.SetField(id, wmc_);
  }
}

// get wmc based on neighbor solid fraction gradient
// template <typename T, typename LatSet>
// void ZhuStefanescu3D<T, LatSet>::Calc_wmc_() {
//   // WMC = -div(grad Anisotropy)
//   CA.UpdatenGradFs(Epsilon);
//   FDM3D<T> FDM3D_anisotropy(Geo.getVoxels(), CA.getAnisotropy());
//   int id = 0;
//   // ny_x = nx_y, nz_x = nx_z, nz_y = ny_z;
// #pragma omp parallel for private(id) num_threads(Thread_Num)
//   for (GrowingZScell3D<T> &cell : Cells) {
//     id = cell.Id;
//     T wmc_ = -FDM3D_anisotropy.divgrad(id);
//     cell.wmc = wmc_;
//   }
// }
// template <typename T, typename LatSet>
// inline T ZhuStefanescu3D<T, LatSet>::getAveNbrPopRho(int id) {
//   // get average poprho of neighbors with State != -1
//   T sum = lbmSO.getPoprho(id);
//   int count = 1;
//   for (int i = 0; i < LatSet::q; i++) {
//     int idn = id + LatSet::Nbr[i];
//     if (CA.getState(idn) != -1) {
//       sum += lbmSO.getPoprho(idn);
//       count++;
//     }
//   }
//   return sum / count;
// }
// template <typename T, typename LatSet>
// inline T ZhuStefanescu3D<T, LatSet>::getStatisticalPopRho() {
//   // get average poprho of all cells
//   // for cells with State = 1 (interface cells)
//   // Rho = C_Solid * fs + C_Liquid * (1 - fs)
//   // for cells with State = 0 (liquid cells) or State = -1 (solid cells)
//   // Rho = PopRho
//   T sum = T(0);
//   int count = 0;
//   for (const Voxel<T, 3> &vox : Geo.getVoxels()) {
//     int id = vox.getId();
//     if (isState(id, CAState::Active)) {
//       if (Fs.getField()[id] > 0) {
//         sum += C_Solid.getField()[id] * Fs.getField()[id] +
//                lbmSO.getPoprho(id) * (1 - Fs.getField()[id]);
//         count++;
//       } else {
//         sum += lbmSO.getPoprho(id);
//       }
//     } else {
//       sum += lbmSO.getPoprho(id);
//     }
//   }
//   return sum / Geo.getVoxels().size();
// }
// // get wmc based on neighbor solid fraction gradient
// template <typename T, typename LatSet>
// void ZhuStefanescu3D<T, LatSet>::Calc_wmc_() {
//   // wmc is weighted mean curvature, given by:
//   // (use eps = \epsilon, Vec_n = grad(f)/|grad(f)| = (nx, ny, nz))
//   // wmc = (3 eps - 1) (partial_x nx + partial_y ny + partial_z nz)
//   //      - 48 eps (nx^2 partial_x nx + ny^2 partial_y ny + nz^2 partial_z
//   nz)
//   //      + 12 eps Q (partial_x nx + partial_y ny + partial_z nz)
//   //      + 12 eps (nx partial_x Q + ny partial_y Q + nz partial_z Q)
//   // where Q = nx^4 + ny^4 + nz^4
//   // (if not specified, column vector is used)
//   // nabla = (partial_x, partial_y, partial_z)
//   // Inabla = (partial_x, 0, 0; 0, partial_y, 0; 0, 0, partial_z)
//   // gradf = nablaf = (partial_x f, partial_y f, partial_z f)
//   // gradfnorm = |gradf| = sqrt(gradf * gradf)
//   // Inablagradf = Inabla * gradf
//   //             = (partial_xx f, partial_yy f, partial_zz f)
//   // divgradf = partial_xx f + partial_yy f + partial_zz f
//   // n = gradf / |gradf| = (nx, ny, nz)
//   //   = (partial_x f, partial_y f, partial_z f)/ gradfnorm
//   // Inablan = Inabla * n = (partial_x nx, partial_y ny, partial_z nz)
//   //         = (partial_xx f, partial_yy f, partial_zz f)/ gradfnorm
//   // grad' = S(transpose) * grad
//   // H' = S(transpose) * H * S
//   // where H is Hessian matrix:
//   // H(f) = [partial_xx f, partial_xy f, partial_xz f]
//   //        [partial_yx f, partial_yy f, partial_yz f]
//   //        [partial_zx f, partial_zy f, partial_zz f]
//   using FDM = FDM3D<T>;
//   int id = 0;
//   Vector<T, 3> gradf;
//   Vector<T, 3> Inablagradf;
//   Vector<T, 3> n;
//   Vector<T, 3> Inablan;
//   T gradfnorm;
//   T nx3, ny3, nz3;
//   T np_xy, np_xz, np_yz;
// #pragma omp parallel for private(id, wmc_, gradf, Inablagradf, n, Inablan,
//
//                                      gradfnorm, nx3, ny3, nz3, np_xy,
//                                      np_xz, \ np_yz)
//                                      num_threads(Thread_Num)
//   for (int i = 0; i < Cells.size(); i++) {
//     id = Cells[i].Id;
//     gradf = FDM::grad(id);
//     // Vector<T, 3> gradf_ = (S_transpose[0]*gradf, S_transpose[1]*gradf,
//     // S_transpose[2]*gradf);
//     Inablagradf = FDM::Inablagrad(id);
//     // Vector<T, 3> Inablagradf_ =
//     gradfnorm = gradf.getnorm();
//     n = gradf.getnormalize();
//     Inablan = Inablagradf / gradfnorm;
//     np_xy = FDM::p_xy(id) / gradfnorm;
//     np_xz = FDM::p_xz(id) / gradfnorm;
//     np_yz = FDM::p_yz(id) / gradfnorm;
//     nx3 = 4 * pow(n[0], 3);
//     ny3 = 4 * pow(n[1], 3);
//     nz3 = 4 * pow(n[2], 3);
//     Cells[i].wmc = (3 * Epsilon - 1) * (Inablan[0] + Inablan[1] +
//     Inablan[2])
//     -
//                    48 * Epsilon *
//                        (n[0] * n[0] * Inablan[0] + n[1] * n[1] * Inablan[1]
//                        +
//                         n[2] * n[2] * Inablan[2]) +
//                    12 * Epsilon * (pow(n[0], 4) + pow(n[1], 4) + pow(n[2],
//                    4)) *
//                        (Inablan[0] + Inablan[1] + Inablan[2]) +
//                    12 * Epsilon *
//                        (n[0] * (nx3 * Inablan[0] + ny3 * np_xy + nz3 *
//                        np_xz)
//                        +
//                         n[1] * (nx3 * np_xy + ny3 * Inablan[1] + nz3 *
//                         np_yz)
//                         + n[2] * (nx3 * np_xz + ny3 * np_yz + nz3 *
//                         Inablan[2]));
//   }
// }

}  // namespace CA