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

#pragma once

#include "ca/zhu_stefanescu2d.h"
#include "zhu_stefanescu2d.h"

namespace CA {
template <typename T, typename LatSet>
ZhuStefanescu2D<T, LatSet>::ZhuStefanescu2D(ZSConverter<T> &convca, RhoLattice<T> &lbmso,
                                            RhoLattice<T> &lbmth, BasicLattice<T, LatSet> &lbmns,
                                            T delta_, T theta, int siteid, int num)
    : Geo(lbmns.getGeo()), ConvCA(convca), lbmSO(lbmso), lbmTH(lbmth), delta(delta_), Theta(theta),
      SiteId(siteid), GT(convca.Lattice_GT_Coef), C0(lbmso.getLatRhoInit()),
      Tl_eq(convca.get_LatTliq(lbmso.getLatRhoInit())), m_l(convca.Lattice_m_Liq),
      Part_Coef(convca.Part_Coef), _Part_Coef(T(1) - convca.Part_Coef),
      Velocity(lbmns.getVelocityField()), Ni(lbmns.getGeo().getNx()), Nj(lbmns.getGeo().getNy()),
      N(lbmns.getGeo().getVoxelsNum()), State(lbmns.getGeo().getVoxelsNum(), CAType::Fluid),
      Flag(lbmns.getGeo().getVoxelsNum(), CAFlag::None), Fs(lbmns.getGeo().getVoxelsNum(), T(0)),
      Delta_Fs(lbmns.getGeo().getVoxelsNum(), T(0)), Curvature(lbmns.getGeo().getVoxelsNum(), T(0)),
      C_Solids(lbmns.getGeo().getVoxelsNum(), T(0)),
      ExcessC(lbmns.getGeo().getVoxelsNum(), Vector<T, LatSet::q>{}),
      ExcessC_(lbmns.getGeo().getVoxelsNum(), T(0)) {
  Delta_Index =
    make_Array<int, LatSet::q>([&](int i) { return LatSet::c[i] * Geo.getProjection(); });
  Interface.reserve(2 * (Ni + Nj));
  Setup(SiteId, num);
}

template <typename T, typename LatSet>
void ZhuStefanescu2D<T, LatSet>::Setup(int id, int num) {
  // set state field from geometry flag
  Geo.forEachVoxel(Geo.getVoidflag(), [&](int id) { State.SetField(id, CAType::Boundary); });

  std::cout << "[Zhu-Stefanescu 2D CA]" << std::endl;
  // get preferred growth angle to x-axis
  std::cout << "preferred growth angle: " << Theta / M_PI << " Pi" << std::endl;
  // get initial undercooling
  T Tl = lbmTH.getLatRhoInit();
  T deltaT = Tl_eq - Tl;
  std::cout << "Initial undercooling: " << deltaT << " | " << ConvCA.TempConv.getPhysDTemp(deltaT)
            << std::endl;
  const Vector<T, 2> &vox = Geo.getVoxel(id);
  // set to solid phase
  if (util::isFlag(State.get(id), CAType::Boundary)) {
    std::cout << "Warning: Setup at (" << vox[0] << ", " << vox[1] << "), CAType = Boundary"
              << std::endl;
  } else {
    State.SetField(id, CAType::Solid);
    Fs.SetField(id, T(1));
    Velocity.SetField(id, Vector<T, LatSet::d>{});
    lbmSO.getRhoField().SetField(id, lbmSO.getRho(id) * Part_Coef);
    std::cout << "Setup at (" << vox[0] << ", " << vox[1] << "), id = " << id << " succeeded"
              << std::endl;
  }

  // num neighbors
  int count = 0;
  for (int i = 1; i < num; i++) {
    int idn = Geo.template getNeighborId<LatSet>(id, i);
    if (util::isFlag(State.get(idn), CAType::Fluid)) {
      State.SetField(idn, CAType::Interface);
      Velocity.SetField(idn, Vector<T, LatSet::d>{});
      Interface.push_back(idn);
      count++;
    }
  }
  if (count != num - 1) {
    std::cout << "Warning: Incomplete neighbors setup. " << count << "/" << num << std::endl;
  }
}

template <typename T, typename LatSet>
void ZhuStefanescu2D<T, LatSet>::UpdateInterface() {
  Interface.clear();
  for (int id = 0; id < N; ++id) {
    if (util::isFlag(State.get(id), CAType::Interface)) Interface.push_back(id);
  }
  // clear excessC_
  ExcessC_.getField(0).Init(T(0));
}

template <typename T, typename LatSet>
void ZhuStefanescu2D<T, LatSet>::UpdateCurvature(T limit) {
  FDM2D<T> FDM(Ni, Nj, Fs.getField(0));
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : Interface) {
    T px = FDM.p_x(id);
    T py = FDM.p_y(id);
    T curv = px * px + py * py;
    curv = curv > 1e-3 ? pow(curv, -T(1.5)) : 0;
    T K_ = curv * (2 * px * py * FDM.p_xy(id) - px * px * FDM.p_yy(id) - py * py * FDM.p_xx(id));
    if (K_ > T(limit)) K_ = T(limit);
    if (K_ < T(-limit)) K_ = T(-limit);
    // set K
    Curvature.SetField(id, K_);
  }
}

template <typename T, typename LatSet>
void ZhuStefanescu2D<T, LatSet>::UpdateDeltaFs() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : Interface) {
    T C_eq = getC_eq(id);
    T deltaf = (C_eq - lbmSO.getRho(id)) / (C_eq * _Part_Coef);
    //
    if (deltaf > 1) {
      const Vector<T, 2> &vox = Geo.getVoxel(id);
      std::cout << "error: at (" << vox[0] << ", " << vox[1] << "), id = " << id
                << " ,deltaf = " << deltaf << ",\t"
                << "Ceq = " << C_eq << ",\t"
                << "Cl = " << lbmSO.getRho(id) << std::endl;
      exit(-1);
    }
    deltaf = deltaf < 0 ? 0 : deltaf;
    Delta_Fs.SetField(id, deltaf);
  }
}
template <typename T, typename LatSet>
T ZhuStefanescu2D<T, LatSet>::getC_eq(int id) {
  T Ceq = C0 + ((Tl_eq - lbmTH.getRho(id)) - GT * Curvature.get(id) * getanisotropy(id)) / m_l;
  if (Ceq > 1 || Ceq < 0) {
    const Vector<T, 2> &vox = Geo.getVoxel(id);
    std::cout << "error: at (" << vox[0] << ", " << vox[1] << "), id = " << id << " ,Ceq = " << Ceq
              << ",\t"
              << "Tl_eq - Tl = " << Tl_eq - lbmTH.getRho(id) << ",\t"
              << "K = " << Curvature.get(id) << ",\t"
              << "anisotropy = " << getanisotropy(id) << ",\t"
              << "GT * Curvature.get(id) * getanisotropy(id) = "
              << GT * Curvature.get(id) * getanisotropy(id) << std::endl;
    exit(-1);
  }
  return Ceq;
}
template <typename T, typename LatSet>
T ZhuStefanescu2D<T, LatSet>::getanisotropy(int id) {
  // g(phi, theta) = 1 - delta * cos(4*(phi - theta))
  // where: delta is the anisotropy coefficient, manually chosen
  // phi is growth angle to x-axis, which can be calculated by solid fraction:
  // theta is the preferred growth angle to x-axis
  return 1 - delta * cos(4 * (getPhi(id) - Theta));
}
template <typename T, typename LatSet>
T ZhuStefanescu2D<T, LatSet>::getPhi(int id) {
  FDM2D<T> FDM(Ni, Nj, Fs.getField(0));
  return atan2(FDM.p_y(id), FDM.p_x(id));
}

template <typename T, typename LatSet>
void ZhuStefanescu2D<T, LatSet>::Grow() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : Interface) {
    T delta_Fs = Delta_Fs.get(id);
    T Fs_temp = Fs.get(id) + delta_Fs;
    if (Fs_temp >= T(1) ||
        !hasNeighborType(id, static_cast<std::uint8_t>(CAType::Fluid | CAType::Interface))) {
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

template <typename T, typename LatSet>
void ZhuStefanescu2D<T, LatSet>::DistributeExcessC(int id, T excessC) {
  std::vector<int> nbrs;
  std::vector<int> dirs;
  nbrs.reserve(LatSet::q);
  dirs.reserve(LatSet::q);
  T sum = T(0);
  for (int i = 1; i < LatSet::q; i++) {
    int idn = Geo.template getNeighborId<LatSet>(id, i);
    if (util::isFlag(State.get(idn), CAType::Fluid)) {
      nbrs.push_back(idn);
      dirs.push_back(i);
      sum += LatSet::w[i];
    }
  }
  T inv_sum = T(1) / sum;
  if (nbrs.size() == 0) return;
  for (int i = 0; i < nbrs.size(); i++) {
    ExcessC.get(nbrs[i])[dirs[i]] = excessC * LatSet::w[dirs[i]] * inv_sum;
  }
}

template <typename T, typename LatSet>
void ZhuStefanescu2D<T, LatSet>::CollectExcessC() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    for (int i = 1; i < LatSet::q; i++) {
      ExcessC_.get(id) += ExcessC.get(id)[i];
    }
  }
  // clear ExcessC
  ExcessC.getField(0).Init(Vector<T, LatSet::q>{});
}

template <typename T, typename LatSet>
void ZhuStefanescu2D<T, LatSet>::SimpleCapture() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    if (util::isFlag(State.get(id), CAType::Fluid) && hasNeighborFlag(id, CAFlag::toSolid)) {
      Flag.SetField(id, CAFlag::toInterface);
      Velocity.SetField(id, Vector<T, LatSet::d>{});
    }
  }
}

template <typename T, typename LatSet>
void ZhuStefanescu2D<T, LatSet>::TypeConversion() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id = 0; id < N; ++id) {
    if (util::isFlag(Flag.get(id), CAFlag::toInterface)) {
      State.SetField(id, CAType::Interface);
    } else if (util::isFlag(Flag.get(id), CAFlag::toSolid)) {
      State.SetField(id, CAType::Solid);
    }
    // clear flag
    Flag.SetField(id, CAFlag::None);
  }
}

template <typename T, typename LatSet>
T ZhuStefanescu2D<T, LatSet>::getSolidCountFracton() {
  SolidCount = 0;
  for (int id = 0; id < N; ++id) {
    if (util::isFlag(State.get(id), CAType::Solid)) ++SolidCount;
  }
  return T(SolidCount) / T(Ni - 2) / T(Nj - 2);
}

template <typename T, typename LatSet>
bool ZhuStefanescu2D<T, LatSet>::hasNeighborType(int id, std::uint8_t type) const {
  for (int i = 1; i < LatSet::q; ++i) {
    if (util::isFlag(State.get(id + Delta_Index[i]), type)) return true;
  }
  return false;
}

template <typename T, typename LatSet>
bool ZhuStefanescu2D<T, LatSet>::hasNeighborFlag(int id, std::uint8_t flag) const {
  for (int i = 1; i < LatSet::q; ++i) {
    if (util::isFlag(Flag.get(id + Delta_Index[i]), flag)) return true;
  }
  return false;
}

// ------------------------------------------------------------------
// ------------------Block Zhu-Stefanescu 2D-------------------------
// ------------------------------------------------------------------

template <typename T, typename LatSet>
BlockZhuStefanescu2D<T, LatSet>::BlockZhuStefanescu2D(VectorFieldAOS<T, 2> &velo,
                                                      ZSConverter<T> &convca, RhoLattice<T> &lbmso,
                                                      RhoLattice<T> &lbmth, Block2D<T> &geo,
                                                      T delta_, T theta, std::size_t siteid,
                                                      bool setup, int num)
    : Geo(geo), ConvCA(convca), lbmSO(lbmso), lbmTH(lbmth), delta(delta_), Theta(theta),
      SiteId(siteid), GT(convca.Lattice_GT_Coef * pow(2, int(geo.getLevel()))),
      C0(lbmso.getLatRhoInit()), Tl_eq(convca.get_LatTliq(lbmso.getLatRhoInit())),
      m_l(convca.Lattice_m_Liq), Part_Coef(convca.Part_Coef), _Part_Coef(T(1) - convca.Part_Coef),
      Velocity(velo), SolidCount(0), Ni(geo.getNx()), Nj(geo.getNy()), N(geo.getN()),
      State(geo.getN(), CAType::Boundary), Flag(geo.getN(), CAFlag::None), Fs(geo.getN(), T(0)),
      Delta_Fs(geo.getN(), T(0)), Curvature(geo.getN(), T(0)), C_Solids(geo.getN(), T(0)),
      ExcessC(geo.getN(), T(0)), ExcessC_(geo.getN(), T(0)) {
  Delta_Index =
    make_Array<int, LatSet::q>([&](int i) { return LatSet::c[i] * Geo.getProjection(); });
  Interface.reserve(2 * (Ni + Nj));
  // set state field from geometry flag
  Geo.forEach(Geo.getGeoFlagField().getField(), Geo.getAABBflag(),
              [&](std::size_t id) { State.SetField(id, CAType::Fluid); });
  if (setup) {
    Setup(SiteId, num);
  }
}

template <typename T, typename LatSet>
void BlockZhuStefanescu2D<T, LatSet>::Setup(std::size_t id, int num) {
  std::cout << "[Zhu-Stefanescu 2D CA]" << std::endl;
  // get preferred growth angle to x-axis
  std::cout << "preferred growth angle: " << Theta / M_PI << " Pi" << std::endl;
  // get initial undercooling
  T Tl = lbmTH.getLatRhoInit();
  T deltaT = Tl_eq - Tl;
  std::cout << "Initial undercooling: " << deltaT << " | " << ConvCA.TempConv.getPhysDTemp(deltaT)
            << std::endl;
  Vector<int, 2> vox = Geo.getLoc(id);
  // set to solid phase
  if (util::isFlag(State.get(id), CAType::Boundary)) {
    std::cout << "Warning: Setup at (" << vox[0] << ", " << vox[1] << "), CAType = Boundary"
              << std::endl;
  } else {
    State.SetField(id, CAType::Solid);
    Fs.SetField(id, T(1));
    Velocity.SetField(id, Vector<T, 2>{});
    lbmSO.getRhoField().SetField(id, lbmSO.getRho(id) * Part_Coef);
    std::cout << "Setup at (" << vox[0] << ", " << vox[1] << "), id = " << id << " succeeded"
              << std::endl;
  }

  // num neighbors
  int count = 0;
  for (int i = 0; i < num; i++) {
    std::size_t idn = id + Delta_Index[i];
    Vector<T, 2> loc_t = Geo.getLoc_t(idn);
    if (Geo.isInside(loc_t)) {
      if (util::isFlag(State.get(idn), CAType::Fluid)) {
        State.SetField(idn, CAType::Interface);
        Velocity.SetField(idn, Vector<T, 2>{});
        Interface.push_back(idn);
        count++;
      }
    }
  }
  if (count != num) {
    std::cout << "Warning: Incomplete neighbors setup. " << count << "/" << num << std::endl;
  }
}

template <typename T, typename LatSet>
void BlockZhuStefanescu2D<T, LatSet>::UpdateInterface() {
  Interface.clear();
  // only update inner cells
  for (int j = Geo.getOverlap(); j < Nj - Geo.getOverlap(); ++j) {
    for (int i = Geo.getOverlap(); i < Ni - Geo.getOverlap(); ++i) {
      std::size_t id = i + j * Ni;
      if (util::isFlag(State.get(id), CAType::Interface)) Interface.push_back(id);
    }
  }
  // clear excessC_
  ExcessC_.getField(0).Init(T(0));
}

template <typename T, typename LatSet>
void BlockZhuStefanescu2D<T, LatSet>::UpdateCurvature(T limit) {
  FDM2D<T> FDM(Ni, Nj, Fs.getField(0));
  for (std::size_t id : Interface) {
    T px = FDM.p_x(id);
    T py = FDM.p_y(id);
    T curv = px * px + py * py;
    curv = curv > 1e-3 ? pow(curv, -T(1.5)) : 0;
    T K_ = curv * (2 * px * py * FDM.p_xy(id) - px * px * FDM.p_yy(id) - py * py * FDM.p_xx(id));
    if (K_ > T(limit)) K_ = T(limit);
    if (K_ < T(-limit)) K_ = T(-limit);
    // set K
    Curvature.SetField(id, K_);
  }
}

template <typename T, typename LatSet>
void BlockZhuStefanescu2D<T, LatSet>::UpdateDeltaFs() {
  for (std::size_t id : Interface) {
    T C_eq = getC_eq(id);
    T deltaf = (C_eq - lbmSO.getRho(id)) / (C_eq * _Part_Coef);
    //
    if (deltaf > 1) {
      Vector<int, 2> vox = Geo.getLoc(id);
      std::cout << "error: at (" << vox[0] << ", " << vox[1] << "), id = " << id
                << " ,deltaf = " << deltaf << ",\t"
                << "Ceq = " << C_eq << ",\t"
                << "Cl = " << lbmSO.getRho(id) << std::endl;
      exit(-1);
    }
    deltaf = deltaf < 0 ? 0 : deltaf;
    Delta_Fs.SetField(id, deltaf);
  }
}

template <typename T, typename LatSet>
T BlockZhuStefanescu2D<T, LatSet>::getC_eq(std::size_t id) {
  T Ceq = C0 + ((Tl_eq - lbmTH.getRho(id)) - GT * Curvature.get(id) * getanisotropy(id)) / m_l;
  if (Ceq > 1 || Ceq < 0) {
    Vector<int, 2> vox = Geo.getLoc(id);
    std::cout << "error: at (" << vox[0] << ", " << vox[1] << "), id = " << id << " ,Ceq = " << Ceq
              << ",\t"
              << "Tl_eq - Tl = " << Tl_eq - lbmTH.getRho(id) << ",\t"
              << "K = " << Curvature.get(id) << ",\t"
              << "anisotropy = " << getanisotropy(id) << ",\t"
              << "GT * Curvature.get(id) * getanisotropy(id) = "
              << GT * Curvature.get(id) * getanisotropy(id) << std::endl;
    exit(-1);
  }
  return Ceq;
}

template <typename T, typename LatSet>
T BlockZhuStefanescu2D<T, LatSet>::getanisotropy(std::size_t id) {
  // g(phi, theta) = 1 - delta * cos(4*(phi - theta))
  // where: delta is the anisotropy coefficient, manually chosen
  // phi is growth angle to x-axis, which can be calculated by solid fraction:
  // theta is the preferred growth angle to x-axis
  return 1 - delta * cos(4 * (getPhi(id) - Theta));
}
template <typename T, typename LatSet>
T BlockZhuStefanescu2D<T, LatSet>::getPhi(std::size_t id) {
  // TODO: efficiency may be improved?
  FDM2D<T> FDM(Ni, Nj, Fs.getField(0));
  return atan2(FDM.p_y(id), FDM.p_x(id));
}

template <typename T, typename LatSet>
void BlockZhuStefanescu2D<T, LatSet>::Grow() {
  for (std::size_t id : Interface) {
    T delta_Fs = Delta_Fs.get(id);
    T Fs_temp = Fs.get(id) + delta_Fs;
    if (Fs_temp >= T(1) ||
        !hasNeighborType(id, static_cast<std::uint8_t>(CAType::Fluid | CAType::Interface))) {
      // get modified delta_Fs
      delta_Fs = T(1) - Fs.get(id);
      Fs.SetField(id, T(1));
      Delta_Fs.SetField(id, T(0));
      C_Solids.get(id) += Part_Coef * delta_Fs * lbmSO.getRho(id);
      // distribute excess solute to neighbors
      DistributeExcessC(id, _Part_Coef * delta_Fs * lbmSO.getRho(id));
      lbmSO.getRho(id) = C_Solids.get(id);
      // Flag.SetField(id, CAFlag::toSolid);
      State.SetField(id, CAType::Solid);
      ++SolidCount;
    } else {
      Fs.SetField(id, Fs_temp);
      C_Solids.get(id) += Part_Coef * delta_Fs * lbmSO.getRho(id);
      ExcessC_.SetField(id, _Part_Coef * delta_Fs * lbmSO.getRho(id));
    }
  }
}

template <typename T, typename LatSet>
void BlockZhuStefanescu2D<T, LatSet>::DistributeExcessC(std::size_t id, T excessC) {
  // get array
  std::array<T *, LatSet::q> TempExcessC = ExcessC.template getArray<T>(id);
  std::vector<int> dirs;
  dirs.reserve(LatSet::q);
  T sum = T(0);
  for (int i = 0; i < LatSet::q; ++i) {
    std::size_t idn = id + Delta_Index[i];
    if (util::isFlag(State.get(idn), CAType::Fluid)) {
      dirs.push_back(i);
      sum += LatSet::w[i];
    }
  }
  T inv_sum = T(1) / sum;
  if (dirs.size() == 0) return;
  for (int i = 0; i < dirs.size(); ++i) {
    *(TempExcessC[dirs[i]]) = excessC * LatSet::w[dirs[i]] * inv_sum;
  }
}

template <typename T, typename LatSet>
void BlockZhuStefanescu2D<T, LatSet>::CollectExcessC() {
  // stream
  for (int i = 0; i < LatSet::q; ++i) {
    ExcessC.getField(i).rotate(Delta_Index[i]);
  }
  // only collect inner cells
  for (int j = Geo.getOverlap(); j < Nj - Geo.getOverlap(); ++j) {
    for (int i = Geo.getOverlap(); i < Ni - Geo.getOverlap(); ++i) {
      std::size_t id = i + j * Ni;
      std::array<T *, LatSet::q> ExcessC_k = ExcessC.template getArray<T>(id);
      for (int k = 0; k < LatSet::q; ++k) {
        ExcessC_.get(id) += *(ExcessC_k[k]);
      }
    }
  }
  // clear ExcessC
  for (int i = 0; i < LatSet::q; ++i) {
    ExcessC.getField(i).Init(T(0));
  }
}

template <typename T, typename LatSet>
void BlockZhuStefanescu2D<T, LatSet>::SimpleCapture() {
  // only capture inner cells
  for (int j = Geo.getOverlap(); j < Nj - Geo.getOverlap(); ++j) {
    for (int i = Geo.getOverlap(); i < Ni - Geo.getOverlap(); ++i) {
      std::size_t id = i + j * Ni;
      if (util::isFlag(State.get(id), CAType::Fluid) && hasNeighborType(id, CAType::Solid)) {
        // Flag.SetField(id, CAFlag::toInterface);
        State.SetField(id, CAType::Interface);
        Velocity.SetField(id, Vector<T, 2>{});
      }
    }
  }
}

template <typename T, typename LatSet>
void BlockZhuStefanescu2D<T, LatSet>::TypeConversion() {
  for (std::size_t id = 0; id < N; ++id) {
    if (util::isFlag(Flag.get(id), CAFlag::toInterface)) {
      State.SetField(id, CAType::Interface);
    } else if (util::isFlag(Flag.get(id), CAFlag::toSolid)) {
      State.SetField(id, CAType::Solid);
    }
  }
  // clear flag
  Flag.getField(0).Init(CAFlag::None);
}

template <typename T, typename LatSet>
bool BlockZhuStefanescu2D<T, LatSet>::hasNeighborType(std::size_t id, std::uint8_t type) const {
  for (int i = 0; i < LatSet::q; ++i) {
    if (util::isFlag(State.get(id + Delta_Index[i]), type)) return true;
  }
  return false;
}

template <typename T, typename LatSet>
bool BlockZhuStefanescu2D<T, LatSet>::hasNeighborFlag(std::size_t id, std::uint8_t flag) const {
  for (int i = 0; i < LatSet::q; ++i) {
    if (util::isFlag(Flag.get(id + Delta_Index[i]), flag)) return true;
  }
  return false;
}

template <typename T, typename LatSet>
void BlockZhuStefanescu2D<T, LatSet>::communicate() {
  // same level communication
  ScalerField<T> &RhoF = lbmSO.getRhoField();
  for (BlockZSCommStru<T, LatSet> &comm : Communicators) {
    BlockZhuStefanescu<T, LatSet> *nZSBlock = comm.SendBlock;
    // const GenericArray<CAFlag> &nFlag = nZSBlock->getFlag().getField();
    const GenericArray<CAType> &nState = nZSBlock->getState().getField();
    const GenericArray<T> &nFs = nZSBlock->getFs().getField();
    const GenericArray<T> &nRho = nZSBlock->getLBMConc().getRhoField().getField();
    // state field, Fs field， rho field communication
    int size = comm.getRecvs().size();
    for (int i = 0; i < size; ++i) {
      std::size_t idrecv = comm.getRecvs()[i];
      std::size_t idsend = comm.getSends()[i];
      // Flag.SetField(idrecv, nFlag[idsend]);
      State.SetField(idrecv, nState[idsend]);
      Fs.SetField(idrecv, nFs[idsend]);
      RhoF.SetField(idrecv, nRho[idsend]);
    }
    // ExcessC field communication
    for (int k = 0; k < LatSet::q; ++k) {
      const CyclicArray<T> &nExcessCk = nZSBlock->getExcessC().getField(k);
      CyclicArray<T> &ExcessCk = ExcessC.getField(k);
      for (int i = 0; i < size; ++i) {
        std::size_t idrecv = comm.getRecvs()[i];
        std::size_t idsend = comm.getSends()[i];
        ExcessCk.set(idrecv, nExcessCk[idsend]);
      }
    }
  }
}

}  // namespace CA
