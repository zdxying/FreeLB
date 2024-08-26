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

#include "ca/gandin.h"

template <typename T, typename LatStru>
GandinCA2D<T, LatStru>::GandinCA2D(CAGField2D<T> &ca,
                                   GandinConverter<T> &convca,
                                   LBManager2D<T> &lbm)
    : CA(ca),
      ConvCA(convca),
      Geo(ca.Geo),
      LBM(lbm),
      message(4),
      LatStru(ca.Geo.getNi(), ca.Geo.getNj()),
      N(ca.Geo.getN()),
      _rd(),
      GetGCells2D(ca.Geo.getNi(), ca.Geo.getNj(), ca.State) {
  int max_nucs =
      ConvCA.Lattice_NucDens_Surf * ((Geo.getNi() + Geo.getNj() - 4) * 2 - 4) +
      ConvCA.Lattice_NucDens_Bulk * (Geo.getNi() - 4) * (Geo.getNj() - 4);
  Interface.reserve(2 * (Geo.getNi() + Geo.getNj()));
  Cells.reserve(2 * (Geo.getNi() + Geo.getNj()));
#ifdef _OPENMP
  IndexStart_OMP.reserve(Thread_Num);
  IndexEnd_OMP.reserve(Thread_Num);
  OMP_helper::Setup_Thread_Vec<Gcell2D<T>>(Nucleus_OMP, max_nucs);
  OMP_helper::Setup_Thread_Vec<Gcell2D<T>>(Captureds_OMP, max_nucs);
  _Atomic_Visited = new std::atomic_flag[N];
  ResetVisited(N);
#else
  _gen = std::mt19937(_rd());
#endif
  // print
  // initial undercooling = Tliq - Tinit
  T Init_DT = ConvCA.get_Tliq(ConvCA.ConcConv.CInit) - ConvCA.TempConv.TInit;
  T Init_LatDT = ConvCA.get_LatTliq(ConvCA.ConcConv.Lattice_CInit) -
                 ConvCA.TempConv.Lattice_TInit;
  std::cout << "[Gandin CA]: " << std::endl;
  std::cout << "Init Undercooling: " << Init_DT << " K\n"
            << "Init Lattice Undercooling: " << Init_LatDT << std::endl;
  std::cout << "Gandin CA setup done!" << std::endl;
}

template <typename T, typename LatStru>
inline GandinCA2D<T, LatStru>::~GandinCA2D() {
#ifdef _OPENMP
  delete[] _Atomic_Visited;
#endif
}

template <typename T, typename LatStru>
void GandinCA2D<T, LatStru>::Nucleation() {
  // Set up normal distribution (C++11) and random number generator,
  // critical undercooling of PreNucs is generated by normal distribution
  std::normal_distribution<T> Gauss_Bulk(ConvCA.Lattice_DT_Mean_Bulk,
                                         ConvCA.Lattice_DT_Std_Bulk);
  std::normal_distribution<T> Gauss_Surf(ConvCA.Lattice_DT_Mean_Surf,
                                         ConvCA.Lattice_DT_Std_Surf);

  std::vector<int> &RemainBulks = this->Get_Bulks();
  std::vector<int> &RemainSurfs = this->Get_Surfs();

  int RemainBulkNum = RemainBulks.size();
  int RemainsurfNum = RemainSurfs.size();
  // Bulk cell
  if (RemainBulkNum == 0) {
    msg_out(0, "All Bulk Cells are Solified");
  } else {
    int Bulk_Sites = round(ConvCA.Lattice_NucDens_Bulk * RemainBulkNum + 0.5);
    if (Bulk_Sites == 0) {
      msg_out(1, "No NucSite in Bulk!");
      Accum_Bulk_Sites += ConvCA.Lattice_NucDens_Bulk * RemainBulkNum;
      if (Accum_Bulk_Sites >= 1) {
        Accum_Bulk_Sites = 0;
        Nucleation_s(RemainBulks, Gauss_Bulk, 1, Cells);
      }
    } else {
#ifdef _OPENMP
      // visited atomic, the array size is N
      // but only the 0th to (RemainIdx.size() - 1)th elements will be used
      ResetVisited(RemainBulkNum);
      // omp
      int i = 0;
      std::vector<int> Bulk_Sites_OMP;
      OMP_helper::Divide_Num(Bulk_Sites_OMP, Bulk_Sites);
#pragma omp parallel for private(i, Gauss_Bulk) num_threads(Thread_Num)
      for (i = 0; i < Bulk_Sites_OMP.size(); i++) {
        Nucleation_omp(RemainBulks, Gauss_Bulk, Bulk_Sites_OMP[i],
                       Nucleus_OMP[i]);
      }
      // merge Nucleus_OMP to Cells
      OMP_helper::Merge<Gcell2D<T>>(Cells, Nucleus_OMP);
#else
      Nucleation_s(RemainBulks, Gauss_Bulk, Bulk_Sites, Cells);
#endif
    }
  }

  // surface cell
  if (RemainsurfNum == 0) {
    msg_out(2, "All Surface Cells are Solified");
  } else {
    int Surf_Sites = round(ConvCA.Lattice_NucDens_Surf * RemainsurfNum + 0.5);
    if (Surf_Sites == 0) {
      msg_out(3, "No NucSite in Surface!");
      Accum_Surf_Sites += ConvCA.Lattice_NucDens_Surf * RemainsurfNum;
      if (Accum_Surf_Sites >= 1) {
        Accum_Surf_Sites = 0;
        Nucleation_s(RemainSurfs, Gauss_Surf, 1, Cells);
      }
    } else {
      Nucleation_s(RemainSurfs, Gauss_Surf, Surf_Sites, Cells);
    }
  }
}

template <typename T, typename LatStru>
void GandinCA2D<T, LatStru>::Nucleation_s(
    const std::vector<int> &RemainIdx,
    std::normal_distribution<T> &GaussDistrib, int NucNum,
    std::vector<Gcell2D<T>> &Nucs) {
  // uniform distribution
  std::uniform_int_distribution<> dis(0, RemainIdx.size() - 1);

  int Cellid = 0;
  T dT_critic = 0;

  for (int i = 0; i < NucNum; i++) {
    Cellid = RemainIdx[dis(_gen)];
    if (CA.State[Cellid] == 0) {
      dT_critic = GaussDistrib(_gen);
      dT_critic = dT_critic < 0 ? 0 : dT_critic;
      if (isNucleated(Cellid, dT_critic)) {
        CA.Activate(Cellid);
        CA.Orine[Cellid] = (std::rand() % 90) * M_PI / T(180);
        Nucs.emplace_back(Cellid, Geo.getVoxel(Cellid)[0],
                          Geo.getVoxel(Cellid)[1]);
      }
    }
  }
}

template <typename T, typename LatStru>
void GandinCA2D<T, LatStru>::Nucleation_omp(
    const std::vector<int> &RemainIdx,
    std::normal_distribution<T> &GaussDistrib, int NucNum,
    std::vector<Gcell2D<T>> &Nucs) {
  Nucs.clear();
  // random number generator
  std::mt19937 gen = std::mt19937(_rd());
  // uniform distribution
  std::uniform_int_distribution<> dis(0, RemainIdx.size() - 1);

  int index = 0;
  int Cellid = 0;
  T dT_critic = 0;

  for (int i = 0; i < NucNum; i++) {
    index = dis(gen);
    Cellid = RemainIdx[index];
    if (CA.State[Cellid] == 0) {
      dT_critic = GaussDistrib(gen);
      dT_critic = dT_critic < 0 ? 0 : dT_critic;
      if (isNucleated(Cellid, dT_critic)) {
#ifdef _OPENMP
        if (!_Atomic_Visited[index].test_and_set(std::memory_order_acquire))
#endif
        {
          CA.Activate(Cellid);
          CA.Orine[Cellid] = (std::rand() % 90) * M_PI / T(180);
          Nucs.emplace_back(Cellid, Geo.getVoxel(Cellid)[0],
                            Geo.getVoxel(Cellid)[1]);
        }
      }
    }
  }
}

template <typename T, typename LatStru>
inline bool GandinCA2D<T, LatStru>::isNucleated(int Cellid, T dT_critic) {
  return LBM.getlbmTH_ref().getPoprho(Cellid) <
         (ConvCA.get_LatTliq(LBM.getlbmSO_ref().getPoprho(Cellid)) - dT_critic);
}

template <typename T, typename LatStru>
void GandinCA2D<T, LatStru>::Single_Nucleation(int latx, int laty, T orine) {
  int Id = Index2D::GetId(latx, laty, Geo.getNi());
  CA.State[Id] = 1;
  CA.Orine[Id] = orine;
  Geo.getVoxel(Id).getFlag() = 1;
  Cells.emplace_back(Id, Geo.getVoxel(Id)[0], Geo.getVoxel(Id)[1]);
}

//------------------- grow -----------------

template <typename T, typename LatStru>
void GandinCA2D<T, LatStru>::Grow() {
  int size = Cells.size();
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < size; i++) {
    Gcell2D<T> &cell = Cells[i];
    cell.getArm() += Get_GrowthVelo(cell.getCellId());
  }
}

template <typename T, typename LatStru>
inline T GandinCA2D<T, LatStru>::Get_GrowthVelo(int id) {
  T dT = ConvCA.get_LatTliq(LBM.getlbmSO_ref().getPoprho(id)) -
         LBM.getlbmTH_ref().getPoprho(id);
  // dT = dT > T(0) ? dT : T(0);
  if (dT <= T(0))
    return T(0);
  else
    return ConvCA.Lattice_GrowthPara * dT * dT;
}

// ------------------capture------------------

template <typename T, typename LatStru>
void GandinCA2D<T, LatStru>::Capture() {
#ifdef _OPENMP
  if (Cells.size() > 100 * Thread_Num) {
    ResetVisited(N);
    OMP_helper::Divide_Index(IndexStart_OMP, IndexEnd_OMP, Cells.size());
    int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < Captureds_OMP.size(); ++i) {
      Capture_s(Captureds_OMP[i], IndexEnd_OMP[i], IndexStart_OMP[i]);
    }
    // add to cells
    OMP_helper::Merge_and_Clear<Gcell2D<T>>(Cells, Captureds_OMP);
  } else {
    Capture_s(Cells, Cells.size());
  }
#else
  Capture_s(Cells, Cells.size());
#endif
  // post capture, erase cells that all nbr cells are growing
  postCapture();
}

template <typename T, typename LatStru>
void GandinCA2D<T, LatStru>::Capture_s(std::vector<Gcell2D<T>> &captureds,
                                       int end, int start) {
  int id_nbr = 0;

  for (int i = start; i < end; ++i) {
    Gcell2D<T> &cell = Cells[i];
    for (int k = 0; k < LatStru::q; ++k) {
      int id_nbr = cell.getCellId() + LatStru::Nbr[k];
      if (CA.State[id_nbr] == 0) {
        CellCapture(cell, id_nbr, captureds);
      }
    }
  }
}

template <typename T, typename LatStru>
void GandinCA2D<T, LatStru>::CellCapture(Gcell2D<T> &gcell,
                                         std::vector<Gcell2D<T>> &captureds,
                                         int id) {
  // get relative position regarding to growing centre
  T r_pos[2];
  Vect2D<T>::R_Loc(Geo.getVoxel(id)[0], Geo.getVoxel(id)[1], gcell[0], gcell[1],
                   CA.Orine[gcell.getCellId()], r_pos);
  // get quadrant of captured(or not) cell centre
  int quad = Vect2D<T>::quad(r_pos);
  // T ArmLen[2] = {gcell.Arm, gcell.Arm};  // ArmLen > 0
  T ArmLen = gcell.getArm();
  // line equation: x/a0 + y/a1 = 1 -> a1*x + a0*y = a0*a1
  if (fabs(r_pos[0] * ArmLen) + fabs(r_pos[1] * ArmLen) <= ArmLen * ArmLen) {
    // -------------------------captured ！--------------------------
    CA.Orine[id] = CA.Orine[gcell.getCellId()];
    // wrong line equation!
    // capture line equation: ArmLen[1]*x + ArmLen[0]*y = ArmLen[0]*ArmLen[1]
    // oppsite capture line: -ArmLen[1]*x - ArmLen[0]*y = ArmLen[0]*ArmLen[1]
    // get distance from cell centre(x,y) to capture line
    // dist = fabs(Ax + By + C) / sqrt(A^2 + B^2)
    // get signed ArmLen
    T ArmLen_s[2] = {ArmLen * quadrant[quad][0], ArmLen * quadrant[quad][1]};
    T dist0 = fabs(ArmLen_s[1] * r_pos[0] + ArmLen_s[0] * r_pos[1] -
                   ArmLen_s[0] * ArmLen_s[1]) /
              std::sqrt(Vect2D<T>::sqr(ArmLen_s));
    T dist1 = fabs(-ArmLen_s[1] * r_pos[0] + ArmLen_s[0] * r_pos[1] +
                   ArmLen_s[0] * ArmLen_s[1]) /
              std::sqrt(Vect2D<T>::sqr(ArmLen_s));
    // Compute new square size: Gandin's thesis
    T New_ArmLen =
        std::min(dist0 / std::sqrt(T(2)), T(1)) + std::min(dist1 / std::sqrt(T(2)), T(1));
    // ----get new growing centre of the captured cell
    // absolute position of arm0 and arm1
    T Arm0[2] = {gcell[0] + ArmLen_s[0] * std::cos(CA.Orine[id]),
                 gcell[1] + ArmLen_s[0] * std::sin(CA.Orine[id])};
    T Arm1[2] = {gcell[0] + ArmLen_s[1] * std::cos(CA.Orine[id] + M_PI / T(2)),
                 gcell[1] + ArmLen_s[1] * std::sin(CA.Orine[id] + M_PI / T(2))};
    // closest corner
    // T Len = dist0 <= dist1 ? ArmLen[0] : ArmLen[1];
    // new growing centre
    T x = gcell[0] + (ArmLen - New_ArmLen) * std::cos(CA.Orine[id]);
    T y = gcell[1] + (ArmLen - New_ArmLen) * std::sin(CA.Orine[id]);
#ifdef _OPENMP
    if (!_Atomic_Visited[id].test_and_set(std::memory_order_acquire))
#endif
    {
      captureds.emplace_back(Gcell2D<T>(id, x, y, New_ArmLen));
      CA.Activate(id);
    }
  }
}

template <typename T, typename LatStru>
void GandinCA2D<T, LatStru>::CellCapture(Gcell2D<T> &cell, int id,
                                         std::vector<Gcell2D<T>> &captureds) {
  // relative position regarding to growing centre of a square
  // the rotation angle is between 0 and Pi/2 degree(rad),
  // the x and y axis of the new coordinate system are the two perpendicular
  // diagonals of the square.
  Vector<T, 2> RLoc =
      getRLoc(Geo.getVoxel(id), cell, CA.Orine[cell.getCellId()]);
  // Vector<T, 2> RLoc;
  // getRLoc(Geo.getVoxel(id), cell, RLoc, CA.Orine[cell.getCellId()]);

  // in the new coordinate system, get Intercept2D of cell centre to be captured
  // which is used to calculate the intercept of the capture line and the axis
  // interceptx=Intercept2D[quad][0]*ArmLen;
  // intercepty=Intercept2D[quad][1]*ArmLen
  int quad = getQuad2D(RLoc);
  // check if captured by 4 edges of the growing square
  T ArmLen = cell.getArm();
  if (fabs(RLoc[0] * ArmLen) + fabs(RLoc[1] * ArmLen) <= std::pow(ArmLen, 2)) {
    // captured！
    CA.Orine[id] = CA.Orine[cell.getCellId()];
    // get the corresponding capture line of [quad-1] and [quad+1]
    int quad_1 = quad == 0 ? 3 : quad - 1;
    int quad_2 = quad == 3 ? 0 : quad + 1;
    // get distance from cell centre to capture line
    // capture line equation: x/ax + y/ay = 1 -> ay*x + ax*y - ax*ay = 0
    // dist = fabs(Ax + By - A*B) / sqrt(A^2 + B^2), where |A| = |B| = ArmLen
    // dist = fabs(A/ArmLen * x + B/ArmLen * y - A*B/ArmLen) / sqrt(2)
    T Dist1 = fabs(Intercept2D[quad_1][1] * RLoc[0] +
                   Intercept2D[quad_1][0] * RLoc[1] -
                   Intercept2D[quad_1][0] * Intercept2D[quad_1][1] * ArmLen) /
              std::sqrt(T(2));
    T Dist2 = fabs(Intercept2D[quad_2][1] * RLoc[0] +
                   Intercept2D[quad_2][0] * RLoc[1] -
                   Intercept2D[quad_2][0] * Intercept2D[quad_2][1] * ArmLen) /
              std::sqrt(T(2));
    // Compute new square's arm length: Gandin's thesis
    T New_ArmLen =
        std::min(Dist1 / std::sqrt(T(2)), T(1)) + std::min(Dist2 / std::sqrt(T(2)), T(1));
    // get the nearest square corner(vertex) to the cell centre
    T minDist = GetDist2(RLoc, Vertex2D[0] * ArmLen);
    int Vertex2D_Id = 0;
    for (int i = 1; i < 4; i++) {
      T dist = GetDist2(RLoc, Vertex2D[i] * ArmLen);
      if (dist < minDist) {
        minDist = dist;
        Vertex2D_Id = i;
      }
    }
    // get new growing square centre of the captured cell
    // get relative location
    Vector<T, 2> RLoc_new = Vertex2D[Vertex2D_Id] * (ArmLen - New_ArmLen);
    // get absolute(global) location
    Vector<T, 2> Loc_new = getGLoc(RLoc_new, cell, CA.Orine[cell.getCellId()]);
#ifdef _OPENMP
    if (!_Atomic_Visited[id].test_and_set(std::memory_order_acquire))
#endif
    {
      captureds.emplace_back(Gcell2D<T>(id, Loc_new, New_ArmLen));
      CA.Activate(id);
    }
  }
}

template <typename T, typename LatStru>
void GandinCA2D<T, LatStru>::postCapture() {
#ifdef _OPENMP
  if (Cells.size() > 1000 * Thread_Num) {
    std::vector<std::vector<Gcell2D<T>>> Cells_OMP;
    OMP_helper::Divide_and_Setup<Gcell2D<T>>(Cells, Cells_OMP);
    int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < Cells_OMP.size(); i++) {
      postcapture_s(Cells_OMP[i]);
    }
    OMP_helper::Merge<Gcell2D<T>>(Cells, Cells_OMP);
  } else {
    postcapture_s(Cells);
  }
#else
  postcapture_s(Cells);
#endif
}

template <typename T, typename LatStru>
inline void GandinCA2D<T, LatStru>::postcapture_s(
    std::vector<Gcell2D<T>> &cells) {
  auto iter = cells.begin();
  int Id = 0;
  while (iter != cells.end()) {
    Id = iter->getCellId();
    if (all_growing(Id)) {
      iter = cells.erase(iter);
      CA.Deactivate(Id);
    } else {
      iter++;
    }
  }
}

// get nbr state, check if all nbr cells are not growing(flag == 0 )
// if at least one nbr cell is not growing, return false else return true
template <typename T, typename LatStru>
inline bool GandinCA2D<T, LatStru>::all_growing(int id) {
  for (int i = 0; i < LatStru::q; i++) {
    if (CA.State[id + LatStru::Nbr[i]] == 0) return false;
  }
  return true;
}
