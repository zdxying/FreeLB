// grow.hh
#include "grow.h"

template <typename T, typename LatStru>
Grow2D<T, LatStru>::Grow2D(CAGField2D<T>& ca, GandinConverter<T>& convca,
                           LBManager2D<T>& lbm)
    : CA(ca),
      Geo(ca.Geo),
      ConvCA(convca),
      LBM(lbm),
      LatStru(ca.Geo.getNi(), ca.Geo.getNj()) {
  int max_BulkNum = (Geo.getNi() - 4) * (Geo.getNj() - 4);
  int max_nucs = ConvCA.Lattice_NucDens_Surf * ((Geo.getNi() + Geo.getNj() - 4) * 2 - 4) +
                 ConvCA.Lattice_NucDens_Bulk * (Geo.getNi() - 4) * (Geo.getNj() - 4);
  Cells.reserve(max_BulkNum);
  New_Nucs.reserve(max_nucs);
  Captureds.reserve(max_BulkNum);
  New_Growings.reserve(max_BulkNum);
  Lat_GrowthPara = ConvCA.Lattice_GrowthPara;

#ifdef _CAPTURE_P
  visited_ato = new std::atomic_flag[N];
  std::fill_n(visited_ato, N, ATOMIC_FLAG_INIT);
#endif
}

template <typename T, typename LatStru>
Grow2D<T, LatStru>::~Grow2D() {
#ifdef _CAPTURE_P
  delete[] visited_ato;
#endif
}

template <typename T, typename LatStru>
void Grow2D<T, LatStru>::Grow() {
  int Id = 0;
#pragma omp parallel for private(Id) num_threads(Thread_Num)
  for (int i = 0; i < Cells.size(); i++) {
    Id = Cells[i].Id;
    Cells[i].Arm += Get_GrowthVelo(Id);
  }
  // merge new_nucs to cells
  Cells.insert(Cells.end(), New_Nucs.begin(), New_Nucs.end());
}

template <typename T, typename LatStru>
T Grow2D<T, LatStru>::Get_GrowthVelo(int id) {
  T dT = ConvCA.get_LatTliq(LBM.getlbmSO_ref().getPoprho(id)) -
         LBM.getlbmTH_ref().getPoprho(id);
  // dT = dT > T(0) ? dT : T(0);
  if (dT <= T(0))
    return T(0);
  else
    return Lat_GrowthPara * dT * dT;
}

template <typename T, typename LatStru>
void Grow2D<T, LatStru>::Capture() {
  int Id = 0;
  int Id_nbr = 0;
  Captureds.clear();

#ifdef _CAPTURE_P
  // no repeated elements in Captureds
  std::fill_n(visited_ato, N, ATOMIC_FLAG_INIT);
  OMP_helper::Setup_Thread_Vec<gcell2D<T>>(Captureds_Th,
                                           Cell.size() * LatStru::q);
  // divide cells may lead to unbalanced workload
#pragma omp parallel for private(Id, Id_nbr) num_threads(Thread_Num)
  for (int i = 0; i < Cells.size(); i++) {
    Id = Cells[i].Id;
    for (int k = 0; k < LatStru::q; k++) {
      Id_nbr = Id + LatStru::Nbr[k];
      if (CA.State[Id_nbr] == 0 &&
          (!visited_ato[Id_nbr].test_and_set(std::memory_order_acquire)))
        CellCapture(Cells[i], Id_nbr);
    }
  }
  // merge captured cells
  OMP_helper::Merge_and_Clear<gcell2D<T>>(Captureds, Captureds_Th);
#else
  for (int i = 0; i < Cells.size(); i++) {
    Id = Cells[i].Id;
    for (int k = 0; k < LatStru::q; k++) {
      Id_nbr = LatStru::Nbr[k];
      if (CA.State[Id_nbr] == 0) {
        CellCapture(Cells[i], Id_nbr);
      }
    }
  }
#endif
  // merge captured cells to Cells
  Cells.insert(Cells.end(), Captureds.begin(), Captureds.end());
  // post capture, erase cells that all nbr cells are growing
  postCapture();
}

// get nbr state, check if all nbr cells are not growing(flag == 0 )
// if at least one nbr cell is not growing, return false else return true
template <typename T, typename LatStru>
inline bool Grow2D<T, LatStru>::all_growing(int id) {
  int nbrstate = 1;
  for (int i = 0; i < LatStru::q; i++) {
    if (CA.State[id + LatStru::Nbr[i]] == 0) return false;
  }
  return true;
}

template <typename T, typename LatStru>
void Grow2D<T, LatStru>::postCapture() {
#ifdef _CAPTURE_P
  OMP_helper::Divide_and_Setup<gcell2D<T>>(Cells, Cells_Th);
#pragma omp parallel for private(Id) num_threads(Thread_Num)
  for (int i = 0; i < Cells_Th.size(); i++) {
    postcapture_s(Cells_Th[i]);
  }
  OMP_helper::Merge_and_Clear<gcell2D<T>>(Cells, Cells_Th);
#else
  postcapture_s(Cells);
#endif
}

template <typename T, typename LatStru>
inline void Grow2D<T, LatStru>::postcapture_s(std::vector<gcell2D<T>>& cells) {
  auto iter = cells.begin();
  int Id = 0;
  while (iter != cells.end()) {
    Id = iter->Id;
    if (all_growing(Id)) {
      iter = cells.erase(iter);
      CA.Deactivate(Id);
    } else {
      iter++;
    }
  }
}

template <typename T, typename LatStru>
void Grow2D<T, LatStru>::CellCapture(gcell2D<T>& gcell, int id) {
  // get relative position regarding to growing centre
  T r_pos[2];
  Vect2D<T>::R_Loc(Geo.getVoxel(id)[0], Geo.getVoxel(Id)[1], gcell.x, gcell.y, CA.Orine[gcell.Id],
                   r_pos);
  // get quadrant of captured(or not) cell centre
  int quad = Vect2D<T>::quad(r_pos);
  // T ArmLen[2] = {gcell.Arm, gcell.Arm};  // ArmLen > 0
  T ArmLen = gcell.Arm;
  // line equation: x/a0 + y/a1 = 1 -> a1*x + a0*y = a0*a1
  if (fabs(r_pos[0] * ArmLen) + fabs(r_pos[1] * ArmLen) <= ArmLen * ArmLen) {
    // -------------------------captured ！--------------------------
    CA.Orine[id] = CA.Orine[gcell.Id];
    CA.Activate(id);
    // capture line equation: ArmLen[1]*x + ArmLen[0]*y = ArmLen[0]*ArmLen[1]
    // oppsite capture line: -ArmLen[1]*x - ArmLen[0]*y = ArmLen[0]*ArmLen[1]
    // get direction from cell centre(x,y) to capture line
    // dir = fabs(Ax + By + C) / sqrt(A^2 + B^2)
    // get signed ArmLen
    T ArmLen_s[2] = {ArmLen * quadrant[quad][0], ArmLen * quadrant[quad][1]};
    T dir0 = fabs(ArmLen_s[1] * r_pos[0] + ArmLen_s[0] * r_pos[1] -
                  ArmLen_s[0] * ArmLen_s[1]) /
             sqrt(Vect2D<T>::sqr(ArmLen_s));
    T dir1 = fabs(-ArmLen_s[1] * r_pos[0] + ArmLen_s[0] * r_pos[1] +
                  ArmLen_s[0] * ArmLen_s[1]) /
             sqrt(Vect2D<T>::sqr(ArmLen_s));
    // Compute new square size: Gandin's thesis
    T New_ArmLen =
        std::min(dir0 / sqrt(T(2)), T(1)) + std::min(dir1 / sqrt(T(2)), T(1));
    // ----get new growing centre of the captured cell
    // absolute position of arm0 and arm1
    T Arm0[2] = {gcell.x + ArmLen_s[0] * cos(CA.Orine[id]),
                 gcell.y + ArmLen_s[0] * sin(CA.Orine[id])};
    T Arm1[2] = {gcell.x + ArmLen_s[1] * cos(CA.Orine[id] + M_PI / T(2)),
                 gcell.y + ArmLen_s[1] * sin(CA.Orine[id] + M_PI / T(2))};
    // closest corner
    // T Len = dir0 <= dir1 ? ArmLen[0] : ArmLen[1];
    // new growing centre
    T x = gcell.x + (ArmLen - New_ArmLen) * cos(CA.Orine[id]);
    T y = gcell.y + (ArmLen - New_ArmLen) * sin(CA.Orine[id]);
#ifdef _CAPTURE_P
    Captureds_Th[omp_get_thread_num()].emplace_back(
        gcell2D<T>(id, x, y, New_ArmLen));
#else
    Captureds.emplace_back(gcell2D<T>(id, x, y, New_ArmLen));
#endif
    // -----------------captured end--------------
  }
}

// template <typename T>
// void Grow2D<T>::Captured(gcell2D<T>& cell, int id, int quad,
//                                  T* ArmLen) {
//   // CA.State[id] = 1;
//   // CA.Orine[id] = CA.Orine[cell.Id];
//   // Geo.getVoxel(id).getFlag() = 1;
//   // capture line equation: ArmLen[1]*x + ArmLen[0]*y = ArmLen[0]*ArmLen[1]
//   // oppsite capture line: ArmLen[1]*x + ArmLen[0]*y = -ArmLen[0]*ArmLen[1]
//   // get direction from cell centre(x,y) to capture line
//   // dir = (Ax + By + C) / sqrt(A^2 + B^2)

//   // project cell centre to capture line

//   // Captureds.emplace_back
// }

// LEGACY: this may lead to unbalanced workload
// template <typename T>
// template <typename LatStru>
// void Grow2D<T>::Capture(LatStru& Lattice) {
//   int Id = 0;
//   int nbr_id = 0;
//   Captureds.clear();
// #ifdef _CAPTURE_P
//   // no repeated elements in Captureds
//   std::fill_n(visited_ato, N, ATOMIC_FLAG_INIT);
//   // divide cells to sub-vectors for parallelization
//   std::vector<std::vector<gcell2D<T>>> Cells_Thread;
//   OMP_helper::Divide<gcell2D<T>>(Cells, Cells_Thread);
// #pragma omp parallel for private(Id, nbr_id) num_threads(Thread_Num)
//   for (int i = 0; i < Cells_Thread.size(); i++) {
//     auto iter = Cells_Thread[i].begin();
//     while (iter != Cells_Thread[i].end()) {
//       Id = iter->Id;
//       // try to catch neighbours
//       for (int k = 0; k < LatStru::q; k++) {
//         nbr_id = Id + LatStru::Nbr[k];
//         if (CA.State[nbr_id] == 0 &&
//             (!visited_ato[nbr_id].test_and_set(std::memory_order_acquire))) {
//           CellCapture(*iter, nbr_id);
//         }
//       }
//       if (CA.getNbrState(Id, Lattice) != 0) {
//         iter = Cells_Thread[i].erase(iter);
//         CA.Deactivate(Id);
//       } else
//         iter++;
//     }
//   }
//   // merge captured cells
//   OMP_helper::Merge<gcell2D<T>>(Captureds, Captureds_Th);
// #else
//   auto iter = Cells.begin();
//   while (iter != Cells.end()) {
//     Id = iter->Id;
//     // get neighbours
//     for (int k = 0; k < LatStru::q; k++) {
//       nbr_id = LatStru::Nbr[k];
//       if (CA.State[nbr_id] == 0) {
//         CellCapture(*iter, nbr_id);
//       }
//     }
//     // check if all nbr cells are captured
//     if (CA.getNbrState(Id, Lattice) != 0) {
//       iter = Cells.erase(iter);
//       CA.Deactivate(Id);
//     } else
//       iter++;
//   }
// #endif
//   // merge captured cells to Cells
//   Cells.insert(Cells.end(), Captureds.begin(), Captureds.end());
// }

// template <typename T>
// void Grow2D<T>::CellCapture(gcell2D<T>& gcell, int id) {
//   // get relative position regarding to growing centre
//   T[2] r_pos = Vect2D<T>::R_Loc(Geo.getVoxel(id)[0], Geo.getVoxel(Id)[1], gcell.x, gcell.y,
//                                 CA.Orine[gcell.Id]);
//   // get quadrant of captured(or not) cell centre
//   int quad = Vect2D<T>::quad(r_pos);
//   // T ArmLen[2] = {gcell.Arm, gcell.Arm};  // ArmLen > 0
//   T ArmLen = gcell.Arm;
//   // line equation: x/a0 + y/a1 = 1 -> a1*x + a0*y = a0*a1
//   if (fabs(r_pos[0] * ArmLen) + fabs(r_pos[1] * ArmLen) <= ArmLen * ArmLen) {
//     // -------------------------captured--------------------------
//     CA.Orine[id] = CA.Orine[gcell.Id];
//     CA.Activate(id);
//     // capture line equation: ArmLen[1]*x + ArmLen[0]*y = ArmLen[0]*ArmLen[1]
//     // oppsite capture line: -ArmLen[1]*x - ArmLen[0]*y = ArmLen[0]*ArmLen[1]
//     // get direction from cell centre(x,y) to capture line
//     // dir = fabs(Ax + By + C) / sqrt(A^2 + B^2)
//     // get signed ArmLen
//     T ArmLen_s[2] = {ArmLen * quadrant[quad][0], ArmLen * quadrant[quad][1]};
//     T dir0 = fabs(ArmLen_s[1] * r_pos[0] + ArmLen_s[0] * r_pos[1] -
//                   ArmLen_s[0] * ArmLen_s[1]) /
//              sqrt(Vect2D<T>::sqr(ArmLen_s));
//     T dir1 = fabs(-ArmLen_s[1] * r_pos[0] + ArmLen_s[0] * r_pos[1] +
//                   ArmLen_s[0] * ArmLen_s[1]) /
//              sqrt(Vect2D<T>::sqr(ArmLen_s));
//     // Compute new square size: Gandin's thesis
//     T New_ArmLen =
//         std::min(dir0 / sqrt(T(2)), T(1)) + std::min(dir1 / sqrt(T(2)),
//         T(1));
//     // ----get new growing centre of the captured cell
//     // absolute position of arm0 and arm1
//     T Arm0[2] = {gcell.x + ArmLen_s[0] * cos(CA.Orine[id]),
//                  gcell.y + ArmLen_s[0] * sin(CA.Orine[id])};
//     T Arm1[2] = {gcell.x + ArmLen_s[1] * cos(CA.Orine[id] + M_PI / T(2)),
//                  gcell.y + ArmLen_s[1] * sin(CA.Orine[id] + M_PI / T(2))};
//     // closest corner
//     // T Len = dir0 <= dir1 ? ArmLen[0] : ArmLen[1];
//     // new growing centre
//     T x = gcell.x + (ArmLen - New_ArmLen) * cos(CA.Orine[id]);
//     T y = gcell.y + (ArmLen - New_ArmLen) * sin(CA.Orine[id]);
// #ifdef _CAPTURE_P
//     Captureds_Th[omp_get_thread_num()].emplace_back(
//         gcell2D<T>(id, x, y, New_ArmLen));
// #else
//     Captureds.emplace_back(gcell2D<T>(id, x, y, New_ArmLen));
// #endif
//     // -----------------captured end--------------
//   }
// }