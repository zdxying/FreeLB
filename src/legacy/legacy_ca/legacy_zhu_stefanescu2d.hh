// Zhu-Stefanescu (Z-S) model for dendrite growth, implementations
#include "legacy/legacy_ca/legacy_zhu_stefanescu2d.h"
#include "utils/exception.h"

template <typename T, typename LatStru>
inline ZhuStefanescu2D<T, LatStru>::ZhuStefanescu2D(CAZSField2D<T> &ca,
                                                    ZSConverter<T> &convca,
                                                    LBManager2D<T> &lbm,
                                                    T delta_, T theta,
                                                    int siteid, int num)
    : CA(ca),
      ConvCA(convca),
      lbmSO(lbm.getlbmSO_ref()),
      lbmTH(lbm.getlbmTH_ref()),
      Geo(ca.Geo),
      FDM2D<T>(ca.Geo.getNi(), ca.Geo.getNj(), ca.f),
      LatStru(ca.Geo.getNi(), ca.Geo.getNj()),
      delta(delta_),
      Ni(ca.Geo.getNi()),
      Nj(ca.Geo.getNj()),
      N(ca.Geo.getNi() * ca.Geo.getNj()),
      Theta(theta),
      SiteId(siteid),
      GT(convca.Lattice_GT_Coef),
      C0(lbm.getlbmSO_ref().getlatrhoinit()),
      Tl_eq(convca.get_LatTliq(lbm.getlbmSO_ref().getlatrhoinit())),
      m_l(convca.Lattice_m_Liq),
      Part_Coef(convca.Part_Coef),
      _Part_Coef(1 - convca.Part_Coef),
      Field(lbm.getUField()) {
  Interface.reserve(2 * (Ni + Nj));
  Cells.reserve(2 * (Ni + Nj));
  PreCaps.reserve(2 * (Ni + Nj));
  Polygon.reserve(4);
  Nbr_Pos.reserve(LatStru::q - 1);
#ifdef _OPENMP
  IndexStart_OMP.reserve(Thread_Num);
  IndexEnd_OMP.reserve(Thread_Num);
  OMP_helper::Setup_Thread_Vec<int>(PreCaps_OMP, 2 * (Ni + Nj));
  OMP_helper::Setup_Thread_Vec<ZScell2D<T>>(Captureds_OMP, 2 * (Ni + Nj));
  _Atomic_Visited = new std::atomic_flag[N];
  ResetVisited(N);
#endif
  // setup
  Setup(SiteId, num);
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::Setup(int id, int num) {
  std::cout << "[Zhu-Stefanescu CA]" << std::endl;
  // get preferred growth angle to x-axis
  std::cout << "preferred growth angle: " << Theta / M_PI << " Pi" << std::endl;
  // get initial undercooling
  T Tl = lbmTH.getlatrhoinit();
  T deltaT = Tl_eq - Tl;
  std::cout << "Initial undercooling: " << deltaT << " | "
            << ConvCA.TempConv.getPhysDTemp(deltaT) << std::endl;
  int x = id % Ni;
  int y = id / Ni;
  // set to solid phase
  if (CA.State[id] == -1) {
    std::cout << "Warning: Setup at (" << x << ", " << y
              << ") failed, State = -1" << std::endl;
  } else {
    CA.Deactivate(id);
    CA.f[id] = 1;
    lbmSO.setPopRho(id, lbmSO.getPoprho(id) * Part_Coef);
    std::cout << "Setup at (" << x << ", " << y << "), id = " << id
              << " succeeded" << std::endl;
  }

  // num neighbors
  int count = 0;
  if (num > LatStru::q) {
    std::cout << "Warning: Setup at (" << x << ", " << y
              << ") failed, num > LatStru::q" << std::endl;
    num = LatStru::q;
  }
  for (int i = 0; i < num; i++) {
    int idn = id + LatStru::Nbr[i];
    if (CA.State[idn] == 0) {
      CA.Activate(idn);
      Cells.emplace_back(idn);
      count++;
      // if fluid field is calculated, set U to 0
      Field.DeActivate(idn);
    }
  }
  if (count != num) {
    std::cout << "Warning: Incomplete neighbors setup. " << count << "/" << num
              << std::endl;
  }
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::EraseInactive() {
#ifdef _OPENMP
  if (Cells.size() > 1000 * Thread_Num) {
    // parallel, divide Cells into sub-vectors
    std::vector<std::vector<ZScell2D<T>>> Cells_Th;
    OMP_helper::Divide_and_Setup<ZScell2D<T>>(Cells, Cells_Th);
    int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < Cells_Th.size(); i++) {
      EraseInactive_s(Cells_Th[i]);
    }
    // merge
    Cells.clear();
    OMP_helper::Merge<ZScell2D<T>>(Cells, Cells_Th);
  } else {
    EraseInactive_s(Cells);
  }
#else
  EraseInactive_s(Cells);
#endif
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::EraseInactive_s(
    std::vector<ZScell2D<T>> &cells) {
  // erase inactive cells from cells
  cells.erase(std::remove_if(cells.begin(), cells.end(),
                             [this](ZScell2D<T> &cell) {
                               return CA.State[cell.Id] == -1;
                             }),
              cells.end());
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::Grow() {
  // traverse all interface cells
  int id = 0;
  T deltaf = 0;
  Solifieds.clear();
  // get deltaf for all interface cells
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < Cells.size(); i++) {
    getdeltaf(Cells[i]);
  }
  // add deltaf to CAfield.f
  // #pragma omp parallel for private(id, deltaf) collapse(1)
  // num_threads(Thread_Num)
  for (int i = 0; i < Cells.size(); i++) {
    ZScell2D<T> &cell = Cells[i];
    id = cell.Id;
    deltaf = cell.deltaf;
    CA.f[id] += deltaf;
    // if solified, set to solid phase
    if (CA.f[id] >= 1) {
      // revised_deltaf = 1 - CA.f[id]old = 1 - (CA.f[id]new - deltaf)
      deltaf += 1 - CA.f[id];
      CA.f[id] = 1;
      // cell.Addto_Cs(Part_Coef * deltaf * getAveNbrPopRho(id), deltaf);
      cell.Addto_Cs(Part_Coef * deltaf * lbmSO.getPoprho(id), deltaf);
      // set pop in lbmSo
      lbmSO.setPopRho(id, cell.C_Solid());
      // add to newly solifieds
      Solifieds.emplace_back(id, deltaf);
      // // calculate rejected solute
      // RejectSolutetoNbr(id, deltaf);
      // // deactivate
      // CA.Deactivate(id);
      // // if fluid field is calculated, set U to 0
      // Field.DeActivate(id);
    } else {
      // update solid phase composition, given by:
      // sum(k*deltaf*C_l)/sum(deltaf)
      // cell.Addto_Cs(Part_Coef * deltaf * getAveNbrPopRho(id), deltaf);
      cell.Addto_Cs(Part_Coef * deltaf * lbmSO.getPoprho(id), deltaf);
      // calculate rejected solute
      // lbmSO.addtoRho(id, _Part_Coef * deltaf * getAveNbrPopRho(id));
      // lbmSO.addtoRho(id, _Part_Coef * deltaf * lbmSO.getPoprho(id));
      lbmSO.addtoPop(id, _Part_Coef * deltaf * lbmSO.getPoprho(id));
      // RejectSolute(id, deltaf);
    }
  }
  // deal with newly solified cells
  Solified_s(Solifieds);
}

// TODO: Enable parallel
template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::Solified_s(
    const std::vector<SolifiedZScell2D<T>> &solifieds) {
  int id = 0;
  auto it = solifieds.begin();
  while (it != solifieds.end()) {
    id = it->Id;
    RejectSolutetoNbr(id, it->deltaf);
    CA.Deactivate(id);
    // set U to 0
    Field.DeActivate(id);
    it++;
  }
}
template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::Post_Grow() {
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < Cells.size(); i++) {
    lbmSO.compute_rho(Cells[i].Id);
    if (all_NbrSolid(Cells[i].Id)) CA.Deactivate(Cells[i].Id);
  }
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::Grow_erase() {
  if (Cells.size() > 100 * Thread_Num) {
    // parallel, divide Cells into sub-vectors
    std::vector<std::vector<ZScell2D<T>>> Cells_Th;
    OMP_helper::Divide_and_Setup<ZScell2D<T>>(Cells, Cells_Th);
    int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < Cells_Th.size(); i++) {
      Grow_erase_s(Cells_Th[i]);
    }
    // merge
    Cells.clear();
    OMP_helper::Merge<ZScell2D<T>>(Cells, Cells_Th);
  } else {
    Grow_erase_s(Cells);
  }
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::Grow_erase_s(
    std::vector<ZScell2D<T>> &cells) {
  // traverse all elements in cells
  int id = 0;
  T deltaf = 0;
  // get deltaf for all interface cells
  for (int i = 0; i < cells.size(); i++) {
    getdeltaf(cells[i]);
  }
  // iterator
  auto it = cells.begin();
  while (it != cells.end()) {
    ZScell2D<T> &cell = *it;
    id = cell.Id;
    deltaf = cell.deltaf;
    CA.f[id] += deltaf;
    // if solified, set to solid phase
    if (CA.f[id] > 1) {
      // revised deltaf = 1 - CA.f[id]old = 1 - (CA.f[id]new - deltaf)
      deltaf += 1 - CA.f[id];
      CA.f[id] = 1;
      cell.Addto_Cs(Part_Coef * deltaf * lbmSO.getPoprho(id), deltaf);
      // set pop in lbmSo
      lbmSO.setPopRho(id, cell.C_Solid());
      // calculate rejected solute
      RejectSolutetoNbr(id, deltaf);
      // deactivate
      CA.Deactivate(id);
      // remove from cells
      it = cells.erase(it);
    } else {
      // update solid phase composition, given by:
      // sum(k*deltaf*C_l)/sum(deltaf)
      cell.Addto_Cs(Part_Coef * deltaf * lbmSO.getPoprho(id), deltaf);
      // calculate rejected solute
      lbmSO.addtoRho(id, _Part_Coef * deltaf * lbmSO.getPoprho(id));
      it++;
    }
  }
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::SimpleCapture() {
  int idn = 0;
  int size = Cells.size();
  for (int i = 0; i < size; i++) {
    // get neighbors
    for (int j = 0; j < LatStru::q; j++) {
      idn = Cells[i].Id + LatStru::Nbr[j];
      if (CA.State[idn] == 0 && has_SolidNbr(idn)) {
        CA.Activate(idn);
        Cells.emplace_back(idn);
        // if fluid field is calculated, set U to 0
        Field.DeActivate(idn);
      }
    }
  }
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::PreCapture_s(std::vector<int> &precaps) {
  precaps.clear();
  int idn = 0;
  int size = Cells.size();
  for (int i = 0; i < size; i++) {
    // get neighbors
    for (int j = 0; j < LatStru::q; j++) {
      idn = Cells[i].Id + LatStru::Nbr[j];
      if (CA.State[idn] == 0 && has_SolidNbr(idn)) {
        // pesudo capture, if not captured, then set State to 0
        CA.State[idn] = 2;
        precaps.emplace_back(idn);
      }
    }
  }
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::PreCapture_omp(std::vector<int> &precaps,
                                                 int start, int end) {
  precaps.clear();
  int idn = 0;
  for (int i = start; i < end; i++) {
    // get neighbors
    for (int j = 0; j < LatStru::q; j++) {
      idn = Cells[i].Id + LatStru::Nbr[j];
      if (CA.State[idn] == 0 && has_SolidNbr(idn)) {
#ifdef _OPENMP
        if (!_Atomic_Visited[idn].test_and_set(std::memory_order_acquire))
#endif
        {
          // pesudo capture, if not captured, then set State to 0
          CA.State[idn] = 2;
          precaps.emplace_back(idn);
        }
      }
    }
  }
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::SLICapture() {
  // erase of inactive cells from Cells will be done after Capture()
  // DO NOT erase them before Capture()

#ifdef _OPENMP
  if (Cells.size() > 100 * Thread_Num) {
    // parallel
    ResetVisited(N);
    OMP_helper::Divide_Index(IndexStart_OMP, IndexEnd_OMP, Cells.size());
    int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < IndexStart_OMP.size(); i++) {
      PreCapture_omp(PreCaps_OMP[i], IndexStart_OMP[i], IndexEnd_OMP[i]);
    }
    // end of parallel
    ////serial
    // PreCapture_s(PreCaps);
    //// divide PreCaps into sub-vectors
    // OMP_helper::Divide<int>(PreCaps, PreCaps_OMP);
    // end of serial
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < PreCaps_OMP.size(); i++) {
      SLICapture_s(PreCaps_OMP[i], Captureds_OMP[i]);
    }
    // add to Cells
    OMP_helper::Merge_and_Clear<ZScell2D<T>>(Cells, Captureds_OMP);
  } else {
    PreCapture_s(PreCaps);
    SLICapture_s(PreCaps, Cells);
  }
#else
  PreCapture_s(PreCaps);
  SLICapture_s(PreCaps, Cells);
#endif
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::SLICapture_s(
    const std::vector<int> &precaps, std::vector<ZScell2D<T>> &captureds) {
  int id, locid_, locid1, locid2, id1, id2;
  std::vector<int> nbrs_pos;
  nbrs_pos.reserve(4);
  for (int i = 0; i < precaps.size(); i++) {
    id = precaps[i];
    // get position of neighbors with State = -1 (solid cells)
    nbrs_pos.clear();
    getNbrs_pos(id, nbrs_pos);
    bool captured = false;
    // create polygonal on 2 neighbors of solid cells
    for (int j = 0; j < nbrs_pos.size(); j++) {
      // local id of solid cell: id_ = id + LatStru::Nbr[locid_]
      // get local index of solid cell
      locid_ = nbrs_pos[j];
      // get local index of 2 neighbors of solid cell
      locid1 = (locid_ + 1) % 8;
      locid2 = (locid_ + 7) % 8;
      // get global index of 2 neighbors of solid cell
      id1 = id + LatStru::Nbr[locid1];
      id2 = id + LatStru::Nbr[locid2];
      // (create quadrilaterals): Cell[i] -> Vert[i] ->Vert[j] ->Cell[j]
      Polygon.emplace_back(Geo.getVoxel(id1)[0], Geo.getVoxel(id1)[1]);
      Polygon.emplace_back(getVertice(id1));
      Polygon.emplace_back(getVertice(id2));
      Polygon.emplace_back(Geo.getVoxel(id2)[0], Geo.getVoxel(id2)[1]);
      // check if centre of cell id is inside the polygonal
      Vector2D<T> point(Geo.getVoxel(id)[0], Geo.getVoxel(id)[1]);
      if (InPolygon(point, Polygon)) {
        captured = true;
        break;
      }
    }
    if (captured) {
      CA.Activate(id);
      captureds.emplace_back(id);
      // if fluid field is calculated, set U to 0
      // Field.DeActivate(id);
    } else {
      // not captured, set State to 0
      CA.State[id] = 0;
    }
  }
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::slicapture(ZScell2D<T> &cells) {
  int id, locid_, locid1, locid2, id1, id2;
  id = cells.Id;
  // get position of neighbors with State = -1 (solid cells)
  Nbr_Pos.clear();
  getNbrs_pos(id, Nbr_Pos);
  bool captured = false;
  // create polygonal on 2 neighbors of solid cells
  for (int j = 0; j < Nbr_Pos.size(); j++) {
    // local id of solid cell: id_ = id + LatStru::Nbr[locid_]
    // get local index of solid cell
    locid_ = Nbr_Pos[j];
    // get local index of 2 neighbors of solid cell
    locid1 = (locid_ + 1) % 8;
    locid2 = (locid_ + 7) % 8;
    // get global index of 2 neighbors of solid cell
    id1 = id + LatStru::Nbr[locid1];
    id2 = id + LatStru::Nbr[locid2];
    // (create quadrilaterals): Cell[i] -> Vert[i] ->Vert[j] ->Cell[j]
    Polygon.clear();
    Polygon.emplace_back(Geo.getVoxel(id1)[0], Geo.getVoxel(id1)[1]);
    Polygon.emplace_back(getVertice(id1));
    Polygon.emplace_back(getVertice(id2));
    Polygon.emplace_back(Geo.getVoxel(id2)[0], Geo.getVoxel(id2)[1]);
    // check if centre of cell id is inside the polygonal
    Vector2D<T> point(Geo.getVoxel(id)[0], Geo.getVoxel(id)[1]);
    if (InPolygon(point, Polygon)) {
      captured = true;
      break;
    }
  }
  if (captured) {
    CA.Activate(id);
    cells.emplace_back(id);
    // if fluid field is calculated, set U to 0
    // Field.DeActivate(id);
  } else {
    // not captured, set State to 0
    CA.State[id] = 0;
  }
}

// template <typename T, typename LatStru>
// void ZhuStefanescu2D<T, LatStru>::MixedCapture() {
//   int id = 0;
//   int idn = 0;
//   // add to PreCaps and capture nearest neighbors
//   PreCaps.clear();
//   for (int i = 0; i < Cells.size(); i++) {
//     id = Cells[i].Id;
//     // get neighbors
//     for (int j = 0; j < LatStru::q; j++) {
//       idn = id + LatStru::Nbr[j];
//       if (CA.State[idn] == 0 && has_SolidNbr(idn)) {
//         // pesudo capture, if not captured, then set State to 0
//         if (has_SolidNbr(idn, 4)) {
//           CA.Activate(idn);
//           Cells.emplace_back(idn);
//         } else {
//           CA.State[idn] = 2;
//           PreCaps.emplace_back(idn);
//         }
//       }
//     }
//   }
//   // capture
//   // slicapture_s(PreCaps, Cells);
// }

// get K based on neighbor solid fraction gradient
template <typename T, typename LatStru>
inline T ZhuStefanescu2D<T, LatStru>::getK(int id) {
  using FDM = FDM2D<T>;
  // note that pow(FDM::p_x(id), 2) + pow(FDM::p_y(id), 2) may < 1e-6
  // causing K very large

  // T curv = pow(FDM::p_x(id), 2) + pow(FDM::p_y(id), 2);
  // curv = curv > 1e-3 ? pow(curv, -T(1.5)) : 0;
  // curv = pow(curv, -T(1.5));
  // T K_ = curv * (2 * FDM::p_x(id) * FDM::p_y(id) * FDM::p_xy(id) -
  //                pow(FDM::p_x(id), 2) * FDM::p_yy(id) -
  //                pow(FDM::p_y(id), 2) * FDM::p_xx(id));

  T K_ = pow(pow(FDM::p_x(id), 2) + pow(FDM::p_y(id), 2), -T(1.5)) *
         (2 * FDM::p_x(id) * FDM::p_y(id) * FDM::p_xy(id) -
          pow(FDM::p_x(id), 2) * FDM::p_yy(id) -
          pow(FDM::p_y(id), 2) * FDM::p_xx(id));
  K_ = K_ > T(2) ? T(2) : K_;
  K_ = K_ < T(-2) ? T(-2) : K_;
#ifdef _FLB_DEBUG
  CA.K[id] = K_;
#endif
  return K_;
}

// get K based on neighbor solid fraction
// K = 1 - 2/(n+1) * (fs + sum(fs_nbrs))
template <typename T, typename LatStru>
inline T ZhuStefanescu2D<T, LatStru>::getK_(int id) {
  int idn = 0;
  T sumfs = CA.f[id];
  for (int j = 0; j < LatStru::q; j++) {
    idn = id + LatStru::Nbr[j];
    sumfs += CA.f[idn];
  }
  return 1 - 2 / (LatStru::q + 1) * sumfs;
}

// get K based on neighbor solid fraction
// K = 1 - 2/(n+1) * (fs + sum(fs_nbrs))
template <typename T, typename LatStru>
inline void ZhuStefanescu2D<T, LatStru>::Calc_K_() {
  // traverse all interface cells
  int id = 0;
  int idn = 0;
#pragma omp parallel for private(id, idn) num_threads(Thread_Num)
  for (int i = 0; i < Cells.size(); i++) {
    id = Cells[i].Id;
    // get neighbors
    T sumfs = CA.f[id];
    for (int j = 0; j < LatStru::q; j++) {
      idn = id + LatStru::Nbr[j];
      sumfs += CA.f[idn];
    }
    Cells[i].K = 1 - 2 / (LatStru::q + 1) * sumfs;
  }
}

// get K based on neighbor solid fraction gradient
template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::Calc_K(T restrict_) {
  //  K = ((df/dx)^2 + (df/dy)^2)^(-3/2) *
  // (2*df/dx*df/dy*d^2f/dxdy - (df/dx)^2*d^2f/dy^2 - (df/dy)^2*d^2f/dx^2)
  // f(x,y) is the solid fraction
  // traverse all interface cells
  using FDM = FDM2D<T>;
  int id = 0;
#pragma omp parallel for private(id) num_threads(Thread_Num)
  for (int i = 0; i < Cells.size(); i++) {
    id = Cells[i].Id;
    // note that pow(FDM::p_x(id), 2) + pow(FDM::p_y(id), 2) may < 1e-6
    // causing K very large
    // check if pow(FDM::p_x(id), 2) + pow(FDM::p_y(id), 2) < 1e-6
    // if true, set K = 0
    T curv = pow(FDM::p_x(id), 2) + pow(FDM::p_y(id), 2);
    curv = curv > 1e-3 ? pow(curv, -T(1.5)) : 0;
    T K_ = curv * (2 * FDM::p_x(id) * FDM::p_y(id) * FDM::p_xy(id) -
                   pow(FDM::p_x(id), 2) * FDM::p_yy(id) -
                   pow(FDM::p_y(id), 2) * FDM::p_xx(id));
    if (K_ > T(restrict_)) K_ = T(restrict_);
    if (K_ < T(-restrict_)) K_ = T(-restrict_);
    // K[id] = K_;
    Cells[i].K = K_;
  }
}

template <typename T, typename LatStru>
inline T ZhuStefanescu2D<T, LatStru>::getPhi(int id) {
  // phi is growth angle to x-axis, which can be calculated by solid fraction:
  // phi = arctan(df/dy / df/dx), i.e., tan(phi) = df/dy / df/dx
  // atan() function returns in radians (-pi/2, pi/2)

  // more general form: cos(phi) = df/dx / sqrt((df/dx)^2 + (df/dy)^2)
  // phi = arccos(df/dx / sqrt((df/dx)^2 + (df/dy)^2))

  // T p_y_ = FDM2D<T>::p_y(id);
  // T p_x_ = FDM2D<T>::p_x(id);
  // T phi;
  // if (p_x_ == 0 && p_y_ == 0) {
  //   phi = 0;
  // } else {
  //   phi = atan(p_y_ / p_x_);
  // }
  // return phi;
  // return acos(FDM2D<T>::p_x(id) /
  //             sqrt(pow(FDM2D<T>::p_x(id), 2) + pow(FDM2D<T>::p_y(id), 2)));

  // atan(0/0) may result in undefined behavior.
  // use atan2() function, returns in radians (-pi, pi)
  return atan2(FDM2D<T>::p_y(id), FDM2D<T>::p_x(id));
}

template <typename T, typename LatStru>
inline T ZhuStefanescu2D<T, LatStru>::getanisotropy(int id) {
  // g(phi, theta) = 1 - delta * cos(4*(phi - theta))
  // where: delta is the anisotropy coefficient, manually chosen
  // phi is growth angle to x-axis, which can be calculated by solid fraction:
  // theta is the preferred growth angle to x-axis
  return 1 - delta * cos(4 * (getPhi(id) - Theta));
}

template <typename T, typename LatStru>
inline T ZhuStefanescu2D<T, LatStru>::getC_eq(ZScell2D<T> &cell) {
  // C_eq is the interface equilibrium composition, given by:
  // C_eq = C0 + [(Tl_eq-T_interface) - GT*K*g(theta)] / m  (m>0)
  // or   = C0 + [(T_interface-Tl_eq) + GT*K*g(theta)] / m  (m<0)
  // C0 is the initial composition,
  // Cl is the actual liquid composition
  // k is the partition coefficient
  // NOTE that m_l here is > 0
  int id = cell.Id;
  T Ceq =
      C0 +
      ((Tl_eq - lbmTH.getPoprho(id)) - GT * cell.K * getanisotropy(id)) / m_l;
  if (Ceq > 1 || Ceq < 0) {
    int x = id % Ni;
    int y = id / Ni;
    std::cout << "error: at (" << x << ", " << y << "), id = " << id
              << " ,Ceq = " << Ceq << ",\t"
              << "Tl_eq - Tl = " << Tl_eq - lbmTH.getPoprho(id) << ",\t"
              << "K = " << cell.K << ",\t"
              << "curv = "
              << pow(FDM2D<T>::p_x(id), 2) + pow(FDM2D<T>::p_y(id), 2) << ",\t"
              << "anisotropy = " << getanisotropy(id) << ",\t"
              << "GT * cell.K * getanisotropy(id) = "
              << GT * cell.K * getanisotropy(id) << std::endl;
    exit(-1);
  }
#ifdef _FLB_DEBUG
  CA.aniso[id] = GT * cell.K * getanisotropy(id);
  // CA.Ceq_Cl[id] = C_eq - getAveNbrPopRho(id);
#endif
  return Ceq;
}

template <typename T, typename LatStru>
void ZhuStefanescu2D<T, LatStru>::getdeltaf(ZScell2D<T> &cell) {
  // deltaf is the increased solid fraction during one time step, given by:
  // deltaf = (C_eq - Cl)/(C_eq*(1-k))
  int id = cell.Id;
  T C_eq = getC_eq(cell);
  T deltaf = (C_eq - lbmSO.getPoprho(id)) / (C_eq * _Part_Coef);
  // getAveNbrPopRho(id)
  //  T deltaf = (C_eq - getAveNbrPopRho(id)) / (C_eq * _Part_Coef);
  if (deltaf > 1) {
    int x = id % Ni;
    int y = id / Ni;
    std::cout << "error: at (" << x << ", " << y << "), id = " << id
              << " ,deltaf = " << deltaf << ",\t"
              << "Ceq = " << C_eq << ",\t"
              << "Cl = " << lbmSO.getPoprho(id) << std::endl;
    exit(-1);
  }
#ifdef _FLB_DEBUG
  CA.Ceq[id] = C_eq;
  CA.Ceq_Cl[id] = C_eq - lbmSO.getPoprho(id);
  // CA.Ceq_Cl[id] = C_eq - getAveNbrPopRho(id);
#endif
  deltaf = deltaf < 0 ? 0 : deltaf;
  cell.deltaf = deltaf;
}

template <typename T, typename LatStru>
inline T ZhuStefanescu2D<T, LatStru>::getLphi(T phi) {
  // Lphi = cellsize/max(|cos(phi)|, |sin(phi)|), where cellsize = 1
  return 1 / std::max(std::abs(cos(phi)), std::abs(sin(phi)));
}

template <typename T, typename LatStru>
inline Vector2D<T> ZhuStefanescu2D<T, LatStru>::getVertice(int id) {
  // Ip = Lphi * fs
  T Phi = getPhi(id);
  T Ip = getLphi(Phi) * CA.f[id];
  // get vertice (x[id]+Ipx, y[id]+Ipy)
  return Vector2D<T>(Geo.getVoxel(id)[0] + Ip * cos(Phi),
                     Geo.getVoxel(id)[1] + Ip * sin(Phi));
}

template <typename T, typename LatStru>
inline bool ZhuStefanescu2D<T, LatStru>::InPolygon(
    Vector2D<T> &point, std::vector<Vector2D<T>> &polygon) {
  bool flag = false;
  // vertices of polygon
  // i assigned to j, then i++, iterate all adjacent vertices
  for (int i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
    Vector2D<T> &p1 = polygon[i];
    Vector2D<T> &p2 = polygon[j];
    if (OnSegment(p1, p2, point)) return true;
    if (((sign(p1.y - point.y) > 0) != (sign(p2.y - point.y) > 0)) &&
        sign((point.x - p1.x) * (p1.y - p2.y) -
             (point.y - p1.y) * (p1.x - p2.x)) < 0)
      flag = !flag;
  }
  return flag;
}

template <typename T, typename LatStru>
inline void ZhuStefanescu2D<T, LatStru>::SortPolygon() {
  // clockwise sort vertices of polygon
}

template <typename T, typename LatStru>
inline bool ZhuStefanescu2D<T, LatStru>::OnSegment(Vector2D<T> &p1,
                                                   Vector2D<T> &p2,
                                                   Vector2D<T> &point) {
  return sign((p1 - point) ^ (p2 - point)) == 0 &&
         sign((p1 - point) * (p2 - point)) <= 0;
}

template <typename T, typename LatStru>
inline bool ZhuStefanescu2D<T, LatStru>::OnSegment_(Vector2D<T> &p1,
                                                    Vector2D<T> &p2,
                                                    Vector2D<T> &point) {
  // Check if point is collinear with p1 and p2
  if (point.x <= std::max(p1.x, p2.x) && point.x >= std::min(p1.x, p2.x) &&
      point.y <= std::max(p1.y, p2.y) && point.y >= std::min(p1.y, p2.y)) {
    // Check if point is on the line defined by p1 and p2
    if (std::abs((p2.y - p1.y) * (point.x - p1.x) -
                 (point.y - p1.y) * (p2.x - p1.x)) < 1e-6) {
      return true;
    }
  }
  return false;
}

template <typename T, typename LatStru>
inline int ZhuStefanescu2D<T, LatStru>::sign(T x) {
  if (std::abs(x) < 1e-6)
    return 0;
  else
    return (x > 0) ? 1 : -1;
}

template <typename T, typename LatStru>
inline void ZhuStefanescu2D<T, LatStru>::RejectSolutetoNbr(int id, T deltaf) {
  // calculate rejected solute
  // get neighbor cells with State = 0
  std::vector<int> nbrs;
  nbrs.reserve(LatStru::q);
  for (int j = 0; j < LatStru::q; j++) {
    int idn = id + LatStru::Nbr[j];
    if (CA.State[idn] == 0) {
      nbrs.emplace_back(idn);
    }
  }
  if (nbrs.size() == 0) return;
  T rejected_each = _Part_Coef * deltaf * lbmSO.getPoprho(id) / nbrs.size();
  // add rejected solute to neighbor cells
  for (int j = 0; j < nbrs.size(); j++) {
    // lbmSO.addtoRho(nbrs[j], rejected_each);
    lbmSO.addtoPop(nbrs[j], rejected_each);
  }
}

template <typename T, typename LatStru>
inline void ZhuStefanescu2D<T, LatStru>::RejectSolute(int id, T deltaf) {
  // calculate rejected solute
  // get neighbor cells with State = 0
  std::vector<int> nbrs;
  nbrs.reserve(LatStru::q + 1);
  nbrs.emplace_back(id);
  for (int j = 0; j < LatStru::q; j++) {
    int idn = id + LatStru::Nbr[j];
    if (CA.State[idn] == 0) {
      nbrs.emplace_back(idn);
    }
  }
  T rejected_each = _Part_Coef * deltaf * lbmSO.getPoprho(id) / nbrs.size();
  // add rejected solute to neighbor cells
  for (int j = 0; j < nbrs.size(); j++) {
    // lbmSO.addtoRho(nbrs[j], rejected_each);
    lbmSO.addtoPop(nbrs[j], rejected_each);
  }
}

template <typename T, typename LatStru>
inline T ZhuStefanescu2D<T, LatStru>::getAveNbrPopRho(int id) {
  // get average poprho of neighbors with State != -1
  T sum = lbmSO.getPoprho(id);
  int count = 1;
  for (int i = 0; i < LatStru::q; i++) {
    int idn = id + LatStru::Nbr[i];
    if (CA.State[idn] != -1) {
      sum += lbmSO.getPoprho(idn);
      count++;
    }
  }
  return sum / count;
}

template <typename T, typename LatStru>
inline T ZhuStefanescu2D<T, LatStru>::getStatisticalPopRho() {
  // get average poprho of all cells
  // for cells with State = 1 (interface cells)
  // Rho = C_Solid * fs + C_Liquid * (1 - fs)
  // for cells with State = 0 (liquid cells) or State = -1 (solid cells)
  // Rho = PopRho
  T sum = 0;
  int count = 0;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      int id = Index2D::GetId(i, j, Ni);
      if (CA.State[id] != 1) {
        sum += lbmSO.getPoprho(id);
      }
    }
  }
  // interface cells
  for (int i = 0; i < Cells.size(); i++) {
    if (CA.f[Cells[i].Id] > 0) {
      sum += Cells[i].C_Solid() * CA.f[Cells[i].Id] +
             lbmSO.getPoprho(Cells[i].Id) * (1 - CA.f[Cells[i].Id]);
    } else
      sum += lbmSO.getPoprho(Cells[i].Id);
  }
  return sum / T(Ni - 2) / T(Nj - 2);
}

template <typename T, typename LatStru>
template <int start>
inline void ZhuStefanescu2D<T, LatStru>::getNbrs_pos(int id,
                                                     std::vector<int> &nbrs) {
  // get position of neighbors which State = -1
  int idn = 0;
  for (int i = start; i < LatStru::q; i++) {
    idn = id + LatStru::Nbr[i];
    if (CA.State[idn] == -1) nbrs.emplace_back(i);
  }
}

// template <typename T, typename LatStru>
// void ZhuStefanescu2D<T, LatStru>::SimpleCapture() {
// #ifdef _OPENMP
//   if (Cells.size() > 20 * Thread_Num) {
//     OMP_helper::Divide_Index(IndexStart_OMP, IndexEnd_OMP, Cells.size());
//     int i = 0;
// #pragma omp parallel for private(i) num_threads(Thread_Num)
//     for (i = 0; i < IndexStart_OMP.size(); i++) {
//       SimpleCapture_omp(Captureds_OMP[i], IndexStart_OMP[i],
//       IndexEnd_OMP[i]);
//     }
//     // merge to Cells
//     for (int i = 0; i < Captureds_OMP.size(); i++) {
//       Cells.insert(Cells.end(), Captureds_OMP[i].begin(),
//                    Captureds_OMP[i].end());
//     }
//   } else {
//     SimpleCapture_s(Cells);
//   }
// #else
//   SimpleCapture_s(Cells);
// #endif
//   // add activated cells to interfaces
//   // Cells.insert(Cells.end(), Captureds.begin(), Captureds.end());
// }
// template <typename T, typename LatStru>
// void ZhuStefanescu2D<T, LatStru>::SimpleCapture_s(
//     std::vector<ZScell2D<T>> &captureds) {
//   int idn = 0;
//   int size = Cells.size();
//   for (int i = 0; i < size; i++) {
//     // get neighbors
//     for (int j = 0; j < LatStru::q; j++) {
//       idn = Cells[i].Id + LatStru::Nbr[j];
//       if (CA.State[idn] == 0 && has_SolidNbr(idn)) {
//         CA.Activate(idn);
//         captureds.emplace_back(idn);
//         // if fluid field is calculated, set U to 0
//         Field.DeActivate(idn);
//       }
//     }
//   }
// }
// template <typename T, typename LatStru>
// void ZhuStefanescu2D<T, LatStru>::SimpleCapture_omp(
//     std::vector<ZScell2D<T>> &captureds, int start, int end) {
//   captureds.clear();
//   int idn = 0;
//   for (int i = start; i < end; i++) {
//     // get neighbors
//     for (int j = 0; j < LatStru::q; j++) {
//       idn = Cells[i].Id + LatStru::Nbr[j];
//       if (CA.State[idn] == 0 && has_SolidNbr(idn)) {
// #ifdef _OPENMP
//         if (!_Atomic_Visited[idn].test_and_set(std::memory_order_acquire))
// #endif
//         {
//           CA.Activate(idn);
//           captureds.emplace_back(idn);
//           // if fluid field is calculated, set U to 0
//           Field.DeActivate(idn);
//         }
//       }
//     }
//   }
// }