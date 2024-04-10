#pragma once

#include "legacy/legacy_lattice/cellidx.h"

IndexManager2D::IndexManager2D(int ni, int nj, std::vector<int>& Idx_)
    : Ni(ni), Nj(nj), N(ni * nj) {
  Idx.reserve((Ni - 2) * (Nj - 2));
  // copy cells to Cells
  Idx.insert(Idx.end(), Idx_.begin(), Idx_.end());
}

template <typename LatSet, typename LatStru>
CellIndexManager2D<LatSet, LatStru>::CellIndexManager2D(
    int ni, int nj, int* flag, std::vector<direction<LatStru::q>>& bound)
    : Ni(ni),
      Nj(nj),
      N(ni * nj),
      Flag(flag),
      LatStru(ni, nj),
      BBcelldirs(bound) {
  Cells.reserve((Ni - 2) * (Nj - 2));
  InCells.reserve((Ni - 4) * (Nj - 4));
  BdCells.reserve((Ni + Nj) * 4);

  isbound_ = new bool[N];
  std::fill_n(isbound_, N, false);
  for (int i = 0; i < bound.size(); i++) {
    isbound_[bound[i].Id] = true;
  }
  // initialise:
  Setup();
}

// low efficiency, can be used to initialise
template <typename LatSet, typename LatStru>
inline void CellIndexManager2D<LatSet, LatStru>::Setup() {
  // Cells
  Cells.clear();
  Traverse_Bulk(Ni, Nj, 1, [this](int id) {
    if (Flag[id] != -1) Cells.emplace_back(id);
  });
  // InCells and BdCells
  InCells.clear();
  BdCells.clear();
  Traverse_Vector<int>(Cells, [this](int id) {
    if (isbound(id))
      BdCells.emplace_back(id);
    else
      InCells.emplace_back(id);
  });
}

template <typename LatSet, typename LatStru>
void CellIndexManager2D<LatSet, LatStru>::Get(std::vector<int>& interface) {
  // Cells
  Cells.clear();
  Traverse_Bulk(Ni, Nj, 1, [this](int id) {
    if (Flag[id] != -1) Cells.emplace_back(id);
  });
  // erase inactivate cells from BBcelldirs
  auto iter = BBcelldirs.begin();
  while (iter != BBcelldirs.end()) {
    if (Flag[iter->Id] == -1) {
      iter = BBcelldirs.erase(iter);
      isbound_[iter->Id] = false;
    } else
      iter++;
  }
  // add to BBcelldirs
  Traverse_Vector<int>(interface, [this](int id) {
    int idn = id;
    for (int i = 0; i < LatStru::q; i++) {
      idn = id + LatStru::Nbr[i];
      if (Flag[idn] != -1 && isbound(idn) && !isbound_[idn]) {
        BBcelldirs.emplace_back(idn);
        isbound_[idn] = true;
        // set BBcelldirs after adding
        // DO NOT set inside the loop
        // BBcelldirs.back().template set<LatSet, LatStru>(Flag, LatStru::Nbr);
      }
    }
  });
  // set celldirs
  int i = 0;
  int num = BBcelldirs.size();
#pragma omp parallel for private(i) collapse(1) num_threads(Thread_Num)
  for (i = 0; i < num; i++) {
    BBcelldirs[i].template set<LatSet, LatStru>(Flag, LatStru::Nbr);
  }

  // InCells
  InCells.clear();
  Traverse_Vector<int>(Cells, [this](int id) {
    if (!isbound_[id]) InCells.emplace_back(id);
  });
}

template <typename LatSet, typename LatStru>
inline void CellIndexManager2D<LatSet, LatStru>::Erase() {
  // Cells
  if (Cells.size() > 1000 * Thread_Num) {
    // erase Cells in parallel
    std::vector<std::vector<int>> Cells_Th;
    OMP_helper::Divide_and_Setup<int>(Cells, Cells_Th);
    int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < Cells_Th.size(); i++) {
      erase_s(Cells_Th[i]);
    }
    // merge
    Cells.clear();
    OMP_helper::Merge<int>(Cells, Cells_Th);
  } else {
    // erase Cells serially
    erase_s(Cells);
  }
  // InCells
  if (InCells.size() > 1e5) {
    // erase InCells in parallel
    std::vector<std::vector<int>> InCells_Th;
    OMP_helper::Divide_and_Setup<int>(InCells, InCells_Th);
    int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < InCells_Th.size(); i++) {
      erase_sin(InCells_Th[i]);
    }
    // merge
    InCells.clear();
    OMP_helper::Merge<int>(InCells, InCells_Th);
  } else {
    // erase InCells serially
    erase_sin(InCells);
  }
}

template <typename LatSet, typename LatStru>
inline void CellIndexManager2D<LatSet, LatStru>::erase_s(
    std::vector<int>& cells) {
  // serially erase Cells
  std::vector<int>::iterator iter = cells.begin();
  while (iter != cells.end()) {
    if (Flag[*iter] == -1)
      iter = cells.erase(iter);
    else
      iter++;
  }
}
template <typename LatSet, typename LatStru>
inline void CellIndexManager2D<LatSet, LatStru>::erase_sin(
    std::vector<int>& incells) {
  // serially erase InCells
  std::vector<int>::iterator iter = incells.begin();
  while (iter != incells.end()) {
    if (Flag[*iter] != 0)
      iter = incells.erase(iter);
    else
      iter++;
  }
}

template <typename LatSet, typename LatStru>
inline bool CellIndexManager2D<LatSet, LatStru>::isbound(int id) {
  for (int i = 1; i < LatStru::q; i++) {
    if (Flag[id + LatStru::Nbr[i]] == -1) return true;
  }
  return false;
}
// template <typename LatSet, typename LatStru>
// inline bool CellIndexManager2D<LatSet, LatStru>::isbound_simd(int id) {
//   // OR(||) operation: the final value of isbound will be true if any of the
//   // loop iterations set it to true
//   bool isbound = false;
// #pragma omp simd reduction(|| : isbound)
//   for (int i = 1; i < LatStru::q; i++) {
//     if (Flag[Id + LatStru::Nbr[i]] == -1) isbound = true;
//   }
//   return isbound;
// }

//////////////////////////////////////
GetGCells2D::GetGCells2D(int ni, int nj, int* flag, int get_method)
    : Ni(ni), Nj(nj), N(ni * nj), Flag(flag), Get_Method(get_method) {
  int max_bulk = (ni - 4) * (nj - 4);
  int max_surf = (ni + nj - 4) * 2 - 4;

  Bulks.reserve(max_bulk);
  Surfs.reserve(max_surf);
  // Cells.reserve(max_bulk + max_surf);
  if (get_method == 1) {
    OMP_helper::Setup_Thread_Vec<int>(Bulks_Thread, max_bulk);
  }
  // initialise:
  Get_pushback();
}

// get remain cells using push_back
void GetGCells2D::Get_pushback() {
  // clear
  Bulks.clear();
  Surfs.clear();

  // using serialised push_back
  Traverse_Bulk(Ni, Nj, 2, [this](int id) {
    if (Flag[id] != -1) Bulks.emplace_back(id);
  });
  Traverse_Peripheral(Ni, Nj, 1, [this](int id) {
    if (Flag[id] != -1) Surfs.emplace_back(id);
  });
}

// get remain cells using erase
void GetGCells2D::Get_erase() {
#ifdef _OPENMP
  if (Bulks.size() > 10000) {  // divide bulks and allocate to threads
    OMP_helper::Divide<int>(Bulks, Bulks_Thread);
    // parallel erase bulk cells
    int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < Bulks_Thread.size(); i++) {
      Get_erase_s(Bulks_Thread[i]);
    }
    // merge
    Bulks.clear();
    OMP_helper::Merge_and_Clear<int>(Bulks, Bulks_Thread);
  } else {
    // serialised erase bulk cells
    Get_erase_s(Bulks);
  }
#else
  // serialised erase bulk cells
  Get_erase_s(Bulks);
#endif
  // serialised erase surf cells
  Get_erase_s(Surfs);
}

inline void GetGCells2D::Get_erase_s(std::vector<int>& cells) {
  std::vector<int>::iterator iter = cells.begin();
  while (iter != cells.end()) {
    if (Flag[*iter] == -1)
      iter = cells.erase(iter);
    else
      iter++;
  }
}
