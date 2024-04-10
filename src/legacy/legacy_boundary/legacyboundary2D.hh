// boundary2D
#include "legacy/legacy_boundary/legacyboundary2D.h"

template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::Setup_Cell_Dirs(int *flags) {
  int i = 0;
#pragma omp parallel for private(i) collapse(1) num_threads(Thread_Num)
  for (i = 0; i < GetBdCellNum(); i++) {
    SetCell_Dirs(i, flags);
  }
}

template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::Setup_Cell_Dirs(Geometry2DLegacy<T> &Geo) {
  int i = 0;
#pragma omp parallel for private(i) collapse(1) num_threads(Thread_Num)
  for (i = 0; i < GetBdCellNum(); i++) {
    SetCell_Dirs(i, Geo);
  }
}

template <typename T, typename LatSet, typename LatStru>
inline void BasicBoundary<T, LatSet, LatStru>::Setup_Cell_Dirs_(int *flags) {
  int i, num = GetBdCellNum();
#pragma omp parallel for private(i) collapse(1) num_threads(Thread_Num)
  for (i = 0; i < num; i++) {
    Cell_Dirs[i].template set<LatSet, LatStru>(flags, Lattice.Nbr);
  }
}

template <typename T, typename LatSet, typename LatStru>
inline void BasicBoundary<T, LatSet, LatStru>::SetCell_Dirs(int i, int *flags) {
  Cell_Dirs[i].Clear();
  // cell index
  int Id = Cell_Dirs[i].Id;
  for (int k = 1; k < LatSet::q; k++) {
    if (flags[Id + Lattice.Nbr[k]] != -1) {
      Cell_Dirs[i].add_inflow(LatSet::opp[k]);
    }
    if (flags[Id + Lattice.Nbr[k]] == -1 &&
        flags[Id + Lattice.Nbr[LatSet::opp[k]]] != -1) {
      Cell_Dirs[i].add_outflow(LatSet::opp[k]);
    }
  }
}

template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::SetCell_Dirs(int i,
                                                     Geometry2DLegacy<T> &Geo) {
  Cell_Dirs[i].Clear();
  // cell index
  int Id = Cell_Dirs[i].Id;
  for (int k = 1; k < LatSet::q; k++) {
    if (Geo.getVoxel(Id + Lattice.Nbr[k]).getFlag() != -1) {
      Cell_Dirs[i].add_inflow(LatSet::opp[k]);
    }
    if (Geo.getVoxel(Id + Lattice.Nbr[k]).getFlag() == -1 &&
        Geo.getVoxel(Id + Lattice.Nbr[LatSet::opp[k]]).getFlag() != -1) {
      Cell_Dirs[i].add_outflow(LatSet::opp[k]);
    }
  }
}

template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::SetFromGeo(int flag,
                                                   Geometry2DLegacy<T> &geo) {
  for (int i = 0; i < geo.getVoxels().size(); ++i) {
    if (geo.getVoxel(i).getFlag() == flag) {
      Add(i);
    }
  }
  Setup_Cell_Dirs(geo);
}
template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::SetFromGeo(int flag, Geometry2DLegacy<T> &geo,
                                                   T latrho) {
  for (int i = 0; i < geo.getVoxels().size(); ++i) {
    if (geo.getVoxel(i).getFlag() == flag) {
      Add(i, latrho);
    }
  }
  Setup_Cell_Dirs(geo);
}
template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::SetFromFlag(int flag, int *flags,
                                                    int ni, int nj) {
  for (int j = 1; j < nj - 1; j++) {
    for (int i = 1; i < ni - 1; i++) {
      int Id = i + j * ni;
      if (flags[Id] == flag) {
        Add(Id);
      }
    }
  }
  Setup_Cell_Dirs(flags);
}
template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::SetFromFlag(int flag, int *flags,
                                                    int ni, int nj, T latrho) {
  for (int j = 1; j < nj - 1; j++) {
    for (int i = 1; i < ni - 1; i++) {
      int Id = i + j * ni;
      if (flags[Id] == flag) {
        Add(Id, latrho);
      }
    }
  }
  Setup_Cell_Dirs(flags);
}

template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::SetFromVector(
    std::vector<int> &cell_ids, int *flags) {
  for (int i = 0; i < cell_ids.size(); i++) {
    Add(cell_ids[i]);
  }
  Setup_Cell_Dirs(flags);
}
template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::SetFromVector(
    std::vector<int> &cell_ids, int *flags, T latrho) {
  for (int i = 0; i < cell_ids.size(); i++) {
    Add(cell_ids[i], latrho);
  }
  Setup_Cell_Dirs(flags);
}
template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::SetFromVector(
    std::vector<int> &cell_ids, int *flags, std::vector<T> &latrhos) {
  if (cell_ids.size() != latrhos.size()) {
    std::cout << "unmatched size: cell_ids and latrhos" << std::endl;
    exit(-1);
  }
  for (int i = 0; i < cell_ids.size(); i++) {
    Add(cell_ids[i], latrhos[i]);
  }
  Setup_Cell_Dirs(flags);
}

template <typename T, typename LatSet, typename LatStru>
void BasicBoundary<T, LatSet, LatStru>::Stream(population<T, LatSet> *pop) {
  int Id = 0;
  int Dir = 0;
#pragma omp parallel for private(Id, Dir) num_threads(Thread_Num)
  for (int i = 0; i < GetBdCellNum(); i++) {
    Id = Cell_Dirs[i].Id;
    pop[Id].f[0] = pop[Id].fpostcol[0];  // IMPROTANT!!!
    for (int k = 0; k < Cell_Dirs[i].inflow.size(); k++) {
      Dir = Cell_Dirs[i].inflow[k];
      pop[Id].f[Dir] = pop[Id + Lattice.Nbr[LatSet::opp[Dir]]].fpostcol[Dir];
    }
  }
}

template <typename T, typename LatSet, typename LatStru>
inline void BasicBoundary<T, LatSet, LatStru>::statistics() {
  std::cout << name << "\t" << GetBdCellNum() << "\n";
}
