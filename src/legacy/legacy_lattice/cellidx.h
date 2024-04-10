// cellidx.h
//
#pragma once

#include "legacy/legacy_boundary/legacyboundary2D.h"
#include "utils/util.h"

class IndexManager2D : public Index2D {
 private:
  int Ni;
  int Nj;
  int N;
  std::vector<int> Idx;

 public:
  IndexManager2D(int ni, int nj, std::vector<int>& Idx_);
  // get
  std::vector<int>& Get() { return Idx; }
  int Get(int i, int j) { return Idx[i + j * Ni]; }
  int Get(int id) { return Idx[id]; }
};

template <typename LatSet, typename LatStru>
class CellIndexManager2D : public Index2D, public LatStru {
 private:
  int Ni;
  int Nj;
  int N;
  // CAFIELD flag(State)
  int* Flag;
  // cells with flag != -1
  std::vector<int> Cells;
  // inner cells with flag == 0
  std::vector<int> InCells;
  // boundary cells, Cells = InCells + BdCells
  std::vector<int> BdCells;
  // refer to celldirs in BB BCs
  std::vector<direction<LatStru::q>>& BBcelldirs;

  // boundary flag
  bool* isbound_ = nullptr;

 public:
  // moving boundaries
  CellIndexManager2D(int ni, int nj, int* flag,
                     std::vector<direction<LatStru::q>>& bound);

  ~CellIndexManager2D() { delete[] isbound_; }

  // get using emplace_back
  // need to set Cell_Dirs in BB after communication
  inline void Setup();
  // call: (CA.GetInterface())
  // set Cell_Dirs inside CellCommunicator
  void Get(std::vector<int>& interface);
  // get using erase
  inline void Erase();
  inline void erase_s(std::vector<int>& cells);
  inline void erase_sin(std::vector<int>& incells);

  // BCs
  inline bool isbound(int id);

  // get
  std::vector<int>& Get_Cells() { return Cells; }
  std::vector<int>& Get_InCells() { return InCells; }
  std::vector<int>& Get_BdCells() { return BdCells; }
  bool* Get_isbound() { return isbound_; }
  // getSolidFraction of the whole domain
  template <typename T = double>
  T getSolidFraction() {
    int count = 0;
    Traverse_Bulk(Ni, Nj, 1, [this, &count](int id) {
      if (Flag[id] == -1) count++;
    });
    return T(count) / T((Ni - 2) * (Nj - 2));
  }
};

// utils
void AddtoIdx(std::vector<int>& idx, int* flags, int flag, int size) {
  for (int i = 0; i < size; i++) {
    if (flags[i] == flag) idx.emplace_back(i);
  }
}
void AddtoIdx(std::vector<int>& idx, std::vector<int>& vec) {
  idx.insert(idx.end(), vec.begin(), vec.end());
}

// make sure Geo is fully set before instantiating this class
// Cells will be set from Geo if (Flag[id] != -1)
class GetGCells2D : public Index2D {
 private:
  int Ni;
  int Nj;
  int N;
  // CAFIELD flag(State)
  int* Flag;

  int Get_Method;

  // store remain cells
  std::vector<int> Bulks;
  std::vector<int> Surfs;
  // merge Remain_Bulk and Remain_Surf
  // std::vector<int> Cells;
  // parallel erase bulk cells，enabled by default (get_method = 1)
  std::vector<std::vector<int>> Bulks_Thread;

 public:
  // get_method: 0: push_back, 1: erase
  GetGCells2D(int ni, int nj, int* flag, int get_method = 0);
  // get remain cells using push_back
  // clear Bulks and Surfs first, then push_back respectively
  // finally merge to Cells(cleared before merge)
  void Get_pushback();
  // get remain cells using erase
  // divide bulks and allocate to threads, perform erase in parallel
  // then merge to Bulks(cleared before merge)
  // serialised erase surf cells(for num of surf cells is small)
  // finally merge to Cells(cleared before merge)
  void Get_erase();
  void Get_erase_s(std::vector<int>& cells);
  // get remain cells
  void Get() {
    if (Get_Method == 0)
      Get_pushback();
    else
      Get_erase();
  }
  // get
  std::vector<int>& Get_Bulks() { return Bulks; }
  std::vector<int>& Get_Surfs() { return Surfs; }
  // std::vector<int>& Get_Cells() { return Cells; }
};
