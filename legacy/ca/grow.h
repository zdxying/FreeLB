// grow.h
#pragma once

#include <atomic>
#include <vector>

#include "lbm/lattice_set.h"
#include "legacy/legacy_lattice/legacyfield.h"
#include "legacy/legacy_lbm/legacy_lbmanager2d.h"
#include "lbm/unit_converter.h"
#include "utils/util.h"

// info about single growing(nucleated or captured) cell in CA
template <typename T>
struct gcell2D {
  int Id;  // corresponding Cell Id
  T x;     // growing centre position, may not be cell centre
  T y;     // growing centre position, may not be cell centre
  T Arm;   // length of arm in 4 directions
  // set to a small value for initialize mush zone
  // T f;  // solid fraction
  gcell2D(int id, T x, T y) : Id(id), x(x), y(y), Arm(T(0.001)) {}
  gcell2D(int id, T x, T y, T arm) : Id(id), x(x), y(y), Arm(arm) {}
};

// store info about growing(nucleated or captured) cells in CA mesh
template <typename T, typename LatStru>
class Grow2D : public LatStru {
 private:
  int Ni;
  int Nj;
  int quadrant[4][2] = {{1, 1}, {1, -1}, {-1, -1}, {-1, 1}};  // 2D quadrant
  T Lat_GrowthPara;
  // growing cells
  std::vector<gcell2D<T>> Cells;
  // newly nucleated cells
  std::vector<gcell2D<T>> New_Nucs;
  // newly captured cells
  std::vector<gcell2D<T>> Captureds;
  // new growings for cell communication
  std::vector<int> New_Growings;
  // new solidified cells for source term calculation
  std::vector<int> New_Solidifieds;

#ifdef _CAPTURE_P
  std::atomic_flag* visited_ato;
  std::vector<std::vector<gcell2D<T>>> Captureds_Th;
  std::vector<std::vector<gcell2D<T>>> Cells_Th;
#endif

  CAGField2D<T>& CA;
  Geometry2DLegacy<T>& Geo;
  GandinConverter<T>& ConvCA;
  LBManager2D<T>& LBM;

 public:
  Grow2D(CAGField2D<T>& ca, GandinConverter<T>& convca, LBManager2D<T>& lbm);
  ~Grow2D();

  // std::vector<gcell2D<T>>& cells = x.getCells();
  std::vector<gcell2D<T>>& getCells() { return Cells; }
  std::vector<gcell2D<T>>& getNewNucs() { return New_Nucs; }
  std::vector<int>& getNewGrowings() { return New_Growings; }

  // newly nucleated cells will be added to Cells
  // after Grow step, inside this function
  // i.e., newly nucleated cells will grow until next time step
  void Grow();

  T Get_GrowthVelo(int id);

  bool all_growing(int id);
  // newly captured cells will be added to Cells,
  // then perform post capture inside this function
  void Capture();

  void postCapture();
  void postcapture_s(std::vector<gcell2D<T>>& cells);

  // cell: growing cell, id: cell id to be captured or not
  void CellCapture(gcell2D<T>& gcell, int id);

  void addto_NewGrowings() {
    New_Growings.clear();
    // insert ids from New_Nucs and Captureds to New_Growings
    int id = 0;
    for (int i = 0; i < New_Nucs.size(); i++) {
      id = New_Nucs[i].Id;
      New_Growings.emplace_back(id);
    }
    for (int i = 0; i < Captureds.size(); i++) {
      id = Captureds[i].Id;
      New_Growings.emplace_back(id);
    }
  }
};