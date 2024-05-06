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

// gandin.h
// this file presents the gandin's model for nucleation and growth

#pragma once

#include <atomic>
#include <random>
#include <vector>

#include "data_struct/Vector.h"
#include "lbm/lattice_set.h"
#include "legacy/legacy_lattice/cellidx.h"
#include "legacy/legacy_lattice/legacyfield.h"
#include "legacy/legacy_lbm/legacy_lbmanager2d.h"
#include "lbm/unit_converter.h"
#include "legacy/legacy_utils/legacy_message.h"
#include "utils/util.h"

// create pre-nuc sites and perform nucleation for ca: nuc2D
// growing(nucleated or captured) cells in CA mesh: Grow2D
// growing cells: Gcell2D

template <typename T>
class Gcell2D : public Vector<T, 2> {
 private:
  int _CellId;
  T _Arm;

 public:
  Gcell2D(int cellid, const Vector<T, 2>& vec, T arm = T(0.001))
      : Vector<T, 2>(vec), _CellId(cellid), _Arm(arm) {}
  Gcell2D(int cellid, T x, T y, T arm = T(0.001))
      : Vector<T, 2>(x, y), _CellId(cellid), _Arm(arm) {}

  int getCellId() const { return _CellId; }
  const T& getArm() const { return _Arm; }
  T& getArm() { return _Arm; }
};

template <typename T, typename LatStru>
class GandinCA2D : public LatStru, public message, public GetGCells2D {
 private:
  int N;
  // nucleation
  T Accum_Bulk_Sites;
  T Accum_Surf_Sites;
  // use lattice unit, need to be converted by UnitConverter first
  int quadrant[4][2] = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};  
  // 2D quadrant
  std::vector<Vector<int, 2>> Intercept2D = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
  std::vector<Vector<int, 2>> Vertex2D = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};

  // seed for random number generator
  // use random_device to generate seed, avoid same series of random numbers
  // call: rd() or rd = std::random_device()
  std::random_device _rd;
  // random number generator
  // call: _gen(rd()) or _gen = std::mt19937(rd())
  std::mt19937 _gen;

  std::vector<int> Interface;
  // growing cells
  std::vector<Gcell2D<T>> Cells;
  // newly solified cells
  std::vector<Gcell2D<T>> Solifieds;

#ifdef _OPENMP
  std::vector<int> IndexStart_OMP;
  std::vector<int> IndexEnd_OMP;
  std::vector<std::vector<Gcell2D<T>>> Nucleus_OMP;
  std::vector<std::vector<Gcell2D<T>>> Captureds_OMP;
  std::atomic_flag* _Atomic_Visited;
  void ResetVisited(int Num) {
#pragma omp parallel for
    for (int i = 0; i < Num; ++i) {
      _Atomic_Visited[i].clear();
    }
  }
#endif

  CAGField2D<T>& CA;
  Geometry2DLegacy<T>& Geo;
  GandinConverter<T>& ConvCA;
  LBManager2D<T>& LBM;

 public:
  GandinCA2D(CAGField2D<T>& ca, GandinConverter<T>& convca,
             LBManager2D<T>& lbm);
  ~GandinCA2D();
  //---------2023/10/11----------------
  // perform nucleation
  void Nucleation();
  void Nucleation_s(const std::vector<int>& RemainIdx,
                    std::normal_distribution<T>& GaussDistrib, int NucNum,
                    std::vector<Gcell2D<T>>& Nucs);
  void Nucleation_omp(const std::vector<int>& RemainIdx,
                      std::normal_distribution<T>& GaussDistrib, int NucNum,
                      std::vector<Gcell2D<T>>& Nucs);

  inline bool isNucleated(int Cellid, T dT_critic);
  void Single_Nucleation(int latx, int laty, T orine);

  // grow
  void Grow();
  inline T Get_GrowthVelo(int id);

  // capture
  void Capture();
  // newly captured cells will be added to Cells,
  // then perform post capture inside this function
  void Capture_s(std::vector<Gcell2D<T>>& captureds, int end, int start = 0);

  void postCapture();
  void postcapture_s(std::vector<Gcell2D<T>>& cells);

  // cell: growing cell, id: cell id to be captured or not
  void CellCapture(Gcell2D<T>& gcell, std::vector<Gcell2D<T>>& captureds,
                   int id);
  void CellCapture(Gcell2D<T>& cell, int id,
                   std::vector<Gcell2D<T>>& captureds);

  bool all_growing(int id);

  // get, transform assumes that the output range has the same size as the input
  void TransformInterface() {
    Interface.clear();
    // get cell id from Cells
    for (auto& cell : Cells) {
      Interface.emplace_back(cell.getCellId());
    }
  }
  std::vector<int>& getInterface() { return Interface; }
  std::vector<int>& Trans_GetInterface() {
    TransformInterface();
    return Interface;
  }

  // remember to call Get_remain cells before lbm
  void apply() {
    this->Get();
    Nucleation();
    Grow();
    Capture();
  }

  // remember to call Get_remain cells  before lbm
  void SingleNuc(int latx, int laty,
                 T orine = (std::rand() % 90) * M_PI / T(180)) {
    Single_Nucleation(latx, laty, orine);
  }

  void apply_SingleNuc() {
    this->Get();
    Grow();
    Capture();
  }
};