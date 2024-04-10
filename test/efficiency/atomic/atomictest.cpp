// atomictest.cpp
// this file test efficiency of std::atomic<bool> and std::atomic_flag
// efficiency of parrallel and serial operation will be tested, too

#include <omp.h>

#include <algorithm>
#include <atomic>
#include <iostream>
#include <vector>

#include "utils/timer.h"
#include "utils/util.h"

int Thread_Num = 16;

// generate a 2D mesh
class MESH {
 public:
  int _Ni;
  int _Nj;
  int _N;
  int* _Flag;
  int* flag;
  int* flag1;
  int* flag2;
  int _Nbr[8];
  std::vector<int> _Cells;

  MESH(int Ni, int Nj)
      : _Ni(Ni), _Nj(Nj), _N(Ni * Nj), _Nbr{1,      Ni,      -1,      -Ni,
                                            1 + Ni, -1 + Ni, -1 - Ni, 1 - Ni} {
    _Flag = new int[_N];
    flag = new int[_N];
    flag1 = new int[_N];
    flag2 = new int[_N];
    for (int i = 0; i < _N; i++) {
      _Flag[i] = 0;
      flag[i] = 0;
      flag1[i] = 0;
      flag2[i] = 0;
    }
    // set flags
    // centre bulk part set to 1
    Index2D::Traverse_Bulk(_Ni, _Nj, int(_Ni / 5),
                           [this](int id) { _Flag[id] = 1; });
    // peripheral part set to -1
    Index2D::Traverse_Peripheral(_Ni, _Nj, 0,
                                 [this](int id) { _Flag[id] = -1; });
    // set index
    _Cells.reserve(_N);
    Index2D::Traverse_Bulk(_Ni, _Nj, 1,
                           [this](int id) { _Cells.emplace_back(id); });
  }
  ~MESH() {
    delete[] _Flag;
    delete[] flag;
    delete[] flag1;
    delete[] flag2;
  }

  inline bool has_Nbr(int id) {
    for (int i = 0; i < 8; ++i) {
      if (_Flag[id + _Nbr[i]] == 1) return true;
    }
    return false;
  }
};

class ATOMIC {
 public:
  int _N;
  std::atomic<bool>* _Atomic_Visited;
  MESH& _Mesh;
  std::vector<int> IndexStart_OMP;
  std::vector<int> IndexEnd_OMP;
  std::vector<int> Cells;
  std::vector<std::vector<int>> Cells_OMP;

  ATOMIC(MESH& mesh) : _Mesh(mesh), _N(mesh._N) {
    _Atomic_Visited = new std::atomic<bool>[_N];
    Reset();
    IndexStart_OMP.reserve(Thread_Num);
    IndexEnd_OMP.reserve(Thread_Num);
    Cells.reserve(_N);
    OMP_helper::Setup_Thread_Vec(Cells_OMP, _N);
  }
  ~ATOMIC() { delete[] _Atomic_Visited; }
  void Reset() {
    // std::fill(_Atomic_Visited, _Atomic_Visited + _N, false);
    std::fill_n(_Atomic_Visited, _N, false);
  }
  int getCellSize() { return Cells.size(); }
  void get() {
    Reset();
    OMP_helper::Divide_Index(IndexStart_OMP, IndexEnd_OMP, _Mesh._Cells.size());
    int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < Cells_OMP.size(); ++i) {
      Get(Cells_OMP[i], IndexStart_OMP[i], IndexEnd_OMP[i]);
    }
    Cells.clear();
    // merge
    OMP_helper::Merge(Cells, Cells_OMP);
  }

  void Get(std::vector<int>& cells, int start, int end) {
    cells.clear();
    int id = 0;
    int idn = 0;
    for (int i = start; i < end; ++i) {
      id = _Mesh._Cells[i];
      if (_Mesh._Flag[id] == 1) {
        // get neighbours
        for (int j = 0; j < 8; j++) {
          idn = id + _Mesh._Nbr[j];
          if (_Mesh.has_Nbr(idn)) {
            // exchange set to true and return the old value
            // if flase, then emplace_back
            if (!_Atomic_Visited[idn].exchange(true,
                                               std::memory_order_acquire)) {
              cells.emplace_back(idn);
              _Mesh.flag[idn] = 1;
              _Mesh.flag1[idn] = 1;
              _Mesh.flag2[idn] = 1;
            }
          }
        }
      }
    }
  }
};

class ATOMICFLAG {
 public:
  int _N;
  std::atomic_flag* _Atomic_Visited;
  MESH& _Mesh;
  std::vector<int> IndexStart_OMP;
  std::vector<int> IndexEnd_OMP;
  std::vector<int> Cells;
  std::vector<std::vector<int>> Cells_OMP;

  ATOMICFLAG(MESH& mesh) : _Mesh(mesh), _N(mesh._N) {
    _Atomic_Visited = new std::atomic_flag[_N];
    Reset();
    IndexStart_OMP.reserve(Thread_Num);
    IndexEnd_OMP.reserve(Thread_Num);
    Cells.reserve(_N);
    OMP_helper::Setup_Thread_Vec(Cells_OMP, _N);
  }
  ~ATOMICFLAG() { delete[] _Atomic_Visited; }
  void Reset() {
#pragma omp parallel for num_threads(Thread_Num)
    for (int i = 0; i < _N; i++) {
      _Atomic_Visited[i].clear();
    }
  }
  int getCellSize() { return Cells.size(); }
  void get() {
    Reset();
    OMP_helper::Divide_Index(IndexStart_OMP, IndexEnd_OMP, _Mesh._Cells.size());
    int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
    for (i = 0; i < Cells_OMP.size(); ++i) {
      Get(Cells_OMP[i], IndexStart_OMP[i], IndexEnd_OMP[i]);
    }
    Cells.clear();
    // merge
    OMP_helper::Merge(Cells, Cells_OMP);
  }
  void Get(std::vector<int>& cells, int start, int end) {
    cells.clear();
    int id = 0;
    int idn = 0;
    for (int i = start; i < end; ++i) {
      id = _Mesh._Cells[i];
      if (_Mesh._Flag[id] == 1) {
        // get neighbours
        for (int j = 0; j < 8; j++) {
          idn = id + _Mesh._Nbr[j];
          if (_Mesh.has_Nbr(idn)) {
            if (!_Atomic_Visited[idn].test_and_set(std::memory_order_acquire)) {
              cells.emplace_back(idn);
              _Mesh.flag[idn] = 1;
              _Mesh.flag1[idn] = 1;
              _Mesh.flag2[idn] = 1;
            }
          }
        }
      }
    }
  }
};

class SERIAL {
 public:
  int _N;
  bool* _Visited;
  MESH& _Mesh;
  std::vector<int> Cells;

  SERIAL(MESH& mesh) : _Mesh(mesh), _N(mesh._N) {
    _Visited = new bool[_N];
    Reset();
    Cells.reserve(_N);
  }
  ~SERIAL() { delete[] _Visited; }
  void Reset() { std::fill_n(_Visited, _N, false); }
  int getCellSize() { return Cells.size(); }
  void Get() {
    Reset();
    Cells.clear();
    int id = 0;
    int idn = 0;
    for (int i = 0; i < _Mesh._Cells.size(); ++i) {
      id = _Mesh._Cells[i];
      if (_Mesh._Flag[id] == 1) {
        // get neighbours
        for (int j = 0; j < 8; j++) {
          idn = id + _Mesh._Nbr[j];
          if (_Mesh.has_Nbr(idn)) {
            if (!_Visited[idn]) {
              Cells.emplace_back(idn);
              _Visited[idn] = true;
              _Mesh.flag[idn] = 1;
              _Mesh.flag1[idn] = 1;
              _Mesh.flag2[idn] = 1;
            }
          }
        }
      }
    }
  }
};

int main() {
  // 1e4, 9e4, 10e4
  // Ni = Nj
  int Ni_List[4] = {102, 302, 1002, 3002};
  int repeat = 20;

  // MESH mesh1(Ni_List[0], Ni_List[0]);
  // MESH mesh2(Ni_List[1], Ni_List[1]);
  // MESH mesh3(Ni_List[2], Ni_List[2]);
  std::vector<MESH> meshes;
  meshes.reserve(4);
  meshes.emplace_back(Ni_List[0], Ni_List[0]);
  meshes.emplace_back(Ni_List[1], Ni_List[1]);
  meshes.emplace_back(Ni_List[2], Ni_List[2]);
  meshes.emplace_back(Ni_List[3], Ni_List[3]);

  // timer
  Timer ATOMIC_TIMER;
  Timer ATOMICFLAG_TIMER;
  Timer SERIAL_TIMER;

  std::vector<double> ATOMIC_Result;
  std::vector<double> ATOMICFLAG_Result;
  std::vector<double> SERIAL_Result;

  for (int i = 0; i < meshes.size(); ++i) {
    MESH& mesh = meshes[i];

    ATOMIC atomic(mesh);
    ATOMICFLAG atomicflag(mesh);
    SERIAL serial(mesh);

    // atomic
    ATOMIC_TIMER.START_TIMER();
    for (int j = 0; j < repeat; ++j) {
      atomic.get();
    }
    ATOMIC_TIMER.END_TIMER();
    ATOMIC_Result.push_back(ATOMIC_TIMER.GetTimeElapsed_Only() / repeat);

    // atomicflag
    ATOMICFLAG_TIMER.START_TIMER();
    for (int j = 0; j < repeat; ++j) {
      atomicflag.get();
    }
    ATOMICFLAG_TIMER.END_TIMER();
    ATOMICFLAG_Result.push_back(ATOMICFLAG_TIMER.GetTimeElapsed_Only() /
                                repeat);

    // serial
    SERIAL_TIMER.START_TIMER();
    for (int j = 0; j < repeat; ++j) {
      serial.Get();
    }
    SERIAL_TIMER.END_TIMER();
    SERIAL_Result.push_back(SERIAL_TIMER.GetTimeElapsed_Only() / repeat);

    // size:
    std::cout << "ATOMIC: " << atomic.getCellSize() << "\n";
    std::cout << "ATOMICFLAG: " << atomicflag.getCellSize() << "\n";
    std::cout << "SERIAL: " << serial.getCellSize() << "\n";
  }

  // print time elapsed
  std::cout << "ATOMIC: "
            << "\n";
  for (int i = 0; i < ATOMIC_Result.size(); ++i) {
    std::cout << "N: " << meshes[i]._Cells.size()
              << "  Time: " << ATOMIC_Result[i] << "\n";
  }
  std::cout << "ATOMICFLAG: "
            << "\n";
  for (int i = 0; i < ATOMICFLAG_Result.size(); ++i) {
    std::cout << "N: " << meshes[i]._Cells.size()
              << "  Time: " << ATOMICFLAG_Result[i] << "\n";
  }
  std::cout << "SERIAL: "
            << "\n";
  for (int i = 0; i < SERIAL_Result.size(); ++i) {
    std::cout << "N: " << meshes[i]._Cells.size()
              << "  Time: " << SERIAL_Result[i] << "\n";
  }

  return 0;
}