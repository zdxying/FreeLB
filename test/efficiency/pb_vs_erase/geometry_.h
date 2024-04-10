// geometry class
// legacy 
#pragma once
#include "utils/util.h"

template <typename T>
struct Geometry2D_legacy : public Index2D {
  /* mesh data, use lattice unit */
  T *x;
  T *y;
  int *Flag;  // geometry flag; -1 means inactive
  /*mesh param*/
  int Ni;
  int Nj;
  int N;

  Geometry2D_legacy(int Ni_, int Nj_);
  ~Geometry2D_legacy();                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
   
  void Setup();
  void Print();
  void PrintActiveCellNum();

  // Set geometry
  void get_rectangle(T x0_, T x1_, T y0_, T y1_, int &i0, int &i1, int &j0,
                     int &j1);
  // set flag = -1, unable to handle other flags
  void circle(T radius, T centre_x, T centre_y, int flag = -1);
  void check_circle(T radius, T centre_x, T centre_y);
  // set flag = -1, unable to handle other flags
  void rectangle(T x0, T x1, T y0, T y1, int flag = -1);
  void check_rectangle(T x0, T x1, T y0, T y1);
  // set flag = -1, unable to handle other flags
  void triangle(T x0, T x1, T x2, T y0, T y1, T y2, int flag = -1);
  void check_triangle(T x0, T x1, T x2, T y0, T y1, T y2);

  // post process, need typename LatSet and LatStru
  // declared in lbm2D.h

  // get/add indices of cells with flag using idx.emplace_back()
  // no other operations on idx
  void get_Idx(std::vector<int> &idx, int flag);
  // add to indices of cells with flag
  // void add_Idx(std::vector<int> &idx, int flag);

  // quick set
  // set boundary cell flag to
  void set_bottom(int flag = 1)  // j = 1
  {
    for (int i = 1; i < Ni - 1; i++) {
      Flag[GetId(i, 1, Ni)] = flag;
    }
  }
  void set_top(int flag = 1) {
    for (int i = 1; i < Ni - 1; i++) {
      Flag[GetId(i, Nj - 2, Ni)] = flag;
    }
  }
  void set_left(int flag = 1) {
    for (int j = 1; j < Nj - 1; j++) {
      Flag[GetId(1, j, Ni)] = flag;
    }
  }
  void set_right(int flag = 1) {
    for (int j = 1; j < Nj - 1; j++) {
      Flag[GetId(Ni - 2, j, Ni)] = flag;
    }
  }
  void quickset_boundary(int flag = 1) {
    set_bottom(flag);
    set_top(flag);
    set_left(flag);
    set_right(flag);
  }
};

template <typename T>
struct Geometry3D : public Index3D {
  /* mesh data */
  T *x;
  T *y;
  T *z;
  int *Flag;
  /*mesh param*/
  int Ni;
  int Nj;
  int Nk;
  int NiNj;
  int N;

  Geometry3D(int Ni_, int Nj_, int Nk_) : Ni(Ni_), Nj(Nj_), Nk(Nk_) {
    N = Ni * Nj * Nk;
    NiNj = Ni * Nj;
    x = new T[N];
    y = new T[N];
    z = new T[N];
    Flag = new int[N];  // 0: fluid, -1: solid, 1: mush zone
    Setup();
    Print();
  }
  ~Geometry3D() {
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] Flag;
  }
  void Setup() {
    // init mesh
    for (int k = 0; k < Nk; k++) {
      for (int j = 0; j < Nj; j++) {
        for (int i = 0; i < Ni; i++) {
          int id = GetId(i, j, k, Ni, NiNj);
          // use lattice unit: cell length = 1.0f
          x[id] = i * T(1);
          y[id] = j * T(1);
          z[id] = k * T(1);
          Flag[id] = 0;
        }
      }
    }
    // periphery cells assumed to be solid, convenient for boundary neighbor
    // search k = 0, k = Nk - 1
    for (int j = 0; j < Nj; j++) {
      for (int i = 0; i < Ni; i++) {
        Flag[GetId(i, j, 0, Ni, NiNj)] = -1;
        Flag[GetId(i, j, Nk - 1, Ni, NiNj)] = -1;
      }
    }
    // j = 0, j = Nj -1
    for (int k = 0; k < Nk; k++) {
      for (int i = 0; i < Ni; i++) {
        Flag[GetId(i, 0, k, Ni, NiNj)] = -1;
        Flag[GetId(i, Nj - 1, k, Ni, NiNj)] = -1;
      }
    }
    // i = 0, i = Ni - 1
    for (int k = 0; k < Nk; k++) {
      for (int j = 0; j < Nj; j++) {
        Flag[GetId(0, j, k, Ni, NiNj)] = -1;
        Flag[GetId(Ni - 1, j, k, Ni, NiNj)] = -1;
      }
    }
  }
  void Print() {
    std::cout << "[Mesh3D info]:\n"
              << "Ni:       " << Ni << "\n"
              << "Nj:       " << Nj << "\n"
              << "Nk:       " << Nk << "\n"
              << "N:        " << N << "\n"
              << std::endl;
  }
  void PrintActiveCellNum() {
    int count = 0;
    for (int k = 1; k < Nk - 1; k++) {
      for (int j = 1; j < Nj - 1; j++) {
        for (int i = 1; i < Ni - 1; i++) {
          int id = GetId(i, j, k, Ni, NiNj);
          if (Flag[id] != -1) count++;
        }
      }
    }
    std::cout << "ComputeCells: " << count << std::endl;
  }

  void circle(T radius, T centre_x, T centre_y, T centre_z, int flag) {
    T x_, y_, z_;
    for (int k = 1; k < Nk - 1; k++) {
      for (int j = 1; j < Nj - 1; j++) {
        for (int i = 1; i < Ni - 1; i++) {
          int id = GetId(i, j, k, Ni, NiNj);
          x_ = x[id];
          y_ = y[id];
          z_ = z[id];
          if (((x_ - centre_x) * (x_ - centre_x) +
               (y_ - centre_y) * (y_ - centre_y) +
               (z_ - centre_z) * (z_ - centre_z)) < radius * radius) {
            Flag[id] = flag;
          }
        }
      }
    }
  }
  void cube(T x0, T x1, T y0, T y1, T z0, T z1, int flag) {
    T x_, y_, z_;
    for (int k = 0; k < Nk; k++) {
      for (int j = 1; j < Nj - 1; j++) {
        for (int i = 1; i < Ni - 1; i++) {
          int id = GetId(i, j, k, Ni, NiNj);
          x_ = x[id];
          y_ = y[id];
          z_ = z[id];
          if (x_ > x0 && x_ < x1 && y_ > y0 && y_ < y1 && z_ > z0 && z_ < z1) {
            Flag[id] = flag;
          }
        }
      }
    }
  }
};
