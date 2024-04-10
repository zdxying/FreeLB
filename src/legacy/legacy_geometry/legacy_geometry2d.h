
// legacy_geometry2d.h

#pragma once

#include <vector>
#include "data_struct/Vector.h"

////////////////////////////////////////////////////////////////////////
// legacy code
////////////////////////////////////////////////////////////////////////

// 22.10.2023: use Vector<T, 2> instead of x, y
// this is a structured mesh

// this is the old version of Geometry2D
template <typename T>
class Geometry2DLegacy {
 private:
  // Geometry2D data: T _data[2]; int flag_;
  std::vector<Voxel<T, 2>> voxels_;
  // Geometry2D info
  int Ni_;
  int Nj_;
  int N_;

 public:
  Geometry2DLegacy(int Ni, int Nj);

  int getNi() const { return Ni_; }
  int getNj() const { return Nj_; }
  int getN() const { return N_; }
  // get position(Vector<T, 2>)
  Vector<T, 2>& getLoc(int id) { return voxels_[id].getLoc(); }
  // get voxel
  Voxel<T, 2>& getVoxel(int id) { return voxels_[id]; }
  // get voxels
  std::vector<Voxel<T, 2>>& getVoxels() { return voxels_; }

  void Setup();
  void Print();
  void PrintActiveCellNum();

  // Set geometry
  void get_rectangle(T x0_, T x1_, T y0_, T y1_, int& i0, int& i1, int& j0,
                     int& j1);
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
  void get_Idx(std::vector<int>& idx, int flag);
  // add to indices of cells with flag
  // void add_Idx(std::vector<int> &idx, int flag);

  // quick set
  // set boundary cell flag to
  void set_bottom(int flag = 1)  // j = 1
  {
    for (int i = 1; i < Ni_ - 1; i++) {
      voxels_[Index2D::GetId(i, 1, Ni_)].setFlag(flag);
    }
  }
  void set_top(int flag = 1) {
    for (int i = 1; i < Ni_ - 1; i++) {
      voxels_[Index2D::GetId(i, Nj_ - 2, Ni_)].setFlag(flag);
    }
  }
  void set_left(int flag = 1) {
    for (int j = 1; j < Nj_ - 1; j++) {
      voxels_[Index2D::GetId(1, j, Ni_)].setFlag(flag);
    }
  }
  void set_right(int flag = 1) {
    for (int j = 1; j < Nj_ - 1; j++) {
      voxels_[Index2D::GetId(Ni_ - 2, j, Ni_)].setFlag(flag);
    }
  }
  void quickset_boundary(int flag = 1) {
    set_bottom(flag);
    set_top(flag);
    set_left(flag);
    set_right(flag);
  }
};