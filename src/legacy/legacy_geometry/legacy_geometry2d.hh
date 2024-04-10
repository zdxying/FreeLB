// legacy_geometry2d.h


#pragma once

#include "legacy/legacy_geometry/legacy_geometry2d.h"

//////////////////////////////////////////////////
// legacy code
///////////////////////////////////////////////////

template <typename T>
Geometry2DLegacy<T>::Geometry2DLegacy(int Ni, int Nj) : Ni_(Ni), Nj_(Nj) {
  N_ = Ni_ * Ni_;
  //
  voxels_.reserve(N_);
  Setup();
  Print();
}

template <typename T>
void Geometry2DLegacy<T>::Setup() {
  // init mesh
  int id = 0;
  for (int j = 0; j < Nj_; j++) {
    for (int i = 0; i < Ni_; i++) {
      id = Index2D::GetId(i, j, Ni_);
      Voxel<T, 2> vox(voxels_.size(), 0, i * T(1), j * T(1));
      // cell length = 1, Voxel(int flag, Args... args)
      voxels_.emplace_back(vox);
    }
  }
  Index2D::Traverse_Peripheral(Ni_, Nj_, 0,
                               [this](int id) { voxels_[id].setFlag(-1); });
}
template <typename T>
void Geometry2DLegacy<T>::Print() {
  std::cout << "[Geometry2DLegacy Statistics]:\n"
            << "Ni:       " << Ni_ << "\n"
            << "Nj:       " << Nj_ << "\n"
            << "N:        " << N_ << "\n"
            << std::endl;
}
template <typename T>
void Geometry2DLegacy<T>::PrintActiveCellNum() {
  int count = 0;
  for (int j = 1; j < Nj_ - 1; j++) {
    for (int i = 1; i < Ni_ - 1; i++) {
      int id = Index2D::GetId(i, j, Ni_);
      if (voxels_[id].getFlag() != -1) count++;
    }
  }
  std::cout << "ActiveCells: " << count << std::endl;
}
template <typename T>
void Geometry2DLegacy<T>::get_rectangle(T x0, T x1, T y0, T y1, int &i0,
                                        int &i1, int &j0, int &j1) {
  if (x0 < 1)
    i0 = 1;
  else
    i0 = int(x0) + 1;  // i.e. 1.5 -> 2

  if (x1 > Ni_ - 2)
    i1 = Ni_ - 2;
  else
    i1 = int(x1);  // i.e. 5.5 -> 5, when used in for loop, need to +1

  if (y0 < 1)
    j0 = 1;
  else
    j0 = int(y0) + 1;  // i.e. 1.5 -> 2

  if (y1 > Nj_ - 2)
    j1 = Nj_ - 2;
  else
    j1 = int(y1);
}

template <typename T>
void Geometry2DLegacy<T>::check_circle(T radius, T centre_x, T centre_y) {
  std::cout << "[checking circle]: ";
  int status = 0;
  // check if the centre is in the domain
  if (centre_x < 1 || centre_x > Ni_ - 2 || centre_y < 1 ||
      centre_y > Nj_ - 2) {
    std::cout << "circle centre is out of the domain" << std::endl;
    status = 1;
  }
  // check if the whole circle is in the domain
  T left = centre_x - radius;
  T right = centre_x + radius;
  T bottom = centre_y - radius;
  T top = centre_y + radius;
  if (left < 1 || right > Ni_ - 2 || bottom < 1 || top > Nj_ - 2) {
    std::cout << "circle is out of the domain" << std::endl;
    status = 1;
  }
  if (status == 0) {
    std::cout << "Correct" << std::endl;
  }
}
template <typename T>
void Geometry2DLegacy<T>::circle(T radius, T centre_x, T centre_y, int flag) {
  check_circle(radius, centre_x, centre_y);
  int i0, i1, j0, j1;  // get rectangle of the circle
  T left = centre_x - radius;
  T right = centre_x + radius;
  T bottom = centre_y - radius;
  T top = centre_y + radius;
  get_rectangle(left, right, bottom, top, i0, i1, j0, j1);
  T x_, y_;
  int id;
#pragma omp parallel for private(x_, y_, id) num_threads(Thread_Num)
  for (int j = j0; j < j1 + 1; j++) {
    for (int i = i0; i < i1 + 1; i++) {
      id = Index2D::GetId(i, j, Ni_);
      x_ = getVoxel(id)[0];
      y_ = getVoxel(id)[1];
      if (((x_ - centre_x) * (x_ - centre_x) +
           (y_ - centre_y) * (y_ - centre_y)) < radius * radius) {
        voxels_[id].getFlag() = flag;
      }
    }
  }
}

template <typename T>
void Geometry2DLegacy<T>::check_rectangle(T x0, T x1, T y0, T y1) {
  std::cout << "[checking rectangle]: " << std::endl;
  int status = 0;
  if (x0 < 1 || x1 > Ni_ - 2 || y0 < 1 || y1 > Nj_ - 2) {
    std::cout << "rectangle is out of the domain" << std::endl;
    status = 1;
  }
  if (status == 0) {
    std::cout << "Correct" << std::endl;
  }
}
template <typename T>
void Geometry2DLegacy<T>::rectangle(T x0, T x1, T y0, T y1, int flag) {
  check_rectangle(x0, x1, y0, y1);
  int i0, i1, j0, j1;  // get rectangle
  get_rectangle(x0, x1, y0, y1, i0, i1, j0, j1);
  T x_, y_;
  int id;
#pragma omp parallel for private(x_, y_, id) num_threads(Thread_Num)
  for (int j = j0; j < j1 + 1; j++) {
    for (int i = i0; i < i1 + 1; i++) {
      id = Index2D::GetId(i, j, Ni_);
      x_ = getVoxel(id)[0];
      y_ = getVoxel(id)[1];
      if (x_ > x0 && x_ < x1 && y_ > y0 && y_ < y1) {
        voxels_[id].getFlag() = flag;
      }
    }
  }
}

template <typename T>
void Geometry2DLegacy<T>::check_triangle(T x0, T x1, T x2, T y0, T y1, T y2) {
  std::cout << "[checking triangle]: " << std::endl;
  int status = 0;
  if (x0 < 1 || x0 > Ni_ - 2 || y0 < 1 || y0 > Nj_ - 2) status = 1;
  if (x1 < 1 || x1 > Ni_ - 2 || y1 < 1 || y1 > Nj_ - 2) status = 1;
  if (x2 < 1 || x2 > Ni_ - 2 || y2 < 1 || y2 > Nj_ - 2) status = 1;
  if (status == 0) {
    std::cout << "Correct" << std::endl;
  }
  if (status == 1) {
    std::cout << "triangle is out of the domain" << std::endl;
  }
}
template <typename T>
inline void Geometry2DLegacy<T>::get_Idx(std::vector<int> &idx, int flag) {
  // traverse all cells and find cells with flag
  Index2D::Traverse_Bulk(Ni_, Nj_, 0, [this, &idx, &flag](int id) {
    if (voxels_[id].getFlag() == flag) idx.emplace_back(id);
  });
}

template <typename T>
void Geometry2DLegacy<T>::triangle(T x0, T x1, T x2, T y0, T y1, T y2,
                                   int flag) {
  check_triangle(x0, x1, x2, y0, y1, y2);
  int i0, i1, j0, j1;
  T left = x0;
  if (x1 < left) left = x1;
  if (x2 < left) left = x2;
  T right = x0;
  if (x1 > right) right = x1;
  if (x2 > right) right = x2;
  T bottom = x0;
  if (x1 < bottom) bottom = x1;
  if (x2 < bottom) bottom = x2;
  T top = x0;
  if (x1 > top) top = x1;
  if (x2 > top) top = x2;
  get_rectangle(left, right, bottom, top, i0, i1, j0, j1);
  T x_, y_, t1, t2, t3;
  int id;
#pragma omp parallel for private(x_, y_, t1, t2, t3, id) num_threads(Thread_Num)
  for (int j = j0; j < j1 + 1; j++) {
    for (int i = i0; i < i1 + 1; i++) {
      id = Index2D::GetId(i, j, Ni_);
      x_ = getVoxel(id)[0];
      y_ = getVoxel(id)[1];
      // check if the point is inside the triangle
      t1 = (x_ - x0) * (y_ - y1) - (x_ - x1) * (y_ - y0);
      t2 = (x_ - x1) * (y_ - y2) - (x_ - x2) * (y_ - y1);
      t3 = (x_ - x2) * (y_ - y0) - (x_ - x0) * (y_ - y2);
      if (t1 * t2 >= 0 && t2 * t3 >= 0) {
        voxels_[id].getFlag() = flag;
      }
    }
  }
}