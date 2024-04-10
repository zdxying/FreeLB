// geometry class
// legacy 
#pragma once
#include "geometry_.h"

template <typename T>
Geometry2D_legacy<T>::Geometry2D_legacy(int Ni_, int Nj_) : Ni(Ni_), Nj(Nj_) {
  N = Ni * Nj;
  x = new T[N];
  y = new T[N];
  Flag = new int[N];  // 0: fluid, -1: solid, 1: mush zone
  Setup();
  Print();
}
template <typename T>
Geometry2D_legacy<T>::~Geometry2D_legacy() {
  delete[] x;
  delete[] y;
  delete[] Flag;
}

template <typename T>
void Geometry2D_legacy<T>::Setup() {
  // init mesh
  for (int j = 0; j < Nj; j++) {
    for (int i = 0; i < Ni; i++) {
      int id = GetId(i, j, Ni);
      // use lattice unit: cell length = 1
      x[id] = i * T(1);
      y[id] = j * T(1);
      Flag[id] = 0;
    }
  }
  // periphery cells assumed to be solid, convenient for boundary neighbor
  // search
  Traverse_Peripheral(Ni, Nj, 0, [this](int id) { Flag[id] = -1; });
}
template <typename T>
void Geometry2D_legacy<T>::Print() {
  std::cout << "[Geometry2D_legacy Statistics]:\n"
            << "Ni:       " << Ni << "\n"
            << "Nj:       " << Nj << "\n"
            << "N:        " << N << "\n"
            << std::endl;
}
template <typename T>
void Geometry2D_legacy<T>::PrintActiveCellNum() {
  int count = 0;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      int id = GetId(i, j, Ni);
      if (Flag[id] != -1) count++;
    }
  }
  std::cout << "ActiveCells: " << count << std::endl;
}
template <typename T>
void Geometry2D_legacy<T>::get_rectangle(T x0, T x1, T y0, T y1, int &i0, int &i1,
                                  int &j0, int &j1) {
  if (x0 < 1)
    i0 = 1;
  else
    i0 = int(x0) + 1;  // i.e. 1.5 -> 2

  if (x1 > Ni - 2)
    i1 = Ni - 2;
  else
    i1 = int(x1);  // i.e. 5.5 -> 5, when used in for loop, need to +1

  if (y0 < 1)
    j0 = 1;
  else
    j0 = int(y0) + 1;  // i.e. 1.5 -> 2

  if (y1 > Nj - 2)
    j1 = Nj - 2;
  else
    j1 = int(y1);
}

template <typename T>
void Geometry2D_legacy<T>::check_circle(T radius, T centre_x, T centre_y) {
  std::cout << "[checking circle]: ";
  int status = 0;
  // check if the centre is in the domain
  if (centre_x < 1 || centre_x > Ni - 2 || centre_y < 1 || centre_y > Nj - 2) {
    std::cout << "circle centre is out of the domain" << std::endl;
    status = 1;
  }
  // check if the whole circle is in the domain
  T left = centre_x - radius;
  T right = centre_x + radius;
  T bottom = centre_y - radius;
  T top = centre_y + radius;
  if (left < 1 || right > Ni - 2 || bottom < 1 || top > Nj - 2) {
    std::cout << "circle is out of the domain" << std::endl;
    status = 1;
  }
  if (status == 0) {
    std::cout << "Correct" << std::endl;
  }
}
template <typename T>
void Geometry2D_legacy<T>::circle(T radius, T centre_x, T centre_y, int flag) {
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
      id = GetId(i, j, Ni);
      x_ = x[id];
      y_ = y[id];
      if (((x_ - centre_x) * (x_ - centre_x) +
           (y_ - centre_y) * (y_ - centre_y)) < radius * radius) {
        Flag[id] = flag;
      }
    }
  }
}

template <typename T>
void Geometry2D_legacy<T>::check_rectangle(T x0, T x1, T y0, T y1) {
  std::cout << "[checking rectangle]: " << std::endl;
  int status = 0;
  if (x0 < 1 || x1 > Ni - 2 || y0 < 1 || y1 > Nj - 2) {
    std::cout << "rectangle is out of the domain" << std::endl;
    status = 1;
  }
  if (status == 0) {
    std::cout << "Correct" << std::endl;
  }
}
template <typename T>
void Geometry2D_legacy<T>::rectangle(T x0, T x1, T y0, T y1, int flag) {
  check_rectangle(x0, x1, y0, y1);
  int i0, i1, j0, j1;  // get rectangle
  get_rectangle(x0, x1, y0, y1, i0, i1, j0, j1);
  T x_, y_;
  int id;
#pragma omp parallel for private(x_, y_, id) num_threads(Thread_Num)
  for (int j = j0; j < j1 + 1; j++) {
    for (int i = i0; i < i1 + 1; i++) {
      id = GetId(i, j, Ni);
      x_ = x[id];
      y_ = y[id];
      if (x_ > x0 && x_ < x1 && y_ > y0 && y_ < y1) {
        Flag[id] = flag;
      }
    }
  }
}

template <typename T>
void Geometry2D_legacy<T>::check_triangle(T x0, T x1, T x2, T y0, T y1, T y2) {
  std::cout << "[checking triangle]: " << std::endl;
  int status = 0;
  if (x0 < 1 || x0 > Ni - 2 || y0 < 1 || y0 > Nj - 2) status = 1;
  if (x1 < 1 || x1 > Ni - 2 || y1 < 1 || y1 > Nj - 2) status = 1;
  if (x2 < 1 || x2 > Ni - 2 || y2 < 1 || y2 > Nj - 2) status = 1;
  if (status == 0) {
    std::cout << "Correct" << std::endl;
  }
  if (status == 1) {
    std::cout << "triangle is out of the domain" << std::endl;
  }
}
template <typename T>
inline void Geometry2D_legacy<T>::get_Idx(std::vector<int> &idx, int flag) {
  // traverse all cells and find cells with flag
  Traverse_Bulk(Ni, Nj, 0, [this, &idx, &flag](int id) {
    if (Flag[id] == flag) idx.emplace_back(id);
  });
}

template <typename T>
void Geometry2D_legacy<T>::triangle(T x0, T x1, T x2, T y0, T y1, T y2, int flag) {
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
      id = GetId(i, j, Ni);
      x_ = x[id];
      y_ = y[id];
      // check if the point is inside the triangle
      t1 = (x_ - x0) * (y_ - y1) - (x_ - x1) * (y_ - y0);
      t2 = (x_ - x1) * (y_ - y2) - (x_ - x2) * (y_ - y1);
      t3 = (x_ - x2) * (y_ - y0) - (x_ - x0) * (y_ - y2);
      if (t1 * t2 >= 0 && t2 * t3 >= 0) {
        Flag[id] = flag;
      }
    }
  }
}