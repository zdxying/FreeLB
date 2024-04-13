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

// 2D finite difference solver
#pragma once

#include "data_struct/field_struct.h"
#include "util.h"
// taylor expansion of f(x,y) at (x0,y0):
// f(x,y) = f(x0,y0) +
//          (x-x0)*df/dx + (y-y0)*df/dy +
//          1/2*(x-x0)^2*d^2f/dx^2 + 1/2*(y-y0)^2*d^2f/dy^2 +
//          (x-x0)*(y-y0)*d^2f/dxdy + ...

// df/dx = (f(x0+dx,y0) - f(x0-dx,y0)) / (2*dx)
// d^2f/dx^2 = (f(x0+dx,y0) - 2*f(x0,y0) + f(x0-dx,y0)) / (dx^2)
// d^2f/dxdy = ((f(x0+dx,y0+dy) - f(x0+dx,y0-dy)) - (f(x0-dx,y0+dy) -
// f(x0-dx,y0-dy))) / (4*dx*dy)

template <typename T>
class FDM2D {
 private:
  int Ni;
  int Nj;
  // neighbor index {1, Ni, -1, -Ni, Ni + 1, Ni - 1, -Ni + 1, -Ni - 1}
  int Nbr[8];
  // cell length is assumed to be 1
  GenericArray<T> &f;

 public:
  FDM2D(int ni, int nj, GenericArray<T> &f_)
      : Ni(ni), Nj(nj), Nbr{1, Ni, -1, -Ni, Ni + 1, Ni - 1, -Ni + 1, -Ni - 1}, f(f_) {}
  ~FDM2D() {}
  // FDM parser
  // partial x, central difference
  inline T p_x(std::size_t id) const { return (f[id + Nbr[0]] - f[id + Nbr[2]]) * T(0.5); }
  // partial y, central difference
  inline T p_y(std::size_t id) const { return (f[id + Nbr[1]] - f[id + Nbr[3]]) * T(0.5); }
  // partial xx
  inline T p_xx(std::size_t id) const { return (f[id + Nbr[0]] - 2 * f[id] + f[id + Nbr[2]]); }
  // partial yy
  inline T p_yy(std::size_t id) const { return (f[id + Nbr[1]] - 2 * f[id] + f[id + Nbr[3]]); }
  // partial xy
  // corrected at 2023-11-16, it should be f[id + Nbr[7]] - f[id + Nbr[6]]
  // instead of f[id + Nbr[6]] - f[id + Nbr[7]]
  // (∂^2f)/∂x∂y=(((f(x+1,y+1)−f(x−1,y+1))−(f(x+1,y−1)−f(x−1,y−1))))/4
  inline T p_xy(std::size_t id) const {
    return (f[id + Nbr[7]] - f[id + Nbr[6]] - f[id + Nbr[5]] + f[id + Nbr[4]]) * T(0.25);
  }
  // get gradient: grad = {partial x, partial y}
  inline Vector<T, 2> grad(std::size_t id) const { return Vector<T, 2>{p_x(id), p_y(id)}; }
  // get Normalizedgradient = grad / |grad|
  inline Vector<T, 2> ngrad(std::size_t id) const {
    Vector<T, 2> grad_ = grad(id);
    grad_.normalize();
    return grad_;
  }
  // get suqare norm of gradient
  inline T gradnorm2(std::size_t id) const {
    const T x = p_x(id);
    const T y = p_y(id);
    return x * x + y * y;
  }
  // get norm of gradient
  inline T gradnorm(std::size_t id) const { return std::sqrt(gradnorm2(id)); }
};

template <typename T>
class Geometry3D;


// TODO: more generic implementation
// DO NOT SUPPORT D3Q15
template <typename T>
class FDM3D {
 private:
  int Ni;
  int Nj;
  int Nk;
  Geometry3D<T> &Geo;

  //   {0, 0, 0},  {1, 0, 0},    {-1, 0, 0}, {0, 1, 0},
  // {0, -1, 0}, {0, 0, 1},    {0, 0, -1},  // 0-6
  // {1, 1, 0},  {-1, -1, 0},  {1, 0, 1},  {-1, 0, -1},
  // {0, 1, 1},  {0, -1, -1},  {1, -1, 0}, {-1, 1, 0},
  // {1, 0, -1}, {-1, 0, 1},   {0, 1, -1}, {0, -1, 1},  // 7-18
  // {1, 1, 1},  {-1, -1, -1}, {1, 1, -1}, {-1, -1, 1},
  // {1, -1, 1}, {-1, 1, -1},  {-1, 1, 1}, {1, -1, -1}  // 19-26
  std::array<int, 27> DIdx;
  // cell length is assumed to be 1
  GenericArray<T> &f;

 public:
  FDM3D(Geometry3D<T> &geo, GenericArray<T> &f_)
      : Geo(geo), Ni(geo.getNx()), Nj(geo.getNy()), Nk(geo.getNz()), DIdx(geo.getDeltaIndex()),
        f(f_) {}
  ~FDM3D() {}
  // FDM parser
  // partial x, central difference
  inline T p_x(int id) const { return (f[id + DIdx[1]] - f[id + DIdx[2]]) * T(0.5); }
  // partial y, central difference
  inline T p_y(int id) const { return (f[id + DIdx[3]] - f[id + DIdx[4]]) * T(0.5); }
  // partial z, central difference
  inline T p_z(int id) const { return (f[id + DIdx[5]] - f[id + DIdx[6]]) * T(0.5); }
  // partial xx
  inline T p_xx(int id) const { return (f[id + DIdx[1]] + f[id + DIdx[2]] - 2 * f[id]); }
  // partial yy
  inline T p_yy(int id) const { return (f[id + DIdx[3]] + f[id + DIdx[4]] - 2 * f[id]); }
  // partial zz
  inline T p_zz(int id) const { return (f[id + DIdx[5]] + f[id + DIdx[6]] - 2 * f[id]); }
  // partial xy
  inline T p_xy(int id) const {
    return (f[id + DIdx[7]] + f[id + DIdx[8]] - f[id + DIdx[13]] - f[id + DIdx[14]]) * T(0.25);
  }
  // partial xz
  inline T p_xz(int id) const {
    return (f[id + DIdx[9]] + f[id + DIdx[10]] - f[id + DIdx[15]] - f[id + DIdx[16]]) * T(0.25);
  }
  // partial yz
  inline T p_yz(int id) const {
    return (f[id + DIdx[11]] + f[id + DIdx[12]] - f[id + DIdx[17]] - f[id + DIdx[18]]) * T(0.25);
  }
  // get gradient: grad = {partial x, partial y, partial z}
  inline Vector<T, 3> grad(int id) const { return Vector<T, 3>{p_x(id), p_y(id), p_z(id)}; }
  // get Normalizedgradient = grad / |grad|
  inline Vector<T, 3> ngrad(int id) const {
    Vector<T, 3> grad_ = grad(id);
    grad_.normalize();
    return grad_;
  }
  // Inabla(grad) = (partial_xx, partial_xy, partial_xz)^T
  Vector<T, 3> Inablagrad(int id) const { return Vector<T, 3>(p_xx(id), p_yy(id), p_zz(id)); }
};


// // 2D finite difference parser
// // nbr[8] = {1, Ni, -1, -Ni, Ni + 1, Ni - 1, -Ni - 1, -Ni + 1}
// template <typename T>
// struct FDMhelper_2D {
//   static inline void partial_x(T *temp, T *f, T dx, int *nbr, int id) {
//     temp[id] = (f[id + nbr[0]] - f[id + nbr[2]]) / (2 * dx);
//   }
//   static inline T partial_x(T *f, T dx, int *nbr, int id) {
//     return (f[id + nbr[0]] - f[id + nbr[2]]) / (2 * dx);
//   }
//   static inline void partial_y(T *temp, T *f, T dy, int *nbr, int id) {
//     temp[id] = (f[id + nbr[1]] - f[id + nbr[3]]) / (2 * dy);
//   }
//   static inline T partial_y(T *f, T dy, int *nbr, int id) {
//     return (f[id + nbr[1]] - f[id + nbr[3]]) / (2 * dy);
//   }
//   static inline void partial_xx(T *temp, T *f, T dx, int *nbr, int id) {
//     temp[id] = (f[id + nbr[0]] - 2 * f[id] + f[id + nbr[2]]) / (dx * dx);
//   }
//   static inline T partial_xx(T *f, T dx, int *nbr, int id) {
//     return (f[id + nbr[0]] - 2 * f[id] + f[id + nbr[2]]) / (dx * dx);
//   }
//   static inline void partial_yy(T *temp, T *f, T dy, int *nbr, int id) {
//     temp[id] = (f[id + nbr[1]] - 2 * f[id] + f[id + nbr[3]]) / (dy * dy);
//   }
//   static inline T partial_yy(T *f, T dy, int *nbr, int id) {
//     return (f[id + nbr[1]] - 2 * f[id] + f[id + nbr[3]]) / (dy * dy);
//   }
//   static inline void partial_xy(T *temp, T *f, T dx, T dy, int *nbr, int id) {
//     temp[id] = ((f[id + nbr[4]] - f[id + nbr[5]]) -
//                 (f[id + nbr[6]] - f[id + nbr[7]])) /
//                (4 * dx * dy);
//   }
//   static inline T partial_xy(T *f, T dx, T dy, int *nbr, int id) {
//     return ((f[id + nbr[4]] - f[id + nbr[5]]) -
//             (f[id + nbr[6]] - f[id + nbr[7]])) /
//            (4 * dx * dy);
//   }
// };
