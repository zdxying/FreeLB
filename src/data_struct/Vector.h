/* This file is part of FreeLB, modified from OpenLB's vector.h, with the following copyright notice:
 *
 * // start of the original OpenLB's copyright notice
 * 
 * This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Asher Zarth, Mathias J. Krause, Albert Mink
 *                2020 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 * 
 * // end of the original OpenLB's copyright notice
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

// Vector.h

#pragma once

#include <math.h>

#include <array>

// the Vector here is not a container, but a mathematical vector
// class Vector stores an array of type T with size D(2 or 3)

template <typename T, unsigned int D>
class Vector {
 private:
  // T _data[D];
  std::array<T, D> _data;

 public:
  using value_type = T;
  static constexpr unsigned int vector_dim = D;

  // constructors
  constexpr Vector() { _data.fill(T(0)); }

  constexpr Vector(T value) { _data.fill(value); }

  template <typename... Args>
  constexpr Vector(Args... args) : _data{args...} {}

  constexpr Vector(const Vector<T, D> &vec) {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = vec[i];
    }
  }

  constexpr Vector(const T *array) {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = array[i];
    }
  }

  ~Vector() = default;

  // return pointer of the first element
  constexpr T *data() { return _data.data(); }
  constexpr const T* data() const { return _data.data(); }

  // (mathematical) operation
  // return square of the norm of the vector
  constexpr T getnorm2() const {
    T result = 0;
    for (unsigned int i = 0; i < D; ++i) {
      result += _data[i] * _data[i];
    }
    return result;
  }
  // return norm of the vector
  constexpr T getnorm() const { return sqrt(getnorm2()); }
  // normalize the vector
  void normalize() {
    T norm = this->getnorm();
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] /= norm;
    }
  }
  // get normalized vector
  constexpr Vector<T, D> getnormalize() const {
    Vector<T, D> result;
    T norm = this->getnorm();
    for (unsigned int i = 0; i < D; ++i) {
      result[i] = _data[i] / norm;
    }
    return result;
  }
  void abs() {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = std::abs(_data[i]);
    }
  }
  constexpr Vector<T, D> getabs() const {
    Vector<T, D> result;
    for (unsigned int i = 0; i < D; ++i) {
      result[i] = std::abs(_data[i]);
    }
    return result;
  }
  void clear() { _data.fill(T(0)); }

  // operators
  // []: return reference of the i th element
  constexpr T &operator[](unsigned int i) { return _data[i]; }
  // []: return const reference of the i th element
  constexpr const T &operator[](unsigned int i) const { return _data[i]; }

  constexpr Vector &operator=(const Vector<T, D> &vec) {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = vec[i];
    }
    return *this;
  }
  constexpr Vector &operator=(const Vector<T, D> &&vec) {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = vec[i];
    }
    return *this;
  }
};

// operators
// +: return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
constexpr Vector<decltype(T{} + U{}), D> operator+(const Vector<T, D> &a, const Vector<U, D> &b) {
  Vector<decltype(T{} + U{}), D> result;
  for (unsigned int i = 0; i < D; ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}
template <typename T, typename U, unsigned int D>
constexpr Vector<decltype(T{} + U{}), D> operator+(const Vector<T, D> &a, U b) {
  return a + Vector<U, D>(b);
}
template <typename T, typename U, unsigned int D>
Vector<decltype(T{} + U{}), D> operator+(T a, const Vector<U, D> &b) {
  return Vector<T, D>(a) + b;
}
// +=: return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
constexpr Vector<decltype(T{} + U{}), D> &operator+=(Vector<T, D> &a, const Vector<U, D> &b) {
  for (unsigned int i = 0; i < D; ++i) {
    a[i] += b[i];
  }
  return a;
}
template <typename T, typename U, unsigned int D>
constexpr Vector<decltype(T{} + U{}), D> &operator+=(Vector<T, D> &a, U b) {
  for (unsigned int i = 0; i < D; ++i) {
    a[i] += b;
  }
  return a;
}
// -: return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
constexpr Vector<decltype(T{} - U{}), D> operator-(const Vector<T, D> &a, const Vector<U, D> &b) {
  Vector<decltype(T{} - U{}), D> result;
  for (unsigned int i = 0; i < D; ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}
template <typename T, typename U, unsigned int D>
constexpr Vector<decltype(T{} - U{}), D> operator-(const Vector<T, D> &a, U b) {
  return a - Vector<U, D>(b);
}
template <typename T, typename U, unsigned int D>
constexpr Vector<decltype(T{} - U{}), D> operator-(T a, const Vector<U, D> &b) {
  return Vector<T, D>(a) - b;
}
// *: inner(dot) product, return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
constexpr decltype(T{} * U{}) operator*(const Vector<T, D> &a, const Vector<U, D> &b) {
  decltype(T{} * U{}) result{};
  for (unsigned int i = 0; i < D; ++i) {
    result += a[i] * b[i];
  }
  return result;
}
// *: Vector multiplies a scaler, return type will be deduced by
// decltype(c++11)
template <typename T, typename U, unsigned int D>
constexpr Vector<decltype(T{} * U{}), D> operator*(const Vector<T, D> &a, U b) {
  Vector<decltype(T{} * U{}), D> result;
  for (unsigned int i = 0; i < D; ++i) {
    result[i] = a[i] * b;
  }
  return result;
}
// scaler multiplies a vector
template <typename T, typename U, unsigned int D>
constexpr Vector<decltype(T{} * U{}), D> operator*(T a, const Vector<U, D> &b) {
  Vector<decltype(T{} * U{}), D> result;
  for (unsigned int i = 0; i < D; ++i) {
    result[i] = a * b[i];
  }
  return result;
}
// /: Vector divides a scaler, return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
constexpr Vector<decltype(T{} / U{}), D> operator/(const Vector<T, D> &a, U b) {
  Vector<decltype(T{} / U{}), D> result;
  for (unsigned int i = 0; i < D; ++i) {
    result[i] = a[i] / b;
  }
  return result;
}

// other mathematical operations
// cross product 2D
template <typename T, typename U>
constexpr decltype(T{} * U{}) CrossProduct2(const Vector<T, 2> &a, const Vector<U, 2> &b) {
  decltype(T{} * U{}) result = a[0] * b[1] - a[1] * b[0];
  return result;
}
// cross product 3D
template <typename T, typename U>
constexpr Vector<decltype(T{} * U{}), 3> CrossProduct3(const Vector<T, 3> &a,
                                                       const Vector<U, 3> &b) {
  Vector<decltype(T{} * U{}), 3> result(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                                        a[0] * b[1] - a[1] * b[0]);
  // result[0] = a[1] * b[2] - a[2] * b[1];
  // result[1] = a[2] * b[0] - a[0] * b[2];
  // result[2] = a[0] * b[1] - a[1] * b[0];
  return result;
}

template <typename T, typename U, unsigned int D>
constexpr auto CrossProduct(const Vector<T, D> &a, const Vector<U, D> &b) {
  if constexpr (D == 2) {
    return CrossProduct2(a, b);
  } else {
    return CrossProduct3(a, b);
  }
}

// get square of distance between 2 vectors
template <typename T, typename U, unsigned int D>
constexpr decltype(T{} * U{}) GetDist2(const Vector<T, D> &a, const Vector<U, D> &b) {
  decltype(T{} * U{}) result{};
  for (unsigned int i = 0; i < D; ++i) {
    result += pow(a[i] - b[i], 2);
  }
  return result;
}
// get distance between 2 vectors
template <typename T, typename U, unsigned int D>
constexpr decltype(T{} * U{}) GetDist(const Vector<T, D> &a, const Vector<U, D> &b) {
  return sqrt(GetDist2(a, b));
}

template <typename T, unsigned int D>
constexpr Vector<T, D> getnormalized(const Vector<T, D> &a) {
  Vector<T, D> result = a;
  T norm = a.getnorm();
  for (unsigned int i = 0; i < D; ++i) {
    result[i] /= norm;
  }
  return result;
}

// get relative location of (x, y) to (x0, y0) in a rotated
// coordinate system，theta is the angle of rotation counterclockwise
// x' = (x - x0) * cos(theta) + (y - y0) * sin(theta)
// y' = -(x - x0) * sin(theta) + (y - y0) * cos(theta)
template <typename T>
constexpr Vector<T, 2> getRLoc(const Vector<T, 2> &loc, const Vector<T, 2> &loc0, T theta = T(0)) {
  T dx = loc[0] - loc0[0];
  T dy = loc[1] - loc0[1];
  return Vector<T, 2>(dx * cos(theta) + dy * sin(theta), -dx * sin(theta) + dy * cos(theta));
}
template <typename T>
void getRLoc(const Vector<T, 2> &loc, const Vector<T, 2> &loc0, Vector<T, 2> &result,
             T theta = T(0)) {
  T dx = loc[0] - loc0[0];
  T dy = loc[1] - loc0[1];
  result[0] = dx * cos(theta) + dy * sin(theta);
  result[1] = -dx * sin(theta) + dy * cos(theta);
  // return Vector<T, 2>(dx * cos(theta) + dy * sin(theta),
  //                     -dx * sin(theta) + dy * cos(theta));
}
template <typename T>
constexpr Vector<T, 2> getGLoc(const Vector<T, 2> &loc, const Vector<T, 2> &glob0, T theta = T(0)) {
  T dx = loc[0];
  T dy = loc[1];
  return Vector<T, 2>(dx * cos(theta) - dy * sin(theta), dx * sin(theta) + dy * cos(theta)) + glob0;
}
// template <typename T>
// Vector<T, 2> getRLoc(T x, T y, const Vector<T, 2> &loc0, T theta) {
//   T dx = x - loc0[0];
//   T dy = y - loc0[1];
//   return Vector<T, 2>(dx * cos(theta) + dy * sin(theta),
//                       -dx * sin(theta) + dy * cos(theta));
// }

template <typename T>
int getQuad2D(const Vector<T, 2> &a) {
  // 1 | 0
  // -----
  // 2 | 3
  if (a[0] >= 0 && a[1] >= 0)
    return 0;
  else if (a[0] < 0 && a[1] >= 0)
    return 1;
  else if (a[0] < 0 && a[1] < 0)
    return 2;
  else
    return 3;
}
