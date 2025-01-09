/* This file is part of FreeLB, modified from OpenLB's vector.h, with the following
 * copyright notice:
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
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// Vector.h

#pragma once

// sin cos sqrt
#include <cmath>

#include <array>

#include "head.h"

// the Vector here is not a container, but a mathematical vector
// class Vector stores an array of type T with size D(2 or 3)

template <typename T, unsigned int D>
class Vector {
 private:
  T _data[D];
  // std::array<T, D> _data;

 public:
  using value_type = T;
  static constexpr unsigned int vector_dim = D;

  // constructors
  // __any__ constexpr Vector() : _data{} {}
  __any__ constexpr Vector() {
    for (unsigned int i = 0; i < D; ++i) _data[i] = T{};
  }
  // __any__ constexpr Vector(T value) { _data.fill(value); }
  __any__ constexpr Vector(T value) {
    for (unsigned int i = 0; i < D; ++i) _data[i] = value;
  }

  template <typename... Args>
  __any__ constexpr Vector(Args... args) : _data{args...} {}
  // Copy constructor
  __any__ constexpr Vector(const Vector<T, D> &vec) {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = vec[i];
    }
  }
  // Move constructor
  __any__ constexpr Vector(Vector<T, D> &&vec) {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = vec[i];
    }
  }

  __any__ constexpr Vector(const T *array) {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = array[i];
    }
  }

  // Conversion constructor
  template <typename U>
  __any__ constexpr Vector(const Vector<U, D> &vec) {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = static_cast<T>(vec[i]);
    }
  }

  ~Vector() = default;

  // return pointer of the first element
  __any__ constexpr T *data() { return _data; }
  __any__ constexpr const T *data() const { return _data; }

  // (mathematical) operation
  // return square of the norm of the vector
  __any__ constexpr T getnorm2() const {
    T result = 0;
    for (unsigned int i = 0; i < D; ++i) {
      result += _data[i] * _data[i];
    }
    return result;
  }
  // return norm of the vector
  __any__ constexpr T getnorm() const { return std::sqrt(getnorm2()); }
  // normalize the vector
  __any__ void normalize() {
    T norm = this->getnorm();
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] /= norm;
    }
  }
  // get normalized vector
  __any__ constexpr Vector<T, D> getnormalize() const {
    Vector<T, D> result;
    T norm = this->getnorm();
    for (unsigned int i = 0; i < D; ++i) {
      result[i] = _data[i] / norm;
    }
    return result;
  }
  __any__ void abs() {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = std::abs(_data[i]);
    }
  }
  __any__ constexpr Vector<T, D> getabs() const {
    Vector<T, D> result;
    for (unsigned int i = 0; i < D; ++i) {
      result[i] = std::abs(_data[i]);
    }
    return result;
  }
  __any__ void clear() { for (unsigned int i = 0; i < D; ++i) _data[i] = T{}; }

  // operators
  // []: return reference of the i th element
  __any__ constexpr T &operator[](unsigned int i) { return _data[i]; }
  // []: return const reference of the i th element
  __any__ constexpr const T &operator[](unsigned int i) const { return _data[i]; }
  // Copy assignment operator
  __any__ constexpr Vector &operator=(const Vector<T, D> &vec) {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = vec[i];
    }
    return *this;
  }
  // Move assignment operator
  __any__ constexpr Vector &operator=(const Vector<T, D> &&vec) {
    for (unsigned int i = 0; i < D; ++i) {
      _data[i] = vec[i];
    }
    return *this;
  }
};

// operators
// +: return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} + U{}), D> operator+(const Vector<T, D> &a,
                                                           const Vector<U, D> &b) {
  Vector<decltype(T{} + U{}), D> result;
  for (unsigned int i = 0; i < D; ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} + U{}), D> operator+(const Vector<T, D> &a, U b) {
  return a + Vector<U, D>(b);
}
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} + U{}), D> operator+(T a, const Vector<U, D> &b) {
  return Vector<T, D>(a) + b;
}
// +=: return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} + U{}), D> &operator+=(Vector<T, D> &a,
                                                             const Vector<U, D> &b) {
  for (unsigned int i = 0; i < D; ++i) {
    a[i] += b[i];
  }
  return a;
}
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} + U{}), D> &operator+=(Vector<T, D> &a, U b) {
  for (unsigned int i = 0; i < D; ++i) {
    a[i] += b;
  }
  return a;
}
// -: return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} - U{}), D> operator-(const Vector<T, D> &a,
                                                           const Vector<U, D> &b) {
  Vector<decltype(T{} - U{}), D> result;
  for (unsigned int i = 0; i < D; ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} - U{}), D> operator-(const Vector<T, D> &a, U b) {
  return a - Vector<U, D>(b);
}
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} - U{}), D> operator-(T a, const Vector<U, D> &b) {
  return Vector<T, D>(a) - b;
}
// -=: return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} - U{}), D> &operator-=(Vector<T, D> &a,
                                                             const Vector<U, D> &b) {
  for (unsigned int i = 0; i < D; ++i) {
    a[i] -= b[i];
  }
  return a;
}
// *: inner(dot) product, return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
__any__ constexpr decltype(T{} * U{}) operator*(const Vector<T, D> &a,
                                                const Vector<U, D> &b) {
  decltype(T{} * U{}) result{};
  for (unsigned int i = 0; i < D; ++i) {
    result += a[i] * b[i];
  }
  return result;
}
// *: Vector multiplies a scalar, return type will be deduced by
// decltype(c++11)
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} * U{}), D> operator*(const Vector<T, D> &a, U b) {
  Vector<decltype(T{} * U{}), D> result;
  for (unsigned int i = 0; i < D; ++i) {
    result[i] = a[i] * b;
  }
  return result;
}
// scalar multiplies a vector
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} * U{}), D> operator*(T a, const Vector<U, D> &b) {
  Vector<decltype(T{} * U{}), D> result;
  for (unsigned int i = 0; i < D; ++i) {
    result[i] = a * b[i];
  }
  return result;
}
// /: Vector divides a scalar, return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} / U{}), D> operator/(const Vector<T, D> &a, U b) {
  Vector<decltype(T{} / U{}), D> result;
  for (unsigned int i = 0; i < D; ++i) {
    result[i] = a[i] / b;
  }
  return result;
}
// /=: Vector divides a scalar, return type will be deduced by decltype(c++11)
template <typename T, typename U, unsigned int D>
__any__ constexpr Vector<decltype(T{} / U{}), D> &operator/=(Vector<T, D> &a, U b) {
  for (unsigned int i = 0; i < D; ++i) {
    a[i] /= b;
  }
  return a;
}

// other mathematical operations
// cross product 2D
template <typename T, typename U>
__any__ constexpr decltype(T{} * U{}) CrossProduct2(const Vector<T, 2> &a,
                                                    const Vector<U, 2> &b) {
  decltype(T{} * U{}) result = a[0] * b[1] - a[1] * b[0];
  return result;
}
// cross product 3D
template <typename T, typename U>
__any__ constexpr Vector<decltype(T{} * U{}), 3> CrossProduct3(const Vector<T, 3> &a,
                                                               const Vector<U, 3> &b) {
  Vector<decltype(T{} * U{}), 3> result(
    a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
  // result[0] = a[1] * b[2] - a[2] * b[1];
  // result[1] = a[2] * b[0] - a[0] * b[2];
  // result[2] = a[0] * b[1] - a[1] * b[0];
  return result;
}

template <typename T, typename U, unsigned int D>
__any__ constexpr auto CrossProduct(const Vector<T, D> &a, const Vector<U, D> &b) {
  if constexpr (D == 2) {
    return CrossProduct2(a, b);
  } else {
    return CrossProduct3(a, b);
  }
}

// get square of distance between 2 vectors
template <typename T, typename U, unsigned int D>
__any__ constexpr decltype(T{} * U{}) GetDist2(const Vector<T, D> &a,
                                               const Vector<U, D> &b) {
  decltype(T{} * U{}) result{};
  for (unsigned int i = 0; i < D; ++i) {
    result += std::pow(a[i] - b[i], 2);
  }
  return result;
}
// get distance between 2 vectors
template <typename T, typename U, unsigned int D>
__any__ constexpr decltype(T{} * U{}) GetDist(const Vector<T, D> &a,
                                              const Vector<U, D> &b) {
  return std::sqrt(GetDist2(a, b));
}

template <typename T, unsigned int D>
__any__ constexpr Vector<T, D> getnormalized(const Vector<T, D> &a) {
  Vector<T, D> result = a;
  T norm = a.getnorm();
  for (unsigned int i = 0; i < D; ++i) {
    result[i] /= norm;
  }
  return result;
}

// get relative location of (x, y) to (x0, y0) in a rotated
// coordinate system，theta is the angle of rotation counterclockwise
// x' = (x - x0) * cos(theta) + (y - y0) * std::sin(theta)
// y' = -(x - x0) * std::sin(theta) + (y - y0) * cos(theta)
template <typename T>
__any__ constexpr Vector<T, 2> getRLoc(const Vector<T, 2> &loc, const Vector<T, 2> &loc0,
                                       T theta = T(0)) {
  T dx = loc[0] - loc0[0];
  T dy = loc[1] - loc0[1];
  return Vector<T, 2>(dx * std::cos(theta) + dy * std::sin(theta),
                      -dx * std::sin(theta) + dy * std::cos(theta));
}
template <typename T>
__any__ void getRLoc(const Vector<T, 2> &loc, const Vector<T, 2> &loc0,
                     Vector<T, 2> &result, T theta = T(0)) {
  T dx = loc[0] - loc0[0];
  T dy = loc[1] - loc0[1];
  result[0] = dx * std::cos(theta) + dy * std::sin(theta);
  result[1] = -dx * std::sin(theta) + dy * std::cos(theta);
  // return Vector<T, 2>(dx * std::cos(theta) + dy * std::sin(theta),
  //                     -dx * std::sin(theta) + dy * std::cos(theta));
}
template <typename T>
__any__ constexpr Vector<T, 2> getGLoc(const Vector<T, 2> &loc, const Vector<T, 2> &glob0,
                                       T theta = T(0)) {
  T dx = loc[0];
  T dy = loc[1];
  return Vector<T, 2>(dx * std::cos(theta) - dy * std::sin(theta),
                      dx * std::sin(theta) + dy * std::cos(theta)) +
         glob0;
}
// template <typename T>
// Vector<T, 2> getRLoc(T x, T y, const Vector<T, 2> &loc0, T theta) {
//   T dx = x - loc0[0];
//   T dy = y - loc0[1];
//   return Vector<T, 2>(dx * std::cos(theta) + dy * std::sin(theta),
//                       -dx * std::sin(theta) + dy * std::cos(theta));
// }

template <typename T>
__any__ int getQuad2D(const Vector<T, 2> &a) {
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
