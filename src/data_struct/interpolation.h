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
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 * 
 */

// interpolation.h
// Coarse/fine interpolation weights and the (IntpSource-based) average /
// interpolation kernels.  Lives in its own header because two unrelated
// layers need them:
//   * data_struct/field.h      (Interpolation() packing overloads)
//   * parallel/communicator.h  (MPI send/recv bookkeeping)
// Keeping them here lets communicator.h stay independent of the (large)
// field.h -- it only ever needed the weight table.

#pragma once

#include <array>
#include <cstddef>

#include "utils/alias.h"  // IntpSource

// get predefined interp weight
template <typename T>
constexpr std::array<std::array<T, 4>, 4> getIntpWeight2D() {
  return 
  {{{T{0.0625}, T{0.1875}, T{0.1875}, T{0.5625}},
    {T{0.1875}, T{0.0625}, T{0.5625}, T{0.1875}},
    {T{0.1875}, T{0.5625}, T{0.0625}, T{0.1875}},
    {T{0.5625}, T{0.1875}, T{0.1875}, T{0.0625}}}};
}

template <typename T>
constexpr std::array<std::array<T, 8>, 8> getIntpWeight3D() {
  return
{{{T{0.015625}, T{0.046875}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.140625}, T{0.421875}},
  {T{0.046875}, T{0.015625}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.421875}, T{0.140625}},
  {T{0.046875}, T{0.140625}, T{0.015625}, T{0.046875}, T{0.140625}, T{0.421875}, T{0.046875}, T{0.140625}},
  {T{0.140625}, T{0.046875}, T{0.046875}, T{0.015625}, T{0.421875}, T{0.140625}, T{0.140625}, T{0.046875}},
  {T{0.046875}, T{0.140625}, T{0.140625}, T{0.421875}, T{0.015625}, T{0.046875}, T{0.046875}, T{0.140625}},
  {T{0.140625}, T{0.046875}, T{0.421875}, T{0.140625}, T{0.046875}, T{0.015625}, T{0.140625}, T{0.046875}},
  {T{0.140625}, T{0.421875}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.015625}, T{0.046875}},
  {T{0.421875}, T{0.140625}, T{0.140625}, T{0.046875}, T{0.140625}, T{0.046875}, T{0.046875}, T{0.015625}}}};
}

template <typename T, unsigned int D>
constexpr auto getIntpWeight() {
  if constexpr (D == 2) {
    return getIntpWeight2D<T>();
  } else {
    return getIntpWeight3D<T>();
  }
}

// simple average over the neighbouring fine cells (4 in 2D, 8 in 3D)
template <typename FloatType, unsigned int Dim, typename ArrayType>
typename ArrayType::value_type getAverage(const ArrayType& Arr,
                                          const IntpSource<Dim>& src) {
  using datatype = typename ArrayType::value_type;
  datatype Aver = datatype{};
  if constexpr (Dim == 2) {
    Aver = (Arr[src[0]] + Arr[src[1]] + Arr[src[2]] + Arr[src[3]]) * FloatType(0.25);
  } else if constexpr (Dim == 3) {
    Aver = (Arr[src[0]] + Arr[src[1]] + Arr[src[2]] + Arr[src[3]] + Arr[src[4]] +
            Arr[src[5]] + Arr[src[6]] + Arr[src[7]]) *
           FloatType(0.125);
  }
  return Aver;
}

// interpolation from neighboring fine cells using precomputed weights
template <unsigned int D, typename FloatType, unsigned int Dim, typename ArrayType>
typename ArrayType::value_type getInterpolation(const ArrayType& Arr,
                                                const IntpSource<Dim>& src) {
  using datatype = typename ArrayType::value_type;
  datatype Intp = datatype{};
  if constexpr (Dim == 2) {
    Intp = Arr[src[0]] * getIntpWeight<FloatType, Dim>()[D][0] +
           Arr[src[1]] * getIntpWeight<FloatType, Dim>()[D][1] +
           Arr[src[2]] * getIntpWeight<FloatType, Dim>()[D][2] +
           Arr[src[3]] * getIntpWeight<FloatType, Dim>()[D][3];
  } else if constexpr (Dim == 3) {
    Intp = Arr[src[0]] * getIntpWeight<FloatType, Dim>()[D][0] +
           Arr[src[1]] * getIntpWeight<FloatType, Dim>()[D][1] +
           Arr[src[2]] * getIntpWeight<FloatType, Dim>()[D][2] +
           Arr[src[3]] * getIntpWeight<FloatType, Dim>()[D][3] +
           Arr[src[4]] * getIntpWeight<FloatType, Dim>()[D][4] +
           Arr[src[5]] * getIntpWeight<FloatType, Dim>()[D][5] +
           Arr[src[6]] * getIntpWeight<FloatType, Dim>()[D][6] +
           Arr[src[7]] * getIntpWeight<FloatType, Dim>()[D][7];
  }
  return Intp;
}
