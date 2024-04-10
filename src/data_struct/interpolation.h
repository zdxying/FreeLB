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

// interpolation.h

#pragma once

#include <array>

template <unsigned int D>
struct BasicInterp {
  // target id
  std::size_t id;
  // source
  std::array<std::size_t, D> src;

  BasicInterp(std::size_t id) : id(id), src({0}) {}
  BasicInterp(std::size_t id, const std::array<std::size_t, D>& src) : id(id), src(src) {}
  template <typename... Args>
  BasicInterp(std::size_t id, Args... args) : id(id), src({args...}) {}

  std::size_t operator[](int i) const { return src[i]; }
};