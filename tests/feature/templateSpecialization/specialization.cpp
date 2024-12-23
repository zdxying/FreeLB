/* This file is part of FreeLB
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

#include "freelb.h"
#include "freelb.hh"

#include "template.h"

using T = FLOAT;



int main() {
  std::vector<T> cellsD2Q9(D2Q9<T>::q, T{1});
  std::vector<T> cellsD3Q19(D3Q19<T>::q, T{1});

  using namespace tempspectest;

  std::cout << rho<D2Q9<T>>::get(cellsD2Q9) << std::endl;
  std::cout << rho<D3Q19<T>>::get(cellsD3Q19) << std::endl;

}
