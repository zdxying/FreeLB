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

#include <vector>
#include <iostream>
#include "lbm/lattice_set.h"

namespace tempspectest {

template <typename LatSet>
struct rhoImpl {
	using T = typename LatSet::FloatType;

  static inline void apply(const std::vector<T>& cell, T& rho_value) {
    rho_value = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) rho_value += cell[i];
		std::cout << "rho::apply" << std::endl;
  }
};

template <typename LatSet>
struct rho {
	using T = typename LatSet::FloatType;

	static inline T get(const std::vector<T>& cell) {
    T rho_value{};
    rhoImpl<LatSet>::apply(cell, rho_value);
    return rho_value;
  }
	static inline void apply(const std::vector<T>& cell, T& rho_value) {
    rhoImpl<LatSet>::apply(cell, rho_value);
  }
};



// Specialization for D2Q9
template <typename T>
struct rhoImpl<D2Q9<T>> {
  static void apply(const std::vector<T>& cell, T& rho_value) {
    rho_value = T{};
    for (unsigned int i = 0; i < D2Q9<T>::q; ++i) rho_value += cell[i];
    std::cout << "rho_apply::apply D2Q9" << std::endl;
  }
};

template <typename T>
struct rhoImpl<D3Q19<T>> {
  static void apply(const std::vector<T>& cell, T& rho_value) {
    rho_value = T{};
    for (unsigned int i = 0; i < D3Q19<T>::q; ++i) rho_value += cell[i];
    std::cout << "rho_apply::apply D3Q19" << std::endl;
  }
};


}  // namespace tempspectest