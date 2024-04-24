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

// force.h
#pragma once

#include "lbm/moment.h"
#include "utils/util.h"

namespace force {

template <typename T, typename LatSet>
struct Force {
  // calculate discrete force
  // c * (c * u) * LatSet::InvCs4
  // Vector<T, LatSet::d> ccu = (LatSet::c[i] * u * LatSet::InvCs4) *
  // LatSet::c[i]; (c - u) * LatSet::InvCs2 Vector<T, LatSet::d> c_u =
  // (LatSet::c[i] - u) * LatSet::InvCs2;
  static std::array<T, LatSet::q> getForcePop(const Vector<T, LatSet::d> &u,
                                              const Vector<T, LatSet::d> &F) {
    std::array<T, LatSet::q> Fi{};
    ComputeForcePop(Fi, u, F);
    return Fi;
  }

  static void ComputeForcePop(std::array<T, LatSet::q> &Fi, const Vector<T, LatSet::d> &u,
                              const Vector<T, LatSet::d> &F) {
    for (int i = 0; i < LatSet::q; ++i) {
      Fi[i] = LatSet::w[i] * F *
              ((LatSet::c[i] - u) * LatSet::InvCs2 +
               (LatSet::c[i] * u * LatSet::InvCs4) * LatSet::c[i]);
    }
  }

  static void ComputeForcePop(std::array<T, LatSet::q> &Fi, const Vector<T, LatSet::d> &u,
                              const T F) {
    constexpr unsigned int d = LatSet::d - 1;
    for (int i = 0; i < LatSet::q; ++i) {
      const T v1 = (LatSet::c[i][d] - u[d]) * LatSet::InvCs2;
      const T v2 = (LatSet::c[i] * u * LatSet::InvCs4) * LatSet::c[i][d];
      Fi[i] = LatSet::w[i] * F * (v1 + v2);
    }
  }
};

}  // namespace force