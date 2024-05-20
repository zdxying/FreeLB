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

// cell_dynamics.h
#pragma once

#include "lbm/collision.h"
#include "lbm/moment.h"

// a template for cell dynamics
// you could use it like:
// using MyDynamics = CellDynamics<task1, task2, task3, ...>;
template <typename... Tasks>
struct CellDynamics {
  template <typename CELL, typename... Args>
  static void Apply(CELL& cell, Args... args) {
    (applyHelper<Tasks>(cell, args...), ...);
  }

  template <typename Task, typename CELL, typename FirstArg, typename... Args>
  static auto applyHelper(CELL& cell, FirstArg firstArg, Args... args)
    -> std::enable_if_t<
      std::is_invocable_v<decltype(&Task::apply), CELL&, FirstArg, Args...>, void> {
    Task::apply(cell, firstArg, args...);
  }
  template <typename Task, typename CELL, typename FirstArg, typename... Args>
  static auto applyHelper(CELL& cell, FirstArg firstArg, Args... args)
    -> std::enable_if_t<
      !std::is_invocable_v<decltype(&Task::apply), CELL&, FirstArg, Args...>, void> {
    applyHelper<Task>(cell, args...);
  }
  template <typename Task, typename CELL>
  static auto applyHelper(CELL& cell)
    -> std::enable_if_t<std::is_invocable_v<decltype(&Task::apply), CELL&>, void> {
    Task::apply(cell);
  }
};

namespace dynamics {
  
}  // namespace dynamics