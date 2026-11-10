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

// mpi_macros.h
//
// Leaf header with the rank-guard logging macros.  The macros expand to calls
// of the flb:: helper functions declared here (defined inline in
// parallel/mpi_manager.h), so a translation unit that only *logs* on rank 0
// (utils/timer.h, utils/directories.h, lbm/unit_converter.h) can *compile*
// without pulling in the full MpiManager definition.  Code that calls mpi()
// methods directly still needs parallel/mpi_manager.h, and every TU that
// expands these macros must link one TU that includes parallel/mpi_manager.hh
// (which defines mpi() and therefore the flb:: helpers).

#pragma once

namespace flb {
// true if the calling process has rank r (MPI builds only)
bool is_rank(int r);
// number of MPI processes
int size();
// MPI_DEBUG_WAIT helper: spin-wait on rank 0, then barrier
void debug_wait();
}  // namespace flb

#ifdef MPI_ENABLED

// guard: run only on the given rank (usually 0), by early return
#define MPI_RANK(x)                     \
  if (!flb::is_rank(x)) {               \
    return;                             \
  }
// followed by a enclosing brace
#define IF_MPI_RANK(x) if (flb::is_rank(x))

#ifdef MPI_DEBUG
#define MPI_DEBUG_WAIT     \
  do {                     \
    flb::debug_wait();     \
  } while (0)
#else
#define MPI_DEBUG_WAIT
#endif  // #ifdef MPI_DEBUG

#else  // #ifdef MPI_ENABLED

#define MPI_RANK(x)
#define IF_MPI_RANK(x)
#define MPI_DEBUG_WAIT

#endif  // #ifdef MPI_ENABLED
