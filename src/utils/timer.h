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

#pragma once

#include <chrono>
#include <cstdint>
#include <iostream>

#include "head.h"

// timer and counters

struct Counter {
  std::int64_t count;
  Counter() { count = 0; }
  ~Counter() = default;
  // reset count
  void reset() { count = 0; }
  // define ++ operator, use: Counter c; ++c; or c++;
  void operator++() { count++; }
  // get count
  // define () operator, use: Counter c; c(); or c.operator()();
  std::int64_t operator()() { return count; }
};

struct Timer : public Counter {
  std::chrono::high_resolution_clock::time_point START;
  std::chrono::high_resolution_clock::time_point END;

  Timer() { START_TIMER(); }
  ~Timer() = default;

  void START_TIMER() { START = std::chrono::high_resolution_clock::now(); }
  void END_TIMER() { END = std::chrono::high_resolution_clock::now(); }

  std::chrono::milliseconds GetDuration() {
    END_TIMER();
    return std::chrono::duration_cast<std::chrono::milliseconds>(END - START);
  }
  // get time elapsed in milliseconds(ms)
  std::int64_t GetDurationCount() {
    END_TIMER();
    return std::chrono::duration_cast<std::chrono::milliseconds>(END - START).count();
  }
  // get time elapsed in milliseconds(ms) but no END_TIMER()
  std::int64_t GetDurationCount_Only() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(END - START).count();
  }
  // get time elapsed in seconds(s)
  double GetTimeElapsed() { return GetDurationCount() / double(1000); }
  // get time elapsed in seconds(s) but no END_TIMER()
  double GetTimeElapsed_Only() { return GetDurationCount_Only() / double(1000); }

  void Print_MainLoopPerformance(std::size_t ni, std::size_t nj, std::size_t nk = 1) {
    Print_MainLoopPerformance(ni * nj * nk);
  }
  void Print_MainLoopPerformance(std::size_t n) {
    MPI_RANK(0)
    std::int64_t count = this->operator()();  // this->operator()() or (*this)()
    double MLUPs = static_cast<double>(count) * static_cast<double>(n) /
                   GetDurationCount() / double(1000);
    std::cout << "[Main_Loop Performance]:"
              << "\n"
              << "Time Elapsed: " << GetTimeElapsed() << " s"
              << "\n"
              << "Total Step:   " << count << "\n"
              << "Ave_MLUPs:    " << MLUPs << "\n";
#ifdef _OPENMP
    std::cout << "MLUPs/Thread: " << MLUPs / Thread_Num << std::endl;
#endif
  }
  void Print_InnerLoopPerformance(std::size_t steps, std::size_t ni, std::size_t nj,
                                  std::size_t nk = 1) {
    Print_InnerLoopPerformance(steps, ni * nj * nk);
  }
  void Print_InnerLoopPerformance(std::size_t steps, std::size_t n) {
    MPI_RANK(0)
    std::int64_t count = this->operator()();  // this->operator()() or (*this)()
    double MLUPs = static_cast<double>(steps) * static_cast<double>(n) /
                   GetDurationCount() / double(1000);
    std::cout << "\n[Step: " << count << "]  "
              << "MLUPs: " << MLUPs << "  ";
    // this->reset();
    START_TIMER();
  }
};
