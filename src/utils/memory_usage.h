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

// memory_usage.h

#pragma once

#include <cstddef>
#include <iostream>

#if defined(__linux__)
#include <unistd.h>
#include <fstream>
#endif
#if defined(__linux__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

// Resident-set measurements for the complex-geometry size model.  They are
// diagnostics only: a platform without the required call returns 0.
namespace flbmem {

inline std::size_t currentRssBytes() {
#if defined(__linux__)
  std::ifstream statm("/proc/self/statm");
  std::size_t totalPages{};
  std::size_t residentPages{};
  if (!(statm >> totalPages >> residentPages)) return 0;
  const long pageSize = sysconf(_SC_PAGESIZE);
  if (pageSize <= 0) return 0;
  return residentPages * static_cast<std::size_t>(pageSize);
#else
  return 0;
#endif
}

inline std::size_t peakRssBytes() {
#if defined(__linux__) || defined(__APPLE__)
  struct rusage usage {};
  if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#if defined(__APPLE__)
  return static_cast<std::size_t>(usage.ru_maxrss);
#else
  return static_cast<std::size_t>(usage.ru_maxrss) * 1024;
#endif
#else
  return 0;
#endif
}

inline void printMemoryUsage(const char* tag = nullptr) {
  std::cout << "[Memory]";
  if (tag != nullptr && tag[0] != '\0') std::cout << " " << tag;
  std::cout << ": current RSS " << currentRssBytes() << " B, peak RSS "
            << peakRssBytes() << " B" << std::endl;
}

}  // namespace flbmem
