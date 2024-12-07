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

// basic_writer.h

#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>
// std::setw
#include <iomanip>

#include "utils/directories.h"

// return the string of the type for vtk output
template <typename T>
void getVTKTypeString(std::string& type) {
  if constexpr (std::is_same<T, double>::value) {
    type = "Float64";
  } else if constexpr (std::is_same<T, float>::value) {
    type = "Float32";
  } else if constexpr (std::is_same<T, std::uint8_t>::value) {
    type = "UInt8";
  } else if constexpr (std::is_same<T, std::int8_t>::value) {
    type = "Int8";
  } else if constexpr (std::is_same<T, std::uint16_t>::value) {
    type = "UInt16";
  } else if constexpr (std::is_same<T, std::int16_t>::value) {
    type = "Int16";
  } else if constexpr (std::is_same<T, std::uint32_t>::value) {
    type = "UInt32";
  } else if constexpr (std::is_same<T, std::int32_t>::value) {
    type = "Int32";
  } else if constexpr (std::is_same<T, std::uint64_t>::value) {
    type = "UInt64";
  } else if constexpr (std::is_same<T, std::int64_t>::value) {
    type = "Int64";
  } else {
    if constexpr (std::is_enum<T>::value) {
      if constexpr (std::is_same<std::underlying_type_t<T>, std::uint8_t>::value) {
        type = "UInt8";
      }
    } else {
      std::cout << "Error: unsupported type" << std::endl;
      exit(-1);
    }
  }
}

class Printer {
 public:
  static void Print_BigBanner(std::string banner) {
    MPI_RANK(0)
    std::cout << "\n------------------------------------------------\n"
              << "                " << banner << "\n"
              << "------------------------------------------------" << std::endl;
  }
  static void Print_Banner(std::string banner) {
    MPI_RANK(0)
    std::cout << "\n----------------" << banner << "----------------" << std::endl;
  }
  template <typename T>
  static void Print_Res(T Res) {
    MPI_RANK(0)
    std::cout << "Res: " << Res << "  ";
  }
  template <typename T>
  static void Print(std::string title, T value) {
    MPI_RANK(0)
    std::cout << title << ": " << value << "  ";
  }
  template <typename T>
  static void Print(std::string title, T value, std::string title2) {
    MPI_RANK(0)
    std::cout << title << ": " << value << title2 << "  ";
  }
  static void PrintTitle(std::string title) {
    MPI_RANK(0)
    std::cout << "[" << title << "]:\n";
  }
  template <typename T>
  static void Print_SolidFraction(T NucPer) {
    MPI_RANK(0)
    std::cout << "SolidFraction: " << T(100) * NucPer << "%"
              << " ";
  }
  static void Endl() {
    MPI_RANK(0)
    std::cout << std::endl;
  }
  static void Flush() {
    MPI_RANK(0)
    std::cout << std::flush;
  }
};