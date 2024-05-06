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

// plot_writer.h

#pragma once

#include "io/basic_writer.h"

template <typename T>
class FLBplot {
 private:
  // x-axis: step
  // y-axis: data
  std::string name;
  std::string xvar;
  std::string yvar;
  std::string work_dir;

 public:
  FLBplot(std::string work_dir_, std::string name = "/FLBplot",
          std::string xvar = "step", std::string yvar = "data")
      : work_dir(work_dir_), name(name), xvar(xvar), yvar(yvar) {
    plotHeader();
  }
  void plotHeader() {
    std::ofstream write;
    write.open(work_dir + name + ".dat");
    write << xvar << "\t" << yvar << "\n";
  }
  // each time call this function, a step-data pair is written
  template <typename X = int, typename U = T>
  void write(X xvalue, U yvalue) {
    std::ofstream write;
    write.open(work_dir + name + ".dat", std::ios::app);
    write << xvalue << "\t" << yvalue << "\n";
  }
  template <typename U = T>
  void write(Counter &counter, U yvalue) {
    write<int, U>(counter(), yvalue);
  }
};