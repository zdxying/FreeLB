/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 */

// io/util: dependency-free IO helpers, merged from the former io/ini and io/dir
// targets.  Both exercise helpers that the vtk/vtm/vtu writers rely on
// (iniReader for parameter files, DirCreator for output directories), so they
// share one binary.

#include "freelb.h"
#include "freelb.hh"
#include "io/ini_reader.h"

#include <filesystem>

#include "../../test_util.h"

using T = float;

int main() {
  // ------------------------------- iniReader -------------------------------
  iniReader reader("test.ini");

  CHECK_EQ(reader.getValue<int>("int", "intx"), 25);
  CHECK_EQ(reader.getValue<std::string>("string", "stringx"), std::string("./output"));
  CHECK_NEAR(reader.getValue<T>("real", "realx"), 0.04, 1e-9);
  CHECK_NEAR(reader.getValue<T>("real", "realy"), 1.5e-4, 1e-9);

  std::vector<int> vec;
  reader.getVector<int>("vector", vec);
  CHECK_EQ(vec.size(), std::size_t(3));
  if (vec.size() == 3) {
    CHECK_EQ(vec[0], 1);
    CHECK_EQ(vec[1], 2);
    CHECK_EQ(vec[2], 3);
  }

  // ------------------------------- DirCreator ------------------------------
  // Make the section deterministic: a previous run without `make clean` may
  // have left these directories behind, which would make the error-path check
  // below (which needs ./test to be absent) fail spuriously.
  for (const char* d : {"./test", "./test0", "./test1", "./test2"}) {
    std::filesystem::remove_all(d);
  }

  // parent ./test does not exist yet -> must fail and create nothing
  DirCreator::Create_Dir("./test/error");
  CHECK(!std::filesystem::exists("./test/error"));

  // multiple sibling directories in one call
  DirCreator::Create_Dir("./test0/", "./test1/", "./test2/");
  CHECK(std::filesystem::exists("./test0/"));
  CHECK(std::filesystem::exists("./test1/"));
  CHECK(std::filesystem::exists("./test2/"));

  // nested directories, created in the given order
  DirCreator::Create_Dir("./test/", "./test/sub/", "./test/sub/sub0/", "./test/sub1/");
  CHECK(std::filesystem::exists("./test/"));
  CHECK(std::filesystem::exists("./test/sub/"));
  CHECK(std::filesystem::exists("./test/sub/sub0/"));
  CHECK(std::filesystem::exists("./test/sub1/"));

  return flbtest::summary("io/util");
}
