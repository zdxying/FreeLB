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

#include <fstream>
#include <string>

#include "offLattice/triangle.h"
#include "utils/directories.h"

namespace offlat {

// only for 3D case

template <typename T>
class TriangleSet {
 private:
  // triangles of the whole domain
  std::vector<Triangle<T>> _triangles;
  // triangle-related index stored like "block strcture"
  std::vector<std::vector<TriangleIdx<T>>> _triIdxs;

 public:
  // UINT8[80] – Header
  // UINT32 – Number of triangles
  // foreach triangle
  // REAL32[3] – Normal vector
  // REAL32[3] – Vertex 1
  // REAL32[3] – Vertex 2
  // REAL32[3] – Vertex 3
  // UINT16 – Attribute byte count
  // end
  TriangleSet() = default;
  TriangleSet(std::size_t size) : _triangles(size), _triIdxs(size) {}

  std::vector<Triangle<T>> &getTriangles() { return _triangles; }
  const std::vector<Triangle<T>> &getTriangles() const { return _triangles; }

  std::vector<std::vector<TriangleIdx<T>>> &getTriangleIdxs() { return _triIdxs; }
  const std::vector<std::vector<TriangleIdx<T>>> &getTriangleIdxs() const {
    return _triIdxs;
  }

  void writeBinarySTL(const std::string &fName, T scale = T{1},
    Vector<T, 3> offset = Vector<T, 3>{}) const {
    mpi().barrier();
    std::size_t TotalTriNum{};
    std::size_t TriNum = _triangles.size();
#ifdef MPI_ENABLED
    mpi().reduce(TriNum, TotalTriNum, MPI_SUM);
#else
    TotalTriNum = TriNum;
#endif

    std::string fullName = "./STLoutput/" + fName + ".stl";
    // trim fullName to less than 80 characters
    if (fullName.size() > 79) {
      fullName = fullName.substr(0, 79);
      // set the last character to '\0'
      fullName[79] = '\0';
    }
    BinarySTLHeader(fullName, static_cast<std::uint32_t>(TotalTriNum));
#ifdef MPI_ENABLED
    // serialize the output, each time only one rank writes to the file
    for (int iRank = 0; iRank < mpi().getSize(); ++iRank) {
      mpi().barrier();
      BinarySTLContent(fullName, scale, offset, iRank);
    }
#else
    BinarySTLContent(fullName, scale, offset, 0);
#endif
    
    // output info
    MPI_RANK(0)
    std::cout << "[TriangleSet]: Write " << TotalTriNum << " triangles" << std::endl;
  }

 private:
  void BinarySTLHeader(const std::string &fullName, std::uint32_t TotalTriNum) const {
    MPI_RANK(0)
    DirCreator::Create_Dir("./STLoutput");
    // UINT8[80] – Header
    char buf[80] = {'\0'};
    // copy the string to buf
    std::copy(fullName.begin(), fullName.end(), buf);
    // write header
    std::ofstream f(fullName, std::ios::out | std::ios::binary);
    f.write(buf, 80);
    // write number of triangles
    f.write(reinterpret_cast<const char *>(&TotalTriNum), sizeof(std::uint32_t));

    f.close();
  }

  void BinarySTLContent(
    const std::string &fullName, T scale, Vector<T, 3> offset, int rank = 0) const {
    MPI_RANK(rank)

    const std::uint32_t num = _triangles.size();
    const std::uint16_t abc{};

    std::ofstream f(fullName, std::ios::out | std::ios::app | std::ios::binary);

    for (std::uint32_t i = 0; i < num; ++i) {
      Vector<T, 3> v0 = scale * _triangles[i][0] + offset;
      Vector<T, 3> v1 = scale * _triangles[i][1] + offset;
      Vector<T, 3> v2 = scale * _triangles[i][2] + offset;
      Vector<T, 3> normal = getTriangleNormal(v0, v1, v2);

      float n[3];
      n[0] = static_cast<float>(normal[0]);
      n[1] = static_cast<float>(normal[1]);
      n[2] = static_cast<float>(normal[2]);
      f.write(reinterpret_cast<const char *>(n), sizeof(float) * 3);
      float v[3];
      v[0] = static_cast<float>(v0[0]);
      v[1] = static_cast<float>(v0[1]);
      v[2] = static_cast<float>(v0[2]);
      f.write(reinterpret_cast<const char *>(v), sizeof(float) * 3);
      v[0] = static_cast<float>(v1[0]);
      v[1] = static_cast<float>(v1[1]);
      v[2] = static_cast<float>(v1[2]);
      f.write(reinterpret_cast<const char *>(v), sizeof(float) * 3);
      v[0] = static_cast<float>(v2[0]);
      v[1] = static_cast<float>(v2[1]);
      v[2] = static_cast<float>(v2[2]);
      f.write(reinterpret_cast<const char *>(v), sizeof(float) * 3);
      f.write(reinterpret_cast<const char *>(&abc), sizeof(std::uint16_t));
    }

    f.close();
  }

};


}  // namespace offlat
