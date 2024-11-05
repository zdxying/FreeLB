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
  std::vector<Triangle<T>> _triangles;

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
	TriangleSet(std::size_t size) : _triangles(size) {}
	
  std::vector<Triangle<T>> &getTriangles() { return _triangles; }

  void writeBinarySTL(std::string fName, T scale = T{1},
                      Vector<T, 3> offset = Vector<T, 3>{}) const {
    writeBinarySTL(fName, _triangles, scale, offset);
  }

  void writeBinarySTL(std::string fName, const std::vector<Triangle<T>> &triangles,
                      T scale = T{1}, Vector<T, 3> offset = Vector<T, 3>{}) const {
    DirCreator::Create_Dir("./STLoutput");
    std::string fullName = "./STLoutput/" + fName + ".stl";
    std::ofstream f(fullName, std::ios::binary);

    const std::int32_t num = triangles.size();
    const std::uint16_t abc = 0;
    char buf[80] = {'\0'};

    // trim fullName to less than 80 characters
    if (fullName.size() > 79) {
      fullName = fullName.substr(0, 79);
      // set the last character to '\0'
      fullName[79] = '\0';
    }
    // copy the string to buf
    std::copy(fullName.begin(), fullName.end(), buf);

    // write header
    f.write(buf, 80);
    // write number of triangles
    f.write(reinterpret_cast<const char *>(&num), sizeof(std::int32_t));
    for (std::int32_t i = 0; i < num; ++i) {
      Vector<T, 3> v0 = scale * triangles[i][0] + offset;
      Vector<T, 3> v1 = scale * triangles[i][1] + offset;
      Vector<T, 3> v2 = scale * triangles[i][2] + offset;
      Vector<T, 3> normal = getTriangleNormal(v0, v1, v2);

      float n[3];
      n[0] = normal[0];
      n[1] = normal[1];
      n[2] = normal[2];
      f.write(reinterpret_cast<const char *>(n), sizeof(float) * 3);
      float v[3];
      v[0] = v0[0];
      v[1] = v0[1];
      v[2] = v0[2];
      f.write(reinterpret_cast<const char *>(v), sizeof(float) * 3);
      v[0] = v1[0];
      v[1] = v1[1];
      v[2] = v1[2];
      f.write(reinterpret_cast<const char *>(v), sizeof(float) * 3);
      v[0] = v2[0];
      v[1] = v2[1];
      v[2] = v2[2];
      f.write(reinterpret_cast<const char *>(v), sizeof(float) * 3);
      f.write(reinterpret_cast<const char *>(&abc), sizeof(std::uint16_t));
    }
		f.close();
		// out put info
		std::cout << "[TriangleSet]: Write " << num << " triangles" << std::endl;
  }
};


}  // namespace offlat
