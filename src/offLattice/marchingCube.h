/* This file is part of FreeLB, modified from marchingCube.h and marchingCube.hh in palabos
 *
 * // start of the original palabos's copyright notice
 * 
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 * 
 * // end of the original palabos's copyright notice
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

#include <vector>

#include "data_struct/field_struct.h"
#include "offLattice/triangleSet.h"


namespace offlat {

struct MarchingCubeConstants {
  static const int edgeTable[256];
  static const int triTable[256][16];
  static const int edgeNeighb[12][3];
  static const int edgeOrient[12];
};

template <typename T, typename FieldType>
class ScalarFieldIsoSurface {
 public:
  ScalarFieldIsoSurface(const BlockField<FieldType, T, 3>& scalarfield, T isoValue)
      : scalarfield(scalarfield), isoValue(isoValue) {}

  bool isInside(std::size_t id) const { return scalarfield.get(id) < isoValue; }

  Vector<T, 3> getSurface(std::size_t id1, std::size_t id2) const {
    T valp1 = scalarfield.get(id1);
    T valp2 = scalarfield.get(id2);

    Vector<T, 3> p1 = scalarfield.getBlock().getLoc_t(id1);
    Vector<T, 3> p2 = scalarfield.getBlock().getLoc_t(id2);

    if (util::fpequal(isoValue, valp1) || util::fpequal(valp1, valp2)) return p1;
    if (util::fpequal(isoValue, valp2)) return p2;
    T mu = (isoValue - valp1) / (valp2 - valp1);
    return Vector<T, 3>{p1[0] + mu * (p2[0] - p1[0]), p1[1] + mu * (p2[1] - p1[1]),
                        p1[2] + mu * (p2[2] - p1[2])};
  }

  int getNx() const { return scalarfield.getBlock().getNx(); }
  int getNy() const { return scalarfield.getBlock().getNy(); }
  int getNz() const { return scalarfield.getBlock().getNz(); }

  const auto& getBlock() const { return scalarfield.getBlock(); }

 private:
  const BlockField<FieldType, T, 3>& scalarfield;
  T isoValue;
};


template <typename T, typename FieldType>
class MarchingCubeSurface {
 public:
  MarchingCubeSurface(const BlockFieldManager<FieldType, T, 3>& scalarBFM, T isoValue)
      : scalarBFM(scalarBFM), isoValue(isoValue) {}

  using mcc = MarchingCubeConstants;

  void generateIsoSurface(std::vector<Triangle<T>>& tvecs) {
    tvecs.clear();
    for (const BlockField<FieldType, T, 3>& blockfield : scalarBFM.getBlockFields()) {
      // const auto& field = blockfield.getFieldType().getField(0);
      const auto& blockxd = blockfield.getBlock();
      const int overlap = blockxd.getOverlap();
      const int Nx = blockxd.getNx();
      const int Ny = blockxd.getNy();
      const int Nz = blockxd.getNz();

      ScalarFieldIsoSurface<T, FieldType> isoSurface(blockfield, isoValue);

      std::size_t NxNy = Nx * Ny;
      std::size_t id{};
      // since the left-bottom-back corner is used to construct the cube
      // the end of the loop should be Nx - overlap - 1 ...
      // however, this will lead to broken(segmented) stl surface when multi-blocks are used
      for (int k = overlap; k < Nz - overlap; ++k) {
        for (int j = overlap; j < Ny - overlap; ++j) {
          id = k * NxNy + j * Nx + overlap;
          for (int i = overlap; i < Nx - overlap; ++i) {
            polygonize(id, tvecs, isoSurface);
            ++id;
          }
        }
      }
    }
  }

  void polygonize(std::size_t id, std::vector<Triangle<T>>& tvecs,
                  const ScalarFieldIsoSurface<T, FieldType>& isoSurface) {
    static constexpr T epsilon = std::numeric_limits<T>::epsilon() * 1.e4;

    int cubeindex{};
    // vertex list
    std::vector<Vector<T, 3>> vertlist(12);
    marchingCubeImpl(id, tvecs, cubeindex, vertlist, isoSurface);
    // Create the triangle
    for (int i = 0; mcc::triTable[cubeindex][i] != -1; i += 3) {
      int edge1 = mcc::triTable[cubeindex][i];
      int edge2 = mcc::triTable[cubeindex][i + 1];
      int edge3 = mcc::triTable[cubeindex][i + 2];

      Triangle<T> triangle;
      triangle[0] = vertlist[edge1];
      triangle[1] = vertlist[edge2];
      triangle[2] = vertlist[edge3];
      if (getTriangleArea(triangle[0], triangle[1], triangle[2]) > epsilon) {
        tvecs.push_back(triangle);
      }
    }
  }

  void marchingCubeImpl(std::size_t id, std::vector<Triangle<T>>& tvecs, int& cubeindex,
                        std::vector<Vector<T, 3>>& vertlist,
                        const ScalarFieldIsoSurface<T, FieldType>& isoSurface) {
    // p0(iX, iY + 1, iZ);            id + Nx
    // p1(iX + 1, iY + 1, iZ);				id + Nx + 1
    // p2(iX + 1, iY, iZ);						id + 1
    // p3(iX, iY, iZ);                id
    // p4(iX, iY + 1, iZ + 1);				id + Nx * Ny + Nx
    // p5(iX + 1, iY + 1, iZ + 1);		id + Nx * Ny + Nx + 1
    // p6(iX + 1, iY, iZ + 1);				id + Nx * Ny + 1
    // p7(iX, iY, iZ + 1);						id + Nx * Ny

    const int Nx = isoSurface.getNx();
    const int Ny = isoSurface.getNy();
    // const int Nz = isoSurface.getNz();
    const auto& block = isoSurface.getBlock();

    std::size_t id0 = id + Nx;
    std::size_t id1 = id + Nx + 1;
    std::size_t id2 = id + 1;
    std::size_t id3 = id;
    std::size_t id4 = id + Nx * Ny + Nx;
    std::size_t id5 = id + Nx * Ny + Nx + 1;
    std::size_t id6 = id + Nx * Ny + 1;
    std::size_t id7 = id + Nx * Ny;

    // IdxLoc
    Vector<int, 3> p3i = block.getLoc(id);
    Vector<int, 3> p0i = p3i + Vector<int, 3>{0, 1, 0};
    Vector<int, 3> p1i = p3i + Vector<int, 3>{1, 1, 0};
    Vector<int, 3> p2i = p3i + Vector<int, 3>{1, 0, 0};
    Vector<int, 3> p4i = p3i + Vector<int, 3>{0, 1, 1};
    Vector<int, 3> p5i = p3i + Vector<int, 3>{1, 1, 1};
    Vector<int, 3> p6i = p3i + Vector<int, 3>{1, 0, 1};
    Vector<int, 3> p7i = p3i + Vector<int, 3>{0, 0, 1};

    Vector<T, 3> p3 = block.getLoc_t(p3i);
    Vector<T, 3> p0 = block.getLoc_t(p0i);
    Vector<T, 3> p1 = block.getLoc_t(p1i);
    Vector<T, 3> p2 = block.getLoc_t(p2i);
    Vector<T, 3> p4 = block.getLoc_t(p4i);
    Vector<T, 3> p5 = block.getLoc_t(p5i);
    Vector<T, 3> p6 = block.getLoc_t(p6i);
    Vector<T, 3> p7 = block.getLoc_t(p7i);

    if (isoSurface.isInside(id0)) cubeindex |= 1;    // Point 0
    if (isoSurface.isInside(id1)) cubeindex |= 2;    // Point 1
    if (isoSurface.isInside(id2)) cubeindex |= 4;    // Point 2
    if (isoSurface.isInside(id3)) cubeindex |= 8;    // Point 3
    if (isoSurface.isInside(id4)) cubeindex |= 16;   // Point 4
    if (isoSurface.isInside(id5)) cubeindex |= 32;   // Point 5
    if (isoSurface.isInside(id6)) cubeindex |= 64;   // Point 6
    if (isoSurface.isInside(id7)) cubeindex |= 128;  // Point 7

    vertlist.resize(12);
    // Cube is entirely in/out of the surface
    if (mcc::edgeTable[cubeindex] == 0) return;

    // Find the vertices where the surface intersects the cube
    if (mcc::edgeTable[cubeindex] & 1) {
      vertlist[0] = isoSurface.getSurface(id0, id1);  // x-edge of y-neighbor.
      removeFromVertex(p0, p1, vertlist[0]);
    }
    if (mcc::edgeTable[cubeindex] & 2) {
      vertlist[1] = isoSurface.getSurface(id1, id2);  // y-edge of x-neighbor.
      removeFromVertex(p1, p2, vertlist[1]);
    }
    if (mcc::edgeTable[cubeindex] & 4) {
      vertlist[2] = isoSurface.getSurface(id2, id3);  // x-edge of current cell.
      removeFromVertex(p2, p3, vertlist[2]);
    }
    if (mcc::edgeTable[cubeindex] & 8) {
      vertlist[3] = isoSurface.getSurface(id3, id0);  // y-edge of current cell.
      removeFromVertex(p3, p0, vertlist[3]);
    }
    if (mcc::edgeTable[cubeindex] & 16) {
      vertlist[4] = isoSurface.getSurface(id4, id5);  // x-edge of y-z-neighbor.
      removeFromVertex(p4, p5, vertlist[4]);
    }
    if (mcc::edgeTable[cubeindex] & 32) {
      vertlist[5] = isoSurface.getSurface(id5, id6);  // y-edge of x-z-neighbor.
      removeFromVertex(p5, p6, vertlist[5]);
    }
    if (mcc::edgeTable[cubeindex] & 64) {
      vertlist[6] = isoSurface.getSurface(id6, id7);  // x-edge of z-neighbor.
      removeFromVertex(p6, p7, vertlist[6]);
    }
    if (mcc::edgeTable[cubeindex] & 128) {
      vertlist[7] = isoSurface.getSurface(id7, id4);  // y-edge of z-neighbor.
      removeFromVertex(p7, p4, vertlist[7]);
    }
    if (mcc::edgeTable[cubeindex] & 256) {
      vertlist[8] = isoSurface.getSurface(id0, id4);  // z-edge of y-neighbor.
      removeFromVertex(p0, p4, vertlist[8]);
    }
    if (mcc::edgeTable[cubeindex] & 512) {
      vertlist[9] = isoSurface.getSurface(id1, id5);  // z-edge of x-y-neighbor.
      removeFromVertex(p1, p5, vertlist[9]);
    }
    if (mcc::edgeTable[cubeindex] & 1024) {
      vertlist[10] = isoSurface.getSurface(id2, id6);  // z-edge of x-neighbor.
      removeFromVertex(p2, p6, vertlist[10]);
    }
    if (mcc::edgeTable[cubeindex] & 2048) {
      vertlist[11] = isoSurface.getSurface(id3, id7);  // z-edge of current cell.
      removeFromVertex(p3, p7, vertlist[11]);
    }
  }

  void removeFromVertex(const Vector<T, 3>& p0, const Vector<T, 3>& p1,
                        Vector<T, 3>& intersection) {
    static const T triangleEpsilon = 1.e-3;
    static const T triangleEpsilonSqr = triangleEpsilon * triangleEpsilon;
    if ((p0 - intersection).getnorm2() < triangleEpsilonSqr) {
      intersection = p0 + triangleEpsilon * (p1 - p0);
    } else if ((p1 - intersection).getnorm2() < triangleEpsilonSqr) {
      intersection = p1 - triangleEpsilon * (p1 - p0);
    }
  }

 private:
  const BlockFieldManager<FieldType, T, 3>& scalarBFM;
  T isoValue;
};

}  // namespace offlat
