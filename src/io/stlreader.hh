/* This file is part of FreeLB, modified from OpenLB's stlReader.hh, with the following copyright notice:
 *
 * // start of the original OpenLB's copyright notice
 * 
 * This file is part of the OpenLB library
 *
 *  Copyright (C) 2010-2015 Thomas Henn, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 * 
 * // end of the original OpenLB's copyright notice
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

// stlreader.hh

#pragma once

#include <fstream>
#include <iostream>
// #include <stdexcept>

#include "data_struct/octree.hh"
#include "io/stlreader.h"
#include "utils/directories.h"

template <typename T>
void Triangle<T>::Init() {
  Vector<T, 3> A = vertex[0];
  Vector<T, 3> B = vertex[1];
  Vector<T, 3> C = vertex[2];
  Vector<T, 3> b, c;
  T bb = 0, cc = 0, bc = 0;

  // for (int i = 0; i < 3; i++) {
  //   b[i] = B[i] - A[i];
  //   c[i] = C[i] - A[i];
  //   bb += b[i] * b[i];
  //   cc += c[i] * c[i];
  //   bc += b[i] * c[i];
  // }
  // normal[0] = b[1] * c[2] - b[2] * c[1];
  // normal[1] = b[2] * c[0] - b[0] * c[2];
  // normal[2] = b[0] * c[1] - b[1] * c[0];
  // T norm = std::sqrt(std::pow(normal[0], 2) + std::pow(normal[1], 2) + std::pow(normal[2], 2));
  // normal[0] /= norm;
  // normal[1] /= norm;
  // normal[2] /= norm;

  b = B - A;
  c = C - A;
  bb = b.getnorm2();
  cc = c.getnorm2();
  bc = b * c;

  normal = CrossProduct(b, c);
  normal.normalize();

  T D = 1.0 / (cc * bb - bc * bc);
  T bbD = bb * D;
  T bcD = bc * D;
  T ccD = cc * D;

  kBeta = 0.;
  kGamma = 0.;
  d = 0.;

  for (int i = 0; i < 3; i++) {
    uBeta[i] = b[i] * ccD - c[i] * bcD;
    uGamma[i] = c[i] * bbD - b[i] * bcD;
    kBeta -= A[i] * uBeta[i];
    kGamma -= A[i] * uGamma[i];
    d += A[i] * normal[i];
  }
}

template <typename T>
Vector<T, 3> Triangle<T>::getCenter() {
  Vector<T, 3> center(T(0));

  center[0] = (vertex[0][0] + vertex[1][0] + vertex[2][0]) / 3.;
  center[1] = (vertex[0][1] + vertex[1][1] + vertex[2][1]) / 3.;
  center[2] = (vertex[0][2] + vertex[1][2] + vertex[2][2]) / 3.;

  return center;
}

template <typename T>
std::vector<T> Triangle<T>::getE0() {
  Vector<T, 3> vec;
  vec[0] = vertex[0][0] - vertex[1][0];
  vec[1] = vertex[0][1] - vertex[1][1];
  vec[2] = vertex[0][2] - vertex[1][2];
  return vec;
}

template <typename T>
std::vector<T> Triangle<T>::getE1() {
  Vector<T, 3> vec;
  vec[0] = vertex[0][0] - vertex[2][0];
  vec[1] = vertex[0][1] - vertex[2][1];
  vec[2] = vertex[0][2] - vertex[2][2];
  return vec;
}

template <typename T>
bool Triangle<T>::IsInside(const Vector<T, 3> &pt) const {
  constexpr T epsilon = std::numeric_limits<T>::epsilon() * T(10);

  const T beta = pt * uBeta + kBeta;
  const T gamma = pt * uGamma + kGamma;

  // check if approximately equal
  Vector<T, 3> vec = pt - (vertex[0] + beta * (vertex[1] - vertex[0]) +
                           gamma * (vertex[2] - vertex[0]));
  if (util::nearZero(vec.getnorm(), epsilon)) {
    const T alpha = T(1) - beta - gamma;
    return (beta >= T(0) || util::nearZero(beta, epsilon)) &&
           (gamma >= T(0) || util::nearZero(gamma, epsilon)) &&
           (alpha >= T(0) || util::nearZero(alpha, epsilon));
  }
  return false;
}

template <typename T>
bool Triangle<T>::testRayIntersect(const Vector<T, 3> &pt,
                                   const Vector<T, 3> &dir, Vector<T, 3> &q,
                                   T &alpha, const T &rad, bool print) {
  T rn = 0.;
  Vector<T, 3> testPt = pt + rad * normal;

  rn = dir * normal;

  // check if the ray is parallel to the triangle
  // Schnitttest Flugrichtung -> Ebene
  if (std::fabs(rn) < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  // calculates the distance(along the ray) from the starting point of the ray
  // to the plane of the triangle
  alpha =
      d - testPt[0] * normal[0] - testPt[1] * normal[1] - testPt[2] * normal[2];
  //  alpha -= testPt[i] * normal[i];
  alpha /= rn;

  // Abstand Partikel Ebene
  // If alpha is negative, the triangle is behind the starting point of the ray,
  // which does not intersect the triangle
  if (alpha < -std::numeric_limits<T>::epsilon()) {
    return false;
  }
  // If alpha is non-negative, calculates the intersection point q
  // moving along the ray from the starting point by a distance alpha
  for (int i = 0; i < 3; i++) {
    q[i] = testPt[i] + alpha * dir[i];
  }
  // checks if the intersection point q is inside the triangle.
  T beta = kBeta;
  for (int i = 0; i < 3; i++) {
    beta += uBeta[i] * q[i];
  }

  // Schnittpunkt q in der Ebene?
  if (beta < -std::numeric_limits<T>::epsilon()) {
    return false;
  }
  T gamma = kGamma;
  for (int i = 0; i < 3; i++) {
    gamma += uGamma[i] * q[i];
  }
  if (gamma < -std::numeric_limits<T>::epsilon()) {
    return false;
  }
  if (1. - beta - gamma < -std::numeric_limits<T>::epsilon()) {
    return false;
  }
  return true;
}

template <typename T>
Vector<T, 3> Triangle<T>::closestPtPointTriangle(const Vector<T, 3> &pt) const {
  const T nEps = -std::numeric_limits<T>::epsilon();
  const T Eps = std::numeric_limits<T>::epsilon();

  Vector<T, 3> ab = vertex[1] - vertex[0];
  Vector<T, 3> ac = vertex[2] - vertex[0];
  Vector<T, 3> bc = vertex[2] - vertex[1];

  T snom = (pt - vertex[0]) * ab;
  T sdenom = (pt - vertex[1]) * (vertex[0] - vertex[1]);

  T tnom = (pt - vertex[0]) * ac;
  T tdenom = (pt - vertex[2]) * (vertex[0] - vertex[2]);

  if (snom < nEps && tnom < nEps) {
    return vertex[0];
  }

  T unom = (pt - vertex[1]) * bc;
  T udenom = (pt - vertex[2]) * (vertex[1] - vertex[2]);

  if (sdenom < nEps && unom < nEps) {
    return vertex[1];
  }
  if (tdenom < nEps && udenom < nEps) {
    return vertex[2];
  }

  T vc = normal * CrossProduct3(vertex[0] - pt, vertex[1] - pt);

  if (vc < nEps && snom > Eps && sdenom > Eps) {
    return vertex[0] + snom / (snom + sdenom) * ab;
  }

  T va = normal * CrossProduct3(vertex[1] - pt, vertex[2] - pt);

  if (va < nEps && unom > Eps && udenom > Eps) {
    return vertex[1] + unom / (unom + udenom) * bc;
  }

  T vb = normal * CrossProduct3(vertex[2] - pt, vertex[0] - pt);

  if (vb < nEps && tnom > Eps && tdenom > Eps) {
    return vertex[0] + tnom / (tnom + tdenom) * ac;
  }

  T u = va / (va + vb + vc);
  T v = vb / (va + vb + vc);
  T w = 1. - u - v;

  return u * vertex[0] + v * vertex[1] + w * vertex[2];
}

// The ASCII STL file:

// solid _name
// facet normal ni nj nk
//    outer loop
//        vertex v1x v1y v1z
//        vertex v2x v2y v2z
//        vertex v3x v3y v3z
//    endloop
// endfacet
// endsolid name

// The binary file:

// UINT8[80] – Header
// UINT32 – Number of triangles
// foreach triangle
// REAL32[3] – Normal vector
// REAL32[3] – Vertex 1
// REAL32[3] – Vertex 2
// REAL32[3] – Vertex 3
// UINT16 – Attribute byte count
// end

template <typename T>
StlMesh<T>::StlMesh(std::string filename, T stlSize)
    : _fName(filename), _min(), _max(), _maxDist2(0), _stlSize(stlSize) {
  // try to open as ASCII STL
  std::ifstream f(_fName.c_str(), std::ios::in);
  Triangles.reserve(10000);
  //
  if (!f.good()) {
    // throw std::runtime_error("STL File not valid.");
    std::cout << "STL File not valid." << std::endl;
    exit(-1);
  }

  char buf[6];
  // the last set to 0
  buf[5] = 0;
  f.read(buf, 5);
  if (std::string(buf) == "solid") {
    // ASCII STL
    // set the file pointer to the very beginning
    f.seekg(0, std::ios::beg);
    if (f.good()) {
      std::string s0, s1;
      int i = 0;
      while (!f.eof()) {
        f >> s0;
        if (s0 == "facet") {
          Triangle<T> tri;
          f >> s1 >> tri.normal[0] >> tri.normal[1] >> tri.normal[2];
          f >> s0 >> s1;
          f >> s0 >> tri.vertex[0][0] >> tri.vertex[0][1] >> tri.vertex[0][2];
          f >> s0 >> tri.vertex[1][0] >> tri.vertex[1][1] >> tri.vertex[1][2];
          f >> s0 >> tri.vertex[2][0] >> tri.vertex[2][1] >> tri.vertex[2][2];
          // skip over the "endfacet" and "endloop" keyword
          f >> s0;
          f >> s0;
          for (int j = 0; j < 3; j++) {
            tri.vertex[0][j] *= stlSize;
            tri.vertex[1][j] *= stlSize;
            tri.vertex[2][j] *= stlSize;
          }
          // get the min and max vertex
          if (i == 0) {
            Getminmax0(tri);
          } else {
            Getminmax(tri);
          }

          i++;
          tri.Init();
          Triangles.push_back(tri);

          _maxDist2 =
              std::max(_maxDist2, GetDist2(tri.vertex[0], tri.vertex[1]));
          _maxDist2 =
              std::max(_maxDist2, GetDist2(tri.vertex[0], tri.vertex[2]));
          _maxDist2 =
              std::max(_maxDist2, GetDist2(tri.vertex[1], tri.vertex[2]));
        } else if (s0 == "endsolid") {
          break;
        }
      }
    }
  } else {
    f.close();
    // Binary STL
    f.open(_fName.c_str(), std::ios::in | std::ios::binary);
    // read the header
    char header[80];
    f.read(header, 80);

    if (!f.good()) {
      // throw std::runtime_error("STL File not valid.");
      std::cout << "STL File not valid." << std::endl;
      exit(-1);
    }

    header[79] = 0;
    int32_t nFacets;
    f.read(reinterpret_cast<char *>(&nFacets), sizeof(int32_t));

    if (!f.good()) {
      // throw std::runtime_error("STL File not valid.");
      std::cout << "STL File not valid." << std::endl;
      exit(-1);
    }

    float v[12];
    std::uint16_t uint16;

    for (int32_t i = 0; i < nFacets; ++i) {
      for (unsigned int j = 0; j < 12; ++j) {
        f.read(reinterpret_cast<char *>(&v[j]), sizeof(float));
      }
      f.read(reinterpret_cast<char *>(&uint16), sizeof(std::uint16_t));
      Triangle<T> tri;
      tri.normal[0] = v[0];
      tri.normal[1] = v[1];
      tri.normal[2] = v[2];
      tri.vertex[0][0] = v[3];
      tri.vertex[0][1] = v[4];
      tri.vertex[0][2] = v[5];
      tri.vertex[1][0] = v[6];
      tri.vertex[1][1] = v[7];
      tri.vertex[1][2] = v[8];
      tri.vertex[2][0] = v[9];
      tri.vertex[2][1] = v[10];
      tri.vertex[2][2] = v[11];

      for (int k = 0; k < 3; k++) {
        tri.vertex[0][k] *= stlSize;
        tri.vertex[1][k] *= stlSize;
        tri.vertex[2][k] *= stlSize;
      }
      // get the min and max vertex
      if (i == 0) {
        Getminmax0(tri);
      } else {
        Getminmax(tri);
      }
      tri.Init();
      Triangles.push_back(tri);

      _maxDist2 = std::max(_maxDist2, GetDist2(tri.vertex[0], tri.vertex[1]));
      _maxDist2 = std::max(_maxDist2, GetDist2(tri.vertex[0], tri.vertex[2]));
      _maxDist2 = std::max(_maxDist2, GetDist2(tri.vertex[1], tri.vertex[2]));
    }
  }
  f.close();
}

template <typename T>
inline void StlMesh<T>::Getminmax0(Triangle<T> &tri) {
  _min[0] = tri.vertex[0][0];
  _min[1] = tri.vertex[0][1];
  _min[2] = tri.vertex[0][2];

  _max[0] = tri.vertex[0][0];
  _max[1] = tri.vertex[0][1];
  _max[2] = tri.vertex[0][2];

  _min[0] = std::min(_min[0], tri.vertex[1][0]);
  _min[1] = std::min(_min[1], tri.vertex[1][1]);
  _min[2] = std::min(_min[2], tri.vertex[1][2]);

  _max[0] = std::max(_max[0], tri.vertex[1][0]);
  _max[1] = std::max(_max[1], tri.vertex[1][1]);
  _max[2] = std::max(_max[2], tri.vertex[1][2]);

  _min[0] = std::min(_min[0], tri.vertex[2][0]);
  _min[1] = std::min(_min[1], tri.vertex[2][1]);
  _min[2] = std::min(_min[2], tri.vertex[2][2]);

  _max[0] = std::max(_max[0], tri.vertex[2][0]);
  _max[1] = std::max(_max[1], tri.vertex[2][1]);
  _max[2] = std::max(_max[2], tri.vertex[2][2]);
}

template <typename T>
inline void StlMesh<T>::Getminmax(Triangle<T> &tri) {
  _min[0] = std::min(_min[0], tri.vertex[0][0]);
  _min[1] = std::min(_min[1], tri.vertex[0][1]);
  _min[2] = std::min(_min[2], tri.vertex[0][2]);

  _max[0] = std::max(_max[0], tri.vertex[0][0]);
  _max[1] = std::max(_max[1], tri.vertex[0][1]);
  _max[2] = std::max(_max[2], tri.vertex[0][2]);

  _min[0] = std::min(_min[0], tri.vertex[1][0]);
  _min[1] = std::min(_min[1], tri.vertex[1][1]);
  _min[2] = std::min(_min[2], tri.vertex[1][2]);

  _max[0] = std::max(_max[0], tri.vertex[1][0]);
  _max[1] = std::max(_max[1], tri.vertex[1][1]);
  _max[2] = std::max(_max[2], tri.vertex[1][2]);

  _min[0] = std::min(_min[0], tri.vertex[2][0]);
  _min[1] = std::min(_min[1], tri.vertex[2][1]);
  _min[2] = std::min(_min[2], tri.vertex[2][2]);

  _max[0] = std::max(_max[0], tri.vertex[2][0]);
  _max[1] = std::max(_max[1], tri.vertex[2][1]);
  _max[2] = std::max(_max[2], tri.vertex[2][2]);
}

template <typename T>
void StlMesh<T>::print(bool full) {
  MPI_RANK(0)
  std::cout << "[StlMesh " << _fName << "]: " << std::endl;
  if (full || Triangles.size() < 20) {
    int i = 1;
    std::cout << "---Triangles---" << std::endl;
    // iterate over all triangles
    for (auto &tri : Triangles) {
      std::cout << i << ": " << tri.vertex[0][0] << " " << tri.vertex[0][1]
                << " " << tri.vertex[0][2] << " | " << tri.vertex[1][0] << " "
                << tri.vertex[1][1] << " " << tri.vertex[1][2] << " | "
                << tri.vertex[2][0] << " " << tri.vertex[2][1] << " "
                << tri.vertex[2][2] << std::endl;
      ++i;
    }
    std::cout << "----end of Triangles----" << std::endl;
  }
  std::cout << "SizeofTriangles = " << Triangles.size()
            << "; maxDist^2 of triangles = " << _maxDist2 << std::endl;
  std::cout << "minPoint = (" << getMin()[0] << "," << getMin()[1] << ","
            << getMin()[2] << "); "
            << "maxPoint = (" << getMax()[0] << "," << getMax()[1] << ","
            << getMax()[2] << ")" << std::endl;
}

template <typename T>
void StlMesh<T>::write(std::string fName) {
  DirCreator::Create_Dir("./STLoutput");
  std::string fullName = "./STLoutput/" + fName + ".stl";
  std::ofstream f(fullName.c_str());
  f << "solid ascii " << fullName << "\n";
  for (unsigned int i = 0; i < Triangles.size(); ++i) {
    Triangle<T> &tri = Triangles[i];
    f << "facet normal " << tri.normal[0] << " " << tri.normal[1] << " "
      << tri.normal[2] << "\n";
    f << "    outer loop\n";
    f << "        vertex " << tri.vertex[0][0] << " " << tri.vertex[0][1] << " "
      << tri.vertex[0][2] << "\n";
    f << "        vertex " << tri.vertex[1][0] << " " << tri.vertex[1][1] << " "
      << tri.vertex[1][2] << "\n";
    f << "        vertex " << tri.vertex[2][0] << " " << tri.vertex[2][1] << " "
      << tri.vertex[2][2] << "\n";
    f << "    endloop\n";
    f << "endfacet\n";
  }
  f << "endsolid\n";
  f.close();
  //
  std::cout << "[StlMesh]: Write to " << fullName << "complete." << std::endl;
}

template <typename T>
bool StlMesh<T>::testRayIntersect(const std::set<unsigned int> &tris,
                                  const Vector<T, 3> &pt,
                                  const Vector<T, 3> &dir, Vector<T, 3> &q,
                                  T &alpha) {
  std::set<unsigned int>::iterator it = tris.begin();
  for (; it != tris.end(); ++it) {
    if (Triangles[*it].testRayIntersect(pt, dir, q, alpha) && alpha < 1) {
      return true;
    }
  }
  return false;
}

template <typename T>
StlReader<T>::StlReader(std::string fName, T voxelSize, T stlSize, int method,
                        bool verbose, T overlap, T max)
    : _voxelSize(voxelSize),
      _stlSize(stlSize),
      _overlap(overlap),
      _fName(fName),
      _mesh(fName, stlSize),
      _verbose(verbose),
      _myMin(),
      _myMax() {
  Vector<T, 3> extension = _mesh.getMax() - _mesh.getMin();
  if (util::nearZero(max)) {
    max = std::max(extension[0], std::max(extension[1], extension[2])) +
          _voxelSize;
  }
  // max depth
  int depth = 0;
  // 2 ^ depth > max, radius = 2^(depth-1) * voxelSize(1)
  for (; _voxelSize * std::pow(2, depth) < max; depth++)
    ;
  Vector<T, 3> center;
  T radius = _voxelSize * std::pow(2, depth - 1);

  /// Find center of tree and move by _voxelSize/4.
  for (unsigned i = 0; i < 3; i++) {
    // center[i] = (_mesh.getMin()[i] + _mesh.getMax()[i]) / 2. - _voxelSize
    // / 4.;
    center[i] = (_mesh.getMin()[i] + _mesh.getMax()[i]) / 2.;
  }

  /// Create tree
  // the first octree: radius = 2^(depth-1) * voxelSize(1) with maxdepth = depth
  // minimal maxdepth = 0, min radius = 2^(-1) *voxelSize
  _tree = new Octree<T>(center, radius, &_mesh, depth, _overlap);

  /// Compute _myMin, _myMax such that they are the smallest (greatest) Voxel
  /// inside the STL.
  for (int i = 0; i < 3; i++) {
    _myMin[i] = center[i] + _voxelSize / 2.;
    _myMax[i] = center[i] - _voxelSize / 2.;
  }
  for (int i = 0; i < 3; i++) {
    while (_myMin[i] > _mesh.getMin()[i]) {
      _myMin[i] -= _voxelSize;
    }
    while (_myMax[i] < _mesh.getMax()[i]) {
      _myMax[i] += _voxelSize;
    }
    _myMax[i] -= _voxelSize;
    _myMin[i] += _voxelSize;
  }

  /// Indicate nodes of the tree. (Inside/Outside)
  switch (method) {
    case 1:
      indicate1();
      break;
    case 3:
      indicate3();
      break;
    case 4:
      indicate2_Xray();
      break;
    case 5:
      indicate2_Yray();
      break;
    default:
      indicate2();
      break;
  }

  // if (_verbose) {
  // print();
  // }
  _mesh.print();
  std::vector<Octree<T>*> leafs;
  _tree->getLeafs(leafs);
  _tree->print(leafs);
}

template <typename T>
void StlReader<T>::indicate1() {
  std::vector<Octree<T> *> leafs;
  _tree->getLeafs(leafs);
  typename std::vector<Octree<T> *>::iterator it = leafs.begin();
  Vector<T, 3> dir, pt, s;

  int intersections = 0;
  int inside = 0;
  Octree<T> *node = nullptr;
  T step = 1. / 1000. * _voxelSize;
  for (; it != leafs.end(); ++it) {
    inside = 0;

    pt = (*it)->getCenter();
    intersections = 0;
    s = pt;  // + step;

    /// X+ dir
    dir[0] = 1;
    dir[1] = 0;
    dir[2] = 0;
    while (s[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections += node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
    }
    inside += (intersections % 2);

    /// Y+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 1;
    dir[2] = 0;
    while (s[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections += node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
    }
    inside += (intersections % 2);

    /// Z+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 0;
    dir[2] = 1;
    while (s[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections += node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
    }
    inside += (intersections % 2);
    (*it)->setInside(inside > 1);
  }
}

template <typename T>
void StlReader<T>::indicate2() {
  T rad = _tree->getRadius();
  Vector<T, 3> rayPt = _tree->getCenter() - rad + .5 * _voxelSize;
  // starting point of the ray
  Vector<T, 3> pt = rayPt;
  // direction of the ray: along the positive z-axis
  Vector<T, 3> rayDir;
  rayDir[0] = 0.;
  rayDir[1] = 0.;
  rayDir[2] = 1.;
  // Vector<T,3> maxEdge = _tree->getCenter() + rad;

  T step = 1. / 1000. * _voxelSize;

  Octree<T> *node = nullptr;
  unsigned short rayInside = 0;
  Vector<T, 3> nodeInters;
  while (pt[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
    node = _tree->find(pt);
    nodeInters = pt;
    nodeInters[2] = node->getCenter()[2] - node->getRadius();
    rayInside = 0;
    while (pt[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(pt);
      nodeInters = pt;
      nodeInters[2] = node->getCenter()[2] - node->getRadius();
      rayInside = 0;
      while (pt[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
        node = _tree->find(pt);
        node->checkRay(nodeInters, rayDir, rayInside);
        node->intersectRayNode(pt, rayDir, nodeInters);
        pt = nodeInters + step * rayDir;
      }
      pt[2] = rayPt[2];
      pt[1] += _voxelSize;
    }
    pt[1] = rayPt[1];
    pt[0] += _voxelSize;
  }
}

template <typename T>
void StlReader<T>::indicate2_Xray() {
  T rad = _tree->getRadius();
  Vector<T, 3> rayPt = _tree->getCenter() - rad + .5 * _voxelSize;
  Vector<T, 3> pt = rayPt;
  Vector<T, 3> rayDir;
  rayDir[0] = 1.;
  rayDir[1] = 0.;
  rayDir[2] = 0.;
  // Vector<T,3> maxEdge = _tree->getCenter() + rad;

  T step = 1. / 1000. * _voxelSize;

  Octree<T> *node = nullptr;
  unsigned short rayInside = 0;
  Vector<T, 3> nodeInters;
  while (pt[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
    node = _tree->find(pt);
    nodeInters = pt;
    nodeInters[0] = node->getCenter()[0] - node->getRadius();
    rayInside = 0;
    while (pt[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(pt);
      nodeInters = pt;
      nodeInters[0] = node->getCenter()[0] - node->getRadius();
      rayInside = 0;
      while (pt[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
        node = _tree->find(pt);
        node->checkRay(nodeInters, rayDir, rayInside);
        node->intersectRayNode(pt, rayDir, nodeInters);
        pt = nodeInters + step * rayDir;
      }
      pt[0] = rayPt[0];
      pt[1] += _voxelSize;
    }
    pt[1] = rayPt[1];
    pt[2] += _voxelSize;
  }
}

template <typename T>
void StlReader<T>::indicate2_Yray() {
  T rad = _tree->getRadius();
  Vector<T, 3> rayPt = _tree->getCenter() - rad + .5 * _voxelSize;
  Vector<T, 3> pt = rayPt;
  Vector<T, 3> rayDir;
  rayDir[0] = 0.;
  rayDir[1] = 1.;
  rayDir[2] = 0.;
  // Vector<T,3> maxEdge = _tree->getCenter() + rad;

  T step = 1. / 1000. * _voxelSize;

  Octree<T> *node = nullptr;
  unsigned short rayInside = 0;
  Vector<T, 3> nodeInters;
  while (pt[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
    node = _tree->find(pt);
    nodeInters = pt;
    nodeInters[1] = node->getCenter()[1] - node->getRadius();
    rayInside = 0;
    while (pt[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(pt);
      nodeInters = pt;
      nodeInters[1] = node->getCenter()[1] - node->getRadius();
      rayInside = 0;
      while (pt[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
        node = _tree->find(pt);
        node->checkRay(nodeInters, rayDir, rayInside);
        node->intersectRayNode(pt, rayDir, nodeInters);
        pt = nodeInters + step * rayDir;
      }
      pt[1] = rayPt[1];
      pt[0] += _voxelSize;
    }
    pt[0] = rayPt[0];
    pt[2] += _voxelSize;
  }
}

template <typename T>
void StlReader<T>::indicate3() {
  std::vector<Octree<T> *> leafs;
  _tree->getLeafs(leafs);
  typename std::vector<Octree<T> *>::iterator it = leafs.begin();

  Vector<T, 3> dir, pt, s;
  Octree<T> *node = nullptr;
  T step = 1. / 1000. * _voxelSize;
  int intersections;
  int sum_intersections;

  for (; it != leafs.end(); ++it) {
    pt = (*it)->getCenter();
    intersections = 0;
    sum_intersections = 0;
    s = pt;  // + step;

    /// X+ dir
    dir[0] = 1;
    dir[1] = 0;
    dir[2] = 0;
    while (s[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Y+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 1;
    dir[2] = 0;
    while (s[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Z+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 0;
    dir[2] = 1;
    while (s[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// X- dir
    intersections = 0;
    s = pt;  // + step;
    dir[0] = -1;
    dir[1] = 0;
    dir[2] = 0;
    while (s[0] > _mesh.getMin()[0] - std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Y- Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = -1;
    dir[2] = 0;
    while (s[1] > _mesh.getMin()[1] - std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Z- Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 0;
    dir[2] = -1;
    while (s[2] > _mesh.getMin()[2] - std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }
    (*it)->setInside(sum_intersections > 5);
  }
}

template <typename T>
bool StlReader<T>::operator()(bool output[], const T input[]) {
  output[0] = false;
  T coords = _tree->getRadius();
  Vector<T, 3> c(_tree->getCenter());
  if (c[0] - coords < input[0] && input[0] < c[0] + coords &&
      c[1] - coords < input[1] && input[1] < c[1] + coords &&
      c[2] - coords < input[2] && input[2] < c[2] + coords) {
    std::vector<T> tmp(input, input + 3);
    output[0] = _tree->find(tmp)->getInside();
  }
  return true;
}

template <typename T>
bool StlReader<T>::distance(T &distance, const Vector<T, 3> &origin,
                            const Vector<T, 3> &direction, int iC) {
  Octree<T> *node = nullptr;
  Vector<T, 3> dir(direction);
  dir.normalize();
  Vector<T, 3> extends = _mesh.getMax() - _mesh.getMin();
  Vector<T, 3> pt(origin);
  Vector<T, 3> q;
  Vector<T, 3> s;
  Vector<T, 3> center = _mesh.getMin() + 1 / 2. * extends;
  T step = _voxelSize / 1000., a = 0;

  for (int i = 0; i < 3; i++) {
    extends[i] /= 2.;
  }

  if (!(_mesh.getMin()[0] < origin[0] && origin[0] < _mesh.getMax()[0] &&
        _mesh.getMin()[1] < origin[1] && origin[1] < _mesh.getMax()[1] &&
        _mesh.getMin()[2] < origin[2] && origin[2] < _mesh.getMax()[2])) {
    T t = T(), d = T();
    bool foundQ = false;

    if (dir[0] > 0) {
      d = _mesh.getMin()[0];
      t = (d - origin[0]) / dir[0];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];

      if (_mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1] &&
          _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    } else if (dir[0] < 0) {
      d = _mesh.getMax()[0];
      t = (d - origin[0]) / dir[0];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1] &&
          _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    }

    if (dir[1] > 0 && !foundQ) {
      d = _mesh.getMin()[1];
      t = (d - origin[1]) / dir[1];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0] &&
          _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    } else if (dir[1] < 0 && !foundQ) {
      d = _mesh.getMax()[1];
      t = (d - origin[1]) / dir[1];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0] &&
          _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    }

    if (dir[2] > 0 && !foundQ) {
      d = _mesh.getMin()[2];
      t = (d - origin[2]) / dir[2];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0] &&
          _mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]) {
        foundQ = true;
      }
    } else if (dir[2] < 0 && !foundQ) {
      d = _mesh.getMax()[2];
      t = (d - origin[2]) / dir[2];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0] &&
          _mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]) {
        foundQ = true;
      }
    }

    if (!foundQ) {
      return false;
    }
  }

  while ((std::fabs(pt[0] - center[0]) < extends[0]) &&
         (std::fabs(pt[1] - center[1]) < extends[1]) &&
         (std::fabs(pt[2] - center[2]) < extends[2])) {
    node = _tree->find(pt);
    if (node->closestIntersection(Vector<T, 3>(origin), dir, q, a)) {
      Vector<T, 3> vek(q - Vector<T, 3>(origin));
      distance = vek.getnorm();
      return true;
    } else {
      Octree<T> *tmpNode = _tree->find(pt);
      tmpNode->intersectRayNode(pt, dir, s);
      for (int i = 0; i < 3; i++) {
        pt[i] = s[i] + step * dir[i];
      }
    }
  }

  return false;
}

template <typename T>
Vector<T, 3> StlReader<T>::findNormalOnSurface(const Vector<T, 3> &pt) {
  // Check if the position is on the corner of a triangle
  unsigned countTriangles = 0;
  Vector<T, 3> normal(T(0));

  for (const Triangle<T> &triangle : _mesh.getTriangles()) {
    if (triangle.IsInside(pt)) {
      ++countTriangles;
      normal = normal + triangle.getNormal();
    }
  }
  if (countTriangles > 0) {
    return normal / countTriangles;
  }

  // if nothing was found return (0,0,0) to indicate that nothing was found
  return normal;
}

template <typename T>
Vector<T, 3> StlReader<T>::evalSurfaceNormal(const Vector<T, 3> &origin) {
  Vector<T, 3> normal(0.);
  Vector<T, 3> closestPoint;
  T distance = std::numeric_limits<T>::max();
  for (const Triangle<T> &triangle : _mesh.getTriangles()) {
    Vector<T, 3> const pointOnTriangle =
        triangle.closestPtPointTriangle(origin);
    Vector<T, 3> const currDistance = pointOnTriangle - origin;
    T currDistanceNorm = currDistance.getnorm();
    if (util::nearZero(currDistanceNorm)) {
      return findNormalOnSurface(origin);
    } else if (distance > currDistanceNorm) {
      normal = currDistance;
      distance = currDistanceNorm;
      closestPoint = pointOnTriangle;
    }
  }

  if (!util::nearZero(normal.getnorm())) {
    normal = normal.normalize();
  } else {
    return normal;
  }

  // Possible change of the sign so that the normal fits to the SDF logic
  if (distance < _voxelSize) {
    bool isInsideInnerPoint;
    this->operator()(&isInsideInnerPoint,
                     (closestPoint - _voxelSize * normal).data());
    bool isInsideOuterPoint;
    this->operator()(&isInsideOuterPoint,
                     (closestPoint + _voxelSize * normal).data());
    normal = normal*(isInsideInnerPoint && !isInsideOuterPoint ? 1 : -1);
  } else {
    bool isInside;
    this->operator()(&isInside, origin.data());
    normal = normal*(isInside ? 1 : -1);
  }
  return normal;
}

template <typename T>
inline T StlReader<T>::signedDistance(const Vector<T, 3> &input) {
  bool isInside;
  this->operator()(&isInside, input.data());
  const short sign = (isInside ? -1 : 1);

  T distance = std::numeric_limits<T>::max();
  for (const Triangle<T> &triangle : _mesh.getTriangles()) {
    Vector<T, 3> const pointOnTriangle = triangle.closestPtPointTriangle(input);
    T currDistance = (pointOnTriangle - input).norm();
    distance = std::min(distance, currDistance);
  }

  return distance * sign;
}

template <typename T>
void StlReader<T>::print() {
  _mesh.print();
  _tree->print();
}

template <typename T>
void StlReader<T>::setNormalsOutside() {
  unsigned int noTris = _mesh.triangleSize();
  Vector<T, 3> center;
  // Octree<T>* node = nullptr;
  for (unsigned int i = 0; i < noTris; i++) {
    center = _mesh.getTri(i).getCenter();
    if (_tree
            ->find(center + _mesh.getTri(i).normal * std::sqrt(3.) * _voxelSize)
            ->getInside()) {
      //      cout << "Wrong direction" << std::endl;
      Vector<T, 3> pt(_mesh.getTri(i).point[0].coords);
      _mesh.getTri(i).point[0].coords = _mesh.getTri(i).point[2].coords;
      _mesh.getTri(i).point[2].coords = pt;
      _mesh.getTri(i).Init();
      //      _mesh.getTri(i).getNormal()[0] *= -1.;
      //      _mesh.getTri(i).getNormal()[1] *= -1.;
      //      _mesh.getTri(i).getNormal()[2] *= -1.;
    }
  }
}

template <typename T>
void StlReader<T>::setBoundaryInsideNodes() {
  std::vector<Octree<T> *> leafs;
  _tree->getLeafs(leafs);
  for (auto it = leafs.begin(); it != leafs.end(); ++it) {
    if ((*it)->getBoundaryNode()) {
      (*it)->setInside(true);
    }
  }
}
