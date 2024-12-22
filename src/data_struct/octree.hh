/* This file is part of FreeLB, modified from OpenLB's octree.hh, with the following
 * copyright notice:
 *
 * // start of the original OpenLB's copyright notice
 *
 * This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Thomas Henn
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

// octree.hh
#pragma once

#include <iomanip>
#include <sstream>

#include "data_struct/octree.h"
#include "utils/directories.h"

template <typename T>
Octree<T>::Octree(Vector<T, 3> center, T rad, StlMesh<T>* mesh, short maxDepth, T overlap,
                  Octree<T>* parent)
    : _center(center), _radius(rad), _mesh(mesh), _maxDepth(maxDepth), _isLeaf(false),
      _boundaryNode(false), _inside(false), _parent(parent), _child(nullptr) {
  findTriangles(overlap);
  //  cout << _triangles.size() << std::endl;
  if (_triangles.size() > 0 && 0 < _maxDepth) {
    _child = new Octree<T>*[8];

    Vector<T, 3> tmpCenter = _center;
    T tmpRad = _radius / 2.;
    tmpCenter[0] = _center[0] - tmpRad;
    tmpCenter[1] = _center[1] - tmpRad;
    tmpCenter[2] = _center[2] + tmpRad;
    _child[0] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

    tmpCenter[0] = _center[0] + tmpRad;
    tmpCenter[1] = _center[1] - tmpRad;
    tmpCenter[2] = _center[2] + tmpRad;
    _child[1] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

    tmpCenter[0] = _center[0] - tmpRad;
    tmpCenter[1] = _center[1] - tmpRad;
    tmpCenter[2] = _center[2] - tmpRad;
    _child[2] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

    tmpCenter[0] = _center[0] + tmpRad;
    tmpCenter[1] = _center[1] - tmpRad;
    tmpCenter[2] = _center[2] - tmpRad;
    _child[3] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

    tmpCenter[0] = _center[0] - tmpRad;
    tmpCenter[1] = _center[1] + tmpRad;
    tmpCenter[2] = _center[2] + tmpRad;
    _child[4] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

    tmpCenter[0] = _center[0] + tmpRad;
    tmpCenter[1] = _center[1] + tmpRad;
    tmpCenter[2] = _center[2] + tmpRad;
    _child[5] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

    tmpCenter[0] = _center[0] - tmpRad;
    tmpCenter[1] = _center[1] + tmpRad;
    tmpCenter[2] = _center[2] - tmpRad;
    _child[6] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

    tmpCenter[0] = _center[0] + tmpRad;
    tmpCenter[1] = _center[1] + tmpRad;
    tmpCenter[2] = _center[2] - tmpRad;
    _child[7] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth - 1, overlap, this);

  } else {
    _isLeaf = true;
    if (_triangles.size() > 0) {
      _boundaryNode = true;
    }
  }
}

template <typename T>
Octree<T>::~Octree() {
  if (_maxDepth != 0 && !_isLeaf) {
    for (int i = 0; i < 8; i++) {
      delete _child[i];
    }
    delete[] _child;
  }
}

template <typename T>
void Octree<T>::findTriangles(T overlap) {
  // inline Triangle<T> &getTri(unsigned int i) { return Triangles[i]; }
  // inline unsigned int triangleSize() const { return Triangles.size(); }
  if (_parent == nullptr) {
    _triangles.reserve(_mesh->triangleSize());
    for (unsigned int i = 0; i < _mesh->triangleSize(); ++i) {
      if (AABBTri(_mesh->getTri(i))) {
        _triangles.push_back(i);
      }
    }
  } else {
    std::vector<unsigned int>::iterator it;
    for (it = _parent->_triangles.begin(); it != _parent->_triangles.end(); ++it) {
      if (AABBTri(_mesh->getTri(*it), overlap)) {
        _triangles.push_back(*it);
      }
    }
  }
}

// check if a triangle intersects with an axis-aligned bounding box (AABB)

// The intersection test is then performed using the separating axis theorem
// (SAT), which states that if there exists an axis along which the projections
// of the AABB and the triangle do not overlap, then the AABB and the triangle
// do not intersect. The code checks for this condition along various axes.

// The p0, p1, and r variables are used to store the projections of the triangle
// vertices and the AABB extents along the current axis. The if statement checks
// if the maximum projection of the triangle is less than the minimum projection
// of the AABB (or vice versa), which would mean that the AABB and the triangle
// do not intersect along the current axis.

template <typename T>
bool Octree<T>::AABBTri(const Triangle<T>& tri, T overlap) {
  std::vector<T> v0(3, T()), v1(3, T()), v2(3, T()), f0(3, T()), f1(3, T()), f2(3, T()),
    e(3, T());

  /* Test intersection cuboids - triangle
   * Intersection test after Christer Ericson - Real time Collision Detection p.
   * TestTriangleAABB p.171 */
  Vector<T, 3> c(_center);
  constexpr T eps = std::numeric_limits<T>::epsilon();

  for (int j = 0; j < 3; j++) {
    v0[j] = tri.vertex[0][j] - _center[j];
    v1[j] = tri.vertex[1][j] - _center[j];
    v2[j] = tri.vertex[2][j] - _center[j];
    e[j] = _radius * 1.01 + overlap;  // + std::numeric_limits<T>::epsilon(); // *1.01;
  }
  for (int j = 0; j < 3; j++) {
    f0[j] = v1[j] - v0[j];
    f1[j] = v2[j] - v1[j];
    f2[j] = v0[j] - v2[j];
  }
  T p0 = T(), p1 = T(), r = T();

  // result[0] = a[1] * b[2] - a[2] * b[1];
  // result[1] = a[2] * b[0] - a[0] * b[2];
  // result[2] = a[0] * b[1] - a[1] * b[0];

  // a00 = u0 x f0 = (1, 0, 0) x f0 = (0, −f0z, f0y)
  // a01 = u0 x f1 = (1, 0, 0) x f1 = (0, −f1z, f1y)
  // a02 = u0 x f2 = (1, 0, 0) x f2 = (0, −f2z, f2y)
  // a10 = u1 x f0 = (0, 1, 0) x f0 = ( f0z, 0, −f0x )
  // a11 = u1 x f1 = (0, 1, 0) x f1 = ( f1z, 0, −f1x )
  // a12 = u1 x f2 = (0, 1, 0) x f2 = ( f2z, 0, −f2x )
  // a20 = u2 x f0 = (0, 0, 1) x f0 = (−f0y, f0x , 0)
  // a21 = u2 x f1 = (0, 0, 1) x f1 = (−f1y, f1x , 0)
  // a22 = u2 x f2 = (0, 0, 1) x f2 = (−f2y, f2x , 0)

  // the projection radius of a box with respect to an axis n is given by
  // r = e0 |u0 · n| + e1 |u1 · n| + e2 |u2 · n| .

  // In the case of n = a00, this simplifies to
  // r = e0 |u0 · a00| + e1 |u1 · a00| + e2 |u2 · a00| ⇔
  // r = e0 |0| + e1 |−a00z| + e2

  // test a00
  // checking for separation along an axis perpendicular to one of the edges of
  // the triangle x - component of cross product/ projection on x-axis
  // p0 = v1 x v0
  // p1 = f0 x v2 = (v1 - v0) x v2
  // projection of the AABB onto the axis
  // finding the maximum extent of the triangle's projection along the axis
  p0 = v0[2] * v1[1] - v0[1] * v1[2];
  p1 = v2[2] * v1[1] - v2[2] * v0[1] + v0[2] * v2[1] - v1[2] * v2[1];
  r = e[1] * std::fabs(f0[2]) + e[2] * std::fabs(f0[1]);
  T mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r + eps) {
    return false;
  }

  // test a01
  // x - component of cross product/ projection on x-axis
  // p0 = f1 x v0 = (v2 - v1) x v0
  // p1 = v2 x v1
  p0 = v0[1] * v1[2] - v0[1] * v2[2] - v0[2] * v1[1] + v0[2] * v2[1];
  p1 = -v1[1] * v2[2] + v1[2] * v2[1];
  r = e[1] * std::fabs(f1[2]) + e[2] * std::fabs(f1[1]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r + eps) {
    return false;
  }

  // test a02
  // x - component of cross product/ projection on x-axis
  // p0 = v0 x v2
  // p1 = f2 x v1 = (v0 - v2) x v1
  p0 = v0[1] * v2[2] - v0[2] * v2[1];
  p1 = v0[1] * v1[2] - v0[2] * v1[1] + v1[1] * v2[2] - v1[2] * v2[1];
  r = e[1] * std::fabs(f2[2]) + e[2] * std::fabs(f2[1]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r + eps) {
    return false;
  }

  // test a10
  p0 = v0[0] * v1[2] - v0[2] * v1[0];
  p1 = v0[0] * v2[2] - v0[2] * v2[0] - v1[0] * v2[2] + v1[2] * v2[0];
  r = e[0] * std::fabs(f0[2]) + e[2] * std::fabs(f0[0]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r + eps) {
    return false;
  }

  // test a11
  p0 = -v0[0] * v1[2] + v0[0] * v2[2] + v0[2] * v1[0] - v0[2] * v2[0];
  p1 = v1[0] * v2[2] - v1[2] * v2[0];
  r = (T)(e[0] * std::fabs(f1[2]) + e[2] * std::fabs(f1[0]));
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r + eps) {
    return false;
  }

  // test a12
  p0 = -v0[0] * v2[2] + v0[2] * v2[0];
  p1 = -v0[0] * v1[2] + v0[2] * v1[0] - v1[0] * v2[2] + v1[2] * v2[0];
  r = e[0] * std::fabs(f2[2]) + e[2] * std::fabs(f2[0]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r + eps) {
    return false;
  }

  // test a20
  p0 = -v0[0] * v1[1] + v0[1] * v1[0];
  p1 = -v0[0] * v2[1] + v0[1] * v2[0] + v1[0] * v2[1] - v1[1] * v2[0];
  r = e[0] * std::fabs(f0[1]) + e[1] * std::fabs(f0[0]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r + eps) {
    return false;
  }

  // test a21
  p0 = v0[0] * v1[1] - v0[0] * v2[1] - v0[1] * v1[0] + v0[1] * v2[0];
  p1 = -v1[0] * v2[1] + v1[1] * v2[0];
  r = e[0] * std::fabs(f1[1]) + e[1] * std::fabs(f1[0]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r + eps) {
    return false;
  }

  // test a22
  p0 = v0[0] * v2[1] - v0[1] * v2[0];
  p1 = v0[0] * v1[1] - v0[1] * v1[0] + v1[0] * v2[1] - v1[1] * v2[0];
  r = e[0] * std::fabs(f2[1]) + e[1] * std::fabs(f2[0]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r + eps) {
    return false;
  }

  // Test the three axes corresponding to the face normals of AABB b (category
  // 1)

  if (std::max(std::max(v0[0], v1[0]), v2[0]) < -e[0] ||
      std::min(std::min(v0[0], v1[0]), v2[0]) > e[0]) {
    return false;
  }
  if (std::max(std::max(v0[1], v1[1]), v2[1]) < -e[1] ||
      std::min(std::min(v0[1], v1[1]), v2[1]) > e[1]) {
    return false;
  }
  if (std::max(std::max(v0[2], v1[2]), v2[2]) < -e[2] ||
      std::min(std::min(v0[2], v1[2]), v2[2]) > e[2]) {
    return false;
  }

  //  Test separating axis corresponding to triangle face normal (category 2)

  /* Test intersection cuboids - triangle plane*/
  r = e[0] * std::fabs(tri.normal[0]) + e[1] * std::fabs(tri.normal[1]) +
      e[2] * std::fabs(tri.normal[2]);
  T s = tri.normal[0] * c[0] + tri.normal[1] * c[1] + tri.normal[2] * c[2] - tri.d;
  return (std::fabs(s) <= r);
}

template <typename T>
Octree<T>* Octree<T>::find(const Vector<T, 3>& pt, const int& maxDepth) {
  //  clout << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
  // modified < to <= to include points on the boundary
  if (_isLeaf || maxDepth == _maxDepth) {
    if (std::abs(_center[0] - pt[0]) <= _radius + std::numeric_limits<T>::epsilon() &&
        std::abs(_center[1] - pt[1]) <= _radius + std::numeric_limits<T>::epsilon() &&
        std::abs(_center[2] - pt[2]) <= _radius + std::numeric_limits<T>::epsilon()) {
      //       clout << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
      return this;
    } else {
      // OstreamManager clout(std::cout, "Octree");

      // std::cout << "Point: " << std::setprecision(10) << pt[0] << " " <<
      // pt[1]
      //           << " " << pt[2] << " " << std::endl;
      // std::cout << "Center: " << std::setprecision(10) << _center[0] << " "
      //           << _center[1] << " " << _center[2] << " " << std::endl;
      // std::cout << "Radius: " << std::setprecision(10) << _radius <<
      // std::endl;

      // throw std::runtime_error("[Octree->find] Point outside of geometry.");
      return nullptr;
    }
  } else {
    if (pt[0] < _center[0]) {
      if (pt[1] < _center[1]) {
        if (pt[2] < _center[2]) {
          return _child[2]->find(pt, maxDepth);
        } else {
          return _child[0]->find(pt, maxDepth);
        }
      } else {
        if (pt[2] < _center[2]) {
          return _child[6]->find(pt, maxDepth);
        } else {
          return _child[4]->find(pt, maxDepth);
        }
      }
    } else {
      if (pt[1] < _center[1]) {
        if (pt[2] < _center[2]) {
          return _child[3]->find(pt, maxDepth);
        } else {
          return _child[1]->find(pt, maxDepth);
        }
      } else {
        if (pt[2] < _center[2]) {
          return _child[7]->find(pt, maxDepth);
        } else {
          return _child[5]->find(pt, maxDepth);
        }
      }
    }
  }
}

template <typename T>
int Octree<T>::testIntersection(const Vector<T, 3>& pt, const Vector<T, 3>& dir,
                                bool print) {
  int intersections = 0;
  Vector<T, 3> q;
  std::vector<Vector<T, 3> > qs;
  T a;

  for (unsigned k = 0; k < _triangles.size(); ++k) {
    if (_mesh->getTri(_triangles[k]).testRayIntersect(pt, dir, q, a)) {
      if (std::fabs(_center[0] - q[0]) <=
            _radius + std::numeric_limits<T>::epsilon() + 1 / 1000. * _radius &&
          std::fabs(_center[1] - q[1]) <=
            _radius + std::numeric_limits<T>::epsilon() + 1 / 1000. * _radius &&
          std::fabs(_center[2] - q[2]) <=
            _radius + std::numeric_limits<T>::epsilon() + 1 / 1000. * _radius) {
        bool newpoint = true;
        for (unsigned i = 0; i < qs.size(); i++) {
          newpoint =
            (!util::nearZero(q[0] - qs[i][0]) || !util::nearZero(q[1] - qs[i][1]) ||
             !util::nearZero(q[2] - qs[i][2]));
        }
        if (newpoint) {
          qs.push_back(q);
          intersections++;
        }
      }
    }
  }
  return intersections;
}

template <typename T>
void Octree<T>::checkRay(const Vector<T, 3>& pt, const Vector<T, 3>& dir,
                         unsigned short& rayInside) {
  unsigned short left = 0, right = 0;
  Vector<T, 3> dirNormed(dir);
  dirNormed.normalize();
  dirNormed = dirNormed * (_radius * 2.);
  Vector<T, 3> q;
  std::vector<Vector<T, 3> > qs;
  T a = 1.;

  for (unsigned int k = 0; k < _triangles.size(); ++k) {
    if (_mesh->getTri(_triangles[k]).testRayIntersect(pt, dirNormed, q, a, 0.) &&
        a < 1.) {
      bool newpoint = true;
      for (unsigned int i = 0; i < qs.size(); i++) {
        newpoint &=
          (!util::nearZero(q[0] - qs[i][0]) || !util::nearZero(q[1] - qs[i][1]) ||
           !util::nearZero(q[2] - qs[i][2]));
      }
      if (newpoint) {
        qs.push_back(q);
        if (a < .5) {
          left++;
        } else {
          right++;
        }
      }
    }
  }
  rayInside += left;
  rayInside %= 2;
  setInside(rayInside);
  rayInside += right;
  rayInside %= 2;
}

template <typename T>
void Octree<T>::getCenterpoints(std::vector<std::vector<T> >& pts) {
  if (_isLeaf) {
    pts.push_back(_center);
  } else {
    for (int i = 0; i < 8; i++) {
      _child[i]->getCenterpoints(pts);
    }
  }
}

template <typename T>
void Octree<T>::getLeafs(std::vector<Octree<T>*>& pts) {
  if (_isLeaf) {
    pts.push_back(this);
  } else {
    for (int i = 0; i < 8; i++) {
      _child[i]->getLeafs(pts);
    }
  }
}

template <typename T>
T Octree<T>::getMinRadius() {
  std::vector<Octree<T>*> leafs;
  getLeafs(leafs);
  T minRadius = _radius;
  typename std::vector<Octree<T>*>::iterator it = leafs.begin();
  for (; it != leafs.end(); ++it) {
    if ((*it)->getMaxdepth() == 0) {
      minRadius = std::min(minRadius, (*it)->getRadius());
      break;
    }
  }
  return minRadius;
}

template <typename T>
T Octree<T>::getMinRadius(std::vector<Octree<T>*>& leafs) {
  T minRadius = _radius;
  typename std::vector<Octree<T>*>::iterator it = leafs.begin();
  for (; it != leafs.end(); ++it) {
    if ((*it)->getMaxdepth() == 0) {
      minRadius = std::min(minRadius, (*it)->getRadius());
      break;
    }
  }
  return minRadius;
}

template <typename T>
void Octree<T>::write(const Vector<T, 3>& pt, const std::string no) {
  if (_triangles.size() > 0 && (std::fabs(pt[0] - _center[0]) < _radius &&
                                std::fabs(pt[1] - _center[1]) < _radius &&
                                std::fabs(pt[2] - _center[2]) < _radius)) {
    DirCreator::Create_Dir("./Octree");
    std::string fullName = "./Octree/Octree_" + no + ".stl";
    std::ofstream f(fullName.c_str());
    if (!f) {
      std::cerr << "[Octree] could not open file: " << fullName << std::endl;
    }
    f << "solid ascii" << std::endl;
    std::vector<unsigned int>::iterator it = _triangles.begin();
    for (; it != _triangles.end(); ++it) {
      f << "facet normal" << _mesh->getTri(*it).normal[0] << " "
        << _mesh->getTri(*it).normal[1] << " " << _mesh->getTri(*it).normal[2] << " "
        << std::endl;
      f << "    outer loop\n";
      f << "        vertex " << _mesh->getTri(*it).vertex[0][0] << " "
        << _mesh->getTri(*it).vertex[0][1] << " " << _mesh->getTri(*it).vertex[0][2]
        << "\n";
      f << "        vertex " << _mesh->getTri(*it).vertex[1][0] << " "
        << _mesh->getTri(*it).vertex[1][1] << " " << _mesh->getTri(*it).vertex[1][2]
        << "\n";
      f << "        vertex " << _mesh->getTri(*it).vertex[2][0] << " "
        << _mesh->getTri(*it).vertex[2][1] << " " << _mesh->getTri(*it).vertex[2][2]
        << "\n";
      f << "    endloop\n";
      f << "endfacet\n";
    }
    f.close();
  }
  if (!_isLeaf) {
    for (int i = 0; i < 8; i++) {
      std::stringstream istr;
      istr << i;
      _child[i]->write(pt, no + istr.str());
    }
  }
}

template <typename T>
void Octree<T>::write(const int depth, const std::string no) {
  if (_triangles.size() > 0 && _maxDepth == depth) {
    DirCreator::Create_Dir("./Octree");
    std::string fullName = "./Octree/Octree_" + no + ".stl";
    std::ofstream f(fullName.c_str());
    if (!f) {
      std::cerr << "[Octree] could not open file: " << fullName << std::endl;
    }
    f << "solid ascii" << std::endl;
    std::vector<unsigned int>::iterator it = _triangles.begin();
    for (; it != _triangles.end(); ++it) {
      f << "facet normal" << _mesh->getTri(*it).normal[0] << " "
        << _mesh->getTri(*it).normal[1] << " " << _mesh->getTri(*it).normal[2] << " "
        << std::endl;
      f << "    outer loop\n";
      f << "        vertex " << _mesh->getTri(*it).vertex[0][0] << " "
        << _mesh->getTri(*it).vertex[0][1] << " " << _mesh->getTri(*it).vertex[0][2]
        << "\n";
      f << "        vertex " << _mesh->getTri(*it).vertex[1][0] << " "
        << _mesh->getTri(*it).vertex[1][1] << " " << _mesh->getTri(*it).vertex[1][2]
        << "\n";
      f << "        vertex " << _mesh->getTri(*it).vertex[2][0] << " "
        << _mesh->getTri(*it).vertex[2][1] << " " << _mesh->getTri(*it).vertex[2][2]
        << "\n";
      f << "    endloop\n";
      f << "endfacet\n";
    }
    f.close();
  }
  if (!_isLeaf) {
    for (int i = 0; i < 8; i++) {
      std::stringstream istr;
      istr << i;
      _child[i]->write(depth, no + istr.str());
    }
  }
}
//
// vtk file format:
// # vtk DataFile Version 2.0
// Octree
// ASCII
// DATASET UNSTRUCTURED_GRID
// POINTS n type    leafs.size()*8  float
// Cells n size     leafs.size()   leafs.size()*9
// CELL_TYPES n
// CELL_DATA n
// Dataset Attribute DataName type
// e.g.:
// SCALARS dataName dataType
// VECTORS dataName dataType
// LOOKUP_TABLE default
// the number of cells n and the size of the cell list size. The cell list size
// is the total number of integer values required to represent the list (i.e.,
// sum of numPoints and connectivity indices over each cell: 8+1)
template <typename T>
void Octree<T>::write(const std::string fName) {
  MPI_RANK(0)
  std::vector<Octree<T>*> leafs;
  getLeafs(leafs);
  if (leafs.size() > 100000) {
    std::cout << "Octree too large:" << leafs.size() << " ,Press any key to continue..."
              << std::endl;
    std::cin.get();
  }
  typename std::vector<Octree<T>*>::iterator it = leafs.begin();
  DirCreator::Create_Dir("./vtkoutput");
  std::string fullName = "./vtkoutput/" + fName + ".vtk";
  std::ofstream f(fullName.c_str());
  if (!f) {
    std::cerr << "[Octree write] could not open file: " << fullName << std::endl;
  }
  std::cout << "Writing Octree(" << leafs.size() << " leafs) to vtk file...";

  f << "# vtk DataFile Version 2.0" << std::endl;
  f << "Octree" << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET UNSTRUCTURED_GRID" << std::endl;
  std::stringstream points;
  std::stringstream cells;
  std::stringstream cell_types;
  // std::stringstream point_data;
  std::stringstream cell_data;
  std::stringstream cell_leaf;
  std::stringstream cell_boundary;

  points << "POINTS " << leafs.size() * 8 << " float" << std::endl;
  cells << "CELLS " << leafs.size() << " " << leafs.size() * 9 << std::endl;
  cell_types << "CELL_TYPES " << leafs.size() << std::endl;
  cell_data << "CELL_DATA " << leafs.size() << std::endl;
  cell_data << "SCALARS insideout int" << std::endl;
  cell_data << "LOOKUP_TABLE default" << std::endl;
  cell_leaf << "SCALARS leaf int" << std::endl;
  cell_leaf << "LOOKUP_TABLE default" << std::endl;
  cell_boundary << "SCALARS boundary int" << std::endl;
  cell_boundary << "LOOKUP_TABLE default" << std::endl;

  Vector<T, 3> center;
  Vector<T, 3> pt;

  T rad;
  int i = 0;
  // /// Gets centerpoint
  // inline const Vector<T, 3>& getCenter() const { return _center; };
  // /// Gets radius
  // inline const T getRadius() const { return _radius; };
  for (; it != leafs.end(); ++it) {
    center = (*it)->getCenter();
    rad = (*it)->getRadius();

    pt[0] = center[0] - rad;
    pt[1] = center[1] - rad;
    pt[2] = center[2] - rad;
    points << pt[0] << " " << pt[1] << " " << pt[2] << " ";

    pt[0] = center[0] + rad;
    pt[1] = center[1] - rad;
    pt[2] = center[2] - rad;
    points << pt[0] << " " << pt[1] << " " << pt[2] << " ";

    pt[0] = center[0] - rad;
    pt[1] = center[1] + rad;
    pt[2] = center[2] - rad;
    points << pt[0] << " " << pt[1] << " " << pt[2] << " ";

    pt[0] = center[0] + rad;
    pt[1] = center[1] + rad;
    pt[2] = center[2] - rad;
    points << pt[0] << " " << pt[1] << " " << pt[2] << " ";

    pt[0] = center[0] - rad;
    pt[1] = center[1] - rad;
    pt[2] = center[2] + rad;
    points << pt[0] << " " << pt[1] << " " << pt[2] << " ";

    pt[0] = center[0] + rad;
    pt[1] = center[1] - rad;
    pt[2] = center[2] + rad;
    points << pt[0] << " " << pt[1] << " " << pt[2] << " ";

    pt[0] = center[0] - rad;
    pt[1] = center[1] + rad;
    pt[2] = center[2] + rad;
    points << pt[0] << " " << pt[1] << " " << pt[2] << " ";

    pt[0] = center[0] + rad;
    pt[1] = center[1] + rad;
    pt[2] = center[2] + rad;
    points << pt[0] << " " << pt[1] << " " << pt[2] << " " << std::endl;

    cells << "8 ";
    for (int j = 0; j < 8; j++) {
      cells << i + j << " ";
    }
    i += 8;
    cells << std::endl;

    cell_types << 11 << std::endl;

    cell_data << (*it)->getInside() << " " << std::endl;
    cell_leaf << (*it)->getMaxdepth() << " " << std::endl;
    cell_boundary << (*it)->getBoundaryNode() << " " << std::endl;
  }

  f << points.str() << cells.str() << cell_types.str() << cell_data.str()
    << cell_leaf.str() << cell_boundary.str();

  // f << "POINT_DATA 0\nCELL_DATA 0\n" << std::endl;
  f.close();
  std::cout << " Done!" << std::endl;
}

template <typename T>
bool Octree<T>::closestIntersectionSphere(const Vector<T, 3>& pt, const T& rad,
                                          const Vector<T, 3>& direction, Vector<T, 3>& q,
                                          T& a, Triangle<T>& tri) {
  a = std::numeric_limits<T>::infinity();
  T alpha = T();
  std::vector<T> qtmp(3, T());
  bool found = false;
  for (unsigned int k = 0; k < _triangles.size(); ++k) {
    if (_mesh->getTri(_triangles[k])
          .testMovingSphereIntersect(pt, rad, direction, qtmp, alpha)) {
      if (alpha < a) {
        a = alpha;
        q = qtmp;
        found = true;
        tri = _mesh->getTri(_triangles[k]);
      }
    }
  }
  return found;
}

template <typename T>
bool Octree<T>::closestIntersection(const Vector<T, 3>& pt, const Vector<T, 3>& direction,
                                    Vector<T, 3>& q, T& a, Triangle<T>& tri, const T& rad,
                                    bool print) {
  a = std::numeric_limits<T>::infinity();
  T alpha = T();
  Vector<T, 3> qtmp;
  bool found = false;

  for (unsigned int k = 0; k < _triangles.size(); ++k) {
    if (_mesh->getTri(_triangles[k]).testRayIntersect(pt, direction, qtmp, alpha, rad)) {
      if (print) {
        std::cout << "Found intersection!" << std::endl;
      }
      if (alpha < a) {
        a = alpha;
        q = qtmp;
        found = true;
        tri = _mesh->getTri(_triangles[k]);
      }
    }
  }
  //  std::cout << a << std::endl;
  return found;
}

template <typename T>
bool Octree<T>::closestIntersection(const Vector<T, 3>& pt, const Vector<T, 3>& direction,
                                    Vector<T, 3>& q, T& a) {
  Triangle<T> tri;
  return closestIntersection(pt, direction, q, a, tri, 0.);
}

template <typename T>
void Octree<T>::intersectRayNode(const Vector<T, 3>& pt, const Vector<T, 3>& dir,
                                 Vector<T, 3>& s) {
  T t, d;
  s = s * T();
  // Plane Normals outside

  // calculates the intersection point of the ray with the plane of the bounding
  // box that is facing towards the ray.
  // then checks if the intersection point is inside the bounding box

  if (dir[0] > 0.) {
    // n = {1, 0, 0}
    d = _center[0] + _radius;
    t = (d - pt[0]) / dir[0];
    s = pt + t * dir;
    if (std::fabs(s[1] - _center[1]) < _radius &&
        std::fabs(s[2] - _center[2]) < _radius) {
      return;
    }
  } else if (dir[0] < 0.) {
    // n = {-1, 0, 0}
    d = _center[0] - _radius;
    t = (d - pt[0]) / dir[0];
    s = pt + t * dir;
    if (std::fabs(s[1] - _center[1]) < _radius &&
        std::fabs(s[2] - _center[2]) < _radius) {
      return;
    }
  }

  if (dir[1] > 0.) {
    d = _center[1] + _radius;
    t = (d - pt[1]) / dir[1];
    s = pt + t * dir;
    if (std::fabs(s[0] - _center[0]) < _radius &&
        std::fabs(s[2] - _center[2]) < _radius) {
      return;
    }
  } else if (dir[1] < 0.) {
    // n = {0, 0, -1}
    d = _center[1] - _radius;
    t = (d - pt[1]) / dir[1];
    s = pt + t * dir;
    if (std::fabs(s[0] - _center[0]) < _radius &&
        std::fabs(s[2] - _center[2]) < _radius) {
      return;
    }
  }

  if (dir[2] > 0.) {
    // n = {0, 0, 1}
    d = _center[2] + _radius;
    t = (d - pt[2]) / dir[2];
    s = pt + t * dir;
    if (std::fabs(s[0] - _center[0]) < _radius &&
        std::fabs(s[1] - _center[1]) < _radius) {
      return;
    }
  } else if (dir[2] < 0.) {
    // n = {0, 0, -1}
    d = _center[2] - _radius;
    t = (d - pt[2]) / dir[2];
    s = pt + t * dir;
    if (std::fabs(s[0] - _center[0]) < _radius &&
        std::fabs(s[1] - _center[1]) < _radius) {
      return;
    }
  }
}

template <typename T>
void Octree<T>::print(std::vector<Octree<T>*>& leafs) {
  MPI_RANK(0)
  std::cout << "[Octree]: " << std::endl;
  std::cout << "TreeRadius = " << _radius << "; Center = (" << _center[0] << ","
            << _center[1] << "," << _center[2] << "); Leafs = " << leafs.size()
            << "; minRadius = " << getMinRadius(leafs) << std::endl;
}

template <typename T>
void Octree<T>::trianglesOnLine(const Vector<T, 3>& pt1, const Vector<T, 3>& pt2,
                                std::set<unsigned int>& tris) {
  tris.clear();
  std::vector<T> line = pt2 - pt1;
  std::vector<T> s = pt1;
  T lineNorm2 = line[0] * line[0] + line[1] * line[1] + line[2] * line[2];
  T dist2 = T();
  Octree<T>* node = NULL;
  int it = 0;
  while (dist2 < lineNorm2 && it < 50) {
    node = find(s);
    tris.insert(node->_triangles.begin(), node->_triangles.end());
    node->intersectRayNode(s, line, s);
    for (int i = 0; i < 3; i++) {
      s[i] = s[i] + line[i] * _radius * 0.001 /* *node->getRadius()*/;
    }
    it++;
    dist2 = (pt1[0] - s[0]) * (pt1[0] - s[0]) + (pt1[1] - s[1]) * (pt1[1] - s[1]) +
            (pt1[2] - s[2]) * (pt1[2] - s[2]);
    // dist2 = GetDist2(pt1, s);
  }
}
