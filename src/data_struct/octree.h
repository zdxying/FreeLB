/* This file is part of FreeLB, modified from OpenLB's octree.h, with the following copyright notice:
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


// octree.h

#pragma once

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "data_struct/Vector.h"
#include "utils/util.h"

template <typename T>
class StlMesh;

template <typename T>
struct Triangle;

template <typename T>
class Octree {
 protected:
  /// _vector _triangles contains number of triangles
  std::vector<unsigned int> _triangles;
  Vector<T, 3> _center;
  T _radius;
  StlMesh<T>* _mesh;
  short _maxDepth;
  bool _isLeaf;
  bool _boundaryNode;
  bool _inside;
  Octree<T>* _parent;
  Octree<T>** _child;
  //
  void findTriangles(T overlap = 0.);
  bool AABBTri(const Triangle<T>& tri, T overlap = 0.);

 public:
  /*
   * Constructs Octree containing triangles of an StlMesh.
   * \param center Centerpoint
   * \param rad Radius
   * \param mesh StlMesh
   * \param maxDepth Maximal depth of tree
   * \param overlap Triangles within rad+overlap are added to this Octree
   */
  Octree(Vector<T, 3> center, T rad, StlMesh<T>* mesh, short maxDepth,
         T overlap = 0., Octree<T>* parent = nullptr);
  /// Destructor destructs
  ~Octree();
  /// Find the node containing the first param with remaining maxDepth
  Octree<T>* find(const Vector<T, 3>&, const int& maxDepth = 0);
  /// Write Octree
  void write(const Vector<T, 3>& pt, const std::string no);
  /// Write Octree
  void write(const int, const std::string);
  /// Write Octree
  void write(const std::string);
  /// Test intersection of ray with all triangles in Octree
  /// returns number of intersections
  int testIntersection(const Vector<T, 3>& pt, const Vector<T, 3>& dir,
                       bool print = false);
  /// Test intersection of ray with all triangles in Octree
  /// q contains point of closest intersection to pt in direction direction
  bool closestIntersection(const Vector<T, 3>& pt,
                           const Vector<T, 3>& direction, Vector<T, 3>& q,
                           T& a);
  /// Test intersection of ray with all triangles in Octree
  /// q contains point of closest intersection to pt in direction direction
  /// tri contains triangle with closest intersection
  bool closestIntersection(const Vector<T, 3>& pt,
                           const Vector<T, 3>& direction, Vector<T, 3>& q, T& a,
                           Triangle<T>& tri, const T& rad = 0.,
                           bool print = false);
  /// Test intersection of sphere moving along ray with radius rad
  /// q contains point of closest intersection to pt in direction direction
  /// tri contains triangle with closest intersection
  bool closestIntersectionSphere(const Vector<T, 3>& pt, const T& rad,
                                 const Vector<T, 3>& direction, Vector<T, 3>& q,
                                 T& a, Triangle<T>& tri);
  /// It's complicated. Computes intersections of a ray with triangles inside
  /// this Octree. Sets _inside depending on value of rayInside and changes
  /// rayInside depending on the number of intersections. Also takes into
  /// account if the intersections happen before or after the center.
  void checkRay(const Vector<T, 3>& pt, const Vector<T, 3>& dir,
                unsigned short& rayInside);
  /// Computes intersection of ray with Octree boundaries
  void intersectRayNode(const Vector<T, 3>& pt, const Vector<T, 3>& dir,
                        Vector<T, 3>& s);
  /// Computes all centerpoints of Octree
  void getCenterpoints(std::vector<std::vector<T> >& pts);
  /// Collectes all leafs
  void getLeafs(std::vector<Octree<T>*>& pts);
  /// Return status of _isLeaf;
  bool isLeaf() { return _isLeaf; }
  // get min radius of leafs
  T getMinRadius();
  T getMinRadius(std::vector<Octree<T>*> &leafs);

  /// Sets Inside
  inline void setInside(bool ins) { _inside = ins; };
  /// Gets Inside
  inline bool getInside() { return _inside; };
  /// Gets _boundarNode
  inline bool getBoundaryNode() { return _boundaryNode; };
  /// Gets Maxdepth
  inline int getMaxdepth() const { return _maxDepth; };
  /// Gets numbers of triangles contained by this Octree
  inline const std::vector<unsigned int>& getTriangles() const {
    return _triangles;
  };
  /// Gets centerpoint
  inline const Vector<T, 3>& getCenter() const { return _center; };
  /// Gets radius
  inline const T getRadius() const { return _radius; };
  /// Prints console output
  void print(std::vector<Octree<T>*> &leafs);
  /// Returns set of indices of all triangles in nodes containing a line.
  void trianglesOnLine(const Vector<T, 3>& pt1, const Vector<T, 3>& pt2,
                       std::set<unsigned int>& tris);
  /// Returns reference to _mesh
  inline StlMesh<T>* getMesh() { return _mesh; }

  ////////////
  // get current depth
  int getDepth();
  //
};
