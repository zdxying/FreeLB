/* This file is part of FreeLB, modified from OpenLB's stlReader.h, with the following copyright notice:
 *
 * // start of the original OpenLB's copyright notice
 * 
 * This file is part of the OpenLB library
 *
 *  Copyright (C) 2010-2015 Thomas Henn, Mathias J. Krause
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

// stlreader.h

// this file read stl file, store and generate a mesh, and voxelise the mesh

#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "data_struct/Vector.h"
#include "data_struct/octree.h"
#include "utils/util.h"

template <typename T>
struct Triangle {
  // 3 verteices of the triangle
  // this could be initialized like: vertex(3) or vertex(3, Vector<T, 3>(value))
  std::vector<Vector<T, 3>> vertex;

  // normal vector of the triangle
  Vector<T, 3> normal;

  /// variables explained in
  /// http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf
  /// page 7-12
  /// precomputed for speedup
  Vector<T, 3> uBeta, uGamma;
  T d, kBeta, kGamma;

  Triangle()
      : vertex(3), normal(), uBeta(), uGamma(), d(0), kBeta(0), kGamma(0) {}

  // Initialize triangle and precompute
  void Init();
  /// Return write access to normal
  inline Vector<T, 3> &getNormal() { return normal; }
  /// Return read access to normal
  inline const Vector<T, 3> &getNormal() const { return normal; }
  /// Returns center
  Vector<T, 3> getCenter();
  /// Returns Pt0-Pt1
  std::vector<T> getE0();
  /// Returns Pt0-Pt2
  std::vector<T> getE1();
  // Check whether a point is inside a triangle
  bool IsInside(const Vector<T, 3> &pt) const;
  /// Compute intersection between Ray and set of triangles; returns true if
  /// intersection is found
  bool testRayIntersect(const Vector<T, 3> &pt, const Vector<T, 3> &dir,
                        Vector<T, 3> &q, T &alpha, const T &rad = T(),
                        bool print = false);
  Vector<T, 3> closestPtPointTriangle(const Vector<T, 3> &pt) const;
};

// this class read stl file and store the mesh in std::vector<Triangle<T>>
// Triangles the min and max points of axis aligned bounding box coordinate is
// computed when reading the stl file, note that the min and max points are not
// neccessarily the actual points of the mesh
template <typename T>
class StlMesh {
 private:
  std::string _fName;
  // Vector of Triangles
  std::vector<Triangle<T>> Triangles;
  /// Min and Max points of axis aligned bounding box coordinate in SI units
  Vector<T, 3> _min, _max;
  /// largest squared length of edge of all triangles
  T _maxDist2;

 public:
  // stl size: vertex *= stlSize
  T _stlSize;
  // get the min and max vertex if i = 0
  inline void Getminmax0(Triangle<T> &tri);
  // get the min and max vertex if i != 0
  inline void Getminmax(Triangle<T> &tri);
  // construct from stl file: vertex *= stlSize
  StlMesh(std::string filename, T stlSize = T(1));
  /// Returns reference to a triangle
  inline Triangle<T> &getTri(unsigned int i) { return Triangles[i]; }
  /// Returns reference to all triangles
  inline std::vector<Triangle<T>> &getTriangles() { return Triangles; }
  /// Returns number of triangles
  inline unsigned int triangleSize() const { return Triangles.size(); }
  /// Returns _min
  inline Vector<T, 3> &getMin() { return _min; }
  inline const Vector<T, 3> &getMin() const { return _min; }
  /// Returns _max
  inline Vector<T, 3> &getMax() { return _max; }
  inline const Vector<T, 3> &getMax() const { return _max; }
  // return AABB
  inline AABB<T, 3> getAABB() const { return AABB<T, 3>{_min, _max}; }
  // return Max - Min
  inline const Vector<T, 3> getMax_Min() const { return _max - _min; }
  /// Returns maxDist squared
  inline float maxDist2() const { return _maxDist2; }
  /// Prints console output
  void print(bool full = false);
  /// Writes STL mesh in Si units
  void write(std::string fName);
  /// Compute intersection between Ray and set of triangles; returns true if
  /// intersection is found
  bool testRayIntersect(const std::set<unsigned int> &tris,
                        const Vector<T, 3> &pt, const Vector<T, 3> &dir,
                        Vector<T, 3> &q, T &alpha);
};

template <typename T>
class StlReader {
 private:
  /// Size of the smallest voxel
  T _voxelSize;
  /// Factor to get Si unit (m), i.e. "0.001" means mm
  T _stlSize;
  /// Overlap increases Octree radius by _overlap
  T _overlap;
  /// Pointer to tree
  Octree<T> *_tree;
  /// The filename
  std::string _fName;
  /// The mesh
  StlMesh<T> _mesh;
  /// Variable for output
  bool _verbose;
  /// Method to indicate inside nodes
  Vector<T, 3> _myMin, _myMax;

 public:
  /// Constructor
  // vertex *= stlSize
  StlReader(std::string fName, T voxelSize, T stlSize = 1, int method = 2,
            bool verbose = false, T overlap = 0., T max = 0.);
  /// Destructor
  ~StlReader() { delete _tree; }

  // get _stlSize
  T getStlSize() const { return _stlSize; }
  // get _voxelSize
  T getVoxelSize() const { return _voxelSize; }

  //////// method
  /*
   *  Old indicate function (slower, more stable)
   *  Define three rays (X-, Y-, Z-direction) for each leaf and count
   * intersections with STL for each ray. Odd number of intersection means
   * inside (Majority vote).
   */
  void indicate1();
  /*
   *  New indicate function (faster, less stable)
   *  Define ray in Z-direction for each Voxel in XY-layer. Indicate all nodes
   * on the fly.
   */
  void indicate2();
  /*
   *  New indicate function (faster, less stable)
   *  Define ray in X-direction for each Voxel in YZ-layer. Indicate all nodes
   * on the fly.
   */
  void indicate2_Xray();
  /*
   *  New indicate function (faster, less stable)
   *  Define ray in Y-direction for each Voxel in XZ-layer. Indicate all nodes
   * on the fly.
   */
  void indicate2_Yray();
  /*
   *  Double ray approach: two times (X-, Y-, Z-direction) for each leaf.
   *  Could be use to deal with double layer triangles and face intersections.
   */
  void indicate3();
  /////////
  /// Returns whether node is inside or not.
  bool operator()(bool output[], const T input[]);

  /// Computes distance to closest triangle intersection
  bool distance(T &distance, const Vector<T, 3> &origin,
                const Vector<T, 3> &direction, int iC = -1);

  /// Finds normal for points on the surface (do not use for points that aren't
  /// on the surface!)
  Vector<T, 3> findNormalOnSurface(const Vector<T, 3> &pt);

  /// Finds surface normal
  Vector<T, 3> evalSurfaceNormal(const Vector<T, 3> &origin);

  /// Computes signed distance to closest triangle in direction of the surface
  /// normal
  T signedDistance(const Vector<T, 3> &input);

  /// Finds and returns normal of the closest surface (triangle)
  Vector<T, 3> surfaceNormal(const Vector<T, 3> &pos, const T meshSize = 0) {
    return evalSurfaceNormal(pos);
  }

  /// Prints console output
  void print();

  /// Writes STL mesh in Si units
  void writeSTL(std::string stlName = "") {
    if (stlName == "") {
      _mesh.write(_fName);
    } else {
      _mesh.write(stlName);
    }
  }

  /// Writes Octree
  void writeOctree() { _tree->write(_fName); }

  /// Rearranges normals of triangles to point outside of geometry
  void setNormalsOutside();

  /// Every octree leaf intersected by the STL will be part of the inside nodes.
  /// Artificially enlarges all details that would otherwise be cut off by the
  /// voxelSize.
  void setBoundaryInsideNodes();
  /// Returns tree
  inline Octree<T> *getTree() const { return _tree; };
  /// Returns reference of mesh
  inline StlMesh<T> &getMesh() { return _mesh; };
  inline const StlMesh<T> &getMesh() const { return _mesh; };
  //
  inline const Vector<T, 3> &getMin() const { return _mesh.getMin(); };
  inline const Vector<T, 3> &getMax() const { return _mesh.getMax(); };
};
