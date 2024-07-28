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

// quadtree.h
#pragma once

#include <omp.h>

#include "data_struct/Vector.h"
#include "head.h"


// namespace tree {

// }  // namespace tree

template <typename T, unsigned int D>
class AbstractTree {
 protected:
  // vox
  std::array<Vector<T, D>*, D> voxels;
  // coarse grid id of uppper level
  std::size_t _id;

 public:
  template <typename... Args>
  AbstractTree(std::size_t id, Args... args) : _id(id), voxels{args...} {}
  std::size_t getId() const { return _id; }
  Vector<T, D>& getVoxel(int i) { return *voxels[i]; }
  const Vector<T, D>& getVoxel(int i) const { return *voxels[i]; }
  std::array<Vector<T, D>*, D>& getVoxels() { return voxels; }
  const std::array<Vector<T, D>*, D>& getVoxels() const { return voxels; }
};

template <typename T, unsigned int D>
class BasicTree : public Vector<T, D> {
 protected:
  // Vector<T, D> _center;
  T _radius;
  std::size_t _id;
  // voexl flag
  std::uint8_t _flag;

 public:
  BasicTree(const Vector<T, D>& center, T radius, std::size_t id, std::uint8_t flag)
      : Vector<T, D>(center), _radius(radius), _id(id), _flag(flag) {}
  Vector<T, D>& getCenter() { return *this; }
  const Vector<T, D>& getCenter() const { return *this; }
  T getRadius() const { return _radius; }
  std::size_t getId() const { return _id; }
  std::uint8_t getFlag() const { return _flag; }
  void setFlag(std::uint8_t flag) { _flag = flag; }
};

// quadtree without neighbors
template <typename T>
class QuadTree : public BasicTree<T, 2> {
 protected:
  bool _isLeaf;
  // refinement level: 1, 2, ...
  std::uint8_t _level;

  // parent and children nodes
  QuadTree<T>* _parent;
  /*
          2 | 3
          0 | 1
  */
  QuadTree<T>** _child;

 public:
  QuadTree(std::size_t id, const Vector<T, 2>& center, std::uint8_t level, std::uint8_t maxlevel,
           std::uint8_t flag, QuadTree<T>* parent = nullptr, T voxsize = T(1))
      : BasicTree<T, 2>(center, pow(T(0.5), level + 1) * voxsize, id, flag), _level(level),
        _parent(parent) {
    if (_level < maxlevel) {
      _isLeaf = false;
      _child = new QuadTree<T>*[4];
      Vector<T, 2> tmpCenter = this->getCenter();
      T halfradius = this->_radius * T(0.5);
      tmpCenter[0] = this->getCenter()[0] - halfradius;
      tmpCenter[1] = this->getCenter()[1] - halfradius;
      _child[0] = new QuadTree<T>(id, tmpCenter, _level + 1, maxlevel, flag, this, voxsize);
      tmpCenter[0] = this->getCenter()[0] + halfradius;
      tmpCenter[1] = this->getCenter()[1] - halfradius;
      _child[1] = new QuadTree<T>(id, tmpCenter, _level + 1, maxlevel, flag, this, voxsize);
      tmpCenter[0] = this->getCenter()[0] - halfradius;
      tmpCenter[1] = this->getCenter()[1] + halfradius;
      _child[2] = new QuadTree<T>(id, tmpCenter, _level + 1, maxlevel, flag, this, voxsize);
      tmpCenter[0] = this->getCenter()[0] + halfradius;
      tmpCenter[1] = this->getCenter()[1] + halfradius;
      _child[3] = new QuadTree<T>(id, tmpCenter, _level + 1, maxlevel, flag, this, voxsize);
    } else {
      _isLeaf = true;
    }
  }
  QuadTree(std::size_t id, const Vector<T, 2>& center, std::uint8_t level, std::uint8_t maxlevel,
           std::uint8_t flag, T rad, QuadTree<T>* parent)
      : BasicTree<T, 2>(center, rad, id, flag), _level(level), _parent(parent) {
    if (_level < maxlevel) {
      _isLeaf = false;
      _child = new QuadTree<T>*[4];
      Vector<T, 2> tmpCenter = this->getCenter();
      T halfradius = this->_radius * T(0.5);
      tmpCenter[0] = this->getCenter()[0] - halfradius;
      tmpCenter[1] = this->getCenter()[1] - halfradius;
      _child[0] = new QuadTree<T>(id, tmpCenter, _level + 1, maxlevel, flag, halfradius, this);
      tmpCenter[0] = this->getCenter()[0] + halfradius;
      tmpCenter[1] = this->getCenter()[1] - halfradius;
      _child[1] = new QuadTree<T>(id, tmpCenter, _level + 1, maxlevel, flag, halfradius, this);
      tmpCenter[0] = this->getCenter()[0] - halfradius;
      tmpCenter[1] = this->getCenter()[1] + halfradius;
      _child[2] = new QuadTree<T>(id, tmpCenter, _level + 1, maxlevel, flag, halfradius, this);
      tmpCenter[0] = this->getCenter()[0] + halfradius;
      tmpCenter[1] = this->getCenter()[1] + halfradius;
      _child[3] = new QuadTree<T>(id, tmpCenter, _level + 1, maxlevel, flag, halfradius, this);
    } else {
      _isLeaf = true;
    }
  }
  QuadTree<T>* getParent() { return _parent; }
  QuadTree<T>** getChildren() { return _child; }
  QuadTree<T>* getChild(int i) { return _child[i]; }

  inline bool isInside(const Vector<T, 2>& pt) {
    return (std::abs(this->getCenter()[0] - pt[0]) < this->_radius &&
            std::abs(this->getCenter()[1] - pt[1]) < this->_radius);
  }
  void getLeafs(std::vector<BasicTree<T, 2>*>& leafs) {
    if (_isLeaf) {
      leafs.push_back(this);
    } else {
      _child[0]->getLeafs(leafs);
      _child[1]->getLeafs(leafs);
      _child[2]->getLeafs(leafs);
      _child[3]->getLeafs(leafs);
    }
  }
  void Refine() {
    _isLeaf = false;
    _child = new QuadTree<T>*[4];
    Vector<T, 2> tmpCenter = this->getCenter();
    T halfradius = this->_radius * T(0.5);
    tmpCenter[0] = this->getCenter()[0] - halfradius;
    tmpCenter[1] = this->getCenter()[1] - halfradius;
    _child[0] = new QuadTree<T>(this->getId(), tmpCenter, _level + 1, _level + 1, this->_flag,
                                halfradius, this);
    tmpCenter[0] = this->getCenter()[0] + halfradius;
    tmpCenter[1] = this->getCenter()[1] - halfradius;
    _child[1] = new QuadTree<T>(this->getId(), tmpCenter, _level + 1, _level + 1, this->_flag,
                                halfradius, this);
    tmpCenter[0] = this->getCenter()[0] - halfradius;
    tmpCenter[1] = this->getCenter()[1] + halfradius;
    _child[2] = new QuadTree<T>(this->getId(), tmpCenter, _level + 1, _level + 1, this->_flag,
                                halfradius, this);
    tmpCenter[0] = this->getCenter()[0] + halfradius;
    tmpCenter[1] = this->getCenter()[1] + halfradius;
    _child[3] = new QuadTree<T>(this->getId(), tmpCenter, _level + 1, _level + 1, this->_flag,
                                halfradius, this);
  }

  // lambda functions
  template <typename Func>
  void forEachLeaf(const Func& func) {
    if (_isLeaf) {
      func(this);
    } else {
      _child[0]->forEachLeaf(func);
      _child[1]->forEachLeaf(func);
      _child[2]->forEachLeaf(func);
      _child[3]->forEachLeaf(func);
    }
  }

  BasicTree<T, 2>* getAbstractTree() { return this; }
  std::uint8_t getLevel() const { return _level; }
};

// // quadtree with neighbors
// namespace quadtree {
// template <typename T, typename LatSet>
// class QuadTree : public BasicTree<T, 2> {
//  protected:
//   bool _isLeaf;
//   // refinement level: 1, 2, ...
//   short _level;
//   // voexl flag
//   std::uint8_t _flag;
//   // coarse cell index
//   int _id;
//   // fine cell index
//   int id;
//   // parent and children nodes
//   QuadTree<T, LatSet>* _parent;
//   /*
//           2 | 3
//           0 | 1
//   */
//   QuadTree<T, LatSet>** _child;

//   // neighbors
//   std::array<QuadTree<T, LatSet>*, LatSet::q> _neighbors;
//   // neighbor count
//   short _neighbor_count;

//  public:
//   QuadTree(int Cid, const Vector<T, 2>& center, short level, std::uint8_t
//   flag,
//            bool isLeaf = true, QuadTree<T, LatSet>* parent = nullptr)
//       : _id(Cid),
//         id(0),
//         BasicTree<T, 2>(center, pow(T(0.5), level + 1)),
//         _level(level),
//         _flag(flag),
//         _parent(parent),
//         _isLeaf(isLeaf) {
//     if (!_isLeaf) {
//       _child = new QuadTree<T, LatSet>*[4];
//       Vector<T, 2> tmpCenter = this->_center;
//       T halfradius = this->_radius * T(0.5);
//       tmpCenter[0] = this->_center[0] - halfradius;
//       tmpCenter[1] = this->_center[1] - halfradius;
//       _child[0] = new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag,
//                                           true, this);
//       tmpCenter[0] = this->_center[0] + halfradius;
//       tmpCenter[1] = this->_center[1] - halfradius;
//       _child[1] = new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag,
//                                           true, this);
//       tmpCenter[0] = this->_center[0] - halfradius;
//       tmpCenter[1] = this->_center[1] + halfradius;
//       _child[2] = new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag,
//                                           true, this);
//       tmpCenter[0] = this->_center[0] + halfradius;
//       tmpCenter[1] = this->_center[1] + halfradius;
//       _child[3] = new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag,
//                                           true, this);
//     }
//   }
//   QuadTree(int Cid, const Vector<T, 2>& center, short level, short maxlevel,
//            std::uint8_t flag, QuadTree<T, LatSet>* parent = nullptr)
//       : _id(Cid),
//         id(0),
//         BasicTree<T, 2>(center, pow(T(0.5), level + 1)),
//         _level(level),
//         _flag(flag),
//         _parent(parent) {
//     if (level < maxlevel) {
//       _isLeaf = false;
//       _child = new QuadTree<T, LatSet>*[4];
//       Vector<T, 2> tmpCenter = this->_center;
//       T halfradius = this->_radius * T(0.5);
//       tmpCenter[0] = this->_center[0] - halfradius;
//       tmpCenter[1] = this->_center[1] - halfradius;
//       _child[0] = new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag,
//                                           true, this);
//       tmpCenter[0] = this->_center[0] + halfradius;
//       tmpCenter[1] = this->_center[1] - halfradius;
//       _child[1] = new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag,
//                                           true, this);
//       tmpCenter[0] = this->_center[0] - halfradius;
//       tmpCenter[1] = this->_center[1] + halfradius;
//       _child[2] = new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag,
//                                           true, this);
//       tmpCenter[0] = this->_center[0] + halfradius;
//       tmpCenter[1] = this->_center[1] + halfradius;
//       _child[3] = new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag,
//                                           true, this);
//     } else {
//       _isLeaf = true;
//     }
//   }
//   ~QuadTree() {
//     if (!_isLeaf) {
//       for (int i = 0; i < 4; ++i) delete _child[i];
//     } else {
//       delete[] _child;
//     }
//   }
//   int getCId() const { return _id; }
//   int getFId() const { return id; }
//   bool isLeaf() { return _isLeaf; }
//   QuadTree<T, LatSet>* getParent() { return _parent; }
//   QuadTree<T, LatSet>** getChildren() { return _child; }
//   QuadTree<T, LatSet>* getChild(int i) { return _child[i]; }

//   // find leaf node
//   QuadTree<T, LatSet>* findLeaf(const Vector<T, 2>& pt) {
//     if (isInside(pt)) {
//       if (_isLeaf) {
//         return this;
//       } else {
//         if (pt[0] < this->_center[0]) {
//           if (pt[1] < this->_center[1]) {
//             return _child[0]->findLeaf(pt);
//           } else {
//             return _child[2]->findLeaf(pt);
//           }
//         } else {
//           if (pt[1] < this->_center[1]) {
//             return _child[1]->findLeaf(pt);
//           } else {
//             return _child[3]->findLeaf(pt);
//           }
//         }
//       }
//     } else {
//       return nullptr;
//     }
//   }
//   void getLeafs(std::vector<QuadTree<T, LatSet>*>& leafs) {
//     if (_isLeaf) {
//       leafs.push_back(this);
//     } else {
//       _child[0]->getLeafs(leafs);
//       _child[1]->getLeafs(leafs);
//       _child[2]->getLeafs(leafs);
//       _child[3]->getLeafs(leafs);
//     }
//   }
//   // find child node within the current node
//   QuadTree<T, LatSet>* findChild(const Vector<T, 2>& pt) {
//     if (isInside(pt)) {
//       if (pt[0] < this->_center[0]) {
//         if (pt[1] < this->_center[1]) {
//           return _child[0];
//         } else {
//           return _child[2];
//         }
//       } else {
//         if (pt[1] < this->_center[1]) {
//           return _child[1];
//         } else {
//           return _child[3];
//         }
//       }
//     } else {
//       return nullptr;
//     }
//   }
//   inline bool isInside(const Vector<T, 2>& pt) {
//     return (std::abs(this->_center[0] - pt[0]) < this->_radius &&
//             std::abs(this->_center[1] - pt[1]) < this->_radius);
//   }
//   void SetFId(int fid) { id = fid; }
//   inline void setNeighbor(QuadTree<T, LatSet>* neighbor, int i) {
//     if (_neighbors[i] == nullptr && neighbor != nullptr) ++_neighbor_count;
//     _neighbors[i] = neighbor;
//   }
//   void SetupChildNbr() {
//     if constexpr (LatSet::q == 5) {
//       // {0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}}
//       // _child[0]
//       _child[0].setNeighbor(_child[1], 1);
//       _child[0].setNeighbor(_child[2], 2);
//       if (_neighbors[3] != nullptr && !_neighbors[3]->isLeaf())
//         _child[0].setNeighbor(_neighbors[3]->_child[1], 3);
//       if (_neighbors[4] != nullptr && !_neighbors[4]->isLeaf())
//         _child[0].setNeighbor(_neighbors[4]->_child[2], 4);

//       // _child[1]
//       _child[1].setNeighbor(_child[3], 2);
//       _child[1].setNeighbor(_child[0], 3);
//       if (_neighbors[1] != nullptr && !_neighbors[1]->isLeaf())
//         _child[1].setNeighbor(_neighbors[1]->_child[0], 1);
//       if (_neighbors[4] != nullptr && !_neighbors[4]->isLeaf())
//         _child[1].setNeighbor(_neighbors[4]->_child[3], 4);

//       // _child[2]
//       _child[2].setNeighbor(_child[3], 1);
//       _child[2].setNeighbor(_child[0], 4);
//       if (_neighbors[2] != nullptr && !_neighbors[2]->isLeaf())
//         _child[2].setNeighbor(_neighbors[2]->child[0], 2);
//       if (_neighbors[3] != nullptr && !_neighbors[3]->isLeaf())
//         _child[2].setNeighbor(_neighbors[3]->_child[3], 3);

//       // _child[3]
//       _child[3].setNeighbor(_child[2], 3);
//       _child[3].setNeighbor(_child[1], 4);
//       if (_neighbors[1] != nullptr && !_neighbors[1]->isLeaf())
//         _child[3].setNeighbor(_neighbors[1]->_child[2], 1);
//       if (_neighbors[2] != nullptr && !_neighbors[2]->isLeaf())
//         _child[3].setNeighbor(_neighbors[2]->_child[1], 2);

//     } else if constexpr (LatSet::q == 9) {
//       // {0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}}
//       // _child[0]
//       _child[0].setNeighbor(_child[1], 1);
//       _child[0].setNeighbor(_child[2], 2);
//       _child[0].setNeighbor(_child[3], 5);
//       if (_neighbors[3] != nullptr && !_neighbors[3]->isLeaf()) {
//         _child[0].setNeighbor(_neighbors[3]->_child[1], 3);
//         _child[0].setNeighbor(_neighbors[3]->_child[3], 6);
//       }
//       if (_neighbors[4] != nullptr && !_neighbors[4]->isLeaf()) {
//         _child[0].setNeighbor(_neighbors[4]->_child[2], 4);
//         _child[0].setNeighbor(_neighbors[4]->_child[3], 8);
//       }
//       if (_neighbors[7] != nullptr && !_neighbors[7]->isLeaf())
//         _child[0].setNeighbor(_neighbors[7]->_child[3], 7);

//       // _child[1]
//       _child[1].setNeighbor(_child[3], 2);
//       _child[1].setNeighbor(_child[0], 3);
//       _child[1].setNeighbor(_child[2], 6);
//       if (_neighbors[1] != nullptr && !_neighbors[1]->isLeaf()) {
//         _child[1].setNeighbor(_neighbors[1]->_child[0], 1);
//         _child[1].setNeighbor(_neighbors[1]->_child[2], 5);
//       }
//       if (_neighbors[4] != nullptr && !_neighbors[4]->isLeaf()) {
//         _child[1].setNeighbor(_neighbors[4]->_child[3], 4);
//         _child[1].setNeighbor(_neighbors[4]->_child[2], 7);
//       }
//       if (_neighbors[8] != nullptr && !_neighbors[8]->isLeaf())
//         _child[1].setNeighbor(_neighbors[8]->_child[2], 8);

//       // _child[2]
//       _child[2].setNeighbor(_child[3], 1);
//       _child[2].setNeighbor(_child[0], 4);
//       _child[2].setNeighbor(_child[1], 8);
//       if (_neighbors[2] != nullptr && !_neighbors[2]->isLeaf()) {
//         _child[2].setNeighbor(_neighbors[2]->_child[0], 2);
//         _child[2].setNeighbor(_neighbors[2]->_child[1], 5);
//       }
//       if (_neighbors[3] != nullptr && !_neighbors[3]->isLeaf()) {
//         _child[2].setNeighbor(_neighbors[3]->_child[3], 3);
//         _child[2].setNeighbor(_neighbors[3]->_child[1], 7);
//       }
//       if (_neighbors[6] != nullptr && !_neighbors[6]->isLeaf())
//         _child[2].setNeighbor(_neighbors[6]->_child[1], 6);

//       // _child[3]
//       _child[3].setNeighbor(_child[2], 3);
//       _child[3].setNeighbor(_child[1], 4);
//       _child[3].setNeighbor(_child[0], 7);
//       if (_neighbors[1] != nullptr && !_neighbors[1]->isLeaf()) {
//         _child[3].setNeighbor(_neighbors[1]->_child[2], 1);
//         _child[3].setNeighbor(_neighbors[1]->_child[0], 8);
//       }
//       if (_neighbors[2] != nullptr && !_neighbors[2]->isLeaf()) {
//         _child[3].setNeighbor(_neighbors[2]->_child[1], 2);
//         _child[3].setNeighbor(_neighbors[2]->_child[0], 6);
//       }
//       if (_neighbors[5] != nullptr && !_neighbors[5]->isLeaf())
//         _child[3].setNeighbor(_neighbors[5]->_child[0], 5);
//     }
//   }
//   void Setup() {
//     for (int i = 1; i < LatSet::q; ++i) {
//       const Vector<T, 2> nbrcenter =
//           this->_center + latset::c<LatSet>(i) * (T(2) * this->_radius);
//       if (_parent->_child[i]->isInside(nbrcenter)) {
//         setNeighbor(_parent->_child[i], i);
//       } else {
//         for (int j = 1; j < LatSet::q; ++j) {
//           if (_parent->_neighbors[j] != nullptr &&
//               _parent->_neighbors[j]->isInside(nbrcenter))
//             setNeighbor(_parent->_neighbors[j]->findChild(nbrcenter), i);
//         }
//       }
//     }
//   }
//   void Refine() {
//     _isLeaf = false;
//     _child = new QuadTree<T, LatSet>*[4];
//     Vector<T, 2> tmpCenter = this->_center;
//     T halfradius = this->_radius * T(0.5);
//     tmpCenter[0] = this->_center[0] - halfradius;
//     tmpCenter[1] = this->_center[1] - halfradius;
//     _child[0] =
//         new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag, true,
//         this);
//     tmpCenter[0] = this->_center[0] + halfradius;
//     tmpCenter[1] = this->_center[1] - halfradius;
//     _child[1] =
//         new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag, true,
//         this);
//     tmpCenter[0] = this->_center[0] - halfradius;
//     tmpCenter[1] = this->_center[1] + halfradius;
//     _child[2] =
//         new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag, true,
//         this);
//     tmpCenter[0] = this->_center[0] + halfradius;
//     tmpCenter[1] = this->_center[1] + halfradius;
//     _child[3] =
//         new QuadTree<T, LatSet>(_id, tmpCenter, _level + 1, _flag, true,
//         this);
//   }
//   const Vector<T, 2>& getCenter() const { return this->_center; }
//   BasicTree<T, 2>* getAbstractTree() { return this; }
//   short getLevel() const { return _level; }
// };

// }  // namespace quadtree
