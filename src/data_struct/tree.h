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

// tree.h

#pragma once

#include <cstdint>

template <typename T, unsigned int D>
struct Tree {
  bool _isLeaf;
  // refinement _level, non-refinement _level = 0
  std::uint8_t _level;
  // index in block
  std::size_t _id;

  Tree<T, D>* _parent;
  Tree<T, D>** _children;

  Tree(std::size_t id, std::uint8_t level, Tree<T, D>* parent = nullptr, bool isLeaf = true)
      : _id(id), _level(level), _parent(parent), _isLeaf(isLeaf) {}
  ~Tree() {
    if (!_isLeaf) {
      for (unsigned int i = 0; i < D; ++i) {
        delete _children[i];
      }
      delete[] _children;
    }
  }
  // void refine() {
  //   if (_isLeaf) {
  //     _isLeaf = false;
  //     _children = new Tree<T, D>*[D];
  //     for (unsigned int i = 0; i < D; ++i) {
  //       _children[i] = new Tree<T, D>(D * _id + i, _level + 1, this);
  //     }
  //   }
  // }
};