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

// refined_geometry2d.h

#pragma once

#include "geometry/geometry2d.h"

// (local) refined geometry
template <typename T>
class RefinedGeometry2D {
 private:
  // mesh data:
  // id_: coarse grid cell global index
  std::vector<QuadTree<T>> _QuadTrees;
  // global index of refined grid cell(_QuadTrees[x])
  std::vector<int> _Index;
  // index of refined grid cells of all levels
  std::vector<std::vector<QuadTree<T>*>> _AllIdx;
  Geometry2D<T>& CoarseGeo;
  // max level
  std::uint8_t maxLevel;
  // refine level flag field
  FlagField RefineLevels;
  // refine level flag field for vtk writer(each cell)
  std::vector<std::uint8_t> RefineLevelsVtk;
  // flags for each cell
  std::vector<std::uint8_t> FlagsVtk;
  // abstract tree for coarse grid
  std::vector<BasicTree<T, 2>> _CTree;
  // abstract tree for all grid
  std::vector<BasicTree<T, 2>*> _Tree;
  // block geometry
  std::vector<std::vector<AABB<int, 2>>> _BlockAABBs;
  // block geometry for all grid
  std::vector<int> _BlockAABBsVtk;
  // blockindex
  ScalerField<int> _BlockIndex;
  // coarse tree in AABBs
  // std::vector<BasicTree<T, 2>*> _CTree_inAABBs;

 public:
  RefinedGeometry2D(Geometry2D<T>& coarseGeo, std::uint8_t maxlevel)
      : CoarseGeo(coarseGeo),
        maxLevel(maxlevel),
        RefineLevels(coarseGeo.getVoxelsNum(), std::uint8_t(0)),
        _BlockIndex(coarseGeo.getVoxelsNum(), 0) {
    _Index.resize(CoarseGeo.getVoxelsNum(), 0);
    for (int id = 0; id < CoarseGeo.getVoxelsNum(); ++id) {
      const Vector<T, 2>& vox = CoarseGeo.getVoxel(id);
      _CTree.emplace_back(vox, T(0.5), id, CoarseGeo.getGeoFlag(id));
    }
  }

  Geometry2D<T>& getCoarseGeo() { return CoarseGeo; }
  FlagField& getRefineLevels() { return RefineLevels; }
  const FlagField& getRefineLevels() const { return RefineLevels; }
  ScalerField<int>& getBlockIndex() { return _BlockIndex; }
  const ScalerField<int>& getBlockIndex() const { return _BlockIndex; }
  std::vector<std::uint8_t>& getRefineLevelsVtk() { return RefineLevelsVtk; }
  const std::vector<std::uint8_t>& getRefineLevelsVtk() const {
    return RefineLevelsVtk;
  }
  std::vector<std::uint8_t>& getFlagsVtk() { return FlagsVtk; }
  const std::vector<std::uint8_t>& getFlagsVtk() const { return FlagsVtk; }
  std::vector<BasicTree<T, 2>*>& getAbstractTree() { return _Tree; }
  const std::vector<BasicTree<T, 2>*>& getAbstractTree() const { return _Tree; }

  // setup from AABBs and flagfield
  template <template <typename> class ArrayType, typename FlagType>
  void Setup(const ArrayType<FlagType>& flagarr, const FlagType voidflag,
             const AABB<T, 2>& AABBs, int layer) {
    CoarseGeo.forEachVoxel(AABBs, [&](int id) {
      if (flagarr[id] != voidflag) {
        RefineLevels.SetField(id, maxLevel);
        const Vector<T, 2>& vox = CoarseGeo.getVoxel(id);
        _Index[id] = _QuadTrees.size();
        _QuadTrees.emplace_back(id, vox, 0, maxLevel,
                                CoarseGeo.getGeoFlag(id));
      }
    });
    if (layer > 1) {
      AddLayer(flagarr, voidflag, layer - 1, maxLevel);
    }
    //
    std::uint8_t level = maxLevel - 1;
    while (level > std::uint8_t(0)) {
      AddLayer(flagarr, voidflag, layer, level);
      --level;
    }
  }
  // setup from flagfield, e.g.: GenericArray<std::uint_8t>
  template <template <typename> class ArrayType, typename FlagType>
  void Setup(const ArrayType<FlagType>& flagarr, const FlagType flag,
             const FlagType voidflag, int layer) {
    CoarseGeo.forEachVoxel(flag, [&](int id) {
      if (flagarr[id] != voidflag) {
        RefineLevels.SetField(id, maxLevel);
        const Vector<T, 2>& vox = CoarseGeo.getVoxel(id);
        _Index[id] = _QuadTrees.size();
        _QuadTrees.emplace_back(id, vox, 0, maxLevel,
                                CoarseGeo.getGeoFlag(id));
      }
    });
    if (layer > 1) {
      AddLayer(flagarr, voidflag, layer - 1, maxLevel);
    }
    //
    std::uint8_t level = maxLevel - 1;
    while (level > std::uint8_t(0)) {
      AddLayer(flagarr, voidflag, layer, level);
      --level;
    }
  }

  // recursively add layer
  template <template <typename> class ArrayType, typename FlagType>
  void AddLayer(const ArrayType<FlagType>& flagarr, const FlagType voidflag,
                int layer, std::uint8_t level) {
    std::size_t PrevSize = _QuadTrees.size();
    std::vector<QuadTree<T>> NewQuadTrees;
    for (const QuadTree<T>& tree : _QuadTrees) {
      std::size_t id = tree.getId();
      // find neighbors with flag != voidflag && RefineLevel == 0
      for (int k = 1; k < 9; ++k) {
        int nbrid = CoarseGeo.getNbrId(id, k);
        if (nbrid >= 0 && nbrid < flagarr.size()) {
          if (flagarr[nbrid] != voidflag &&
              RefineLevels.get(nbrid) == std::uint8_t(0)) {
            RefineLevels.SetField(nbrid, level);
            const Vector<T, 2>& vox = CoarseGeo.getVoxel(nbrid);
            _Index[nbrid] = PrevSize + NewQuadTrees.size();
            NewQuadTrees.emplace_back(nbrid, vox, 0, level,
                                      CoarseGeo.getGeoFlag(nbrid));
          }
        }
      }
    }
    _QuadTrees.insert(_QuadTrees.end(), NewQuadTrees.begin(),
                      NewQuadTrees.end());
    --layer;
    if (layer > 0) {
      AddLayer(flagarr, voidflag, layer, level);
    }
  }

  // update geometry flag based on AABBs
  // Abstract Tree must be updated before calling this function
  void UpdateFlag(const AABB<T, 2>& AABBs, std::uint8_t flag) {
#pragma omp parallel for num_threads(Thread_Num)
    for (QuadTree<T>& tree : _QuadTrees) {
      tree.forEachLeaf([&](QuadTree<T>* leaf) {
        if (AABBs.isInside(leaf->getCenter())) {
          leaf->setFlag(flag);
        }
      });
    }
  }
  void UpdateAbstractTree() {
    _Tree.clear();
    // get simple tree from coarse grid
    for (int id = 0; id < CoarseGeo.getVoxelsNum(); ++id) {
      if (RefineLevels.get(id) == std::uint8_t(0)) {
        _Tree.push_back(&_CTree[id]);
      } else {
        _QuadTrees[_Index[id]].getLeafs(_Tree);
      }
    }
  }
  std::vector<BasicTree<T, 2>*>& getUpdatedAbstractTree() {
    UpdateAbstractTree();
    return _Tree;
  }
  // Abstract Tree must be updated before calling this function
  void UpdateRefineLevelsVtk() {
    RefineLevelsVtk.resize(_Tree.size(), 0);
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id = 0; id < _Tree.size(); ++id) {
      RefineLevelsVtk[id] = RefineLevels.get(_Tree[id]->getId());
    }
  }
  // Abstract Tree must be updated before calling this function
  void UpdateFlagsVtk() {
    FlagsVtk.resize(_Tree.size(), 0);
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id = 0; id < _Tree.size(); ++id) {
      FlagsVtk[id] = _Tree[id]->getFlag();
    }
  }
  // Abstract Tree must be updated before calling this function
  void UpdateBlockAABBVtk() {
    // _BlockAABBsVtk.resize(_CTree.size(), 0);
    int count = 0;
    for (const AABB<int, 2>& aabb : _BlockAABBs[0]) {
      count++;
      CoarseGeo.forEachVoxelint(
          aabb, [&](int id) { _BlockIndex.SetField(id, count); });
    }
  }

  void UpdateBlockAABB(std::uint8_t flag, T limfrac = T(0.7)) {
    // split domain into blocks for each refine level
    // level 1
    _BlockAABBs.clear();
    std::vector<AABB<int, 2>> BlockAABBs;
    // split using a rectangle structure
    // find the min and max of the tree
    T minx = CoarseGeo.getMax()[0];
    T miny = CoarseGeo.getMax()[1];
    T maxx = CoarseGeo.getMin()[0];
    T maxy = CoarseGeo.getMin()[1];
    for (const QuadTree<T>& tree : _QuadTrees) {
      const Vector<T, 2>& centre = tree.getCenter();
      if (centre[0] < minx) minx = centre[0];
      if (centre[0] > maxx) maxx = centre[0];
      if (centre[1] < miny) miny = centre[1];
      if (centre[1] > maxy) maxy = centre[1];
    }
    int xmin = CoarseGeo.getGridx(minx);
    int xmax = CoarseGeo.getGridx(maxx);
    int ymin = CoarseGeo.getGridy(miny);
    int ymax = CoarseGeo.getGridy(maxy);
    AABB<int, 2> LargeRec(Vector<int, 2>(xmin, ymin),
                          Vector<int, 2>(xmax, ymax));
    // get all coarse trees in aabb
    // for (const BasicTree<T, 2>& tree : _CTree) {
    //   if (LargeRec.isInside(tree.getCenter()))
    //   _CTree_inAABBs.push_back(&tree);
    // }
    // split the domain into blocks
    SplitAABB(flag, LargeRec, BlockAABBs, limfrac);
    _BlockAABBs.push_back(BlockAABBs);
  }

  // if limit fraction is not reached, split the AABBs
  void SplitAABB(std::uint8_t flag, const AABB<int, 2>& AABBs,
                 std::vector<AABB<int, 2>>& BlockAABBs, T limfrac) {
    std::size_t Total = 0;
    std::size_t Num = 0;
    CoarseGeo.forEachVoxelint(AABBs, [&](int id) {
      ++Total;
      if (static_cast<bool>(_CTree[id].getFlag() & flag)) ++Num;
    });
    if (Total < 16) return;
    // for (const BasicTree<T, 2>& tree : _CTree) {
    //   if (AABBs.isInside(tree.getCenter())) {
    //     ++Total;
    //     if (static_cast<bool>(tree.getFlag() & flag)) ++Num;
    //   }
    // }
    if (T(Num) / T(Total) < limfrac) {
      // split AABBs
      // get longest side
      if (AABBs.getExtension()[0] < AABBs.getExtension()[1]) {
        // split along y
        int midy = (AABBs.getMin()[1] + AABBs.getMax()[1]) / T(2);
        AABB<int, 2> subAABB1(AABBs.getMin(),
                              Vector<int, 2>(AABBs.getMax()[0], midy));
        SplitAABB(flag, subAABB1, BlockAABBs, limfrac);
        AABB<int, 2> subAABB2(Vector<int, 2>(AABBs.getMin()[0], midy + 1),
                              AABBs.getMax());
        SplitAABB(flag, subAABB2, BlockAABBs, limfrac);
      } else {
        // split along x
        int midx = (AABBs.getMin()[0] + AABBs.getMax()[0]) / T(2);
        AABB<int, 2> subAABB1(AABBs.getMin(),
                              Vector<int, 2>(midx, AABBs.getMax()[1]));
        SplitAABB(flag, subAABB1, BlockAABBs, limfrac);
        AABB<int, 2> subAABB2(Vector<int, 2>(midx + 1, AABBs.getMin()[1]),
                              AABBs.getMax());
        SplitAABB(flag, subAABB2, BlockAABBs, limfrac);
      }
    } else {
      BlockAABBs.push_back(AABBs);
    }
  }

  //   void SetupNbr() {
  //     SetupLevel0();
  //     int level = 0;
  //     while (level < maxLevel) {
  //       ++level;
  //       SetupNbrofLevel(level);
  //     }
  //   }
  //   // init fine grid of level = 1 setup neighbors
  //   void SetupLevel0() {
  // #pragma omp parallel for num_threads(Thread_Num) schedule(static)
  //     for (QuadTree<T, LatSet>& vox : _QuadTrees) {
  //       int id = vox.getCId();
  //       for (int i = 1; i < LatSet::q; ++i) {
  //         int nbrid = CoarseGeo.getNbrId(id, i);
  //         int pos_in_vector = _Index[nbrid];
  //         if (pos_in_vector != 0) vox.setNeighbor(&_QuadTrees[pos_in_vector],
  //         i);
  //       }
  //     }
  //   }
  //   // set neighbor of level x > 1
  //   void SetupNbrofLevel(short level) {
  //     const std::vector<QuadTree<T, LatSet>*>& idx = _AllIdx[level - 1];
  // #pragma omp parallel for num_threads(Thread_Num) schedule(static)
  //     for (QuadTree<T, LatSet>* vox : idx) {
  //       vox->SetupChildNbr();
  //     }
  //   }
  // void GetVoxelofLevel(std::vector<std::vector<QuadTree<T, LatSet>*>>& IDX) {
  //   IDX.clear();
  //   std::vector<QuadTree<T, LatSet>*> idx_level_c;
  //   GetVoxelofLevel1(idx_level_c);
  //   IDX.push_back(idx_level_c);
  //   int level = 1;
  //   while (level < maxLevel) {
  //     std::vector<QuadTree<T, LatSet>*> idx_level_f;
  //     GetVoxelofNextLevel(idx_level_c, idx_level_f);
  //     IDX.push_back(idx_level_f);
  //     ++level;
  //     idx_level_c = idx_level_f;
  //   }
  // }
  // void GetVoxelofLevel1(std::vector<QuadTree<T, LatSet>*>& idx) {
  //   idx.clear();
  //   int count = 0;
  //   for (QuadTree<T, LatSet>& vox : _QuadTrees) {
  //     vox.SetFId(count);
  //     idx.push_back(&vox);
  //     ++count;
  //   }
  // }
  // void GetVoxelofNextLevel(std::vector<QuadTree<T, LatSet>*>& Cidx,
  //                          std::vector<QuadTree<T, LatSet>*>& Fidx) {
  //   Fidx.clear();
  //   int count = 0;
  //   for (QuadTree<T, LatSet>* vox : Cidx) {
  //     if (!vox->isLeaf()) {
  //       QuadTree<T, LatSet>* voxchild0 = vox->getChild(0);
  //       voxchild0->SetFId(count);
  //       Fidx.push_back(voxchild0);
  //       ++count;
  //       QuadTree<T, LatSet>* voxchild1 = vox->getChild(1);
  //       voxchild1->SetFId(count);
  //       Fidx.push_back(voxchild1);
  //       ++count;
  //       QuadTree<T, LatSet>* voxchild2 = vox->getChild(2);
  //       voxchild2->SetFId(count);
  //       Fidx.push_back(voxchild2);
  //       ++count;
  //       QuadTree<T, LatSet>* voxchild3 = vox->getChild(3);
  //       voxchild3->SetFId(count);
  //       Fidx.push_back(voxchild3);
  //       ++count;
  //     }
  //   }
  // }
};


