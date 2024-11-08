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

// basic_geometry.hh

#include "geometry/basic_geometry.h"

template <typename T, unsigned int D>
template <typename S>
bool AABB<T, D>::IsInside(const Vector<S, D>& pt) const {
  T eps = std::numeric_limits<decltype(T{} * S{})>::epsilon();
  if constexpr (D == 2) {
    return (pt[0] >= _min[0] - eps && pt[0] <= _max[0] + eps && pt[1] >= _min[1] - eps &&
            pt[1] <= _max[1] + eps);
  } else {
    return (pt[0] >= _min[0] - eps && pt[0] <= _max[0] + eps && pt[1] >= _min[1] - eps &&
            pt[1] <= _max[1] + eps && pt[2] >= _min[2] - eps && pt[2] <= _max[2] + eps);
  }
}

template <typename T, unsigned int D>
bool AABB<T, D>::isInside(const Vector<T, D>& pt) const {
  constexpr T eps = std::numeric_limits<T>::epsilon();
  if constexpr (D == 2) {
    return (pt[0] >= _min[0] - eps && pt[0] <= _max[0] + eps && pt[1] >= _min[1] - eps &&
            pt[1] <= _max[1] + eps);
  } else {
    return (pt[0] >= _min[0] - eps && pt[0] <= _max[0] + eps && pt[1] >= _min[1] - eps &&
            pt[1] <= _max[1] + eps && pt[2] >= _min[2] - eps && pt[2] <= _max[2] + eps);
  }
}

template <typename T, unsigned int D>
AABB<T, D> AABB<T, D>::getExtended(T extension) const {
  return getExtended(Vector<T, D>{extension});
}

template <typename T, unsigned int D>
AABB<T, D> AABB<T, D>::getExtended(const Vector<T, D>& extension) const {
  Vector<T, D> min = _min - extension;
  Vector<T, D> max = _max + extension;
  return AABB<T, D>(min, max);
}

template <typename T, unsigned int D>
AABB<T, D> AABB<T, D>::getResized(T factor) const {
  Vector<T, D> extension = _extension * factor;
  Vector<T, D> min = _center - extension / T(2);
  Vector<T, D> max = _center + extension / T(2);
  return AABB<T, D>(min, max);
}

// TODO: by shifting intsec in advance, getLocIdx of each for loop can be simplified
// this may be removed in the future, use BasicBlock<T, D>::getCellIdx instead
// beacuse this can't handle blocks with different refine levels
template <typename T, unsigned int D>
void AABB<T, D>::getOverlappedCellIdx(const AABB<int, D>& aabb, std::vector<int>& cellIdx) const {
  if (isOverlapped(*this, aabb)) {
    if constexpr (D == 2) {
      const AABB<int, D> intsec = getIntersection(*this, aabb);
      for (int j = intsec.getMin()[1]; j <= intsec.getMax()[1]; j++) {
        for (int i = intsec.getMin()[0]; i <= intsec.getMax()[0]; i++) {
          cellIdx.push_back(getLocIdx(Vector<int, 2>{i, j}));
        }
      }
    } else if constexpr (D == 3) {
      const AABB<int, D> intsec = getIntersection(*this, aabb);
      for (int k = intsec.getMin()[2]; k <= intsec.getMax()[2]; k++) {
        for (int j = intsec.getMin()[1]; j <= intsec.getMax()[1]; j++) {
          for (int i = intsec.getMin()[0]; i <= intsec.getMax()[0]; i++) {
            cellIdx.push_back(getLocIdx(Vector<int, 3>{i, j, k}));
          }
        }
      }
    }
  }
}

template <typename T, unsigned int D>
std::size_t AABB<T, D>::getLocIdx(const Vector<int, D>& pt) const {
  // 2d projection: {1, extension[0]+1}
  // 3d projection: {1, extension[0]+1, (extension[0]+1)*(extension[1]+1)}
  if constexpr (D == 2) {
    return pt[0] - _min[0] + (_extension[0] + 1) * (pt[1] - _min[1]);
  } else if constexpr (D == 3) {
    return pt[0] - _min[0] + (_extension[0] + 1) * (pt[1] - _min[1]) +
           (_extension[0] + 1) * (_extension[1] + 1) * (pt[2] - _min[2]);
  }
}

template <typename T, unsigned int D>
void AABB<T, D>::getOverlappedCellIdx(const AABB<int, D>& aabb, std::vector<int>& cellIdx,
                                      const Vector<int, D>& shift) const {
  if (isOverlapped(*this, aabb)) {
    if constexpr (D == 2) {
      const AABB<int, D> intsec = getIntersection(*this, aabb);
      for (int j = intsec.getMin()[1]; j <= intsec.getMax()[1]; j++) {
        for (int i = intsec.getMin()[0]; i <= intsec.getMax()[0]; i++) {
          cellIdx.push_back(getLocIdx(Vector<int, 2>{i, j}, shift));
        }
      }
    } else if constexpr (D == 3) {
      const AABB<int, D> intsec = getIntersection(*this, aabb);
      for (int k = intsec.getMin()[2]; k <= intsec.getMax()[2]; k++) {
        for (int j = intsec.getMin()[1]; j <= intsec.getMax()[1]; j++) {
          for (int i = intsec.getMin()[0]; i <= intsec.getMax()[0]; i++) {
            cellIdx.push_back(getLocIdx(Vector<int, 3>{i, j, k}, shift));
          }
        }
      }
    }
  }
}

template <typename T, unsigned int D>
std::size_t AABB<T, D>::getLocIdx(const Vector<int, D>& pt, const Vector<int, D>& shift) const {
  // 2d projection: {1, extension[0]+1}
  // 3d projection: {1, extension[0]+1, (extension[0]+1)*(extension[1]+1)}
  Vector<int, D> base = _min + shift;
  Vector<int, D> ext = _extension - (shift * 2);
  if constexpr (D == 2) {
    return pt[0] - base[0] + (ext[0] + 1) * (pt[1] - base[1]);
  } else if constexpr (D == 3) {
    return pt[0] - base[0] + (ext[0] + 1) * (pt[1] - base[1]) +
           (ext[0] + 1) * (ext[1] + 1) * (pt[2] - base[2]);
  }
}

template <typename T, unsigned int D>
void AABB<T, D>::divide(int Nx, int Ny, std::vector<AABB<int, 2>>& subAABBs) const {
  static_assert(D == 2, "Only 2D AABB can be divided into Nx * Ny parts.");
  if (Nx == 0 || Ny == 0){
    return;
  }
  int Nx_child = 0;
  int Ny_child = 0;

  T minx_child = _min[0];
  T miny_child = _min[1];

  // Calculate the size of each part
  int partX = (_extension[0] + 1) / Nx;
  int partY = (_extension[1] + 1) / Ny;
  int partX_ = partX + 1;
  int partY_ = partY + 1;

  // Calculate the remainder
  int remainderX = (_extension[0] + 1) % Nx;
  int remainderY = (_extension[1] + 1) % Ny;
  int rX = Nx - remainderX;
  int rY = Ny - remainderY;

  for (int iY = 0; iY < Ny; iY++) {
    if (iY < rY) {
      Ny_child = partY;
    } else {
      Ny_child = partY_;
    }
    for (int iX = 0; iX < Nx; iX++) {
      if (iX < rX) {
        Nx_child = partX;
      } else {
        Nx_child = partX_;
      }
      subAABBs.emplace_back(Vector<int, 2>{minx_child, miny_child},
                            Vector<int, 2>{minx_child + Nx_child - 1, miny_child + Ny_child - 1});
      minx_child += Nx_child;
    }
    minx_child = _min[0];
    miny_child += Ny_child;
  }
}

template <typename T, unsigned int D>
void AABB<T, D>::divide(int Nx, int Ny, int Nz, std::vector<AABB<int, 3>>& subAABBs) const {
  static_assert(D == 3, "Only 3D AABB can be divided into Nx * Ny * Nz parts.");
  if (Nx == 0 || Ny == 0 || Nz == 0){
    return;
  }
  int Nx_child = 0;
  int Ny_child = 0;
  int Nz_child = 0;
  T minx_child = _min[0];
  T miny_child = _min[1];
  T minz_child = _min[2];
  // Calculate the size of each part
  int partX = (_extension[0] + 1) / Nx;
  int partY = (_extension[1] + 1) / Ny;
  int partZ = (_extension[2] + 1) / Nz;
  int partX_ = partX + 1;
  int partY_ = partY + 1;
  int partZ_ = partZ + 1;
  // Calculate the remainder
  int remainderX = (_extension[0] + 1) % Nx;
  int remainderY = (_extension[1] + 1) % Ny;
  int remainderZ = (_extension[2] + 1) % Nz;
  int rX = Nx - remainderX;
  int rY = Ny - remainderY;
  int rZ = Nz - remainderZ;

  for (int iZ = 0; iZ < Nz; iZ++) {
    if (iZ < rZ) {
      Nz_child = partZ;
    } else {
      Nz_child = partZ_;
    }
    for (int iY = 0; iY < Ny; iY++) {
      if (iY < rY) {
        Ny_child = partY;
      } else {
        Ny_child = partY_;
      }
      for (int iX = 0; iX < Nx; iX++) {
        if (iX < rX) {
          Nx_child = partX;
        } else {
          Nx_child = partX_;
        }
        subAABBs.emplace_back(Vector<int, 3>{minx_child, miny_child, minz_child},
                              Vector<int, 3>{minx_child + Nx_child - 1, miny_child + Ny_child - 1,
                                             minz_child + Nz_child - 1});
        minx_child += Nx_child;
      }
      minx_child = _min[0];
      miny_child += Ny_child;
    }
    miny_child = _min[1];
    minz_child += Nz_child;
  }
}

// ---------------------basicblock----------------------

template <typename T, unsigned int D>
BasicBlock<T, D> BasicBlock<T, D>::getExtBlock(int ext) const {
  return getExtBlock(Vector<int, D>{ext});
}

template <typename T, unsigned int D>
BasicBlock<T, D> BasicBlock<T, D>::getExtBlock(const Vector<int, D>& extension) const {
  const Vector<T, D> extension_t = extension * VoxelSize;
  Vector<int, D> extension_idx;
  for (unsigned int i = 0; i < D; i++) {
    extension_idx[i] = std::ceil(extension[i] / std::pow(2, _level));
  }
  const Vector<int, D> extmesh = Mesh + 2 * extension;
  const AABB<T, D> extaabb = AABB<T, D>::getExtended(extension_t);
  const AABB<int, D> extidxblock = IndexBlock.getExtended(extension_idx);
  return BasicBlock<T, D>(_level, VoxelSize, BlockId, extaabb, extidxblock, extmesh);
}

template <typename T, unsigned int D>
BasicBlock<T, D> BasicBlock<T, D>::getRefinedBlock(std::uint8_t deltalevel) const {
  int ratio = std::pow(2, int(deltalevel));
  T newvoxsize = VoxelSize / T(ratio);
  const Vector<int, D> refmesh = Mesh * ratio;
  return BasicBlock<T, D>(deltalevel + _level, newvoxsize, BlockId, *this, IndexBlock, refmesh);
}


template <typename T, unsigned int D>
void BasicBlock<T, D>::refine(std::uint8_t deltalevel) {
  int ratio = std::pow(2, int(deltalevel));
  Mesh = Mesh * ratio;
  VoxelSize /= T(ratio);
  _level += deltalevel;
  MinCenter = AABB<T, D>::_min + Vector<T, D>{T(0.5) * VoxelSize};
  if constexpr (D == 2) {
    N = Mesh[0] * Mesh[1];
    Projection = Vector<int, 2>{1, Mesh[0]};
  } else if constexpr (D == 3) {
    N = Mesh[0] * Mesh[1] * Mesh[2];
    Projection = Vector<int, 3>{1, Mesh[0], Mesh[0] * Mesh[1]};
  }
}

template <typename T, unsigned int D>
BasicBlock<T, D> BasicBlock<T, D>::getCoasenedBlock(std::uint8_t _deltalevel) const {
  // check if the block can be coarsened
  if (_level == 0) {
    std::cerr << "The block is already the coarsest block" << std::endl;
  }
  int ratio = std::pow(2, int(_deltalevel));
  T newvoxsize = VoxelSize * T(ratio);
  const Vector<int, D> coarsemesh = Mesh / ratio;
  return BasicBlock<T, D>(_level - _deltalevel, newvoxsize, BlockId, *this, IndexBlock, coarsemesh);
}

template <typename T, unsigned int D>
void BasicBlock<T, D>::coarsen(std::uint8_t deltalevel) {
  int ratio = std::pow(2, int(deltalevel));
  Mesh = Mesh / ratio;
  VoxelSize *= T(ratio);
  _level -= deltalevel;
  MinCenter = AABB<T, D>::_min + Vector<T, D>{T(0.5) * VoxelSize};
  if constexpr (D == 2) {
    N = Mesh[0] * Mesh[1];
    Projection = Vector<int, 2>{1, Mesh[0]};
  } else if constexpr (D == 3) {
    N = Mesh[0] * Mesh[1] * Mesh[2];
    Projection = Vector<int, 3>{1, Mesh[0], Mesh[0] * Mesh[1]};
  }
}

template <typename T, unsigned int D>
std::size_t BasicBlock<T, D>::getIndex_t(const Vector<T, D>& loc) const {
  const Vector<T, D> ext = loc - MinCenter;
  const int x = static_cast<int>(std::round(ext[0] / VoxelSize));
  const int y = static_cast<int>(std::round(ext[1] / VoxelSize));
  return getIndex(Vector<int, D>{x, y});
}

template <typename T, unsigned int D>
std::size_t BasicBlock<T, D>::getIndex(const Vector<int, D>& locidx) const {
  if constexpr (D == 2) {
    return locidx[0] + locidx[1] * Mesh[0];
  } else if constexpr (D == 3) {
    return locidx[0] + locidx[1] * Projection[1] + locidx[2] * Projection[2];
  }
}

template <typename T, unsigned int D>
void BasicBlock<T, D>::getLoc(std::size_t id, Vector<int, D>& locidx) const {
  if constexpr (D == 2) {
    locidx[1] = id / Mesh[0];
    locidx[0] = id - locidx[1] * Mesh[0];
  } else if constexpr (D == 3) {
    locidx[2] = id / Projection[2];
    int temp = id - locidx[2] * Projection[2];
    locidx[1] = temp / Projection[1];
    locidx[0] = temp - locidx[1] * Projection[1];
  }
}

template <typename T, unsigned int D>
Vector<int, D> BasicBlock<T, D>::getLoc(std::size_t id) const {
  Vector<int, D> locidx;
  getLoc(id, locidx);
  return locidx;
}

template <typename T, unsigned int D>
void BasicBlock<T, D>::getLoc_t(std::size_t id, Vector<T, D>& loc) const {
  Vector<int, D> locidx;
  getLoc(id, locidx);
  loc = MinCenter + (VoxelSize * locidx);
}

template <typename T, unsigned int D>
Vector<T, D> BasicBlock<T, D>::getLoc_t(std::size_t id) const {
  Vector<T, D> loc;
  getLoc_t(id, loc);
  return loc;
}

template <typename T, unsigned int D>
Vector<T, D> BasicBlock<T, D>::getLoc_t(Vector<int, D>& locidx) const {
  return MinCenter + (VoxelSize * locidx);
}
//TODO: delete -1 in idxmax calculation
template <typename T, unsigned int D>
void BasicBlock<T, D>::getLocIdxRange(const AABB<T, D>& AABBs, Vector<int, D>& idx_min,
                                      Vector<int, D>& idx_max) const {
  const Vector<T, D>& min = AABBs.getMin();
  const Vector<T, D>& max = AABBs.getMax();
  if constexpr (std::is_same<T, int>::value) {
    for (unsigned int i = 0; i < D; ++i) {
      int idxmin = min[i] - IndexBlock.getMin()[i];
      idx_min[i] = idxmin < 0 ? 0 : idxmin;
    }
    for (unsigned int i = 0; i < D; ++i) {
      int idxmax = max[i] - IndexBlock.getMin()[i];
      idx_max[i] = idxmax > (Mesh[i] - 1) ? (Mesh[i] - 1) : idxmax;
    }
  } else {
    // get index of min and max
    for (unsigned int i = 0; i < D; ++i) {
      int idxmin = int(std::round((min[i] - AABB<T, D>::_min[i]) / VoxelSize));
      idx_min[i] = idxmin < 0 ? 0 : idxmin;
    }
    for (unsigned int i = 0; i < D; ++i) {
      int idxmax = int(std::round((max[i] - AABB<T, D>::_min[i]) / VoxelSize)) - 1;
      idx_max[i] = idxmax > (Mesh[i] - 1) ? (Mesh[i] - 1) : idxmax;
    }
  }
}

template <typename T, unsigned int D>
template <typename Func>
void BasicBlock<T, D>::forEach(const Func& func) {
  for (std::size_t id = 0; id < N; ++id) {
    func(id);
  }
}

template <typename T, unsigned int D>
template <typename Func>
void BasicBlock<T, D>::forEach(const AABB<T, D>& AABBs, const Func& func) {
  Vector<int, D> idx_min;
  Vector<int, D> idx_max;
  getLocIdxRange(AABBs, idx_min, idx_max);
  if constexpr (D == 2) {
    for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
      for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
        const Vector<T, 2> pt = MinCenter + (VoxelSize * Vector<T, 2>{T(i), T(j)});
        if (AABBs.isInside(pt)) func(getIndex(Vector<int, 2>{i, j}));
      }
    }
  } else if constexpr (D == 3) {
    for (int k = idx_min[2]; k <= idx_max[2]; ++k) {
      for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
        for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
          const Vector<T, 3> pt = MinCenter + (VoxelSize * Vector<T, 3>{T(i), T(j), T(k)});
          if (AABBs.isInside(pt)) func(getIndex(Vector<int, 3>{i, j, k}));
        }
      }
    }
  }
}

template <typename T, unsigned int D>
template <typename ArrayType, typename Func>
void BasicBlock<T, D>::forEach(const AABB<T, D>& AABBs, const ArrayType& flag,
                               std::uint8_t fromflag, const Func& func) {
  Vector<int, D> idx_min;
  Vector<int, D> idx_max;
  getLocIdxRange(AABBs, idx_min, idx_max);
  if constexpr (D == 2) {
    for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
      for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
        const Vector<T, 2> pt = MinCenter + (VoxelSize * Vector<T, 2>{T(i), T(j)});
        std::size_t id = getIndex(Vector<int, 2>{i, j});
        if (AABBs.isInside(pt) && util::isFlag(flag[id], fromflag)) func(id);
      }
    }
  } else if constexpr (D == 3) {
    for (int k = idx_min[2]; k <= idx_max[2]; ++k) {
      for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
        for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
          const Vector<T, 3> pt = MinCenter + (VoxelSize * Vector<T, 3>{T(i), T(j), T(k)});
          std::size_t id = getIndex(Vector<int, 3>{i, j, k});
          if (AABBs.isInside(pt) && util::isFlag(flag[id], fromflag)) func(id);
        }
      }
    }
  }
}

template <typename T, unsigned int D>
template <typename ArrayType, typename Func>
void BasicBlock<T, D>::forEach(const ArrayType& flag, std::uint8_t fromflag, const Func& func) {
  for (std::size_t id = 0; id < N; ++id) {
    if (util::isFlag(flag[id], fromflag)) func(id);
  }
  // if constexpr (D == 2) {
  //   for (int j = 0; j < Mesh[1]; ++j) {
  //     for (int i = 0; i < Mesh[0]; ++i) {
  //       std::size_t id = getIndex(Vector<int, 2>{i, j});
  //       if (util::isFlag(flag[id], fromflag)) func(id);
  //     }
  //   }
  // } else if constexpr (D == 3) {
  //   for (int k = 0; k < Mesh[2]; ++k) {
  //     for (int j = 0; j < Mesh[1]; ++j) {
  //       for (int i = 0; i < Mesh[0]; ++i) {
  //         std::size_t id = getIndex(Vector<int, 3>{i, j, k});
  //         if (util::isFlag(flag[id], fromflag)) func(id);
  //       }
  //     }
  //   }
  // }
}

template <typename T, unsigned int D>
void BasicBlock<T, D>::getCellIdx(const AABB<T, D>& AABB0, const AABB<T, D>& AABB1,
                                  std::vector<std::size_t>& cellIdx) const {
  // get intersection
  const AABB<T, D> intsec = getIntersection(AABB0, AABB1);
  // get Mesh index range
  Vector<int, D> idx_min;
  Vector<int, D> idx_max;
  getLocIdxRange(intsec, idx_min, idx_max);
  cellIdx.clear();
  if constexpr (D == 2) {
    cellIdx.reserve((idx_max[0] - idx_min[0] + 1) * (idx_max[1] - idx_min[1] + 1));
    for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
      for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
        const Vector<int, 2> idx = Vector<int, 2>{i, j};
        cellIdx.push_back(getIndex(idx));
      }
    }
  } else if constexpr (D == 3) {
    cellIdx.reserve((idx_max[0] - idx_min[0] + 1) * (idx_max[1] - idx_min[1] + 1) *
                    (idx_max[2] - idx_min[2] + 1));
    for (int k = idx_min[2]; k <= idx_max[2]; ++k) {
      for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
        for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
          const Vector<int, 3> idx = Vector<int, 3>{i, j, k};
          cellIdx.push_back(getIndex(idx));
        }
      }
    }
  }
}


template <typename T, unsigned int D>
bool BasicBlock<T, D>::ExcludeCornerIdx(
  std::vector<std::size_t>& cellIdxbase, std::vector<std::size_t>& cellIdxnbr, 
  std::vector<std::size_t>& excellIdxbase, std::vector<std::size_t>& excellIdxnbr) const {
  
  bool hasCorner = false;

  excellIdxbase.clear();
  excellIdxnbr.clear();

  auto it = cellIdxbase.begin();
  auto it_nbr = cellIdxnbr.begin();

  while (it != cellIdxbase.end()) {
    if (whichCorner(*it) != -1) {
      // if is corner cell
      // add to excluded cell list
      excellIdxbase.push_back(*it);
      excellIdxnbr.push_back(*it_nbr);
      // remove from original list
      it = cellIdxbase.erase(it);
      it_nbr = cellIdxnbr.erase(it_nbr);
      hasCorner = true;
    } else {
      ++it;
      ++it_nbr;
    }
  }
  return hasCorner;
}

template <typename T, unsigned int D>
bool BasicBlock<T, D>::ExcludeEdgeIdx(
  std::vector<std::size_t>& cellIdxbase, std::vector<std::size_t>& cellIdxnbr,
  std::vector<std::vector<std::size_t>>& excellIdxbase, std::vector<std::vector<std::size_t>>& excellIdxnbr) const {
  static_assert(D == 3, "ExcludeEdgeIdx is only for 3D block");

  bool hasEdge = false;

  excellIdxbase.clear();
  excellIdxnbr.clear();
  // 12 edges for 3D block
  excellIdxbase.resize(12);
  excellIdxnbr.resize(12);

  auto it = cellIdxbase.begin();
  auto it_nbr = cellIdxnbr.begin();

  while (it != cellIdxbase.end()) {
    if (whichEdge(*it) != -1) {
      for (int i = 0; i < 12; ++i) {
        // if is edge cell, add to excluded cell list
        if (whichEdge(*it) == i) {
          excellIdxbase[i].push_back(*it);
          excellIdxnbr[i].push_back(*it_nbr);
          break;
        }
      }
      // remove from original list
      it = cellIdxbase.erase(it);
      it_nbr = cellIdxnbr.erase(it_nbr);
      hasEdge = true;
    } else {
      ++it;
      ++it_nbr;
    }
  }
  return hasEdge;
}

template <typename T, unsigned int D>
bool BasicBlock<T, D>::ExcludeInnerIdx(
  std::vector<std::size_t>& cellIdxbase, std::vector<std::size_t>& cellIdxnbr,
  std::vector<std::size_t>& excellIdxbase, std::vector<std::size_t>& excellIdxnbr) const {
  
  bool hasInner = false;

  excellIdxbase.clear();
  excellIdxnbr.clear();

  auto it = cellIdxbase.begin();
  auto it_nbr = cellIdxnbr.begin();

  if constexpr (D == 2) {
    while (it != cellIdxbase.end()) {
      if (whichEdge(*it) != -1) {
        // if is edge cell, continue
        ++it;
        ++it_nbr;
        continue;
      } 
      // else if (whichCorner(*it) != -1) {
      //   // if is corner cell, continue
      //   ++it;
      //   ++it_nbr;
      //   continue;
      // } 
      else {
        // if is inner cell
        // add to excluded cell list
        excellIdxbase.push_back(*it);
        excellIdxnbr.push_back(*it_nbr);
        // remove from original list
        it = cellIdxbase.erase(it);
        it_nbr = cellIdxnbr.erase(it_nbr);
        hasInner = true;
      }
    }
  } else if constexpr (D == 3) {
    while (it != cellIdxbase.end()) {
      if (whichFace(*it) != -1) {
        // if is face cell, continue
        ++it;
        ++it_nbr;
        continue;
      } 
      // else if (whichEdge(*it) != -1) {
      //   // if is edge cell, continue
      //   ++it;
      //   ++it_nbr;
      //   continue;
      // } else if (whichCorner(*it) != -1) {
      //   // if is corner cell, continue
      //   ++it;
      //   ++it_nbr;
      //   continue;
      // }
       else {
        // if is inner cell
        // add to excluded cell list
        excellIdxbase.push_back(*it);
        excellIdxnbr.push_back(*it_nbr);
        // remove from original list
        it = cellIdxbase.erase(it);
        it_nbr = cellIdxnbr.erase(it_nbr);
        hasInner = true;
      }
    }
  }
  return hasInner;
}

template <typename T, unsigned int D>
void BasicBlock<T, D>::getCornerIdx(std::vector<std::size_t>& cornerIdx) const {
  if constexpr (D == 2) {
    // (0,0) (Nx,0) (0,Ny) (Nx,Ny)
    cornerIdx = {0, Mesh[0] - 1, Mesh[0] * (Mesh[1] - 1), Mesh[0] * Mesh[1] - 1};
  } else if constexpr (D == 3) {
    // (0,0,0) (Nx,0,0) (0,Ny,0) (Nx,Ny,0) (0,0,Nz) (Nx,0,Nz) (0,Ny,Nz) (Nx,Ny,Nz)
    cornerIdx = {0, Mesh[0] - 1, Mesh[0] * (Mesh[1] - 1), Mesh[0] * Mesh[1] - 1,
                 Mesh[0] * Mesh[1] * (Mesh[2] - 1), Mesh[0] * Mesh[1] * (Mesh[2] - 1) + Mesh[0] - 1,
                 Mesh[0] * Mesh[1] * Mesh[2] - Mesh[0], Mesh[0] * Mesh[1] * Mesh[2] - 1};
  }
}

template <typename T, unsigned int D>
void BasicBlock<T, D>::getEdgeIdx(std::vector<std::size_t>& edgeIdx, const Vector<int, D>& dir) const {
  edgeIdx.clear();
  if constexpr (D == 2) {
    if (dir[0] == 0) {
      // y direction
      edgeIdx.reserve(Mesh[0]);
      if (dir[1] > 0) {
        for (int i = 0; i < Mesh[0]; i++) {
          edgeIdx.push_back(getIndex(Vector<int, 2>{i, Mesh[1] - 1}));
        }
      } else {
        for (int i = 0; i < Mesh[0]; i++) {
          edgeIdx.push_back(getIndex(Vector<int, 2>{i, 0}));
        }
      }
    } else if (dir[1] == 0) {
      // x direction
      edgeIdx.reserve(Mesh[1]);
      if (dir[0] > 0) {
        for (int j = 0; j < Mesh[1]; j++) {
          edgeIdx.push_back(getIndex(Vector<int, 2>{Mesh[0] - 1, j}));
        }
      } else {
        for (int j = 0; j < Mesh[1]; j++) {
          edgeIdx.push_back(getIndex(Vector<int, 2>{0, j}));
        }
      }
    }
  } else if constexpr (D == 3) {
    if (dir[0] == 0) {
      // yz plane
      edgeIdx.reserve(Mesh[0]);
      if (dir[1] > 0) {
        if (dir[2] > 0) {
          for (int i = 0; i < Mesh[0]; ++i) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{i, Mesh[1] - 1, Mesh[2] - 1}));
          }
        } else {
          for (int i = 0; i < Mesh[0]; ++i) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{i, Mesh[1] - 1, 0}));
          }
        }
      } else {
        if (dir[2] > 0) {
          for (int i = 0; i < Mesh[0]; ++i) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{i, 0, Mesh[2] - 1}));
          }
        } else {
          for (int i = 0; i < Mesh[0]; ++i) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{i, 0, 0}));
          }
        }
      }
    } else if (dir[1] == 0) {
      // xz plane
      edgeIdx.reserve(Mesh[1]);
      if (dir[0] > 0) {
        if (dir[2] > 0) {
          for (int j = 0; j < Mesh[1]; ++j) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{Mesh[0] - 1, j, Mesh[2] - 1}));
          }
        } else {
          for (int j = 0; j < Mesh[1]; ++j) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{Mesh[0] - 1, j, 0}));
          }
        }
      } else {
        if (dir[2] > 0) {
          for (int j = 0; j < Mesh[1]; ++j) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{0, j, Mesh[2] - 1}));
          }
        } else {
          for (int j = 0; j < Mesh[1]; ++j) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{0, j, 0}));
          }
        }
      }
    } else if (dir[2] == 0) {
      // xy plane
      edgeIdx.reserve(Mesh[2]);
      if (dir[0] > 0) {
        if (dir[1] > 0) {
          for (int k = 0; k < Mesh[2]; ++k) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{Mesh[0] - 1, Mesh[1] - 1, k}));
          }
        } else {
          for (int k = 0; k < Mesh[2]; ++k) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{Mesh[0] - 1, 0, k}));
          }
        }
      } else {
        if (dir[1] > 0) {
          for (int k = 0; k < Mesh[2]; ++k) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{0, Mesh[1] - 1, k}));
          }
        } else {
          for (int k = 0; k < Mesh[2]; ++k) {
            edgeIdx.push_back(getIndex(Vector<int, 3>{0, 0, k}));
          }
        }
      }
    }
  }
}


template <typename T, unsigned int D>
void BasicBlock<T, D>::getFaceIdx(std::vector<std::size_t>& faceIdx, const Vector<int, D>& dir) const {
  static_assert(D == 3, "Only 3D block has face index.");
  faceIdx.clear();
  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
    std::cerr << "[BasicBlock<T, D>::getFaceIdx]: Invalid direction" << std::endl;
    exit(1);
  }
  if (dir[0] == 0 && dir[1] == 0) {
    // z direction
    faceIdx.reserve(Mesh[0] * Mesh[1]);
    if (dir[2] > 0) {
      for (int j = 0; j < Mesh[1]; ++j) {
        for (int i = 0; i < Mesh[0]; ++i) {
          faceIdx.push_back(getIndex(Vector<int, 3>{i, j, Mesh[2] - 1}));
        }
      }
    } else {
      for (int j = 0; j < Mesh[1]; ++j) {
        for (int i = 0; i < Mesh[0]; ++i) {
          faceIdx.push_back(getIndex(Vector<int, 3>{i, j, 0}));
        }
      }
    }
  } else if (dir[0] == 0 && dir[2] == 0) {
    // y direction
    faceIdx.reserve(Mesh[0] * Mesh[2]);
    if (dir[1] > 0) {
      for (int k = 0; k < Mesh[2]; ++k) {
        for (int i = 0; i < Mesh[0]; ++i) {
          faceIdx.push_back(getIndex(Vector<int, 3>{i, Mesh[1] - 1, k}));
        }
      }
    } else {
      for (int k = 0; k < Mesh[2]; ++k) {
        for (int i = 0; i < Mesh[0]; ++i) {
          faceIdx.push_back(getIndex(Vector<int, 3>{i, 0, k}));
        }
      }
    }
  } else if (dir[1] == 0 && dir[2] == 0) {
    // x direction
    faceIdx.reserve(Mesh[1] * Mesh[2]);
    if (dir[0] > 0) {
      for (int k = 0; k < Mesh[2]; ++k) {
        for (int j = 0; j < Mesh[1]; ++j) {
          faceIdx.push_back(getIndex(Vector<int, 3>{Mesh[0] - 1, j, k}));
        }
      }
    } else {
      for (int k = 0; k < Mesh[2]; ++k) {
        for (int j = 0; j < Mesh[1]; ++j) {
          faceIdx.push_back(getIndex(Vector<int, 3>{0, j, k}));
        }
      }
    }
  }
}

// 2d
// +y
// 2-----3
// |     |
// |     |
// 0-----1 +x
//
// 3d
//   +z
//   4-----6
//  /|    /|
// 5-----7 |
// | 0---|-2  +y
// |/    |/
// 1-----3
// +x

template <typename T, unsigned int D>
int BasicBlock<T, D>::whichCorner(const Vector<int, D>& pt) const {
  if constexpr (D == 2) {
    if (pt[0] == 0 && pt[1] == 0) {
      return 0;
    } else if (pt[0] == Mesh[0] - 1 && pt[1] == 0) {
      return 1;
    } else if (pt[0] == 0 && pt[1] == Mesh[1] - 1) {
      return 2;
    } else if (pt[0] == Mesh[0] - 1 && pt[1] == Mesh[1] - 1) {
      return 3;
    } else {
      return -1;
    }
  } else if constexpr (D == 3) {
    if (pt[0] == 0 && pt[1] == 0 && pt[2] == 0) {
      return 0;
    } else if (pt[0] == Mesh[0] - 1 && pt[1] == 0 && pt[2] == 0) {
      return 1;
    } else if (pt[0] == 0 && pt[1] == Mesh[1] - 1 && pt[2] == 0) {
      return 2;
    } else if (pt[0] == Mesh[0] - 1 && pt[1] == Mesh[1] - 1 && pt[2] == 0) {
      return 3;
    } else if (pt[0] == 0 && pt[1] == 0 && pt[2] == Mesh[2] - 1) {
      return 4;
    } else if (pt[0] == Mesh[0] - 1 && pt[1] == 0 && pt[2] == Mesh[2] - 1) {
      return 5;
    } else if (pt[0] == 0 && pt[1] == Mesh[1] - 1 && pt[2] == Mesh[2] - 1) {
      return 6;
    } else if (pt[0] == Mesh[0] - 1 && pt[1] == Mesh[1] - 1 && pt[2] == Mesh[2] - 1) {
      return 7;
    } else {
      return -1;
    }
  }
}
template <typename T, unsigned int D>
int BasicBlock<T, D>::whichCorner(std::size_t idx) const {
  const Vector<int, D> pt = getLoc(idx);
  return whichCorner(pt);
}

// 2d
// +y  3
//  -------
// 0|     |1
//  |     |
//  -------  +x
//     2
//
// 3d
//   +z
//      ---6---
//   10/|0 11/|2
//    ---7--  |
//  1| ---4|3-   +y
//   |/8   |/9
//   ---5---
// +x  

template <typename T, unsigned int D>
int BasicBlock<T, D>::whichEdge(const Vector<int, D>& pt) const {
  if constexpr (D == 2) {
    if (pt[0] == 0) {
      return 0;
    } else if (pt[0] == Mesh[0] - 1) {
      return 1;
    } else if (pt[1] == 0) {
      return 2;
    } else if (pt[1] == Mesh[1] - 1) {
      return 3;
    } else {
      return -1;
    }
  } else if constexpr (D == 3) {
    if (pt[0] == 0 && pt[1] == 0) {
      return 0;
    } else if (pt[0] == Mesh[0] - 1 && pt[1] == 0) {
      return 1;
    } else if (pt[0] == 0 && pt[1] == Mesh[1] - 1) {
      return 2;
    } else if (pt[0] == Mesh[0] - 1 && pt[1] == Mesh[1] - 1) {
      return 3;
    } else if (pt[0] == 0 && pt[2] == 0) {
      return 4;
    } else if (pt[0] == Mesh[0] - 1 && pt[2] == 0) {
      return 5;
    } else if (pt[0] == 0 && pt[2] == Mesh[2] - 1) {
      return 6;
    } else if (pt[0] == Mesh[0] - 1 && pt[2] == Mesh[2] - 1) {
      return 7;
    } else if (pt[1] == 0 && pt[2] == 0) {
      return 8;
    } else if (pt[1] == Mesh[1] - 1 && pt[2] == 0) {
      return 9;
    } else if (pt[1] == 0 && pt[2] == Mesh[2] - 1) {
      return 10;
    } else if (pt[1] == Mesh[1] - 1 && pt[2] == Mesh[2] - 1) {
      return 11;
    } else {
      return -1;
    }
  }
}
template <typename T, unsigned int D>
int BasicBlock<T, D>::whichEdge(std::size_t idx) const {
  const Vector<int, D> pt = getLoc(idx);
  return whichEdge(pt);
}


// 3d
//   +z
//      ---5---
//     /| 0  /|
//    ------  |
//  2| ----|--3 +y
//   |/ 1  |/
//   ---4---
// +x  

template <typename T, unsigned int D>
int BasicBlock<T, D>::whichFace(const Vector<int, D>& pt) const {
  static_assert(D == 3, "Only 3D block has face index.");
  // x face
  if (pt[0] == 0) {
    return 0;
  } else if (pt[0] == Mesh[0] - 1) {
    return 1;
  } 
  // y face
    else if (pt[1] == 0) {
    return 2;
  } else if (pt[1] == Mesh[1] - 1) {
    return 3;
  }
  // z face 
  else if (pt[2] == 0) {
    return 4;
  } else if (pt[2] == Mesh[2] - 1) {
    return 5;
  } else {
    return -1;
  }
}
template <typename T, unsigned int D>
int BasicBlock<T, D>::whichFace(std::size_t idx) const {
  const Vector<int, D> pt = getLoc(idx);
  return whichFace(pt);
}