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

#include "basic_geometry.h"
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
  T eps = std::numeric_limits<T>::epsilon();
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


// ---------------------basicblock----------------------

template <typename T, unsigned int D>
BasicBlock<T, D> BasicBlock<T, D>::getExtBlock(int ext) const {
  return getExtBlock(Vector<int, D>{ext});
}

template <typename T, unsigned int D>
BasicBlock<T, D> BasicBlock<T, D>::getExtBlock(const Vector<int, D>& extension) const {
  const Vector<T, D> extension_t = extension * VoxelSize;
  Vector<int, D> extension_idx;
  for (int i = 0; i < D; i++) {
    extension_idx[i] = std::ceil(extension[i] / pow(2, _level));
  }
  const Vector<int, D> extmesh = Mesh + 2 * extension;
  const AABB<T, D> extaabb = AABB<T, D>::getExtended(extension_t);
  const AABB<int, D> extidxblock = IndexBlock.getExtended(extension_idx);
  return BasicBlock<T, D>(_level, VoxelSize, BlockId, extaabb, extidxblock, extmesh);
}

template <typename T, unsigned int D>
BasicBlock<T, D> BasicBlock<T, D>::getRefinedBlock(std::uint8_t deltalevel) const {
  int ratio = pow(2, int(deltalevel));
  T newvoxsize = VoxelSize / T(ratio);
  const Vector<int, D> refmesh = Mesh * ratio;
  return BasicBlock<T, D>(deltalevel + _level, newvoxsize, BlockId, *this, IndexBlock, refmesh);
}


template <typename T, unsigned int D>
void BasicBlock<T, D>::refine(std::uint8_t deltalevel) {
  int ratio = pow(2, int(deltalevel));
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
  int ratio = pow(2, int(_deltalevel));
  T newvoxsize = VoxelSize * T(ratio);
  const Vector<int, D> coarsemesh = Mesh / ratio;
  return BasicBlock<T, D>(_level - _deltalevel, newvoxsize, BlockId, *this, IndexBlock, coarsemesh);
}

template <typename T, unsigned int D>
void BasicBlock<T, D>::coarsen(std::uint8_t deltalevel) {
  int ratio = pow(2, int(deltalevel));
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
void BasicBlock<T, D>::forEach(const AABB<T, D>& AABBs, Func func) {
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
template <typename flagtype, typename Func>
void BasicBlock<T, D>::forEach(const AABB<T, D>& AABBs, const GenericArray<flagtype>& flag,
                               std::uint8_t fromflag, Func func) {
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
template <typename flagtype, typename Func>
void BasicBlock<T, D>::forEach(const GenericArray<flagtype>& flag, std::uint8_t fromflag, Func func) {
  if constexpr (D == 2) {
    for (int j = 0; j < Mesh[1]; ++j) {
      for (int i = 0; i < Mesh[0]; ++i) {
        std::size_t id = getIndex(Vector<int, 2>{i, j});
        if (util::isFlag(flag[id], fromflag)) func(id);
      }
    }
  } else if constexpr (D == 3) {
    for (int k = 0; k < Mesh[2]; ++k) {
      for (int j = 0; j < Mesh[1]; ++j) {
        for (int i = 0; i < Mesh[0]; ++i) {
          std::size_t id = getIndex(Vector<int, 3>{i, j, k});
          if (util::isFlag(flag[id], fromflag)) func(id);
        }
      }
    }
  }
}

template <typename T, unsigned int D>
void BasicBlock<T, D>::getCellIdx(const AABB<T, D>& base, const AABB<T, D>& AABBs,
                                  std::vector<std::size_t>& cellIdx) const {
  // get intersection of AABBs and IndexBlock
  const AABB<T, D> intsec = getIntersection(base, AABBs);
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
        const Vector<T, 2> pt = MinCenter + (VoxelSize * idx);
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
          const Vector<T, 3> pt = MinCenter + (VoxelSize * idx);
          cellIdx.push_back(getIndex(idx));
        }
      }
    }
  }
}
