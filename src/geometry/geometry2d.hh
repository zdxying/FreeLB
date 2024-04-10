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

// geometry class
#pragma once

#include "geometry/geometry2d.h"

///////////////////////////////////////
// Geometry2D
////////////////////////////////////////

// peripheral cells will be automatically added to the geometry
template <typename T>
Geometry2D<T>::Geometry2D(int Nx, int Ny, const AABB<T, 2> &AABBs, T voxelSize,
                          const Vector<T, 2> &min, std::uint8_t AABBflag,
                          std::uint8_t voidflag)
    : _Nx(Nx + 2),
      _Ny(Ny + 2),
      N((Nx + 2) * (Ny + 2)),
      _voxelSize(voxelSize),
      AABB<T, 2>(
          min - Vector<T, 2>{voxelSize, voxelSize},
          min + Vector<T, 2>{(Nx + 1) * voxelSize, (Ny + 1) * voxelSize}),
      GeometryFlag((Nx + 2) * (Ny + 2), voidflag),
      Projection{1, Nx + 2},
      _AABBflag(AABBflag),
      _voidflag(voidflag),
      _MinCenter(min - Vector<T, 2>{voxelSize * T(0.5), voxelSize * T(0.5)}) {
  _Voxels.reserve(N);
  print();
  // read from AABBs
  ReadAABBs(AABBs, AABBflag);
}

template <typename T>
void Geometry2D<T>::ReadAABBs(const AABB<T, 2> &AABBs, std::uint8_t AABBflag) {
  // loop over all sites in the AABB defined by (0,0,0) and (Nx,Ny,Nz)
  for (int j = 0; j < _Ny; ++j) {
    for (int i = 0; i < _Nx; ++i) {
      // loc in 2d array
      Vector<int, 2> locidx{i, j};
      // get the location of the voxel
      Vector<T, 2> loc = _voxelSize * locidx + _MinCenter;
      // add to _Voxels
      _Voxels.push_back(loc);
      // check if it is inside
      if (AABBs.isInside(loc))
        GeometryFlag.SetField(getIndex(locidx), AABBflag);
    }
  }
  std::cout << "[Geometry2D]: ReadAABBs Done!"
            << "\n"
            << "Voxels size: " << _Voxels.size() << std::endl;
}

template <typename T>
template <typename LatSet>
void Geometry2D<T>::SetupBoundary(std::uint8_t AABBflag,
                                  std::uint8_t boundaryflag) {
  FlagField TransFlag(N, std::uint8_t(0));
  for (int id = 0; id < N; ++id) {
    if (GeometryFlag.get(id) == AABBflag) {
      const Vector<T, 2> &voxel = _Voxels[id];
      int count = 0;
      for (int i = 1; i < LatSet::q; ++i) {
        Vector<T, 2> pt = voxel + LatSet::c[i] * _voxelSize;
        if (this->isInside(pt)) {
          int idx = LatSet::c[i] * Projection + id;
          if (GeometryFlag.get(idx) == AABBflag) count++;
        }
      }
      if (count < LatSet::q - 1) TransFlag.SetField(id, std::uint8_t(1));
    }
  }
  for (int id = 0; id < N; ++id) {
    if (static_cast<bool>(TransFlag.get(id)))
      GeometryFlag.SetField(id, boundaryflag);
  }
}

template <typename T>
void Geometry2D<T>::setFlag(const AABB<T, 2> &AABBs, std::uint8_t flag) {
  // // get cuboid defined by AABBs
  // Vector<T, 2> min = AABBs.getMin();
  // Vector<T, 2> max = AABBs.getMax();
  // // get index of min and max
  // Vector<int, 2> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
  //                        int((min[1] - this->getMin()[1]) / _voxelSize));
  // for (int i = 0; i < 2; ++i) {
  //   if (idx_min[i] < 0) idx_min[i] = 0;
  // }
  // Vector<int, 2> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
  //                        int((max[1] - this->getMin()[1]) / _voxelSize));
  // if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  // if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  // // set flag
  // for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
  //   for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
  //     int id = getIndex(Vector<int, 2>(i, j));
  //     // check if there is a voxel
  //     const Vector<T, 2> &voxel = _Voxels[id];
  //     // check if voxel is inside AABBs
  //     if (AABBs.isInside(voxel)) GeometryFlag.SetField(id, flag);
  //   }
  // }
  // use lambda function
  forEachVoxel(AABBs,
               [this, flag](int id) { GeometryFlag.SetField(id, flag); });
}

template <typename T>
void Geometry2D<T>::setFlag(const AABB<T, 2> &AABBs, std::uint8_t fromflag,
                            std::uint8_t flag) {
  // // get cuboid defined by AABBs
  // Vector<T, 2> min = AABBs.getMin();
  // Vector<T, 2> max = AABBs.getMax();
  // // get index of min and max
  // Vector<int, 2> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
  //                        int((min[1] - this->getMin()[1]) / _voxelSize));
  // for (int i = 0; i < 2; ++i) {
  //   if (idx_min[i] < 0) idx_min[i] = 0;
  // }
  // Vector<int, 2> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
  //                        int((max[1] - this->getMin()[1]) / _voxelSize));
  // if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  // if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  // // set flag
  // for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
  //   for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
  //     int id = getIndex(Vector<int, 2>(i, j));
  //     const Vector<T, 2> &voxel = _Voxels[id];
  //     // check if voxel is inside AABBs
  //     if (AABBs.isInside(voxel) && GeometryFlag.get(id) == fromflag)
  //       GeometryFlag.SetField(id, flag);
  //   }
  // }
  // use lambda function
  forEachVoxel(AABBs, fromflag,
               [this, flag](int id) { GeometryFlag.SetField(id, flag); });
}

// lambda function
template <typename T>
template <typename Func>
void Geometry2D<T>::forEachVoxel(const AABB<T, 2> &AABBs, Func func) {
  // get cuboid defined by AABBs
  Vector<T, 2> min = AABBs.getMin();
  Vector<T, 2> max = AABBs.getMax();
  // get index of min and max
  Vector<int, 2> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize));
  if (idx_min[0] < 0) idx_min[0] = 0;
  if (idx_min[1] < 0) idx_min[1] = 0;
  Vector<int, 2> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  // set flag
  for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
    for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
      int id = getIndex(Vector<int, 2>(i, j));
      const Vector<T, 2> &voxel = _Voxels[id];
      // check if voxel is inside AABBs
      if (AABBs.isInside(voxel)) {
        func(id);
      }
    }
  }
}
template <typename T>
template <typename Func>
void Geometry2D<T>::forEachVoxelint(const AABB<int, 2> &AABBs, Func func) {
  if constexpr (std::is_same<T, int>::value) {
    // get index of min and max
    Vector<int, 2> idx_min = AABBs.getMin();
    if (idx_min[0] < 0) idx_min[0] = 0;
    if (idx_min[1] < 0) idx_min[1] = 0;
    Vector<int, 2> idx_max = AABBs.getMax();
    if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
    if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
    // set flag
    for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
      for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
        const Vector<int, 2> pos{i, j};
        int id = getIndex(pos);
        // check if voxel is inside AABBs
        if (AABBs.isInside(pos)) {
          func(id);
        }
      }
    }
  }
}
template <typename T>
template <typename Func>
void Geometry2D<T>::forEachVoxel(const AABB<T, 2> &AABBs, std::uint8_t fromflag,
                                 Func func) {
  // get cuboid defined by AABBs
  Vector<T, 2> min = AABBs.getMin();
  Vector<T, 2> max = AABBs.getMax();
  // get index of min and max
  Vector<int, 2> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize));
  if (idx_min[0] < 0) idx_min[0] = 0;
  if (idx_min[1] < 0) idx_min[1] = 0;
  Vector<int, 2> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  // set flag
  for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
    for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
      int id = getIndex(Vector<int, 2>(i, j));
      const Vector<T, 2> &voxel = _Voxels[id];
      // check if voxel is inside AABBs
      if (AABBs.isInside(voxel) && GeometryFlag.get(id) == fromflag) {
        func(id);
      }
    }
  }
}
template <typename T>
template <typename Func>
void Geometry2D<T>::forEachVoxel(std::uint8_t fromflag, Func func) {
  for (int id = 0; id < N; ++id) {
    if (GeometryFlag.get(id) == fromflag) {
      func(id);
    }
  }
}

