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

// geometry3d.hh
#pragma once

#include "geometry/geometry3d.h"

//////////////////////////////////////////////////////////////
/////// Implementation of Geometry3D
//////////////////////////////////////////////////////////////

template <typename T>
Geometry3D<T>::Geometry3D(const StlReader<T>& reader, std::uint8_t stlflag)
    : _voxelSize(reader.getVoxelSize()),
      AABB<T, 3>(reader.getMesh().getMin() - reader.getVoxelSize(),
                 reader.getMesh().getMax() + reader.getVoxelSize()),
      _Nx(int(reader.getMesh().getMax_Min()[0] / reader.getStlSize()) + 2),
      _Ny(int(reader.getMesh().getMax_Min()[1] / reader.getStlSize()) + 2),
      _Nz(int(reader.getMesh().getMax_Min()[2] / reader.getStlSize()) + 2),
      GeometryFlag(
          (int(reader.getMesh().getMax_Min()[0] / reader.getStlSize()) + 2) *
              (int(reader.getMesh().getMax_Min()[1] / reader.getStlSize()) +
               2) *
              (int(reader.getMesh().getMax_Min()[2] / reader.getStlSize()) + 2),
          std::uint8_t(0)),
      Projection{
          1, int(reader.getMesh().getMax_Min()[0] / reader.getStlSize()) + 2,
          (int(reader.getMesh().getMax_Min()[0] / reader.getStlSize()) + 2) *
              (int(reader.getMesh().getMax_Min()[1] / reader.getStlSize()) +
               2)} {
  N = _Nx * _Ny * _Nz;
  _Voxels.reserve(N);
  print();
  // read from octree
  ReadOctree(reader.getTree(), stlflag);
}

template <typename T>
Geometry3D<T>::Geometry3D(int Nx, int Ny, int Nz, const AABB<T, 3>& AABBs,
                          T voxelSize, const Vector<T, 3>& min,
                          std::uint8_t AABBflag, std::uint8_t voidflag)
    :AABB<T, 3>(min - voxelSize,
                 min + Vector<T, 3>{(Nx + 1) * voxelSize, (Ny + 1) * voxelSize,
                                    (Nz + 1) * voxelSize}), 
      _Nx(Nx + 2),
      _Ny(Ny + 2),
      _Nz(Nz + 2),
      N((Nx + 2) * (Ny + 2) * (Nz + 2)),
      _voxelSize(voxelSize),
      _AABBflag(AABBflag),
      _voidflag(voidflag),
      GeometryFlag((Nx + 2) * (Ny + 2) * (Nz + 2), voidflag) {
  _Voxels.reserve(N);
  print();
  // read from AABBs
  ReadAABBs(AABBs, AABBflag);
}

template <typename T>
Geometry3D<T>::Geometry3D(int Nx, int Ny, int Nz, const AABB<T, 3>& AABBs,
                          const StlReader<T>& reader, std::uint8_t AABBflag,
                          std::uint8_t stlflag, const Vector<T, 3>& min)
    : _Nx(Nx + 2),
      _Ny(Ny + 2),
      _Nz(Nz + 2),
      N((Nx + 2) * (Ny + 2) * (Nz + 2)),
      _voxelSize(reader.getVoxelSize()),
      AABB<T, 3>(min - reader.getVoxelSize(),
                 min + Vector<T, 3>{(Nx + 1) * reader.getVoxelSize(),
                                    (Ny + 1) * reader.getVoxelSize(),
                                    (Nz + 1) * reader.getVoxelSize()}),
      GeometryFlag((Nx + 2) * (Ny + 2) * (Nz + 2), std::uint8_t(0)),
      Projection{1, Nx + 2, (Nx + 2) * (Ny + 2)} {
  _Voxels.reserve(N);
  print();
  ReadAABBs(AABBs, AABBflag);
  // set flag from stlreader
  setFlag(AABBs, reader, stlflag);
}

template <typename T>
void Geometry3D<T>::ReadOctree(Octree<T>* tree, std::uint8_t stlflag) {
  // loop over all sites in the AABB defined by (0,0,0) and (Nx,Ny,Nz)
  for (int k = 0; k < _Nz; ++k) {
    for (int j = 0; j < _Ny; ++j) {
      for (int i = 0; i < _Nx; ++i) {
        // loc in 3d array
        Vector<int, 3> locarray{i, j, k};
        // get the location of the voxel
        Vector<T, 3> loc =
            _voxelSize * locarray + this->_min + _voxelSize * T(0.5);
        // add to _Voxels
        _Voxels.push_back(loc);
        // get the node containing the voxel
        Octree<T>* node = tree->find(loc);
        if (node != nullptr) {
          // check if it is a [leaf] node and if it is [inside]
          if (node->isLeaf() && node->getInside())
            GeometryFlag.SetField(getIndex(locarray), stlflag);
        }
      }
    }
  }
  std::cout << "[Geometry3D]: ReadOctree Done!"
            << "\n"
            << "Voxels size: " << _Voxels.size() << std::endl;
}

template <typename T>
void Geometry3D<T>::ReadAABBs(const AABB<T, 3>& AABBs, std::uint8_t AABBflag) {
  // loop over all sites in the AABB defined by (0,0,0) and (Nx,Ny,Nz)
  for (int k = 0; k < _Nz; ++k) {
    for (int j = 0; j < _Ny; ++j) {
      for (int i = 0; i < _Nx; ++i) {
        // loc in 3d array
        Vector<int, 3> locarray{i, j, k};
        // get the location of the voxel
        Vector<T, 3> loc =
            _voxelSize * locarray + this->_min + _voxelSize * T(0.5);
        // add to _Voxels
        _Voxels.push_back(loc);
        // check if it is inside
        if (AABBs.isInside(loc))
          GeometryFlag.SetField(getIndex(locarray), AABBflag);
      }
    }
  }
  std::cout << "[Geometry3D]: ReadAABBs Done!"
            << "\n"
            << "Voxels size: " << _Voxels.size() << std::endl;
}

template <typename T>
template <typename LatSet>
void Geometry3D<T>::SetupBoundary(std::uint8_t AABBflag,
                                  std::uint8_t boundaryflag) {
  FlagField TransFlag(N, std::uint8_t(0));
  for (int id = 0; id < N; ++id) {
    if (GeometryFlag.get(id) == AABBflag) {
      Vector<T, 3>& voxel = _Voxels[id];
      unsigned int count = 0;
      for (unsigned int i = 1; i < LatSet::q; ++i) {
        Vector<T, 3> pt = voxel + LatSet::c[i] * _voxelSize;
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
void Geometry3D<T>::setFlag(const AABB<T, 3>& AABBs, std::uint8_t flag) {
  forEachVoxel(AABBs,
               [this, flag](int id) { GeometryFlag.SetField(id, flag); });
}

template <typename T>
void Geometry3D<T>::setFlag(const AABB<T, 3>& AABBs, std::uint8_t fromflag,
                            std::uint8_t flag) {
  forEachVoxel(AABBs, fromflag,
               [this, flag](int id) { GeometryFlag.SetField(id, flag); });
}

template <typename T>
void Geometry3D<T>::setFlag(const AABB<T, 3>& AABBs, const StlReader<T>& reader,
                            std::uint8_t flag) {
  // get cuboid defined by AABBs
  Vector<T, 3> min = AABBs.getMin();
  Vector<T, 3> max = AABBs.getMax();
  // get index of min and max in
  Vector<int, 3> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize),
                         int((min[2] - this->getMin()[2]) / _voxelSize));
  for (int i = 0; i < 3; ++i) {
    if (idx_min[i] < 0) idx_min[i] = 0;
  }
  Vector<int, 3> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize),
                         int((max[2] - this->getMin()[2]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  if (idx_max[2] > _Nz - 1) idx_max[2] = _Nz - 1;
  // set flag
  for (int k = idx_min[2]; k <= idx_max[2]; ++k) {
    for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
      for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
        int id = getIndex(Vector<int, 3>(i, j, k));
        Vector<T, 3>& voxel = _Voxels[id];
        Octree<T>* node = reader.getTree()->find(voxel);
        if (node != nullptr) {
          // check if it is a [leaf] node and if it is [inside]
          if (node->isLeaf() && node->getInside())
            GeometryFlag.SetField(id, flag);
        }
      }
    }
  }
}

template <typename T>
void Geometry3D<T>::setFlag(const AABB<T, 3>& AABBs, const StlReader<T>& reader,
                            std::uint8_t fromflag, std::uint8_t flag) {
  // get cuboid defined by AABBs
  Vector<T, 3> min = AABBs.getMin();
  Vector<T, 3> max = AABBs.getMax();
  // get index of min and max in
  Vector<int, 3> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize),
                         int((min[2] - this->getMin()[2]) / _voxelSize));
  for (int i = 0; i < 3; ++i) {
    if (idx_min[i] < 0) idx_min[i] = 0;
  }
  Vector<int, 3> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize),
                         int((max[2] - this->getMin()[2]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  if (idx_max[2] > _Nz - 1) idx_max[2] = _Nz - 1;
  // set flag
  for (int k = idx_min[2]; k <= idx_max[2]; ++k) {
    for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
      for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
        int id = getIndex(Vector<int, 3>(i, j, k));
        Vector<T, 3>& voxel = _Voxels[id];
        if (GeometryFlag.get(id) == fromflag) {
          Octree<T>* node = reader.getTree()->find(voxel);
          if (node != nullptr) {
            // check if it is a [leaf] node and if it is [inside]
            if (node->isLeaf() && node->getInside())
              GeometryFlag.SetField(id, flag);
          }
        }
      }
    }
  }
}

template <typename T>
template <typename Func>
void Geometry3D<T>::forEachVoxel(const AABB<T, 3>& AABBs, Func func) {
  // get cuboid defined by AABBs
  Vector<T, 3> min = AABBs.getMin();
  Vector<T, 3> max = AABBs.getMax();
  // get index of min and max in
  Vector<int, 3> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize),
                         int((min[2] - this->getMin()[2]) / _voxelSize));
  for (int i = 0; i < 3; ++i) {
    if (idx_min[i] < 0) idx_min[i] = 0;
  }
  Vector<int, 3> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize),
                         int((max[2] - this->getMin()[2]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  if (idx_max[2] > _Nz - 1) idx_max[2] = _Nz - 1;
  // set flag
  for (int k = idx_min[2]; k <= idx_max[2]; ++k) {
    for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
      for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
        int id = getIndex(Vector<int, 3>(i, j, k));
        const Vector<T, 3>& voxel = _Voxels[id];
        // check if voxel is inside AABBs
        if (AABBs.isInside(voxel)) {
          func(id);
        }
      }
    }
  }
}

template <typename T>
template <typename Func>
void Geometry3D<T>::forEachVoxel(const AABB<T, 3>& AABBs, std::uint8_t fromflag,
                                 Func func) {
  // get cuboid defined by AABBs
  Vector<T, 3> min = AABBs.getMin();
  Vector<T, 3> max = AABBs.getMax();
  // get index of min and max in
  Vector<int, 3> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize),
                         int((min[2] - this->getMin()[2]) / _voxelSize));
  for (int i = 0; i < 3; ++i) {
    if (idx_min[i] < 0) idx_min[i] = 0;
  }
  Vector<int, 3> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize),
                         int((max[2] - this->getMin()[2]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  if (idx_max[2] > _Nz - 1) idx_max[2] = _Nz - 1;
  // set flag
  for (int k = idx_min[2]; k <= idx_max[2]; ++k) {
    for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
      for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
        int id = getIndex(Vector<int, 3>(i, j, k));
        const Vector<T, 3>& voxel = _Voxels[id];
        // check if voxel is inside AABBs
        if (AABBs.isInside(voxel) && GeometryFlag.get(id) == fromflag) {
          func(id);
        }
      }
    }
  }
}

template <typename T>
template <typename Func>
void Geometry3D<T>::forEachVoxel(std::uint8_t fromflag, Func func) {
  for (int id = 0; id < N; ++id) {
    if (GeometryFlag.get(id) == fromflag) {
      func(id);
    }
  }
}