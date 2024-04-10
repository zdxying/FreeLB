#pragma once

#include "legacy/legacy_geometry/voxel_geometry3d.h"


template <typename T>
VoxelGeometry3D<T>::VoxelGeometry3D(const StlReader<T>& reader)
    : _voxelSize(reader.getVoxelSize()),
      AABB<T, 3>(reader.getMesh().getMin(), reader.getMesh().getMax()) {
  _Nx = int(this->_extension[0] / reader.getStlSize());
  _Ny = int(this->_extension[1] / reader.getStlSize());
  _Nz = int(this->_extension[2] / reader.getStlSize());

  _idx = Vector<int, 3>(1, _Nx, _Nx * _Ny);
  _GlobalIdx.resize(_Nx * _Ny * _Nz, -1);
  print();
  // read from octree
  ReadOctree(reader.getTree());
  // setNeighbor(*(this));
}

template <typename T>
VoxelGeometry3D<T>::VoxelGeometry3D(int Nx, int Ny, int Nz,
                                    const AABB<T, 3>& AABBs, T voxelSize,
                                    const Vector<T, 3>& min)
    : _Nx(Nx),
      _Ny(Ny),
      _Nz(Nz),
      _voxelSize(voxelSize),
      AABB<T, 3>(min, min + Vector<T, 3>{Nx * voxelSize, Ny * voxelSize,
                                         Nz * voxelSize}) {
  _idx = Vector<int, 3>(1, _Nx, _Nx * _Ny);
  _GlobalIdx.resize(_Nx * _Ny * _Nz, -1);
  print();
  // read from AABBs
  ReadAABBs(AABBs);
  // setNeighbor(AABBs);
}

template <typename T>
VoxelGeometry3D<T>::VoxelGeometry3D(int Nx, int Ny, int Nz,
                                    const AABB<T, 3>& AABBs,
                                    const StlReader<T>& reader, int flag,
                                    const Vector<T, 3>& min)
    : _Nx(Nx),
      _Ny(Ny),
      _Nz(Nz),
      _voxelSize(reader.getVoxelSize()),
      AABB<T, 3>(min, min + Vector<T, 3>{Nx * reader.getVoxelSize(),
                                         Ny * reader.getVoxelSize(),
                                         Nz * reader.getVoxelSize()}) {
  _idx = Vector<int, 3>(1, _Nx, _Nx * _Ny);
  _GlobalIdx.resize(_Nx * _Ny * _Nz, -1);
  print();
  ReadAABBs(AABBs);
  // set flag from stlreader
  setFlag(AABBs, reader, flag);
}

template <typename T>
void VoxelGeometry3D<T>::ReadOctree(Octree<T>* tree) {
  // loop over all sites in the AABB defined by (0,0,0) and (Nx,Ny,Nz)
  for (int k = 0; k < _Nz; ++k) {
    for (int j = 0; j < _Ny; ++j) {
      for (int i = 0; i < _Nx; ++i) {
        // get the location of the voxel
        Vector<T, 3> loc = _voxelSize * Vector<T, 3>{T(i), T(j), T(k)} +
                           this->_min + _voxelSize * T(0.5);
        // get the node containing the voxel
        Octree<T>* node = tree->find(loc);
        if (node != nullptr) {
          // check if it is a [leaf] node and if it is [inside]
          if (node->isLeaf() && node->getInside()) {
            // create a voxel at (i,j,k)
            Voxel<T, 3> voxel(loc, _Voxels.size());
            // set _GlobalIdx
            _GlobalIdx[i + j * _Nx + k * _Nx * _Ny] = _Voxels.size();
            // add to _Voxels
            _Voxels.push_back(voxel);
          }
        }
      }
    }
  }
  std::cout << "[VoxelGeometry3D]: ReadOctree Done!"
            << "\n"
            << "_Voxels size: " << _Voxels.size() << std::endl;
}

template <typename T>
void VoxelGeometry3D<T>::ReadAABBs(const AABB<T, 3>& AABBs) {
  // loop over all sites in the AABB defined by (0,0,0) and (Nx,Ny,Nz)
  for (int k = 0; k < _Nz; ++k) {
    for (int j = 0; j < _Ny; ++j) {
      for (int i = 0; i < _Nx; ++i) {
        // get the location of the voxel
        Vector<T, 3> loc = _voxelSize * Vector<T, 3>{T(i), T(j), T(k)} +
                           this->_min + _voxelSize * T(0.5);
        // check if it is inside
        if (AABBs.isInside(loc)) {
          // create a voxel at (i,j,k)
          Voxel<T, 3> voxel(loc, _Voxels.size());
          // set _GlobalIdx
          _GlobalIdx[i + j * _Nx + k * _Nx * _Ny] = _Voxels.size();
          // add to _Voxels
          _Voxels.push_back(voxel);
        }
      }
    }
  }
  std::cout << "[VoxelGeometry3D]: ReadAABBs Done!"
            << "\n"
            << "_Voxels size: " << _Voxels.size() << std::endl;
}

template <typename T>
template <typename LatSet>
void VoxelGeometry3D<T>::Setup() {
  // traverse _GlobalIdx
  int size = _GlobalIdx.size();
  for (int id = 0; id < size; ++id) {
    if (_GlobalIdx[id] != -1) {
      // get voxel
      Voxel<T, 3>& voxel = _Voxels[_GlobalIdx[id]];
      // init neighbor
      voxel.initNeighborList(LatSet::q);
      // get neighbor
      for (int i = 1; i < LatSet::q; ++i) {
        // check if neighbor is inside
        Vector<T, 3> pt = voxel.getLoc() + LatSet::c[i] * _voxelSize;
        if (this->isInside(pt)) {
          // check if neighbor index is inside
          int idx = LatSet::c[i] * _idx + id;
          if (idx >= 0 && idx < size) {
            // check if neighbor is initialized with index of _Voxels
            if (_GlobalIdx[idx] != -1) {
              voxel.setNeighbor(i, &_Voxels[_GlobalIdx[idx]]);
            }
          }
        }
      }
      // set voxel flag
      // is boundary voxel, set flag = 1
      // is inner voxel, set flag = 0
      if (voxel.getNeighborCount() < LatSet::q - 1)
        voxel.setFlag(1);
      else
        voxel.setFlag(0);
    }
  }
}

template <typename T>
template <typename LatSet>
void VoxelGeometry3D<T>::Setup(int AABBflag, int AABBBdflag, int stlflag,
                               int stlBdflag) {
  // traverse _GlobalIdx
  int size = _GlobalIdx.size();
  for (int id = 0; id < size; ++id) {
    if (_GlobalIdx[id] != -1) {
      // get voxel
      Voxel<T, 3>& voxel = _Voxels[_GlobalIdx[id]];
      // init neighbor
      voxel.initNeighborList(LatSet::q);
      // get neighbor
      for (int i = 1; i < LatSet::q; ++i) {
        // check if neighbor is inside
        Vector<T, 3> pt = voxel.getLoc() + LatSet::c[i] * _voxelSize;
        if (this->isInside(pt)) {
          // check if neighbor index is inside
          int idx = LatSet::c[i] * _idx + id;
          if (idx >= 0 && idx < size) {
            // check if neighbor is initialized with index of _Voxels
            if (_GlobalIdx[idx] != -1) {
              // set neighbor if has same flag
              if (voxel.getFlag() == _Voxels[_GlobalIdx[idx]].getFlag())
                voxel.setNeighbor(i, &_Voxels[_GlobalIdx[idx]]);
            }
          }
        }
      }
    }
  }
  for (int id = 0; id < size; ++id) {
    if (_GlobalIdx[id] != -1) {
      // get voxel
      Voxel<T, 3>& voxel = _Voxels[_GlobalIdx[id]];
      // set voxel flag
      if (voxel.getNeighborCount() < LatSet::q - 1) {
        if (voxel.getFlag() == AABBflag) {
          voxel.setFlag(AABBBdflag);
        } else if (voxel.getFlag() == stlflag) {
          voxel.setFlag(stlBdflag);
        }
      }
    }
  }
}
// this is an example of vtk file from :
// https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
// # vtk DataFile Version 2.0
// Volume example
// ASCII
// DATASET STRUCTURED_POINTS
// DIMENSIONS 3 4 6
// ASPECT_RATIO 1 1 1
// ORIGIN 0 0 0
// POINT_DATA 72
// SCALARS volume_scalars char 1
// LOOKUP_TABLE default
// ...

template <typename T>
template <bool boundaryinfo>
void VoxelGeometry3D<T>::WriteStruPoints() {
  // write to file
  DirCreator::Create_Dir("./vtkoutput");
  std::string fullName = "./vtkoutput/voxels_StructPoint.vtk";
  std::ofstream f(fullName.c_str());
  std::cout << "Writing " << _Voxels.size() << " Voxels to vtk file...";
  f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
  f << "Voxels" << std::endl << std::flush;
  f << "ASCII" << std::endl << std::flush;

  f << "DATASET STRUCTURED_POINTS" << std::endl << std::flush;
  f << "DIMENSIONS " << _Nx << " " << _Ny << " " << _Nz << std::endl
    << std::flush;
  f << "ORIGIN " << this->_min[0] + _voxelSize * T(0.5) << " "
    << this->_min[1] + _voxelSize * T(0.5) << " "
    << this->_min[2] + _voxelSize * T(0.5) << std::endl
    << std::flush;
  f << "ASPECT_RATIO " << _voxelSize << " " << _voxelSize << " " << _voxelSize
    << std::endl
    << std::flush;
  f << "POINT_DATA " << _Nx * _Ny * _Nz << std::endl << std::flush;

  // stream
  std::stringstream flags;
  flags << "SCALARS Flag int" << std::endl << std::flush;
  flags << "LOOKUP_TABLE default" << std::endl << std::flush;
  for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
    if (_GlobalIdx[i] != -1)
      flags << _Voxels[_GlobalIdx[i]].getFlag() << " ";
    else
      flags << -1 << " ";
  }
  flags << std::endl;
  f << flags.str();
  // boundary info
  if constexpr (boundaryinfo) {
    std::stringstream boundary;
    boundary << "SCALARS BoundaryNum int" << std::endl << std::flush;
    boundary << "LOOKUP_TABLE default" << std::endl << std::flush;
    for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
      if (_GlobalIdx[i] != -1)
        boundary << _Voxels[_GlobalIdx[i]].getNeighborCount() << " ";
      else
        boundary << -1 << " ";
    }
    boundary << std::endl;
    f << boundary.str();
  }
  f.close();
  std::cout << " Done!" << std::endl;
}

template <typename T>
void VoxelGeometry3D<T>::WriteStruGrid() {
  // write to file
  DirCreator::Create_Dir("./vtkoutput");
  std::string fullName = "./vtkoutput/voxels_StructGrid.vtk";
  std::ofstream f(fullName.c_str());
  std::cout << "Writing " << _Voxels.size() << " Voxels to vtk file...";
  f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
  f << "Voxels" << std::endl << std::flush;
  f << "ASCII" << std::endl << std::flush;

  f << "DATASET STRUCTURED_GRID" << std::endl << std::flush;
  f << "DIMENSIONS " << _Nx << " " << _Ny << " " << _Nz << std::endl
    << std::flush;
  f << "ORIGIN " << this->_min[0] + _voxelSize * T(0.5) << " "
    << this->_min[1] + _voxelSize * T(0.5) << " "
    << this->_min[2] + _voxelSize * T(0.5) << std::endl
    << std::flush;
  f << "ASPECT_RATIO " << _voxelSize << " " << _voxelSize << " " << _voxelSize
    << std::endl
    << std::flush;
  f << "POINT_DATA " << _Nx * _Ny * _Nz << std::endl << std::flush;
  f << "SCALARS Flag int" << std::endl << std::flush;
  f << "LOOKUP_TABLE default" << std::endl << std::flush;

  std::stringstream flags;
  for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
    if (_GlobalIdx[i] != -1) {
      flags << _Voxels[_GlobalIdx[i]].getFlag() << " ";
    } else {
      flags << -1 << " ";
    }
  }
  f << flags.str() << std::endl;
  f.close();
  std::cout << " Done!" << std::endl;
}

template <typename T>
void VoxelGeometry3D<T>::setFlag(const AABB<T, 3>& AABBs, int flag) {
  // get cuboid defined by AABBs
  Vector<T, 3> min = AABBs.getMin();
  Vector<T, 3> max = AABBs.getMax();
  // get index of min and max in _GlobalIdx
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
        int id = Vector<int, 3>(i, j, k) * _idx;
        // check if there is a voxel
        if (_GlobalIdx[id] != -1) {
          Voxel<T, 3>& voxel = _Voxels[_GlobalIdx[id]];
          // check if voxel is inside AABBs
          if (AABBs.isInside(voxel.getLoc())) voxel.setFlag(flag);
        }
      }
    }
  }
}

template <typename T>
void VoxelGeometry3D<T>::setFlag(const AABB<T, 3>& AABBs, int fromflag,
                                 int flag) {
  // get cuboid defined by AABBs
  Vector<T, 3> min = AABBs.getMin();
  Vector<T, 3> max = AABBs.getMax();
  // get index of min and max in _GlobalIdx
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
        int id = Vector<int, 3>(i, j, k) * _idx;
        if (_GlobalIdx[id] != -1) {
          Voxel<T, 3>& voxel = _Voxels[_GlobalIdx[id]];
          // check if voxel is inside AABBs
          if (AABBs.isInside(voxel.getLoc()) && voxel.getFlag() == fromflag)
            voxel.setFlag(flag);
        }
      }
    }
  }
}

template <typename T>
void VoxelGeometry3D<T>::setFlag(const AABB<T, 3>& AABBs,
                                 const StlReader<T>& reader, int flag) {
  // get cuboid defined by AABBs
  Vector<T, 3> min = AABBs.getMin();
  Vector<T, 3> max = AABBs.getMax();
  // get index of min and max in _GlobalIdx
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
        int id = Vector<int, 3>(i, j, k) * _idx;
        if (_GlobalIdx[id] != -1) {
          Voxel<T, 3>& voxel = _Voxels[_GlobalIdx[id]];
          Octree<T>* node = reader.getTree()->find(voxel.getLoc());
          if (node != nullptr) {
            // check if it is a [leaf] node and if it is [inside]
            if (node->isLeaf() && node->getInside()) voxel.setFlag(flag);
          }
        }
      }
    }
  }
}

template <typename T>
void VoxelGeometry3D<T>::setFlag(const AABB<T, 3>& AABBs,
                                 const StlReader<T>& reader, int fromflag,
                                 int flag) {
  // get cuboid defined by AABBs
  Vector<T, 3> min = AABBs.getMin();
  Vector<T, 3> max = AABBs.getMax();
  // get index of min and max in _GlobalIdx
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
        int id = Vector<int, 3>(i, j, k) * _idx;
        if (_GlobalIdx[id] != -1) {
          Voxel<T, 3>& voxel = _Voxels[_GlobalIdx[id]];
          if (voxel.getFlag() == fromflag) {
            Octree<T>* node = reader.getTree()->find(voxel.getLoc());
            if (node != nullptr) {
              // check if it is a [leaf] node and if it is [inside]
              if (node->isLeaf() && node->getInside()) voxel.setFlag(flag);
            }
          }
        }
      }
    }
  }
}