
// voxel_geometry2d.hh

#pragma once

#include "legacy/legacy_geometry/voxel_geometry2d.h"

template <typename T>
VoxelGeometry2D<T>::VoxelGeometry2D(int Nx, int Ny, const AABB<T, 2> &AABBs,
                                    T voxelSize, const Vector<T, 2> &min)
    : _Nx(Nx),
      _Ny(Ny),
      _voxelSize(voxelSize),
      AABB<T, 2>(min, min + Vector<T, 2>{Nx * voxelSize, Ny * voxelSize}) {
  _idx = Vector<int, 2>(1, _Nx);
  _GlobalIdx.resize(_Nx * _Ny, -1);
  print();
  // read from AABBs
  ReadAABBs(AABBs);
  // setNeighbor(AABBs);
}

template <typename T>
void VoxelGeometry2D<T>::ReadAABBs(const AABB<T, 2> &AABBs) {
  // loop over all sites in the AABB
  for (int j = 0; j < _Ny; ++j) {
    for (int i = 0; i < _Nx; ++i) {
      // get the location of the voxel
      Vector<T, 2> loc = _voxelSize * Vector<T, 2>{T(i), T(j)} + this->_min +
                         _voxelSize * T(0.5);
      // check if it is inside
      if (AABBs.isInside(loc)) {
        // create a voxel at (i,j,k)
        Voxel<T, 2> voxel(loc, _Voxels.size());
        // set _GlobalIdx
        _GlobalIdx[i + j * _Nx] = _Voxels.size();
        // add to _Voxels
        _Voxels.push_back(voxel);
      }
    }
  }
  std::cout << "[VoxelGeometry2D]: ReadAABBs Done!"
            << "\n"
            << "_Voxels size: " << _Voxels.size() << std::endl;
}

template <typename T>
template <typename LatSet>
void VoxelGeometry2D<T>::Setup(int innerflag, int boundaryflag) {
  InnerFlag = innerflag;
  BoundaryFlag = boundaryflag;
  // traverse _GlobalIdx
  int size = _GlobalIdx.size();
  for (int id = 0; id < size; ++id) {
    if (_GlobalIdx[id] != -1) {
      // get voxel
      Voxel<T, 2> &voxel = _Voxels[_GlobalIdx[id]];
      // init neighbor
      voxel.initNeighborList(LatSet::q);
      // get neighbor
      for (int i = 1; i < LatSet::q; ++i) {
        // check if neighbor is inside
        Vector<T, 2> pt = voxel.getLoc() + LatSet::c[i] * _voxelSize;
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
      if (voxel.getNeighborCount() < LatSet::q - 1)
        voxel.setFlag(BoundaryFlag);
      else
        voxel.setFlag(InnerFlag);
    }
  }
}

// template <typename T>
// template <typename LatSet>
// void VoxelGeometry2D<T>::Setup(int fromflag, int toBdflag) {
//   // traverse _GlobalIdx
//   int size = _GlobalIdx.size();
//   for (int id = 0; id < size; ++id) {
//     if (_GlobalIdx[id] != -1) {
//       // get voxel
//       Voxel<T, 2> &voxel = _Voxels[_GlobalIdx[id]];
//       // init neighbor
//       voxel.initNeighborList(LatSet::q);
//       // get neighbor
//       for (int i = 1; i < LatSet::q; ++i) {
//         // check if neighbor is inside
//         Vector<T, 2> pt = voxel.getLoc() + LatSet::c[i] * _voxelSize;
//         if (this->isInside(pt)) {
//           // check if neighbor index is inside
//           int idx = LatSet::c[i] * _idx + id;
//           if (idx >= 0 && idx < size) {
//             // check if neighbor is initialized with index of _Voxels
//             if (_GlobalIdx[idx] != -1) {
//               // set neighbor if has same flag // NOT A GOOD IDEA
//               if (voxel.getFlag() == _Voxels[_GlobalIdx[idx]].getFlag())
//                 voxel.setNeighbor(i, &_Voxels[_GlobalIdx[idx]]);
//             }
//           }
//         }
//       }
//     }
//   }
//   for (int id = 0; id < size; ++id) {
//     if (_GlobalIdx[id] != -1) {
//       // get voxel
//       Voxel<T, 2> &voxel = _Voxels[_GlobalIdx[id]];
//       // set voxel flag
//       if (voxel.getNeighborCount() < LatSet::q - 1) {
//         if (voxel.getFlag() == fromflag) {
//           voxel.setFlag(toBdflag);
//         }
//       }
//     }
//   }
// }

template <typename T>
void VoxelGeometry2D<T>::setFlag(const AABB<T, 2> &AABBs, int flag) {
  // get cuboid defined by AABBs
  Vector<T, 2> min = AABBs.getMin();
  Vector<T, 2> max = AABBs.getMax();
  // get index of min and max in _GlobalIdx
  Vector<int, 2> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize));
  for (int i = 0; i < 2; ++i) {
    if (idx_min[i] < 0) idx_min[i] = 0;
  }
  Vector<int, 2> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  // set flag
  for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
    for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
      int id = Vector<int, 2>(i, j) * _idx;
      // check if there is a voxel
      if (_GlobalIdx[id] != -1) {
        Voxel<T, 2> &voxel = _Voxels[_GlobalIdx[id]];
        // check if voxel is inside AABBs
        if (AABBs.isInside(voxel.getLoc())) voxel.setFlag(flag);
      }
    }
  }
}

template <typename T>
void VoxelGeometry2D<T>::setFlag(const AABB<T, 2> &AABBs, int fromflag,
                                 int toflag) {
  // get cuboid defined by AABBs
  Vector<T, 2> min = AABBs.getMin();
  Vector<T, 2> max = AABBs.getMax();
  // get index of min and max in _GlobalIdx
  Vector<int, 2> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize));
  for (int i = 0; i < 2; ++i) {
    if (idx_min[i] < 0) idx_min[i] = 0;
  }
  Vector<int, 2> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  // set flag
  for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
    for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
      int id = Vector<int, 2>(i, j) * _idx;
      if (_GlobalIdx[id] != -1) {
        Voxel<T, 2> &voxel = _Voxels[_GlobalIdx[id]];
        // check if voxel is inside AABBs
        if (AABBs.isInside(voxel.getLoc()) && voxel.getFlag() == fromflag)
          voxel.setFlag(toflag);
      }
    }
  }
}

template <typename T>
template <typename U>
void VoxelGeometry2D<T>::setFieldFlag(const AABB<T, 2> &AABBs,
                                      std::vector<U> &field, const U &flag) {
  // get cuboid defined by AABBs
  Vector<T, 2> min = AABBs.getMin();
  Vector<T, 2> max = AABBs.getMax();
  // get index of min and max in _GlobalIdx
  Vector<int, 2> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize));
  for (int i = 0; i < 2; ++i) {
    if (idx_min[i] < 0) idx_min[i] = 0;
  }
  Vector<int, 2> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
// set flag
#pragma omp parallel for num_threads(Thread_Num) collapse(2) schedule(static)
  for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
    for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
      int id = Vector<int, 2>(i, j) * _idx;
      if (_GlobalIdx[id] != -1) {
        Voxel<T, 2> &voxel = _Voxels[_GlobalIdx[id]];
        // check if voxel is inside AABBs
        if (AABBs.isInside(voxel.getLoc())) field[voxel.getId()] = flag;
      }
    }
  }
}

template <typename T>
template <typename U>
void VoxelGeometry2D<T>::setFieldFlag(const AABB<T, 2> &AABBs,
                                      std::vector<U> &field, const U &fromflag,
                                      const U &toflag) {
  // get cuboid defined by AABBs
  Vector<T, 2> min = AABBs.getMin();
  Vector<T, 2> max = AABBs.getMax();
  // get index of min and max in _GlobalIdx
  Vector<int, 2> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize));
  for (int i = 0; i < 2; ++i) {
    if (idx_min[i] < 0) idx_min[i] = 0;
  }
  Vector<int, 2> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  // set flag
  for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
    for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
      int id = Vector<int, 2>(i, j) * _idx;
      if (_GlobalIdx[id] != -1) {
        Voxel<T, 2> &voxel = _Voxels[_GlobalIdx[id]];
        // check if voxel is inside AABBs
        if (AABBs.isInside(voxel.getLoc()) && field[voxel.getId()] == fromflag)
          field[voxel.getId()] = toflag;
      }
    }
  }
}

template <typename T>
template <typename U>
void VoxelGeometry2D<T>::setFieldFlag(const AABB<T, 2> &AABBs,
                                      std::vector<U> &field, const U &flag,
                                      int voxflag) {
  // get cuboid defined by AABBs
  Vector<T, 2> min = AABBs.getMin();
  Vector<T, 2> max = AABBs.getMax();
  // get index of min and max in _GlobalIdx
  Vector<int, 2> idx_min(int((min[0] - this->getMin()[0]) / _voxelSize),
                         int((min[1] - this->getMin()[1]) / _voxelSize));
  for (int i = 0; i < 2; ++i) {
    if (idx_min[i] < 0) idx_min[i] = 0;
  }
  Vector<int, 2> idx_max(int((max[0] - this->getMin()[0]) / _voxelSize),
                         int((max[1] - this->getMin()[1]) / _voxelSize));
  if (idx_max[0] > _Nx - 1) idx_max[0] = _Nx - 1;
  if (idx_max[1] > _Ny - 1) idx_max[1] = _Ny - 1;
  // set flag
  for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
    for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
      int id = Vector<int, 2>(i, j) * _idx;
      if (_GlobalIdx[id] != -1) {
        Voxel<T, 2> &voxel = _Voxels[_GlobalIdx[id]];
        // check if voxel is inside AABBs
        if (AABBs.isInside(voxel.getLoc()) && voxel.getFlag() == voxflag)
          field[voxel.getId()] = flag;
      }
    }
  }
}

template <typename T>
template <bool boundaryinfo>
void VoxelGeometry2D<T>::WriteStruPoints() {
  // write to file
  DirCreator::Create_Dir("./vtkoutput");
  std::string fullName = "./vtkoutput/voxels_StructPoint.vtk";
  std::ofstream f(fullName.c_str());
  std::cout << "Writing " << _Voxels.size() << " Voxels to vtk file...";
  f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
  f << "Voxels" << std::endl << std::flush;
  f << "ASCII" << std::endl << std::flush;

  f << "DATASET STRUCTURED_POINTS" << std::endl << std::flush;
  f << "DIMENSIONS " << _Nx << " " << _Ny << " 1" << std::endl << std::flush;
  f << "ORIGIN " << this->_min[0] + _voxelSize * T(0.5) << " "
    << this->_min[1] + _voxelSize * T(0.5) << " 0" << std::endl
    << std::flush;
  f << "ASPECT_RATIO " << _voxelSize << " " << _voxelSize << " 1" << std::endl
    << std::flush;
  f << "POINT_DATA " << _Nx * _Ny << std::endl << std::flush;

  // stream
  std::stringstream flags;
  flags << "SCALARS Flag int" << std::endl << std::flush;
  flags << "LOOKUP_TABLE default" << std::endl << std::flush;
  for (int i = 0; i < _Nx * _Ny; ++i) {
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
    for (int i = 0; i < _Nx * _Ny; ++i) {
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