// legacyfield3D.h
#pragma once

#include <vector>

#include "geometry/geometry3d.h"
#include "lbm/unit_converter.h"
#include "utils/fdm_solver.h"
#include "utils/directories.h"
#include "utils/util.h"

template <typename T>
class VelocityField3D {
 private:
  // velocity field
  std::vector<Vector<T, 3>> _Velocity;
  // velocity field init, lattice unit
  Vector<T, 3> _Velocity_Init;
  // geometry
  VoxelGeometry3D<T> &_Geo;
  // converter
  BaseConverter<T> &_Conv;
  // index to get global index of Velocity
  // std::vector<int> &_GlobalIdx;
 public:
  // constructor
  VelocityField3D(const Vector<T, 3> &latu_init, BaseConverter<T> &conv,
                  VoxelGeometry3D<T> &Geo3D)
      : _Velocity_Init(latu_init), _Geo(Geo3D), _Conv(conv) {
    // set all velocity to init
    _Velocity.clear();
    _Velocity.resize(_Geo.getVoxels().size(), _Velocity_Init);
    // set boundary velocity to (0,0,0)
    setVelocity(1, Vector<T, 3>(T(0), T(0), T(0)));
    std::cout << "[VelocityField3D]: "
              << "U_Init = (" << _Velocity_Init[0] << ", " << _Velocity_Init[1]
              << ", " << _Velocity_Init[2] << ")" << std::endl;
  }
  VelocityField3D(BaseConverter<T> &conv, const Vector<T, 3> &physu_init,
                  VoxelGeometry3D<T> &Geo3D)
      : VelocityField3D(conv.getLatticeU(physu_init), conv, Geo3D) {}
  // get
  VoxelGeometry3D<T> &getGeo() { return _Geo; }
  // get velocity with lattice unit
  Vector<T, 3> &getVelocity(int id) { return _Velocity[id]; }
  const Vector<T, 3> &getVelocity(int id) const { return _Velocity[id]; }
  Vector<T, 3> *getVelocityPtr(int id) { return &_Velocity[id]; }
  // get velocity with physical unit
  T getPhysVelocityU0(int id) const { return _Conv.getPhysU(_Velocity[id][0]); }
  T getPhysVelocityU1(int id) const { return _Conv.getPhysU(_Velocity[id][1]); }
  T getPhysVelocityU2(int id) const { return _Conv.getPhysU(_Velocity[id][2]); }

  // set
  void reset(int id) { _Velocity[id] = _Velocity_Init; }
  void setCellVelocity(int Id, const Vector<T, 3> &latu) {
    _Velocity[Id] = latu;
  }
  void setCellVelocity(const Vector<T, 3> &physu, int Id) {
    _Velocity[Id] = _Conv.getLatticeU(physu);
  }
  // set all velocity with GeoFlag = flag to vel
  void setVelocity(const Vector<T, 3> &physu, int flag) {
    setVelocity(flag, _Conv.getLatticeU(physu));
  }
  void setVelocity(int flag, const Vector<T, 3> &latu) {
    for (int i = 0; i < _Velocity.size(); ++i) {
      if (_Geo.getVoxel(i).getFlag() == flag) {
        _Velocity[i] = latu;
      }
    }
  }
  void setVelocity(const Vector<T, 3> &physu, const AABB<T, 3> &AABBs) {
    setVelocity(AABBs, _Conv.getLatticeU(physu));
  }
  void setVelocity(const AABB<T, 3> &AABBs, const Vector<T, 3> &latu) {
    // get cuboid defined by AABBs
    Vector<T, 3> min = AABBs.getMin();
    Vector<T, 3> max = AABBs.getMax();
    // get index of min and max in _GlobalIdx
    Vector<int, 3> idx_min(
        int((min[0] - _Geo.getMin()[0]) / _Geo.getVoxelSize()),
        int((min[1] - _Geo.getMin()[1]) / _Geo.getVoxelSize()),
        int((min[2] - _Geo.getMin()[2]) / _Geo.getVoxelSize()));
    for (int i = 0; i < 3; ++i) {
      if (idx_min[i] < 0) idx_min[i] = 0;
    }
    Vector<int, 3> idx_max(
        int((max[0] - _Geo.getMin()[0]) / _Geo.getVoxelSize()),
        int((max[1] - _Geo.getMin()[1]) / _Geo.getVoxelSize()),
        int((max[2] - _Geo.getMin()[2]) / _Geo.getVoxelSize()));
    if (idx_max[0] > _Geo.getNx() - 1) idx_max[0] = _Geo.getNx() - 1;
    if (idx_max[1] > _Geo.getNy() - 1) idx_max[1] = _Geo.getNy() - 1;
    if (idx_max[2] > _Geo.getNz() - 1) idx_max[2] = _Geo.getNz() - 1;
    // set velocity
    for (int k = idx_min[2]; k <= idx_max[2]; ++k) {
      for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
        for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
          int id = Vector<int, 3>(i, j, k) * _Geo.getIdx();
          if (_Geo.getGlobalIdx()[id] != -1) {
            Voxel<T, 3> &voxel = _Geo.getVoxel(id);
            // check if voxel is inside AABBs
            if (AABBs.isInside(voxel.getLoc()))
              _Velocity[_Geo.getGlobalIdx()[id]] = latu;
          }
        }
      }
    }
  }
  void setVelocity(const Vector<T, 3> &physu, const AABB<T, 3> &AABBs,
                   int flag) {
    setVelocity(AABBs, _Conv.getLatticeU(physu), flag);
  }
  void setVelocity(const AABB<T, 3> &AABBs, const Vector<T, 3> &latu,
                   int flag) {
    // get cuboid defined by AABBs
    Vector<T, 3> min = AABBs.getMin();
    Vector<T, 3> max = AABBs.getMax();
    // get index of min and max in _GlobalIdx
    Vector<int, 3> idx_min(
        int((min[0] - _Geo.getMin()[0]) / _Geo.getVoxelSize()),
        int((min[1] - _Geo.getMin()[1]) / _Geo.getVoxelSize()),
        int((min[2] - _Geo.getMin()[2]) / _Geo.getVoxelSize()));
    for (int i = 0; i < 3; ++i) {
      if (idx_min[i] < 0) idx_min[i] = 0;
    }
    Vector<int, 3> idx_max(
        int((max[0] - _Geo.getMin()[0]) / _Geo.getVoxelSize()),
        int((max[1] - _Geo.getMin()[1]) / _Geo.getVoxelSize()),
        int((max[2] - _Geo.getMin()[2]) / _Geo.getVoxelSize()));
    if (idx_max[0] > _Geo.getNx() - 1) idx_max[0] = _Geo.getNx() - 1;
    if (idx_max[1] > _Geo.getNy() - 1) idx_max[1] = _Geo.getNy() - 1;
    if (idx_max[2] > _Geo.getNz() - 1) idx_max[2] = _Geo.getNz() - 1;
    // set velocity
    for (int k = idx_min[2]; k <= idx_max[2]; ++k) {
      for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
        for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
          int id = Vector<int, 3>(i, j, k) * _Geo.getIdx();
          if (_Geo.getGlobalIdx()[id] != -1) {
            Voxel<T, 3> &voxel = _Geo.getVoxel(id);
            // check if voxel is inside AABBs
            if (AABBs.isInside(voxel.getLoc()) && voxel.getFlag() == flag)
              _Velocity[_Geo.getGlobalIdx()[id]] = latu;
          }
        }
      }
    }
  }
};

// CA field
template <typename T>
class CAZSField3D {
 private:
  VoxelGeometry3D<T> &_Geo;
  VelocityField3D<T> &Velocity;
  // ZS field data
  std::vector<int> State;
  std::vector<T> Fs;
  // normalized gradient of solid fraction: nx, ny, nz
  std::vector<T> nGradFs_x;
  std::vector<T> nGradFs_y;
  std::vector<T> nGradFs_z;
  // anisotropy funvtion:
  // norm(gradFs) * (1 - 3eps + 4eps * (nx^4 + ny^4 + nz^4))
  std::vector<T> Anisotropy;

 public:
  CAZSField3D(VelocityField3D<T> &velocity)
      : Velocity(velocity), _Geo(velocity.getGeo()) {
    State.resize(_Geo.getVoxels().size(), 0);
    Fs.resize(_Geo.getVoxels().size(), T(0));
    nGradFs_x.resize(_Geo.getVoxels().size(), T(0));
    nGradFs_y.resize(_Geo.getVoxels().size(), T(0));
    nGradFs_z.resize(_Geo.getVoxels().size(), T(0));
    Anisotropy.resize(_Geo.getVoxels().size(), T(0));
  }
  void SetState(int flag, int toState = -1) {
    const std::vector<Voxel<T, 3>> &Voxels = _Geo.getVoxels();
    for (int i = 0; i < Voxels.size(); ++i) {
      if (Voxels[i].getFlag() == flag) State[i] = toState;
    }
  }
  inline void Activate(int id) {
    State[id] = 1;
    // _Geo.getVoxel(id).setFlag(1);
    // set velocity to 0
    Velocity.getVelocity(id).clear();
  }
  inline void Deactivate(int id) {
    State[id] = -1;
    // _Geo.getVoxel(id).setFlag(-1);
    // set velocity to 0
    // Velocity.getVelocity(id).clear();
  }
  void AddSolidFraction(int id, T sf) { Fs[id] += sf; }
  void SetSolidFraction(int id, T sf) { Fs[id] = sf; }
  void UpdatenGradFs(T Epsilon = T(0.03)) {
    const std::vector<Voxel<T, 3>> &Voxels = _Geo.getVoxels();
    // gradient of solid fraction
    Vector<T, 3> gradFs = Vector<T, 3>(T(0), T(0), T(0));
    T normgradFs = T(0);
    int id = 0;
// traverse all voxels
// handle /0
#pragma omp parallel for private(gradFs, normgradFs, id) num_threads(Thread_Num)
    for (const Voxel<T, 3> &voxel : Voxels) {
      id = voxel.getId();
      if (voxel.getFlag() == 0) {
        gradFs[0] = (Fs[voxel.getNeighbor(1)->getId()] -
                     Fs[voxel.getNeighbor(2)->getId()]) /
                    2;
        gradFs[1] = (Fs[voxel.getNeighbor(3)->getId()] -
                     Fs[voxel.getNeighbor(4)->getId()]) /
                    2;
        gradFs[2] = (Fs[voxel.getNeighbor(5)->getId()] -
                     Fs[voxel.getNeighbor(6)->getId()]) /
                    2;
        // get normalized gradient of solid fraction
        normgradFs = gradFs.getnorm();
        if (normgradFs < std::numeric_limits<T>::epsilon()) {
          nGradFs_x[id] = T(0);
          nGradFs_y[id] = T(0);
          nGradFs_z[id] = T(0);
        } else {
          nGradFs_x[id] = gradFs[0] / normgradFs;
          nGradFs_y[id] = gradFs[1] / normgradFs;
          nGradFs_z[id] = gradFs[2] / normgradFs;
        }
      } else {
        nGradFs_x[id] = T(0);
        nGradFs_y[id] = T(0);
        nGradFs_z[id] = T(0);
      }
      // calculate anisotropy function
      // norm(gradFs) * (1 - 3eps + 4eps * (nx^4 + ny^4 + nz^4))
      Anisotropy[id] =
          normgradFs * (T(1) - 3 * Epsilon +
                        4 * Epsilon *
                            (pow(nGradFs_x[id], 4) + pow(nGradFs_y[id], 4) +
                             pow(nGradFs_z[id], 4)));
    }
  }
  int getState(int id) const { return State[id]; }
  T getSolidFraction(int id) const { return Fs[id]; }
  T getnGradFs_x(int id) const { return nGradFs_x[id]; }
  T getnGradFs_y(int id) const { return nGradFs_y[id]; }
  T getnGradFs_z(int id) const { return nGradFs_z[id]; }

  VoxelGeometry3D<T> &getGeo() { return _Geo; }
  VelocityField3D<T> &getVelocity() { return Velocity; }
  std::vector<int> &getStates() { return State; }
  std::vector<T> &getSolidFraction() { return Fs; }
  std::vector<T> &getnGradFs_x() { return nGradFs_x; }
  std::vector<T> &getnGradFs_y() { return nGradFs_y; }
  std::vector<T> &getnGradFs_z() { return nGradFs_z; }
  std::vector<T> &getAnisotropy() { return Anisotropy; }

  void WriteStates(int step) const {
    int _Nx = _Geo.getNx();
    int _Ny = _Geo.getNy();
    int _Nz = _Geo.getNz();
    const std::vector<int> &_GlobalIdx = _Geo.getGlobalIdx();
    // write to file
    DirCreator::Create_Dir("./vtkoutput");
    std::string fullName = "./vtkoutput/State" + std::to_string(step) + ".vtk";
    std::ofstream f(fullName.c_str());
    f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
    f << "Voxels" << std::endl << std::flush;
    f << "ASCII" << std::endl << std::flush;

    f << "DATASET STRUCTURED_POINTS" << std::endl << std::flush;
    f << "DIMENSIONS " << _Nx << " " << _Ny << " " << _Nz << std::endl
      << std::flush;
    f << "ORIGIN " << _Geo.getMin()[0] + _Geo.getVoxelSize() * T(0.5) << " "
      << _Geo.getMin()[1] + _Geo.getVoxelSize() * T(0.5) << " "
      << _Geo.getMin()[2] + _Geo.getVoxelSize() * T(0.5) << std::endl
      << std::flush;
    f << "ASPECT_RATIO " << _Geo.getVoxelSize() << " " << _Geo.getVoxelSize()
      << " " << _Geo.getVoxelSize() << std::endl
      << std::flush;
    f << "POINT_DATA " << _Nx * _Ny * _Nz << std::endl << std::flush;

    std::stringstream state;
    state << "SCALARS state int" << std::endl << std::flush;
    state << "LOOKUP_TABLE default" << std::endl << std::flush;

    for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
      if (_GlobalIdx[i] != -1) {
        state << State[i] << " ";
      } else {
        state << -2 << " ";
      }
    }
    state << std::endl;
    f << state.str();
    f.close();
  }
  void WriteFs(int step) const {
    int _Nx = _Geo.getNx();
    int _Ny = _Geo.getNy();
    int _Nz = _Geo.getNz();
    const std::vector<int> &_GlobalIdx = _Geo.getGlobalIdx();
    // write to file
    DirCreator::Create_Dir("./vtkoutput");
    std::string fullName = "./vtkoutput/Fs" + std::to_string(step) + ".vtk";
    std::ofstream f(fullName.c_str());
    f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
    f << "Voxels" << std::endl << std::flush;
    f << "ASCII" << std::endl << std::flush;

    f << "DATASET STRUCTURED_POINTS" << std::endl << std::flush;
    f << "DIMENSIONS " << _Nx << " " << _Ny << " " << _Nz << std::endl
      << std::flush;
    f << "ORIGIN " << _Geo.getMin()[0] + _Geo.getVoxelSize() * T(0.5) << " "
      << _Geo.getMin()[1] + _Geo.getVoxelSize() * T(0.5) << " "
      << _Geo.getMin()[2] + _Geo.getVoxelSize() * T(0.5) << std::endl
      << std::flush;
    f << "ASPECT_RATIO " << _Geo.getVoxelSize() << " " << _Geo.getVoxelSize()
      << " " << _Geo.getVoxelSize() << std::endl
      << std::flush;
    f << "POINT_DATA " << _Nx * _Ny * _Nz << std::endl << std::flush;

    std::stringstream fs;
    if constexpr (std::is_same<T, double>::value) {
      fs << "SCALARS Fs double" << std::endl << std::flush;
    } else {
      fs << "SCALARS Fs float" << std::endl << std::flush;
    }
    fs << "LOOKUP_TABLE default" << std::endl << std::flush;
    for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
      if (_GlobalIdx[i] != -1) {
        fs << Fs[i] << " ";
      } else {
        fs << 0 << " ";
      }
    }
    fs << std::endl;
    f << fs.str();
    f.close();
  }
};