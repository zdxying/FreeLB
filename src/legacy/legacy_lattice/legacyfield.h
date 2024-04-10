// lattice
#pragma once

#include <vector>

#include "geometry/geometry2d.h"
#include "lbm/unit_converter.h"
#include "utils/util.h"

template <typename T>
class VelocityField2D {
 private:
  // velocity field
  std::vector<Vector<T, 2>> _Velocity;
  // velocity field init, lattice unit
  Vector<T, 2> _Velocity_Init;
  // geometry
  VoxelGeometry2D<T> &_Geo;
  // converter
  BaseConverter<T> &_Conv;
  // index to get global index of Velocity
  // std::vector<int> &_GlobalIdx;
 public:
  // constructor
  VelocityField2D(const Vector<T, 2> &latu_init, BaseConverter<T> &conv,
                  VoxelGeometry2D<T> &Geo2D)
      : _Velocity_Init(latu_init), _Geo(Geo2D), _Conv(conv) {
    // set all velocity to init
    _Velocity.clear();
    _Velocity.resize(_Geo.getVoxels().size(), _Velocity_Init);
    // set boundary velocity to (0,0,0)
    setVelocity(1, Vector<T, 2>{T(0), T(0)});
    std::cout << "[VelocityField2D]: "
              << "U_Init = (" << _Velocity_Init[0] << ", " << _Velocity_Init[1]
              << ")" << std::endl;
  }
  VelocityField2D(BaseConverter<T> &conv, const Vector<T, 2> &physu_init,
                  VoxelGeometry2D<T> &Geo2D)
      : VelocityField2D(conv.getLatticeU(physu_init), conv, Geo2D) {}
  // get
  VoxelGeometry2D<T> &getGeo() { return _Geo; }
  // get velocity with lattice unit
  Vector<T, 2> &getVelocity(int id) { return _Velocity[id]; }
  const Vector<T, 2> &getVelocity(int id) const { return _Velocity[id]; }
  Vector<T, 2> *getVelocityPtr(int id) { return &_Velocity[id]; }
  std::vector<Vector<T, 2>> &getVelocity() { return _Velocity; }
  const std::vector<Vector<T, 2>> &getVelocity() const { return _Velocity; }
  // get velocity with physical unit
  T getPhysVelocityU0(int id) const { return _Conv.getPhysU(_Velocity[id][0]); }
  T getPhysVelocityU1(int id) const { return _Conv.getPhysU(_Velocity[id][1]); }

  // set
  void reset(int id) { _Velocity[id] = _Velocity_Init; }
  void setCellVelocity(int Id, const Vector<T, 2> &latu) {
    _Velocity[Id] = latu;
  }
  void setCellVelocity(const Vector<T, 2> &physu, int Id) {
    _Velocity[Id] = _Conv.getLatticeU(physu);
  }
  // set all velocity with GeoFlag = flag to vel
  void setVelocity(const Vector<T, 2> &physu, int flag) {
    setVelocity(flag, _Conv.getLatticeU(physu));
  }
  void setVelocity(int flag, const Vector<T, 2> &latu) {
    for (int i = 0; i < _Velocity.size(); ++i) {
      if (_Geo.getVoxel(i).getFlag() == flag) {
        _Velocity[i] = latu;
      }
    }
  }
  void setVelocity(const Vector<T, 2> &physu, const AABB<T, 2> &AABBs) {
    setVelocity(AABBs, _Conv.getLatticeU(physu));
  }
  void setVelocity(const AABB<T, 2> &AABBs, const Vector<T, 2> &latu) {
    // get cuboid defined by AABBs
    Vector<T, 2> min = AABBs.getMin();
    Vector<T, 2> max = AABBs.getMax();
    // get index of min and max in _GlobalIdx
    Vector<int, 2> idx_min(
        int((min[0] - _Geo.getMin()[0]) / _Geo.getVoxelSize()),
        int((min[1] - _Geo.getMin()[1]) / _Geo.getVoxelSize()));
    for (int i = 0; i < 2; ++i) {
      if (idx_min[i] < 0) idx_min[i] = 0;
    }
    Vector<int, 2> idx_max(
        int((max[0] - _Geo.getMin()[0]) / _Geo.getVoxelSize()),
        int((max[1] - _Geo.getMin()[1]) / _Geo.getVoxelSize()));
    if (idx_max[0] > _Geo.getNx() - 1) idx_max[0] = _Geo.getNx() - 1;
    if (idx_max[1] > _Geo.getNy() - 1) idx_max[1] = _Geo.getNy() - 1;
    // set velocity
    for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
      for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
        int id = Vector<int, 2>(i, j) * _Geo.getIdx();
        if (_Geo.getGlobalIdx()[id] != -1) {
          Voxel<T, 2> &voxel = _Geo.getVoxel(id);
          // check if voxel is inside AABBs
          if (AABBs.isInside(voxel.getLoc()))
            _Velocity[_Geo.getGlobalIdx()[id]] = latu;
        }
      }
    }
  }
  void setVelocity(const Vector<T, 2> &physu, const AABB<T, 2> &AABBs,
                   int flag) {
    setVelocity(AABBs, _Conv.getLatticeU(physu), flag);
  }
  void setVelocity(const AABB<T, 2> &AABBs, const Vector<T, 2> &latu,
                   int flag) {
    // get cuboid defined by AABBs
    Vector<T, 2> min = AABBs.getMin();
    Vector<T, 2> max = AABBs.getMax();
    // get index of min and max in _GlobalIdx
    Vector<int, 2> idx_min(
        int((min[0] - _Geo.getMin()[0]) / _Geo.getVoxelSize()),
        int((min[1] - _Geo.getMin()[1]) / _Geo.getVoxelSize()));
    for (int i = 0; i < 2; ++i) {
      if (idx_min[i] < 0) idx_min[i] = 0;
    }
    Vector<int, 2> idx_max(
        int((max[0] - _Geo.getMin()[0]) / _Geo.getVoxelSize()),
        int((max[1] - _Geo.getMin()[1]) / _Geo.getVoxelSize()));
    if (idx_max[0] > _Geo.getNx() - 1) idx_max[0] = _Geo.getNx() - 1;
    if (idx_max[1] > _Geo.getNy() - 1) idx_max[1] = _Geo.getNy() - 1;
    // set velocity
    for (int j = idx_min[1]; j <= idx_max[1]; ++j) {
      for (int i = idx_min[0]; i <= idx_max[0]; ++i) {
        int id = Vector<int, 2>(i, j) * _Geo.getIdx();
        if (_Geo.getGlobalIdx()[id] != -1) {
          Voxel<T, 2> &voxel = _Geo.getVoxel(id);
          // check if voxel is inside AABBs
          if (AABBs.isInside(voxel.getLoc()) && voxel.getFlag() == flag)
            _Velocity[_Geo.getGlobalIdx()[id]] = latu;
        }
      }
    }
  }
};

// velocity field, receive lattice initial velocity
template <typename T>
struct Velocity2D : public Index2D {
  BaseConverter<T> &Conv;
  Geometry2DLegacy<T> &Geo;
  T **U;
  T U_Init[2];

  // receives physical velocity
  Velocity2D(T physu_init[2], BaseConverter<T> &conv, Geometry2DLegacy<T> &Geo_)
      : Velocity2D(conv, Geo_, conv.getLatticeU(physu_init[0]),
                   conv.getLatticeU(physu_init[1])) {}
  // this constructor receives lattice velocity
  Velocity2D(BaseConverter<T> &conv, Geometry2DLegacy<T> &Geo_, T latu_init[2])
      : Velocity2D(conv, Geo_, latu_init[0], latu_init[1]) {}
  // this constructor receives lattice velocity
  Velocity2D(BaseConverter<T> &conv, Geometry2DLegacy<T> &Geo_, T latu_init0 = 0,
             T latu_init1 = 0)
      : Conv(conv), Geo(Geo_), U_Init{latu_init0, latu_init1} {
    U = new T *[Geo.getN()];
    for (int i = 0; i < Geo.getN(); i++) {
      U[i] = new T[2];
      U[i][0] = T(0);
      U[i][1] = T(0);
    }
    Traverse_Bulk(Geo.getNi(), Geo.getNj(), 1, [this](int id) {
      U[id][0] = U_Init[0];
      U[id][1] = U_Init[1];
    });
    // print log
    std::cout << "[Velocity2D]: "
              << "\n"
              << "using lattice unit"
              << "\n"
              << "U_Init = (" << U_Init[0] << ", " << U_Init[1] << ")"
              << std::endl;
  }

  ~Velocity2D() {
    for (int i = 0; i < Geo.getN(); i++) {
      delete[] U[i];
    }
    delete[] U;
  }

  // Set field
  inline void SetFieldU0(int id, T u0) { U[id][0] = u0; }
  inline void SetFieldU1(int id, T u1) { U[id][1] = u1; }
  inline void SetFieldU(int id, T u0, T u1) {
    U[id][0] = u0;
    U[id][1] = u1;
  }
  void DeActivate(int id) {
    U[id][0] = T(0);
    U[id][1] = T(0);
  }

  void Setby_LatU(std::vector<int> &idx, T LatU[2]) {
    Setby_LatU(idx, LatU[0], LatU[1]);
  }
  void Setby_LatU(std::vector<int> &idx, T LatU0, T LatU1) {
    for (int i = 0; i < idx.size(); i++) {
      SetFieldU(idx[i], LatU0, LatU1);
    }
  }
  void Setby_LatU(int flag, T LatU[2]) { Setby_LatU(flag, LatU[0], LatU[1]); }
  void Setby_LatU(int flag, T LatU0, T LatU1) {
    Traverse_Bulk(Geo.getNi(), Geo.getNj(), 1,
                  [this, flag, LatU0, LatU1](int id) {
                    if (Geo.getVoxel(id).getFlag() == flag) {
                      SetFieldU(id, LatU0, LatU1);
                    }
                  });
  }

  void Setby_PhysU(std::vector<int> &idx, T PhysU[2]) {
    Setby_PhysU(idx, PhysU[0], PhysU[1]);
  }
  void Setby_PhysU(std::vector<int> &idx, T PhysU0, T PhysU1) {
    T LatU0 = Conv.getLatticeU(PhysU0);
    T LatU1 = Conv.getLatticeU(PhysU1);
    Setby_LatU(idx, LatU0, LatU1);
  }
  void Setby_PhysU(int flag, T PhysU[2]) {
    Setby_PhysU(flag, PhysU[0], PhysU[1]);
  }
  void Setby_PhysU(int flag, T PhysU0, T PhysU1) {
    T LatU0 = Conv.getLatticeU(PhysU0);
    T LatU1 = Conv.getLatticeU(PhysU1);
    Setby_LatU(flag, LatU0, LatU1);
  }
};

// CA field
template <typename T>
struct StateField2D : public Index2D {
  Geometry2DLegacy<T> &Geo;
  int *State;
  StateField2D(Geometry2DLegacy<T> &geo) : Geo(geo) {
    State = new int[Geo.getN()];
    Traverse_Peripheral(Geo.getNi(), Geo.getNj(), 0,
                        [this](int id) { State[id] = -1; });

    Traverse_Bulk_OMP(Geo.getNi(), Geo.getNj(), 1, [this](int id) {
      if (Geo.getVoxel(id).getFlag() == -1) {
        State[id] = -1;
      } else {
        State[id] = 0;
      }
    });
  }
  ~StateField2D() { delete[] State; }
  inline void activate(int id) {
    State[id] = 1;
    Geo.getVoxel(id).setFlag(1);
  }
  inline void deactivate(int id) {
    State[id] = -1;
    Geo.getVoxel(id).setFlag(-1);
  }
  virtual T getSolidFraction(int id) = 0;
  virtual T getOrine(int id) = 0;
};

template <typename T>
struct CAGField2D final : public StateField2D<T> {
  using StateField2D<T>::Geo;
  T *Orine;
  T *f;  // solid fraction

  CAGField2D(Geometry2DLegacy<T> &geo_) : StateField2D<T>(geo_) {
    Orine = new T[Geo.getN()];
    f = new T[Geo.getN()];
    Index2D::Traverse_Peripheral(Geo.getNi(), Geo.getNj(), 0, [this](int id) {
      Orine[id] = T(-1);
      f[id] = T(-1);
    });

    Index2D::Traverse_Bulk_OMP(Geo.getNi(), Geo.getNj(), 1, [this](int id) {
      if (Geo.getVoxel(id).getFlag() == -1) {
        Orine[id] = T(-1);
        f[id] = T(-1);
      } else {
        Orine[id] = T(0);
        f[id] = T(0);
      }
    });
  }
  ~CAGField2D() {
    delete[] Orine;
    delete[] f;
  }
  inline void Activate(int id) { StateField2D<T>::activate(id); }
  // set state and geoFlag to -1
  inline void Deactivate(int id) { StateField2D<T>::deactivate(id); }
  T getSolidFraction(int id) override { return f[id]; }
  T getOrine(int id) override { return Orine[id]; }
};

template <typename T>
struct CAZSField2D final : public StateField2D<T> {
  using StateField2D<T>::Geo;
  T *f;  // solid fraction

#ifdef _FLB_DEBUG
  // Ceq - Cl
  T *Ceq;
  T *Ceq_Cl;
  T *aniso;
  // curvature of Solid-Liquid interface
  T *K;
#endif

  CAZSField2D(Geometry2DLegacy<T> &geo) : StateField2D<T>(geo) {
    f = new T[Geo.getN()];
    Index2D::Traverse_Peripheral(Geo.getNi(), Geo.getNj(), 0,
                                 [this](int id) { f[id] = T(-1); });

    Index2D::Traverse_Bulk_OMP(Geo.getNi(), Geo.getNj(), 1, [this](int id) {
      if (Geo.getVoxel(id).getFlag() == -1) {
        f[id] = T(1);
      } else {
        f[id] = T(0);
      }
    });
#ifdef _FLB_DEBUG
    Ceq_Cl = new T[Geo.getN()];
    std::fill_n(Ceq_Cl, Geo.getN(), T(0));
    Ceq = new T[Geo.getN()];
    std::fill_n(Ceq, Geo.getN(), T(0));
    aniso = new T[Geo.getN()];
    std::fill_n(aniso, Geo.getN(), T(0));
    K = new T[Geo.getN()];
    std::fill_n(K, Geo.getN(), T(0));
#endif
  }
  ~CAZSField2D() {
    delete[] f;
#ifdef _FLB_DEBUG
    delete Ceq_Cl;
    delete Ceq;
    delete aniso;
    delete[] K;
#endif
  }
  inline void Activate(int id) { StateField2D<T>::activate(id); }
  // set state and geoFlag to -1
  inline void Deactivate(int id) { StateField2D<T>::deactivate(id); }
  T getSolidFraction(int id) override { return f[id]; }
  T getOrine(int id) override { return 0; }
};
