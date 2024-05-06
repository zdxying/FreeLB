// lattice
#pragma once

#include <vector>

#include "geometry/geometry2d.h"
#include "lbm/unit_converter.h"
#include "utils/util.h"

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
