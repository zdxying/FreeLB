// legacy_equilibrium.h

#pragma once

#include "legacy/legacy_lbm/legacy_populations.h"
#include "utils/util.h"

template <typename T, typename LatSet>
struct LegacyEquilibrium {

  static inline T Order1(int k, const T *u, T rho) {
    return LatSet::w[k] * rho *
           (T(1) + LatSet::InvCs2 * Vect2D<T>::dot(u, LatSet::c[k]));
    // return LatSet::w[k] * pop.rho * (T(1) + LatSet::InvCs2 * uc);
    // openlb:
    // return LatSet::w[k] * pop.rho * (1 + LatSet::InvCs2 * uc) - LatSet::w[k];
    // = LatSet::w[k] * (pop.rho + LatSet::InvCs2 * uc * pop.rho) - LatSet::w[k]
    // = LatSet::w[k] * (pop.rho + LatSet::InvCs2 * uc * pop.rho - 1)
  }

  static inline T Order1_Incompresible(int k, const T *u, T rho) {
    return LatSet::w[k] *
           (rho + LatSet::InvCs2 * Vect2D<T>::dot(u, LatSet::c[k]));
    // return LatSet::w[k] * (pop.rho + LatSet::InvCs2 * uc);
  }

  static inline T Order2(int k, const T *u, T rho, T u2) {
    // return LatSet::w[k] * rho *
    //        (T(1) + LatSet::InvCs2 * Vect2D<T>::dot(u, LatSet::c[k]) +
    //         LatSet::InvCs2 * LatSet::InvCs2 * Vect2D<T>::dot(u, LatSet::c[k])
    //         *
    //             Vect2D<T>::dot(u, LatSet::c[k]) * T(0.5) -
    //         LatSet::InvCs2 * u2 * T(0.5));
    return LatSet::w[k] * rho *
           (T(1) + LatSet::InvCs2 * Vect2D<T>::dot(u, LatSet::c[k]) +
            pow(LatSet::InvCs2 * Vect2D<T>::dot(u, LatSet::c[k]), 2) * T(0.5) -
            LatSet::InvCs2 * u2 * T(0.5));
    // return LatSet::w[k] * pop.rho * (1 + LatSet::InvCs2 * uc
    // + LatSet::InvCs2 * LatSet::InvCs2 * uc * uc * T(0.5)
    // - LatSet::InvCs2 * u2 * T(0.5));
  }

  static inline void Feq_firstOrder(T *feq, const T *u, T rho) {
    // #pragma omp simd
    for (int k = 0; k < LatSet::q; k++) {
      feq[k] = Order1(k, u, rho);
    }
  }

  static inline void Feq_firstOrder_Incompresible(T *feq, const T *u, T rho) {
    // #pragma omp simd
    for (int k = 0; k < LatSet::q; k++) {
      feq[k] = Order1_Incompresible(k, u, rho);
    }
  }

  static inline void Feq_secondOrder(T *feq, const T *u, T rho) {
    T u2 = Vect2D<T>::sqr(u);
    // #pragma omp simd
    for (int k = 0; k < LatSet::q; k++) {
      feq[k] = Order2(k, u, rho, u2);
    }
  }
};

template <typename T, typename LatSet>
struct Equilibrium3D {
  T _OMEGA;
  T _1_OMEGA;
  // constructor
  Equilibrium3D(T OMEGA) : _OMEGA(OMEGA), _1_OMEGA(T(1) - OMEGA) {}
  // METHOD
  static inline T Order1(int k, const AbstractPopulation<T, LatSet> &pop) {
    return LatSet::w[k] * pop.getRho() *
           (T(1) + LatSet::InvCs2 * (pop.getVelocity() * LatSet::c[k]));
  }
  static inline T Order1_Incompresible(
      int k, const AbstractPopulation<T, LatSet> &pop) {
    return LatSet::w[k] *
           (pop.getRho() + LatSet::InvCs2 * (pop.getVelocity() * LatSet::c[k]));
  }
  static inline T Order2(int k, const AbstractPopulation<T, LatSet> &pop) {
    const T uc = pop.getVelocity() * LatSet::c[k];
    return LatSet::w[k] * pop.getRho() *
           (T(1) + LatSet::InvCs2 * (uc) +
            pow(uc, 2) * T(0.5) * LatSet::InvCs4 -
            LatSet::InvCs2 * pop.getVelocity().getnorm2() * T(0.5));
  }
  static inline T Order2(int k, const AbstractPopulation<T, LatSet> &pop,
                         T u2) {
    const T uc = pop.getVelocity() * LatSet::c[k];
    return LatSet::w[k] * pop.getRho() *
           (T(1) + LatSet::InvCs2 * (uc) +
            pow(uc, 2) * T(0.5) * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T(0.5));
  }
  // bgk: pop.fpostcol[k] = _OMEGA * feq[k] + _1_OMEGA * pop.f[k];
  void Feq_firstOrder(AbstractPopulation<T, LatSet> &pop) {
    // #pragma omp simd
    for (int k = 0; k < LatSet::q; ++k) {
      pop.getDDFpostcol(k) = _OMEGA * Order1(k, pop) + _1_OMEGA * pop.getDDF(k);
    }
  }
  void Feq_firstOrder_Incompresible(AbstractPopulation<T, LatSet> &pop) {
    // #pragma omp simd
    for (int k = 0; k < LatSet::q; ++k) {
      pop.getDDFpostcol(k) =
          _OMEGA * Order1_Incompresible(k, pop) + _1_OMEGA * pop.getDDF(k);
    }
  }
  void Feq_secondOrder(AbstractPopulation<T, LatSet> &pop) {
    T u2 = pop.getVelocity().getnorm2();
    // #pragma omp simd
    for (int k = 0; k < LatSet::q; ++k) {
      pop.getDDFpostcol(k) =
          _OMEGA * Order2(k, pop, u2) + _1_OMEGA * pop.getDDF(k);
    }
  }
};

// void NEE(population<T, LatSet> *LB, int id1, int id2, int k)
// {
//     // T feq1 = feqO2(k, id, &LB1, Lat.Field[id].dotU(c[k]),
//     Lat.Field[id].ScalarU2());
//     // T feq2 = feqO2(k, id2, &LB2, Lat.Field[id2].dotU(c[k]),
//     Lat.Field[id2].ScalarU2()); T feq1 = feqO1_Incompresible(k, LB[id1],
//     Lat.Field[id1].dotU(c[k])); T feq2 = feqO1_Incompresible(k, LB[id2],
//     Lat.Field[id2].dotU(c[k]));
//     // LB1.f[k] = LB2.f[k] + feq1 - feq2;
//     LB[id1].f[k] = LB[id2].f[k] + feq1 - feq2;
// }