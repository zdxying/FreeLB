// force.h

#include "utils/util.h"

template <typename T>
struct Force2D {
  T F[2];
  Force2D() : F{T(0), T(0)} {}
  void Reset() {
    F[0] = T(0);
    F[1] = T(0);
  }
};
template <typename T, typename LatSet>
struct ForceMethod2D {
  static T VelocityO1_SpatialO1(
      int k, T *f, T *u, T omega)  // velocity of order 1, spatial order 1
  {
    return LatSet::w[k] * Vect2D<T>::dot(f, LatSet::c[k]) * LatSet::InvCs2;
  }
  static T VelocityO1_SpatialO2(
      int k, T *f, T *u, T omega)  // velocity of order 1, spatial order 2
  {
    return LatSet::w[k] * Vect2D<T>::dot(f, LatSet::c[k]) * LatSet::InvCs2 *
           (T(1) - omega * T(0.5));
  }
  static T VelocityO2_SpatialO1(
      int k, T *f, T *u, T omega)  // velocity of order 2, spatial order 1
  {
    return LatSet::w[k] *
           (f[0] * (LatSet::InvCs2 * (LatSet::c[k][0] - u[0]) +
                    LatSet::InvCs4 * Vect2D<T>::dot(u, LatSet::c[k]) *
                        LatSet::c[k][0]) +
            f[1] * (LatSet::InvCs2 * (LatSet::c[k][1] - u[1]) +
                    LatSet::InvCs4 * Vect2D<T>::dot(u, LatSet::c[k]) *
                        LatSet::c[k][1]));
  }
  static T VelocityO2_SpatialO2(
      int k, T *f, T *u, T omega)  // velocity of order 2, spatial order 2
  {
    return LatSet::w[k] *
           (f[0] * (LatSet::InvCs2 * (LatSet::c[k][0] - u[0]) +
                    LatSet::InvCs4 * Vect2D<T>::dot(u, LatSet::c[k]) *
                        LatSet::c[k][0]) +
            f[1] * (LatSet::InvCs2 * (LatSet::c[k][1] - u[1]) +
                    LatSet::InvCs4 * Vect2D<T>::dot(u, LatSet::c[k]) *
                        LatSet::c[k][1])) *
           (T(1) - omega * T(0.5));
  }
};

template <typename T, typename LatSet, typename lbmType>
class ForceManager {
 private:
  T _OMEGA;
  const int dim = LatSet::d - 1;
  std::vector<Vector<T, LatSet::d>> _Force;
  lbmType &_LbmNS;
  // source
  std::vector<Basiclbm<T> *> _LbmSource;

 public:
  ForceManager(lbmType &LbmNS,
               const Vector<T, LatSet::d> &force = Vector<T, LatSet::d>{})
      : _LbmNS(LbmNS), _OMEGA(LbmNS.getOmega()) {
    _Force.resize(_LbmNS.getPopSize(), force);
  }
  template <typename... Args>
  void AddSource(Basiclbm<T> *Lbm, Args*... args) {
    _LbmSource.push_back(Lbm);
    AddSource(args...);
  }
  void AddSource(Basiclbm<T> *Lbm) { _LbmSource.push_back(Lbm); }
  // get
  const std::vector<Vector<T, LatSet::d>> &getForce() const { return _Force; }
  std::vector<Vector<T, LatSet::d>> &getForce() { return _Force; }

  // call virtual function may cause performance issue
  void AddtoBuoyancy() {
    for (Basiclbm<T> *lbmSource : _LbmSource) {
      const std::vector<int> *idx = lbmSource->getIdx();
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
      for (int id : *idx) {
        _Force[id][dim] += lbmSource->getBuoyancy(id);
      }
    }
  }
  // call virtual function may cause performance issue
  void AddForcetoPop() {
    const std::vector<int> *idx = _LbmNS.getIdx();
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id : *idx) {
      // calculate force
      _LbmNS.addtoPoppostcol(id, _Force[id][dim]);
    }
  }
  void ResetForce() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (Vector<T, LatSet::d> &force : _Force) {
      force.clear();
    }
  }
  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T),
            bool InnerRho = true, bool InnerU = true>
  void applyBuoyancy() {
    _LbmNS.template Collide_BGK<get_feq>();
    _LbmNS.Stream();
    _LbmNS.BCs();
    AddtoBuoyancy();
    _LbmNS.addForce(getForce());
    _LbmNS.template PostStreamRho<InnerRho>();
    _LbmNS.template PostStreamU<InnerU>();
    ResetForce();
  }
  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T),
            bool InnerRho = true, bool InnerU = true>
  void applyConstBulkForce() {
    _LbmNS.template Collide_BGK<get_feq>();
    _LbmNS.Stream();
    _LbmNS.BCs();
    _LbmNS.addForce(getForce());
    _LbmNS.template PostStreamRho<InnerRho>();
    _LbmNS.template PostStreamU<InnerU>();
  }
};
