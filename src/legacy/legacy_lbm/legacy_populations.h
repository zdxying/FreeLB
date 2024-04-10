#pragma once
// populations:
// density distribution functions in one single cell
#include <array>
#include <string>
#include <vector>

#include "data_struct/Vector.h"

template <typename T>
class Basiclbm {
 public:
  // get
  virtual T getPoprho(int Id) const = 0;
  virtual T getPhysrho(int Id) const = 0;
  virtual T getlatrhoinit() const = 0;
  virtual std::string getname() const = 0;
  virtual std::string getnamerho() const = 0;
  virtual T getOmega() const = 0;
  virtual T getPopSize() const = 0;
  virtual T getBuoyancy(int id) const = 0;
  virtual std::vector<int> *getIdx() = 0;
  virtual std::vector<int> *getInIdx() = 0;
  virtual std::vector<int> *getPostRhoCells() = 0;
  virtual std::vector<int> *getPostUCells() = 0;
  // virtual Force2D<T> *getForce() = 0;
  // set
  virtual void setPopRho(int Id, T rho) = 0;
  // virtual void setPopRho(int ni, int nj, T rho) = 0;
  // virtual void setpoprho_From_Vector(std::vector<int> &popidx, T rho_) = 0;
  // virtual void setPopRho_From_VoxelFlag(int flag, T rho) = 0;
  // add
  // virtual void addtoBuoyancy(std::vector<int> &popidx, Force2D<T> *Force) =
  // 0; virtual void resetBuoyancy() = 0;
  virtual void addtoRho(int id, T value) = 0;
  virtual void addtoPop(int id, T value) = 0;
  // virtual void addtoPoppostcol(int id, T value);// compute
  // virtual void compute_rho(int id) = 0;
};

template <typename T, typename LatSet>
struct PopulationMethod {
  static inline void Compute_rho(const T *f, T &rho) {
    rho = T(0);
    for (int i = 0; i < LatSet::q; i++) {
      rho += f[i];
    }
  }
  static inline T Get_U(const T *f, const T rho, const int k) {  // k: 0, 1, 2
    T rhoU = T(0);
    for (int i = 0; i < LatSet::q; i++) {
      rhoU += f[i] * LatSet::c[i][k];
    }
    return rhoU / rho;
  }
  static inline void Compute_U(const T *f, const T rho, T *U) {
    for (int i = 0; i < LatSet::d; i++) {
      U[i] = Get_U(f, rho, i);
    }
  }
};
template <typename T, typename LatSet>
struct BasicPopulation {
  T f[LatSet::q];
  T rho;
  // Pointers to neighboring cells
  std::vector<BasicPopulation<T, LatSet> *> nbrs;
  BasicPopulation() : rho(0) {
    nbrs.resize(LatSet::q, nullptr);
    for (int i = 0; i < LatSet::q; i++) {
      f[i] = 0;
    }
  }
  BasicPopulation(T rho_) : rho(rho_) {
    nbrs.resize(LatSet::q, nullptr);
    for (int i = 0; i < LatSet::q; i++) {
      f[i] = rho * LatSet::w[i];
    }
  }
  void Compute_rho() { PopulationMethod<T, LatSet>::Compute_rho(f, rho); }
  T Get_U(int k) { return PopulationMethod<T, LatSet>::Get_U(f, rho, k); }
  void Compute_U(T *U) { PopulationMethod<T, LatSet>::Compute_U(f, rho, U); }
};

// distribution function with a post-collision distribution function (buffer)
// fpostcol[] is used to avoid race condition when streaming
// streaming sstep: fpostcol[] -> f[]
template <typename T, typename LatSet>
struct population {
  T f[LatSet::q];         // density distribution function
  T fpostcol[LatSet::q];  // f after collision
  T rho;
  // constructor
  population(T rho_ = T(0)) : rho(rho_) {
    for (int i = 0; i < LatSet::q; i++) {
      f[i] = rho * LatSet::w[i];
      fpostcol[i] = f[i];
    }
  }
  inline void Compute_rho() {
    PopulationMethod<T, LatSet>::Compute_rho(f, rho);
  }
  inline T Get_U(const int k) {
    return PopulationMethod<T, LatSet>::Get_U(f, rho, k);
  }
  inline void Compute_U(T *U) {
    PopulationMethod<T, LatSet>::Compute_U(f, rho, U);
  }
  inline void addtoRho(const T value) { rho += value; }
  inline void addtoPop(const T value) {
    for (int i = 0; i < LatSet::q; i++) {
      f[i] += value * LatSet::w[i];
    }
  }

  inline void Compute_Boundaryrho(const std::vector<int> &outflow,
                                  const int *opp) {
    rho = T(0);
    Compute_rho();
    for (int i = 0; i < outflow.size(); i++) {
      rho += f[opp[outflow[i]]] - f[outflow[i]];
    }
  }
};

template <typename T, typename LatSet>
struct AbstractPopulation {
  virtual void Compute_rho() = 0;
  virtual void Compute_U() = 0;
  virtual void addtoRho(const T value) {}
  virtual void addtoPop(const T value) {}
  // get density distribution function
  virtual inline T &getDDF(int k) = 0;
  virtual inline const T &getDDF(int k) const = 0;
  // get density distribution function - post collision
  virtual inline T &getDDFpostcol(int k) {
    static T dummy;
    return dummy;
  }
  virtual inline const T &getDDFpostcol(int k) const {
    static T dummy;
    return dummy;
  }
  virtual inline T getRho() const = 0;
  virtual inline AbstractPopulation<T, LatSet> *getNeighbor(int k) {
    static AbstractPopulation<T, LatSet> *dummy;
    return dummy;
  }
  virtual inline Voxel<T, LatSet::d> &getVoxel() {
    static Voxel<T, LatSet::d> dummy;
    return dummy;
  }
  virtual inline const Voxel<T, LatSet::d> &getVoxel() const {
    static Voxel<T, LatSet::d> dummy;
    return dummy;
  }
  virtual inline Vector<T, LatSet::d> &getVelocity() {
    static Vector<T, LatSet::d> dummy;
    return dummy;
  }
  virtual inline const Vector<T, LatSet::d> &getVelocity() const {
    static Vector<T, LatSet::d> dummy;
    return dummy;
  }
};

// population with neighbor information
// the most memory consuming population: including postcol and nbrs
template <typename T, typename LatSet>
struct PopulationNbr final : public AbstractPopulation<T, LatSet> {
  // T f[LatSet::q];
  // T fpostcol[LatSet::q];
  T rho;
  std::array<T, LatSet::q> f;
  std::array<T, LatSet::q> fpostcol;
  // Pointer to neighboring cells
  std::array<PopulationNbr<T, LatSet> *, LatSet::q> nbrs;
  // pointer to voxel
  Voxel<T, LatSet::d> &voxel;
  // pointer to velocity field
  Vector<T, LatSet::d> &velocity;

  // constructor
  PopulationNbr(Voxel<T, LatSet::d> &voxel_, Vector<T, LatSet::d> &velocity_,
                T rho_ = T(0))
      : voxel(voxel_), velocity(velocity_), rho(rho_) {
    for (int i = 0; i < LatSet::q; ++i) {
      nbrs[i] = nullptr;
      f[i] = rho * LatSet::w[i];
      fpostcol[i] = f[i];
    }
  }
  void Initialize(T rho_) {
    rho = rho_;
    for (int i = 0; i < LatSet::q; ++i) {
      f[i] = rho * LatSet::w[i];
      fpostcol[i] = f[i];
    }
  }
  inline T &getDDF(int k) override { return f[k]; }
  inline const T &getDDF(int k) const override { return f[k]; }
  inline T &getDDFpostcol(int k) override { return fpostcol[k]; }
  inline const T &getDDFpostcol(int k) const override { return fpostcol[k]; }
  inline T getRho() const override { return rho; }
  // get pointer to neighbor population
  inline PopulationNbr<T, LatSet> *getNeighbor(int i) override {
    return nbrs[i];
  }
  inline Voxel<T, LatSet::d> &getVoxel() override { return voxel; }
  inline const Voxel<T, LatSet::d> &getVoxel() const override { return voxel; }
  inline Vector<T, LatSet::d> &getVelocity() override { return velocity; }
  inline const Vector<T, LatSet::d> &getVelocity() const override {
    return velocity;
  }
  T getVelocityx() const { return velocity[0]; }
  T getVelocityy() const { return velocity[1]; }
  T getVelocityz() const { return velocity[2]; }
  // set neighbor
  void setNeighbor(int i, PopulationNbr<T, LatSet> *nbr) { nbrs[i] = nbr; }
  void resetNeighbor(int i) { nbrs[i] = nullptr; }
  void resetNeighbor() {
    for (int i = 0; i < LatSet::q; ++i) nbrs[i] = nullptr;
  }

  void Compute_rho() override {
    rho = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho += f[i];
  }
  void Compute_rho(T &rho_) {
    rho_ = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho_ += f[i];
  }
  void Compute_U() override {
    for (int k = 0; k < LatSet::d; ++k) {
      T rhoU = T(0);
      for (int i = 0; i < LatSet::q; ++i) {
        rhoU += f[i] * LatSet::c[i][k];
      }
      velocity[k] = rhoU / rho;
    }
  }
  void Compute_U(Vector<T, LatSet::d> &velocity_) {
    for (int k = 0; k < LatSet::d; ++k) {
      T rhoU = T(0);
      for (int i = 0; i < LatSet::q; ++i) {
        rhoU += f[i] * LatSet::c[i][k];
      }
      velocity_[k] = rhoU / rho;
    }
  }
  void Compute_rhoU() {
    rho = T(0);
    for (int j = 0; j < LatSet::q; ++j) rho += f[j];
    for (int k = 0; k < LatSet::d; ++k) {
      T rhoU = T(0);
      for (int i = 0; i < LatSet::q; ++i) {
        rhoU += f[i] * LatSet::c[i][k];
      }
      velocity[k] = rhoU / rho;
    }
  }
  void Compute_rhoU(T &rho_, Vector<T, LatSet::d> &velocity_) {
    rho_ = T(0);
    for (int i = 0; i < LatSet::q; ++i) rho_ += f[i];
    for (int k = 0; k < LatSet::d; ++k) {
      T rhoU = T(0);
      for (int i = 0; i < LatSet::q; ++i) {
        rhoU += f[i] * LatSet::c[i][k];
      }
      velocity_[k] = rhoU / rho_;
    }
  }
  inline void BGK_O1(const T omega, const T _omega) {
    for (int i = 0; i < LatSet::q; ++i) {
      fpostcol[i] =
          omega * (LatSet::w[i] * rho *
                   (T(1) + LatSet::InvCs2 * (velocity * LatSet::c[i]))) +
          _omega * f[i];
    }
  }
  // incompressible
  inline void BGK_O1incomp(const T omega, const T _omega) {
    for (int i = 0; i < LatSet::q; ++i) {
      fpostcol[i] =
          omega * (LatSet::w[i] *
                   (rho + LatSet::InvCs2 * (velocity * LatSet::c[i]))) +
          _omega * f[i];
    }
  }
  inline void BGK_O2(const T omega, const T _omega) {
    // get u2
    T u2 = velocity.getnorm2();
    for (int i = 0; i < LatSet::q; ++i) {
      T uc = velocity * LatSet::c[i];
      fpostcol[i] = omega * (LatSet::w[i] * rho *
                             (T(1) + LatSet::InvCs2 * uc +
                              uc * uc * T(0.5) * LatSet::InvCs4 -
                              LatSet::InvCs2 * u2 * T(0.5))) +
                    _omega * f[i];
    }
  }
  void addForceO2(const Vector<T, LatSet::d> &Force) {
    for (int k = 0; k < LatSet::q; ++k) {
      // get c * u * LatSet::InvCs4
      T cu = LatSet::c[k] * velocity * LatSet::InvCs4;
      // (c - u) * LatSet::InvCs2
      Vector<T, LatSet::d> c_u = (LatSet::c[k] - velocity) * LatSet::InvCs2;
      fpostcol[k] += Force * (c_u + cu * LatSet::c[k]) * LatSet::w[k];
    }
  }
  void addForceO2(const Vector<T, LatSet::d> &Force, const T omega) {
    for (int k = 0; k < LatSet::q; ++k) {
      // get c * u * LatSet::InvCs4
      T cu = LatSet::c[k] * velocity * LatSet::InvCs4;
      // (c - u) * LatSet::InvCs2
      Vector<T, LatSet::d> c_u = (LatSet::c[k] - velocity) * LatSet::InvCs2;
      T force = Force * (c_u + cu * LatSet::c[k]);
      force *= LatSet::w[k] * (T(1) - omega * T(0.5));
      fpostcol[k] += force;
    }
  }
  void addForceO1(const Vector<T, LatSet::d> &Force, const T omega) {
    const T cons = (T(1) - omega * T(0.5)) * LatSet::InvCs2;
    for (int k = 0; k < LatSet::q; ++k) {
      fpostcol[k] += Force * LatSet::c[k] * LatSet::w[k] * cons;
    }
  }
  void addForceO1(const Vector<T, LatSet::d> &Force) {
    for (int k = 0; k < LatSet::q; ++k) {
      fpostcol[k] += Force * LatSet::c[k] * LatSet::w[k] * LatSet::InvCs2;
    }
  }
};
//////////////////////////////////////////////////
// distribution function using a memory saving scheme
// only one copy of density distribution function is used
// known as AA-Pattern, by Bailey, et al. in 2009
template <typename T, typename LatSet>
struct population_AA : public BasicPopulation<T, LatSet> {
  population_AA() : BasicPopulation<T, LatSet>() {}
  population_AA(T rho_) : BasicPopulation<T, LatSet>(rho_) {}
  // TODO: may be optimized by using smaller ftemp[]
  // e.g., ftemp[LatSet::q/2]
  // by doing so, the lattice set need to be rearranged
  void exchange() {
    T ftemp[LatSet::q];
    for (int i = 0; i < LatSet::q; i++) {
      ftemp[i] = this->f[i];
    }
    for (int i = 0; i < LatSet::q; i++) {
      this->f[i] = ftemp[LatSet::opp[i]];
    }
  }
};
