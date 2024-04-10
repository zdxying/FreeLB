// bcs 2d
#pragma once

#include "legacy/legacy_boundary/legacy_basic_boundary.h"
#include "legacy/legacy_lattice/legacyfield3D.h"
#include "legacy/legacy_lbm/legacy_populations.h"

template <typename T, typename LatSet>
class BoundaryCell3D_ {
 private:
  // outflow directions
  std::vector<int> outflows;
  // inflow directions
  std::vector<int> inflows;
  // pointer to the population
  AbstractPopulation<T, LatSet> &BdPop;

 public:
  BoundaryCell3D_(AbstractPopulation<T, LatSet> &pop) : BdPop(pop) { set(); }
  // get
  AbstractPopulation<T, LatSet> &getBdPop() { return BdPop; }
  std::vector<int> &getoutflows() { return outflows; }
  const std::vector<int> &getoutflows() const { return outflows; }
  std::vector<int> &getinflows() { return inflows; }
  const std::vector<int> &getinflows() const { return inflows; }
  const int getoutflow(int i) const { return outflows[i]; }
  const int getinflow(int i) const { return inflows[i]; }
  // set
  void set() {
    inflows.clear();
    outflows.clear();
    const std::vector<Voxel<T, LatSet::d> *> &voxnbrs =
        BdPop.getVoxel().getNeighborList();
    for (int k = 1; k < voxnbrs.size(); ++k) {
      // inflow
      if (voxnbrs[k] != nullptr) {
        if (voxnbrs[k]->getFlag() != -1) inflows.push_back(LatSet::opp[k]);
      }
      // outflow
      if (voxnbrs[k] == nullptr && voxnbrs[LatSet::opp[k]] != nullptr) {
        if (voxnbrs[LatSet::opp[k]]->getFlag() != -1)
          outflows.push_back(LatSet::opp[k]);
      } else if (voxnbrs[k] != nullptr) {
        if (voxnbrs[k]->getFlag() == -1 && voxnbrs[LatSet::opp[k]] != nullptr) {
          if (voxnbrs[LatSet::opp[k]]->getFlag() != -1)
            outflows.push_back(LatSet::opp[k]);
        }
      }
    }
  }
};
template <typename T, typename LatSet>
class BoundaryCell3D {
 private:
  // outflow directions
  std::vector<int> outflows;
  // inflow directions
  std::vector<int> inflows;
  // pointer to the population
  PopulationNbr<T, LatSet> *BdPop;
  // Population id
  int Id;

 public:
  BoundaryCell3D(PopulationNbr<T, LatSet> *pop)
      : BdPop(pop), Id(pop->getVoxel().getId()) {
    set();
  }
  BoundaryCell3D(Voxel<T, LatSet::d> &vox, std::vector<int> &State)
      : Id(vox.getId()), BdPop(nullptr) {
    set(vox, State);
  }
  // get
  PopulationNbr<T, LatSet> *getBdPop() { return BdPop; }
  std::vector<int> &getoutflows() { return outflows; }
  const std::vector<int> &getoutflows() const { return outflows; }
  std::vector<int> &getinflows() { return inflows; }
  const std::vector<int> &getinflows() const { return inflows; }
  int getoutflow(int i) const { return outflows[i]; }
  int getinflow(int i) const { return inflows[i]; }
  int getId() const { return Id; }
  // set
  void set() {
    inflows.clear();
    outflows.clear();
    const std::vector<Voxel<T, LatSet::d> *> &voxnbrs =
        BdPop->getVoxel().getNeighborList();
    for (int k = 1; k < voxnbrs.size(); ++k) {
      // inflow
      if (voxnbrs[k] != nullptr) {
        if (voxnbrs[k]->getFlag() != -1) inflows.push_back(LatSet::opp[k]);
      }
      // outflow
      if (voxnbrs[k] == nullptr && voxnbrs[LatSet::opp[k]] != nullptr) {
        if (voxnbrs[LatSet::opp[k]]->getFlag() != -1)
          outflows.push_back(LatSet::opp[k]);
      } else if (voxnbrs[k] != nullptr) {
        if (voxnbrs[k]->getFlag() == -1 && voxnbrs[LatSet::opp[k]] != nullptr) {
          if (voxnbrs[LatSet::opp[k]]->getFlag() != -1)
            outflows.push_back(LatSet::opp[k]);
        }
      }
    }
  }
  void set(Voxel<T, LatSet::d> &vox, std::vector<int> &State) {
    inflows.clear();
    outflows.clear();
    const std::vector<Voxel<T, LatSet::d> *> &voxnbrs = vox.getNeighborList();
    for (int k = 1; k < LatSet::q; ++k) {
      // TODO: may be improved
      if (voxnbrs[k] == nullptr) {
        std::cout << "Reached Boundary of the Computational Domain"
                  << std::endl;
        exit(-2);
      }
      int id = voxnbrs[k]->getId();
      // inflow
      if (State[id] != -1) inflows.push_back(LatSet::opp[k]);
      // outflow
      if (State[id] == -1 && State[voxnbrs[LatSet::opp[k]]->getId()] != -1)
        outflows.push_back(LatSet::opp[k]);
    }
  }
};

template <typename T, typename LatSet>
class BasicBoundary3D {
 protected:
  std::vector<BoundaryCell3D<T, LatSet>> BdPops;
  // geometry flag
  int _flag;
  // name for output info
  std::string _name;
  // boundary values
  // std::vector<T> BoundaryValues;

 public:
  BasicBoundary3D(std::string name, int flag = 1) : _flag(flag), _name(name) {}

  const int &getFlag() const { return _flag; }
  std::vector<BoundaryCell3D<T, LatSet>> &getBdPops() { return BdPops; }
  const std::vector<BoundaryCell3D<T, LatSet>> &getBdPops() const {
    return BdPops;
  }
  void clear() { BdPops.clear(); }
  void addtoBd(BoundaryCell3D<T, LatSet> &bdpop) { BdPops.emplace_back(bdpop); }
  void addtoBd(Voxel<T, LatSet::d> &vox, std::vector<int> &State) {
    BdPops.emplace_back(BoundaryCell3D<T, LatSet>(vox, State));
  }
  // this will be called in lbm3D
  void Setup(std::vector<PopulationNbr<T, LatSet>> &pops) {
    for (PopulationNbr<T, LatSet> &pop : pops) {
      if (pop.getVoxel().getFlag() == _flag)
        BdPops.emplace_back(BoundaryCell3D<T, LatSet>(&pop));
    }
  }
  // booundary cell streaming
  void Stream() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (auto &flowdir : BdPops) {
      auto *pop = flowdir.getBdPop();
      pop->f[0] = pop->fpostcol[0];  // IMPROTANT!!!
      for (auto dir : flowdir.getinflows()) {
        pop->f[dir] = pop->getNeighbor(dir)->fpostcol[dir];
      }
    }
  }
  void Stream(std::vector<PopulationNbr<T, LatSet>> &pops) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (auto &flowdir : BdPops) {
      auto &pop = pops[flowdir.getId()];
      pop.f[0] = pop.fpostcol[0];  // IMPROTANT!!!
      for (auto dir : flowdir.getinflows()) {
        pop.f[dir] = pop.getNeighbor(dir)->fpostcol[dir];
      }
    }
  }

  virtual void Apply() = 0;
  virtual void Apply(std::vector<PopulationNbr<T, LatSet>> &pops) = 0;
  void getinfo() {
    std::cout << std::setw(18) << std::left << _name << std::setw(10)
              << std::left << BdPops.size() << std::endl;
  }
};

template <typename T, typename LatSet>
struct BBlikemethod3D_ {
  // pop.f[k] = pop.fpostcol[LatSet::opp[k]];
  static inline void normal_bounceback(AbstractPopulation<T, LatSet> &pop,
                                       int k) {
    pop.getDDF(k) = pop.getDDFpostcol(LatSet::opp[k]);
  }
  // pop.f[k] = 2 * rho * LatSet::w[k] - pop.fpostcol[LatSet::opp[k]];
  static inline void anti_bounceback_simplified(
      AbstractPopulation<T, LatSet> &pop, int k) {
    pop.getDDF(k) =
        2 * pop.getRho() * LatSet::w[k] - pop.getDDFpostcol(LatSet::opp[k]);
  }
  // pop.f[k] = pop.fpostcol[LatSet::opp[k]] + 2 * LatSet::w[k] * pop.rho * uc *
  // LatSet::InvCs2; call: (pop[Id], Dir, Field.U[Id])
  static inline void movingwall_bounceback(AbstractPopulation<T, LatSet> &pop,
                                           int k) {
    pop.getDDF(k) = pop.getDDFpostcol(LatSet::opp[k]) +
                    2 * LatSet::w[k] * pop.getRho() *
                        (pop.getVelocity() * LatSet::c[k]) * LatSet::InvCs2;
  }
  // pop.f[k] = 2 * feq - pop.fpostcol[LatSet::opp[k]]; first order
  // call: (pop[Id], Dir, Field.U[Id])
  static inline void anti_bounceback_O1(AbstractPopulation<T, LatSet> &pop,
                                        int k) {
    pop.getDDF(k) =
        2 * Equilibrium<T, LatSet>::Order1(k, pop.getVelocity(), pop.getRho()) -
        pop.getDDFpostcol(LatSet::opp[k]);
  }
  // pop.f[k] = 2 * feq - pop.fpostcol[LatSet::opp[k]]; second order
  // call: (pop[Id], Dir, Field.U[Id])
  static inline void anti_bounceback_O2(AbstractPopulation<T, LatSet> &pop,
                                        int k) {
    pop.getDDF(k) =
        2 * Equilibrium<T, LatSet>::Order2(k, pop.getVelocity(), pop.getRho(),
                                           pop.getVelocity().getnorm2()) -
        pop.getDDFpostcol(LatSet::opp[k]);
  }
  // pressure boundary condition: anti bounceback scheme
  // pop.f[k] = 2 * LatSet::w[k] * pop.rho * (1 + LatSet::InvCs2 *
  // LatSet::InvCs2 * uc * uc * T(0.5) - LatSet::InvCs2 * u2 * T(0.5))
  static inline void anti_bounceback_pressure(
      AbstractPopulation<T, LatSet> &pop, int k) {
    // get the interpolated velocity
    const Vector<T, LatSet::d> uwall =
        pop.getVelocity() +
        T(0.5) * (pop.getVelocity() -
                  pop.getNeighbor(LatSet::opp[k])->getVelocity());
    pop.getDDF(k) =
        2 * pop.getRho() * LatSet::w[k] *
            (T(1) + pow((uwall * LatSet::c[k]), 2) * T(0.5) * LatSet::InvCs4 -
             uwall.getnorm2() * T(0.5) * LatSet::InvCs2) -
        pop.getDDFpostcol(LatSet::opp[k]);
  }
};
template <typename T, typename LatSet,
          void (*BBLikemethod)(AbstractPopulation<T, LatSet> &, int)>
class BounceBackLike3D_ final : public BasicBoundary3D<T, LatSet> {
 public:
  // methods
  BounceBackLike3D_(int flag = 1, std::string name = "Bounce-Back-Like")
      : BasicBoundary3D<T, LatSet>(name, flag) {}
  void Apply() override {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (auto &flowdir : this->BdPops) {
      auto &pop = flowdir.getBdPop();
      for (const auto &dir : flowdir.getoutflows()) {
        BBLikemethod(pop, dir);
      }
    }
  }
  void Apply(std::vector<PopulationNbr<T, LatSet>> &pops) override {
    #pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (auto &flowdir : this->BdPops) {
      auto &pop = pops[flowdir.getId()];
      for (const auto &dir : flowdir.getoutflows()) {
        BBLikemethod(pop, dir);
      }
    }
  }
};


template <typename T, typename LatSet,
          void (*BBLikemethod)(PopulationNbr<T, LatSet> &, int)>
class BounceBackLike3D final : public BasicBoundary3D<T, LatSet> {
 public:
  // methods
  BounceBackLike3D(int flag = 1, std::string name = "Bounce-Back-Like")
      : BasicBoundary3D<T, LatSet>(name, flag) {}
  void Apply() override {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (auto &flowdir : this->BdPops) {
      auto *pop = flowdir.getBdPop();
      for (const auto &dir : flowdir.getoutflows()) {
        BBLikemethod(*pop, dir);
      }
    }
  }
  void Apply(std::vector<PopulationNbr<T, LatSet>> &pops) override {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (auto &flowdir : this->BdPops) {
      auto &pop = pops[flowdir.getId()];
      for (const auto &dir : flowdir.getoutflows()) {
        BBLikemethod(pop, dir);
      }
    }
  }
};

template <typename T, typename LatSet>
class BoundaryManager3D {
 private:
  std::vector<BasicBoundary3D<T, LatSet> *> _Boundaries;

 public:
  BoundaryManager3D(std::vector<BasicBoundary3D<T, LatSet> *> &boundaries)
      : _Boundaries(boundaries) {}
  template <typename... Args>
  BoundaryManager3D(Args... args) : _Boundaries{args...} {}

  void Setup(std::vector<PopulationNbr<T, LatSet>> &pop) {
    for (BasicBoundary3D<T, LatSet> *boundary : _Boundaries)
      boundary->Setup(pop);
  }
  template <typename PopulationType>
  void Setup_(std::vector<PopulationType> &pop) {
    for (BasicBoundary3D<T, LatSet> *boundary : _Boundaries)
      boundary->Setup(pop);
  }
  void Stream() {
    for (BasicBoundary3D<T, LatSet> *boundary : _Boundaries) boundary->Stream();
  }
  void Apply() {
    for (BasicBoundary3D<T, LatSet> *boundary : _Boundaries) boundary->Apply();
  }
  void printinfo() {
    std::cout << "[Boundary Statistics]: "
              << "\n"
              << "Boundary Type  |  Number of Boundary Cells" << std::endl;
    for (BasicBoundary3D<T, LatSet> *boundary : _Boundaries)
      boundary->getinfo();
  }
};

template <typename T, typename LatSet>
class MovingBoundaryManager3D {
 private:
  BasicBoundary3D<T, LatSet> &MovingBd;
  VoxelGeometry3D<T> &Geo;
  T *StateField;

  std::vector<PopulationNbr<T, LatSet>> *Pop;

 public:
  MovingBoundaryManager3D(BasicBoundary3D<T, LatSet> &movingbd,
                          VoxelGeometry3D<T> &geo, std::vector<int> &statefield,
                          std::vector<PopulationNbr<T, LatSet>> *pop = nullptr)
      : MovingBd(movingbd), Geo(geo), StateField(statefield), Pop(pop) {}

  // communicate
  void Communicate(std::vector<int> &interface) {
    MovingBd.clear();
    for (int id : interface) {
      MovingBd.addtoBd(Geo.getVoxel(id), StateField);
    }
  }
};