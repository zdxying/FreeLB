// legacy_basic_boundary.h

#include "legacy/legacy_lbm/legacy_populations.h"

#pragma once


// -----------------------------
// old version of boundary
// -----------------------------

template <typename T, typename LatSet>
class GenericBoundaryCell {
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
  GenericBoundaryCell(PopulationNbr<T, LatSet> *pop)
      : BdPop(pop), Id(pop->getVoxel().getId()) {
    outflows.reserve(LatSet::q);
    inflows.reserve(LatSet::q);
    set();
  }
  GenericBoundaryCell(Voxel<T, LatSet::d> &vox, std::vector<int> &State)
      : Id(vox.getId()), BdPop(nullptr) {
    outflows.reserve(LatSet::q);
    inflows.reserve(LatSet::q);
    set(vox, State);
  }
  GenericBoundaryCell(Voxel<T, LatSet::d> &vox)
      : Id(vox.getId()), BdPop(nullptr) {
    outflows.reserve(LatSet::q);
    inflows.reserve(LatSet::q);
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
    for (int k = 1; k < LatSet::q; ++k) {
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
  void set(const Voxel<T, LatSet::d> &vox, const std::vector<int> &State) {
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
  template <typename FlagType>
  void set(const Voxel<T, LatSet::d> &vox, const std::vector<FlagType> &State,
           FlagType flag) {
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
      if (State[id] != flag) inflows.push_back(LatSet::opp[k]);
      // outflow
      if (State[id] == flag && State[voxnbrs[LatSet::opp[k]]->getId()] != flag)
        outflows.push_back(LatSet::opp[k]);
    }
  }
};

template <typename T, typename LatSet>
class GenericBoundary {
 protected:
  std::vector<GenericBoundaryCell<T, LatSet>> BdPops;
  // geometry flag
  int _flag;
  // name for output info
  std::string _name;
  // boundary values
  // std::vector<T> BoundaryValues;

 public:
  GenericBoundary(std::string name, int flag = 1) : _flag(flag), _name(name) {}

  const int &getFlag() const { return _flag; }
  std::vector<GenericBoundaryCell<T, LatSet>> &getBdPops() { return BdPops; }
  const std::vector<GenericBoundaryCell<T, LatSet>> &getBdPops() const {
    return BdPops;
  }
  void clear() { BdPops.clear(); }
  void addtoBd(GenericBoundaryCell<T, LatSet> &bdpop) {
    BdPops.emplace_back(bdpop);
  }
  void addtoBd(Voxel<T, LatSet::d> &vox, std::vector<int> &State) {
    BdPops.emplace_back(GenericBoundaryCell<T, LatSet>(vox, State));
  }
  void addtoBd(Voxel<T, LatSet::d> &vox) {
    BdPops.emplace_back(GenericBoundaryCell<T, LatSet>(vox));
  }
  // this will be called in lbm3D
  void Setup(std::vector<PopulationNbr<T, LatSet>> &pops) {
    for (PopulationNbr<T, LatSet> &pop : pops) {
      if (pop.getVoxel().getFlag() == _flag)
        BdPops.emplace_back(GenericBoundaryCell<T, LatSet>(&pop));
    }
  }
  // set based on field flag
  template <typename FIELDTYPE>
  void Setup(const std::vector<FIELDTYPE> &field, const FIELDTYPE &flag,
             const std::vector<PopulationNbr<T, LatSet>> &pops) {
    for (const PopulationNbr<T, LatSet> &pop : pops) {
      if (field[pop.getVoxel().getId()] == flag)
        BdPops.emplace_back(GenericBoundaryCell<T, LatSet>(&pop));
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
struct BBlikemethod {
  // pop.f[k] = pop.fpostcol[LatSet::opp[k]];
  // it is more recommended to use bounceback<T, LatSet>::apply
  static inline void normal_bounceback(PopulationNbr<T, LatSet> &pop, int k) {
    pop.f[k] = pop.fpostcol[LatSet::opp[k]];
  }
  // pop.f[k] = 2 * rho * LatSet::w[k] - pop.fpostcol[LatSet::opp[k]];
  // it is more recommended to use antibounceback<T, LatSet>::apply
  static inline void anti_bounceback_simplified(PopulationNbr<T, LatSet> &pop,
                                                int k) {
    pop.f[k] = 2 * pop.rho * LatSet::w[k] - pop.fpostcol[LatSet::opp[k]];
  }

  // pop.f[k] = pop.fpostcol[LatSet::opp[k]] + 2 * LatSet::w[k] * pop.rho * uc *
  // LatSet::InvCs2;
  // call: (pop[Id], Dir, Field.U[Id])
  static inline void movingwall_bounceback(PopulationNbr<T, LatSet> &pop,
                                           int k) {
    pop.f[k] = pop.fpostcol[LatSet::opp[k]] +
               2 * LatSet::w[k] * pop.rho * (pop.getVelocity() * LatSet::c[k]) *
                   LatSet::InvCs2;
  }
  // pop.f[k] = 2 * feq - pop.fpostcol[LatSet::opp[k]]; first order
  // call: (pop[Id], Dir, Field.U[Id])
  static inline void anti_bounceback_O1(PopulationNbr<T, LatSet> &pop, int k) {
    pop.f[k] =
        2 * Equilibrium<T, LatSet>::Order1(k, pop.getVelocity(), pop.rho) -
        pop.fpostcol[LatSet::opp[k]];
  }
  // pop.f[k] = 2 * feq - pop.fpostcol[LatSet::opp[k]]; second order
  // call: (pop[Id], Dir, Field.U[Id])
  static inline void anti_bounceback_O2(PopulationNbr<T, LatSet> &pop, int k) {
    pop.f[k] =
        2 * Equilibrium<T, LatSet>::Order2(k, pop.getVelocity(), pop.rho,
                                           pop.getVelocity().getnorm2()) -
        pop.fpostcol[LatSet::opp[k]];
  }
  // pressure boundary condition: anti bounceback scheme
  // pop.f[k] = 2 * LatSet::w[k] * pop.rho * (1 + LatSet::InvCs2 *
  // LatSet::InvCs2 * uc * uc * T(0.5) - LatSet::InvCs2 * u2 * T(0.5))
  static inline void anti_bounceback_pressure(PopulationNbr<T, LatSet> &pop,
                                              int k) {
    // get the interpolated velocity
    const Vector<T, LatSet::d> uwall =
        pop.getVelocity() +
        T(0.5) * (pop.getVelocity() -
                  pop.getNeighbor(LatSet::opp[k])->getVelocity());
    pop.f[k] =
        2 * pop.rho * LatSet::w[k] *
            (T(1) + pow((uwall * LatSet::c[k]), 2) * T(0.5) * LatSet::InvCs4 -
             uwall.getnorm2() * T(0.5) * LatSet::InvCs2) -
        pop.fpostcol[LatSet::opp[k]];
  }
};

template <typename T, typename LatSet,
          void (*BBLikemethod)(PopulationNbr<T, LatSet> &, int)>
class GenericBounceBackLike final : public GenericBoundary<T, LatSet> {
 public:
  // methods
  GenericBounceBackLike(int flag = 1, std::string name = "Bounce-Back-Like")
      : GenericBoundary<T, LatSet>(name, flag) {}
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
class GenericBoundaryManager {
 private:
  std::vector<GenericBoundary<T, LatSet> *> _Boundaries;

 public:
  GenericBoundaryManager(std::vector<GenericBoundary<T, LatSet> *> &boundaries)
      : _Boundaries(boundaries) {}
  template <typename... Args>
  GenericBoundaryManager(Args... args) : _Boundaries{args...} {}

  void Setup(std::vector<PopulationNbr<T, LatSet>> &pop) {
    for (GenericBoundary<T, LatSet> *boundary : _Boundaries)
      boundary->Setup(pop);
  }
  template <typename PopulationType>
  void Setup_(std::vector<PopulationType> &pop) {
    for (GenericBoundary<T, LatSet> *boundary : _Boundaries)
      boundary->Setup(pop);
  }
  void Stream() {
    for (GenericBoundary<T, LatSet> *boundary : _Boundaries) boundary->Stream();
  }
  void Apply() {
    for (GenericBoundary<T, LatSet> *boundary : _Boundaries) boundary->Apply();
  }
  void printinfo() {
    std::cout << "[Boundary Statistics]: "
              << "\n"
              << "Boundary Type  |  Number of Boundary Cells" << std::endl;
    for (GenericBoundary<T, LatSet> *boundary : _Boundaries)
      boundary->getinfo();
  }
};

template <typename T, typename LatSet>
class GenericMovingBoundary {
 private:
  GenericBoundary<T, LatSet> &MovingBd;
  VoxelGeometry<T, LatSet::d> &Geo;
  std::vector<int> *StateField;

  std::vector<PopulationNbr<T, LatSet>> *Pops;

 public:
  GenericMovingBoundary(GenericBoundary<T, LatSet> &movingbd,
                        VoxelGeometry<T, LatSet::d> &geo,
                        std::vector<int> *statefield,
                        std::vector<PopulationNbr<T, LatSet>> *pop = nullptr)
      : MovingBd(movingbd), Geo(geo), StateField(statefield), Pops(pop) {}

  // communicate
  void Communicate(const std::vector<int> &interface) {
    MovingBd.clear();
    for (int id : interface) {
      MovingBd.addtoBd(Geo.getVoxel(id), *StateField);
    }
  }
  template <typename FlagType>
  void Communicate(const std::vector<int> &interface,
                   const std::vector<FlagType> &State, FlagType flag) {
    MovingBd.clear();
    for (int id : interface) {
      MovingBd.addtoBd(Geo.getVoxel(id));
    }
    // set boundary cell
    for (GenericBoundaryCell<T, LatSet> &bdpop : MovingBd.getBdPops()) {
      bdpop.template set<FlagType>(Geo.getVoxel(bdpop.getId()), State, flag);
    }
  }
};

// direct perform boundary condition
template <typename T, typename LatSet, typename FlagType,
          void (*BBLikemethod)(PopulationNbr<T, LatSet> &, int)>
class DirectGenericBoundary {
 private:
  std::vector<int> &Interface;
  std::vector<FlagType> &State;
  FlagType Flag;
  std::vector<PopulationNbr<T, LatSet>> &Pops;
  T _OMEGA;
  T _1_OMEGA;

 public:
  DirectGenericBoundary(std::vector<int> &interface,
                        std::vector<FlagType> &state, FlagType flag,
                        std::vector<PopulationNbr<T, LatSet>> &pops, T omega)
      : Interface(interface),
        State(state),
        Flag(flag),
        Pops(pops),
        _OMEGA(omega),
        _1_OMEGA(1 - omega) {}

  void Collide_BGKO2() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id : Interface) {
      Pops[id].BGK_O2(_OMEGA, _1_OMEGA);
    }
  }
  void Collide_BGKO1() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id : Interface) {
      Pops[id].BGK_O1(_OMEGA, _1_OMEGA);
    }
  }
  void Stream() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id : Interface) {
      PopulationNbr<T, LatSet> &pop = Pops[id];
      pop.f[0] = pop.fpostcol[0];
      for (int k = 1; k < LatSet::q; ++k) {
        PopulationNbr<T, LatSet> *popnbr = pop.getNeighbor(k);
        if (State[popnbr->getVoxel().getId()] != Flag)
          pop.f[LatSet::opp[k]] = popnbr->fpostcol[LatSet::opp[k]];
      }
    }
  }
  void Stream(FlagType additionalFlag) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id : Interface) {
      PopulationNbr<T, LatSet> &pop = Pops[id];
      pop.f[0] = pop.fpostcol[0];
      for (int k = 1; k < LatSet::q; ++k) {
        PopulationNbr<T, LatSet> *popnbr = pop.getNeighbor(k);
        int idn = popnbr->getVoxel().getId();
        if (State[idn] == additionalFlag)
          pop.f[LatSet::opp[k]] = popnbr->fpostcol[LatSet::opp[k]];
      }
    }
  }
  void ApplyBCs() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id : Interface) {
      PopulationNbr<T, LatSet> &pop = Pops[id];
      for (int k = 1; k < LatSet::q; ++k) {
        int idn = pop.getVoxel().getNeighborId(k);
        int idnopp = pop.getVoxel().getNeighborId(LatSet::opp[k]);
        if (State[idn] == Flag && State[idnopp] != Flag)
          BBLikemethod(pop, LatSet::opp[k]);
      }
    }
  }
  void ApplyBCs(FlagType additionalFlag) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int id : Interface) {
      PopulationNbr<T, LatSet> &pop = Pops[id];
      for (int k = 1; k < LatSet::q; ++k) {
        int idn = pop.getVoxel().getNeighborId(k);
        int idnopp = pop.getVoxel().getNeighborId(LatSet::opp[k]);
        if (State[idn] == Flag && State[idnopp] != Flag)
          BBLikemethod(pop, LatSet::opp[k]);
        else if (State[idn] == additionalFlag &&
                 State[idnopp] != additionalFlag)
          BBLikemethod(pop, LatSet::opp[k]);
      }
    }
  }
};