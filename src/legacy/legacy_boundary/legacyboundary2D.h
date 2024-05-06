// bcs 2d
#pragma once

#include "legacy/legacy_lattice/legacyfield.h"

// flow direction, used to handle boundary cells' streaming and bcs
template <int q>
struct direction {
  int Id;  // cell id
  // used in BCs
  std::vector<int> outflow;
  // inflow from direction opp[x]
  // used in streaming step for boundary cells
  std::vector<int> inflow;

  direction() : Id(-1) {}
  direction(int Id_) : Id(Id_) {
    outflow.reserve(q);
    inflow.reserve(q);
  }
  void Clear() {
    outflow.clear();
    inflow.clear();
  }
  inline void add_outflow(int k) { outflow.emplace_back(k); }
  inline void add_inflow(int k) { inflow.emplace_back(k); }

  template <typename LatSet, typename LatStru>
  inline void set(const int *flags, int (&Nbr)[LatSet::q]) {
    Clear();
    for (int k = 1; k < q; k++) {
      if (flags[Id + Nbr[k]] != -1) {
        add_inflow(LatSet::opp[k]);
      }
      if (flags[Id + Nbr[k]] == -1 && flags[Id + Nbr[LatSet::opp[k]]] != -1) {
        add_outflow(LatSet::opp[k]);
      }
    }
  }
};

// todo: add fixed boundary cells not participating in cell communication

// template <typename T, typename LatSet, typename LatInfo, typename
// BoundaryType> public BoundaryType<T, LatSet>
template <typename T, typename LatSet, typename LatStru>
class BasicBoundary {
 protected:
  std::vector<direction<LatSet::q>> Cell_Dirs;
  std::vector<int> BdCells;
  std::string name;
  // boundary values
  // std::vector<T> BoundaryValues;
  LatStru &Lattice;

  BasicBoundary(LatStru &lattice, std::string name_)
      : Lattice(lattice), name(name_) {
    Cell_Dirs.reserve((Lattice.Ni + Lattice.Nj) * 2);
    BdCells.reserve((Lattice.Ni + Lattice.Nj) * 2);
  }

 public:
  inline void Add(int id) { Cell_Dirs.emplace_back(id); }
  inline void Add(int id, T latrho) {
    Cell_Dirs.emplace_back(id);
    // BoundaryValues.emplace_back(latrho);
  }
  inline void Clear() { Cell_Dirs.clear(); }
  inline int GetBdCellNum() { return Cell_Dirs.size(); }
  // get
  inline std::vector<direction<LatSet::q>> &GetBdCell_Dirs() {
    return Cell_Dirs;
  }
  // just return reference of BdCells, transform is NOT performed
  inline std::vector<int> &GetBdCell_Ids() { return BdCells; }
  // get ids from Cell_Dirs, transform is performed
  inline std::vector<int> &GetBd_Ids() {
    BdCells.clear();
    std::transform(Cell_Dirs.begin(), Cell_Dirs.end(),
                   std::back_inserter(BdCells),
                   [](const direction<LatSet::q> &dir) { return dir.Id; });
    return BdCells;
  }

  void Setup_Cell_Dirs(int *flags);
  void Setup_Cell_Dirs(Geometry2DLegacy<T> &Geo);
  void Setup_Cell_Dirs_(int *flags);
  // set Cell_Dirs[i], i is not cell index
  inline void SetCell_Dirs(int i, int *flags);
  inline void SetCell_Dirs(int i, Geometry2DLegacy<T> &Geo);

  void SetFromGeo(int flag, Geometry2DLegacy<T> &geo);
  void SetFromGeo(int flag, Geometry2DLegacy<T> &geo, T latrho);
  void SetFromFlag(int flag, int *flags, int ni, int nj);
  void SetFromFlag(int flag, int *flags, int ni, int nj, T latrho);
  // user-defined in .cpp
  void SetFromVector(std::vector<int> &cell_ids, int *flags);
  void SetFromVector(std::vector<int> &cell_ids, int *flags, T latrho);
  void SetFromVector(std::vector<int> &cell_ids, int *flags,
                     std::vector<T> &latrhos);

  // apply bcs
  virtual void Apply(population<T, LatSet> *pop) = 0;
  // boundary cell streaming
  void Stream(population<T, LatSet> *pop);
  // print statistics
  inline void statistics();
};

// local bcs: bb, ab, bbmw
template <typename T, typename LatSet>
struct bounceback {
  static void apply(population<T, LatSet> &pop, int k) {
    pop.f[k] = pop.fpostcol[LatSet::opp[k]];
  }
};

template <typename T, typename LatSet, typename LatStru>
struct BounceBack final: public BasicBoundary<T, LatSet, LatStru>,
                    public bounceback<T, LatSet> {
  BounceBack(LatStru &lattice)
      : BasicBoundary<T, LatSet, LatStru>(lattice, "Bounce-Back") {}
  void Apply(population<T, LatSet> *pop) override {
    int Id = 0;
    int Dir = 0;
#pragma omp parallel for private(Id, Dir) num_threads(Thread_Num)
    for (int i = 0; i < this->GetBdCellNum(); i++) {
      Id = this->Cell_Dirs[i].Id;
      for (int k = 0; k < this->Cell_Dirs[i].outflow.size(); k++) {
        Dir = this->Cell_Dirs[i].outflow[k];
        this->apply(pop[Id], Dir);
      }
    }
  }
};

template <typename T, typename LatSet>
struct bouncebackmethod {
  // pop.f[k] = pop.fpostcol[LatSet::opp[k]];
  // it is more recommended to use bounceback<T, LatSet>::apply
  static inline void normal_bounceback(population<T, LatSet> &pop, int k,
                                       const T *u) {
    normal_bounceback(pop, k);
  }
  static inline void normal_bounceback(population<T, LatSet> &pop, int k) {
    pop.f[k] = pop.fpostcol[LatSet::opp[k]];
  }
  // pop.f[k] = 2 * rho * LatSet::w[k] - pop.fpostcol[LatSet::opp[k]];
  // it is more recommended to use antibounceback<T, LatSet>::apply
  static inline void anti_bounceback_simplified(population<T, LatSet> &pop,
                                                int k, const T *u) {
    anti_bounceback_simplified(pop, k);
  }
  static inline void anti_bounceback_simplified(population<T, LatSet> &pop,
                                                int k) {
    pop.f[k] = 2 * pop.rho * LatSet::w[k] - pop.fpostcol[LatSet::opp[k]];
  }

  // pop.f[k] = pop.fpostcol[LatSet::opp[k]] + 2 * LatSet::w[k] * pop.rho * uc *
  // LatSet::InvCs2;
  // call: (pop[Id], Dir, Field.U[Id])
  static inline void movingwall_bounceback(population<T, LatSet> &pop, int k,
                                           const T *u) {
    pop.f[k] = pop.fpostcol[LatSet::opp[k]] +
               2 * LatSet::w[k] * pop.rho * Vect2D<T>::dot(u, LatSet::c[k]) *
                   LatSet::InvCs2;
  }
  // pop.f[k] = 2 * feq - pop.fpostcol[LatSet::opp[k]]; first order
  // call: (pop[Id], Dir, Field.U[Id])
  static inline void anti_bounceback_O1(population<T, LatSet> &pop, int k,
                                        const T *u) {
    pop.f[k] = 2 * Equilibrium<T, LatSet>::Order1(k, u, pop.rho) -
               pop.fpostcol[LatSet::opp[k]];
  }
  // pop.f[k] = 2 * feq - pop.fpostcol[LatSet::opp[k]]; second order
  // call: (pop[Id], Dir, Field.U[Id])
  static inline void anti_bounceback_O2(population<T, LatSet> &pop, int k,
                                        const T *u) {
    pop.f[k] =
        2 * Equilibrium<T, LatSet>::Order2(k, u, pop.rho, Vect2D<T>::sqr(u)) -
        pop.fpostcol[LatSet::opp[k]];
  }
  // pressure boundary condition: anti bounceback scheme
  // pop.f[k] = 2 * LatSet::w[k] * pop.rho * (1 + LatSet::InvCs2 *
  // LatSet::InvCs2 * uc * uc * T(0.5) - LatSet::InvCs2 * u2 * T(0.5))
  static inline void anti_bounceback_pressure(population<T, LatSet> &pop, int k,
                                              const T *u) {
    pop.f[k] =
        2 * pop.rho * LatSet::w[k] *
            (T(1) +
             pow(LatSet::InvCs2 * Vect2D<T>::dot(u, LatSet::c[k]), 2) * T(0.5) -
             LatSet::InvCs2 * Vect2D<T>::sqr(u) * T(0.5)) -
        pop.fpostcol[LatSet::opp[k]];
  }
};

template <typename T, typename LatSet, typename LatStru,
          void (*bcmethod)(population<T, LatSet> &, int, const T *)>
struct BounceBackLike final: public BasicBoundary<T, LatSet, LatStru>,
                        public bounceback<T, LatSet> {
  // members
  Velocity2D<T> &Field;
  // methods
  BounceBackLike(LatStru &lattice, Velocity2D<T> &field,
                 std::string name = "Bounce-Back-Like")
      : BasicBoundary<T, LatSet, LatStru>(lattice, name), Field(field) {}
  void Apply(population<T, LatSet> *pop) override {
    int Id = 0;
    int Dir = 0;
#pragma omp parallel for private(Id, Dir) num_threads(Thread_Num)
    for (int i = 0; i < this->GetBdCellNum(); i++) {
      Id = this->Cell_Dirs[i].Id;
      for (int k = 0; k < this->Cell_Dirs[i].outflow.size(); k++) {
        Dir = this->Cell_Dirs[i].outflow[k];
        bcmethod(pop[Id], Dir, Field.U[Id]);
      }
    }
  }
};

template <typename T, typename LatSet>
struct antibounceback {
  // simplified version when wall velocity is zero
  static void apply(population<T, LatSet> &pop, int k, T rho) {
    pop.f[k] = 2 * rho * LatSet::w[k] - pop.fpostcol[LatSet::opp[k]];
  }
};
template <typename T, typename LatSet, typename LatStru>
struct AntiBounceBack final: public BasicBoundary<T, LatSet, LatStru>,
                        public antibounceback<T, LatSet> {
  AntiBounceBack(LatStru &lattice)
      : BasicBoundary<T, LatSet, LatStru>(lattice, "Anti-Bounce-Back") {}
  void Apply(population<T, LatSet> *pop) override {
    int Id = 0;
    int Dir = 0;
#pragma omp parallel for private(Id, Dir) num_threads(Thread_Num)
    for (int i = 0; i < this->GetBdCellNum(); i++) {
      Id = this->Cell_Dirs[i].Id;
      for (int k = 0; k < this->Cell_Dirs[i].outflow.size(); k++) {
        Dir = this->Cell_Dirs[i].outflow[k];
        this->apply(pop[Id], Dir, pop[Id].rho);
      }
    }
  }
};

template <typename T, typename LatSet>
struct bouncebackmovingwall {
  static void apply(population<T, LatSet> &pop, int k, T uc) {
    pop.f[k] = pop.fpostcol[LatSet::opp[k]] +
               2 * LatSet::w[k] * pop.rho * uc * LatSet::InvCs2;
    // pop.f[i] = pop.fpostcol[opp9[i]] + 2 * w[i] * pop.rho *
    // (Lattice_U_Wall[0] * c[i][0] + Lattice_U_Wall[1] * c[i][1]) *
    // LatSet::InvCs2 / (1 - Lattice_U_Wall[0]);
  }
};
template <typename T, typename LatSet, typename LatStru>
struct BounceBackMovingWall final: public BasicBoundary<T, LatSet, LatStru>,
                              public bouncebackmovingwall<T, LatSet>,
                              public Vect2D<T> {
  Velocity2D<T> &Field;
  BounceBackMovingWall(LatStru &lattice, Velocity2D<T> &field)
      : BasicBoundary<T, LatSet, LatStru>(lattice, "Bounce-Back-Moving-Wall"),
        Field(field) {}
  void Apply(population<T, LatSet> *pop) override {
    int Id = 0;
    int Dir = 0;
#pragma omp parallel for private(Id, Dir) num_threads(Thread_Num)
    for (int i = 0; i < this->GetBdCellNum(); i++) {
      Id = this->Cell_Dirs[i].Id;
      for (int k = 0; k < this->Cell_Dirs[i].outflow.size(); k++) {
        Dir = this->Cell_Dirs[i].outflow[k];
        this->apply(pop[Id], Dir, Vect2D<T>::dot(Field.U[Id], LatSet::c[Dir]));
      }
    }
  }
};

/*using neighbor cells*/
template <typename T, typename LatSet>
struct periodic {
  static void apply(population<T, LatSet> *pop, int id1, int id2, int k) {
    pop[id1].f[k] = pop[id2].fpostcol[k];
  }
};
template <typename T, typename LatSet, typename LatStru>
struct Periodic final: public BasicBoundary<T, LatSet, LatStru>,
                  public periodic<T, LatSet> {
  Periodic(LatStru &lattice)
      : BasicBoundary<T, LatSet, LatStru>(lattice, "Periodic") {}
  void Apply(population<T, LatSet> *pop) override {
    int Id = 0;
    int Dir = 0;
#pragma omp parallel for private(Id, Dir) num_threads(Thread_Num)
    for (int i = 0; i < this->GetBdCellNum(); i++) {
      Id = this->Cell_Dirs[i].Id;
      for (int k = 0; k < this->Cell_Dirs[i].outflow.size(); k++) {
        Dir = this->Cell_Dirs[i].outflow[k];
        this->apply(pop, Id, Id + this->Lattice.Per[Dir], Dir);
      }
    }
  }
};

// should be a post-collision operator
template <typename T, typename LatSet>
struct generalisedperiodic {
  static void apply(population<T, LatSet> *pop, int id1, int id2, int k,
                    Velocity2D<T> &field) {
    // use id2's Velocity
    T feq1 = Equilibrium<T, LatSet>::Order1_Incompresible(k, field.U[id2],
                                                          pop[id1].rho);
    T feq2 = Equilibrium<T, LatSet>::Order1_Incompresible(k, field.U[id2],
                                                          pop[id2].rho);
    pop[id1].f[k] = feq1 + pop[id2].fpostcol[k] - feq2;
  }
};
template <typename T, typename LatSet, typename LatStru>
struct GeneralisedPeriodic final: public BasicBoundary<T, LatSet, LatStru>,
                             public generalisedperiodic<T, LatSet> {
  Velocity2D<T> &Field;
  T boundary_rho;
  GeneralisedPeriodic(LatStru &lattice, Velocity2D<T> &field,
                      T boundary_rho_ = 1)
      : BasicBoundary<T, LatSet, LatStru>(lattice, "Generalised-Periodic"),
        Field(field),
        boundary_rho(boundary_rho_) {}
  void Apply(population<T, LatSet> *pop) override {
    int Id = 0;
    int Dir = 0;
#pragma omp parallel for private(Id, Dir) num_threads(Thread_Num)
    for (int i = 0; i < this->GetBdCellNum(); i++) {
      Id = this->Cell_Dirs[i].Id;
      for (int k = 0; k < this->Cell_Dirs[i].outflow.size(); k++) {
        Dir = this->Cell_Dirs[i].outflow[k];
        this->apply(pop, Id, Id + this->Lattice.Per[Dir], Dir, Field);
      }
    }
  }
};

template <typename T, typename LatSet, typename LatStru>
struct BoundaryManager2D {
  BounceBack<T, LatSet, LatStru> *BB = nullptr;
  AntiBounceBack<T, LatSet, LatStru> *ABB = nullptr;
  BounceBackMovingWall<T, LatSet, LatStru> *BBMW = nullptr;
  Periodic<T, LatSet, LatStru> *Per = nullptr;
  GeneralisedPeriodic<T, LatSet, LatStru> *GP = nullptr;

  std::vector<BasicBoundary<T, LatSet, LatStru> *> Bds;

  BoundaryManager2D(BounceBack<T, LatSet, LatStru> *bb = nullptr,
                    AntiBounceBack<T, LatSet, LatStru> *abb = nullptr,
                    BounceBackMovingWall<T, LatSet, LatStru> *bbmw = nullptr,
                    Periodic<T, LatSet, LatStru> *per = nullptr,
                    GeneralisedPeriodic<T, LatSet, LatStru> *gp = nullptr)
      : BB(bb), ABB(abb), BBMW(bbmw), Per(per), GP(gp), Bds() {
    if (BB != nullptr) Bds.push_back(BB);
    if (ABB != nullptr) Bds.push_back(ABB);
    if (BBMW != nullptr) Bds.push_back(BBMW);
    if (Per != nullptr) Bds.push_back(Per);
    if (GP != nullptr) Bds.push_back(GP);
    statistic();
  }
  BoundaryManager2D(std::vector<BasicBoundary<T, LatSet, LatStru> *> &boundary)
      : Bds(boundary) {
    Statistics();
  }

  void Apply(population<T, LatSet> *pop) {
    for (int i = 0; i < Bds.size(); i++) {
      Bds[i]->Apply(pop);
    }
  }
  void Stream(population<T, LatSet> *pop) {
    for (int i = 0; i < Bds.size(); i++) {
      Bds[i]->Stream(pop);
    }
  }

  void statistic() {
    std::cout << "[Boundary Statistics]: " << std::endl;
    if (BB != nullptr) std::cout << "BB: " << BB->GetBdCellNum() << std::endl;
    if (ABB != nullptr) std::cout << "ABB: " << ABB->GetBdCellNum() << std::endl;
    if (BBMW != nullptr) std::cout << "BBMW: " << BBMW->GetBdCellNum() << std::endl;
    if (Per != nullptr) std::cout << "Per: " << Per->GetBdCellNum() << std::endl;
    if (GP != nullptr) std::cout << "GP: " << GP->GetBdCellNum() << std::endl;
  }

  void Statistics() {
    std::cout << "[Boundary Statistics]: "
              << "\n"
              << "Boundary Type  |  Number of Boundary Cells"
              << "\n";
    for (int i = 0; i < Bds.size(); i++) {
      Bds[i]->statistics();
    }
  }
};

// other bcs method
template <typename T, typename LatSet>
T getRhoInlet(BaseConverter<T> &Conv, int Ni, int Nj, T deltaP = 0) {
  if (deltaP == 0) {
    deltaP = 8 * Conv.Lattice_VisKine * Conv.Lattice_charU * Ni / Nj / Nj;
  }
  return deltaP * LatSet::InvCs2 + 1;
}