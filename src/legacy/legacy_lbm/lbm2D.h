/*Lattice2D Boltzmann method 2D implementations*/
#pragma once

#include <string>
#include <vector>

#include "legacy/legacy_boundary/legacyboundary2D.h"
#include "boundary/periodic_boundary.h"
#include "boundary/bounce_back_boundary.h"
#include "legacy/legacy_lbm/force.h"


template <typename T, typename LatSet>
class BasicLBM{
  private:
  BasicLattice<T, LatSet> &Lattice;
  BoundaryManager &BDM;

  public:
  BasicLBM(BasicLattice<T, LatSet> &lattice, BoundaryManager &bdm)
      : Lattice(lattice), BDM(bdm) {}

  template <void (*get_feq)(T*, const Vector<T, LatSet::d>&, T)>
  void Apply() {
    // Lattice.UpdateRhoU(Lattice.getInnerIndex());
    Lattice.template BGK<get_feq>();
    Lattice.Stream();
    BDM.Apply();
  }
  template <void (*get_feq)(T*, const Vector<T, LatSet::d>&, T)>
  void Apply_ColIdx() {
    Lattice.template BGK<get_feq>(Lattice.getIndex());
    Lattice.Stream();
    BDM.Apply();
  }      
};

template <typename T, typename LatSet>
class Genericlbm2D final : public Basiclbm<T> {
 private:
  std::vector<PopulationNbr<T, LatSet>> _Pop;
  // init
  T _Lattice_Rho_Init;
  // Omega, 1 - OMEGA
  T _OMEGA;
  T _1_OMEGA;

  T _Lattice_gbeta;
  /*index*/
  // built-in index, set from Geo, if (Geo.getVoxel(id).getFlag() != -1)
  std::vector<int> _PopIdx;  // pop index = inner + boundary
  // inner pop index, set from Geo, if (Geo.getVoxel(id).getFlag() == 0)
  std::vector<int> _PopInIdx;
  // post-streaming index to update rho
  std::vector<int> _PostRhoIdx;
  // post-streaming index to update u
  std::vector<int> _PostUIdx;
  // index
  std::vector<int> *_Cells;
  std::vector<int> *_InCells;
  std::vector<int> *_PostRhoCells;
  std::vector<int> *_PostUCells;

  std::string _name;
  std::string _namerho;

  VoxelGeometry2D<T> &_Geo;
  VelocityField2D<T> &_Velocity;
  AbstractConverter<T> &_Conv;
  GenericBoundaryManager<T, LatSet> &_BdManager;
  GenericBoundary<T, LatSet> *_MovingBd;

  // tolerance
  std::vector<T> _RhoOld;
  T _RhoRes;
  std::vector<Vector<T, LatSet::d>> _UOld;
  T _URes;

 public:
  Genericlbm2D(VelocityField2D<T> &velocity, AbstractConverter<T> &conv,
               GenericBoundaryManager<T, LatSet> &bdmanager, std::string name,
               std::string rhoname,
               GenericBoundary<T, LatSet> *movingbd = nullptr,
               std::vector<int> *cells = nullptr,
               std::vector<int> *incells = nullptr,
               std::vector<int> *postrhocells = nullptr,
               std::vector<int> *postucells = nullptr);

  /*-------- get ----------*/
  // if at least one nbr cell flag == -1, return -1, else return 0
  T getPoprho(int Id) const override { return _Pop[Id].rho; }
  T getPhysrho(int Id) const override { return _Conv.getPhysRho(_Pop[Id].rho); }
  T getlatrhoinit() const override { return _Lattice_Rho_Init; }
  std::vector<PopulationNbr<T, LatSet>> &getPops() { return _Pop; }
  const std::vector<PopulationNbr<T, LatSet>> &getPops() const { return _Pop; }
  PopulationNbr<T, LatSet> &getPop(int Id) { return _Pop[Id]; }
  const PopulationNbr<T, LatSet> &getPop(int Id) const { return _Pop[Id]; }
  T getOmega() const override { return _OMEGA; }
  T getPopSize() const override { return _Pop.size(); }

  std::string getname() const override { return _name; }
  std::string getnamerho() const override { return _namerho; }

  std::vector<int> *getIdx() override { return _Cells; }
  std::vector<int> *getInIdx() override { return _InCells; }
  std::vector<int> *getPostRhoCells() override { return _PostRhoCells; }
  std::vector<int> *getPostUCells() override { return _PostUCells; }
  void setIdx(std::vector<int> *cells) { _Cells = cells; }
  void setInIdx(std::vector<int> *cells) { _InCells = cells; }
  void setPostRhoCells(std::vector<int> *cells) { _PostRhoCells = cells; }
  void setPostUCells(std::vector<int> *cells) { _PostUCells = cells; }

  VoxelGeometry2D<T> &getGeo() { return _Geo; }
  VelocityField2D<T> &getVelocityField() { return _Velocity; }

  /*---------set---------*/
  void DefaultSetupIndex();
  void defaultSetupBCs() {
    _BdManager.Setup(_Pop);
    _BdManager.printinfo();
  }
  void addtoIndex_From_VoxelFlag(std::vector<int> *idx, int flag) {
    for (const auto &voxel : _Geo.getVoxels()) {
      if (voxel.getFlag() == flag) idx->push_back(voxel.getId());
    }
  }
  void setPopRho(int Id, T rho) override { _Pop[Id].rho = rho; }
  // void setpoprho_From_Vector(std::vector<int> &popidx, T rho_) override;
  void setPopRho_From_VoxelFlag(int flag, T rho);
  template <typename U = int>
  void InitPop_From_Field(const std::vector<U> &idx, const U &flag, T rho);
  // post process of Geometry2D, set flags of boundary cells

  // BGK collision
  // T *feq, T *u, T rho
  // call Collide_BGK: Collide_BGK<Equilibrium<T,LatSet>::Feq_secondOrder>()
  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T)>
  void Collide_BGK();
  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T)>
  void Collide_BGK(const std::vector<int> &idx);
  void Collide_BGKO2();

  void Stream();
  void InnerStream();

  void BCs() { _BdManager.Apply(); }

  // poststream
  template <bool InnerRho = true>
  void PostStreamRho() {
    if constexpr (InnerRho)
      Compute_Rho(*_PostRhoCells);
    else
      Compute_Rho();
  }
  template <bool InnerU = true>
  void PostStreamU() {
    if constexpr (InnerU)
      Compute_U(*_PostUCells);
    else
      Compute_U();
  }
  void Compute_Rho();
  void Compute_Rho(const std::vector<int> &RhoIdx);
  void Compute_U();
  void Compute_U(const std::vector<int> &UIdx);

  // source term
  // for lbmTH and lbmSO
  // add thermal buoyancy or solutal buoyancy to Force in lbmNS
  // call addtoBuoyancy: lbmTH->addtoBuoyancy(lbmTH->get_Idx(),
  // lbmNS->getForce());
  T getBuoyancy(int id) const override {
    return (_Pop[id].rho - _Lattice_Rho_Init) * _Lattice_gbeta;
  }
  void addtoBuoyancy(std::vector<Vector<T, LatSet::d>> &Force) const;
  void addForce(const std::vector<Vector<T, LatSet::d>> &Force);
  void addForce(const std::vector<Vector<T, LatSet::d>> &Force,
                const std::vector<T> &Field, const std::vector<int> &idx);

  // add to popRho
  void addtoRho(int id, T value) override { _Pop[id].rho += value; }
  // add to pop in each directrion
  void addtoPop(int id, T value) override {
    for (int k = 0; k < LatSet::q; ++k) {
      _Pop[id].f[k] += value * LatSet::w[k];
    }
  }
  void addtoPoppostcol(int id, T value) {
    for (int k = 0; k < LatSet::q; ++k) {
      _Pop[id].fpostcol[k] += value * LatSet::w[k];
    }
  }

  ////////////////////
  // apply
  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T),
            bool InnerRho = true, bool InnerU = true>
  void Run() {
    Collide_BGK<get_feq>();
    Stream();
    BCs();
    PostStreamRho<InnerRho>();
    PostStreamU<InnerU>();
  }
  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T),
            bool InnerRho = true>
  void Run_Rho() {
    Collide_BGK<get_feq>();
    Stream();
    BCs();
    PostStreamRho<InnerRho>();
  }

  // apply with moment update first
  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T),
            bool InnerRho = true, bool InnerU = true>
  void apply_() {
    PostStreamRho<InnerRho>();
    PostStreamU<InnerU>();
    Collide_BGK<get_feq>();
    Stream();
    BCs();
  }
  // test
  template <bool InnerRho = true, bool InnerU = true>
  void Run() {
    Collide_BGKO2();
    Stream();
    BCs();
    PostStreamRho<InnerRho>();
    PostStreamU<InnerU>();
  }

  // tolerance
  void EnableToleranceRho() {
    _RhoOld.reserve(_Pop.size());
    for (const auto &pop : _Pop) _RhoOld.push_back(pop.rho);
  }
  void EnableToleranceU() {
    _UOld.reserve(_Pop.size());
    for (const auto &pop : _Pop) _UOld.push_back(pop.getVelocity());
  }

  T getToleranceRho();
  T getToleranceU();
  void WriteStruPoints(int step) const;
  void WriteRho(int step) const;
  void WritePop(int step) const;
};

// basic data structure for LBM performing collision and streaming

// LB inflow and outflow direction for each cell
// BCs handling is performed by setting the corresponding inflow and outflow
// directions
template <typename T>
class Basiclbm2D {
 public:
  // get

  virtual T getPoprho(int Id) = 0;
  virtual T getPhysrho(int Id) = 0;
  virtual T getlatrhoinit() = 0;
  virtual std::string getname() = 0;
  virtual std::string getnamerho() = 0;
  virtual std::vector<int> &get_Idx() = 0;
  virtual std::vector<int> &get_InIdx() = 0;
  virtual Force2D<T> *getForce() = 0;
  // set

  virtual void setPopRho(int Id, T rho) = 0;
  // add

  virtual void addtoBuoyancy(std::vector<int> &popidx, Force2D<T> *Force) = 0;
  virtual void resetBuoyancy() = 0;
  virtual void addtoRho(int id, T value) = 0;
  virtual void addtoPop(int id, T value) = 0;

  // compute
  virtual void compute_rho(int id) = 0;
};


// e.g.: LatSet = D2Q9, lat::LatSet need to be instantiated
template <typename T, typename LatSet, typename LatStru>
class lbm2D final : public Basiclbm2D<T>, public Index2D, public LatStru {
 private:
  int Ni;
  int Nj;
  int N;
  // init
  T Lattice_Rho_Init;
  T OMEGA;
  // 1 - OMEGA
  T _OMEGA;
  T Lattice_gbeta;
  /*index*/
  // built-in index, set from Geo, if (Geo.getVoxel(id).getFlag() != -1)
  std::vector<int> popIdx;  // pop index = inner + boundary
  // inner pop index, set from Geo, if (Geo.getVoxel(id).getFlag() == 0)
  std::vector<int> popInIdx;
  // Index from CellIndexManager2D.Cells, Flag[id] != -1
  std::vector<int> &Cells;
  // Index from CellIndexManager2D.InCells, Flag[id] == 0
  std::vector<int> &InCells;

  // if enable _COLLISION_P, built-in index must be used
  // legacy
#ifdef _COLLISION_P
  std::vector<std::vector<int>> popIdx_Th;
  std::vector<std::vector<int>> popInIdx_Th;
#endif

  std::string name;
  std::string namerho;

  Velocity2D<T> &Field;
  Geometry2DLegacy<T> &Geo;
  AbstractConverter<T> &Conv;
  BoundaryManager2D<T, LatSet, LatStru> &BDM;

  population<T, LatSet> *pop;
  Force2D<T> *Force;
  bool force_enabled = false;

 public:
  lbm2D(Velocity2D<T> &field, AbstractConverter<T> &Conv_,
        BoundaryManager2D<T, LatSet, LatStru> &boundaryman, std::string name_,
        std::string rhoname_);
  lbm2D(Velocity2D<T> &field, AbstractConverter<T> &Conv_,
        BoundaryManager2D<T, LatSet, LatStru> &boundaryman, std::string name_,
        std::string rhoname_, std::vector<int> &cells,
        std::vector<int> &incells);

  void Setup();
  ~lbm2D();

  void EnableForce();
  // add to InCells by Geo.Flag
  void AddtoInIdx(int flag);

  /*-------- get ----------*/
  // if at least one nbr cell flag == -1, return -1, else return 0
  T getPoprho(int Id) override { return pop[Id].rho; }
  T getPhysrho(int Id) override { return Conv.getPhysRho(pop[Id].rho); }
  T getlatrhoinit() override { return Lattice_Rho_Init; }
  std::string getname() override { return name; }
  std::string getnamerho() override { return namerho; }
  // return Cells
  std::vector<int> &get_Idx() override { return Cells; }
  // return InCells
  std::vector<int> &get_InIdx() override { return InCells; }
  Force2D<T> *getForce() override { return Force; }
  population<T, LatSet> *getPop() { return pop; }
  /*---------set---------*/
  void setPopRho(int Id, T rho) override { pop[Id].rho = rho; }
  void setPopRho(int ni, int nj, T rho) {
    pop[Index2D::GetId(ni, nj, Ni)] = rho;
  }
  void setpoprho_From_Vector(std::vector<int> &popidx, T rho_);
  void setPopRho_From_VoxelFlag(int flag, T rho);
  // post process of Geometry2D, set flags of boundary cells

  // BGK collision
  // T *feq, T *u, T rho
  // call Collide_BGK: Collide_BGK<Equilibrium<T,LatSet>::Feq_secondOrder>()
  template <void (*get_feq)(T *, const T *, T)>
  void Collide_BGK(const std::vector<int> &idx);
  template <void (*get_feq)(T *, const T *, T)>
  void Collide_BGK_p(const std::vector<std::vector<int>> &idx_th);
  template <void (*get_feq)(T *, const T *, T)>
  void Collide_BGK_Partial(const std::vector<direction<LatSet::q>> &BdCell_Dirs,
                           T *fraction);

  void Stream(const std::vector<int> &idx);

  template <void (lbm2D<T, LatSet, LatStru>::*rho_method)(),
            void (lbm2D<T, LatSet, LatStru>::*u_method)()>
  void PostStreamrhoU();

  void BCsManager();

  // poststream
  void compute_rho(const std::vector<int> &idx);
  void compute_rho(int id) override;
  void compute_u(const std::vector<int> &idx);
  void poststream_rho();
  void poststream_innerrho();
  void poststream_u();
  void poststream_inneru();

  // source term
  // for lbmTH and lbmSO
  // add thermal buoyancy or solutal buoyancy to Force in lbmNS
  // call addtoBuoyancy: lbmTH->addtoBuoyancy(lbmTH->get_Idx(),
  // lbmNS->getForce());
  void addtoBuoyancy(std::vector<int> &popidx, Force2D<T> *Force) override;
  void resetBuoyancy() override;
  // int k, T f, T *u, T omega
  template <T (*get_force)(int, T *, T *, T)>
  void AddBuoyancytoPop(std::vector<int> &idx);

  // source
  // add to popRho
  void addtoRho(int id, T value) override;
  // add to pop in each directrion
  void addtoPop(int id, T value) override;

  // statistics of lattice field
  // compute average rho, if partially filled with liquid, DO NOT use this
  T getAverRho();

  ////////////////////
  // apply
  // normal apply
  template <void (*get_feq)(T *, const T *, T),
            void (lbm2D<T, LatSet, LatStru>::*rho_method)() =
                &lbm2D<T, LatSet, LatStru>::poststream_rho,
            void (lbm2D<T, LatSet, LatStru>::*u_method)() =
                &lbm2D<T, LatSet, LatStru>::poststream_inneru>
  void apply() {
#ifdef _COLLISION_P
    Collide_BGK_p<get_feq>(popIdx_Th);
#else
    Collide_BGK<get_feq>(Cells);
#endif
    Stream(InCells);
    BCsManager();
    PostStreamrhoU<rho_method, u_method>();
  }

  // apply with idxrho and idxu
  template <void (*get_feq)(T *, const T *, T)>
  void apply(std::vector<int> &idxrho, std::vector<int> &idxu) {
    Collide_BGK<get_feq>(Cells);
    Stream(InCells);
    BCsManager();
    compute_rho(idxrho);
    compute_u(idxu);
  }

  // apply with moment update first
  template <void (*get_feq)(T *, const T *, T),
            void (lbm2D<T, LatSet, LatStru>::*rho_method)() =
                &lbm2D<T, LatSet, LatStru>::poststream_rho,
            void (lbm2D<T, LatSet, LatStru>::*u_method)() =
                &lbm2D<T, LatSet, LatStru>::poststream_inneru>
  void apply_() {
    PostStreamrhoU<rho_method, u_method>();
#ifdef _COLLISION_P
    Collide_BGK_p<get_feq>(popIdx_Th);
#else
    Collide_BGK<get_feq>(Cells);
#endif
    Stream(InCells);
    BCsManager();
  }

  // apply with force
  template <void (*get_feq)(T *, const T *, T),
            T (*get_force)(int, T *, T *, T),
            void (lbm2D<T, LatSet, LatStru>::*rho_method)() =
                &lbm2D<T, LatSet, LatStru>::poststream_rho,
            void (lbm2D<T, LatSet, LatStru>::*u_method)() =
                &lbm2D<T, LatSet, LatStru>::poststream_inneru>
  void applyF() {
    Collide_BGK<get_feq>(Cells);
    AddBuoyancytoPop<get_force>(Cells);
    Stream(InCells);
    BCsManager();
    PostStreamrhoU<rho_method, u_method>();
  }

  // apply with FORCE with idxrho and idxu
  template <void (*get_feq)(T *, const T *, T),
            T (*get_force)(int, T *, T *, T)>
  void applyF(std::vector<int> &idxrho, std::vector<int> &idxu) {
    Collide_BGK<get_feq>(Cells);
    AddBuoyancytoPop<get_force>(Cells);
    Stream(InCells);
    BCsManager();
    compute_rho(idxrho);
    compute_u(idxu);
  }

  // apply without poststreamU
  template <void (*get_feq)(T *, const T *, T),
            void (lbm2D<T, LatSet, LatStru>::*rho_method)() =
                &lbm2D<T, LatSet, LatStru>::poststream_innerrho>
  void applyRho() {
    Collide_BGK<get_feq>(Cells);
    Stream(InCells);
    BCsManager();
    (this->*rho_method)();
  }

  // apply without poststreamU, with cells partially filled with liquid
  template <void (*get_feq)(T *, const T *, T),
            void (lbm2D<T, LatSet, LatStru>::*rho_method)() =
                &lbm2D<T, LatSet, LatStru>::poststream_innerrho>
  void applyRho_Partial(T *fraction) {
    Collide_BGK<get_feq>(Cells);
    Collide_BGK_Partial<get_feq>(BDM.BB->GetBdCell_Dirs(), fraction);
    Stream(InCells);
    BCsManager();
    (this->*rho_method)();
  }
};
