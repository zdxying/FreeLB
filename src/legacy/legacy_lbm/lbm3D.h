// lbm3D
/*Lattice3D Boltzmann method 2D implementations*/
#pragma once

#include <string>
#include <vector>

#include "legacy/legacy_boundary/legacyboundary3D.h"
#include "legacy/legacy_lattice/legacyfield3D.h"
#include "legacy/legacy_lbm/legacy_populations.h"

template <typename T, typename LatSet>
class lbm3D final : public Equilibrium3D<T, LatSet>, public Basiclbm<T> {
 private:
  // density distribution function
  std::vector<PopulationNbr<T, LatSet>> _Pop;
  // init
  T _Lattice_Rho_Init;
  // collision
  T _OMEGA;
  T _1_OMEGA;
  // buoyancy
  T _Lattice_gbeta;
  // name for writing files
  std::string _name;
  std::string _namerho;
  // reference to geometry and velocity field
  VoxelGeometry3D<T> &_Geo;
  VelocityField3D<T> &_Velocity;
  AbstractConverter<T> &_Conv;
  BoundaryManager3D<T, LatSet> &_BdManager;
  // /*index*/
  // built-in index, set from Geo, if (Geo.getVoxel(id).getFlag() != -1)
  std::vector<int> PopInIdx;

  // cell states
  std::vector<int> *State;
  // moving boundary
  BasicBoundary3D<T, LatSet> *MovingBd;

  // post-streaming index to update rho
  std::vector<int> PostRhoIdx;
  // post-streaming index to update u
  std::vector<int> PostUIdx;
  // index
  std::vector<int> *_Cells;
  std::vector<int> *_InCells;
  std::vector<int> *_PostRhoCells;
  std::vector<int> *_PostUCells;

  // tolerance
  std::vector<T> _RhoOld;
  T _RhoRes;
  std::vector<Vector<T, LatSet::d>> _UOld;
  T _URes;

 public:
  // constructor
  lbm3D(VelocityField3D<T> &velocity, AbstractConverter<T> &conv,
        BoundaryManager3D<T, LatSet> &bdmanager, std::string name,
        std::string rhoname, std::vector<int> *state = nullptr,
        BasicBoundary3D<T, LatSet> *movingbd = nullptr, bool bdsetup = true);

  void BdSetup() {
    _BdManager.Setup(_Pop);
    _BdManager.printinfo();
  }
  T getPoprho(int Id) const override { return _Pop[Id].rho; }
  T getPhysrho(int Id) const override { return _Conv.getPhysRho(_Pop[Id].rho); }
  T getlatrhoinit() const override { return _Lattice_Rho_Init; }
  std::vector<PopulationNbr<T, LatSet>> &getPops() { return _Pop; }
  const std::vector<PopulationNbr<T, LatSet>> &getPops() const { return _Pop; }
  PopulationNbr<T, LatSet> &getPop(int Id) { return _Pop[Id]; }
  const PopulationNbr<T, LatSet> &getPop(int Id) const { return _Pop[Id]; }
  T getOmega() const override { return _OMEGA; }
  T getPopSize() const override { return _Pop.size(); }
  T getBuoyancy(int id) const override {
    return (_Pop[id].rho - _Lattice_Rho_Init) * _Lattice_gbeta;
  }
  std::vector<int> *getIdx() override { return _Cells; }
  std::vector<int> *getInIdx() override { return _InCells; }
  std::vector<int> *getPostRhoCells() override { return _PostRhoCells; }
  std::vector<int> *getPostUCells() override { return _PostUCells; }

  std::string getname() const override { return _name; }
  std::string getnamerho() const override { return _namerho; }

  // set
  void setPopRho_From_VoxelFlag(int flag, T rho);
  void addtoPostRhoIdx(int flag);
  void addtoPostUIdx(int flag);
  void setPopRho(int Id, T rho) override { _Pop[Id].rho = rho; }
  void addtoPop(int id, T value) {
    for (int k = 0; k < LatSet::q; ++k) {
      _Pop[id].f[k] += value * LatSet::w[k];
    }
  }
  void addtoRho(int id, T value) { _Pop[id].rho += value; }

  // ---------LB methods------------
  // template <void (*get_feq)(T *, const T *, const T)>
  // void Collide_BGK(const std::vector<int> &idx);

  template <void (Equilibrium3D<T, LatSet>::*Collide)(
      AbstractPopulation<T, LatSet> &)>
  void Collide_BGK();

  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T)>
  void Collide_BGK();
  // check cell state before collide
  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T)>
  void Collide_BGK_StateCheck();

  void Stream();
  // check cell state before stream
  void Stream_StateCheck();
  inline void BCs();

  template <bool InnerRho = true, bool InnerU = true>
  void PostStream() {
    if constexpr (InnerRho)
      Compute_Rho(PostRhoIdx);
    else
      Compute_Rho();
    if constexpr (InnerU)
      Compute_U(PostUIdx);
    else
      Compute_U();
  }
  void Compute_Rho();
  void Compute_Rho(std::vector<int> &RhoIdx);
  void Compute_Rho_Commun();
  void Compute_U();
  void Compute_U(std::vector<int> &UIdx);
  void Compute_U_Commun();

  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T),
            bool InnerRho = true, bool InnerU = true>
  void Run() {
    Collide_BGK<get_feq>();
    Stream();
    BCs();
    PostStream<InnerRho, InnerU>();
  }
  // with cell communication
  template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T),
            bool InnerRho = true, bool InnerU = true>
  void Run_Commun() {
    Collide_BGK_StateCheck<get_feq>();
    Stream_StateCheck();
    MovingBd->Stream(getPops());
    BCs();
    MovingBd->Apply(getPops());
    Compute_Rho_Commun();
    Compute_U_Commun();
    // PostStream<InnerRho, InnerU>();
  }
  template <void (Equilibrium3D<T, LatSet>::*Collide)(
                AbstractPopulation<T, LatSet> &),
            bool InnerRho = true, bool InnerU = true>
  void Run() {
    Collide_BGK<Collide>();
    Stream();
    BCs();
    PostStream<InnerRho, InnerU>();
  }

  // tolerance
  void EnableToleranceRho() {
    _RhoOld.reserve(_Pop.size());
    for (const auto &pop : _Pop) _RhoOld.push_back(pop.rho);
    // _RhoOld.reserve(_Pop.size());
    // for (int i = 0; i < _Pop.size(); ++i) {
    //   T rho = _Pop[i].rho;
    //   _RhoOld.push_back(rho);
    // }
  }
  void EnableToleranceU() {
    _UOld.reserve(_Pop.size());
    for (const auto &pop : _Pop) _UOld.push_back(pop.getVelocity());
    // _UOld.reserve(_Pop.size());
    // for (int i = 0; i < _Pop.size(); ++i) {
    //   T Uold0 = _Pop[i].getVelocity()[0];
    //   T Uold1 = _Pop[i].getVelocity()[1];
    //   T Uold2 = _Pop[i].getVelocity()[2];
    //   _UOld.emplace_back(Uold0, Uold1, Uold2);
    // }
  }
  T getToleranceRho();
  T getToleranceU();
  void WriteStruPoints(int step) const;
  void WriteStruPoints_(int step) const;
  void WriteRho(int step) const;
};
