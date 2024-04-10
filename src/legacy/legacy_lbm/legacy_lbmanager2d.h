/*Lattice2D Boltzmann method Manager*/
#pragma once

#include "legacy/legacy_lbm/lbm2D.h"

// TODO: Compute_Boundaryrho()

// manage multiple lb2d instances
template <typename T>
class LBManager2D : public Index2D {
 private:
  int Ni;
  int Nj;
  int N;
  // tolerance
  T tol;
  T** Uold = nullptr;
  T* Rhoold = nullptr;

  Velocity2D<T>& Field;
  Basiclbm2D<T>* lbmNS;
  Basiclbm2D<T>* lbmTH;
  Basiclbm2D<T>* lbmSO;
  std::vector<Basiclbm2D<T>*> lbms;

  // souce term: thermal and solutal buoyancy
  T Lattice_gbetaT;
  T Lattice_gbetaC;

 public:
  LBManager2D(Velocity2D<T>& field, Basiclbm2D<T>* lbm0_ = nullptr,
              Basiclbm2D<T>* lbm1_ = nullptr, Basiclbm2D<T>* lbm2_ = nullptr);
  LBManager2D(Velocity2D<T>& field, std::vector<Basiclbm2D<T>*>& lbms_);
  ~LBManager2D();
  void ReOrder_lbms(std::string NS = "NS", std::string TH = "TH",
                    std::string SO = "SO");

  // get pointers
  std::vector<Basiclbm2D<T>*>* getlbms() { return &lbms; }
  Basiclbm2D<T>* getlbmNS() { return lbmNS; }
  Basiclbm2D<T>* getlbmTH() { return lbmTH; }
  Basiclbm2D<T>* getlbmSO() { return lbmSO; }
  // get reference
  Basiclbm2D<T>& getlbmNS_ref() { return *lbmNS; }
  Basiclbm2D<T>& getlbmTH_ref() { return *lbmTH; }
  Basiclbm2D<T>& getlbmSO_ref() { return *lbmSO; }
  Velocity2D<T>& getUField() { return Field; }


  // tolerance
  T Current_Res = 1;
  void SetupToleranceU(T tol_);
  void SetupToleranceRho();
  void calcToleranceU();
  bool NOTConvergedU = true;

  void calcToleranceRho(Basiclbm2D<T>* lb);
  bool NOTConvergedRho = true;

// FORCE
// only add force, no update
  void add_ThermalBuoyancy();
  void add_SolutalBuoyancy();
  // set used to update and add force
  void set_ThermalBuoyancy();
  void set_SolutalBuoyancy();
  void set_ThermalandSolutalBuoyancy();
// Source term
  void add_Conc(PhaseDiagramConverter<T> &PDConv);
  void add_Temp();

// test
  void Simple_Solification();
};
