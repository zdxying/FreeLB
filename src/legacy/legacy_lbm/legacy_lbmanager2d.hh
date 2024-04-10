/*Lattice2D Boltzmann method Manager*/

#include "legacy/legacy_lbm/legacy_lbmanager2d.h"

template <typename T>
LBManager2D<T>::LBManager2D(Velocity2D<T>& field, Basiclbm2D<T>* lbm0_,
                            Basiclbm2D<T>* lbm1_, Basiclbm2D<T>* lbm2_)
    : Field(field),
      Ni(field.Geo.getNi()),
      Nj(field.Geo.getNj()),
      N(field.Geo.getN()),
      lbms() {
  if (lbm0_ != nullptr) lbms.push_back(lbm0_);
  if (lbm1_ != nullptr) lbms.push_back(lbm1_);
  if (lbm2_ != nullptr) lbms.push_back(lbm2_);
  if (lbms.size() == 0) {
    std::cout << "No lb method!" << std::endl;
    exit(-1);
  }
  ReOrder_lbms();
}
template <typename T>
LBManager2D<T>::LBManager2D(Velocity2D<T>& field,
                            std::vector<Basiclbm2D<T>*>& lbms_)
    : Field(field),
      lbms(lbms_),
      Ni(field.Geo.getNi()),
      Nj(field.Geo.getNj()),
      N(field.Geo.getN()) {
  // lbmNS = lbms[0];
  // if (lbms.size() > 1) lbmTH = lbms[1];
  // if (lbms.size() > 2) lbmSO = lbms[2];
  if (lbms.size() == 0) {
    std::cout << "No lb method!" << std::endl;
    exit(-1);
  }
  ReOrder_lbms();
}

template <typename T>
LBManager2D<T>::~LBManager2D() {
  if (Uold != nullptr) {
    for (int i = 0; i < 2; i++) {
      delete[] Uold[i];
    }
    delete[] Uold;
  }
  if (Rhoold != nullptr) {
    delete[] Rhoold;
  }
}

template <typename T>
void LBManager2D<T>::ReOrder_lbms(std::string NS, std::string TH,
                                  std::string SO) {
  for (int i = 0; i < lbms.size(); i++) {
    if (lbms[i]->getname() == NS)
      lbmNS = lbms[i];
    else if (lbms[i]->getname() == TH)
      lbmTH = lbms[i];
    else if (lbms[i]->getname() == SO)
      lbmSO = lbms[i];
  }
}

/*tolerance*/
template <typename T>
void LBManager2D<T>::SetupToleranceU(T tol_) {
  tol = tol_;
  Uold = new T*[N];
#pragma omp parallel for num_threads(Thread_Num)
  for (int i = 0; i < N; i++) {
    Uold[i] = new T[2];
    Uold[i][0] = T(0);
    Uold[i][1] = T(0);
  }
}
template <typename T>
void LBManager2D<T>::SetupToleranceRho() {
  Rhoold = new T[N];
  Traverse_Bulk_OMP(Ni, Nj, 1, [this](int id) { Rhoold[id] = T(0); });
}

template <typename T>
void LBManager2D<T>::calcToleranceU() {
  T res0, res1, res;
  T maxres = T(0);
#ifdef _OPENMP
  T thread_maxres[Thread_Num] = {T(0)};  // for each thread, get the max res
#endif
  int Id = 0;
#pragma omp parallel for private(Id, res0, res1, res) num_threads(Thread_Num)
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      Id = GetId(i, j, Ni);
      res0 = std::abs(Field.U[Id][0] - Uold[Id][0]);
      res1 = std::abs(Field.U[Id][1] - Uold[Id][1]);
      res = std::max(res0, res1);
// if define omp parallel, then use thread_maxres
#ifdef _OPENMP
      if (res > thread_maxres[omp_get_thread_num()]) {
        thread_maxres[omp_get_thread_num()] = res;
      }
#else
      if (res > maxres) {
        maxres = res;
      }
#endif
    }
  }

#ifdef _OPENMP
  for (int i = 0; i < Thread_Num; i++) {
    if (thread_maxres[i] > maxres) {
      maxres = thread_maxres[i];
    }
  }
#endif

  Current_Res = maxres;
  if (Current_Res < tol) {
    std::cout << "Converged! Res: " << Current_Res << std::endl;
    NOTConvergedU = false;
  } else {
    if (Current_Res > 10) {
      std::cout << "Diverged! Res: " << Current_Res << std::endl;
      exit(-1);
    }
    Traverse_Bulk_OMP(Ni, Nj, 1, [this](int id) {
      Uold[id][0] = Field.U[id][0];
      Uold[id][1] = Field.U[id][1];
    });
  }
}

template <typename T>
void LBManager2D<T>::calcToleranceRho(Basiclbm2D<T>* lb) {
  T res;
  T maxres = T(0);
#ifdef _OPENMP
  T thread_maxres[Thread_Num] = {T(0)};  // for each thread, get the max res
#endif
  int Id = 0;
#pragma omp parallel for private(Id, res) num_threads(Thread_Num)
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      Id = GetId(i, j, Ni);
      res = std::abs(lb->getPoprho() - Rhoold[Id]);
#ifdef _OPENMP
      if (res > thread_maxres[omp_get_thread_num()]) {
        thread_maxres[omp_get_thread_num()] = res;
      }
#else
      if (res > maxres) {
        maxres = res;
      }
#endif
    }
  }
#ifdef _OPENMP
  for (int i = 0; i < Thread_Num; i++) {
    if (thread_maxres[i] > maxres) {
      maxres = thread_maxres[i];
    }
  }
#endif

  Current_Res = maxres;
  if (Current_Res < tol) {
    std::cout << "Converged! Res: " << Current_Res << std::endl;
    NOTConvergedRho = false;
  } else {
    if (Current_Res > 10) {
      std::cout << "Diverged! Res: " << Current_Res << std::endl;
      exit(-1);
    }
#pragma omp parallel for private(Id) num_threads(Thread_Num)
    for (int j = 1; j < Nj - 1; j++) {
      for (int i = 1; i < Ni - 1; i++) {
        Id = GetId(i, j, Ni);
        Rhoold[Id] = lb->getPoprho();
      }
    }
  }
}

template <typename T>
void LBManager2D<T>::add_ThermalBuoyancy() {
  if (lbmTH == nullptr) {
    std::cout << "No thermal lbm!" << std::endl;
    exit(-1);
  }
  lbmTH->addtoBuoyancy(lbmTH->get_Idx(), lbmNS->getForce());
}
template <typename T>
void LBManager2D<T>::add_SolutalBuoyancy() {
  if (lbmSO == nullptr) {
    std::cout << "No solutal lbm!" << std::endl;
    exit(-1);
  }
  lbmSO->addtoBuoyancy(lbmSO->get_Idx(), lbmNS->getForce());
}
template <typename T>
void LBManager2D<T>::set_ThermalBuoyancy() {
  lbmNS->resetBuoyancy();
  add_ThermalBuoyancy();
}
template <typename T>
void LBManager2D<T>::set_SolutalBuoyancy() {
  lbmNS->resetBuoyancy();
  add_SolutalBuoyancy();
}
template <typename T>
void LBManager2D<T>::set_ThermalandSolutalBuoyancy() {
  lbmNS->resetBuoyancy();
  add_ThermalBuoyancy();
  add_SolutalBuoyancy();
}

template <typename T>
void LBManager2D<T>::Simple_Solification() {}
