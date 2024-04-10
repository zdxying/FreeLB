#include "legacy/legacy_lbm/lbm2D.h"

template <typename T, typename LatSet>
Genericlbm2D<T, LatSet>::Genericlbm2D(
    VelocityField2D<T> &velocity, AbstractConverter<T> &conv,
    GenericBoundaryManager<T, LatSet> &bdmanager, std::string name,
    std::string rhoname, GenericBoundary<T, LatSet> *movingbd,
    std::vector<int> *cells, std::vector<int> *incells,
    std::vector<int> *postrhocells, std::vector<int> *postucells)
    : _Velocity(velocity),
      _Geo(velocity.getGeo()),
      _Conv(conv),
      _BdManager(bdmanager),
      _name(name),
      _namerho(rhoname),
      _OMEGA(conv.GetOMEGA()),
      _1_OMEGA(T(1) - conv.GetOMEGA()),
      _Lattice_Rho_Init(conv.getLatRhoInit()),
      _Lattice_gbeta(conv.GetLattice_gbeta()),
      _MovingBd(movingbd),
      _Cells(cells),
      _InCells(incells),
      _PostRhoCells(postrhocells),
      _PostUCells(postucells) {
  _Pop.reserve(_Geo.getVoxels().size());
  for (int i = 0; i < _Geo.getVoxels().size(); ++i) {
    _Pop.emplace_back(_Geo.getVoxel(i), _Velocity.getVelocity(i),
                      _Lattice_Rho_Init);
    // set PopInIdx
    // if (_Geo.getVoxel(i).getFlag() == 0) _PopInIdx.push_back(i);
  }
  // set pop neighbor
  for (int i = 0; i < _Geo.getVoxels().size(); ++i) {
    const std::vector<Voxel<T, LatSet::d> *> &voxnbrs =
        _Geo.getVoxel(i).getNeighborList();
    for (int k = 1; k < LatSet::q; ++k) {
      // NOTE THAT:
      // add INVERSE direction(inflow source) to facilitate streaming
      // e.g., direction 1 {1, 0} (right hand side) != nullptr
      //       set to direction 3 {-1, 0} (left hand side)
      if (voxnbrs[k] != nullptr) {
        if (voxnbrs[k]->getFlag() != -1)
          _Pop[i].setNeighbor(LatSet::opp[k], &_Pop[voxnbrs[k]->getId()]);
      }
    }
  }
}

template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::DefaultSetupIndex() {
  // set post-streaming index
  if (_Cells == nullptr) {
    _PopIdx.reserve(_Pop.size());
    for (int i = 0; i < _Pop.size(); ++i) {
      if (_Geo.getVoxel(i).getFlag() != -1) _PopIdx.push_back(i);
    }
    _Cells = &_PopIdx;
  }
  if (_InCells == nullptr) {
    _PopInIdx.reserve(_Pop.size());
    for (int i = 0; i < _Pop.size(); ++i) {
      if (_Geo.getVoxel(i).getFlag() == 0) _PopInIdx.push_back(i);
    }
    _InCells = &_PopInIdx;
  }
  if (_PostRhoCells == nullptr) {
    _PostRhoIdx.reserve(_Pop.size());
    for (int i = 0; i < _Pop.size(); ++i) {
      if (_Geo.getVoxel(i).getFlag() == 0) _PostRhoIdx.push_back(i);
    }
    _PostRhoCells = &_PostRhoIdx;
  }
  if (_PostUCells == nullptr) {
    _PostUIdx.reserve(_Pop.size());
    for (int i = 0; i < _Pop.size(); ++i) {
      if (_Geo.getVoxel(i).getFlag() == 0) _PostUIdx.push_back(i);
    }
    _PostUCells = &_PostUIdx;
  }
}

template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::setPopRho_From_VoxelFlag(int flag, T rho) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (auto &pop : _Pop) {
    if (pop.getVoxel().getFlag() == flag) pop.Initialize(rho);
  }
}
template <typename T, typename LatSet>
template <typename U>
void Genericlbm2D<T, LatSet>::InitPop_From_Field(const std::vector<U> &idx,
                                                   const U &flag, T rho) {
  for (auto &pop : _Pop) {
    if (idx[pop.getVoxel().getId()] == flag) pop.Initialize(rho);
  }
}

template <typename T, typename LatSet>
template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T)>
void Genericlbm2D<T, LatSet>::Collide_BGK() {
  T feq[LatSet::q] = {T(0)};
#pragma omp parallel for private(feq) num_threads(Thread_Num) schedule(static)
  for (int id : *_Cells) {
    auto &pop = _Pop[id];
    get_feq(feq, pop.getVelocity(), pop.rho);
    for (int k = 0; k < LatSet::q; ++k) {
      pop.fpostcol[k] = _OMEGA * feq[k] + _1_OMEGA * pop.f[k];
    }
  }
}
template <typename T, typename LatSet>
template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T)>
void Genericlbm2D<T, LatSet>::Collide_BGK(const std::vector<int> &idx) {
  T feq[LatSet::q] = {T(0)};
#pragma omp parallel for private(feq) num_threads(Thread_Num) schedule(static)
  for (int id : idx) {
    auto &pop = _Pop[id];
    get_feq(feq, pop.getVelocity(), pop.rho);
    for (int k = 0; k < LatSet::q; ++k) {
      pop.fpostcol[k] = _OMEGA * feq[k] + _1_OMEGA * pop.f[k];
    }
  }
}
template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::Collide_BGKO2() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : *_Cells) {
    _Pop[id].BGK_O2(_OMEGA, _1_OMEGA);
  }
}

template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::Stream() {
// inner cell streaming
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : *_InCells) {
    auto &pop = _Pop[id];
    pop.f[0] = pop.fpostcol[0];  // IMPROTANT!!!
    for (int k = 1; k < LatSet::q; ++k) {
      pop.f[k] = pop.getNeighbor(k)->fpostcol[k];
    }
  }
  // boundary cell streaming
  _BdManager.Stream();
}

template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::InnerStream() {
  // inner cell streaming
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : *_InCells) {
    auto &pop = _Pop[id];
    pop.f[0] = pop.fpostcol[0];  // IMPROTANT!!!
    for (int k = 1; k < LatSet::q; ++k) {
      pop.f[k] = pop.getNeighbor(k)->fpostcol[k];
    }
  }
}

template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::Compute_Rho() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (auto &pop : _Pop) {
    pop.Compute_rho();
  }
}
template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::Compute_Rho(const std::vector<int> &RhoIdx) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : RhoIdx) {
    _Pop[id].Compute_rho();
  }
}
template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::Compute_U() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (auto &pop : _Pop) {
    pop.Compute_U();
  }
}
template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::Compute_U(const std::vector<int> &UIdx) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : UIdx) {
    _Pop[id].Compute_U();
  }
}

template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::addtoBuoyancy(
    std::vector<Vector<T, LatSet::d>> &Force) const {
  const int dim = LatSet::d - 1;
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : *_Cells) {
    Force[id][dim] += (_Pop[id].rho - _Lattice_Rho_Init) * _Lattice_gbeta;
  }
}

template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::addForce(
    const std::vector<Vector<T, LatSet::d>> &Force) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : *_Cells) {
    for (int k = 0; k < LatSet::q; ++k) {
      // get c * u * LatSet::InvCs4
      T cu = LatSet::c[k] * _Pop[id].getVelocity() * LatSet::InvCs4;
      // (c - u) * LatSet::InvCs2
      Vector<T, LatSet::d> c_u =
          (LatSet::c[k] - _Pop[id].getVelocity()) * LatSet::InvCs2;
      T force = Force[id] * (c_u + cu * LatSet::c[k]);
      force *= LatSet::w[k] * (T(1) - _OMEGA * T(0.5));
      _Pop[id].f[k] += force;
    }
  }
}

template <typename T, typename LatSet>
inline void Genericlbm2D<T, LatSet>::addForce(
    const std::vector<Vector<T, LatSet::d>> &Force,
    const std::vector<T> &Field, const std::vector<int> &idx) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (int id : idx) {
    for (int k = 0; k < LatSet::q; ++k) {
      // get c * u * LatSet::InvCs4
      T cu = LatSet::c[k] * _Pop[id].getVelocity() * LatSet::InvCs4;
      // (c - u) * LatSet::InvCs2
      Vector<T, LatSet::d> c_u =
          (LatSet::c[k] - _Pop[id].getVelocity()) * LatSet::InvCs2;
      T force = Force[id] * (c_u + cu * LatSet::c[k]);
      force *= LatSet::w[k] * (T(1) - _OMEGA * T(0.5));
      _Pop[id].fpostcol[k] += force * Field[id];
    }
  }
}

// force
// template <typename T, typename LatSet>
// void Genericlbm2D<T, LatSet>::addtoBuoyancy(std::vector<int> &popidx,
//                                             Force2D<T> *force) {
//   int Id = 0;
// #pragma omp parallel for private(Id) num_threads(Thread_Num)
//   for (int i = 0; i < popidx.size(); i++) {
//     Id = popidx.at(i);
//     force[Id].F[1] += (pop[Id].rho - Lattice_Rho_Init) * Lattice_gbeta;
//   }
// }
// template <typename T, typename LatSet>
// void Genericlbm2D<T, LatSet>::resetBuoyancy() {
// #pragma omp parallel for private(Id) num_threads(Thread_Num)
//   for (Id = 0; Id < N; Id++) {
//     Force[Id].clear();
//   }
// }

template <typename T, typename LatSet>
T Genericlbm2D<T, LatSet>::getToleranceRho() {
  T res;
  T maxres = T(0);
  for (int i = 0; i < _Pop.size(); ++i) {
    res = std::abs(_Pop[i].rho - _RhoOld[i]);
    maxres = std::max(res, maxres);
    // set _RhoOld
    _RhoOld[i] = _Pop[i].rho;
  }
  return maxres;
}
template <typename T, typename LatSet>
T Genericlbm2D<T, LatSet>::getToleranceU() {
  T res0, res1, res;
  T maxres = T(0);
  for (int i = 0; i < _Pop.size(); ++i) {
    res0 = std::abs(_Pop[i].getVelocity()[0] - _UOld[i][0]);
    res1 = std::abs(_Pop[i].getVelocity()[1] - _UOld[i][1]);
    res = std::max(res0, res1);
    maxres = std::max(res, maxres);
    // set UOld
    _UOld[i][0] = _Pop[i].getVelocity()[0];
    _UOld[i][1] = _Pop[i].getVelocity()[1];
  }
  return maxres;
}

template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::WriteRho(int step) const {
  int _Nx = _Geo.getNx();
  int _Ny = _Geo.getNy();
  const std::vector<int> &_GlobalIdx = _Geo.getGlobalIdx();
  // write to file
  DirCreator::Create_Dir("./vtkoutput");
  std::string fullName =
      "./vtkoutput/" + _namerho + std::to_string(step) + ".vtk";
  std::ofstream f(fullName.c_str());
  f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
  f << "Voxels" << std::endl << std::flush;
  f << "ASCII" << std::endl << std::flush;

  f << "DATASET STRUCTURED_POINTS" << std::endl << std::flush;
  f << "DIMENSIONS " << _Nx << " " << _Ny << " 1" << std::endl << std::flush;
  f << "ORIGIN " << _Geo.getMin()[0] + _Geo.getVoxelSize() * T(0.5) << " "
    << _Geo.getMin()[1] + _Geo.getVoxelSize() * T(0.5) << " 0" << std::endl
    << std::flush;
  f << "ASPECT_RATIO " << _Geo.getVoxelSize() << " " << _Geo.getVoxelSize()
    << " 1" << std::endl
    << std::flush;
  f << "POINT_DATA " << _Nx * _Ny << std::endl << std::flush;

  std::stringstream rho;
  if constexpr (std::is_same<T, double>::value) {
    rho << "SCALARS rho double" << std::endl << std::flush;
  } else {
    rho << "SCALARS rho float" << std::endl << std::flush;
  }
  rho << "LOOKUP_TABLE default" << std::endl << std::flush;
  for (int i = 0; i < _Nx * _Ny; ++i) {
    if (_GlobalIdx[i] != -1) {
      const PopulationNbr<T, LatSet> &pop = _Pop[_GlobalIdx[i]];
      rho << pop.rho << " ";
    } else {
      rho << 0 << " ";
    }
  }
  rho << std::endl;
  f << rho.str();
  f.close();
}

template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::WritePop(int step) const {
  int _Nx = _Geo.getNx();
  int _Ny = _Geo.getNy();
  const std::vector<int> &_GlobalIdx = _Geo.getGlobalIdx();
  // write to file
  DirCreator::Create_Dir("./vtkoutput");
  std::string fullName = "./vtkoutput/Pop" + std::to_string(step) + ".vtk";
  std::ofstream f(fullName.c_str());
  f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
  f << "Voxels" << std::endl << std::flush;
  f << "ASCII" << std::endl << std::flush;

  f << "DATASET STRUCTURED_POINTS" << std::endl << std::flush;
  f << "DIMENSIONS " << _Nx << " " << _Ny << " 1" << std::endl << std::flush;
  f << "ORIGIN " << _Geo.getMin()[0] + _Geo.getVoxelSize() * T(0.5) << " "
    << _Geo.getMin()[1] + _Geo.getVoxelSize() * T(0.5) << " 0" << std::endl
    << std::flush;
  f << "ASPECT_RATIO " << _Geo.getVoxelSize() << " " << _Geo.getVoxelSize()
    << " 1" << std::endl
    << std::flush;
  f << "POINT_DATA " << _Nx * _Ny << std::endl << std::flush;

  // write pops
  for (int k = 0; k < LatSet::q; ++k) {
    std::stringstream rho;
    std::string popname = "pop" + std::to_string(k);
    if constexpr (std::is_same<T, double>::value) {
      rho << "SCALARS " << popname << " double" << std::endl << std::flush;
    } else {
      rho << "SCALARS " << popname << " float" << std::endl << std::flush;
    }
    rho << "LOOKUP_TABLE default" << std::endl << std::flush;
    for (int i = 0; i < _Nx * _Ny; ++i) {
      if (_GlobalIdx[i] != -1) {
        const PopulationNbr<T, LatSet> &pop = _Pop[_GlobalIdx[i]];
        rho << pop.f[k] << " ";
      } else {
        rho << 0 << " ";
      }
    }
    rho << std::endl;
    f << rho.str();
  }

  f.close();
}

template <typename T, typename LatSet>
void Genericlbm2D<T, LatSet>::WriteStruPoints(int step) const {
  int _Nx = _Geo.getNx();
  int _Ny = _Geo.getNy();
  const std::vector<int> &_GlobalIdx = _Geo.getGlobalIdx();
  // write to file
  DirCreator::Create_Dir("./vtkoutput");
  std::string fullName = "./vtkoutput/" + _name + std::to_string(step) + ".vtk";
  std::ofstream f(fullName.c_str());
  f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
  f << "Voxels" << std::endl << std::flush;
  f << "ASCII" << std::endl << std::flush;

  f << "DATASET STRUCTURED_POINTS" << std::endl << std::flush;
  f << "DIMENSIONS " << _Nx << " " << _Ny << " 1" << std::endl << std::flush;
  f << "ORIGIN " << _Geo.getMin()[0] + _Geo.getVoxelSize() * T(0.5) << " "
    << _Geo.getMin()[1] + _Geo.getVoxelSize() * T(0.5) << " 0" << std::endl
    << std::flush;
  f << "ASPECT_RATIO " << _Geo.getVoxelSize() << " " << _Geo.getVoxelSize()
    << " 1" << std::endl
    << std::flush;
  f << "POINT_DATA " << _Nx * _Ny << std::endl << std::flush;

  std::stringstream rho;
  std::stringstream velocity;
  if constexpr (std::is_same<T, double>::value) {
    rho << "SCALARS rho double" << std::endl << std::flush;
    velocity << "VECTORS velocity double" << std::endl << std::flush;
  } else {
    rho << "SCALARS rho float" << std::endl << std::flush;
    velocity << "VECTORS velocity float" << std::endl << std::flush;
  }
  rho << "LOOKUP_TABLE default" << std::endl << std::flush;
  // velocity << "LOOKUP_TABLE default" << std::endl << std::flush;

  for (int i = 0; i < _Nx * _Ny; ++i) {
    if (_GlobalIdx[i] != -1) {
      const PopulationNbr<T, LatSet> &pop = _Pop[_GlobalIdx[i]];
      rho << pop.rho << " ";
      velocity << _Conv.getPhysU(pop.getVelocityx()) << " "
               << _Conv.getPhysU(pop.getVelocityy()) << " " << 0 << " ";
    } else {
      rho << 0 << " ";
      velocity << 0 << " " << 0 << " " << 0 << " ";
    }
  }
  rho << std::endl;
  velocity << std::endl;
  f << rho.str() << velocity.str();
  f.close();
}

//=============================================
// legacy lbm2D class
//=============================================
///////////////////////

template <typename T, typename LatSet, typename LatStru>
lbm2D<T, LatSet, LatStru>::lbm2D(
    Velocity2D<T> &field, AbstractConverter<T> &Conv_,
    BoundaryManager2D<T, LatSet, LatStru> &boundaryman, std::string name_,
    std::string rhoname_)
    : lbm2D(field, Conv_, boundaryman, name_, rhoname_, popIdx, popInIdx) {}

// if delegated constructor is used, Do Not use initialization list to
// initialize other members

template <typename T, typename LatSet, typename LatStru>
lbm2D<T, LatSet, LatStru>::lbm2D(
    Velocity2D<T> &field, AbstractConverter<T> &Conv_,
    BoundaryManager2D<T, LatSet, LatStru> &boundaryman, std::string name_,
    std::string rhoname_, std::vector<int> &cells, std::vector<int> &incells)
    : Lattice_Rho_Init(Conv_.getLatRhoInit()),
      Field(field),
      Geo(field.Geo),
      Conv(Conv_),
      BDM(boundaryman),
      Ni(field.Geo.getNi()),
      Nj(field.Geo.getNj()),
      N(field.Geo.getN()),
      name(name_),
      namerho(rhoname_),
      LatStru(field.Geo.getNi(), field.Geo.getNj()),
      OMEGA(Conv_.GetOMEGA()),
      _OMEGA(T(1) - Conv_.GetOMEGA()),
      Lattice_gbeta(Conv_.GetLattice_gbeta()),
      Cells(cells),
      InCells(incells) {
  pop = new population<T, LatSet>[N];
  // Init
  Traverse_Bulk(Ni, Nj, 1, [this](int id) {
    if (Geo.getVoxel(id).getFlag() != -1) {
      pop[id] = population<T, LatSet>(Lattice_Rho_Init);
    }
  });

  Setup();
}

template <typename T, typename LatSet, typename LatStru>
inline void lbm2D<T, LatSet, LatStru>::Setup() {
  std::cout << "[" << name << ":]"
            << "\n";
  if (&Cells == &popIdx) {
    std::cout << "use built-in popIndex" << std::endl;
    popIdx.reserve(N);
    Traverse_Bulk(Ni, Nj, 1, [this](int id) {
      if (Geo.getVoxel(id).getFlag() != -1) popIdx.emplace_back(id);
    });
  } else {
    std::cout << "use index from CellComm" << std::endl;
  }
  if (&InCells == &popInIdx) {
    std::cout << "use built-in popinIndex" << std::endl;
    popInIdx.reserve(N);
    Traverse_Bulk(Ni, Nj, 1, [this](int id) {
      if (Geo.getVoxel(id).getFlag() == 0) popInIdx.emplace_back(id);
    });
  } else {
    std::cout << "use index from CellComm" << std::endl;
  }
  // print index statistics
  std::cout << "CellIndex: " << Cells.size() << "\n"
            << "InCellIndex: " << InCells.size() << std::endl;
}

template <typename T, typename LatSet, typename LatStru>
lbm2D<T, LatSet, LatStru>::~lbm2D() {
  delete[] pop;
  if (force_enabled) {
    delete[] Force;
  }
}
template <typename T, typename LatSet, typename LatStru>
void lbm2D<T, LatSet, LatStru>::EnableForce() {
  Force = new Force2D<T>[N];
  force_enabled = true;
}
template <typename T, typename LatSet, typename LatStru>
void lbm2D<T, LatSet, LatStru>::AddtoInIdx(int flag) {
  Traverse_Bulk(Ni, Nj, 1, [this, &flag](int id) {
    if (Geo.getVoxel(id).getFlag() == flag) InCells.emplace_back(id);
  });
}

/*-------------Collision and Streaming----------------*/
// bgk collision:
// fp = f + (feq - f) * OMEGA = OMEGA * feq + _OMEGA * f
// where OMEGA = 1 / Tau, Tau > 1/2 -> 0 < OMEGA < 2
// Tau > 1, 0 < OMEGA < 1, f decays exponentially towards feq
// Tau < 1, 1 < OMEGA < 2, f oscillates around feq with an exponentially
// decreasing amplitude
// Tau = 1, OMEGA = 1, f decays directly towards feq
template <typename T, typename LatSet, typename LatStru>
template <void (*get_feq)(T *, const T *, T)>
void lbm2D<T, LatSet, LatStru>::Collide_BGK(const std::vector<int> &idx) {
  int Id = 0;
  T feq[LatSet::q] = {T(0)};
#pragma omp parallel for private(Id, feq) collapse(1) num_threads(Thread_Num)
  for (int i = 0; i < idx.size(); ++i) {
    Id = idx[i];
    get_feq(feq, Field.U[Id], pop[Id].rho);
    for (int k = 0; k < LatSet::q; ++k) {
      // pop[Id].fpostcol[k] = pop[Id].f[k] + (feq[k] - pop[Id].f[k]) * OMEGA;
      pop[Id].fpostcol[k] = OMEGA * feq[k] + _OMEGA * pop[Id].f[k];
    }
  }
}

// template <typename T, typename LatSet, typename LatStru>
// template <void (*get_feq)(T *, const T *, T)>
// void lbm2D<T, LatSet, LatStru>::Collide_BGK_p(
//     const std::vector<std::vector<int>> &idx_th) {
//   int Id = 0;
//   T feq[LatSet::q] = {T(0)};
//   int thn = 0;
// #pragma omp parallel private(Id, feq, thn) num_threads(Thread_Num)
//   {
//     thn = omp_get_thread_num();
//     for (int i = 0; i < idx_th[thn].size(); i++) {
//       Id = idx_th[thn][i];
//       get_feq(feq, Field.U[Id], pop[Id].rho);
//       for (int k = 0; k < LatSet::q; k++) {
//         // pop[Id].fpostcol[k] = pop[Id].f[k] + (feq[k] - pop[Id].f[k]) *
//         OMEGA; pop[Id].fpostcol[k] = OMEGA * feq[k] + _OMEGA * pop[Id].f[k];
//       }
//     }
//   }
// }

template <typename T, typename LatSet, typename LatStru>
template <void (*get_feq)(T *, const T *, T)>
inline void lbm2D<T, LatSet, LatStru>::Collide_BGK_Partial(
    const std::vector<direction<LatSet::q>> &BdCell_Dirs, T *fraction) {
  int Id = 0;
  int Dir = 0;
  T feq[LatSet::q] = {T(0)};
#pragma omp parallel for private(Id, Dir, feq) collapse(1) \
    num_threads(Thread_Num)
  for (int i = 0; i < BdCell_Dirs.size(); i++) {
    Id = BdCell_Dirs[i].Id;
    get_feq(feq, Field.U[Id], pop[Id].rho);
    for (int k = 0; k < BdCell_Dirs[i].inflow.size(); k++) {
      Dir = LatSet::opp[BdCell_Dirs[i].inflow[k]];
      pop[Id].fpostcol[Dir] = pop[Id].f[Dir] + (feq[Dir] - pop[Id].f[Dir]) *
                                                   OMEGA * (1 - fraction[Id]);
      // pop[Id].fpostcol[k] = OMEGA * feq[k] + _OMEGA * pop[Id].f[k];
    }
    // for (int k = 0; k < LatSet::q; k++) {
    //   pop[Id].fpostcol[k] = pop[Id].f[k] + (feq[k] - pop[Id].f[k]) * OMEGA *
    //   (1 - fraction[Id]);
    // }
  }
}

template <typename T, typename LatSet, typename LatStru>
void lbm2D<T, LatSet, LatStru>::Stream(const std::vector<int> &idx) {
  // boundary cells can be handled by direction class (inflow = stream)
  int Id = 0;
// inner cell streaming
#pragma omp parallel for private(Id) num_threads(Thread_Num)
  for (int i = 0; i < idx.size(); ++i) {
    Id = idx[i];
    pop[Id].f[0] = pop[Id].fpostcol[0];  // IMPROTANT!!!
    for (int k = 1; k < LatSet::q; ++k) {
      pop[Id].f[k] = pop[Id + LatStru::Nbr[LatSet::opp[k]]].fpostcol[k];
      // pop[Id - Lattice::Nbr[direction]].fpostcol[direction];
    }
  }
  // boundary cell streaming
  BDM.Stream(pop);
}

/*-----------------Boundary Conditions-----------------*/
template <typename T, typename LatSet, typename LatStru>
void lbm2D<T, LatSet, LatStru>::BCsManager() {
  // iterate over std::vector<BasicBoundary*> Bds
  for (int i = 0; i < BDM.Bds.size(); i++) {
    BDM.Bds[i]->Apply(pop);
  }
}

/*---------------post processing----------------*/
template <typename T, typename LatSet, typename LatStru>
template <void (lbm2D<T, LatSet, LatStru>::*rho_method)(),
          void (lbm2D<T, LatSet, LatStru>::*u_method)()>
void lbm2D<T, LatSet, LatStru>::PostStreamrhoU() {
  (this->*rho_method)();
  (this->*u_method)();
}

template <typename T, typename LatSet, typename LatStru>
inline void lbm2D<T, LatSet, LatStru>::compute_rho(
    const std::vector<int> &idx) {
  int i = 0;
#pragma omp parallel for private(i) num_threads(Thread_Num)
  for (i = 0; i < idx.size(); i++) {
    pop[idx[i]].Compute_rho();
  }
}

template <typename T, typename LatSet, typename LatStru>
inline void lbm2D<T, LatSet, LatStru>::compute_rho(int id) {
  pop[id].Compute_rho();
}

template <typename T, typename LatSet, typename LatStru>
inline void lbm2D<T, LatSet, LatStru>::compute_u(const std::vector<int> &idx) {
  int i = 0;
  int Id = 0;
#pragma omp parallel for private(i, Id) num_threads(Thread_Num)
  for (i = 0; i < idx.size(); i++) {
    Id = idx[i];
    pop[Id].Compute_U(Field.U[Id]);
  }
}

template <typename T, typename LatSet, typename LatStru>
inline void lbm2D<T, LatSet, LatStru>::poststream_rho() {
  compute_rho(Cells);
}
template <typename T, typename LatSet, typename LatStru>
inline void lbm2D<T, LatSet, LatStru>::poststream_innerrho() {
  compute_rho(InCells);
}
template <typename T, typename LatSet, typename LatStru>
inline void lbm2D<T, LatSet, LatStru>::poststream_u() {
  compute_u(Cells);
}
template <typename T, typename LatSet, typename LatStru>
inline void lbm2D<T, LatSet, LatStru>::poststream_inneru() {
  compute_u(InCells);
}

// force
template <typename T, typename LatSet, typename LatStru>
void lbm2D<T, LatSet, LatStru>::addtoBuoyancy(std::vector<int> &popidx,
                                              Force2D<T> *force) {
  int Id = 0;
#pragma omp parallel for private(Id) num_threads(Thread_Num)
  for (int i = 0; i < popidx.size(); i++) {
    Id = popidx.at(i);
    force[Id].F[1] += (pop[Id].rho - Lattice_Rho_Init) * Lattice_gbeta;
  }
}
template <typename T, typename LatSet, typename LatStru>
void lbm2D<T, LatSet, LatStru>::resetBuoyancy() {
  int Id = 0;
#pragma omp parallel for private(Id) num_threads(Thread_Num)
  for (Id = 0; Id < N; Id++) {
    Force[Id].Reset();
  }
}
template <typename T, typename LatSet, typename LatStru>
inline void lbm2D<T, LatSet, LatStru>::addtoRho(int id, T value) {
  pop[id].addtoRho(value);
}
template <typename T, typename LatSet, typename LatStru>
inline void lbm2D<T, LatSet, LatStru>::addtoPop(int id, T value) {
  pop[id].addtoPop(value);
}

template <typename T, typename LatSet, typename LatStru>
inline T lbm2D<T, LatSet, LatStru>::getAverRho() {
  T averrho = T(0);
#pragma omp parallel for reduction(+ : averrho) num_threads(Thread_Num)
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      averrho += pop[Index2D::GetId(i, j, Ni)].rho;
    }
  }
  return averrho / T(Ni - 2) / T(Nj - 2);
}

template <typename T, typename LatSet, typename LatStru>
template <T (*get_force)(int, T *, T *, T)>
void lbm2D<T, LatSet, LatStru>::AddBuoyancytoPop(std::vector<int> &idx) {
  int Id = 0;
#pragma omp parallel for private(Id) num_threads(Thread_Num)
  for (int i = 0; i < idx.size(); i++) {
    Id = idx[i];
    for (int k = 1; k < LatSet::q; k++) {
      pop[Id].fpostcol[k] += get_force(k, Force[Id].F, Field.U[Id], OMEGA);
    }
  }
}

template <typename T, typename LatSet, typename LatStru>
void lbm2D<T, LatSet, LatStru>::setpoprho_From_Vector(std::vector<int> &popidx,
                                                      T rho_) {
  int id;
#pragma omp parallel for private(id) num_threads(Thread_Num)
  for (int i = 0; i < popidx.size(); i++) {
    id = popidx.at(i);
    pop[id].rho = rho_;
  }
}

template <typename T, typename LatSet, typename LatStru>
void lbm2D<T, LatSet, LatStru>::setPopRho_From_VoxelFlag(int flag, T rho_) {
  int id;
#pragma omp parallel for private(id) num_threads(Thread_Num)
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      id = Index2D::GetId(i, j, Ni);
      if (Geo.getVoxel(id).getFlag() == flag) {
        pop[id].rho = rho_;
      }
    }
  }
}
