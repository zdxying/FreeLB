#include "legacy/legacy_lbm/lbm3D.h"

template <typename T, typename LatSet>
lbm3D<T, LatSet>::lbm3D(VelocityField3D<T> &velocity,
                        AbstractConverter<T> &conv,
                        BoundaryManager3D<T, LatSet> &bdmanager,
                        std::string name, std::string rhoname,
                        std::vector<int> *state,
                        BasicBoundary3D<T, LatSet> *movingbd, bool bdsetup)
    : Equilibrium3D<T, LatSet>(conv.GetOMEGA()),
      _Velocity(velocity),
      _Geo(velocity.getGeo()),
      _Conv(conv),
      _BdManager(bdmanager),
      _name(name),
      _namerho(rhoname),
      _OMEGA(conv.GetOMEGA()),
      _1_OMEGA(T(1) - conv.GetOMEGA()),
      _Lattice_Rho_Init(conv.getLatRhoInit()),
      _Lattice_gbeta(conv.GetLattice_gbeta()),
      State(state),
      MovingBd(movingbd) {
  // init _Pop based on voxel geometry and setup neighbor
  _Pop.clear();
  _Pop.reserve(_Geo.getVoxels().size());
  for (int i = 0; i < _Geo.getVoxels().size(); ++i) {
    _Pop.emplace_back(_Geo.getVoxel(i), _Velocity.getVelocity(i),
                      _Lattice_Rho_Init);
    // set PopInIdx
    if (_Geo.getVoxel(i).getFlag() == 0) PopInIdx.push_back(i);
  }
  // set pop neighbor
  for (int i = 0; i < _Geo.getVoxels().size(); ++i) {
    const std::vector<Voxel<T, LatSet::d> *> &voxnbrs =
        _Geo.getVoxel(i).getNeighborList();
    for (int k = 1; k < LatSet::q; ++k) {
      // NOTE THAT:
      // add INVERSE direction(inflow source) to facilitate streaming
      // e.g., direction 1 {1, 0, 0} (right hand side) != nullptr
      //       set to direction 3 {-1, 0, 0} (left hand side)
      if (voxnbrs[k] != nullptr) {
        if (voxnbrs[k]->getFlag() != -1)
          _Pop[i].setNeighbor(LatSet::opp[k], &_Pop[voxnbrs[k]->getId()]);
      }
    }
  }
  // set post-streaming index
  PostRhoIdx.clear();
  PostUIdx.clear();
  PostRhoIdx.insert(PostRhoIdx.end(), PopInIdx.begin(), PopInIdx.end());
  PostUIdx.insert(PostUIdx.end(), PopInIdx.begin(), PopInIdx.end());
  if (bdsetup == true) {
    // set bcs
    _BdManager.Setup(_Pop);
    _BdManager.printinfo();
  }
}

template <typename T, typename LatSet>
void lbm3D<T, LatSet>::setPopRho_From_VoxelFlag(int flag, T rho) {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (auto &pop : _Pop) {
    if (pop.getVoxel().getFlag() == flag) pop.rho = rho;
  }
}

template <typename T, typename LatSet>
void lbm3D<T, LatSet>::addtoPostRhoIdx(int flag) {
  for (int i = 0; i < _Pop.size(); ++i) {
    PopulationNbr<T, LatSet> &pop = _Pop[i];
    if (pop.getVoxel().getFlag() == flag) PostRhoIdx.push_back(i);
  }
}

template <typename T, typename LatSet>
void lbm3D<T, LatSet>::addtoPostUIdx(int flag) {
  for (int i = 0; i < _Pop.size(); ++i) {
    PopulationNbr<T, LatSet> &pop = _Pop[i];
    if (pop.getVoxel().getFlag() == flag) PostUIdx.push_back(i);
  }
}

/*-------------Collision and Streaming----------------*/
// BGK collision:
// fp = f + (feq - f) * OMEGA = OMEGA * feq + _OMEGA * f
// where OMEGA = 1 / Tau, Tau > 1/2 -> 0 < OMEGA < 2
// Tau > 1, 0 < OMEGA < 1, f decays exponentially towards feq
// Tau < 1, 1 < OMEGA < 2, f oscillates around feq with an exponentially
// decreasing amplitude
// Tau = 1, OMEGA = 1, f decays directly towards feq

// template <typename T, typename LatSet>
// template <void (*get_feq)(T *, const T *, const T)>
// void lbm3D<T, LatSet>::Collide_BGK(const std::vector<int> &idx) {
//   int Id = 0;
//   T feq[LatSet::q] = {T(0)};
// #pragma omp parallel for private(Id, feq) collapse(1) num_threads(Thread_Num)
//   for (int i = 0; i < idx.size(); i++) {
//     Id = idx[i];
//     get_feq(feq, Field.U[Id], pop[Id].rho);
//     for (int k = 0; k < LatSet::q; k++) {
//       _Pop[Id].fpostcol[k] = _OMEGA * feq[k] + _1_OMEGA * _Pop[Id].f[k];
//     }
//   }
// }
template <typename T, typename LatSet>
template <
    void (Equilibrium3D<T, LatSet>::*Collide)(AbstractPopulation<T, LatSet> &)>
void lbm3D<T, LatSet>::Collide_BGK() {
  int i = 0;
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (PopulationNbr<T, LatSet> &pop : _Pop) (this->*Collide)(pop);
}

template <typename T, typename LatSet>
template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T)>
void lbm3D<T, LatSet>::Collide_BGK() {
  int size = _Pop.size();
  T feq[LatSet::q] = {T(0)};
#pragma omp parallel for private(feq) num_threads(Thread_Num) schedule(static)
  for (PopulationNbr<T, LatSet> &pop : _Pop) {
    get_feq(feq, pop.getVelocity(), pop.rho);
    for (int k = 0; k < LatSet::q; ++k) {
      pop.fpostcol[k] = _OMEGA * feq[k] + _1_OMEGA * pop.f[k];
    }
  }
}

template <typename T, typename LatSet>
template <void (*get_feq)(T *, const Vector<T, LatSet::d> &, T)>
void lbm3D<T, LatSet>::Collide_BGK_StateCheck() {
  T feq[LatSet::q] = {T(0)};
#pragma omp parallel for private(feq) num_threads(Thread_Num) schedule(static)
  for (PopulationNbr<T, LatSet> &pop : _Pop) {
    if (State->operator[](pop.getVoxel().getId()) != -1) {
      get_feq(feq, pop.getVelocity(), pop.rho);
      for (int k = 0; k < LatSet::q; ++k) {
        pop.fpostcol[k] = _OMEGA * feq[k] + _1_OMEGA * pop.f[k];
      }
    }
  }
}

template <typename T, typename LatSet>
void lbm3D<T, LatSet>::Stream() {
  int i = 0;
  int size = PopInIdx.size();
// inner cell streaming
#pragma omp parallel for private(i) num_threads(Thread_Num) schedule(static)
  for (i = 0; i < size; ++i) {
    PopulationNbr<T, LatSet> &pop = _Pop[PopInIdx[i]];
    pop.f[0] = pop.fpostcol[0];  // IMPROTANT!!!
    for (int k = 1; k < LatSet::q; ++k) {
      pop.f[k] = pop.getNeighbor(k)->fpostcol[k];
    }
  }
  // boundary cell streaming
  _BdManager.Stream();
}

template <typename T, typename LatSet>
void lbm3D<T, LatSet>::Stream_StateCheck() {
// inner cell streaming
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (PopulationNbr<T, LatSet> &pop : _Pop) {
    if (State->operator[](pop.getVoxel().getId()) == 0) {
      pop.f[0] = pop.fpostcol[0];  // IMPROTANT!!!
      for (int k = 1; k < LatSet::q; ++k) {
        pop.f[k] = pop.getNeighbor(k)->fpostcol[k];
      }
    }
  }
  // boundary cell streaming
  _BdManager.Stream();
}

template <typename T, typename LatSet>
inline void lbm3D<T, LatSet>::BCs() {
  _BdManager.Apply();
}

template <typename T, typename LatSet>
void lbm3D<T, LatSet>::Compute_Rho() {
  int i = 0;
  int size = _Pop.size();
#pragma omp parallel for private(i) num_threads(Thread_Num) schedule(static)
  for (i = 0; i < size; ++i) {
    _Pop[i].Compute_rho();
  }
}
template <typename T, typename LatSet>
void lbm3D<T, LatSet>::Compute_Rho(std::vector<int> &RhoIdx) {
  int i = 0;
  int size = RhoIdx.size();
#pragma omp parallel for private(i) num_threads(Thread_Num) schedule(static)
  for (i = 0; i < size; ++i) {
    _Pop[RhoIdx[i]].Compute_rho();
  }
}
template <typename T, typename LatSet>
void lbm3D<T, LatSet>::Compute_Rho_Commun() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (PopulationNbr<T, LatSet> &pop : _Pop) {
    if (State->operator[](pop.getVoxel().getId()) == 0) {
      pop.Compute_rho();
    }
  }
}
template <typename T, typename LatSet>
void lbm3D<T, LatSet>::Compute_U() {
  int i = 0;
  int size = _Pop.size();
#pragma omp parallel for private(i) num_threads(Thread_Num) schedule(static)
  for (i = 0; i < size; ++i) {
    _Pop[i].Compute_U();
  }
}
template <typename T, typename LatSet>
void lbm3D<T, LatSet>::Compute_U(std::vector<int> &UIdx) {
  int i = 0;
  int size = UIdx.size();
#pragma omp parallel for private(i) num_threads(Thread_Num) schedule(static)
  for (i = 0; i < size; ++i) {
    _Pop[UIdx[i]].Compute_U();
  }
}
template <typename T, typename LatSet>
void lbm3D<T, LatSet>::Compute_U_Commun() {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
  for (PopulationNbr<T, LatSet> &pop : _Pop) {
    if (State->operator[](pop.getVoxel().getId()) == 0) {
      pop.Compute_U();
    }
  }
}

template <typename T, typename LatSet>
T lbm3D<T, LatSet>::getToleranceRho() {
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
T lbm3D<T, LatSet>::getToleranceU() {
  T res0, res1, res2, res;
  T maxres = T(0);
  for (int i = 0; i < _Pop.size(); ++i) {
    res0 = std::abs(_Pop[i].getVelocity()[0] - _UOld[i][0]);
    res1 = std::abs(_Pop[i].getVelocity()[1] - _UOld[i][1]);
    res2 = std::abs(_Pop[i].getVelocity()[2] - _UOld[i][2]);
    res = std::max(res0, std::max(res1, res2));
    maxres = std::max(res, maxres);
    // set _UOld
    _UOld[i][0] = _Pop[i].getVelocity()[0];
    _UOld[i][1] = _Pop[i].getVelocity()[1];
    _UOld[i][2] = _Pop[i].getVelocity()[2];
  }
  return maxres;
}

// dataType is one of the types bit, unsigned_char, char, unsigned_short, short,
// unsigned_int, int, unsigned_long, long, float, or double
template <typename T, typename LatSet>
void lbm3D<T, LatSet>::WriteStruPoints(int step) const {
  int _Nx = _Geo.getNx();
  int _Ny = _Geo.getNy();
  int _Nz = _Geo.getNz();
  const std::vector<int> &_GlobalIdx = _Geo.getGlobalIdx();
  // write to file
  DirCreator::Create_Dir("./vtkoutput");
  std::string fullName = "./vtkoutput/" + _name + std::to_string(step) + ".vtk";
  std::ofstream f(fullName.c_str());
  f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
  f << "Voxels" << std::endl << std::flush;
  f << "ASCII" << std::endl << std::flush;

  f << "DATASET STRUCTURED_POINTS" << std::endl << std::flush;
  f << "DIMENSIONS " << _Nx << " " << _Ny << " " << _Nz << std::endl
    << std::flush;
  f << "ORIGIN " << _Geo.getMin()[0] + _Geo.getVoxelSize() * T(0.5) << " "
    << _Geo.getMin()[1] + _Geo.getVoxelSize() * T(0.5) << " "
    << _Geo.getMin()[2] + _Geo.getVoxelSize() * T(0.5) << std::endl
    << std::flush;
  f << "ASPECT_RATIO " << _Geo.getVoxelSize() << " " << _Geo.getVoxelSize()
    << " " << _Geo.getVoxelSize() << std::endl
    << std::flush;
  f << "POINT_DATA " << _Nx * _Ny * _Nz << std::endl << std::flush;

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

  for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
    if (_GlobalIdx[i] != -1) {
      const PopulationNbr<T, LatSet> &pop = _Pop[_GlobalIdx[i]];
      rho << pop.rho << " ";
      velocity << _Conv.getPhysU(pop.getVelocityx()) << " "
               << _Conv.getPhysU(pop.getVelocityy()) << " "
               << _Conv.getPhysU(pop.getVelocityz()) << " ";
      // velocity << _Conv.getPhysU(_Velocity.getVelocity(_GlobalIdx[i])[0]) <<
      // " "
      //          << _Conv.getPhysU(_Velocity.getVelocity(_GlobalIdx[i])[1]) <<
      //          " "
      //          << _Conv.getPhysU(_Velocity.getVelocity(_GlobalIdx[i])[2]) <<
      //          " ";
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

template <typename T, typename LatSet>
void lbm3D<T, LatSet>::WriteStruPoints_(int step) const {
  int _Nx = _Geo.getNx();
  int _Ny = _Geo.getNy();
  int _Nz = _Geo.getNz();
  const std::vector<int> &_GlobalIdx = _Geo.getGlobalIdx();
  // write to file
  DirCreator::Create_Dir("./vtkoutput");
  std::string fullName = "./vtkoutput/" + _name + std::to_string(step) + ".vtk";
  std::ofstream f(fullName.c_str());
  f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
  f << "Voxels" << std::endl << std::flush;
  f << "ASCII" << std::endl << std::flush;

  f << "DATASET STRUCTURED_POINTS" << std::endl << std::flush;
  f << "DIMENSIONS " << _Nx << " " << _Ny << " " << _Nz << std::endl
    << std::flush;
  f << "ORIGIN " << _Geo.getMin()[0] + _Geo.getVoxelSize() * T(0.5) << " "
    << _Geo.getMin()[1] + _Geo.getVoxelSize() * T(0.5) << " "
    << _Geo.getMin()[2] + _Geo.getVoxelSize() * T(0.5) << std::endl
    << std::flush;
  f << "ASPECT_RATIO " << _Geo.getVoxelSize() << " " << _Geo.getVoxelSize()
    << " " << _Geo.getVoxelSize() << std::endl
    << std::flush;
  f << "POINT_DATA " << _Nx * _Ny * _Nz << std::endl << std::flush;

  if constexpr (std::is_same<T, double>::value) {
    f << "SCALARS rho double" << std::endl << std::flush;
  } else {
    f << "SCALARS rho float" << std::endl << std::flush;
  }
  f << "LOOKUP_TABLE default" << std::endl << std::flush;
  for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
    if (_GlobalIdx[i] != -1) {
      const PopulationNbr<T, LatSet> &pop = _Pop[_GlobalIdx[i]];
      f << pop.rho << " ";
    } else {
      f << 0 << " ";
    }
  }
  f << std::endl;
  if constexpr (std::is_same<T, double>::value) {
    f << "VECTORS velocity double" << std::endl << std::flush;
  } else {
    f << "VECTORS velocity float" << std::endl << std::flush;
  }
  for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
    if (_GlobalIdx[i] != -1) {
      const PopulationNbr<T, LatSet> &pop = _Pop[_GlobalIdx[i]];
      f << _Conv.getPhysU(pop.getVelocityx()) << " "
        << _Conv.getPhysU(pop.getVelocityy()) << " "
        << _Conv.getPhysU(pop.getVelocityz()) << " ";
    } else {
      f << 0 << " " << 0 << " " << 0 << " ";
    }
  }
  f << std::endl;
  f.close();
}

template <typename T, typename LatSet>
void lbm3D<T, LatSet>::WriteRho(int step) const {
  int _Nx = _Geo.getNx();
  int _Ny = _Geo.getNy();
  int _Nz = _Geo.getNz();
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
  f << "DIMENSIONS " << _Nx << " " << _Ny << " " << _Nz << std::endl
    << std::flush;
  f << "ORIGIN " << _Geo.getMin()[0] + _Geo.getVoxelSize() * T(0.5) << " "
    << _Geo.getMin()[1] + _Geo.getVoxelSize() * T(0.5) << " "
    << _Geo.getMin()[2] + _Geo.getVoxelSize() * T(0.5) << std::endl
    << std::flush;
  f << "ASPECT_RATIO " << _Geo.getVoxelSize() << " " << _Geo.getVoxelSize()
    << " " << _Geo.getVoxelSize() << std::endl
    << std::flush;
  f << "POINT_DATA " << _Nx * _Ny * _Nz << std::endl << std::flush;

  std::stringstream rho;
  if constexpr (std::is_same<T, double>::value) {
    rho << "SCALARS rho double" << std::endl << std::flush;
  } else {
    rho << "SCALARS rho float" << std::endl << std::flush;
  }
  rho << "LOOKUP_TABLE default" << std::endl << std::flush;
  for (int i = 0; i < _Nx * _Ny * _Nz; ++i) {
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
