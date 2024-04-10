
#include "legacy/legacy_io/legacy_writer.h"

template <typename T>
void DataWriter2D<T>::Write_Geometry(std::string name) {
  std::ofstream geo;
  geo.open(work_dir + name + ".dat");
  Write_datHeader(geo, XY, std::string("Flag"), 0);  // set offset to 0
  int id;
  for (int j = 0; j < Nj; j++) {
    for (int i = 0; i < Ni; i++) {
      cord_Lat(geo, i, j);
      id = GetId(i, j, Ni);
      Write_LatFlag(geo, id);
    }
  }
  geo.close();
}

template <typename T>
inline void DataWriter2D<T>::Wrtie_FromFlag(int step, std::vector<int>& CellIds,
                                            std::string name) {
  int N = Ni * Nj;
  int* flag_field = new int[N];
  std::fill_n(flag_field, N, 0);
  // write
  for (int i = 0; i < CellIds.size(); i++) {
    flag_field[CellIds[i]] = 1;
  }
  // write to file
  std::ofstream Flags;
  Flags.open(work_dir + name + std::to_string(step) + ".dat");
  Write_datHeader(Flags, XY, std::string("Flag"), 0);
  int id;
  for (int j = 0; j < Nj; j++) {
    for (int i = 0; i < Ni; i++) {
      cord_Lat(Flags, i, j);
      id = GetId(i, j, Ni);
      Write_Flag(Flags, id, flag_field);
    }
  }
  Flags.close();
  delete[] flag_field;
}

template <typename T>
inline void DataWriter2D<T>::Write_MovingBCs(int step,
                                             std::vector<int>& CellIds,
                                             std::string name) {
  // overwriting flag
  bool overwriting_ = false;
  std::vector<int> overwriting_ids;
  // prepare int *isbound_
  int N = Ni * Nj;
  int* isbound_ = new int[N];
  std::fill_n(isbound_, N, 0);
  // write to isbound_
  // boundary
  for (int i = 0; i < CellIds.size(); i++) {
    isbound_[CellIds[i]] = 1;
  }
  // inavtice cells
  for (int i = 0; i < N; i++) {
    if (field.Geo.getVoxel(i).getFlag() == -1) {
      if (isbound_[i] == 0) {
        isbound_[i] = -1;
      } else {
        // 2 means overwriting
        isbound_[i] = 2;
        overwriting_ = true;
        overwriting_ids.emplace_back(i);
      }
    }
  }
  // write to file
  std::ofstream BCs;
  BCs.open(work_dir + name + std::to_string(step) + ".dat");
  Write_datHeader(BCs, XY, std::string("isbound"), 0);
  int id;
  for (int j = 0; j < Nj; j++) {
    for (int i = 0; i < Ni; i++) {
      cord_Lat(BCs, i, j);
      id = GetId(i, j, Ni);
      Write_Flag(BCs, id, isbound_);
    }
  }
  BCs.close();
  delete[] isbound_;

  // overwriting log
  if (overwriting_) {
    std::cout << "Warning: overwriting" << overwriting_ids.size()
              << " inactive cells with moving boundary cells" << std::endl;
    // std::ofstream overwriting_log;
    // overwriting_log.open(work_dir + name + std::to_string(step) +
    //                      "_overwriting.log");
    // for (int i = 0; i < overwriting_ids.size(); i++) {
    //   overwriting_log << overwriting_ids[i] << "\n";
    // }
    // overwriting_log.close();
  }
}

template <typename T>
template <typename LatSet>
inline void DataWriter2D<T>::Write_BCsDirection(
    int step, std::vector<direction<LatSet::q>>& Cell_Dirs, std::string name) {
  // prepare int *inflow_ and *outflow_
  int N = Ni * Nj;
  int* inflow_ = new int[N];
  int* outflow_ = new int[N];
  std::fill_n(inflow_, N, 0);
  std::fill_n(outflow_, N, 0);
  // write to inflow_ and outflow_
  for (int i = 0; i < Cell_Dirs.size(); i++) {
    int id = Cell_Dirs[i].Id;
    inflow_[id] = Cell_Dirs[i].inflow.size();
    outflow_[id] = Cell_Dirs[i].outflow.size();
  }
  // write to file
  std::ofstream BCs;
  BCs.open(work_dir + name + std::to_string(step) + ".dat");
  Write_datHeader(BCs, XY, std::string("inflow,outflow"), 0);
  int id;
  for (int j = 0; j < Nj; j++) {
    for (int i = 0; i < Ni; i++) {
      cord_Lat(BCs, i, j);
      id = GetId(i, j, Ni);
      Write_Flag(BCs, id, inflow_);
      Write_Flag(BCs, id, outflow_);
    }
  }
  BCs.close();
  delete[] inflow_;
  delete[] outflow_;
}

template <typename T>
template <typename LatSet>
inline void DataWriter2D<T>::Write_Populations(int step,
                                               population<T, LatSet>* pop) {
  std::ofstream popfile;
  popfile.open(work_dir + "/Populations_" + std::to_string(step) + ".dat");
  // generate header, add to name, like: f0, f1, f2, ...
  std::string name = "f" + std::to_string(0);
  for (int i = 1; i < LatSet::q; i++) {
    name += ",f" + std::to_string(i);
  }
  Write_datHeader(popfile, XY, name);
  int id;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      cord_Lat(popfile, i, j);
      id = GetId(i, j, Ni);
      for (int k = 0; k < LatSet::q; k++) {
        write_data(popfile, pop[id].f[k]);
      }
    }
  }
  popfile.close();
}

template <typename T>
template <typename LatSet>
inline void DataWriter2D<T>::Write_Populations(int step, int k,
                                               population<T, LatSet>* pop) {
  std::ofstream popfile;
  popfile.open(work_dir + "/Populations_" + std::to_string(step) + ".dat");
  std::string name = "f" + std::to_string(k);
  Write_datHeader(popfile, XY, name);
  int id;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      cord_Lat(popfile, i, j);
      id = GetId(i, j, Ni);
      write_data(popfile, pop[id].f[k]);
    }
  }
  popfile.close();
}

template <typename T>
template <typename LatSet>
inline void DataWriter2D<T>::Write_Populations(int step, std::vector<int>& dir,
                                               population<T, LatSet>* pop) {
  int size = dir.size();
  if (size == 1) {
    Write_Populations(step, dir[0], pop);
    return;
  } else {
    std::ofstream popfile;
    popfile.open(work_dir + "/Populations_" + std::to_string(step) + ".dat");
    // generate header, add to name, like: f0, f1, f2, ...
    std::string name = "f" + std::to_string(dir[0]);
    for (int i = 1; i < size; i++) {
      name += ",f" + std::to_string(dir[i]);
    }
    Write_datHeader(popfile, XY, name);
    int id;
    for (int j = 1; j < Nj - 1; j++) {
      for (int i = 1; i < Ni - 1; i++) {
        cord_Lat(popfile, i, j);
        id = GetId(i, j, Ni);
        for (int k = 0; k < size; k++) {
          write_data(popfile, pop[id].f[dir[k]]);
        }
      }
    }
    popfile.close();
  }
}

template <typename T>
void DataWriter2D<T>::Write_Fluid_Phys_Only(int step, std::string name) {
  std::ofstream fluid;
  fluid.open(work_dir + name + std::to_string(step) + ".dat");
  Write_datHeader(fluid, XY, Us);
  int id;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      cord_Lat(fluid, i, j);
      id = GetId(i, j, Ni);
      Write_U_Phys(fluid, id);
    }
  }
  fluid.close();
}

template <typename T>
void DataWriter2D<T>::Write_Rho_Phys_Only(int step, std::string name) {
  std::ofstream rho;
  rho.open(work_dir + name + std::to_string(step) + ".dat");
  Write_datHeader(rho, XY, getnamerho());
  int id;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      cord_Lat(rho, i, j);
      id = GetId(i, j, Ni);
      GeneralRho_Phys(rho, id);
    }
  }
  rho.close();
}

template <typename T>
inline void DataWriter2D<T>::Write_ZSCA_Only(int step, std::string name) {
  std::ofstream rho;
  rho.open(work_dir + name + std::to_string(step) + ".dat");
  Write_datHeader(rho, XY, std::string(",Fs,State"));
  int id;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      cord_Lat(rho, i, j);
      id = GetId(i, j, Ni);
      Write_SolidFrac(rho, id);
      Write_CAState(rho, id);
    }
  }
  rho.close();
}

template <typename T>
inline void DataWriter2D<T>::Write_ZSCA(int step, CAZSField2D<T>& zsca,
                                        std::string name) {
  std::ofstream rho;
  rho.open(work_dir + name + std::to_string(step) + ".dat");
  Write_datHeader(rho, XY, std::string(",Ceq,Ceq_Cl,aniso,Curv"));
  int id;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      cord_Lat(rho, i, j);
      id = GetId(i, j, Ni);
      // Write_SolidFrac(rho, id);
      // Write_CAState(rho, id);
      // Write_Curvature(rho, id);
#ifdef _FLB_DEBUG
      write_data(rho, zsca.Ceq[id]);
      write_data(rho, zsca.Ceq_Cl[id]);
      write_data(rho, zsca.aniso[id]);
      write_data(rho, zsca.K[id]);
#endif
    }
  }
  rho.close();
}

template <typename T>
void DataWriter2D<T>::Write_GCA(int step, CAGField2D<T>& gca,
                                std::string name) {
  std::ofstream ca;
  ca.open(work_dir + name + std::to_string(step) + ".dat");
  Write_datHeader(ca, XY, std::string(",Orine,Fs,State"));
  int id;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      cord_Lat(ca, i, j);
      id = GetId(i, j, Ni);
      write_data(ca, gca.Orine[id]);
      write_data(ca, gca.f[id]);
      Write_CAState(ca, id);
    }
  }
  ca.close();
}

template <typename T>
void DataWriter2D<T>::Write_Fluid_Rho_Phys(int step, std::string name) {
  std::ofstream rho;
  rho.open(work_dir + name + std::to_string(step) + ".dat");
  Write_datHeader(rho, XY, Us + getnamerho());
  int id;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      cord_Lat(rho, i, j);
      id = GetId(i, j, Ni);
      Write_U_Phys(rho, id);
      GeneralRho_Phys(rho, id);
    }
  }
  rho.close();
}

template <typename T>
inline void DataWriter2D<T>::Write_Fluid_Rho_Phys_And_ZSCA(int step,
                                                           std::string name) {
  std::ofstream rho;
  rho.open(work_dir + name + std::to_string(step) + ".dat");
  Write_datHeader(rho, XY, Us + getnamerho() + std::string(",Fs"));
  int id;
  for (int j = 1; j < Nj - 1; j++) {
    for (int i = 1; i < Ni - 1; i++) {
      cord_Lat(rho, i, j);
      id = GetId(i, j, Ni);
      Write_U_Phys(rho, id);
      GeneralRho_Phys(rho, id);
      Write_SolidFrac(rho, id);
    }
  }
  rho.close();
}

// void DataWriter2D::Get_BoundaryType(Boundary &bound)
// {
//     // boundary type:
//     // 0: no boundary; 1: bounceback; 2: anti-bounceback 3: bounceback with
//     moving wall 4:periodic 5:Generalized periodic 6:Non equilibrium
//     7:Symmetry
//     // create a new array to store boundary type
//     delBoundType = true;
//     int N = Ni * Nj;
//     BoundType = new int[N];
//     for (int i = 0; i < N; i++)
//     {
//         BoundType[i] = 0;
//     }
//     // write boundary type
//     for (int i = 0; i < bound.BB.size(); i++)
//     {
//         BoundType[bound.BB.at(i)] = 1;
//     }
//     for (int i = 0; i < bound.ABB.size(); i++)
//     {
//         BoundType[bound.ABB.at(i)] = 2;
//     }
//     for (int i = 0; i < bound.BBMW.size(); i++)
//     {
//         BoundType[bound.BBMW.at(i)] = 3;
//     }
//     for (int i = 0; i < bound.P.size(); i++)
//     {
//         BoundType[bound.P.at(i)] = 4;
//     }
//     for (int i = 0; i < bound.GP.size(); i++)
//     {
//         BoundType[bound.GP.at(i)] = 5;
//     }
//     for (int i = 0; i < bound.NEE.size(); i++)
//     {
//         BoundType[bound.NEE.at(i)] = 6;
//     }
//     for (int i = 0; i < bound.Sym.size(); i++)
//     {
//         BoundType[bound.Sym.at(i)] = 7;
//     }
// }

// void DataWriter2D::Write_CA(int step)
// {
//     std::ofstream CA;
//     CA.open(work_dir + "/CA_" + std::to_string(step) + ".dat");
//     Write_datHeader(CA, std::string("X,Y"), std::string("Orine"));
//     for (int j = 1; j < Nj - 1; j++)
//     {
//         for (int i = 1; i < Ni - 1; i++)
//         {
//             int Id = i + Ni * j;
//             Write_cord_Lat(CA, i, j);
//             Write_Orinentation(CA, Id);
//             Write_CAStatus(CA, Id);
//             CA << "\n";
//         }
//     }
//     CA.close();
// }

// void DataWriter2D::Write_Flag()
// {
//     std::ofstream LB;
//     LB.open(work_dir + "/BCs" + ".dat");
//     Write_datHeader(LB, std::string("X,Y"), std::string("inflow,outflow"));
//     for (int j = 1; j < Nj - 1; j++)
//     {
//         for (int i = 1; i < Ni - 1; i++)
//         {
//             int Id = i + Ni * j;
//             // LB << i << "\t" << j << "\t" << sqrt(field.Field[Id].U[0] *
//             field.Field[Id].U[0] + field.Field[Id].U[1] *
//             field.Field[Id].U[1])
//             //    << "\t" << field.Field[Id].U[0] << "\t" <<
//             field.Field[Id].U[1] << "\n"; LB << i << "\t" << j << "\t" <<
//             lb->dir[Id].inflow.size() << "\t" << lb->dir[Id].outflow.size()
//             << std::endl;
//             //    << field.LB_Mesh_T[Id].rho << "\t"
//             //    << field.LB_Mesh_C[Id].rho << "\t"
//         }
//     }
//     LB.close();
// }
