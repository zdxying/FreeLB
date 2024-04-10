/*write to file manager class*/
#pragma once

#include <fstream>

#include "legacy/legacy_lbm/legacy_lbmanager2d.h"
#include "utils/directories.h"
#include "utils/timer.h"

template <typename T>
class DataWriter2D : public Index2D {
 public:
  /*extract data from Mesh*/
  Velocity2D<T> &field;
  LBManager2D<T> *lbman;
  StateField2D<T> *ca;
  UnitConvManager<T> &ConvMan;
  int *BoundType;
  bool delBoundType = false;
  std::string work_dir;

  /*Mesh*/
  int Ni;
  int Nj;

  // string
  std::string XY = "X,Y";
  std::string Us = "Ux,Uy,U";
  // return e.g. ",rho,T,C"
  std::string getnamerho() {
    // get e.g. ",rho,T,C"
    std::string var;
    std::vector<Basiclbm2D<T> *> lbms = *lbman->getlbms();
    for (int i = 0; i < lbms.size(); i++) {
      var += "," + lbms[i]->getnamerho();
    }
    return var;
  }

  DataWriter2D(std::string work_dir_, Velocity2D<T> &field_,
               UnitConvManager<T> &convman, LBManager2D<T> *lb_ = nullptr,
               StateField2D<T> *ca_ = nullptr)
      : work_dir(work_dir_),
        field(field_),
        ConvMan(convman),
        Ni(field.Geo.getNi()),
        Nj(field.Geo.getNj()),
        lbman(lb_),
        ca(ca_) {
    DirCreator::Create_Dir(work_dir_);
  }
  ~DataWriter2D() {
    if (delBoundType) {
      delete[] BoundType;
    }
  }

  // Write_datHeader(LBCA, std::string("X,Y"),
  // std::string("Ux,Uy,U,rho,f1,f2"));
  // call: (write, XY, "Ux,Uy,U,rho,f1,f2"), NO "," before Ux
  void Write_datHeader(
      std::ofstream &write, std::string var1 = "X,Y",
      std::string var2 = "Ux,Uy,U",
      int offset = 2)  // declare offset = 2 as default in defination
  {
    write << "TITLE = \"Data\"\n"
          << "variables = " << var1 << "," << var2 << "\n"
          << "Zone t=\"data\"\n"
          << "I=" << Ni - offset << ",J=" << Nj - offset << ",F=POINT"
          << std::endl;
  }
  void cord_Lat(std::ofstream &write, int i, int j) {
    write << i << "\t" << j << "\t";
  }
  void Write_LatFlag(std::ofstream &write, int Id) {
    write << field.Geo.getVoxel(Id).getFlag() << "\t";
  }

  void Write_Geometry(std::string name = "/Geometry");
  void Wrtie_FromFlag(int step, std::vector<int> &CellIds,
                      std::string name = "/Flags_");
  // write inactive and moving boundary cells, overwriting check is enabled
  void Write_MovingBCs(int step, std::vector<int> &CellIds,
                       std::string name = "/BCs_");
  template <typename LatSet>
  void Write_BCsDirection(int step,
                          std::vector<direction<LatSet::q>> &Cell_Dirs,
                          std::string name = "/BCsDir_");
  void Write_Fluid_Phys_Only(int step, std::string name = "/FluidPhys_");
  void Write_Rho_Phys_Only(int step, std::string name = "/RhoPhys_");
  void Write_ZSCA_Only(int step, std::string name = "/ZSCA_");
  void Write_ZSCA(int step, CAZSField2D<T> &zsca,
                  std::string name = "/ZSCADBG_");
  void Write_GCA(int step, CAGField2D<T> &gca, std::string name = "/GCA_");
  void Write_Fluid_Rho_Phys(int step, std::string name = "/FluidPhys_RhoPhys_");
  void Write_Fluid_Rho_Phys_And_ZSCA(
      int step, std::string name = "/FluidPhys_RhoPhys_ZSCA_");
  void Write_Fluid_Rho_Phys_And_GCA(int step);
  // debug
  template <typename LatSet>
  void Write_Populations(int step, population<T, LatSet> *pop);
  template <typename LatSet>
  void Write_Populations(int step, int k, population<T, LatSet> *pop);
  template <typename LatSet>
  void Write_Populations(int step, std::vector<int> &dir,
                         population<T, LatSet> *pop);

  void Write_U_Lat(std::ofstream &write, int Id) {
    componentU_Lat(write, Id);
    scalerU_Lat(write, Id);
  }
  void Write_U_Phys(std::ofstream &write, int Id) {
    componentU_Phys(write, Id);
    scalerU_Phys(write, Id);
  }

  /*-------------Rho Method---------------*/
  inline void rho_Lat(std::ofstream &write, int Id, Basiclbm2D<T> *lbm) {
    write << lbm->getPoprho(Id) << "\t";
  }
  inline void GeneralRho_Lat(std::ofstream &write, int Id) {
    std::vector<Basiclbm2D<T> *> lbms = *lbman->getlbms();
    for (int i = 0; i < lbms.size(); i++) {
      rho_Lat(write, Id, lbms[i]);
    }
  }
  inline void rho_Phys(std::ofstream &write, int Id, Basiclbm2D<T> *lbm) {
    write << lbm->getPhysrho(Id) << "\t";
  }
  inline void GeneralRho_Phys(std::ofstream &write, int Id) {
    std::vector<Basiclbm2D<T> *> lbms = *lbman->getlbms();
    for (int i = 0; i < lbms.size(); i++) {
      rho_Phys(write, Id, lbms[i]);
    }
  }

  /*-------------U Method---------------*/
  inline void scalerU_Lat(std::ofstream &write, int Id) {
    write << Vect2D<T>::norm(field.U[Id]) << "\t";
  }
  inline void scalerU_Phys(std::ofstream &write, int Id) {
    write << ConvMan.BaseConv->getPhysU(Vect2D<T>::norm(field.U[Id])) << "\t";
  }
  inline void componentU_Lat(std::ofstream &write, int Id) {
    write << field.U[Id][0] << "\t" << field.U[Id][1] << "\t";
  }
  inline void componentU_Phys(std::ofstream &write, int Id) {
    write << ConvMan.BaseConv->getPhysU(field.U[Id][0]) << "\t"
          << ConvMan.BaseConv->getPhysU(field.U[Id][1]) << "\t";
  }

  /*----------CA-------------*/
  // ca->getOrine(Id)
  inline void Write_Orinentation(std::ofstream &write, int Id) {
    write << ca->getOrine(Id) << "\t";
  }
  // ca->State[Id]
  inline void Write_CAState(std::ofstream &write, int Id) {
    write << ca->State[Id] << "\t";
  }
  // ca->getSolidFraction(Id)
  inline void Write_SolidFrac(std::ofstream &write, int Id) {
    write << ca->getSolidFraction(Id) << "\t";
  }

  /*----------Flags-------------*/
  inline void Write_Flag(std::ofstream &write, int Id, int *flag) {
    write << flag[Id] << "\t";
  }

  inline void Write_Step(std::ofstream &write, int step) {
    write << step << "\t";
  }

  // BCS
  inline void SetupBoundType() { BoundType = new int[Ni * Nj]; }
  // 0: no boundary; 1: bounceback; 2: anti-bounceback 3: bounceback with moving
  // wall 4: free slip 5:periodic 6: non-equilibrium-extrapolation
  inline void Write_BoundaryType(std::ofstream &write, int Id) {
    write << BoundType[Id] << "\n";
  }
  template <typename LatSet>
  inline void Write_outflow(std::ofstream &write,
                            direction<LatSet::q> &Cell_Dirs) {
    write << Cell_Dirs.outflow.size() << "\n";
  }
  template <typename LatSet>
  inline void Write_inflow(std::ofstream &write,
                           direction<LatSet::q> &Cell_Dirs) {
    write << Cell_Dirs.inflow.size() << "\n";
  }

  // debug
  // distribution function
  template <typename LatSet>
  inline void Write_distribution_function(std::ofstream &write, int Id, int k,
                                          population<T, LatSet> *pop) {
    write << pop[Id].f[k] << "\n";
  }

  // performance

  void Write_Performance(int step, T MLUPS);

  // experimental
  template <typename U = T>
  inline void write_data(std::ofstream &write_, const U &data) {
    write_ << data << "\n";
  }
};

// freelb plot
template <typename T>
class FLBplot {
 private:
  // x-axis: step
  // y-axis: data
  std::string name;
  std::string xvar;
  std::string yvar;
  std::string work_dir;

 public:
  FLBplot(std::string work_dir_, std::string name = "/FLBplot",
          std::string xvar = "step", std::string yvar = "data")
      : work_dir(work_dir_), name(name), xvar(xvar), yvar(yvar) {
    prepare_plot();
  }
  inline void prepare_plot() {
    std::ofstream write;
    write.open(work_dir + name + ".dat");
    write << xvar << "\t" << yvar << "\n";
  }
  // each time call this function, a step-data pair is written
  template <typename X = int, typename U = T>
  inline void write_plot(X xvalue, U yvalue) {
    std::ofstream write;
    write.open(work_dir + name + ".dat", std::ios::app);
    write << xvalue << "\t" << yvalue << "\n";
  }
  template <typename U = T>
  inline void write_plot(Counter &counter, U yvalue) {
    write_plot<int, U>(counter(), yvalue);
  }
};

// template <typename T>
// class DataWriter3D : {
//  public:

// };

namespace wrt {

// call: (Ni, Nj, Nk, write, "X,Y",  "Ux,Uy,U"), NO "," before Ux
void Write_datHeader(int Ni, int Nj, int Nk, std::ofstream &write,
                     std::string var1, std::string var2) {
  write << "TITLE = \"Data\"\n"
        << "variables = " << var1 << "," << var2 << "\n"
        << "Zone t=\"data\"\n"
        << "I=" << Ni << ",J=" << Nj << ",K=" << Nk << ",F=POINT" << std::endl;
}
}  // namespace wrt