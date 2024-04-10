/*CA method manager class*/
#pragma once

#include "head.h"
#include "legacy/legacy_lattice/legacyfield.h"
#include "nuc.h"
#include "grow.h"

// // #include "LatticeCommunicator.h"
// // a single cell in CA mesh
// template <typename T>
// struct Cell {
//   /*CA Mesh Index*/
//   int Id;  // Identity number of a Cell, also could be used as an index to find
//            // a certain Cell
//   int Grain_Id;  // connection with Grain: Grain.Id
//   /*State Index*/
//   /*
//    * characterize phase and neighborhood of a cell
//    * 0: Liquid
//    * 1: Not Liquid but at least one neigbor is Liquid
//    * -1:Neither itself and neighbors are Liquid
//    */
//   int State;
//   /*
//    * track the development of the mush zone within cell
//    * 0: Liquid (State = 0)
//    * 1: Not Liquid and mush zone is developing (solidfraction < 1)
//    * -1:Not Liquid and mush zone's development is complete (solidfraction = 1)
//    */
//   int Mush;
//   T DT_Critic;  // Critical Undercooling by normal distribution
//   /*Mush Zone Related Parameters*/
//   T A_max;
//   T A_min;
//   T Mush_Frac;       // mush zone fraction (gm)
//   T DMush_Frac;      // mush zone fraction change (gm)
//   T Env_Solid_Frac;  // solid fraction in the envelop (gsm)
//   T Solid_Frac;      // solid fraction (gs)
//   T Orine;
//   Cell(int id_) : Id(id_) {
//     Grain_Id = -1;
//     State = 0;
//     Mush = 0;
//     DT_Critic = 0;
//     A_max = 0;
//     A_min = 0;
//     Mush_Frac = 0;
//     DMush_Frac = 0;
//     Env_Solid_Frac = 0;
//     Solid_Frac = 0;
//     Orine = 0;
//   }
//   Cell() : Id(-1) {
//     Grain_Id = -1;
//     State = 1;
//     Mush = -1;
//     DT_Critic = 0;
//     A_max = 0;
//     A_min = 0;
//     Mush_Frac = 0;
//     DMush_Frac = 0;
//     Env_Solid_Frac = 0;
//     Solid_Frac = 0;
//     Orine = 0;
//   }
// };

// /*characterize Mush zone*/
// template <typename T>
// struct Grain {
//   int Id;           // Identity number of a Grain, sorted by adding order
//   int Cell_Id;      // corresponding Cell Id: Cell.Id
//   int Dendrite_Id;  // corresponding Dendrite Id
//   T x;  // Grain centre position, not cell centre
//   T y;  // Grain centre position, not cell centre
//   T ArmLen[4];  // length of arm in 4 directions
//   T Orien;      // 0-90 degree orientation angle
//   Grain() {}
//   Grain(T x_, T y_, int Grain_id_, int Cell_id_, T Orine_)
//       : x(x_), y(y_), Id(Grain_id_), Cell_Id(Cell_id_), Orien(Orine_) {
//     ArmLen[0] = 0.001f;  // set to a little value for initialize mush zone; if
//                          // single nucleation, can't be 0
//     ArmLen[1] = 0.001f;
//     ArmLen[2] = 0.001f;
//     ArmLen[3] = 0.001f;
//     Dendrite_Id = 0;
//   }
// };

// /*Dendrite*/
// struct Dendrite {
//   int Id;
//   std::vector<int> Cell_Ids;
//   std::vector<int> Grain_Ids;
//   Dendrite() {
//     Cell_Ids.reserve(10);
//     Grain_Ids.reserve(10);
//   }
//   Dendrite(int Id_, int Grain_Id_, int Cell_Id_) : Id(Id_) {
//     Cell_Ids.emplace_back(Cell_Id_);
//     Grain_Ids.emplace_back(Grain_Id_);
//   }
//   ~Dendrite() {}
//   void add(int Cell_Id_, int Grain_Id_) {
//     Cell_Ids.emplace_back(Cell_Id_);
//     Grain_Ids.emplace_back(Grain_Id_);
//   }
// };

template <typename T>
class CAManager {
 public:
  CAManager(Geometry2DLegacy<T> &geo, GandinConverter<T> &ConvCA_,
            UnitConvManager<T> &convman)
      : Geo(geo), ConvCA(ConvCA_), ConvMan(convman) {
    Ni = Geo.getNi();
    Nj = Geo.getNj();
    N = Geo.getN();
    Nuc_Count = 0;
    Accumulated_NucCount_Bulk = 0;
    Accumulated_NucCount_Surf = 0;
    ComputeDomain = N - 2 * Ni - 2 * Nj + 4;
    std::cout << "PartCoeff = " << Part_Coef << std::endl;

    ActiveCells.reserve(ComputeDomain);

    for (int i = 0; i < sizeof(announced) / sizeof(announced[0]); i++) {
      announced[i] = 0;
    }

    CA = new Cell<T>[N];
    for (int j = 1; j < Nj - 1; j++) {
      for (int i = 1; i < Ni - 1; i++) {
        int Id = i + Ni * j;
        /*Init: Id, x, y, Temp, Conc, Enth*/
        CA[Id] = Cell<T>(Id);
        /*peripheral*/
      }
    }

    Lat_DT_Mean_Bulk = ConvCA.Lattice_DT_Mean_Bulk;
    Lat_DT_Std_Bulk = ConvCA.Lattice_DT_Std_Bulk;
    Lat_DT_Mean_Surf = ConvCA.Lattice_DT_Mean_Surf;
    Lat_DT_Std_Surf = ConvCA.Lattice_DT_Std_Surf;
    Lat_Nuc_Dens_Bulk = ConvCA.Lattice_NucDens_Bulk;
    Lat_Nuc_Dens_Surf = ConvCA.Lattice_NucDens_Surf;
    Lat_Growth_Para = ConvCA.Lattice_GrowthPara;
    // phase diagram
    Lat_T_Melt = ConvCA.Lattice_T_Melt;
    Lat_T_Eute = ConvCA.Lattice_T_Eute;
    Lat_m_Solidus = ConvCA.Lattice_m_Sol;
    Lat_m_Liquidus = ConvCA.Lattice_m_Liq;
    Part_Coef = ConvCA.Part_Coef;

    LatHeat_SHeatCap = ConvMan.TempConv->Lattice_LatHeat_SHeatCap;
  }
  ~CAManager() { delete[] CA; }

  void Nucleation();
  void Grain_capture();

  T get_estimated_Max_TimeStep(T CA_TimeStep_Coeff_);

  void init_Mushzone(int Grain_Id_);
  void set_Mushzone(int Grain_Id_);

  /*test and study*/

  void Simple_Phase_Transition();

  /*get*////////////////////////
  int getNucCount() { return Nuc_Count; }
  T getSolifiedPercent() {
    return static_cast<T>(Nuc_Count) / static_cast<T>(ComputeDomain);
  }
  // phase diagram
  T get_Temp_Liquidus(T conc) { return Lat_T_Melt - conc * Lat_m_Liquidus; }
  T get_Temp_Liquidus_Eute(T conc) {
    T Temp_Li = Lat_T_Melt - conc * Lat_m_Liquidus;
    if (Temp_Li < Lat_T_Eute) {
      Temp_Li = Lat_T_Eute;
    }
    return Temp_Li;
  }
  T get_Temp_Solidus(T conc) { return Lat_T_Melt - conc * Lat_m_Solidus; }
  T get_Temp_Solidus_Eute(T conc) {
    T Temp_So = Lat_T_Melt - conc * Lat_m_Solidus;
    if (Temp_So < Lat_T_Eute) {
      Temp_So = Lat_T_Eute;
    }
    return Temp_So;
  }///////////////////////////

  // lattice communicator and Active_Grains Index
  // strore all active cells by Cell Id
  std::vector<int> ActiveCells;
  // store newly mushy cells in one single time step.
  std::vector<int> NewMushyCells;
  // store newly solid cells in one single time step. in latcomm, index to set
  // lat field to 0
  std::vector<int> NewSolidCells;
  // in grain capture step, add to ActiveCells
  std::vector<int> CapturedCells;

  // test
  std::vector<int> FluidCells;

 private:
  UnitConvManager<T> &ConvMan;
  GandinConverter<T> &ConvCA;
  Cell<T> *CA;
  Geometry2DLegacy<T> &Geo;
  // lbm2D<T, LatSet, LatStru> &lbmT;
  // lbm2D<T, LatSet, LatStru> &lbmC;
  Basiclbm2D<T> *lbmT;
  Basiclbm2D<T> *lbmC;
  /* lattice structure */
  int Ni;
  int Nj;
  int N;
  int ComputeDomain;  // ComputeDomain size, = N - 2*Ni - 2*Nj + 4
  /*counters*/
  int Nuc_Count;
  T Accumulated_NucCount_Bulk;
  T Accumulated_NucCount_Surf;

  // CA parameters from GandinConverter
  T Lat_DT_Mean_Bulk;
  T Lat_DT_Std_Bulk;
  T Lat_DT_Mean_Surf;
  T Lat_DT_Std_Surf;
  T Lat_Nuc_Dens_Surf;
  T Lat_Nuc_Dens_Bulk;
  T Lat_Growth_Para;
  // phase diagram
  T Lat_T_Melt;
  T Lat_T_Eute;
  T Lat_m_Solidus;
  T Lat_m_Liquidus;
  T Part_Coef;

  T LatHeat_SHeatCap;

  std::vector<int> Pre_NucSites;  // store nuc SITES by Cell Id(not nucleated yet)
                              // use vector for parallel
  std::vector<Grain<T>> Grains;  // visit grain by Grain Id
  std::vector<Dendrite> Dendrites;

  int trav_in[8] = {-1, -1, -1, 0, 0, 1, 1, 1};  //-1,-1,-1, 0, 0, 1, 1, 1
  int trav_jn[8] = {-1, 0, 1, -1, 1, -1, 0, 1};  //-1, 0, 1,-1, 1,-1, 0, 1
};
