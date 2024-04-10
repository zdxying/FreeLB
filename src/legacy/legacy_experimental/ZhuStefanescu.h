// Zhu-Stefanescu (Z-S) model for dendrite growth
// Virtual front tracking model for the quantitative modeling
// of dendritic growth in solidification of alloys 2007
#pragma once

#include <algorithm>
#include <atomic>
#include <vector>

#include "lbm/lattice_set.h"
#include "legacy/legacy_lattice/legacyfield.h"
#include "legacy/legacy_lbm/legacy_lbmanager2d.h"
#include "lbm/unit_converter.h"
#include "utils/fdm_solver.h"
#include "utils/util.h"

// curvature of Solid-Liquid interface, K is calculated by:
//  K = ((df/dx)^2 + (df/dy)^2)^(-3/2) *
// (2*df/dx*df/dy*d^2f/dxdy - (df/dx)^2*d^2f/dy^2 - (df/dy)^2*d^2f/dx^2)
// f(x,y) is the solid fraction

// function accounting for the anisotropy of the surface tension
// g(phi, theta) = 1 - delta * cos(4*(phi - theta))
// where: delta is the anisotropy coefficient, manually chosen
// phi is growth angle to x-axis, which can be calculated by solid fraction:
// phi = arctan(df/dy / df/dx), i.e., tan(phi) = df/dy / df/dx
// theta is the preferred growth angle to x-axis
// then: g(theta) = 1 - delta * cos(4*arctan(df/dy / df/dx) - 4*theta)

// the driving force for dendritic growth is considered to be controlled by
// the difference between local interface equilibrium composition and
// local actual liquid composition
// deltaf is the increased solid fraction during one time step, given by:
// deltaf = (C_eq - Cl)/(C_eq*(1-k)), where:
// C_eq is the interface equilibrium composition, given by:
// C_eq = C0 + [(T_interface - Tl_eq) + GT*K*g(theta)] / m
// C0 is the initial composition,
// Cl is the actual liquid composition
// k is the partition coefficient

// rejected solute in an interface cell at each timestep
// deltaC = C_l(1-k)*deltaf

// virtual interface tracking scheme:
// Lphi = cellsize/max(|cos(phi)|, |sin(phi)|), where cellsize = 1
// Lphi is measured from cell center along the direction normal to local SL
// interface, position of the sharp SL interface Ip is given by:
// Ip = Lphi * fs, by connecting the point Ip of all interface cells,
// the sharp interface can thus be obtained

// nucleation: use single nucleation
// capture: virtual interface tracking scheme

template <typename T>
struct ZScell2D {
  // cell id
  int Id;
  // delta solid fraction
  T deltaf;
  // curvature of Solid-Liquid interface
  T K;
  // Solid phase composition
  Fraction<T> C_Solid;
  // vertice of polygonal for capture
  // Vector2D<T> Vertice;
  ZScell2D() : Id(0), K(T(0)), deltaf(T(0)), C_Solid() {}
  ZScell2D(int id) : Id(id), K(T(0)), deltaf(T(0)), C_Solid() {}
  // void SetVertice(const Vector2D<T>& vertice) { Vertice = vertice; }
  void Addto_Cs(T num, T denum) { C_Solid.add(num, denum); }
  void set_Cs(T num, T denum) { C_Solid.set(num, denum); }
};

template <typename T>
struct SolifiedZScell2D {
  // cell id
  int Id;
  // revised deltaf in the last timestep
  T deltaf;

  // constructor
  SolifiedZScell2D() : Id(0), deltaf(0) {}
  SolifiedZScell2D(int id, T deltaf_) : Id(id), deltaf(deltaf_) {}
};

// must use D2Q8 as LatStru
// multi-nuclei is NOT supported inside this class, to use multi-nuclei,
// it is recommended to instantiate multiple objects of this class
template <typename T, typename LatStru>
class ZhuStefanescu2D : public FDM2D<T>, public LatStru {
 private:
  int Ni;
  int Nj;
  // curvature of Solid-Liquid interface, stored in CAfield
  // T* K;
  // preferred growth angle to x-axis, manually chosen
  T Theta;
  // nucleation sites
  int SiteId;
  // anisotropy coefficient, manually chosen
  T delta;
  // Gibbs-Thomson coefficient
  T GT;
  // initial composition C0
  T C0;
  // equilibrium liquidus temperature at the initial composition C0
  T Tl_eq;
  // slope of liquidus line
  T m_l;
  // partition coefficient
  T Part_Coef;
  // (1 - partition coefficient)
  T _Part_Coef;

  // interface cells' Idx
  std::vector<int> CellsIdx;

#ifdef _OPENMP
  std::vector<std::vector<ZScell2D<T>>> Cells_OMP;
  std::vector<std::vector<SolifiedZScell2D<T>>> Solifieds_OMP;
  std::vector<std::vector<int>> PreCaps_OMP;
#else
  // newly solified cells
  std::vector<SolifiedZScell2D<T>> Solifieds;
#endif
  // interface cells
  std::vector<ZScell2D<T>> Cells;
  // SLI/mixed Capture:
  std::vector<int> PreCaps;

  // Simple Capture: add newly captured cells serially no need to use omp
  // std::vector<ZScell2D<T>> Captureds;

  CAZSField2D<T>& CA;
  Geometry2DLegacy<T>& Geo;
  ZSConverter<T>& ConvCA;
  LBManager2D<T>& LBM;
  Basiclbm2D<T>& lbmSO;
  Basiclbm2D<T>& lbmTH;
  // velocity field
  Velocity2D<T>& Field;

 public:
  ZhuStefanescu2D(CAZSField2D<T>& ca, ZSConverter<T>& convca,
                  LBManager2D<T>& lbm, T delta_, T theta, int siteid,
                  int num = 4);
  ZhuStefanescu2D() = delete;
  ~ZhuStefanescu2D() {}

  // setup
  // one cell State = -1, 4 or 8 nearest neighbors State = 0
  void Setup(int id, int num);
  // erase inactive cells
  void EraseInactive();
  // erase inactive cells, serially
  inline void EraseInactive_s(std::vector<ZScell2D<T>>& cells);

  // grow, omp parallel is enabled, no erase operation
  // if cells are large, use Grow_erase_s may improve performance
  void Grow();
  // Grow, serially
  inline void Grow_s(std::vector<ZScell2D<T>>& cells,
                     std::vector<SolifiedZScell2D<T>>& solifieds);
  // deal with newly solified cells serially
  inline void Solified_s(const std::vector<SolifiedZScell2D<T>>& solifieds);
  // grow while erase, using while loop and iterator, serially
  // divide cells into sub-vectors and execute in parallel is recommended
  void Grow_erase();
  inline void Grow_erase_s(std::vector<ZScell2D<T>>& cells,
                           std::vector<SolifiedZScell2D<T>>& solifieds);

  // serially get precaptured cells, DO NOT call this function parallelly
  void PreCapture_s(std::vector<int>& precaps, std::vector<ZScell2D<T>>& cells);

  // capture scheme:
  // 1. simple capture:
  // capture cells which have solid neighbors (at least one out of 8)
  void SimpleCapture();
  // simpled capture, serially
  void SimpleCapture_s(std::vector<ZScell2D<T>>& cells);
  // 2. virtual interface tracking scheme Zhu-Stefanescu
  // capture cells which were inside virtual SL interface polygonal
  // DO NOT use grow_erase() with this capture scheme
  void SLICapture();
  inline void slicap_s(int id, std::vector<ZScell2D<T>> &cells);
  inline void SLICapture_s(std::vector<int>& precaps,
                           std::vector<ZScell2D<T>>& cells);
  // 3. mixed capture:
  // capture cells which have solid neighbors(at least one out of 4)
  // and using virtual interface tracking scheme
  // void MixedCapture();

  // get, transform assumes that the output range has the same size as the input
  void TransformInterface() {
    CellsIdx.clear();
#ifdef _OPENMP
    for (auto& cells : Cells_OMP) {
      std::transform(cells.begin(), cells.end(), std::back_inserter(CellsIdx),
                     [](const ZScell2D<T>& cell) { return cell.Id; });
    }
#else
    std::transform(Cells.begin(), Cells.end(), std::back_inserter(CellsIdx),
                   [](const ZScell2D<T>& cell) { return cell.Id; });
#endif
  }
  std::vector<int>& getInterface() { return CellsIdx; }
  std::vector<int>& Trans_GetInterface() {
    TransformInterface();
    return CellsIdx;
  }

  inline T getK(int id);
  inline T getK_(int id);

  // call: &ZhuStefanescu2D<T, LatStru>::getK
  template <T (ZhuStefanescu2D<T, LatStru>::*get_K)(int) =
                &ZhuStefanescu2D<T, LatStru>::getK>
  void Calc_K();
  template <T (ZhuStefanescu2D<T, LatStru>::*get_K)(int)>
  inline void Calc_K_s(std::vector<ZScell2D<T>>& cells);
  // get curvature of Solid-Liquid interface

  // calc increased solid fraction during one time step
  // get interface normal to x-axis
  inline T getPhi(int id);
  // get anisotropy function g(theta)
  inline T getanisotropy(int id);
  // get C_eq, interface equilibrium composition
  inline T getC_eq(ZScell2D<T>& cell);
  // get deltaf, increased solid fraction during one time step
  inline T getdeltaf(ZScell2D<T>& cell);

  // calc position of the vertices of the polygonal for capture
  // get Lphi
  inline T getLphi(T phi);
  // get vertice of polygonal
  inline Vector2D<T> getVertice(int id);
  // get polygonal
  void getPolygon();
  // check if a point is inside a polygon
  bool InPolygon(Vector2D<T>& point, std::vector<Vector2D<T>>& polygon);
  // sort the vertices of the polygonal
  void SortPolygon();
  // check if a point is on a line segment
  inline bool OnSegment(Vector2D<T>& p1, Vector2D<T>& p2, Vector2D<T>& point);
  inline bool OnSegment_(Vector2D<T>& p1, Vector2D<T>& p2, Vector2D<T>& point);
  // return 3 state -1, 0, 1, compared with epsilon and 0
  inline int sign(T x);
  // get position of neighbors with State = -1
  template <int start = 0>
  inline void getNbrs_pos(int id, std::vector<int>& nbrs);
  // detect solid neighbors, return true if has solid neighbors
  inline bool has_SolidNbr(int id, int end = LatStru::q) {
    for (int i = 0; i < end; i++) {
      if (CA.State[id + LatStru::Nbr[i]] == -1) return true;
    }
    return false;
  }
  //   inline bool has_SolidNbr_simd(int id) {
  //     bool inactive = false;
  // #pragma omp simd reduction(|| : inactive)
  //     for (int i = 1; i < LatSet::q; i++) {
  //       if (CA.State[Id + LatStru::Nbr[i]] == -1) inactive = true;
  //     }
  //     return inactive;
  //   }

  // experimental
  // reject solute to neighbors
  inline void RejectSolutetoNbr(int id, T deltaf);
  // reject solute to neighbors and cell itself
  inline void RejectSolute(int id, T deltaf);
  inline T getAveNbrPopRho(int id);
  inline T getStatisticalPopRho();

  void apply_SimpleCapture() {
    Calc_K();
    Grow();
    EraseInactive();
    SimpleCapture();
  }
  void apply_SLICapture() {
    Calc_K();
    Grow();
    SLICapture();
    // erase inactive cells after capture
    EraseInactive();
  }
  // void apply_MixedCapture() {
  //   Calc_K();
  //   Grow();
  //   MixedCapture();
  //   EraseInactive();
  // }
};

struct Curvature2D {};
