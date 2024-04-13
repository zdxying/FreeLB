/* This file is part of FreeLB
 *
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 *
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with FreeLB. If
 * not, see <https://www.gnu.org/licenses/>.
 *
 */

// Zhu-Stefanescu (Z-S) model for dendrite growth
// Virtual front tracking model for the quantitative modeling
// of dendritic growth in solidification of alloys 2007
#pragma once

#include <algorithm>
#include <atomic>
#include <vector>

#include "ca/cazs.h"
#include "data_struct/block_lattice.h"
#include "lbm/lattice_set.h"
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


namespace CA {

// D2Q4
template <typename T>
struct D2Q4 : public Basic_Lattice_Set<2, 4> {
  static constexpr Vector<int, 2> c[q] = {{0, -1}, {-1, 0}, {1, 0}, {0, 1}};
  static constexpr T w[q] = {T(1) / T(4), T(1) / T(4), T(1) / T(4), T(1) / T(4)};
  static constexpr int opp[q] = {3, 2, 1, 0};
};

// D2Q8
template <typename T>
struct D2Q8 : public Basic_Lattice_Set<2, 8> {
  static constexpr Vector<int, 2> c[q] = {{-1, -1}, {0, -1}, {1, -1}, {-1, 0},
                                          {1, 0},   {-1, 1}, {0, 1},  {1, 1}};
  static constexpr T w[q] = {T(1) / T(20), T(1) / T(5),  T(1) / T(20), T(1) / T(5),
                             T(1) / T(5),  T(1) / T(20), T(1) / T(5),  T(1) / T(20)};
  static constexpr int opp[q] = {7, 6, 5, 4, 3, 2, 1, 0};
};

template <typename T, typename LatSet>
class ZhuStefanescu2D {
 private:
  int Ni;
  int Nj;
  int N;
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

  // interface cells
  std::vector<std::size_t> Interface;

  Geometry2D<T>& Geo;
  ZSConverter<T>& ConvCA;
  RhoLattice<T>& lbmSO;
  RhoLattice<T>& lbmTH;
  VectorFieldAOS<T, LatSet::d>& Velocity;

  // state field std::uint8_t
  ScalerField<CAType> State;
  // flag field std::uint8_t
  ScalerField<CAFlag> Flag;
  // solid fraction
  ScalerField<T> Fs;
  // delta solid fraction
  ScalerField<T> Delta_Fs;
  // curvature of Solid-Liquid interface
  ScalerField<T> Curvature;
  // Solid phase composition
  ScalerField<T> C_Solids;
  // excess rho
  VectorFieldAOS<T, LatSet::q> ExcessC;
  // colleted excess rho
  ScalerField<T> ExcessC_;
  // solid count
  std::size_t SolidCount;

  // nbr index
  std::array<int, LatSet::q> Delta_Index;

 public:
  ZhuStefanescu2D(ZSConverter<T>& convca, RhoLattice<T>& lbmso, RhoLattice<T>& lbmth,
                  BasicLattice<T, LatSet>& lbmns, T delta_, T theta, int siteid,
                  int num = LatSet::q);
  ~ZhuStefanescu2D() {}

  // get field data
  std::vector<std::size_t>& getInterface() { return Interface; }
  ScalerField<CAType>& getState() { return State; }
  ScalerField<CAFlag>& getFlag() { return Flag; }
  ScalerField<T>& getFs() { return Fs; }
  ScalerField<T>& getDeltaFs() { return Delta_Fs; }
  ScalerField<T>& getCurvature() { return Curvature; }
  ScalerField<T>& getCSolids() { return C_Solids; }
  VectorFieldAOS<T, LatSet::q>& getExcessC() { return ExcessC; }
  ScalerField<T>& getExcessC_() { return ExcessC_; }
  // setup
  void Setup(int id, int num);
  void UpdateInterface();
  void UpdateCurvature(T limit = 4);
  // get C_eq, interface equilibrium composition
  T getC_eq(int id);
  // get anisotropy function g(theta)
  T getanisotropy(int id);
  // get interface normal to x-axis
  T getPhi(int id);

  void UpdateDeltaFs();
  // grow, omp parallel is enabled, no erase operation
  // if cells are large, use Grow_erase_s may improve performance
  void Grow();
  void DistributeExcessC(int id, T excessC);
  void CollectExcessC();
  void SimpleCapture();

  void TypeConversion();

  void PreCapture_s(std::vector<int>& precaps);
  void PreCapture_omp(std::vector<int>& precaps, int start, int end);
  // void SLICapture();
  // void slicapture(ZScell2D<T>& cells);
  // void SLICapture_s(const std::vector<int>& precaps,
  //                   std::vector<ZScell2D<T>>& cells);
  T getSolidCountFracton();
  // calc position of the vertices of the polygonal for capture
  // get Lphi
  // inline T getLphi(T phi);
  // get vertice of polygonal
  // inline Vector2D<T> getVertice(int id);
  // check if a point is inside a polygon
  // bool InPolygon(Vector2D<T>& point, std::vector<Vector2D<T>>& polygon);
  // sort the vertices of the polygonal
  // void SortPolygon();
  // // check if a point is on a line segment
  // inline bool OnSegment(Vector2D<T>& p1, Vector2D<T>& p2, Vector2D<T>& point);
  // inline bool OnSegment_(Vector2D<T>& p1, Vector2D<T>& p2, Vector2D<T>& point);
  // return 3 state -1, 0, 1, compared with epsilon and 0
  inline int sign(T x);
  // get position of neighbors with State = -1
  // template <int start = 0>
  // inline void getNbrs_pos(int id, std::vector<int>& nbrs);

  // experimental
  inline T getAveNbrPopRho(int id);
  inline T getStatisticalPopRho();

  void apply_SimpleCapture() {
    UpdateInterface();
    UpdateCurvature();
    UpdateDeltaFs();
    Grow();
    CollectExcessC();
    SimpleCapture();
    TypeConversion();
  }
  // void apply_SLICapture() {
  //   UpdateCurvature();
  //   UpdateDeltaFs();
  //   Grow();
  //   CollectExcessC();
  //   TypeConversion();
  // }

  bool hasNeighborType(int id, std::uint8_t type) const;
  bool hasNeighborFlag(int id, std::uint8_t flag) const;
};

// --------------------------------------------------------------------
// -------------------------BlockZhuStefanescu2D-----------------------
// --------------------------------------------------------------------

template <typename T, typename LatSet>
struct BlockZSCommStru {
  BlockZhuStefanescu2D<T, LatSet>* SendBlock;
  BlockComm<T, LatSet::d>* Comm;

  BlockZSCommStru(BlockZhuStefanescu2D<T, LatSet>* sblock,
                  BlockComm<T, LatSet::d>* blockcomm)
      : SendBlock(sblock), Comm(blockcomm) {}

  std::vector<std::size_t>& getSends() { return Comm->SendCells; }
  std::vector<std::size_t>& getRecvs() { return Comm->RecvCells; }
};

template <typename T, typename LatSet>
class BlockZhuStefanescu2D {
 private:
  // preferred growth angle to x-axis, manually chosen
  T Theta;
  // anisotropy coefficient, manually chosen
  T delta;
  // Gibbs-Thomson coefficient
  T GT;
  // initial composition C0
  T C0;
  // low limit of temperature
  T Tl;
  // equilibrium liquidus temperature at the initial composition C0
  T Tl_eq;
  // slope of liquidus line
  T m_l;
  // partition coefficient
  T Part_Coef;
  // (1 - partition coefficient)
  T _Part_Coef;

  // interface cells
  std::vector<std::size_t> Interface;

  Block2D<T>& Geo;
  ZSConverter<T>& ConvCA;

  ScalerField<T>& Conc;
  ScalerField<T>& Temp;

  VectorFieldAOS<T, 2>& Velocity;
  // COMM
  std::vector<BlockZSCommStru<T, LatSet>> Communicators;

  // state field std::uint8_t
  ScalerField<CAType>& State;
  // solid fraction
  ScalerField<T>& Fs;
  // delta solid fraction
  ScalerField<T>& Delta_Fs;
  // curvature of Solid-Liquid interface
  ScalerField<T>& Curvature;
  // Solid phase composition
  ScalerField<T>& C_Solids;
  // excess rho
  // cyclic array is used to prevent race condition
  // PopulationField<T, LatSet::q> ExcessC;
  // pre streamed excess C field
  ScalerField<T>& PreExcessC;
  // colleted excess rho
  ScalerField<T>& ExcessC;
  // solid count
  std::size_t SolidCount;

  // nbr index
  std::array<int, LatSet::q> Delta_Index;

 public:
  BlockZhuStefanescu2D(BlockField<VectorFieldAOS<T, 2>, T, 2>& veloFM,
                       ZSConverter<T>& convca, BlockRhoLattice<T>& latso,
                       BlockRhoLattice<T>& latth, ScalerField<CAType>& state,
                       ScalerField<T>& fs, ScalerField<T>& delta_fs,
                       ScalerField<T>& curvature, ScalerField<T>& csolids,
                       ScalerField<T>& preexcessc, ScalerField<T>& excessc, T delta,
                       T theta, std::size_t siteid, int num = LatSet::q);
  ~BlockZhuStefanescu2D() {}

  Block2D<T>& getGeo() { return Geo; }
  std::vector<BlockZSCommStru<T, LatSet>>& getCommunicators() { return Communicators; }

  std::uint8_t getLevel() const { return Geo.getLevel(); }
  std::size_t getSolidCount() const { return SolidCount; }

  // get field data
  std::vector<std::size_t>& getInterface() { return Interface; }
  ScalerField<CAType>& getState() { return State; }
  ScalerField<T>& getFs() { return Fs; }
  ScalerField<T>& getDeltaFs() { return Delta_Fs; }
  ScalerField<T>& getCurvature() { return Curvature; }
  ScalerField<T>& getCSolids() { return C_Solids; }
  ScalerField<T>& getPreExcessC() { return PreExcessC; }
  ScalerField<T>& getExcessC() { return ExcessC; }
  // setup
  void Setup(std::size_t id, int num);
  void UpdateInterface();
  void UpdateCurvature(T limit = 4);
  // get C_eq, interface equilibrium composition
  T getC_eq(std::size_t id);
  // get anisotropy function g(theta)
  T getanisotropy(std::size_t id);
  // get interface normal to x-axis
  T getPhi(std::size_t id);

  void UpdateDeltaFs();
  // grow, omp parallel is enabled, no erase operation
  // if cells are large, use Grow_erase_s may improve performance
  void Grow();
  // this function should be called after communication
  void DistributeExcessC();
  void SimpleCapture();

  void apply_SimpleCapture() {
    DistributeExcessC();
    SimpleCapture();
  }

  void apply_grow() {
    UpdateInterface();
    UpdateCurvature();
    UpdateDeltaFs();
    Grow();
  }

  void communicate();

  bool hasNeighborType(std::size_t id, std::uint8_t type) const;
};

template <typename T, typename CALatSet>
class BlockZhuStefanescu2DManager {
 private:
  std::vector<BlockZhuStefanescu2D<T, CALatSet>> BlockZS;

  std::vector<std::vector<std::size_t>*> Interfaces;

  BlockGeometry2D<T>& BlockGeo;

  // preferred growth angle to x-axis, manually chosen
  T Theta;
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

  // --- CA Field ---
  // state field
  BlockFieldManager<ScalerField<CAType>, T, 2> StateFM;
  // solid fraction field
  BlockFieldManager<ScalerField<T>, T, 2> FsFM;
  // delta solid fraction field
  BlockFieldManager<ScalerField<T>, T, 2> DeltaFsFM;
  // curvature of Solid-Liquid interface
  BlockFieldManager<ScalerField<T>, T, 2> CurvFM;
  // Solid phase composition
  BlockFieldManager<ScalerField<T>, T, 2> CSolidsFM;
  // pre streamed excess C field
  // BlockFieldManager<PopulationField<T, CALatSet::q> ,T , 2> ExcessCFM;
  // comapared with population field, the excess C field has less data exchange, use
  // scalar field instead
  BlockFieldManager<ScalerField<T>, T, 2> PreExcessCFM;
  // collected excess C field
  BlockFieldManager<ScalerField<T>, T, 2> ExcessCFM;

 public:
  BlockFieldManager<ScalerField<CAType>, T, 2>& getStateFM() { return StateFM; }
  BlockFieldManager<ScalerField<T>, T, 2>& getFsFM() { return FsFM; }
  BlockFieldManager<ScalerField<T>, T, 2>& getDeltaFsFM() { return DeltaFsFM; }
  BlockFieldManager<ScalerField<T>, T, 2>& getCurvFM() { return CurvFM; }
  BlockFieldManager<ScalerField<T>, T, 2>& getCSolidsFM() { return CSolidsFM; }
  BlockFieldManager<ScalerField<T>, T, 2>& getPreExcessCFM() { return PreExcessCFM; }
  BlockFieldManager<ScalerField<T>, T, 2>& getExcessCFM() { return ExcessCFM; }

  template <typename LatSet0, typename LatSet1>
  BlockZhuStefanescu2DManager(BlockFieldManager<VectorFieldAOS<T, 2>, T, 2>& veloFM,
                              ZSConverter<T>& convca,
                              BlockLatticeManager<T, LatSet0>& LatSos,
                              BlockLatticeManager<T, LatSet1>& LatThs, T delta_, T theta,
                              std::size_t SiteId, int num = CALatSet::q)
      : BlockGeo(veloFM.getGeo()), Theta(theta), delta(delta_), GT(0.1), C0(0.1),
        Tl_eq(0.1), m_l(0.1), Part_Coef(0.1), _Part_Coef(1 - Part_Coef),
        StateFM(veloFM.getGeo(), CAType::Boundary), FsFM(veloFM.getGeo(), T(0)),
        DeltaFsFM(veloFM.getGeo(), T(0)), CurvFM(veloFM.getGeo(), T(0)),
        CSolidsFM(veloFM.getGeo(), T(0)), PreExcessCFM(veloFM.getGeo(), T(0)),
        ExcessCFM(veloFM.getGeo(), T(0)) {
    // find block to setup
    std::vector<std::size_t> blockids;
    blockids.resize(BlockGeo.getBlockNum(), std::size_t(0));

    Vector<T, 2> loc_t = BlockGeo.getLoc_t(SiteId);
    for (std::size_t i = 0; i < BlockGeo.getBlockNum(); ++i) {
      if (BlockGeo.getBlock(i).getBaseBlock().isInside(loc_t)) {
        blockids[i] = BlockGeo.getBlock(i).getIndex_t(loc_t);
      }
    }
    // create BlockZhuStefanescu2D
    for (std::size_t i = 0; i < BlockGeo.getBlockNum(); ++i) {
      BlockZS.emplace_back(
        veloFM.getBlockField(i), convca, LatSos.getBlockLat(i), LatThs.getBlockLat(i),
        StateFM.getBlockField(i).getField(), FsFM.getBlockField(i).getField(),
        DeltaFsFM.getBlockField(i).getField(), CurvFM.getBlockField(i).getField(),
        CSolidsFM.getBlockField(i).getField(), PreExcessCFM.getBlockField(i).getField(),
        ExcessCFM.getBlockField(i).getField(), delta_, theta, blockids[i], num);
    }
    // init States, ExcessC_s, Interfaces
    for (auto& zs : BlockZS) {
      Interfaces.push_back(&(zs.getInterface()));
    }
    InitCommunicators();
  }

  // freeze cells based on CA field(if all points are solid, freeze the blockcell)
  void FreezeBlockCells(BlockGeometryHelper2D<T>& GeoHelper) {
#pragma omp parallel for num_threads(Thread_Num)
    for (BasicBlock<T, 2>& block : GeoHelper.getBlockCells()) {
      if (util::isFlag(GeoHelper.getBlockCellTag(block.getBlockId()),
                       BlockCellTag::Solid)) {
        continue;
      }
      Vector<T, 2> centre = block.getCenter();
      // find block
      int blockFid = 0;
      for (auto& zs : BlockZS) {
        if (zs.getGeo().isInside(centre)) {
          blockFid = zs.getGeo().getBlockId();
          break;
        }
      }
      BasicBlock<T, 2>& blockF = BlockGeo.getBlock(blockFid);
      GenericArray<CAType>& StateF = BlockZS[blockFid].getState().getField(0);

      // get start index
      Vector<T, 2> Ext = block.getMinCenter() - blockF.getMinCenter();
      int x = static_cast<int>(Ext[0] / blockF.getCellSize());
      int y = static_cast<int>(Ext[1] / blockF.getCellSize());
      // check if all points within the blockcell are solid
      bool freeze = true;
      for (int iy = y; iy < y + block.getNy(); ++iy) {
        for (int ix = x; ix < x + block.getNx(); ++ix) {
          std::size_t id = ix + iy * blockF.getNx();
          if (!util::isFlag(StateF[id], CAType::Solid)) {
            freeze = false;
            break;
          }
        }
      }
      if (freeze) {
        GeoHelper.getBlockCellTag(block.getBlockId()) = BlockCellTag::Solid;
      }
    }
  }


  std::vector<std::vector<std::size_t>*>& getInterfaces() { return Interfaces; }

  std::size_t getInterfaceNum() {
    std::size_t num = 0;
    for (auto& vec : Interfaces) {
      num += vec->size();
    }
    return num;
  }
  std::size_t getSolidCount() {
    std::size_t count = 0;
    for (auto& zs : BlockZS) {
      count += zs.getSolidCount();
    }
    return count;
  }

  void apply_SimpleCapture() {
#pragma omp parallel for num_threads(Thread_Num)
    for (auto& zs : BlockZS) {
      zs.apply_grow();
    }
    Communicate();
#pragma omp parallel for num_threads(Thread_Num)
    for (auto& zs : BlockZS) {
      zs.apply_SimpleCapture();
    }
  }

  void InitCommunicators() {
    for (auto& zs : BlockZS) {
      Block2D<T>& Geo = zs.getGeo();
      std::vector<BlockZSCommStru<T, CALatSet>>& Communicators = zs.getCommunicators();
      for (BlockComm<T, CALatSet::d>& comm : Geo.getCommunicators()) {
        Communicators.emplace_back(&BlockZS[comm.getSendId()], &comm);
      }
    }
  }

  void Communicate() {
#pragma omp parallel for num_threads(Thread_Num)
    for (auto& zs : BlockZS) {
      zs.communicate();
    }
  }
};
}  // namespace CA
