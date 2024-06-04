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
#include "data_struct/lattice.h"
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
  ScalarField<CAType> State;
  // flag field std::uint8_t
  ScalarField<CAFlag> Flag;
  // solid fraction
  ScalarField<T> Fs;
  // delta solid fraction
  ScalarField<T> Delta_Fs;
  // curvature of Solid-Liquid interface
  ScalarField<T> Curvature;
  // Solid phase composition
  ScalarField<T> C_Solids;
  // excess rho
  VectorFieldAOS<T, LatSet::q> ExcessC;
  // colleted excess rho
  ScalarField<T> ExcessC_;
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
  ScalarField<CAType>& getState() { return State; }
  ScalarField<CAFlag>& getFlag() { return Flag; }
  ScalarField<T>& getFs() { return Fs; }
  ScalarField<T>& getDeltaFs() { return Delta_Fs; }
  ScalarField<T>& getCurvature() { return Curvature; }
  ScalarField<T>& getCSolids() { return C_Solids; }
  VectorFieldAOS<T, LatSet::q>& getExcessC() { return ExcessC; }
  ScalarField<T>& getExcessC_() { return ExcessC_; }
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

// define unique CA Field
struct STATEBase : public FieldBase<1> {};
struct FSBase : public FieldBase<1> {};
struct DELTAFSBase : public FieldBase<1> {};
struct CURVATUREBase : public FieldBase<1> {};
struct CSOLIDSBase : public FieldBase<1> {};
struct PREEXCESSCBase : public FieldBase<1> {};
struct EXCESSCBase : public FieldBase<1> {};


using STATE = GenericField<GenericArray<CAType>, STATEBase>;
template <typename T>
using FS = GenericField<GenericArray<T>, FSBase>;
template <typename T>
using DELTAFS = GenericField<GenericArray<T>, DELTAFSBase>;
template <typename T>
using CURVATURE = GenericField<GenericArray<T>, CURVATUREBase>;
template <typename T>
using CSOLIDS = GenericField<GenericArray<T>, CSOLIDSBase>;
template <typename T>
using PREEXCESSC = GenericField<GenericArray<T>, PREEXCESSCBase>;
template <typename T>
using EXCESSC = GenericField<GenericArray<T>, EXCESSCBase>;

template <typename T>
using CAFIELDS = TypePack<STATE, FS<T>, DELTAFS<T>, CURVATURE<T>, CSOLIDS<T>,
                          PREEXCESSC<T>, EXCESSC<T>, TEMPINIT<T>, CONCINIT<T>>;
template <typename T>
using REFFIELDS = TypePack<VELOCITY<T, 2>, CONC<T>, TEMP<T>>;
template <typename T>
using FIELDPACK = TypePack<CAFIELDS<T>, REFFIELDS<T>>;
template <typename T>
using ALLFIELDS = typename ExtractFieldPack<FIELDPACK<T>>::mergedpack;

template <typename T, typename LatSet>
class BlockZhuStefanescu2D : public BlockLatticeBase<T, LatSet, ALLFIELDS<T>> {
 private:
  // preferred growth angle to x-axis, manually chosen
  T Theta;
  // anisotropy coefficient, manually chosen
  T delta;
  // Gibbs-Thomson coefficient
  T GT;
  // initial composition C0
  // low limit of temperature
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
  // solid count
  std::size_t SolidCount;

  ZSConverter<T>& ConvCA;

 public:
  template <typename... FIELDPTRS>
  BlockZhuStefanescu2D(Block2D<T>& geo, ZSConverter<T>& convca,
                       std::tuple<FIELDPTRS...> fieldptrs, T delta, T theta);
  ~BlockZhuStefanescu2D() {}

  std::uint8_t getLevel() const { return this->BlockGeo.getLevel(); }
  std::size_t getSolidCount() const { return SolidCount; }

  // get field data
  std::vector<std::size_t>& getInterface() { return Interface; }
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

  void apply_grow() {
    UpdateInterface();
    UpdateCurvature();
    UpdateDeltaFs();
    Grow();
    DistributeExcessC();
  }

  bool hasNeighborType(std::size_t id, std::uint8_t type) const;
};

template <typename T, typename CALatSet>
class BlockZhuStefanescu2DManager
    : public BlockLatticeManagerBase<T, CALatSet, FIELDPACK<T>> {
      public:
      
 private:
  std::vector<BlockZhuStefanescu2D<T, CALatSet>> BlockZS;

  std::vector<std::vector<std::size_t>*> Interfaces;

  ZSConverter<T>& ConvCA;

  // preferred growth angle to x-axis, manually chosen
  T Theta;
  // anisotropy coefficient, manually chosen
  T delta;

  // --- CA Field ---
  // communication needed:
  // state field, solid fraction field, collected excess C field

  // block-local field:
  // delta solid fraction field, curvature of Solid-Liquid interface, Solid phase
  // composition, pre streamed excess C field

  // total solid count
  std::size_t SolidCount;

 public:
  template <typename INITVALUEPACK, typename... FIELDPTRTYPES>
  BlockZhuStefanescu2DManager(BlockGeometry2D<T>& blockgeo, ZSConverter<T>& convca, T delta_,
                              T theta, INITVALUEPACK& initvalues,
                              FIELDPTRTYPES*... fieldptrs)
      : BlockLatticeManagerBase<T, CALatSet, FIELDPACK<T>>(blockgeo, initvalues, fieldptrs...),
        ConvCA(convca), Theta(theta), delta(delta_), SolidCount(std::size_t{}) {
    // create BlockZhuStefanescu2D
    for (std::size_t i = 0; i < this->BlockGeo.getBlockNum(); ++i) {
      BlockZS.emplace_back(this->BlockGeo.getBlock(i), ConvCA,
                           ExtractFieldPtrs<T, CALatSet, FIELDPACK<T>>::getFieldPtrTuple(
                             i, this->Fields, this->FieldPtrs),
                           delta, Theta);
    }
    // init Interfaces
    for (auto& zs : BlockZS) {
      Interfaces.push_back(&(zs.getInterface()));
    }
  }

  void CAFieldDataInit(BlockGeometryHelper2D<T>& GeoHelper) {
    // data transfer
    // StateFM.InitCopy(GeoHelper, CAType::Boundary);
    this->template getField<FS<T>>().InitAndComm(GeoHelper);
    this->template getField<CSOLIDS<T>>().InitAndComm(GeoHelper);
    this->template getField<EXCESSC<T>>().InitAndComm(GeoHelper);

    // field init without data transfer
    this->template getField<DELTAFS<T>>().Init(T(0));
    this->template getField<CURVATURE<T>>().Init(T(0));
    this->template getField<PREEXCESSC<T>>().Init(T(0));
  }

  void Init() {
    BlockZS.clear();
    for (std::size_t i = 0; i < this->BlockGeo.getBlockNum(); ++i) {
      BlockZS.emplace_back(this->BlockGeo.getBlock(i), ConvCA,
                           ExtractFieldPtrs<T, CALatSet, FIELDPACK<T>>::getFieldPtrTuple(
                             i, this->Fields, this->FieldPtrs),
                           delta, Theta);
    }
    // init Interfaces
    Interfaces.clear();
    for (auto& zs : BlockZS) {
      Interfaces.push_back(&(zs.getInterface()));
    }
  }

  void Setup(std::size_t SiteId, int num = CALatSet::q) {
    // find block to setup
    std::size_t blockid;
    std::vector<int> nums;
    nums.resize(this->BlockGeo.getBlockNum(), 0);

    Vector<T, 2> loc_t = this->BlockGeo.getLoc_t(SiteId);
    for (int i = 0; i < this->BlockGeo.getBlockNum(); ++i) {
      if (this->BlockGeo.getBlock(i).getBaseBlock().isInside(loc_t)) {
        blockid = this->BlockGeo.getBlock(i).getIndex_t(loc_t);
        nums[i] = num;
        break;
      }
    }

    std::cout << "[Zhu-Stefanescu 2D CA]" << std::endl;
    // get preferred growth angle to x-axis
    std::cout << "preferred growth angle: " << Theta / M_PI << " Pi" << std::endl;
    // get initial undercooling
    // Tl(latth.getLatRhoInit()), Tl_eq(convca.get_LatTliq(latso.getLatRhoInit()))
    T Tl_eq = ConvCA.get_LatTliq(this->template getField<CONCINIT<T>>().getBlockField(0).get());
    T deltaT = Tl_eq - this->template getField<TEMPINIT<T>>().getBlockField(0).get();
    std::cout << "Initial undercooling: " << deltaT << " | "
              << ConvCA.TempConv.getPhysDTemp(deltaT) << std::endl;

    for (int i = 0; i < this->BlockGeo.getBlockNum(); ++i) {
      BlockZS[i].Setup(blockid, nums[i]);
    }

    Communicate();
    for (auto& zs : BlockZS) {
      zs.SimpleCapture();
    }
    Communicate();
  }

  std::vector<std::vector<std::size_t>*>& getInterfaces() { return Interfaces; }


  void Apply_SimpleCapture() {
#pragma omp parallel for num_threads(Thread_Num)
    for (auto& zs : BlockZS) {
      zs.apply_grow();
    }
    Communicate();
#pragma omp parallel for num_threads(Thread_Num)
    for (auto& zs : BlockZS) {
      zs.SimpleCapture();
    }
  }

  void Communicate() {
    // State field
    this->template getField<STATE>().NormalCommunicate();
    // Fs field
    this->template getField<FS<T>>().NormalCommunicate();
    // ExcessC field
    this->template getField<EXCESSC<T>>().NormalCommunicate();
    // ExcessC field post communication

    // Reversed Communicate should only apply to newly solidified cells's neighbors

    // #pragma omp parallel for num_threads(Thread_Num)
    //     for (BlockField<ScalarField<T>, T, 2>& blockF : ExcessCFM.getBlockFields()) {
    //       // data in Recvs is sent from blockF to Sends of nblockF
    //       for (BlockFieldComm<ScalarField<T>, T, 2>& comm : blockF.getComms()) {
    //         BlockField<ScalarField<T>, T, 2>* nblockF = comm.BlockF;
    //         std::size_t size = comm.getRecvs().size();
    //         for (int iArr = 0; iArr < blockF.Size(); ++iArr) {
    //           auto& nArray = nblockF->getField(iArr);
    //           const auto& Array = blockF.getField(iArr);
    //           for (std::size_t id = 0; id < size; ++id) {
    //             nArray[comm.getSend(id)] += Array[comm.getRecv(id)];
    //           }
    //         }
    //       }
    //     }
  }

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
    SolidCount += count;
    return SolidCount;
  }

  bool WillRefineBlockCells(BlockGeometryHelper2D<T>& GeoHelper) {
    const std::uint8_t LevelLimit = GeoHelper.getLevelLimit();
    bool willrefine = false;
    std::vector<bool> hasSolid(GeoHelper.getBlockCells().size(), false);
#pragma omp parallel for num_threads(Thread_Num)
    for (int icell = 0; icell < GeoHelper.getBlockCells().size(); ++icell) {
      // cell block
      const BasicBlock<T, 2>& cellblock = GeoHelper.getBlockCell(icell);
      if (cellblock.getLevel() <= LevelLimit) {
        T voxsize = cellblock.getVoxelSize();
        // find corresponding block
        for (int iblock = 0; iblock < this->BlockGeo.getBlockNum(); ++iblock) {
          const BasicBlock<T, 2>& block = this->BlockGeo.getBlock(iblock);
          const BasicBlock<T, 2>& baseblock =
            this->BlockGeo.getBlock(iblock).getBaseBlock();
          if (isOverlapped(cellblock, baseblock)) {
            // get CA State field
            const GenericArray<CAType>& StateF =
              this->template getField<STATE>().getBlockField(iblock).getField(0);
            const AABB<T, 2> intsec = getIntersection(cellblock, baseblock);
            int Nx = static_cast<int>(std::round(intsec.getExtension()[0] / voxsize));
            int Ny = static_cast<int>(std::round(intsec.getExtension()[1] / voxsize));
            // block start
            Vector<T, 2> blockstart = intsec.getMin() - block.getMin();
            int blockstartx = static_cast<int>(std::round(blockstart[0] / voxsize));
            int blockstarty = static_cast<int>(std::round(blockstart[1] / voxsize));

            for (int iy = 0; iy < Ny; ++iy) {
              for (int ix = 0; ix < Nx; ++ix) {
                std::size_t idblock =
                  (iy + blockstarty) * block.getNx() + ix + blockstartx;
                if (util::isFlag(StateF[idblock], CAType::Solid)) {
                  hasSolid[icell] = true;
                  // exit for loop
                  ix = Nx;
                  iy = Ny;
                }
              }
            }
          }
        }
      }
    }
    GeoHelper.TagRefineLayer(hasSolid, willrefine);

    return willrefine;
  }
};
}  // namespace CA
