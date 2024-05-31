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

// block_lattice.h

#pragma once

#include "data_struct/lattice.h"
#include "utils/fdm_solver.h"

template <typename T, typename LatSet, typename TypePack>
struct BlockLatComm {
  BlockLattice<T, LatSet, TypePack>* SendBlock;
  BlockComm<T, LatSet::d>* Comm;

  BlockLatComm(BlockLattice<T, LatSet, TypePack>* sblock,
               BlockComm<T, LatSet::d>* blockcomm)
      : SendBlock(sblock), Comm(blockcomm) {}

  std::vector<std::size_t>& getSends() { return Comm->SendCells; }
  std::vector<std::size_t>& getRecvs() { return Comm->RecvCells; }
};

template <typename T, typename LatSet, typename TypePack>
struct InterpBlockLatComm {
  BlockLattice<T, LatSet, TypePack>* SendBlock;
  InterpBlockComm<T, LatSet::d>* Comm;

  InterpBlockLatComm(BlockLattice<T, LatSet, TypePack>* sblock,
                     InterpBlockComm<T, LatSet::d>* interpcomm)
      : SendBlock(sblock), Comm(interpcomm) {}

  std::vector<std::size_t>& getRecvs() { return Comm->RecvCells; }
  std::vector<InterpSource<LatSet::d>>& getSends() { return Comm->SendCells; }
  // std::vector<InterpWeight<T, LatSet::d>>& getWeights() { return Comm->InterpWeights; }
};


// // a base lattice for BlockLattice
// template <typename T>
// class BlockRhoLattice {
//  protected:
//   ScalarField<T>& Rho;
//   // converter
//   AbstractConverter<T>& Conv;

//  public:
//   BlockRhoLattice(AbstractConverter<T>& conv, ScalarField<T>& rho)
//       : Conv(conv), Rho(rho) {}

//   ScalarField<T>& getRhoField() { return Rho; }
//   const ScalarField<T>& getRhoField() const { return Rho; }

//   const T& getRho(int i) const { return Rho.get(i); }
//   T& getRho(int i) { return Rho.get(i); }

//   void SetRhoField(std::size_t id, T value) { Rho.SetField(id, value); }

//   T getLatRhoInit() const { return Conv.getLatRhoInit(); }
//   T getLatgBeta(std::uint8_t level) const {
//     return RefineConverter<T>::getLattice_gbetaF(Conv.getLattice_gbeta(), level);
//   }
// };

// block structure for refined lattice
template <typename T, typename LatSet, typename TypePack>
class BlockLattice {
 public:
  using CellType = BCell<T, LatSet, TypePack>;
  using LatticeSet = LatSet;

 protected:
  // nbr index
  std::array<int, LatSet::q> Delta_Index;
  // geometry
  Block<T, LatSet::d>& BlockGeo;

  // field
  FieldPtrCollection<TypePack> Fields;

  // --- lattice communication structure ---
  // conmmunicate with same level block
  std::vector<BlockLatComm<T, LatSet, TypePack>> Communicators;
  // average blocklat comm, get from higher level block
  std::vector<InterpBlockLatComm<T, LatSet, TypePack>> AverageComm;
  // interp blocklat comm, get from lower level block
  std::vector<InterpBlockLatComm<T, LatSet, TypePack>> InterpComm;

  // omega = 1 / tau
  T Omega;
  // 1 - omega
  T _Omega;
  // 1 - omega/2
  T fOmega;

  // tolerance
  T RhoRes;
  std::vector<T> RhoOld;
  T URes;
  std::vector<Vector<T, LatSet::d>> UOld;

 public:
  // BlockLattice(Block<T, LatSet::d>& block, ScalarField<T>& rho,
  //              VectorFieldAOS<T, LatSet::d>& velocity,
  //              PopulationField<T, LatSet::q>& pops, AbstractConverter<T>& conv,
  //              bool initpop = true);
  template <typename... FIELDPTRS>
  BlockLattice(Block<T, LatSet::d>& block, AbstractConverter<T>& conv,
               std::tuple<FIELDPTRS...> fieldptrs);

  template <typename FieldType>
  auto& getField() {
    return Fields.template getField<FieldType>();
  }
  template <typename FieldType>
  const auto& getField() const {
    return Fields.template getField<FieldType>();
  }

  std::array<T*, LatSet::q> getPop(std::size_t id) {
    return getField<POP<T, LatSet::q>>().getArray(id);
  }

  std::size_t getNbrId(std::size_t id, int dir) const { return id + Delta_Index[dir]; }

  Block<T, LatSet::d>& getGeo() { return BlockGeo; }
  const Block<T, LatSet::d>& getGeo() const { return BlockGeo; }
  int getNx() const { return BlockGeo.getNx(); }
  int getNy() const { return BlockGeo.getNy(); }
  int getNz() const { return BlockGeo.getNz(); }
  std::size_t getN() const { return BlockGeo.getN(); }
  int getOverlap() const { return BlockGeo.getOverlap(); }
  inline T getOmega() const { return Omega; }
  inline T get_Omega() const { return _Omega; }
  inline T getfOmega() const { return fOmega; }
  std::uint8_t getLevel() const { return BlockGeo.getLevel(); }
  const Vector<int, LatSet::d>& getProjection() const { return BlockGeo.getProjection(); }
  const std::array<int, LatSet::q>& getDelta_Index() const { return Delta_Index; }

  std::vector<BlockLatComm<T, LatSet, TypePack>>& getCommunicators() {
    return Communicators;
  }
  std::vector<InterpBlockLatComm<T, LatSet, TypePack>>& getAverageComm() {
    return AverageComm;
  }
  std::vector<InterpBlockLatComm<T, LatSet, TypePack>>& getInterpComm() {
    return InterpComm;
  }

  // convert from fine to coarse, call after all average communication(including pops)
  void PopConvFineToCoarse();
  // convert from coarse to fine, call after all interp communication(including pops)
  void PopConvCoarseToFine();
  // normal communication, which can be done using normalcommunicate() in blockFM
  void communicate();
  // average communication
  void avercommunicate();
  // interp communication
  void interpcommunicate();

  void Stream();

  template <typename CELLDYNAMICS, typename ArrayType>
  void ApplyCellDynamics(const ArrayType& flagarr);

  // tolerance
  void EnableToleranceRho(T rhores = T(1e-5));
  void EnableToleranceU(T ures = T(1e-5));
  T getToleranceRho();
  T getToleranceU();
  // get inner block tolerance
  T getTolRho(int shift = 1);
  T getTolU(int shift = 1);
};

// block lattice manager
template <typename T, typename LatSet, typename TypePack>
class BlockLatticeManager {
 public:
  using FIELDS = typename ExtractFieldPack<TypePack>::pack1;
  using FIELDPTRS = typename ExtractFieldPack<TypePack>::pack2;
  using ALLFIELDS = typename ExtractFieldPack<TypePack>::mergedpack;
  using BLOCKLATTICE = BlockLattice<T, LatSet, ALLFIELDS>;

 private:
  std::vector<BlockLattice<T, LatSet, ALLFIELDS>> BlockLats;
  BlockGeometry<T, LatSet::d>& BlockGeo;
  AbstractConverter<T>& Conv;

  // field
  BlockFieldManagerCollection<T, LatSet, FIELDS> Fields;
  BlockFieldManagerPtrCollection<T, LatSet, FIELDPTRS> FieldPtrs;

 public:
  template <typename... FIELDPTRTYPES>
  BlockLatticeManager(BlockGeometry<T, LatSet::d>& blockgeo, AbstractConverter<T>& conv,
                      FIELDPTRTYPES&... fieldptrs);
  template <typename INITVALUEPACK, typename... FIELDPTRTYPES>
  BlockLatticeManager(BlockGeometry<T, LatSet::d>& blockgeo, INITVALUEPACK& initvalues,
                      AbstractConverter<T>& conv, FIELDPTRTYPES&... fieldptrs);

  template <typename FieldType>
  auto& getField() {
    if constexpr (isTypeInTuple<FieldType, FIELDS>::value) {
      return Fields.template getField<FieldType>();
    } else if constexpr (isTypeInTuple<FieldType, FIELDPTRS>::value) {
      return FieldPtrs.template getField<FieldType>();
    }
  }
  template <typename FieldType>
  const auto& getField() const {
    if constexpr (isTypeInTuple<FieldType, FIELDS>::value) {
      return Fields.template getField<FieldType>();
    } else if constexpr (isTypeInTuple<FieldType, FIELDPTRS>::value) {
      return FieldPtrs.template getField<FieldType>();
    }
  }

  void Init();

  AbstractConverter<T>& getConverter() { return Conv; }

  void InitComm();
  void InitAverComm();
  void InitIntpComm();

  // get block lattice with block id
  BlockLattice<T, LatSet, ALLFIELDS>& findBlockLat(int blockid) {
    for (BlockLattice<T, LatSet, ALLFIELDS>& blocklat : BlockLats) {
      if (blocklat.getGeo().getBlockId() == blockid) return blocklat;
    }
    std::cerr << "[BlockLatticeManager]: can't find blocklat with blockid " << blockid
              << std::endl;
    exit(1);
  }

  inline std::uint8_t getMaxLevel() const { return BlockGeo.getMaxLevel(); }

  BlockLattice<T, LatSet, ALLFIELDS>& getBlockLat(int i) { return BlockLats[i]; }
  const BlockLattice<T, LatSet, ALLFIELDS>& getBlockLat(int i) const {
    return BlockLats[i];
  }

  std::vector<BlockLattice<T, LatSet, ALLFIELDS>>& getBlockLats() { return BlockLats; }
  const std::vector<BlockLattice<T, LatSet, ALLFIELDS>>& getBlockLats() const {
    return BlockLats;
  }

  BlockGeometry<T, LatSet::d>& getGeo() { return BlockGeo; }
  const BlockGeometry<T, LatSet::d>& getGeo() const { return BlockGeo; }

  void Stream(std::int64_t count);

  template <typename CELLDYNAMICS, typename FieldType>
  void ApplyCellDynamics(std::int64_t count,
                         const BlockFieldManager<FieldType, T, LatSet::d>& BFM);

#ifdef MPI_ENABLED
  void MPIAverComm(std::int64_t count);
  void MPIInterpComm(std::int64_t count);
#endif

  void Communicate(std::int64_t count);

  // tolerance
  void EnableToleranceRho(T rhores = T(1e-5));
  void EnableToleranceU(T ures = T(1e-5));
  T getToleranceRho(int shift = 0);
  T getToleranceU(int shift = 0);
};

// dynamic block lattice, refine and coarsen based on gradient of rho
template <typename T, typename LatSet, typename TypePack>
class DynamicBlockLatticeHelper2D {
 private:
  BlockLatticeManager<T, LatSet, TypePack>& BlockLatMan;
  BlockGeometry2D<T>& BlockGeo;
  BlockGeometryHelper2D<T>& BlockGeoHelper;
  // square norm of gradient of rho field, each BlockCell corresponds to a scalar field
  std::vector<ScalarField<T>> _GradNorm2s;
  // lattice refine threshold
  std::vector<T> _RefineTholds;
  // lattice coarsen threshold
  std::vector<T> _CoarsenTholds;
  // uer defined max refine level
  int _MaxRefineLevel;

  std::vector<T> _MaxGradNorm2s;
  // ScalarField<T> _GradNorm2F;

  BlockFieldManager<VectorFieldAOS<T, 2>, T, 2>& VelocityFM;


 public:
  DynamicBlockLatticeHelper2D(BlockLatticeManager<T, LatSet, TypePack>& blocklatman,
                              BlockGeometryHelper2D<T>& geohelper,
                              BlockFieldManager<VectorFieldAOS<T, 2>, T, 2>& velocityfm,
                              const std::vector<T>& refineth,
                              const std::vector<T>& coarsenth, int MaxRefineLevel = 2)
      : BlockLatMan(blocklatman), BlockGeo(blocklatman.getGeo()),
        BlockGeoHelper(geohelper), VelocityFM(velocityfm), _RefineTholds(refineth),
        _CoarsenTholds(coarsenth), _MaxRefineLevel(MaxRefineLevel)
  // ,_GradNorm2F(BlockGeo.getBaseBlock().getN(), T(0)) {
  {
    // int minsize = std::min(_RefineTholds.size(), _CoarsenTholds.size());
    // _MaxRefineLevel = std::min(_MaxRefineLevel, minsize);
    // init gradnorm2
    for (BasicBlock<T, 2>& block : BlockGeoHelper.getBlockCells()) {
      _GradNorm2s.emplace_back(block.getN(), T(0));
      _MaxGradNorm2s.push_back(T(0));
    }
  }

  void ComputeGradNorm2();
  void UpdateMaxGradNorm2();

  bool WillRefineOrCoarsen();

  void GeoRefine(int OptProcNum, int MaxProcNum = -1, bool enforce = true);

  // field data transfer based on GeoHelper, before re-init geometry
  // pop date transfer between blocks with different level should be treated separately
  // other field data like rho and velo transfer should use Init() in BlockFieldManager
  // you should transfer pop data after other field data hs been transferred
  void PopFieldInit();


  void PopConversionFineToCoarse(const ScalarField<T>& RhoF,
                                 const VectorFieldAOS<T, LatSet::d>& VelocityF,
                                 PopulationField<T, LatSet::q>& PopsF,
                                 const BasicBlock<T, 2>& FBaseBlock,
                                 const BasicBlock<T, 2>& CBlock,
                                 const BasicBlock<T, 2>& CBaseBlock, T OmegaF);

  void PopConversionCoarseToFine(const ScalarField<T>& RhoF,
                                 const VectorFieldAOS<T, LatSet::d>& VelocityF,
                                 PopulationField<T, LatSet::q>& PopsF,
                                 const BasicBlock<T, 2>& CBaseBlock,
                                 const BasicBlock<T, 2>& FBlock,
                                 const BasicBlock<T, 2>& FBaseBlock, T OmegaC);

  // experimental
  std::vector<T>& getMaxGradNorm2s() { return _MaxGradNorm2s; }

  // gather all ScalarField<T> in _GradNorm2s to _GradNorm2F, only used for testing
  // uinform grid void UpdateGradNorm2F(); ScalarField<T>& getGradNorm2F() { return
  // _GradNorm2F; }
};
