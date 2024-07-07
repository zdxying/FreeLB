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

#include "data_struct/block_lattice_base.h"
// AbstractConverter
#include "lbm/unit_converter.h"

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

// block structure for refined lattice
template <typename T, typename LatSet, typename TypePack>
class BlockLattice : public BlockLatticeBase<T, LatSet, TypePack> {
 public:
  using CellType = Cell<T, LatSet, TypePack>;
  using LatticeSet = LatSet;
  using FloatType = T;

  using GenericRho = typename FindGenericRhoType<T, TypePack>::type;

 protected:
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
  template <typename... FIELDPTRS>
  BlockLattice(Block<T, LatSet::d>& block, AbstractConverter<T>& conv,
               std::tuple<FIELDPTRS...> fieldptrs);

  std::array<T*, LatSet::q> getPop(std::size_t id) {
    return this->template getField<POP<T, LatSet::q>>().getArray(id);
  }

  inline T getOmega() const { return Omega; }
  inline T get_Omega() const { return _Omega; }
  inline T getfOmega() const { return fOmega; }

  std::vector<BlockLatComm<T, LatSet, TypePack>>& getCommunicators() {
    return Communicators;
  }
  std::vector<InterpBlockLatComm<T, LatSet, TypePack>>& getAverageComm() {
    return AverageComm;
  }
  std::vector<InterpBlockLatComm<T, LatSet, TypePack>>& getInterpComm() {
    return InterpComm;
  }

  // normal communication, which can be done using normalcommunicate() in blockFM
  void communicate();
  // average communication
  void avercommunicate();
  // interp communication
  void interpcommunicate();

  void Stream();

  template <typename CELLDYNAMICS, typename ArrayType>
  void ApplyCellDynamics(const ArrayType& flagarr);

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics();

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics(const Genericvector<std::size_t>& Idx);

  template <typename DYNAMICS, typename elementType>
  void ApplyDynamics(const Genericvector<elementType>& Idx);

  template <typename CELLDYNAMICS, typename ArrayType>
  void ApplyInnerCellDynamics(const ArrayType& flagarr);

  template <typename CELLDYNAMICS>
  void ApplyInnerCellDynamics();

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
class BlockLatticeManager : public BlockLatticeManagerBase<T, LatSet, TypePack> {
 public:
  using FIELDS = typename ExtractFieldPack<TypePack>::pack1;
  using FIELDPTRS = typename ExtractFieldPack<TypePack>::pack2;
  using ALLFIELDS = typename ExtractFieldPack<TypePack>::mergedpack;

  using BLOCKLATTICE = BlockLattice<T, LatSet, ALLFIELDS>;
  using CellType = Cell<T, LatSet, ALLFIELDS>;
  using FloatType = T;
  using LatticeSet = LatSet;

  using GenericRho = typename GetGenericRhoType<T, FIELDS>::type;

 private:
  std::vector<BlockLattice<T, LatSet, ALLFIELDS>> BlockLats;
  AbstractConverter<T>& Conv;

 public:
  template <typename... FIELDPTRTYPES>
  BlockLatticeManager(BlockGeometry<T, LatSet::d>& blockgeo, AbstractConverter<T>& conv,
                      FIELDPTRTYPES*... fieldptrs);
  template <typename INITVALUEPACK, typename... FIELDPTRTYPES>
  BlockLatticeManager(BlockGeometry<T, LatSet::d>& blockgeo, INITVALUEPACK& initvalues,
                      AbstractConverter<T>& conv, FIELDPTRTYPES*... fieldptrs);

  void Init();
  template <typename... InitValues>
  void Init(BlockGeometryHelper<T, LatSet::d>& GeoHelper,
            const std::tuple<InitValues...>& initvalues);

  template <typename Func>
  void ForEachField(Func&& func) {
    this->Fields.ForEachField(std::forward<Func>(func));
  }

  template <typename FieldType>
  void addField(BlockFieldManager<FieldType, T, LatSet::d>& field) {
    this->FieldPtrs.template addField<FieldType>(field);
    // constexpr std::size_t ith = this->FieldPtrs.is_at_index<FieldType>() +
    // FIELDS::size;
    int i = 0;
    for (auto& blocklat : BlockLats) {
      // blocklat.template addField<FieldType, ith>(field);
      blocklat.template addField<FieldType>(field.getBlockField(i));
      ++i;
    }
  }

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

  inline std::uint8_t getMaxLevel() const { return this->BlockGeo.getMaxLevel(); }

  BlockLattice<T, LatSet, ALLFIELDS>& getBlockLat(int i) { return BlockLats[i]; }
  const BlockLattice<T, LatSet, ALLFIELDS>& getBlockLat(int i) const {
    return BlockLats[i];
  }
  std::vector<BlockLattice<T, LatSet, ALLFIELDS>>& getBlockLats() { return BlockLats; }
  const std::vector<BlockLattice<T, LatSet, ALLFIELDS>>& getBlockLats() const {
    return BlockLats;
  }

  void Stream(std::int64_t count);
  void Stream();

  template <typename CELLDYNAMICS, typename FieldType>
  void ApplyCellDynamics(std::int64_t count,
                         const BlockFieldManager<FieldType, T, LatSet::d>& BFM);
  template <typename CELLDYNAMICS, typename FieldType>
  void ApplyCellDynamics(const BlockFieldManager<FieldType, T, LatSet::d>& BFM);

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics(std::int64_t count);
  template <typename CELLDYNAMICS>
  void ApplyCellDynamics();

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics(std::int64_t count, const GenericvectorManager<std::size_t>& blockids);
  template <typename CELLDYNAMICS>
  void ApplyCellDynamics(const GenericvectorManager<std::size_t>& blockids);

  template <typename DYNAMICS, typename elementType>
  void ApplyDynamics(std::int64_t count, const GenericvectorManager<elementType>& blockids);
  template <typename DYNAMICS, typename elementType>
  void ApplyDynamics(const GenericvectorManager<elementType>& blockids);


  template <typename CELLDYNAMICS, typename FieldType>
  void ApplyInnerCellDynamics(std::int64_t count,
                              const BlockFieldManager<FieldType, T, LatSet::d>& BFM);
  template <typename CELLDYNAMICS, typename FieldType>
  void ApplyInnerCellDynamics(const BlockFieldManager<FieldType, T, LatSet::d>& BFM);

  template <typename CELLDYNAMICS>
  void ApplyInnerCellDynamics(std::int64_t count);
  template <typename CELLDYNAMICS>
  void ApplyInnerCellDynamics();

  template <typename Func>
  void ForEachBlockLattice(Func&& func) {
    for (auto& blocklat : BlockLats) {
      func(blocklat);
    }
  }

  template <typename Func>
  void ForEachBlockLattice(std::int64_t count, Func&& func) {
    for (auto& blocklat : BlockLats) {
      const int deLevel = static_cast<int>(getMaxLevel() - blocklat.getLevel());
      if (count % (static_cast<int>(pow(2, deLevel))) == 0) {
        func(blocklat);
      }
    }
  }


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

// coupling block lattice manager
template <typename BlockLatManager0, typename BlockLatManager1>
class BlockLatManagerCoupling {
 public:
  using T = typename BlockLatManager0::FloatType;
  using CELL0 = typename BlockLatManager0::CellType;
  using CELL1 = typename BlockLatManager1::CellType;

  BlockLatManagerCoupling(BlockLatManager0& blocklatman0, BlockLatManager1& blocklatman1)
      : BlockLatMan0(blocklatman0), BlockLatMan1(blocklatman1) {}

  template <typename CELLDYNAMICS, typename BLOCKFIELDMANAGER>
  void ApplyCellDynamics(std::int64_t count, const BLOCKFIELDMANAGER& BFM) {
    int size = static_cast<int>(BlockLatMan0.getBlockLats().size());
#pragma omp parallel for num_threads(Thread_Num)
    for (int i = 0; i < size; ++i) {
      auto& blocklat0 = BlockLatMan0.getBlockLat(i);
      auto& blocklat1 = BlockLatMan1.getBlockLat(i);
      const int deLevel =
        static_cast<int>(BlockLatMan0.getMaxLevel() - blocklat0.getLevel());
      if (count % (static_cast<int>(pow(2, deLevel))) == 0) {
        const auto& flagArray = BFM.getBlockField(i).getField(0);
        CELL0 cell0(0, blocklat0);
        CELL1 cell1(0, blocklat1);
        for (std::size_t id = 0; id < blocklat0.getN(); ++id) {
          cell0.setId(id);
          cell1.setId(id);
          CELLDYNAMICS::Execute(flagArray[id], cell0, cell1);
        }
      }
    }
  }

 private:
  BlockLatManager0& BlockLatMan0;
  BlockLatManager1& BlockLatMan1;
};


// dynamic block lattice, refine and coarsen based on gradient of rho
template <typename T, typename LatSet, typename TypePack>
class DynamicBlockLatticeHelper2D {
 public:
  using FIELDS = typename ExtractFieldPack<TypePack>::pack1;
  using GenericRho = typename GetGenericRhoType<T, FIELDS>::type;

  using FIELDPTRS = typename ExtractFieldPack<TypePack>::pack2;
  using ALLFIELDS = typename ExtractFieldPack<TypePack>::mergedpack;

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


 public:
  DynamicBlockLatticeHelper2D(BlockLatticeManager<T, LatSet, TypePack>& blocklatman,
                              BlockGeometryHelper2D<T>& geohelper,
                              const std::vector<T>& refineth,
                              const std::vector<T>& coarsenth, int MaxRefineLevel = 2)
      : BlockLatMan(blocklatman), BlockGeo(blocklatman.getGeo()),
        BlockGeoHelper(geohelper), _RefineTholds(refineth), _CoarsenTholds(coarsenth),
        _MaxRefineLevel(MaxRefineLevel)
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
