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

template <typename T, typename LatSet>
struct BlockLatCommStru {
  BlockLattice<T, LatSet>* SendBlock;
  BlockComm<T, LatSet::d>* Comm;

  BlockLatCommStru(BlockLattice<T, LatSet>* sblock, BlockComm<T, LatSet::d>* blockcomm)
      : SendBlock(sblock), Comm(blockcomm) {}

  std::vector<std::size_t>& getSends() { return Comm->SendCells; }
  std::vector<std::size_t>& getRecvs() { return Comm->RecvCells; }
};

template <typename T, typename LatSet>
struct InterpBlockLatCommStru {
  BlockLattice<T, LatSet>* SendBlock;
  InterpBlockComm<T, LatSet::d>* Comm;

  InterpBlockLatCommStru(BlockLattice<T, LatSet>* sblock,
                         InterpBlockComm<T, LatSet::d>* interpcomm)
      : SendBlock(sblock), Comm(interpcomm) {}

  std::vector<std::size_t>& getRecvs() { return Comm->RecvCells; }
  std::vector<InterpSource<LatSet::d>>& getSends() { return Comm->SendCells; }
  std::vector<InterpWeight<T, LatSet::d>>& getWeights() { return Comm->InterpWeights; }
};


// a base lattice for BlockLattice
template <typename T>
class BlockRhoLattice {
 protected:
  ScalerField<T>& Rho;
  // converter
  AbstractConverter<T>& Conv;
  // rho init
  T Lattice_Rho_Init;
  // buoyancy
  T Lattice_gbeta;

 public:
  BlockRhoLattice(AbstractConverter<T>& conv, ScalerField<T>& rho)
      : Conv(conv), Lattice_Rho_Init(conv.getLatRhoInit()),
        Lattice_gbeta(conv.GetLattice_gbeta()), Rho(rho) {}

  ScalerField<T>& getRhoField() { return Rho; }
  const ScalerField<T>& getRhoField() const { return Rho; }

  const T& getRho(int i) const { return Rho.get(i); }
  T& getRho(int i) { return Rho.get(i); }

  void SetRhoField(std::size_t id, T value) { Rho.SetField(id, value); }

  T getLatRhoInit() const { return Lattice_Rho_Init; }
  T getLatgBeta() const { return Lattice_gbeta; }
};

// block structure for refined lattice
template <typename T, typename LatSet>
class BlockLattice : public BlockRhoLattice<T> {
 protected:
  // nbr index
  std::array<int, LatSet::q> Delta_Index;
  // geometry
  Block<T, LatSet::d>& BlockGeo;
  // populations
  PopulationField<T, LatSet::q>& Pops;
  // velocity field
  VectorFieldAOS<T, LatSet::d>& Velocity;

  // --- lattice communication structure ---
  // conmmunicate with same level block
  std::vector<BlockLatCommStru<T, LatSet>> Communicators;
  // average blocklat comm, get from higher level block
  std::vector<InterpBlockLatCommStru<T, LatSet>> AverageComm;
  // interp blocklat comm, get from lower level block
  std::vector<InterpBlockLatCommStru<T, LatSet>> InterpComm;

  // omega
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
// MPI
#ifdef MPI_ENABLED
  int _Rank;
  bool _NeedMPIComm;
  // buffers
  MPIBlockBuffer<T> MPIBuffer;
#endif

 public:
  BlockLattice(Block<T, LatSet::d>& block, ScalerField<T>& rho,
               VectorFieldAOS<T, LatSet::d>& velocity,
               PopulationField<T, LatSet::q>& pops, AbstractConverter<T>& conv, bool initpop = true);

  void InitPop(int Id, T rho) {
    for (int i = 0; i < LatSet::q; ++i) Pops.getField(i)[Id] = rho * LatSet::w[i];
  }
  std::array<T*, LatSet::q> getPop(std::size_t id) {
    return Pops.template getArray<T>(id);
  }
  T& getPopdir(std::size_t id, int dir) { return Pops.getField(dir)[id]; }
  const T& getPopdir(std::size_t id, int dir) const { return Pops.getField(dir)[id]; }
  BCell<T, LatSet> getNeighbor(const BCell<T, LatSet>& cell, int i) const {
    return BCell<T, LatSet>(cell.getId() + Delta_Index[i], *this);
  }
  BCell<T, LatSet> getNeighbor(const BCell<T, LatSet>& cell,
                               const Vector<int, LatSet::d>& direction) const {
    return BCell<T, LatSet>(cell.getId() + direction * getProjection(), *this);
  }

  std::size_t getNbrId(std::size_t id, int dir) const { return id + Delta_Index[dir]; }
  const std::array<int, LatSet::q>& getDelta_Index() const { return Delta_Index; }

  Block<T, LatSet::d>& getGeo() { return BlockGeo; }
  const Block<T, LatSet::d>& getGeo() const { return BlockGeo; }

  PopulationField<T, LatSet::q>& getPopField() { return Pops; }
  VectorFieldAOS<T, LatSet::d>& getVelocityField() { return Velocity; }

  int getNx() const { return BlockGeo.getNx(); }
  int getNy() const { return BlockGeo.getNy(); }
  int getNz() const { return BlockGeo.getNz(); }
  std::size_t getN() const { return BlockGeo.getN(); }
  inline T getOmega() const { return Omega; }
  inline T get_Omega() const { return _Omega; }
  inline T getfOmega() const { return fOmega; }
  std::uint8_t getLevel() const { return BlockGeo.getLevel(); }
  const Vector<int, LatSet::d>& getProjection() const { return BlockGeo.getProjection(); }
  Vector<T, LatSet::d>& getVelocity(int i) { return Velocity.get(i); }
  std::vector<BlockLatCommStru<T, LatSet>>& getCommunicators() { return Communicators; }
  std::vector<InterpBlockLatCommStru<T, LatSet>>& getAverageComm() { return AverageComm; }
  std::vector<InterpBlockLatCommStru<T, LatSet>>& getInterpComm() { return InterpComm; }
  BlockLatCommStru<T, LatSet>& getCommunicator(int i) { return Communicators[i]; }

  // communication, before streaming
  void communicate();
  // average communication
  void avercommunicate();
  // interp communication
  void interpcommunicate();
  // update rho based on a flag field of flagtype(usually std::uint8_t or enum
  // of std::uint8_t)
  template <typename flagtype>
  void UpdateRho(const GenericArray<flagtype>& flagarr, std::uint8_t flag);

  template <typename flagtype>
  void UpdateRho_Source(const GenericArray<flagtype>& flagarr, std::uint8_t flag,
                        const GenericArray<T>& source);

  template <typename flagtype>
  void UpdateU(const GenericArray<flagtype>& flagarr, std::uint8_t flag);

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
            typename flagtype = std::uint8_t>
  void BGK(const GenericArray<flagtype>& flagarr, std::uint8_t flag);

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
            typename flagtype = std::uint8_t>
  void BGK_Source(const GenericArray<flagtype>& flagarr, std::uint8_t flag,
                  const GenericArray<T>& source);

  void Stream();

  // tolerance
  void EnableToleranceRho(T rhores = T(1e-5));
  void EnableToleranceU(T ures = T(1e-5));
  T getToleranceRho();
  T getToleranceU();
  // get inner block tolerance
  T getTolRho(int shift = 1);
  T getTolU(int shift = 1);


#ifdef MPI_ENABLED
  int getRank() const { return _Rank; }
  bool getNeedMPIComm() const { return _NeedMPIComm; }

  MPIBlockBuffer<T>& getMPIBlockBuffer() { return MPIBuffer; }
  const MPIBlockBuffer<T>& getMPIBlockBuffer() const { return MPIBuffer; }

  void MPIcommunicate(MPIBlockCommStru& MPIComm);
#endif
};

// block lattice manager
template <typename T, typename LatSet>
class BlockLatticeManager {
 private:
  std::vector<BlockLattice<T, LatSet>> BlockLats;
  BlockGeometry<T, LatSet::d>& BlockGeo;
  BlockFieldManager<VectorFieldAOS<T, LatSet::d>, T, LatSet::d>& VelocityFM;
  AbstractConverter<T>& Conv;

  // Rho Field
  BlockFieldManager<ScalerField<T>, T, LatSet::d> RhoFM;
  // Population Field
  BlockFieldManager<PopulationField<T, LatSet::q>, T, LatSet::d> PopsFM;

 public:
  BlockLatticeManager(
    BlockGeometry<T, LatSet::d>& blockgeo, AbstractConverter<T>& conv,
    BlockFieldManager<VectorFieldAOS<T, LatSet::d>, T, LatSet::d>& blockvelocity);
  
  void Init();

  AbstractConverter<T>& getConverter() { return Conv; }

  void InitCommunicators();
  void InitAverComm();
  void InitIntpComm();

  void UpdateMaxLevel() { BlockGeo.UpdateMaxLevel(); }
  inline std::uint8_t getMaxLevel() const { return BlockGeo.getMaxLevel(); }

  BlockFieldManager<ScalerField<T>, T, LatSet::d>& getRhoFM() { return RhoFM; }
  const BlockFieldManager<ScalerField<T>, T, LatSet::d>& getRhoFM() const {
    return RhoFM;
  }

  BlockFieldManager<PopulationField<T, LatSet::q>, T, LatSet::d>& getPopsFM() {
    return PopsFM;
  }
  const BlockFieldManager<PopulationField<T, LatSet::q>, T, LatSet::d>& getPopsFM()
    const {
    return PopsFM;
  }

  BlockFieldManager<VectorFieldAOS<T, LatSet::d>, T, LatSet::d>& getVelocityFM() {
    return VelocityFM;
  }
  const BlockFieldManager<VectorFieldAOS<T, LatSet::d>, T, LatSet::d>& getVelocityFM()
    const {
    return VelocityFM;
  }

  BlockLattice<T, LatSet>& getBlockLat(int i) { return BlockLats[i]; }
  const BlockLattice<T, LatSet>& getBlockLat(int i) const { return BlockLats[i]; }

  std::vector<BlockLattice<T, LatSet>>& getBlockLats() { return BlockLats; }
  const std::vector<BlockLattice<T, LatSet>>& getBlockLats() const { return BlockLats; }

  BlockGeometry<T, LatSet::d>& getGeo() { return BlockGeo; }
  const BlockGeometry<T, LatSet::d>& getGeo() const { return BlockGeo; }


  template <typename flagtype = std::uint8_t>
  void UpdateRho(std::int64_t count, std::uint8_t flag,
                 const BlockFieldManager<ScalerField<flagtype>, T, LatSet::d>& BFM);

  template <typename flagtype = std::uint8_t>
  void UpdateRho_Source(std::int64_t count, std::uint8_t flag,
                        const BlockFieldManager<ScalerField<flagtype>, T, LatSet::d>& BFM,
                        const BlockFieldManager<ScalerField<T>, T, LatSet::d>& source);

  template <typename flagtype = std::uint8_t>
  void UpdateU(std::int64_t count, std::uint8_t flag,
               const BlockFieldManager<ScalerField<flagtype>, T, LatSet::d>& BFM);

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
            typename flagtype = std::uint8_t>
  void BGK(std::int64_t count, std::uint8_t flag,
           const BlockFieldManager<ScalerField<flagtype>, T, LatSet::d>& BFM);

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
            typename flagtype = std::uint8_t>
  void BGK_Source(std::int64_t count, std::uint8_t flag,
                  const BlockFieldManager<ScalerField<flagtype>, T, LatSet::d>& BFM,
                  const BlockFieldManager<ScalerField<T>, T, LatSet::d>& source);

  void Stream(std::int64_t count);

  void Communicate(std::int64_t count);

  // tolerance
  void EnableToleranceRho(T rhores = T(1e-5));
  void EnableToleranceU(T ures = T(1e-5));
  T getToleranceRho(int shift = 0);
  T getToleranceU(int shift = 0);
};

// dynamic block lattice, refine and coarsen based on gradient of rho
template <typename T, typename LatSet>
class DynamicBlockLatticeHelper2D {
 private:
  BlockLatticeManager<T, LatSet>& BlockLatMan;
  BlockGeometry2D<T>& BlockGeo;
  BlockGeometryHelper2D<T>& BlockGeoHelper;
  // square norm of gradient of rho field, each BlockCell corresponds to a scalar field
  std::vector<ScalerField<T>> _GradNorm2s;
  // lattice refine threshold
  std::vector<T> _RefineTholds;
  // lattice coarsen threshold
  std::vector<T> _CoarsenTholds;
  // uer defined max refine level
  int _MaxRefineLevel;

  std::vector<T> _MaxGradNorm2s;
  // ScalerField<T> _GradNorm2F;

  BlockFieldManager<VectorFieldAOS<T, 2>, T, 2>& VelocityFM;


 public:
  DynamicBlockLatticeHelper2D(BlockLatticeManager<T, LatSet>& blocklatman,
                              BlockGeometryHelper2D<T>& geohelper,
                              BlockFieldManager<VectorFieldAOS<T, 2>, T, 2>& velocityfm,
                              const std::vector<T>& refineth,
                              const std::vector<T>& coarsenth, int MaxRefineLevel = 2)
      : BlockLatMan(blocklatman), BlockGeo(blocklatman.getGeo()),
        BlockGeoHelper(geohelper), VelocityFM(velocityfm), _RefineTholds(refineth),
        _CoarsenTholds(coarsenth), _MaxRefineLevel(MaxRefineLevel)
  // ,_GradNorm2F(BlockGeo.getBaseBlock().getN(), T(0)) {
  {
    int minsize = std::min(_RefineTholds.size(), _CoarsenTholds.size());
    _MaxRefineLevel = std::min(_MaxRefineLevel, minsize);
    // init gradnorm2
    for (BasicBlock<T, 2>& block : BlockGeoHelper.getBlockCells()) {
      _GradNorm2s.emplace_back(block.getN(), T(0));
      _MaxGradNorm2s.push_back(T(0));
    }
  }

  void ComputeGradNorm2();
  void UpdateMaxGradNorm2();

  bool WillRefineOrCoarsen();

  void GeoRefineOrCoarsen(int OptProcNum, int MaxProcNum = -1, bool enforce = true);

  // field data transfer based on GeoHelper, before re-init geometry
  // pop date transfer between blocks with different level should be treated separately
  // other field data like rho and velo transfer should use Init() in BlockFieldManager
  // you should transfer pop data after other field data hs been transferred
  void PopFieldInit();


  void PopConversionFineToCoarse(const ScalerField<T>& RhoF,
                                 const VectorFieldAOS<T, LatSet::d>& VelocityF,
                                 PopulationField<T, LatSet::q>& PopsF,
                                 const BasicBlock<T, 2>& FBaseBlock,
                                 const BasicBlock<T, 2>& CBlock,
                                 const BasicBlock<T, 2>& CBaseBlock, T OmegaF);

  void PopConversionCoarseToFine(const ScalerField<T>& RhoF,
                                 const VectorFieldAOS<T, LatSet::d>& VelocityF,
                                 PopulationField<T, LatSet::q>& PopsF,
                                 const BasicBlock<T, 2>& CBaseBlock,
                                 const BasicBlock<T, 2>& FBlock,
                                 const BasicBlock<T, 2>& FBaseBlock, T OmegaC);

  // experimental
  std::vector<T>& getMaxGradNorm2s() { return _MaxGradNorm2s; }

  // gather all ScalerField<T> in _GradNorm2s to _GradNorm2F, only used for testing
  // uinform grid void UpdateGradNorm2F(); ScalerField<T>& getGradNorm2F() { return
  // _GradNorm2F; }
};
