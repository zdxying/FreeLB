/* This file is part of FreeLB
 * 
 * Copyright (C) 2024 Yuan Man
 * E-mail contact: ymmanyuan@outlook.com
 * The most recent progress of FreeLB will be updated at
 * <https://github.com/zdxying/FreeLB>
 * 
 * FreeLB is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * FreeLB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with FreeLB. If not, see
 * <https://www.gnu.org/licenses/>.
 * 
 */

// block_lattice.h

#pragma once

#include "data_struct/lattice.h"
#include "utils/fdm_solver.h"

template <typename T, typename LatSet>
struct BlockLatCommStru {
  BlockLattice<T, LatSet>* SendBlock;
  BlockCommStru<T, LatSet::d>* Comm;

  BlockLatCommStru(BlockLattice<T, LatSet>* sblock, BlockCommStru<T, LatSet::d>* blockcomm)
      : SendBlock(sblock), Comm(blockcomm) {}

  std::vector<std::size_t>& getSends() { return Comm->SendCells; }
  std::vector<std::size_t>& getRecvs() { return Comm->RecvCells; }
};

template <typename T, typename LatSet>
struct InterpBlockLatCommStru {
  BlockLattice<T, LatSet>* SendBlock;
  InterpBlockCommStru<T, LatSet::d>* Comm;

  InterpBlockLatCommStru(BlockLattice<T, LatSet>* sblock,
                         InterpBlockCommStru<T, LatSet::d>* interpcomm)
      : SendBlock(sblock), Comm(interpcomm) {}

  std::vector<std::size_t>& getRecvs() { return Comm->RecvCells; }
  std::vector<InterpSource<LatSet::d>>& getSends() { return Comm->SendCells; }
  std::vector<InterpWeight<T, LatSet::d>>& getWeights() { return Comm->InterpWeights; }
};

// TODO: use a unified geometry reference like BasicBlock

// block structure for refined lattice
template <typename T, typename LatSet>
class BlockLattice : public RhoLattice<T> {
 protected:
  // nbr index
  std::array<int, LatSet::q> Delta_Index;
  // geometry
  Block<T, LatSet::d>& BlockGeometry;
  // populations
  PopulationField<T, LatSet::q> Pops;
  // velocity field
  VectorFieldAOS<T, LatSet::d>& Velocity;
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
  BlockLattice(Block<T, LatSet::d>& blockgeo, AbstractConverter<T>& conv,
               VectorFieldAOS<T, LatSet::d>& velocity);

  void InitPop(int Id, T rho) {
    for (int i = 0; i < LatSet::q; ++i) Pops.getField(i)[Id] = rho * LatSet::w[i];
  }
  std::array<T*, LatSet::q> getPop(std::size_t id) { return Pops.template getArray<T>(id); }
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
  Block<T, LatSet::d>& getGeo() { return BlockGeometry; }
  const Block<T, LatSet::d>& getGeo() const { return BlockGeometry; }
  PopulationField<T, LatSet::q>& getPopField() { return Pops; }
  VectorFieldAOS<T, LatSet::d>& getVelocityField() { return Velocity; }

  int getNx() const { return BlockGeometry.getNx(); }
  int getNy() const { return BlockGeometry.getNy(); }
  int getNz() const { return BlockGeometry.getNz(); }
  std::size_t getN() const { return BlockGeometry.getN(); }
  inline T getOmega() const { return Omega; }
  inline T get_Omega() const { return _Omega; }
  inline T getfOmega() const { return fOmega; }
  std::uint8_t getLevel() const { return BlockGeometry.getLevel(); }
  const Vector<int, LatSet::d>& getProjection() const { return BlockGeometry.getProjection(); }
  Vector<T, LatSet::d>& getVelocity(int i) { return Velocity.get(i); }
  std::vector<BlockLatCommStru<T, LatSet>>& getCommunicators() { return Communicators; }
  std::vector<InterpBlockLatCommStru<T, LatSet>>& getAverageComm() { return AverageComm; }
  std::vector<InterpBlockLatCommStru<T, LatSet>>& getInterpComm() { return InterpComm; }
  BlockLatCommStru<T, LatSet>& getCommunicator(int i) { return Communicators[i]; }
  void Refine(bool refinevelo = true, std::uint8_t deltalevel = std::uint8_t(1));
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

  // EXPERIMENTAL
  template <void (*Dynamics)(BCell<T, LatSet>&), typename flagtype>
  void forBlockCells(const GenericArray<flagtype>& flagarr, std::uint8_t flag);

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
  BlockVectFieldAOS<T, LatSet::d>& BlockVelocity;

 public:
  BlockLatticeManager(BlockGeometry<T, LatSet::d>& blockgeo, AbstractConverter<T>& conv,
                      BlockVectFieldAOS<T, LatSet::d>& blockvelocity);

  void InitCommunicators();
  void InitAverComm();
  void InitIntpComm();

  void UpdateMaxLevel() { BlockGeo.UpdateMaxLevel(); }
  inline std::uint8_t getMaxLevel() const { return BlockGeo.getMaxLevel(); }

  std::vector<RhoLattice<T>*> getRhoLattices();

  BlockLattice<T, LatSet>& getBlockLat(int i) { return BlockLats[i]; }
  const BlockLattice<T, LatSet>& getBlockLat(int i) const { return BlockLats[i]; }

  std::vector<BlockLattice<T, LatSet>>& getBlockLats() { return BlockLats; }
  const std::vector<BlockLattice<T, LatSet>>& getBlockLats() const { return BlockLats; }

  BlockGeometry<T, LatSet::d>& getGeo() { return BlockGeo; }
  const BlockGeometry<T, LatSet::d>& getGeo() const { return BlockGeo; }

  BlockVectFieldAOS<T, LatSet::d>& getBlockVelocity() { return BlockVelocity; }
  const BlockVectFieldAOS<T, LatSet::d>& getBlockVelocity() const { return BlockVelocity; }

  // get all rho fields
  std::vector<ScalerField<T>*> getRhoField();
  std::vector<PopulationField<T, LatSet::q>*> getPopField();

  template <typename flagtype>
  void UpdateRho(std::int64_t count, std::vector<ScalerField<flagtype>*>& field, std::uint8_t flag);

  void UpdateRho(std::int64_t count, std::uint8_t flag);

  template <typename flagtype>
  void UpdateRho_Source(std::int64_t count, std::vector<ScalerField<flagtype>*>& field,
                        std::uint8_t flag, std::vector<ScalerField<T>*>& source);

  template <typename flagtype>
  void UpdateU(std::int64_t count, std::vector<ScalerField<flagtype>*>& field, std::uint8_t flag);

  void UpdateU(std::int64_t count, std::uint8_t flag);

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
            typename flagtype = std::uint8_t>
  void BGK(std::int64_t count, std::vector<ScalerField<flagtype>*>& field, std::uint8_t flag);

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
            typename flagtype = std::uint8_t>
  void BGK(std::int64_t count, std::uint8_t flag);

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T),
            typename flagtype = std::uint8_t>
  void BGK_Source(std::int64_t count, std::vector<ScalerField<flagtype>*>& field, std::uint8_t flag,
                  std::vector<ScalerField<T>*>& source);

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
  std::vector<T> _RefineThresholds;
  // lattice coarsen threshold
  std::vector<T> _CoarsenThresholds;
  // uer defined max refine level
  std::uint8_t _MaxRefineLevel;

  // experimental
  std::vector<T> _MaxGradNorm2s;
  // ScalerField<T> _GradNorm2F;


 public:
  DynamicBlockLatticeHelper2D(BlockLatticeManager<T, LatSet>& blocklatman,
                              BlockGeometryHelper2D<T>& geohelper, const std::vector<T>& refineth,
                              const std::vector<T>& coarsenth)
      : BlockLatMan(blocklatman), BlockGeo(blocklatman.getGeo()), BlockGeoHelper(geohelper),
        _RefineThresholds(refineth), _CoarsenThresholds(coarsenth)
  // ,_GradNorm2F(BlockGeo.getBaseBlock().getN(), T(0)) {
  {
    // init gradnorm2
    for (BasicBlock<T, 2>& block : BlockGeoHelper.getBlockCells()) {
      _GradNorm2s.emplace_back(block.getN(), T(0));
      _MaxGradNorm2s.push_back(T(0));
    }
  }

  void ComputeGradNorm2();
  void UpdateMaxGradNorm2();

  // refine and coarsen based on gradient of rho
  void GeoRefineAndCoarsen(int OptProcNum, int MaxProcNum = -1, bool enforce = true);
  // set lattice based on GeoRefineAndCoarsen
  void LatticeRefineAndCoarsen();

  void LatticeRefine(int blockid) {}

  // experimental
  std::vector<T>& getMaxGradNorm2s() { return _MaxGradNorm2s; }

  // gather all ScalerField<T> in _GradNorm2s to _GradNorm2F, only used for testing uinform grid
  // void UpdateGradNorm2F();
  // ScalerField<T>& getGradNorm2F() { return _GradNorm2F; }
};
