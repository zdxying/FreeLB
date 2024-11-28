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

#pragma once

#include "boundary/basic_boundary.h"
#include "lbm/moment.h"


template <typename T, typename LatSet>
class FixedPeriodicBoundary final : public AbstractBoundary {
 private:
  // virtual cells
  std::vector<std::size_t> VCells;
  // virtual cells get info from Real cells
  std::vector<std::size_t> RCells;
  // boundary cells for macroscopic update
  std::vector<std::size_t> BdCells;
  // lattice
  PopLattice<T, LatSet> &Lat;
  // geometry
  Geometry<T, LatSet::d> &Geo;
  // geometry flag
  GenericArray<std::uint8_t> &Field;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;
  // AABB
  AABB<T, LatSet::d> &Box;
  // from AABB
  AABB<T, LatSet::d> &FromBox;
  std::string _name = "FixedPeriodic";

 public:
  FixedPeriodicBoundary(PopLattice<T, LatSet> &lat, AABB<T, LatSet::d> &box,
                        AABB<T, LatSet::d> &frombox, std::uint8_t cellflag,
                        std::uint8_t voidflag = std::uint8_t(1))
      : Lat(lat), Geo(lat.getGeo()), Field(lat.getGeo().getGeoFlagField().getField()),
        Box(box), FromBox(frombox), BdCellFlag(cellflag), voidFlag(voidflag) {
    Vector<T, LatSet::d> boxext = Box.getExtension();
    Vector<int, 3> SIZE(1);
    if constexpr (LatSet::d == 2) {
      SIZE[0] = int(boxext[0] / Geo.getVoxelSize()) + 2;
      SIZE[1] = int(boxext[1] / Geo.getVoxelSize()) + 2;
    } else if constexpr (LatSet::d == 3) {
      SIZE[0] = int(boxext[0] / Geo.getVoxelSize()) + 2;
      SIZE[1] = int(boxext[1] / Geo.getVoxelSize()) + 2;
      SIZE[2] = int(boxext[2] / Geo.getVoxelSize()) + 2;
    }
    int size = SIZE[0] * SIZE[1] * SIZE[2];
    VCells.reserve(size);
    RCells.reserve(size);
    BdCells.reserve(size);
    Geo.forEachVoxel(Box, BdCellFlag, [this](int id) { BdCells.push_back(id); });
    // dist between two AABB
    Vector<T, LatSet::d> dist = FromBox.getCenter() - Box.getCenter();
    // get corrected dist
    for (int i = 0; i < LatSet::d; ++i) {
      if (dist[i] < T(-(1e-6))) {
        dist[i] -= Geo.getVoxelSize();
      } else if (dist[i] > T(1e-6)) {
        dist[i] += Geo.getVoxelSize();
      }
    }
    // extended box
    AABB<T, LatSet::d> extBox = Box.getExtended(Geo.getVoxelSize());
    // get virtual cells and corresponding real cells
    Geo.forEachVoxel(extBox, voidFlag, [this, &dist](int id) {
      if (Geo.template hasNeighborFlag<LatSet>(id, BdCellFlag)) {
        VCells.push_back(id);
        int idr = Geo.findCellId(Geo.getVoxel(id) + dist);
        RCells.push_back(idr);
      }
    });
  }

  void Apply() override {
    int size = VCells.size();
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (int i = 0; i < size; ++i) {
      BasicPopCell<T, LatSet> vcell(VCells[i], Lat);
      BasicPopCell<T, LatSet> rcell(RCells[i], Lat);
      for (int k = 1; k < LatSet::q; ++k) vcell[k] = rcell[k];
    }
  }
  void getinfo() override {
    std::cout << std::setw(18) << std::left << _name << std::setw(10) << std::left << BdCells.size()
              << std::endl;
  }
  void UpdateRho() override {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (std::size_t id : BdCells) {
      BasicPopCell<T, LatSet> cell(id, Lat);
      moment::Rho<T, LatSet>::apply(cell, Lat.getRho(id));
    }
  }
  void UpdateU() override {
#pragma omp parallel for num_threads(Thread_Num) schedule(static)
    for (std::size_t id : BdCells) {
      BasicPopCell<T, LatSet> cell(id, Lat);
      moment::Velocity<T, LatSet>::apply(cell, Lat.getVelocity(id));
    }
  }
};


// --------- block boundary ---------

template <typename BLOCKLATTICE, typename ArrayType>
class BlockFixedPeriodicBoundary : public NonLocalBoundary<BLOCKLATTICE, ArrayType> {
 public:
  BlockFixedPeriodicBoundary(BLOCKLATTICE &lat, const ArrayType &f, std::uint8_t cellflag,
                           std::uint8_t voidflag)
      : NonLocalBoundary<BLOCKLATTICE, ArrayType>(lat, f, cellflag, voidflag) {}
};

template <typename BLOCKLATTICEMANAGER, typename BLOCKFIELDMANAGER>
class FixedPeriodicBoundaryManager {
 public:
  using BLOCKLATTICE = typename BLOCKLATTICEMANAGER::BLOCKLATTICE;
  using CELL = typename BLOCKLATTICE::CellType;
  using T = typename BLOCKLATTICE::FloatType;
  using LatSet = typename BLOCKLATTICE::LatticeSet;
  static constexpr unsigned int D = BLOCKLATTICE::LatticeSet::d;
  using ArrayType = typename BLOCKFIELDMANAGER::array_type;
  using TypePack = typename BLOCKLATTICE::FieldTypePack;

 private:
  std::string _name;
  std::vector<BlockFixedPeriodicBoundary<BLOCKLATTICE, ArrayType>> BdBlocks;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;

  BLOCKLATTICEMANAGER &LatMan;

  BLOCKFIELDMANAGER &BlockFManager;

 public:
  FixedPeriodicBoundaryManager(std::string name, BLOCKLATTICEMANAGER &lat,
                            BLOCKFIELDMANAGER &BlockFM, std::uint8_t cellflag,
                            std::uint8_t voidflag = std::uint8_t(1))
      : _name(name), BdCellFlag(cellflag), voidFlag(voidflag), LatMan(lat),
        BlockFManager(BlockFM) {
    Init();
  }

  void Init() {
    BdBlocks.clear();
    // for each blocks in blocklat
    for (int i = 0; i < LatMan.getGeo().getBlockNum(); ++i) {
      BdBlocks.emplace_back(LatMan.getBlockLat(i),
                            BlockFManager.getBlockField(i).getField(0), BdCellFlag,
                            voidFlag);
    }
  }

  void Setup(const AABB<T, D>& box0, NbrDirection dir0, 
  const AABB<T, D>& box1, NbrDirection dir1) {
    // set periodic boundary condition for cells within box0 and box1
    for (std::size_t i = 0; i < BdBlocks.size(); ++i) {
      auto &bdBlock = BdBlocks[i];
      // find which block contain the 2 boxes
      if (bdBlock.getBdCells().size() > 0){
        BasicCommSet<T, D>& BaseCommSet = bdBlock.getBaseCommSet();
        LatticeCommSet<LatSet, TypePack>& LatCommSet = bdBlock.getLatticeComm();

        auto &lat = bdBlock.getLat();
        auto &block = bdBlock.getLat().getGeo();
        auto &baseblock = block.getBaseBlock();
        // find which box is in the base block
        bool box0in = isOverlapped(baseblock, box0);
        bool box1in = isOverlapped(baseblock, box1);

        if (box0in && box1in) {
          // 2 boxes are in the same block
          BaseCommSet.Recvs.emplace_back(&block);
          block.getCellIdx(box0, BaseCommSet.Recvs.back().Cells);
          LatCommSet.Recvs.emplace_back(lat, BaseCommSet.Recvs.back());
          getReconstPopDir<LatSet>(dir0, LatCommSet.Recvs.back().StreamDirections);

          BaseCommSet.Sends.emplace_back(&block);
          block.getCellIdx(box1, BaseCommSet.Sends.back().Cells);
          LatCommSet.Sends.emplace_back(lat, BaseCommSet.Sends.back());
          getReconstPopDir<LatSet>(dir0, LatCommSet.Sends.back().StreamDirections);

          BaseCommSet.Recvs.emplace_back(&block);
          block.getCellIdx(box1, BaseCommSet.Recvs.back().Cells);
          LatCommSet.Recvs.emplace_back(lat, BaseCommSet.Recvs.back());
          getReconstPopDir<LatSet>(dir1, LatCommSet.Recvs.back().StreamDirections);

          BaseCommSet.Sends.emplace_back(&block);
          block.getCellIdx(box0, BaseCommSet.Sends.back().Cells);
          LatCommSet.Sends.emplace_back(lat, BaseCommSet.Sends.back());
          getReconstPopDir<LatSet>(dir1, LatCommSet.Sends.back().StreamDirections);

        } else if (box0in) {
          // box0 is in the block
          // find which block contains box1
          for (std::size_t j = i + 1; j < BdBlocks.size(); ++j) {
            auto &bdBlock1 = BdBlocks[j];
            if (bdBlock1.getBdCells().size() > 0){
              BasicCommSet<T, D>& BaseCommSet1 = bdBlock1.getBaseCommSet();
              LatticeCommSet<LatSet, TypePack>& LatCommSet1 = bdBlock1.getLatticeComm();
              auto &lat1 = bdBlock1.getLat();
              auto &block1 = bdBlock1.getLat().getGeo();
              auto &baseblock1 = block1.getBaseBlock();
              if (isOverlapped(baseblock1, box1)) {
                BaseCommSet1.Sends.emplace_back(&block);
                block1.getCellIdx(box1, BaseCommSet1.Sends.back().Cells);
                LatCommSet1.Sends.emplace_back(lat, BaseCommSet1.Sends.back());
                getReconstPopDir<LatSet>(dir0, LatCommSet1.Sends.back().StreamDirections);

                BaseCommSet1.Recvs.emplace_back(&block);
                block1.getCellIdx(box1, BaseCommSet1.Recvs.back().Cells);
                LatCommSet1.Recvs.emplace_back(lat, BaseCommSet1.Recvs.back());
                getReconstPopDir<LatSet>(dir1, LatCommSet1.Recvs.back().StreamDirections);

                BaseCommSet.Sends.emplace_back(&block1);
                block.getCellIdx(box0, BaseCommSet.Sends.back().Cells);
                LatCommSet.Sends.emplace_back(lat1, BaseCommSet.Sends.back());
                getReconstPopDir<LatSet>(dir1, LatCommSet.Sends.back().StreamDirections);

                BaseCommSet.Recvs.emplace_back(&block1);
                block.getCellIdx(box0, BaseCommSet.Recvs.back().Cells);
                LatCommSet.Recvs.emplace_back(lat1, BaseCommSet.Recvs.back());
                getReconstPopDir<LatSet>(dir0, LatCommSet.Recvs.back().StreamDirections);

                i = j;
                break;
              }
            }
          }
        } else if (box1in) {
          // box1 is in the block
          // find which block contains box0
          for (std::size_t j = i + 1; j < BdBlocks.size(); ++j) {
            auto &bdBlock0 = BdBlocks[j];
            if (bdBlock0.getBdCells().size() > 0){
              BasicCommSet<T, D>& BaseCommSet0 = bdBlock0.getBaseCommSet();
              LatticeCommSet<LatSet, TypePack>& LatCommSet0 = bdBlock0.getLatticeComm();
              auto &lat0 = bdBlock0.getLat();
              auto &block0 = bdBlock0.getLat().getGeo();
              auto &baseblock0 = block0.getBaseBlock();
              if (isOverlapped(baseblock0, box1)) {
                BaseCommSet0.Sends.emplace_back(&block);
                block0.getCellIdx(box0, BaseCommSet0.Sends.back().Cells);
                LatCommSet0.Sends.emplace_back(lat, BaseCommSet0.Sends.back());
                getReconstPopDir<LatSet>(dir1, LatCommSet0.Sends.back().StreamDirections);

                BaseCommSet0.Recvs.emplace_back(&block);
                block0.getCellIdx(box0, BaseCommSet0.Recvs.back().Cells);
                LatCommSet0.Recvs.emplace_back(lat, BaseCommSet0.Recvs.back());
                getReconstPopDir<LatSet>(dir0, LatCommSet0.Recvs.back().StreamDirections);

                BaseCommSet.Sends.emplace_back(&block0);
                block.getCellIdx(box1, BaseCommSet.Sends.back().Cells);
                LatCommSet.Sends.emplace_back(lat0, BaseCommSet.Sends.back());
                getReconstPopDir<LatSet>(dir0, LatCommSet.Sends.back().StreamDirections);

                BaseCommSet.Recvs.emplace_back(&block0);
                block.getCellIdx(box1, BaseCommSet.Recvs.back().Cells);
                LatCommSet.Recvs.emplace_back(lat0, BaseCommSet.Recvs.back());
                getReconstPopDir<LatSet>(dir1, LatCommSet.Recvs.back().StreamDirections);

                i = j;
                break;
              }
            }
          }
        }
      }
    }
  }

  void Apply() {
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
    for (auto &bdBlock : BdBlocks) {
      LatticeCommSet<LatSet, TypePack>& latcomm = bdBlock.getLatticeComm();
      BLOCKLATTICE& recvblocklat = bdBlock.getLat();
      for (const LatticeComm<LatSet, TypePack>& recvcomm : latcomm.Recvs) {
        auto &sendbdBlock = BdBlocks[recvcomm.Comm.TargetBlockId];
        BLOCKLATTICE& sendblocklat = sendbdBlock.getLat();
        LatticeCommSet<LatSet, TypePack>& sendlatcomm = sendbdBlock.getLatticeComm();
        const LatticeComm<LatSet, TypePack>& sendcomm = sendlatcomm.getSendComm(recvblocklat.getGeo().getBlockId());
        for (std::size_t id = 0; id < recvcomm.Comm.Cells.size(); ++id) {
          CELL recvcell(recvcomm.Comm.Cells[id], recvblocklat);
          CELL sendcell(sendcomm.Comm.Cells[id], sendblocklat);
          for (unsigned int k : recvcomm.StreamDirections) {
            recvcell[k] = sendcell[k];
          }
        }
      }
    }
  }
};
