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

 private:
  // A periodic copy runs before Stream(), and Stream() pulls: the population
  // that must land on a face cell is read from the neighbour on the far side
  // of the face.  Writing the face cell itself therefore lands the copy one
  // cell inside the domain, so the recv cells are moved one voxel across every
  // periodic axis into the block's own halo, which is where the pull reads.
  // The boxes are periodic domain faces, so that halo cell is padding.
  static void MoveToHalo(const Block<T, D>& cellsFrom, NbrDirection dir,
                         std::vector<std::size_t>& cells) {
    const Vector<int, D>& proj = cellsFrom.getProjection();
    const int overlap = cellsFrom.getOverlap();
    const Vector<int, D>& mesh = cellsFrom.getBaseBlock().getMesh();
    std::ptrdiff_t offset = 0;
    int step[D] = {};
    for (unsigned int a = 0; a < D; ++a) {
      const std::uint8_t neg = std::uint8_t(1) << (2 * a);
      const std::uint8_t pos = std::uint8_t(2) << (2 * a);
      if (util::isFlag(static_cast<std::uint8_t>(dir), neg)) {
        step[a] = -1;
      } else if (util::isFlag(static_cast<std::uint8_t>(dir), pos)) {
        step[a] = 1;
      }
      offset += static_cast<std::ptrdiff_t>(step[a]) * proj[a];
    }
    if (offset == 0) return;
    // The halo written below is the block's own padding.  A box that is not the
    // block's boundary layer would make the copy overwrite a real neighbour
    // cell, and a block without an overlap has no padding at all; both are
    // caller errors, so fail here instead of corrupting the domain.
    if (overlap < 1) {
      std::cerr << "[FixedPeriodicBoundaryManager]: a periodic face needs a "
                   "block overlap"
                << std::endl;
      std::exit(-1);
    }
    for (std::size_t& id : cells) {
      const Vector<int, D> loc = cellsFrom.getLoc(id);
      for (unsigned int a = 0; a < D; ++a) {
        if (step[a] == 0) continue;
        const int boundary = step[a] < 0 ? overlap : overlap + mesh[a] - 1;
        if (loc[a] != boundary) {
          std::cerr << "[FixedPeriodicBoundaryManager]: periodic face box is not "
                       "the block boundary (axis "
                    << a << ", local " << loc[a] << ", expected " << boundary
                    << ")" << std::endl;
          std::exit(-1);
        }
      }
      id += offset;
    }
  }

  std::size_t AddRecv(BasicCommSet<T, D>& base, LatticeCommSet<LatSet, TypePack>& latc,
                      BLOCKLATTICE& lat, Block<T, D>& target, Block<T, D>& cellsFrom,
                      const AABB<T, D>& box, NbrDirection dir) {
    base.Recvs.emplace_back(&target);
    cellsFrom.getCellIdx(box, base.Recvs.back().Cells);
    MoveToHalo(cellsFrom, dir, base.Recvs.back().Cells);
    latc.Recvs.emplace_back(lat, base.Recvs.back());
    getPeriodicPopDir<LatSet>(dir, latc.Recvs.back().StreamDirections);
    return latc.Recvs.size() - 1;
  }

  std::size_t AddSend(BasicCommSet<T, D>& base, LatticeCommSet<LatSet, TypePack>& latc,
                      BLOCKLATTICE& lat, Block<T, D>& target, Block<T, D>& cellsFrom,
                      const AABB<T, D>& box, NbrDirection dir) {
    base.Sends.emplace_back(&target);
    cellsFrom.getCellIdx(box, base.Sends.back().Cells);
    latc.Sends.emplace_back(lat, base.Sends.back());
    getPeriodicPopDir<LatSet>(dir, latc.Sends.back().StreamDirections);
    return latc.Sends.size() - 1;
  }

 public:
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
        auto &block = bdBlock.getLat().getBlock();
        auto &baseblock = block.getBaseBlock();
        // find which box is in the base block
        bool box0in = isOverlapped(baseblock, box0);
        bool box1in = isOverlapped(baseblock, box1);

        if (box0in && box1in) {
          // 2 boxes are in the same block.  Each recv is paired with the send
          // holding its opposite face; a block with two periodic pairs has
          // several sends to the same block id, so the id cannot identify it.
          const std::size_t r0 =
            AddRecv(BaseCommSet, LatCommSet, lat, block, block, box0, dir0);
          const std::size_t s0 =
            AddSend(BaseCommSet, LatCommSet, lat, block, block, box1, dir0);
          LatCommSet.Recvs[r0].SendCommIndex = static_cast<int>(s0);

          const std::size_t r1 =
            AddRecv(BaseCommSet, LatCommSet, lat, block, block, box1, dir1);
          const std::size_t s1 =
            AddSend(BaseCommSet, LatCommSet, lat, block, block, box0, dir1);
          LatCommSet.Recvs[r1].SendCommIndex = static_cast<int>(s1);

        } else if (box0in) {
          // box0 is in the block
          // find which block contains box1
          for (std::size_t j = i + 1; j < BdBlocks.size(); ++j) {
            auto &bdBlock1 = BdBlocks[j];
            if (bdBlock1.getBdCells().size() > 0){
              BasicCommSet<T, D>& BaseCommSet1 = bdBlock1.getBaseCommSet();
              LatticeCommSet<LatSet, TypePack>& LatCommSet1 = bdBlock1.getLatticeComm();
              auto &lat1 = bdBlock1.getLat();
              auto &block1 = bdBlock1.getLat().getBlock();
              auto &baseblock1 = block1.getBaseBlock();
              if (isOverlapped(baseblock1, box1)) {
                const std::size_t s1 =
                  AddSend(BaseCommSet1, LatCommSet1, lat, block, block1, box1, dir0);
                const std::size_t r1 =
                  AddRecv(BaseCommSet1, LatCommSet1, lat, block, block1, box1, dir1);

                const std::size_t s2 =
                  AddSend(BaseCommSet, LatCommSet, lat1, block1, block, box0, dir1);
                const std::size_t r2 =
                  AddRecv(BaseCommSet, LatCommSet, lat1, block1, block, box0, dir0);

                LatCommSet1.Recvs[r1].SendCommIndex = static_cast<int>(s2);
                LatCommSet.Recvs[r2].SendCommIndex = static_cast<int>(s1);

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
              auto &block0 = bdBlock0.getLat().getBlock();
              auto &baseblock0 = block0.getBaseBlock();
              if (isOverlapped(baseblock0, box1)) {
                const std::size_t s0 =
                  AddSend(BaseCommSet0, LatCommSet0, lat, block, block0, box0, dir1);
                const std::size_t r0 =
                  AddRecv(BaseCommSet0, LatCommSet0, lat, block, block0, box0, dir0);

                const std::size_t s1 =
                  AddSend(BaseCommSet, LatCommSet, lat0, block0, block, box1, dir0);
                const std::size_t r1 =
                  AddRecv(BaseCommSet, LatCommSet, lat0, block0, block, box1, dir1);

                LatCommSet0.Recvs[r0].SendCommIndex = static_cast<int>(s1);
                LatCommSet.Recvs[r1].SendCommIndex = static_cast<int>(s0);

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
        if (recvcomm.SendCommIndex < 0) {
          std::cerr << "[FixedPeriodicBoundaryManager]: recv without a paired send"
                    << std::endl;
          std::exit(-1);
        }
        const LatticeComm<LatSet, TypePack>& sendcomm =
          sendlatcomm.Sends[static_cast<std::size_t>(recvcomm.SendCommIndex)];
        // The copy pairs the lists by position, so a mismatch is a broken
        // pairing; without this check it reads past the shorter list.
        if (recvcomm.Comm.Cells.size() != sendcomm.Comm.Cells.size()) {
          std::cerr << "[FixedPeriodicBoundaryManager]: recv/send cell lists "
                       "differ ("
                    << recvcomm.Comm.Cells.size() << " vs "
                    << sendcomm.Comm.Cells.size() << ")" << std::endl;
          std::exit(-1);
        }
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
