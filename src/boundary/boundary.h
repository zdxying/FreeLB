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

// boundary.h

#pragma once

#include "boundary/bounce_back_boundary.h"
#include "boundary/periodic_boundary.h"


// --------------------------------------------------------------------------
// BBLike Fixed BlockBoundary
// --------------------------------------------------------------------------

template <typename BLOCKLATTICE, typename ArrayType>
class BBLikeFixedBlockBoundary : public BlockFixedBoundary<BLOCKLATTICE, ArrayType> {
 public:
  using CELL = typename BLOCKLATTICE::CellType;

  BBLikeFixedBlockBoundary(BLOCKLATTICE &lat, const ArrayType &f, std::uint8_t cellflag,
                           std::uint8_t voidflag)
      : BlockFixedBoundary<BLOCKLATTICE, ArrayType>(lat, f, cellflag, voidflag) {}

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics() {
    CELL cell(0, this->Lat);
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
    for (const auto &bdcell : this->BdCells) {
      cell.setId(bdcell.Id);
      for (unsigned int k : bdcell.outflows) {
        CELLDYNAMICS::apply(cell, k);
      }
    }
  }
};

template <typename CELLDYNAMICS, typename BLOCKLATTICEMANAGER, typename BLOCKFIELDMANAGER>
class BBLikeFixedBlockBdManager final : public AbstractBlockBoundary {
 public:
  using BLOCKLATTICE = typename BLOCKLATTICEMANAGER::BLOCKLATTICE;
  using ArrayType = typename BLOCKFIELDMANAGER::array_type;

 private:
  std::string _name;
  std::vector<BBLikeFixedBlockBoundary<BLOCKLATTICE, ArrayType>> BdBlocks;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;

  BLOCKLATTICEMANAGER &LatMan;

  BLOCKFIELDMANAGER &BlockFManager;

 public:
  BBLikeFixedBlockBdManager(std::string name, BLOCKLATTICEMANAGER &lat,
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
  void Apply(std::int64_t count) override {
    std::uint8_t MaxLevel = LatMan.getMaxLevel();
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
    for (auto &bdBlock : BdBlocks) {
      if (count %
            (static_cast<int>(std::pow(2, int(MaxLevel - bdBlock.getLat().getLevel())))) ==
          0)
        bdBlock.template ApplyCellDynamics<CELLDYNAMICS>();
    }
  }
};


// --------------------------------------------------------------------------
// BBLike Moving BlockBoundary
// --------------------------------------------------------------------------

template <typename BLOCKLATTICE, typename ArrayType>
class BBLikeMovingBlockBoundary : public BlockMovingBoundary<BLOCKLATTICE, ArrayType> {
 public:
  using CELL = typename BLOCKLATTICE::CellType;
  using LatSet = typename CELL::LatticeSet;

  BBLikeMovingBlockBoundary(BLOCKLATTICE &lat, std::vector<std::size_t> *ids,
                            ArrayType &f, std::uint8_t voidflag, std::uint8_t cellflag)
      : BlockMovingBoundary<BLOCKLATTICE, ArrayType>(lat, *ids, f, voidflag, cellflag) {}

  template <typename CELLDYNAMICS>
  void ApplyCellDynamics() {
    CELL cell(0, this->Lat);
#ifdef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num) schedule(static) firstprivate(cell)
#endif
    for (std::size_t id : this->Ids) {
      cell.setId(id);
      for (unsigned int k = 1; k < LatSet::q; ++k) {
        if (util::isFlag(this->Field[this->Lat.getNbrId(id, k)], this->voidFlag)) {
          CELLDYNAMICS::apply(cell, latset::opp<LatSet>(k));
        }
      }
    }
  }
};

template <typename CELLDYNAMICS, typename BLOCKLATTICEMANAGER, typename BLOCKFIELDMANAGER>
class BBLikeMovingBlockBdManager final : public AbstractBlockBoundary {
 public:
  using BLOCKLATTICE = typename BLOCKLATTICEMANAGER::BLOCKLATTICE;
  using ArrayType = typename BLOCKFIELDMANAGER::array_type;

 private:
  std::string _name;
  std::vector<BBLikeMovingBlockBoundary<BLOCKLATTICE, ArrayType>> BdBlocks;
  // boundary cell flag
  std::uint8_t BdCellFlag;
  // boundary flag
  std::uint8_t voidFlag;
  BLOCKLATTICEMANAGER &LatMan;
  // ids
  std::vector<std::vector<std::size_t> *> &IDss;

  BLOCKFIELDMANAGER &BlockFManager;

 public:
  BBLikeMovingBlockBdManager(std::string name, BLOCKLATTICEMANAGER &lat,
                             std::vector<std::vector<std::size_t> *> &idss,
                             BLOCKFIELDMANAGER &BlockFM, std::uint8_t voidflag,
                             std::uint8_t cellflag = std::uint8_t(0))
      : _name(name), BdCellFlag(cellflag), voidFlag(voidflag), LatMan(lat), IDss(idss),
        BlockFManager(BlockFM) {
    Init();
  }

  void Init() {
    BdBlocks.clear();
    int size = static_cast<int>(LatMan.getBlockLats().size());
    for (int i = 0; i < size; ++i) {
      BdBlocks.emplace_back(LatMan.getBlockLat(i), IDss[i],
                            BlockFManager.getBlockField(i).getField(0), voidFlag,
                            BdCellFlag);
    }
  }
  void Apply(std::int64_t count) override {
    std::uint8_t MaxLevel = LatMan.getMaxLevel();
#ifndef SingleBlock_OMP
#pragma omp parallel for num_threads(Thread_Num)
#endif
    for (auto &bdBlock : BdBlocks) {
      if (count %
            (static_cast<int>(std::pow(2, int(MaxLevel - bdBlock.getLat().getLevel())))) ==
          0)
        bdBlock.template ApplyCellDynamics<CELLDYNAMICS>();
    }
  }
};


// fixed boundary helper
template <typename BLOCKLATTICEMANAGER>
struct FixedBoundaryHelper {
  using LatSet = typename BLOCKLATTICEMANAGER::LatticeSet;

  template <typename BLOCKFIELDMANAGER>
  static void Setup(BLOCKLATTICEMANAGER &latman,
                    GenericvectorManager<FixedBdCell> &blockids,
                    BLOCKFIELDMANAGER &BlockFM, std::uint8_t BdCellFlag,
                    std::uint8_t voidFlag) {
    int i = 0;
    for (auto &Lat : latman.getBlockLats()) {
      std::vector<FixedBdCell> &BdCells = blockids.getvector(i).getvector();
      const auto &Field = BlockFM.getBlockField(i).getField();
      //
      std::size_t reserveSize;
      if constexpr (LatSet::d == 2) {
        reserveSize = (Lat.getNx() + Lat.getNy()) * 2;
      } else if constexpr (LatSet::d == 3) {
        reserveSize = (Lat.getNx() * Lat.getNy() + Lat.getNx() * Lat.getNz() +
                       Lat.getNy() * Lat.getNz()) *
                      2;
      }
      BdCells.reserve(reserveSize);
      // add inner cells
      if constexpr (LatSet::d == 2) {
        for (int iy = Lat.getOverlap(); iy < Lat.getNy() - Lat.getOverlap(); ++iy) {
          for (int ix = Lat.getOverlap(); ix < Lat.getNx() - Lat.getOverlap(); ++ix) {
            std::size_t id = ix + iy * Lat.getNx();
            if (util::isFlag(Field[id], BdCellFlag)) {
              BdCells.emplace_back(id, LatSet::q);
              FixedBdCell &fixedbdcell = BdCells.back();
              // get neighbor
              // Attention: if voidFlag is 0, DO NOT use util::isFlag
              for (int k = 1; k < LatSet::q; ++k) {
                // if (util::isFlag(Field[Lat.getNbrId(id, k)], voidFlag) &&
                //     !util::isFlag(Field[Lat.getNbrId(id, latset::opp<LatSet>(k))], voidFlag)) {
                // using util::isFlag requires voidFlag to be non-zero
                if (util::isFlag(Field[Lat.getNbrId(id, k)], voidFlag)) {
                  fixedbdcell.outflows.push_back(latset::opp<LatSet>(k));
                }
              }
            }
          }
        }
      } else if constexpr (LatSet::d == 3) {
        for (int iz = Lat.getOverlap(); iz < Lat.getNz() - Lat.getOverlap(); ++iz) {
          for (int iy = Lat.getOverlap(); iy < Lat.getNy() - Lat.getOverlap(); ++iy) {
            for (int ix = Lat.getOverlap(); ix < Lat.getNx() - Lat.getOverlap(); ++ix) {
              std::size_t id =
                ix + iy * Lat.getProjection()[1] + iz * Lat.getProjection()[2];
              if (util::isFlag(Field[id], BdCellFlag)) {
                BdCells.emplace_back(id, LatSet::q);
                FixedBdCell &fixedbdcell = BdCells.back();
                // get neighbor
                // Attention: if voidFlag is 0, DO NOT use util::isFlag
                for (int k = 1; k < LatSet::q; ++k) {
                  // if (util::isFlag(Field[Lat.getNbrId(id, k)], voidFlag) &&
                  //     !util::isFlag(Field[Lat.getNbrId(id, latset::opp<LatSet>(k))], voidFlag))
                  //     {
                  // using util::isFlag requires voidFlag to be non-zero
                  if (util::isFlag(Field[Lat.getNbrId(id, k)], voidFlag)) {
                    fixedbdcell.outflows.push_back(latset::opp<LatSet>(k));
                  }
                }
              }
            }
          }
        }
      }
      // shrink capacity to actual size
      BdCells.shrink_to_fit();
      //
      ++i;
    }
  }
};
