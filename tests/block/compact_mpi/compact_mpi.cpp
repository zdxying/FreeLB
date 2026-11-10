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

// compact_mpi.cpp
// Two-rank compact-path smoke test: InitVoxelMap prunes the communicators, the
// lattice plan is stamped with the resulting geometry generation, and the
// compact iteration visits exactly the active cells on each rank.
//
// A shared link packs Scell, Rcell, so the send slot belongs to the send block.
// The sparse pattern makes the two blocks' masks disagree at the same local
// index; a receiver-side test would then drop live pairs.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "freelb.h"
#include "freelb.hh"
#include "../../test_util.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

int Thread_Num = 1;

static constexpr std::uint8_t VoidFlag = std::uint8_t(1);
static constexpr std::uint8_t ActiveFlag = std::uint8_t(2);

struct CountDynamics {
  static inline std::size_t count = 0;

  template <typename CELL>
  static void Execute(std::uint8_t, CELL&) {
    ++count;
  }
};

static std::size_t globalSum(std::size_t local) {
  unsigned long long value = local;
  unsigned long long sum = 0;
  MPI_Allreduce(&value, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  return static_cast<std::size_t>(sum);
}

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT;

  constexpr int Ni = 32;
  constexpr int Nj = 32;
  const AABB<T, 2> cavity(Vector<T, 2>{}, Vector<T, 2>{T(Ni), T(Nj)});

  BlockGeometryHelper2D<T> helper(Ni, Nj, cavity, T(1), 8, 2, 1);
  helper.CreateBlocks(2, 2);
  helper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> geo(helper);

  BlockFieldManager<FLAG, T, 2> FlagFM(geo, VoidFlag);
  FlagFM.forEach([&](auto& field, std::size_t id) {
    field.SetField(id, ActiveFlag);
  });

  // The masks are defined in global coordinates, so both ends of a link agree
  // on whether the physical cell is active.
  auto fillMask = [&](const auto& isActive) {
    for (std::size_t iblock = 0; iblock < geo.getBlocks().size(); ++iblock) {
      Block2D<T>& block = geo.getBlock(static_cast<int>(iblock));
      auto& field = FlagFM.getBlockField(static_cast<int>(iblock)).getField(0);
      for (std::size_t id = 0; id < field.size(); ++id) {
        const Vector<T, 2> vox = block.getVoxel(block.getLoc(id));
        field[id] = isActive(vox) ? ActiveFlag : VoidFlag;
      }
    }
  };

  // Count what pruning must remove, using the sender's own field for the send
  // slot.  shiftedLive counts live pairs whose receiver-local cell at the send
  // index is void: exactly the pairs a receiver-side test would mis-prune.
  auto scanPruning = [&](std::size_t& sharedPairs, std::size_t& mpiCells,
                         std::size_t& shiftedLive) {
    sharedPairs = 0;
    mpiCells = 0;
    shiftedLive = 0;
    for (std::size_t iblock = 0; iblock < geo.getBlocks().size(); ++iblock) {
      const auto& recvField =
        FlagFM.getBlockField(static_cast<int>(iblock)).getField(0);
      const Communicator& COMM =
        geo.getBlock(static_cast<int>(iblock)).getCommunicator();
      auto scanShared = [&](const std::vector<SharedComm>& comms) {
        for (const SharedComm& comm : comms) {
          const auto& sendField =
            FlagFM
              .getBlockField(
                static_cast<int>(geo.findBlockIndex(comm.SendBlockId)))
              .getField(0);
          for (std::size_t i = 0; i + 1 < comm.SendRecvCells.size(); i += 2) {
            const std::size_t sendId = comm.SendRecvCells[i];
            const std::size_t recvId = comm.SendRecvCells[i + 1];
            const bool sendLive = sendField[sendId] != VoidFlag;
            const bool recvLive = recvField[recvId] != VoidFlag;
            if (!sendLive || !recvLive) {
              ++sharedPairs;
            } else if (sendId < recvField.size() &&
                       recvField[sendId] == VoidFlag) {
              ++shiftedLive;
            }
          }
        }
      };
      scanShared(COMM.Comm.Comms);
      scanShared(COMM.AllComm.Comms);
#ifdef MPI_ENABLED
      auto scanMpi = [&](const std::vector<DistributedComm>& comms) {
        for (const DistributedComm& comm : comms) {
          for (std::size_t id : comm.Cells) {
            if (recvField[id] == VoidFlag) ++mpiCells;
          }
        }
      };
      scanMpi(COMM.MPIComm.Recvs);
      scanMpi(COMM.MPIComm.Sends);
      scanMpi(COMM.AllMPIComm.Recvs);
      scanMpi(COMM.AllMPIComm.Sends);
      scanMpi(COMM.DirRecvs);
      scanMpi(COMM.DirSends);
#endif
    }
  };

  // Every pair that survived pruning must still have both ends live on their
  // own side.
  auto countDeadSurvivors = [&]() {
    std::size_t dead = 0;
    for (std::size_t iblock = 0; iblock < geo.getBlocks().size(); ++iblock) {
      const auto& recvField =
        FlagFM.getBlockField(static_cast<int>(iblock)).getField(0);
      const Communicator& COMM =
        geo.getBlock(static_cast<int>(iblock)).getCommunicator();
      auto scanShared = [&](const std::vector<SharedComm>& comms) {
        for (const SharedComm& comm : comms) {
          const auto& sendField =
            FlagFM
              .getBlockField(
                static_cast<int>(geo.findBlockIndex(comm.SendBlockId)))
              .getField(0);
          for (std::size_t i = 0; i + 1 < comm.SendRecvCells.size(); i += 2) {
            if (sendField[comm.SendRecvCells[i]] == VoidFlag ||
                recvField[comm.SendRecvCells[i + 1]] == VoidFlag) {
              ++dead;
            }
          }
        }
      };
      scanShared(COMM.Comm.Comms);
      scanShared(COMM.AllComm.Comms);
    }
    return globalSum(dead);
  };

  // --- sparse pattern: the send slot must use the send block's mask --------
  fillMask([](const Vector<T, 2>& vox) {
    const int gx = static_cast<int>(std::floor(static_cast<double>(vox[0])));
    const int gy = static_cast<int>(std::floor(static_cast<double>(vox[1])));
    return (gx + 2 * gy) % 5 == 0;
  });

  std::size_t expectedShared = 0;
  std::size_t expectedMpi = 0;
  std::size_t shiftedLive = 0;
  scanPruning(expectedShared, expectedMpi, shiftedLive);
  const std::size_t globalExpectedShared = globalSum(expectedShared);
  const std::size_t globalExpectedMpi = globalSum(expectedMpi);
  const std::size_t globalShiftedLive = globalSum(shiftedLive);
  // The pattern must actually exercise both kinds of pruning and leave live
  // pairs that only the sender-side test keeps.
  CHECK(globalExpectedShared > 0);
  CHECK(globalExpectedMpi > 0);
  CHECK(globalShiftedLive > 0);

  const std::uint64_t sparseGeneration = geo.getGeometryGeneration();
  geo.InitVoxelMap(FlagFM);
  CHECK(geo.getGeometryGeneration() > sparseGeneration);
  CHECK_EQ(geo.getLastPrunedSharedPairs(), globalExpectedShared);
  CHECK_EQ(geo.getLastPrunedMpiCells(), globalExpectedMpi);
  CHECK_EQ(countDeadSurvivors(), std::size_t(0));
  if (mpi().getRank() == 0) {
    std::cout << "[block/compact_mpi] sparse pruning: shared pairs "
              << geo.getLastPrunedSharedPairs() << ", MPI cells "
              << geo.getLastPrunedMpiCells() << ", shifted live pairs "
              << globalShiftedLive << std::endl;
  }

  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  BaseConverter<T> conv(LatSet::cs2);
  ValuePack initValues(conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(geo, initValues, conv);

  auto checkActiveCount = [&](const char* phase) {
    std::size_t expected = 0;
    for (std::size_t iblock = 0; iblock < geo.getBlocks().size(); ++iblock) {
      const auto& field =
        FlagFM.getBlockField(static_cast<int>(iblock)).getField(0);
      for (std::size_t id = 0; id < field.size(); ++id) {
        if (field[id] != VoidFlag) ++expected;
      }
    }
    CountDynamics::count = 0;
    NSLattice.ApplyCellDynamics<CountDynamics>(FlagFM);
    CHECK_EQ(CountDynamics::count, expected);
    const std::size_t globalCount = globalSum(CountDynamics::count);
    CHECK_EQ(globalCount, globalSum(expected));
    if (mpi().getRank() == 0) {
      std::cout << "[block/compact_mpi] " << phase << ": active cells "
                << globalCount << std::endl;
    }
  };

  // Build the full MPI plan after the map exists; the generation stamp must
  // accept it.  The plan itself is exercised by the existing comm suite.
  NSLattice.NormalFullCommunicate();
  checkActiveCount("sparse");

  // --- slab: keep the original cross-face coverage check -------------------
  fillMask([](const Vector<T, 2>& vox) {
    return vox[0] >= T(15) && vox[0] <= T(17);
  });
  const std::uint64_t slabGeneration = geo.getGeometryGeneration();
  geo.InitVoxelMap(FlagFM);
  CHECK(geo.getGeometryGeneration() > slabGeneration);
  checkActiveCount("slab");

  if (mpi().getRank() == 0) {
    std::cout << "[block/compact_mpi] ranks " << mpi().getSize()
              << ", geometry generation " << geo.getGeometryGeneration()
              << ", pruned shared pairs " << geo.getLastPrunedSharedPairs()
              << ", pruned MPI cells " << geo.getLastPrunedMpiCells()
              << std::endl;
  }

  return flbtest::summary("block/compact_mpi");
}
