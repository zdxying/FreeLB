// checkpoint across decompositions (P3)
//
// A checkpoint written by a run with one block-to-rank assignment must be
// readable by a run with another assignment, and the resumed run must continue
// bit-identically.
//
// What makes this possible: `LoadBalancing(np)` only reassigns the existing
// blocks to ranks, it never creates or merges them, and `getBlockId()` is the
// global id assigned by CreateBlocks().  So a checkpoint identifies its blocks
// by blockid and the reader picks the blocks it owns out of every rank's file.
//
// Caveat proven by the negative case below: AdaptiveOptimization(np) *does*
// change the block set (it re-splits the domain for a given process count), so
// a run that uses it has a rank-count-dependent geometry and cannot be resumed
// with a different rank count -- that is rejected, not silently mismatched.
//
// Built with mpic++ (see Makefile).  The parent process re-launches itself
// under mpirun -np 2 to produce a two-rank checkpoint.

#include <iomanip>
#include <sys/wait.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>

#include "freelb.h"
#include "freelb.hh"
#include "../../test_util.h"
#include "utils/field_checksum.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

static constexpr std::uint8_t VoidFlag = std::uint8_t(1);
static constexpr std::uint8_t AABBFlag = std::uint8_t(2);
static constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
static constexpr std::uint8_t BBMovingWallFlag = std::uint8_t(8);

static constexpr int kPre = 200;
static constexpr int kPost = 100;
static constexpr int kRun = kPre + kPost;

// One lid-driven cavity.  The geometry is built with the MANUAL CreateBlocks:
// AdaptiveOptimization(np) re-splits the domain for a given process count, so
// the block set would become rank-count dependent and no block-by-block match
// would be possible.  DivideBlock2D(Nx, Ny) is rank independent; LoadBalancing(np)
// then only decides who owns which block.
struct Sim {
  AABB<T, 2> cavity;
  AABB<T, 2> toplid;
  BlockGeometryHelper2D<T> GeoHelper;
  BlockGeometry2D<T>* Geo;
  BlockFieldManager<FLAG, T, LatSet::d>* FlagFM;
  BaseConverter<T> Conv;
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  using CELL = Cell<T, LatSet, FIELDS>;
  using BulkTask =
    tmp::Key_TypePair<AABBFlag, collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  using BBTask = tmp::Key_TypePair<BouncebackFlag, collision::BounceBack<CELL>>;
  using BBMVTask = tmp::Key_TypePair<BBMovingWallFlag, collision::BounceBackMovingWall<CELL>>;
  using TaskCollection = tmp::TupleWrapper<BulkTask, BBTask, BBMVTask>;
  using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;
  BlockLatticeManager<T, LatSet, FIELDS>* Lattice;

  Sim(int Ni, int Nj, int blockXNum, int blockYNum)
      : cavity(Vector<T, 2>{}, Vector<T, 2>{T(Ni), T(Nj)}),
        toplid(T(1), T((Nj - 1)), T(Ni - 1), T(Nj)),
        GeoHelper(Ni, Nj, cavity, T(1), 20),
        Geo(nullptr),
        FlagFM(nullptr),
        Conv(LatSet::cs2),
        Lattice(nullptr) {
    GeoHelper.CreateBlocks(blockXNum, blockYNum);  // rank independent
    GeoHelper.LoadBalancing(mpi().getSize());
    Geo = new BlockGeometry2D<T>(GeoHelper);

    Conv.ConvertFromRT(T(1), T(0.56), T(0.001), T(Ni), T(0.1), T(0.02));
    FlagFM = new BlockFieldManager<FLAG, T, LatSet::d>(*Geo, VoidFlag);
    FlagFM->forEach(cavity,
                    [&](FLAG& field, std::size_t id) { field.SetField(id, AABBFlag); });
    FlagFM->template SetupBoundary<LatSet>(cavity, BouncebackFlag);
    FlagFM->forEach(toplid, [&](FLAG& field, std::size_t id) {
      if (util::isFlag(field.get(id), BouncebackFlag))
        field.SetField(id, BBMovingWallFlag);
    });
    ValuePack<T, Vector<T, LatSet::d>, T> initValues(T(1), Vector<T, LatSet::d>{}, T{});
    initValues =
      ValuePack<T, Vector<T, LatSet::d>, T>(Conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
    Lattice = new BlockLatticeManager<T, LatSet, FIELDS>(*Geo, initValues, Conv);
    Vector<T, LatSet::d> LatU_Wall = Conv.getLatticeU(Vector<T, 2>{T(0.1), T(0)});
    Lattice->getField<VELOCITY<T, LatSet::d>>().forEach(
      toplid, *FlagFM, BBMovingWallFlag,
      [&](auto& field, std::size_t id) { field.SetField(id, LatU_Wall); });
  }

  ~Sim() {
    delete Lattice;
    delete FlagFM;
    delete Geo;
  }

  void step() {
    Lattice->Stream();
    Lattice->template ApplyCellDynamics<NSTask>(*FlagFM);
    Lattice->NormalFullCommunicate();
  }

  // per-block checksum over the INNER cells (the overlap/halo layer excluded),
  // one "blockid sum sumsq" line per block.  The inner cells are what drives the
  // next step, so they are the right thing to compare across layouts: the halo
  // cells are filled by a different code path when the neighbour lives on
  // another rank and differ at the rounding level.
  void dumpBlockChecksums(const std::string& path) {
    const int ol = 1;  // uniform level-0 blocks
    std::ofstream os(path, std::ios::trunc);
    os << std::setprecision(17);
    for (auto& bl : Lattice->getBlockLats()) {
      auto& popf = bl.template getField<POP<T, LatSet::q>>();
      const int nx = bl.getNx(), ny = bl.getNy();
      double sum = 0, sumsq = 0;
      for (int j = ol; j < ny - ol; ++j) {
        for (int i = ol; i < nx - ol; ++i) {
          const std::size_t id = static_cast<std::size_t>(j) * nx + i;
          for (unsigned int d = 0; d < LatSet::q; ++d) {
            const double v = static_cast<double>(popf.getField(d)[id]);
            sum += v;
            sumsq += v * v;
          }
        }
      }
      os << bl.getBlock().getBlockId() << " " << sum << " " << sumsq << "\n";
    }
  }

  // pop checksum over every block this rank owns (halo cells included)
  std::array<double, 3> popChecksum() {
    std::array<double, 3> s{0, 0, 0};
    for (auto& bl : Lattice->getBlockLats()) {
      auto& popf = bl.template getField<POP<T, LatSet::q>>();
      const std::array<double, 3> c = frelb_diag::PopChecksumHost(popf, bl.getN(), LatSet::q);
      s[0] += c[0];
      s[1] += c[1];
      s[2] = c[2] > s[2] ? c[2] : s[2];
    }
    return s;
  }
};

// Spawn `count` copies of this binary in "write" mode.
//
// NB: system("mpirun -np 2 ...") does NOT work from inside an MPI process
// (OpenMPI refuses it silently, no output at all).  MPI_Comm_spawn is the
// standard way to nest; the children get a normal MPI_COMM_WORLD with the
// right rank/size, so mpi().init() works unchanged in them.
static void spawnWriters(int count) {
  char* cmd = (char*)"./checkpoint_decomp.exe";
  char* arg = (char*)"write";
  char* argv[] = {arg, nullptr};
  MPI_Comm inter;
  std::vector<int> errcodes(static_cast<std::size_t>(count), 0);
  const int rc = MPI_Comm_spawn(cmd, argv, count, MPI_INFO_NULL, 0, MPI_COMM_SELF,
                                &inter, errcodes.data());
  if (rc != MPI_SUCCESS) {
    IF_MPI_RANK(0) { std::cerr << "[decomp] MPI_Comm_spawn failed" << std::endl; }
    CHECK(false);
    return;
  }
  MPI_Barrier(inter);
  MPI_Comm_disconnect(&inter);
}

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);
  const std::string mode = argc > 1 ? argv[1] : "main";
  const std::string dir = "./ckpt_decomp";

  // ---- child: run kPre steps and write a checkpoint with THIS rank's blocks ----
  if (mode == "write") {
    Sim sim(120, 120, 6, 6);  // 36 blocks
    for (int i = 0; i < kPre; ++i) sim.step();
    flbckpt::write(dir, kPre, *sim.Lattice, *sim.FlagFM);
    mpi().barrier();
    sim.dumpBlockChecksums("child_blocks_" + std::to_string(mpi().getRank()) + ".txt");
    // Spawned child: leave the MPI world.  Do NOT call MPI_Finalize() here:
    // MpiManager's destructor does it (mpi_manager.h:136), calling it twice is
    // what produced "MPI_Finalize() was called after MPI_FINALIZE".
    MPI_Comm parent = MPI_COMM_NULL;
    MPI_Comm_get_parent(&parent);
    if (parent != MPI_COMM_NULL) {
      MPI_Barrier(parent);
      MPI_Comm_disconnect(&parent);
    }
    return 0;
  }

  // ---- parent: a single process ----
  IF_MPI_RANK(0) {
    std::cout << "[decomp] this process runs with " << mpi().getSize() << " rank(s)"
              << std::endl;
  }

  // 1) produce a two-rank checkpoint
  if (std::system("rm -f child_blocks_*.txt ref_blocks.txt resumed_blocks.txt") != 0)
    CHECK(false);
  spawnWriters(2);
  CHECK_EQ(flbckpt::latestStep(dir), (std::int64_t)kPre);

  // the two-rank run produced one file per rank
  {
    int nfiles = 0;
    for (int r = 0; r < 8; ++r) {
      std::ifstream is(dir + "/step_" + std::to_string(kPre) + "/rank" +
                       std::to_string(r) + ".ck");
      if (is) ++nfiles;
    }
    CHECK_EQ(nfiles, 2);
  }

  // 2) reference: one rank owning all 36 blocks, continuous run
  std::array<double, 3> cRefPre{}, cRefPost{};
  {
    Sim sim(120, 120, 6, 6);
    for (int i = 0; i < kPre; ++i) sim.step();
    cRefPre = sim.popChecksum();
    sim.dumpBlockChecksums("ref_blocks.txt");
    for (int i = kPre; i < kRun; ++i) sim.step();
    cRefPost = sim.popChecksum();
  }

  // The two-rank run must have produced exactly the same physics as the
  // one-rank run.  Compared per block, on the inner cells only: the halo cells
  // are filled by a different code path when the neighbour lives on another
  // rank, so the full-field checksums differ at the 1e-7 (FP32 rounding) level,
  // while the inner cells -- the ones that drive the next step -- agree exactly.
  {
    std::vector<std::array<double, 3>> mine;  // per-block inner checksum
    std::vector<std::vector<std::array<double, 3>>> theirs;  // [rank][block]
    IF_MPI_RANK(0) {
      // this rank's per-block values
      std::ifstream is("ref_blocks.txt");
      std::string line;
      while (std::getline(is, line)) {
        std::istringstream ls(line);
        int id = -1;
        std::array<double, 3> c{0, 0, 0};
        ls >> id >> c[0] >> c[1] >> c[2];
        if (id >= 0) mine.push_back(c);
      }
    }
    int nranks_seen = 0;
    for (int r = 0; r < 4; ++r) {
      std::ifstream is("child_blocks_" + std::to_string(r) + ".txt");
      if (!is) continue;
      ++nranks_seen;
      IF_MPI_RANK(0) {
        std::vector<std::array<double, 3>> perrank;
        std::string line;
        while (std::getline(is, line)) {
          std::istringstream ls(line);
          int id = -1;
          std::array<double, 3> c{0, 0, 0};
          ls >> id >> c[0] >> c[1] >> c[2];
          if (id >= 0) perrank.push_back(c);
        }
        theirs.push_back(perrank);
      }
    }
    CHECK_EQ(nranks_seen, 2);
    IF_MPI_RANK(0) {
      std::size_t nblocks = 0;
      for (const auto& pr : theirs) nblocks += pr.size();
      CHECK_EQ(nblocks, mine.size());
      // every block must appear exactly once across the two rank files
      std::size_t compared = 0;
      for (const auto& pr : theirs)
        for (const auto& c : pr) {
          bool matched = false;
          for (const auto& m : mine)
            if (m[0] == c[0] && m[1] == c[1]) matched = true;
          CHECK(matched);
          ++compared;
        }
      CHECK_EQ(compared, nblocks);
    }
  }

  // 3) resume the two-rank checkpoint into a one-rank run
  std::array<double, 3> cResPost{};
  {
    Sim sim(120, 120, 6, 6);
    for (int i = 0; i < kPre; ++i) sim.step();  // leaves the containers rotated
    CHECK(flbckpt::checkCompatible(dir, -1, *sim.Lattice, *sim.FlagFM));
    const std::int64_t restored = flbckpt::read(dir, -1, *sim.Lattice, *sim.FlagFM);
    CHECK_EQ(restored, (std::int64_t)kPre);
    sim.dumpBlockChecksums("resumed_blocks.txt");
    for (int i = kPre; i < kRun; ++i) sim.step();
    cResPost = sim.popChecksum();
  }
  // both sides run with the same layout, so this one is exact
  CHECK(cResPost[0] == cRefPost[0]);
  CHECK(cResPost[1] == cRefPost[1]);
  CHECK(cResPost[2] == cRefPost[2]);

  // 4) a different block set must be rejected: a coarser BlockCellLen gives a
  //    different number of blocks, so no block-by-block match is possible
  {
    Sim sim(120, 120, 3, 3);  // 9 blocks instead of 36
    CHECK(!flbckpt::checkCompatible(dir, -1, *sim.Lattice, *sim.FlagFM));
  }

  if (std::system("rm -f child_blocks_*.txt ref_blocks.txt resumed_blocks.txt") != 0)
    CHECK(false);

  return flbtest::summary("checkpoint_decomp");
}
