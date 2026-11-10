// Coalesced plain communication: results must not depend on the rank count.
//
// This exercises NormalCommunicate(), NormalFullCommunicate() and
// NormalAllCommunicate() -- for both the lattice and a BlockFieldManager --
// post one message per peer instead of one per communicator.  Coalescing only
// changes where the message boundaries are, not which cell receives which
// value, so a fixed decomposition must produce bit-identical fields no matter
// how the blocks are spread over ranks.
//
// The reference is this process running alone (np = 1): every neighbour is
// local, so the MPI path is not involved at all.  The children then run the
// same decomposition on 2 and 4 ranks; their per-block checksums must match the
// reference exactly (compared at 17 digits, i.e. round-trip exact).
//
// This also covers the link order: the Dir plan orders refs by the
// receiver-derived Order and uses the constant DirTag, while Full/All use their
// constants and order by the (receiver, sender) block pair.  A wrong order
// would silently swap two blocks' halos, which shows up here immediately.
//
// The field managers start from distinct per-block values, so a mis-ordered
// exchange cannot hide behind equal data, and the checksum covers the halo
// cells (the ones the exchange writes), not just the interior.
//
// Built with mpic++ (see Makefile).  The parent process re-launches itself with
// MPI_Comm_spawn, the same way tests/block/checkpoint_decomp does.

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "freelb.h"
#include "freelb.hh"
#include "../../test_util.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

static constexpr std::uint8_t VoidFlag = std::uint8_t(1);
static constexpr std::uint8_t AABBFlag = std::uint8_t(2);
static constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
static constexpr std::uint8_t BBMovingWallFlag = std::uint8_t(8);

static constexpr int kSteps = 100;
static constexpr int kBlocksX = 4;
static constexpr int kBlocksY = 4;

using BlockMap = std::map<int, std::array<double, 3>>;

// One scalar out of a field element: the field managers below store double,
// uint8_t and Vector<double,2>, and the checksum only needs a number that
// identifies the value.
static double scalarOf(double v) { return v; }
static double scalarOf(std::uint8_t v) { return static_cast<double>(v); }
template <typename U, unsigned int D>
static double scalarOf(const Vector<U, D>& v) {
  return static_cast<double>(v[0]);
}

struct Sim {
  AABB<T, 2> cavity;
  AABB<T, 2> toplid;
  BlockGeometryHelper2D<T> GeoHelper;
  BlockGeometry2D<T>* Geo;
  BlockFieldManager<FLAG, T, LatSet::d>* FlagFM;
  BlockFieldManager<RHO<T>, T, LatSet::d>* RhoFM;
  BlockFieldManager<VELOCITY<T, LatSet::d>, T, LatSet::d>* VelFM;
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
        RhoFM(nullptr),
        VelFM(nullptr),
        Conv(LatSet::cs2),
        Lattice(nullptr) {
    // Manual decomposition: AdaptiveOptimization(np) would re-split the domain
    // for a given process count, so the runs would not be solving the same
    // problem.  LoadBalancing(np) only decides who owns which block.
    GeoHelper.CreateBlocks(blockXNum, blockYNum);
    GeoHelper.LoadBalancing(mpi().getSize());
    Geo = new BlockGeometry2D<T>(GeoHelper);

    Conv.ConvertFromRT(T(1), T(0.56), T(0.001), T(Ni), T(0.1), T(0.02));
    FlagFM = new BlockFieldManager<FLAG, T, LatSet::d>(*Geo, VoidFlag);
    FlagFM->forEach(cavity,
                    [&](FLAG& field, std::size_t id) { field.SetField(id, AABBFlag); });
    FlagFM->template SetupBoundary<LatSet>(cavity, BouncebackFlag);
    FlagFM->forEach(toplid, [&](FLAG& field, std::size_t id) {
      if (util::isFlag(field.get(id), BouncebackFlag)) field.SetField(id, BBMovingWallFlag);
    });
    ValuePack<T, Vector<T, LatSet::d>, T> initValues(T(1), Vector<T, LatSet::d>{}, T{});
    initValues =
      ValuePack<T, Vector<T, LatSet::d>, T>(Conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
    Lattice = new BlockLatticeManager<T, LatSet, FIELDS>(*Geo, initValues, Conv);
    Vector<T, LatSet::d> LatU_Wall = Conv.getLatticeU(Vector<T, 2>{T(0.1), T(0)});
    Lattice->getField<VELOCITY<T, LatSet::d>>().forEach(
      toplid, *FlagFM, BBMovingWallFlag,
      [&](auto& field, std::size_t id) { field.SetField(id, LatU_Wall); });

    // Field managers for the field communication paths.  Their values depend on
    // the block id, so receiving the wrong block's halo is visible.
    RhoFM = new BlockFieldManager<RHO<T>, T, LatSet::d>(*Geo, Conv.getLatRhoInit());
    VelFM = new BlockFieldManager<VELOCITY<T, LatSet::d>, T, LatSet::d>(
      *Geo, Vector<T, LatSet::d>{});
    for (auto& bf : RhoFM->getBlockFields()) {
      const double b0 = 1.0 + 0.01 * bf.getBlock().getBlockId();
      auto& arr = bf.getField(0);
      for (std::size_t i = 0; i < arr.size(); ++i) arr[i] = b0 + 1e-6 * double(i);
    }
    for (auto& bf : VelFM->getBlockFields()) {
      const double b0 = 0.01 * bf.getBlock().getBlockId();
      // VELOCITY stores a Vector per cell in a single component array, so the
      // buffer element is a Vector and the Vector iSend overload is used.
      auto& arr = bf.getField(0);
      for (std::size_t i = 0; i < arr.size(); ++i)
        arr[i] = Vector<T, LatSet::d>{T(b0 + 1e-3), T(-b0 - 1e-6 * double(i))};
    }
  }

  ~Sim() {
    delete Lattice;
    delete FlagFM;
    delete RhoFM;
    delete VelFM;
    delete Geo;
  }

  void step(const std::string& mode) {
    if (mode == "fieldrho") {
      RhoFM->NormalCommunicate();
      return;
    }
    if (mode == "fieldvel") {
      VelFM->NormalCommunicate();
      return;
    }
    if (mode == "fieldflag") {
      FlagFM->NormalCommunicate();
      return;
    }
    if (mode == "fieldall") {
      RhoFM->AllNormalCommunicate();
      VelFM->AllNormalCommunicate();
      FlagFM->AllNormalCommunicate();
      return;
    }
    Lattice->Stream();
    Lattice->template ApplyCellDynamics<NSTask>(*FlagFM);
    if (mode == "partial") {
      Lattice->NormalCommunicate();
    } else if (mode == "full") {
      Lattice->NormalFullCommunicate();
    } else if (mode == "mixed") {
      // Both plans in one step: they must not share or clobber each other's
      // buffers.  The full exchange is a superset of the partial one, so the
      // field has to end up identical to the full-only run.
      Lattice->NormalCommunicate();
      Lattice->NormalFullCommunicate();
    } else {
      Lattice->NormalAllCommunicate();
    }
  }

  // Population checksum over the inner cells.  The overlap layer is filled by
  // the comm path and legitimately depends on which scheme ran, so it is not
  // compared; the inner cells are what drives the next step.
  void addPopChecksum(BlockMap& out) const {
    const int ol = 1;
    for (const auto& bl : Lattice->getBlockLats()) {
      auto& popf = bl.template getField<POP<T, LatSet::q>>();
      const int nx = bl.getNx(), ny = bl.getNy();
      double sum = 0, sumsq = 0, max = 0;
      for (int j = ol; j < ny - ol; ++j) {
        for (int i = ol; i < nx - ol; ++i) {
          const std::size_t id = static_cast<std::size_t>(j) * nx + i;
          for (unsigned int d = 0; d < LatSet::q; ++d) {
            const double v = static_cast<double>(popf.getField(d)[id]);
            sum += v;
            sumsq += v * v;
            max = std::abs(v) > max ? std::abs(v) : max;
          }
        }
      }
      auto& c = out[bl.getBlock().getBlockId()];
      c[0] += sum;
      c[1] += sumsq;
      c[2] = std::max(c[2], max);
    }
  }

  // Field checksum over EVERY cell, halo included: the exchange only writes the
  // halo, so that is exactly what has to be compared.
  template <typename FM>
  void addFieldChecksum(const FM& fm, BlockMap& out) const {
    for (const auto& bf : fm.getBlockFields()) {
      double sum = 0, sumsq = 0, max = 0;
      for (unsigned int a = 0; a < FM::field_type::array_dim; ++a) {
        const auto& arr = bf.getField(a);
        for (std::size_t i = 0; i < arr.size(); ++i) {
          const double v = scalarOf(arr[i]);
          sum += v;
          sumsq += v * v;
          max = std::abs(v) > max ? std::abs(v) : max;
        }
      }
      auto& c = out[bf.getBlock().getBlockId()];
      c[0] += sum;
      c[1] += sumsq;
      c[2] = std::max(c[2], max);
    }
  }

  BlockMap checksum(const std::string& mode) const {
    BlockMap out;
    if (mode == "fieldrho") {
      addFieldChecksum(*RhoFM, out);
    } else if (mode == "fieldvel") {
      addFieldChecksum(*VelFM, out);
    } else if (mode == "fieldflag") {
      addFieldChecksum(*FlagFM, out);
    } else if (mode == "fieldall") {
      addFieldChecksum(*RhoFM, out);
      addFieldChecksum(*VelFM, out);
      addFieldChecksum(*FlagFM, out);
    } else {
      addPopChecksum(out);
    }
    return out;
  }
};

static bool writeBlocks(const std::string& path, const BlockMap& blocks) {
  std::ofstream os(path, std::ios::trunc);
  if (!os) return false;
  os << std::setprecision(17);
  for (const auto& kv : blocks) {
    os << kv.first << " " << kv.second[0] << " " << kv.second[1] << " " << kv.second[2] << "\n";
  }
  return true;
}

static bool readBlocks(const std::string& path, BlockMap& out) {
  std::ifstream is(path);
  if (!is) return false;
  std::string line;
  while (std::getline(is, line)) {
    std::istringstream ls(line);
    int id = -1;
    std::array<double, 3> c{0, 0, 0};
    ls >> id >> c[0] >> c[1] >> c[2];
    if (id >= 0) out[id] = c;
  }
  return true;
}

static std::string childFile(const std::string& mode, int np, int rank) {
  return "child_" + mode + "_" + std::to_string(np) + "_" + std::to_string(rank) + ".txt";
}

// Spawn `count` copies of this binary in "child" mode for the given comm path.
static void spawnChildren(int count, const std::string& mode) {
  char* cmd = (char*)"./plain_comm_mpi.exe";
  char* arg0 = (char*)"child";
  char* arg1 = (char*)mode.c_str();
  char* argv[] = {arg0, arg1, nullptr};
  MPI_Comm inter;
  std::vector<int> errcodes(static_cast<std::size_t>(count), 0);
  const int rc = MPI_Comm_spawn(cmd, argv, count, MPI_INFO_NULL, 0, MPI_COMM_SELF, &inter,
                                errcodes.data());
  if (rc != MPI_SUCCESS) {
    std::cerr << "[plain_comm_mpi] MPI_Comm_spawn failed" << std::endl;
    CHECK(false);
    return;
  }
  MPI_Barrier(inter);
  MPI_Comm_disconnect(&inter);
}

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);
  const std::string mode = argc > 1 ? argv[1] : "main";

  // ---- child: run the requested comm path on this rank, dump the blocks ----
  if (mode == "child") {
    const std::string comm = argc > 2 ? argv[2] : "partial";
    Sim sim(120, 120, kBlocksX, kBlocksY);
    for (int i = 0; i < kSteps; ++i) sim.step(comm);
    writeBlocks(childFile(comm, mpi().getSize(), mpi().getRank()), sim.checksum(comm));
    // Spawned child: leave the MPI world.  MpiManager's destructor calls
    // MPI_Finalize(), so do not call it here as well.
    MPI_Comm parent = MPI_COMM_NULL;
    MPI_Comm_get_parent(&parent);
    if (parent != MPI_COMM_NULL) {
      MPI_Barrier(parent);
      MPI_Comm_disconnect(&parent);
    }
    return 0;
  }

  // ---- parent: one rank, all blocks local, no MPI in the reference --------
  const char* modes[] = {"partial", "full", "all",    "mixed",
                         "fieldrho", "fieldvel", "fieldflag", "fieldall"};
  const int counts[] = {2, 4};

  if (std::system("rm -f child_*.txt ref_*.txt") != 0) CHECK(false);

  for (const char* m : modes) {
    BlockMap ref;
    {
      Sim sim(120, 120, kBlocksX, kBlocksY);
      for (int i = 0; i < kSteps; ++i) sim.step(m);
      ref = sim.checksum(m);
    }
    CHECK(!ref.empty());

    for (int np : counts) {
      spawnChildren(np, m);
      BlockMap got;
      int files = 0;
      for (int r = 0; r < np; ++r) {
        BlockMap part;
        if (!readBlocks(childFile(m, np, r), part)) continue;
        ++files;
        for (const auto& kv : part) {
          // every block must be owned by exactly one rank
          CHECK(got.find(kv.first) == got.end());
          got[kv.first] = kv.second;
        }
      }
      CHECK_EQ(files, np);
      CHECK_EQ(got.size(), ref.size());
      std::size_t compared = 0;
      for (const auto& kv : ref) {
        auto it = got.find(kv.first);
        if (it == got.end()) {
          CHECK(false);
          continue;
        }
        // bitwise: coalescing must not move a single value
        CHECK(it->second[0] == kv.second[0]);
        CHECK(it->second[1] == kv.second[1]);
        CHECK(it->second[2] == kv.second[2]);
        ++compared;
      }
      CHECK_EQ(compared, ref.size());
      if (mpi().getRank() == 0) {
        std::cout << "[plain_comm_mpi] " << m << " np=" << np << ": " << got.size()
                  << " blocks match the single-rank reference" << std::endl;
      }
    }
  }

  // keep the dumps around when something failed, they are the first thing to
  // look at; a passing run leaves nothing behind
  if (flbtest::failures() == 0) {
    if (std::system("rm -f child_*.txt ref_*.txt") != 0) CHECK(false);
  }
  return flbtest::summary("plain_comm_mpi");
}
