// AMR communication: results must not depend on the rank count.
//
// The exchanges that only exist for an adaptively refined block structure --
// the level-filtered overloads and the coarse<->fine Aver/Intp conversions --
// had no test at all: nothing in tests/ called them.  This target pins them the
// same way tests/block/plain_comm_mpi pins the uniform exchanges: run the same
// decomposition on one rank (every neighbour local, no MPI in the reference)
// and on two ranks, and require the per-block checksums to match exactly.
//
// The geometry is deliberately not refblock's: that one produces blocks with a
// 22-row extent that exercise the bounds check, so it would only test that.
// Here every block is square-ish and the refined region is inset, which keeps
// the target on "does the exchange move the right values".
//
// Modes, each compared across ranks:
//   normal-count / full-count / all-count   level-filtered lattice exchange
//   field-count                             level-filtered field exchange
//   aver-count                              fine -> coarse population conversion
//   intp-count                              coarse -> fine population conversion
//
// intp-count used to be excluded: the MPI Intp exchange built its send and
// recv lists from two different enumerations of "the coarse cells that take
// part" (lengths 264 vs 65 coarse units for the block-0/block-2 pair), the
// posted buffer sizes disagreed and the run hung in the following barrier.
// Both sides enumerate one shared stencil set (src/geometry/intp_stencil.h)
// and the pair invariant send == recv * (4 in 2D / 8 in 3D) holds by
// construction.
//
// Built with mpic++ (see Makefile).  The parent re-launches itself with
// MPI_Comm_spawn, the same way plain_comm_mpi and checkpoint_decomp do.

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

static constexpr int kNi = 64;
static constexpr int kNj = 64;
static constexpr T kCellLen = T(1);
static constexpr int kBlockLen = 16;
static constexpr int kSteps = 6;

using BlockMap = std::map<int, std::array<double, 3>>;

static double scalarOf(double v) { return v; }
static double scalarOf(std::uint8_t v) { return static_cast<double>(v); }
template <typename U, unsigned int D>
static double scalarOf(const Vector<U, D>& v) {
  return static_cast<double>(v[0]);
}

struct Sim {
  AABB<T, 2> cavity;
  BlockGeometryHelper2D<T> GeoHelper;
  BlockGeometry2D<T>* Geo;
  BlockFieldManager<FLAG, T, 2>* FlagFM;
  BlockFieldManager<RHO<T>, T, 2>* RhoFM;
  BlockFieldManager<VELOCITY<T, 2>, T, 2>* VelFM;
  BaseConverter<T> Conv;
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, 2>, POP<T, LatSet::q>>;
  using CELL = Cell<T, LatSet, FIELDS>;
  using BulkTask = tmp::Key_TypePair<
    AABBFlag, collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  using TaskCollection = tmp::TupleWrapper<BulkTask>;
  using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;
  BlockLatticeManager<T, LatSet, FIELDS>* Lattice;

  Sim()
      : cavity(Vector<T, 2>{}, Vector<T, 2>{T(kNi * kCellLen), T(kNj * kCellLen)}),
        GeoHelper(kNi, kNj, cavity, kCellLen, kBlockLen),
        Geo(nullptr),
        FlagFM(nullptr),
        RhoFM(nullptr),
        VelFM(nullptr),
        Conv(LatSet::cs2),
        Lattice(nullptr) {
    // Refine an inset square by two levels and let CheckRefine() add the
    // level-1 transition ring, so getMaxLevel() == 2 and the level filter has
    // three distinct depths (blocks at deLevel 0, 1 and 2).
    AABB<T, 2> inner(Vector<T, 2>{T(24), T(24)}, Vector<T, 2>{T(40), T(40)});
    GeoHelper.forEachBlockCell([&](BasicBlock<T, 2>& block) {
      if (isOverlapped(block, inner)) block.refine(std::uint8_t(2));
    });
    GeoHelper.CheckRefine();
    GeoHelper.CreateBlocks();
    // Manual decomposition plus LoadBalancing: AdaptiveOptimization would
    // re-split the domain for a given rank count, so the runs would not be
    // solving the same problem.
    GeoHelper.LoadBalancing(mpi().getSize());
    Geo = new BlockGeometry2D<T>(GeoHelper);
    CHECK(Geo->getMaxLevel() == 2);

    Conv.ConvertFromRT(T(1), T(0.56), T(0.001), T(kNi), T(0.1), T(0.02));

    FlagFM = new BlockFieldManager<FLAG, T, 2>(*Geo, VoidFlag);
    FlagFM->forEach(cavity,
                    [&](FLAG& field, std::size_t id) { field.SetField(id, AABBFlag); });

    ValuePack<T, Vector<T, 2>, T> initValues(Conv.getLatRhoInit(), Vector<T, 2>{}, T{});
    Lattice = new BlockLatticeManager<T, LatSet, FIELDS>(*Geo, initValues, Conv);

    // Field managers whose values depend on the block id, so receiving the
    // wrong block's halo is visible rather than hidden by equal data.
    RhoFM = new BlockFieldManager<RHO<T>, T, 2>(*Geo, Conv.getLatRhoInit());
    VelFM = new BlockFieldManager<VELOCITY<T, 2>, T, 2>(*Geo, Vector<T, 2>{});
    for (auto& bf : RhoFM->getBlockFields()) {
      const double b0 = 1.0 + 0.01 * bf.getBlock().getBlockId();
      auto& arr = bf.getField(0);
      for (std::size_t i = 0; i < arr.size(); ++i) arr[i] = b0 + 1e-6 * double(i);
    }
    for (auto& bf : VelFM->getBlockFields()) {
      const double b0 = 0.01 * bf.getBlock().getBlockId();
      auto& arr = bf.getField(0);
      for (std::size_t i = 0; i < arr.size(); ++i)
        arr[i] = Vector<T, 2>{T(b0 + 1e-3), T(-b0 - 1e-6 * double(i))};
    }
  }

  ~Sim() {
    delete Lattice;
    delete FlagFM;
    delete RhoFM;
    delete VelFM;
    delete Geo;
  }

  // `count` drives the level filter.  The step deliberately does not stream or
  // collide: this target is about the exchange, and keeping the state fixed
  // makes the Aver/Intp conversion the only thing that moves values.
  void step(const std::string& mode, std::int64_t count) {
    if (mode == "normal-count") {
      Lattice->NormalCommunicate(count);
    } else if (mode == "full-count") {
      Lattice->NormalFullCommunicate(count);
    } else if (mode == "all-count") {
      Lattice->NormalAllCommunicate(count);
    } else if (mode == "field-count") {
      RhoFM->NormalCommunicate(count);
      VelFM->NormalCommunicate(count);
      FlagFM->NormalCommunicate(count);
    } else if (mode == "aver-count") {
      Lattice->AverageCommunicate(count);
    } else if (mode == "intp-count") {
      Lattice->InterpolateCommunicate(count);
    }
  }

  // Population checksum over the inner cells (the overlap layer is what the
  // exchange writes, so it legitimately depends on the path taken).
  void addPopChecksum(BlockMap& out) const {
    for (const auto& bl : Lattice->getBlockLats()) {
      auto& popf = bl.template getField<POP<T, LatSet::q>>();
      const int nx = bl.getNx(), ny = bl.getNy();
      const int ol = bl.getOverlap();
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
    if (mode == "field-count") {
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
    os << kv.first << " " << kv.second[0] << " " << kv.second[1] << " " << kv.second[2]
       << "\n";
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

static std::string childFile(const std::string& mode, int count, int np, int rank) {
  return "child_" + mode + "_" + std::to_string(count) + "_" + std::to_string(np) + "_" +
         std::to_string(rank) + ".txt";
}

static void spawnChildren(int nranks, const std::string& mode, std::int64_t count) {
  char* cmd = (char*)"./amr_comm.exe";
  char* arg0 = (char*)"child";
  char* arg1 = (char*)mode.c_str();
  std::string countStr = std::to_string(count);
  char* arg2 = (char*)countStr.c_str();
  char* argv[] = {arg0, arg1, arg2, nullptr};
  MPI_Comm inter;
  std::vector<int> errcodes(static_cast<std::size_t>(nranks), 0);
  const int rc = MPI_Comm_spawn(cmd, argv, nranks, MPI_INFO_NULL, 0, MPI_COMM_SELF, &inter,
                                errcodes.data());
  if (rc != MPI_SUCCESS) {
    std::cerr << "[amr_comm] MPI_Comm_spawn failed" << std::endl;
    CHECK(false);
    return;
  }
  MPI_Barrier(inter);
  MPI_Comm_disconnect(&inter);
}

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);
  const std::string mode = argc > 1 ? argv[1] : "main";

  if (mode == "child") {
    const std::string comm = argc > 2 ? argv[2] : "normal-count";
    const std::int64_t count = argc > 3 ? std::atoll(argv[3]) : 1;
    Sim sim;
    for (int i = 0; i < kSteps; ++i) sim.step(comm, count + i);
    writeBlocks(childFile(comm, count, mpi().getSize(), mpi().getRank()), sim.checksum(comm));
    MPI_Comm parent = MPI_COMM_NULL;
    MPI_Comm_get_parent(&parent);
    if (parent != MPI_COMM_NULL) {
      MPI_Barrier(parent);
      MPI_Comm_disconnect(&parent);
    }
    return 0;
  }

  // count = 1,2,3,4 drives the level filter: deLevel 0 qualifies for every
  // count, deLevel 1 only for even counts, deLevel 2 only for multiples of 4.
  // Aver only gets 1, 2 and 4 so the run stays cheap.
  const char* modes[] = {"normal-count", "full-count", "all-count", "field-count",
                         "aver-count", "intp-count"};
  const std::int64_t counts[] = {1, 2, 3, 4};
  const int nranksList[] = {2};

  if (std::system("rm -f child_*.txt") != 0) CHECK(false);

  for (const char* m : modes) {
    const bool isAver = std::string(m) == "aver-count";
    for (std::int64_t count : counts) {
      if (isAver && count == 3) continue;
      BlockMap ref;
      {
        Sim sim;
        for (int i = 0; i < kSteps; ++i) sim.step(m, count + i);
        ref = sim.checksum(m);
      }
      CHECK(!ref.empty());

      for (int np : nranksList) {
        spawnChildren(np, m, count);
        BlockMap got;
        int files = 0;
        for (int r = 0; r < np; ++r) {
          BlockMap part;
          if (!readBlocks(childFile(m, count, np, r), part)) continue;
          ++files;
          for (const auto& kv : part) {
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
          // bitwise: the exchange must not move a single value
          CHECK(it->second[0] == kv.second[0]);
          CHECK(it->second[1] == kv.second[1]);
          CHECK(it->second[2] == kv.second[2]);
          ++compared;
        }
        CHECK_EQ(compared, ref.size());
        if (mpi().getRank() == 0) {
          std::cout << "[amr_comm] " << m << " count=" << count << " np=" << np << ": "
                    << got.size() << " blocks match the single-rank reference"
                    << std::endl;
        }
      }
    }
  }

  if (flbtest::failures() == 0) {
    if (std::system("rm -f child_*.txt") != 0) CHECK(false);
  }
  return flbtest::summary("amr_comm");
}
