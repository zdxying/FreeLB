// hidden communication must be equivalent to a normal communicate
//
// ApplyCellDynamicsWithHiddenComm hides the MPI latency behind the collision of
// the interior cells.  Two things are required for it to be equivalent to
// Stream + ApplyCellDynamics + NormalCommunicate, and both were missing:
//
//   1. the shared-memory (same-rank neighbour) halo exchange -- the scheme only
//      ever did the MPI part, so a block whose neighbour lives on the same rank
//      never had its overlap refreshed.  It went unnoticed because every
//      configuration exercised so far had exactly one block per rank;
//   2. the second collision phase must cover every cell except the inner ring.
//      The communicate only writes the directions the neighbour needs; the other
//      directions of an overlap cell keep whatever is in it, and a normal
//      communicate has just collided them.  Skipping that collision made the two
//      paths differ in the overlap, which reached the fluid through the rotation
//      and the bounce-back's getPrevious().
//
// The same schedule now also drives the full (MPIComm) and all (AllMPIComm)
// exchanges, and a BlockFieldManager exchange can be split around the step with
// PostNormalCommunicate() / FinishNormalCommunicate().  All four modes are
// compared here, on one rank (hidden vs plain) and across two ranks (child
// process vs the one-rank reference).
//
// The test builds a 2x2 block grid, so every block has neighbours in more than
// one direction, and compares the two schemes cell by cell.

#include <sys/wait.h>

#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
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

static constexpr int kSteps = 20;

using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
using CELL = Cell<T, LatSet, FIELDS>;
using BulkTask =
  tmp::Key_TypePair<AABBFlag, collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
using BBTask = tmp::Key_TypePair<BouncebackFlag, collision::BounceBack<CELL>>;
using BBMVTask = tmp::Key_TypePair<BBMovingWallFlag, collision::BounceBackMovingWall<CELL>>;
using TaskCollection = tmp::TupleWrapper<BulkTask, BBTask, BBMVTask>;
using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;
using VelField = BlockFieldManager<VELOCITY<T, LatSet::d>, T, LatSet::d>;

enum class Mode { Partial, Full, All, Field };

static const char* modeName(Mode m) {
  switch (m) {
    case Mode::Partial: return "partial";
    case Mode::Full: return "full";
    case Mode::All: return "all";
    case Mode::Field: return "field";
  }
  return "?";
}

static Mode parseMode(const std::string& s) {
  if (s == "full") return Mode::Full;
  if (s == "all") return Mode::All;
  if (s == "field") return Mode::Field;
  return Mode::Partial;
}

// One lid-driven cavity on a fixed 2x2 block grid.  CreateBlocks(bx, by) is
// rank independent, LoadBalancing(np) only decides who owns which block, so the
// same block set can be run with any number of ranks.
struct Sim {
  AABB<T, 2> cavity;
  AABB<T, 2> toplid;
  BlockGeometryHelper2D<T> GeoHelper;
  BlockGeometry2D<T>* Geo;
  BlockFieldManager<FLAG, T, LatSet::d>* FlagFM;
  VelField* VelFM;
  BaseConverter<T> Conv;
  BlockLatticeManager<T, LatSet, FIELDS>* Lattice;

  Sim(int Ni, int Nj, int bx, int by)
      : cavity(Vector<T, 2>{}, Vector<T, 2>{T(Ni), T(Nj)}),
        toplid(T(1), T((Nj - 1)), T(Ni - 1), T(Nj)),
        GeoHelper(Ni, Nj, cavity, T(1), 20),
        Geo(nullptr),
        FlagFM(nullptr),
        VelFM(nullptr),
        Conv(LatSet::cs2),
        Lattice(nullptr) {
    GeoHelper.CreateBlocks(bx, by);
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
    ValuePack<T, Vector<T, LatSet::d>, T> iv(T(1), Vector<T, LatSet::d>{}, T{});
    iv = ValuePack<T, Vector<T, LatSet::d>, T>(Conv.getLatRhoInit(),
                                               Vector<T, LatSet::d>{}, T{});
    Lattice = new BlockLatticeManager<T, LatSet, FIELDS>(*Geo, iv, Conv);
    Vector<T, LatSet::d> LatU_Wall = Conv.getLatticeU(Vector<T, 2>{T(0.1), T(0)});
    Lattice->getField<VELOCITY<T, LatSet::d>>().forEach(
      toplid, *FlagFM, BBMovingWallFlag,
      [&](auto& field, std::size_t id) { field.SetField(id, LatU_Wall); });

    // Independent field manager used by the split-exchange test.  Its values
    // depend on the block id, so receiving the wrong block's halo is visible.
    VelFM = new VelField(*Geo, Vector<T, LatSet::d>{});
    for (auto& bf : VelFM->getBlockFields()) {
      const double b0 = 0.01 * bf.getBlock().getBlockId();
      auto& arr = bf.getField(0);
      for (std::size_t i = 0; i < arr.size(); ++i)
        arr[i] = Vector<T, LatSet::d>{T(b0 + 1e-3), T(-b0 - 1e-6 * double(i))};
    }
  }
  ~Sim() {
    delete Lattice;
    delete VelFM;
    delete FlagFM;
    delete Geo;
  }

  // Local update of the test field.  It reads the halos too, so over several
  // steps a wrong halo refresh propagates into the interior.
  void updateField() {
    for (auto& bf : VelFM->getBlockFields()) {
      auto& arr = bf.getField(0);
      for (std::size_t i = 0; i < arr.size(); ++i) {
        arr[i] = Vector<T, LatSet::d>{T(0.5) * arr[i][0] + T(1e-4),
                                      T(0.5) * arr[i][1] - T(1e-4)};
      }
    }
  }

  template <Mode M>
  void stepHidden() {
    if constexpr (M == Mode::Field) {
      updateField();
      VelFM->PostNormalCommunicate();
      Lattice->Stream();
      Lattice->template ApplyCellDynamicsWithHiddenComm<
        NSTask, flbcomm::CommMode::Full>(*FlagFM);
      VelFM->FinishNormalCommunicate();
    } else {
      Lattice->Stream();
      if constexpr (M == Mode::Partial) {
        // Partial is deliberately explicit: the default is Full.
        Lattice->template ApplyCellDynamicsWithHiddenComm<
          NSTask, flbcomm::CommMode::Partial>(*FlagFM);
      } else if constexpr (M == Mode::Full) {
        Lattice->template ApplyCellDynamicsWithHiddenComm<
          NSTask, flbcomm::CommMode::Full>(*FlagFM);
      } else {
        Lattice->template ApplyCellDynamicsWithHiddenComm<
          NSTask, flbcomm::CommMode::All>(*FlagFM);
      }
    }
  }

  template <Mode M>
  void stepPlain() {
    if constexpr (M == Mode::Field) {
      updateField();
      VelFM->NormalCommunicate();
      Lattice->Stream();
      Lattice->template ApplyCellDynamics<NSTask>(*FlagFM);
      Lattice->NormalFullCommunicate();
    } else {
      Lattice->Stream();
      Lattice->template ApplyCellDynamics<NSTask>(*FlagFM);
      if constexpr (M == Mode::Partial) {
        Lattice->NormalCommunicate();
      } else if constexpr (M == Mode::Full) {
        Lattice->NormalFullCommunicate();
      } else {
        Lattice->NormalAllCommunicate();
      }
    }
  }

  void step(Mode m, bool hidden) {
    switch (m) {
      case Mode::Partial: hidden ? stepHidden<Mode::Partial>() : stepPlain<Mode::Partial>(); break;
      case Mode::Full: hidden ? stepHidden<Mode::Full>() : stepPlain<Mode::Full>(); break;
      case Mode::All: hidden ? stepHidden<Mode::All>() : stepPlain<Mode::All>(); break;
      case Mode::Field: hidden ? stepHidden<Mode::Field>() : stepPlain<Mode::Field>(); break;
    }
  }

  // Full-field checksum of every local block: "blockid sum sumsq maxabs"
  void dumpBlockChecksums(const std::string& path) const {
    std::ofstream os(path, std::ios::trunc);
    os << std::setprecision(17);
    for (auto& bl : Lattice->getBlockLats()) {
      auto& popf = bl.template getField<POP<T, LatSet::q>>();
      const std::size_t n = bl.getN();
      double sum = 0, sumsq = 0, mx = 0;
      for (std::size_t id = 0; id < n; ++id)
        for (unsigned int d = 0; d < LatSet::q; ++d) {
          const double v = static_cast<double>(popf.getField(d)[id]);
          sum += v;
          sumsq += v * v;
          const double a = v < 0 ? -v : v;
          if (a > mx) mx = a;
        }
      os << bl.getBlock().getBlockId() << " " << sum << " " << sumsq << " " << mx << "\n";
    }
  }

  void dumpFieldChecksums(const std::string& path) const {
    std::ofstream os(path, std::ios::trunc);
    os << std::setprecision(17);
    for (const auto& bf : VelFM->getBlockFields()) {
      const auto& arr = bf.getField(0);
      double sum = 0, sumsq = 0, mx = 0;
      for (std::size_t id = 0; id < arr.size(); ++id) {
        const double v = static_cast<double>(arr[id][0]);
        sum += v;
        sumsq += v * v;
        const double a = v < 0 ? -v : v;
        if (a > mx) mx = a;
      }
      os << bf.getBlock().getBlockId() << " " << sum << " " << sumsq << " " << mx << "\n";
    }
  }

  void dump(Mode m, const std::string& path) const {
    if (m == Mode::Field)
      dumpFieldChecksums(path);
    else
      dumpBlockChecksums(path);
  }
};

// Spawn `count` copies of this binary in "write" mode.  system("mpirun ...")
// cannot be used from inside an MPI process; MPI_Comm_spawn is the standard
// nesting mechanism and gives the children a normal MPI_COMM_WORLD.
static void spawnWriters(int count, const std::string& mode) {
  char* cmd = (char*)"./hidden_comm.exe";
  // MPI_Comm_spawn expects a NULL-terminated argv (like execv), otherwise it
  // reads past the end -- which segfaults intermittently under mpirun.
  char* argv[] = {(char*)"write", (char*)mode.c_str(), nullptr};
  MPI_Comm inter;
  std::vector<int> errcodes(static_cast<std::size_t>(count), 0);
  const int rc = MPI_Comm_spawn(cmd, argv, count, MPI_INFO_NULL, 0, MPI_COMM_SELF,
                                &inter, errcodes.data());
  if (rc != MPI_SUCCESS) {
    IF_MPI_RANK(0) { std::cerr << "[hidden_comm] MPI_Comm_spawn failed" << std::endl; }
    CHECK(false);
    return;
  }
  MPI_Barrier(inter);
  MPI_Comm_disconnect(&inter);
}

using Rows = std::vector<std::array<double, 4>>;

static Rows load(const std::string& path) {
  Rows rows;
  std::ifstream is(path);
  std::string line;
  while (std::getline(is, line)) {
    std::istringstream ls(line);
    std::array<double, 4> r{0, 0, 0, 0};
    if (ls >> r[0] >> r[1] >> r[2] >> r[3]) rows.push_back(r);
  }
  return rows;
}

int main(int argc, char** argv) {
  mpi().init(&argc, &argv);
  const std::string modeArg = argc > 1 ? argv[1] : "main";

  // ---- child: run the hidden scheme with THIS rank's blocks ----
  if (modeArg == "write") {
    const Mode m = parseMode(argc > 2 ? argv[2] : "partial");
    Sim sim(60, 60, 2, 2);
    for (int i = 0; i < kSteps; ++i) sim.step(m, true);
    sim.dump(m, "hidden_child_" + std::string(modeName(m)) + "_" +
                  std::to_string(mpi().getRank()) + ".txt");
    MPI_Comm parent = MPI_COMM_NULL;
    MPI_Comm_get_parent(&parent);
    if (parent != MPI_COMM_NULL) {
      MPI_Barrier(parent);
      MPI_Comm_disconnect(&parent);
      // no MPI_Finalize() here: MpiManager's destructor does it, calling it
      // twice aborts with "MPI_Finalize() was called after MPI_FINALIZE"
    }
    return 0;
  }

  // ---- parent: a single process ----
  if (std::system("rm -f hidden_child_*.txt hidden_*_1rank.txt plain_*_1rank.txt") != 0)
    CHECK(false);

  const Mode modes[] = {Mode::Partial, Mode::Full, Mode::All, Mode::Field};
  for (Mode m : modes) {
    const std::string name = modeName(m);

    // 1) one rank owning all four blocks: hidden vs plain, cell by cell
    Rows hidden;
    {
      Sim hid(60, 60, 2, 2);
      Sim pln(60, 60, 2, 2);
      for (int i = 0; i < kSteps; ++i) {
        hid.step(m, true);
        pln.step(m, false);
      }
      const std::string hp = "hidden_" + name + "_1rank.txt";
      const std::string pp = "plain_" + name + "_1rank.txt";
      hid.dump(m, hp);
      pln.dump(m, pp);
      IF_MPI_RANK(0) {
        std::cout << "[hidden_comm] mode=" << name << " ranks=" << mpi().getSize()
                  << " local blocks=" << hid.Lattice->getBlockLats().size() << std::endl;
      }
      hidden = load(hp);
      const Rows plain = load(pp);
      CHECK(hidden.size() == plain.size());
      CHECK(hidden.size() == std::size_t(4));
      for (std::size_t i = 0; i < hidden.size(); ++i) {
        CHECK(hidden[i][0] == plain[i][0]);  // blockid
        CHECK(hidden[i][1] == plain[i][1]);  // sum
        CHECK(hidden[i][2] == plain[i][2]);  // sumsq
        CHECK(hidden[i][3] == plain[i][3]);  // maxabs
      }
    }

    // 2) two ranks, two blocks each: must match the one-rank reference
    spawnWriters(2, name);
    {
      Rows child;
      for (int r = 0; r < 4; ++r) {
        const Rows rows =
          load("hidden_child_" + name + "_" + std::to_string(r) + ".txt");
        if (!rows.empty()) child.insert(child.end(), rows.begin(), rows.end());
      }
      CHECK(child.size() == std::size_t(4));
      for (const auto& c : child) {
        bool matched = false;
        for (const auto& h : hidden)
          if (h[0] == c[0]) {
            CHECK(h[1] == c[1]);
            CHECK(h[2] == c[2]);
            CHECK(h[3] == c[3]);
            matched = true;
          }
        CHECK(matched);  // every child block must be one of the local ones
      }
    }
  }

  if (std::system("rm -f hidden_child_*.txt hidden_*_1rank.txt plain_*_1rank.txt") != 0)
    CHECK(false);

  return flbtest::summary("hidden_comm");
}
