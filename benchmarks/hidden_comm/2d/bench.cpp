// hidden communication benchmark, 2D
//
// Measures the throughput of ApplyCellDynamicsWithHiddenComm and prints a
// phase breakdown.  Build with _HIDDEN_PROFILE (the Makefile does) to get the
// breakdown; the timers are compiled out of normal builds.
//
//   ./hidden_comm_bench2d.exe [Ni] [blocksPerSide] [steps] [mode] [scheme]
//   mpirun -np 2 ./hidden_comm_bench2d.exe 120 10 3000 full hidden
//
// mode:   partial | full | all | field
// scheme: hidden  | plain
//
// The block grid is fixed (CreateBlocks is rank independent); LoadBalancing
// decides which rank owns which block, so the same command works for any -np.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;
using LatSet = D2Q9<T>;

static constexpr std::uint8_t VoidFlag = std::uint8_t(1);
static constexpr std::uint8_t AABBFlag = std::uint8_t(2);
static constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
static constexpr std::uint8_t BBMovingWallFlag = std::uint8_t(8);

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

static Mode parseMode(const char* s) {
  const std::string m = s ? s : "partial";
  if (m == "full") return Mode::Full;
  if (m == "all") return Mode::All;
  if (m == "field") return Mode::Field;
  return Mode::Partial;
}

static const char* modeName(Mode m) {
  switch (m) {
    case Mode::Partial: return "partial";
    case Mode::Full: return "full";
    case Mode::All: return "all";
    case Mode::Field: return "field";
  }
  return "?";
}

// Lid-driven cavity on a fixed block grid.  Only the hidden-communication step
// is timed; the setup is identical to tests/feature/hidden_comm.
struct Sim {
  AABB<T, 2> cavity;
  AABB<T, 2> toplid;
  BlockGeometryHelper2D<T> GeoHelper;
  BlockGeometry2D<T>* Geo;
  BlockFieldManager<FLAG, T, LatSet::d>* FlagFM;
  VelField* VelFM;
  BaseConverter<T> Conv;
  BlockLatticeManager<T, LatSet, FIELDS>* Lattice;

  Sim(int Ni, int Nj, int bx, int by, bool withField)
      : cavity(Vector<T, 2>{}, Vector<T, 2>{T(Ni), T(Nj)}),
        toplid(T(1), T(Nj - 1), T(Ni - 1), T(Nj)),
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

    if (withField) {
      VelFM = new VelField(*Geo, Vector<T, LatSet::d>{});
      for (auto& bf : VelFM->getBlockFields()) {
        const double b0 = 0.01 * bf.getBlock().getBlockId();
        auto& arr = bf.getField(0);
        for (std::size_t i = 0; i < arr.size(); ++i)
          arr[i] = Vector<T, LatSet::d>{T(b0 + 1e-3), T(-b0 - 1e-6 * double(i))};
      }
    }
  }
  ~Sim() {
    delete Lattice;
    delete VelFM;
    delete FlagFM;
    delete Geo;
  }

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
      Lattice->template ApplyCellDynamicsWithHiddenComm<NSTask>(*FlagFM);
      VelFM->FinishNormalCommunicate();
    } else {
      Lattice->Stream();
      if constexpr (M == Mode::Partial) {
        Lattice->template ApplyCellDynamicsWithHiddenComm<NSTask>(*FlagFM);
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
      Lattice->NormalCommunicate();
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

  std::size_t localCells() const {
    std::size_t n = 0;
    for (auto& bl : Lattice->getBlockLats()) n += bl.getN();
    return n;
  }
};

int main(int argc, char** argv) {
  mpi().init(&argc, &argv);
  const int Ni = argc > 1 ? std::atoi(argv[1]) : 120;
  const int bx = argc > 2 ? std::atoi(argv[2]) : 10;
  const int steps = argc > 3 ? std::atoi(argv[3]) : 3000;
  const Mode mode = parseMode(argc > 4 ? argv[4] : nullptr);
  const bool hidden = argc > 5 ? std::string(argv[5]) != "plain" : true;

  Sim sim(Ni, Ni, bx, bx, mode == Mode::Field);

  // Message layout: the communicators the geometry produces, and how they are
  // coalesced into one message per target rank.  Most of the raw communicators
  // carry a single cell (the corner ones).
  std::map<std::size_t, std::size_t> hist;
  std::map<int, std::size_t> perRank;
  std::size_t entries = 0, cells = 0, blocks = 0;
  for (auto& bl : sim.Lattice->getBlockLats()) {
    ++blocks;
    for (auto& comm : bl.getBlock().getCommunicator().DirSends) {
      ++entries;
      cells += comm.Cells.size();
      perRank[comm.TargetRank] += comm.Cells.size();
      ++hist[comm.Cells.size()];
    }
  }
  if (mpi().getRank() == 0) {
    std::printf("[config] Ni=%d blocks=%dx%d ranks=%d steps=%d mode=%s scheme=%s\n", Ni, bx, bx,
                mpi().getSize(), steps, modeName(mode), hidden ? "hidden" : "plain");
    std::printf("[layout] local blocks=%zu  communicators/step=%zu  cells/comm=%.1f\n", blocks,
                entries, entries ? double(cells) / entries : 0.0);
    std::printf("[layout] comm size histogram (cells x count):");
    for (auto& kv : hist) std::printf(" %zux%zu", kv.first, kv.second);
    std::printf("\n");
    std::printf("[layout] coalesced messages/step=%zu  cells/message=%.1f\n", perRank.size(),
                perRank.empty() ? 0.0 : double(cells) / perRank.size());
  }

  for (int i = 0; i < 100; ++i) sim.step(mode, hidden);  // warmup
  mpi().barrier();

  double best = 1e30;
  for (int rep = 0; rep < 3; ++rep) {
    mpi().barrier();
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < steps; ++i) sim.step(mode, hidden);
    const auto t1 = std::chrono::steady_clock::now();
    const double s = std::chrono::duration<double>(t1 - t0).count();
    if (s < best) best = s;
    mpi().barrier();
  }

  if (mpi().getRank() == 0) {
    const double mlups = double(sim.localCells()) * steps / best / 1e6;
    std::printf("[result] best=%.4f s  MLUPS=%.2f\n", best, mlups);
#ifdef _HIDDEN_PROFILE
    const HiddenCommProfile& p = g_hiddenCommProfile;
    const double tot = p.collideFirst + p.send + p.recv + p.collideLast + p.waitall + p.shared + p.set;
    auto pct = [tot](double v) { return tot > 0 ? 100.0 * v / tot : 0.0; };
    std::printf("[profile] steps=%ld  requests/step=%.1f\n", p.calls,
                p.calls ? double(p.sendRequests) / p.calls : 0.0);
    std::printf("[profile] collideFirst=%5.1f%%  send=%5.1f%%  recv=%5.1f%%  collideLast=%5.1f%%\n",
                pct(p.collideFirst), pct(p.send), pct(p.recv), pct(p.collideLast));
    std::printf("[profile] waitall=%5.1f%%  shared=%5.1f%%  set=%5.1f%%\n", pct(p.waitall),
                pct(p.shared), pct(p.set));
    std::printf("[profile] absolute: collideFirst=%.3f send=%.3f recv=%.3f collideLast=%.3f "
                "waitall=%.3f shared=%.3f set=%.3f\n",
                p.collideFirst, p.send, p.recv, p.collideLast, p.waitall, p.shared, p.set);
#endif
  }
  return 0;
}
