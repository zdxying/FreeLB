// plain communication benchmark, 2D
//
// Measures the three lattice exchanges (partial / full / all) and the field
// exchange, each as a complete step and as the exchange on its own.  For every
// path it also prints the message layout: how many communicators the geometry
// produces and how many MPI messages the coalesced plan posts for them.
//
//   ./plain_comm_bench2d.exe [Ni] [blocksPerSide] [steps]
//   mpirun -np 2 ./plain_comm_bench2d.exe 120 10 2000
//
// The block grid is fixed (CreateBlocks is rank independent); LoadBalancing
// decides which rank owns which block, so the same command works for any -np.
// With one rank there is no MPI traffic at all, which makes a useful baseline.

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <set>

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

// Lid-driven cavity on a fixed block grid, same setup as the other benchmarks.
struct Sim {
  AABB<T, 2> cavity;
  AABB<T, 2> toplid;
  BlockGeometryHelper2D<T> GeoHelper;
  BlockGeometry2D<T>* Geo;
  BlockFieldManager<FLAG, T, LatSet::d>* FlagFM;
  BaseConverter<T> Conv;
  BlockLatticeManager<T, LatSet, FIELDS>* Lattice;

  Sim(int Ni, int Nj, int bx, int by)
      : cavity(Vector<T, 2>{}, Vector<T, 2>{T(Ni), T(Nj)}),
        toplid(T(1), T(Nj - 1), T(Ni - 1), T(Nj)),
        GeoHelper(Ni, Nj, cavity, T(1), 20),
        Geo(nullptr),
        FlagFM(nullptr),
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
    ValuePack<T, Vector<T, LatSet::d>, T> iv(Conv.getLatRhoInit(),
                                             Vector<T, LatSet::d>{}, T{});
    Lattice = new BlockLatticeManager<T, LatSet, FIELDS>(*Geo, iv, Conv);
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

  std::size_t localCells() const {
    std::size_t n = 0;
    for (auto& bl : Lattice->getBlockLats()) n += bl.getN();
    return n;
  }
};

// communicators on this rank, distinct peer ranks, and cells per step
struct Layout {
  std::size_t comms{};
  std::size_t messages{};
  std::size_t cells{};
};

template <typename SendList>
static Layout layoutOf(Sim& sim, SendList&& sendsOf) {
  Layout l;
  std::set<int> peers;
  for (auto& bl : sim.Lattice->getBlockLats()) {
    for (const DistributedComm& c : sendsOf(bl.getBlock().getCommunicator())) {
      ++l.comms;
      l.cells += c.Cells.size();
      peers.insert(c.TargetRank);
    }
  }
  l.messages = peers.size();
  return l;
}

static void printLayout(const char* name, const Layout& l) {
  std::printf("[layout] %-8s communicators/step=%-5zu messages/step=%-3zu cells/step=%zu\n", name,
              l.comms, l.messages, l.cells);
}

int main(int argc, char** argv) {
  mpi().init(&argc, &argv);
  const int Ni = argc > 1 ? std::atoi(argv[1]) : 120;
  const int bx = argc > 2 ? std::atoi(argv[2]) : 10;
  const int steps = argc > 3 ? std::atoi(argv[3]) : 2000;

  Sim sim(Ni, Ni, bx, bx);

  if (mpi().getRank() == 0) {
    std::printf("[config] Ni=%d blocks=%dx%d ranks=%d steps=%d\n", Ni, bx, bx, mpi().getSize(),
                steps);
  }

  const Layout dirL = layoutOf(sim, [](const Communicator& c) -> const std::vector<DistributedComm>& {
    return c.DirSends;
  });
  const Layout fullL = layoutOf(sim, [](const Communicator& c) -> const std::vector<DistributedComm>& {
    return c.MPIComm.Sends;
  });
  const Layout allL = layoutOf(sim, [](const Communicator& c) -> const std::vector<DistributedComm>& {
    return c.AllMPIComm.Sends;
  });
  if (mpi().getRank() == 0) {
    printLayout("partial", dirL);
    printLayout("full", fullL);
    printLayout("all", allL);
    printLayout("field", fullL);
  }

  auto run = [&](const char* name, auto&& step, auto&& comm) {
    for (int i = 0; i < 50; ++i) step();  // warmup
    mpi().barrier();
    double bestStep = 1e30;
    for (int rep = 0; rep < 3; ++rep) {
      mpi().barrier();
      const auto t0 = std::chrono::steady_clock::now();
      for (int i = 0; i < steps; ++i) step();
      const auto t1 = std::chrono::steady_clock::now();
      bestStep = std::min(bestStep, std::chrono::duration<double>(t1 - t0).count());
      mpi().barrier();
    }
    double bestComm = 1e30;
    for (int rep = 0; rep < 3; ++rep) {
      mpi().barrier();
      const auto t0 = std::chrono::steady_clock::now();
      for (int i = 0; i < steps; ++i) comm();
      const auto t1 = std::chrono::steady_clock::now();
      bestComm = std::min(bestComm, std::chrono::duration<double>(t1 - t0).count());
      mpi().barrier();
    }
    if (mpi().getRank() == 0) {
      std::printf("[result] %-8s step=%8.1f us  comm=%7.1f us  MLUPS=%6.2f\n", name,
                  1e6 * bestStep / steps, 1e6 * bestComm / steps,
                  double(sim.localCells()) * steps / bestStep / 1e6);
    }
  };

  run("partial",
      [&] {
        sim.Lattice->Stream();
        sim.Lattice->template ApplyCellDynamics<NSTask>(*sim.FlagFM);
        sim.Lattice->NormalCommunicate();
      },
      [&] { sim.Lattice->NormalCommunicate(); });

  run("full",
      [&] {
        sim.Lattice->Stream();
        sim.Lattice->template ApplyCellDynamics<NSTask>(*sim.FlagFM);
        sim.Lattice->NormalFullCommunicate();
      },
      [&] { sim.Lattice->NormalFullCommunicate(); });

  run("all",
      [&] {
        sim.Lattice->Stream();
        sim.Lattice->template ApplyCellDynamics<NSTask>(*sim.FlagFM);
        sim.Lattice->NormalAllCommunicate();
      },
      [&] { sim.Lattice->NormalAllCommunicate(); });

  run("field", [&] { sim.FlagFM->NormalCommunicate(); },
      [&] { sim.FlagFM->NormalCommunicate(); });

  return 0;
}
