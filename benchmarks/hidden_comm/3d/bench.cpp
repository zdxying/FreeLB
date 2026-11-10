// hidden communication benchmark, 3D
//
// Same measurement as the 2D benchmark: throughput of
// ApplyCellDynamicsWithHiddenComm, the message layout, and a per-phase
// breakdown behind _HIDDEN_PROFILE.
//
//   ./hidden_comm_bench3d.exe [N] [blocksPerSide] [steps]
//   mpirun -np 2 ./hidden_comm_bench3d.exe 100 4 300
//
// The domain is a cube of N^3 cells split into blocksPerSide^3 blocks, so the
// block edge is N/blocksPerSide.  The block grid is rank independent;
// LoadBalancing decides which rank owns which block.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <map>

#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;
using LatSet = D3Q19<T>;

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

// Lid-driven cavity on a fixed block grid.  Only the hidden-communication step
// is timed; the setup mirrors examples/cavity3d.
struct Sim {
  AABB<T, 3> cavity;
  AABB<T, 3> toplid;
  BlockGeometryHelper3D<T> GeoHelper;
  BlockGeometry3D<T>* Geo;
  BlockFieldManager<FLAG, T, LatSet::d>* FlagFM;
  BaseConverter<T> Conv;
  BlockLatticeManager<T, LatSet, FIELDS>* Lattice;

  Sim(int N, int b)
      : cavity(Vector<T, 3>{}, Vector<T, 3>{T(N), T(N), T(N)}),
        toplid(Vector<T, 3>{T(1), T(1), T(N - 1)}, Vector<T, 3>{T(N - 1), T(N - 1), T(N)}),
        GeoHelper(N, N, N, cavity, T(1), 20),
        Geo(nullptr),
        FlagFM(nullptr),
        Conv(LatSet::cs2),
        Lattice(nullptr) {
    GeoHelper.CreateBlocks(b, b, b);
    GeoHelper.LoadBalancing(mpi().getSize());
    Geo = new BlockGeometry3D<T>(GeoHelper);
    Conv.ConvertFromRT(T(1), T(0.56), T(0.001), T(N), T(0.1), T(0.02));
    FlagFM = new BlockFieldManager<FLAG, T, LatSet::d>(*Geo, VoidFlag);
    FlagFM->forEach(cavity,
                    [&](FLAG& field, std::size_t id) { field.SetField(id, AABBFlag); });
    FlagFM->template SetupBoundary<LatSet>(cavity, BouncebackFlag);
    FlagFM->forEach(toplid, [&](FLAG& field, std::size_t id) {
      if (util::isFlag(field.get(id), BouncebackFlag))
        field.SetField(id, BBMovingWallFlag);
    });
    ValuePack<T, Vector<T, LatSet::d>, T> iv(Conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
    Lattice = new BlockLatticeManager<T, LatSet, FIELDS>(*Geo, iv, Conv);
    Vector<T, LatSet::d> LatU_Wall = Conv.getLatticeU(Vector<T, 3>{T(0.1), T(0), T(0)});
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
    Lattice->template ApplyCellDynamicsWithHiddenComm<NSTask>(*FlagFM);
  }

  std::size_t localCells() const {
    std::size_t n = 0;
    for (auto& bl : Lattice->getBlockLats()) n += bl.getN();
    return n;
  }
};

int main(int argc, char** argv) {
  mpi().init(&argc, &argv);
  const int N = argc > 1 ? std::atoi(argv[1]) : 100;
  const int b = argc > 2 ? std::atoi(argv[2]) : 4;
  const int steps = argc > 3 ? std::atoi(argv[3]) : 300;

  Sim sim(N, b);

  // Message layout: the communicators the geometry produces, and how they are
  // coalesced into one message per target rank.
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
    std::printf("[config] N=%d blocks=%dx%dx%d ranks=%d steps=%d\n", N, b, b, b, mpi().getSize(),
                steps);
    std::printf("[layout] local blocks=%zu  communicators/step=%zu  cells/comm=%.1f\n", blocks,
                entries, entries ? double(cells) / entries : 0.0);
    std::printf("[layout] comm size histogram (cells x count):");
    for (auto& kv : hist) std::printf(" %zux%zu", kv.first, kv.second);
    std::printf("\n");
    std::printf("[layout] coalesced messages/step=%zu  cells/message=%.1f\n", perRank.size(),
                perRank.empty() ? 0.0 : double(cells) / perRank.size());
  }

  for (int i = 0; i < 20; ++i) sim.step();  // warmup
  mpi().barrier();

  double best = 1e30;
  for (int rep = 0; rep < 3; ++rep) {
    mpi().barrier();
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < steps; ++i) sim.step();
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
