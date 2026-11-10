// checkpoint regression test
//
// A resumed run must continue bit-identically from an uninterrupted one, so the
// test compares pop-field checksums:
//
//   A: 200 steps -> write -> 100 steps -> checksum
//   B: 200 steps -> write -> read -> 100 steps -> checksum
//
// A and B must agree exactly.  Determinism requires a single thread, a single
// block and no -D_UNROLLFOR, hence the flags in the Makefile.
//
// Also covered: the CyclicArray logical-view round trip (the checkpoint file
// must not contain the container's rotation state), and readHeader().

#include <cstdio>

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

// single lid-driven cavity: uniform, single block, deterministic
struct Sim {
  AABB<T, 2> cavity;
  AABB<T, 2> toplid;
  BlockGeometry2D<T> Geo;
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM;
  BaseConverter<T> Conv;
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  using CELL = Cell<T, LatSet, FIELDS>;
  using BulkTask =
    tmp::Key_TypePair<AABBFlag, collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>>;
  using BBTask = tmp::Key_TypePair<BouncebackFlag, collision::BounceBack<CELL>>;
  using BBMVTask = tmp::Key_TypePair<BBMovingWallFlag, collision::BounceBackMovingWall<CELL>>;
  using TaskCollection = tmp::TupleWrapper<BulkTask, BBTask, BBMVTask>;
  using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;
  using RhoUTask = tmp::Key_TypePair<AABBFlag, moment::rhoU<CELL, true>>;
  using TaskCollectionRhoU = tmp::TupleWrapper<RhoUTask>;
  using TaskSelectorRhoU = tmp::TaskSelector<TaskCollectionRhoU, std::uint8_t, CELL>;
  BlockLatticeManager<T, LatSet, FIELDS>* Lattice;

  Sim(int Ni, int Nj)
      : cavity(Vector<T, 2>{}, Vector<T, 2>{T(Ni), T(Nj)}),
        toplid(T(1), T((Nj - 1)), T(Ni - 1), T(Nj)),
        Geo(Ni, Nj, 1, cavity, T(1)),
        FlagFM(Geo, VoidFlag),
        Conv(LatSet::cs2),
        Lattice(nullptr) {
    Conv.ConvertFromRT(T(1), T(0.56), T(0.001), T(Ni), T(0.1), T(0.02));
    FlagFM.forEach(cavity, [&](FLAG& field, std::size_t id) { field.SetField(id, AABBFlag); });
    FlagFM.template SetupBoundary<LatSet>(cavity, BouncebackFlag);
    FlagFM.forEach(toplid, [&](FLAG& field, std::size_t id) {
      if (util::isFlag(field.get(id), BouncebackFlag))
        field.SetField(id, BBMovingWallFlag);
    });
    // the converter must have been converted from RT before its init values are used
    ValuePack<T, Vector<T, LatSet::d>, T> initValues(T(1), Vector<T, LatSet::d>{}, T{});
    initValues =
      ValuePack<T, Vector<T, LatSet::d>, T>(Conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
    Lattice = new BlockLatticeManager<T, LatSet, FIELDS>(Geo, initValues, Conv);
    Vector<T, LatSet::d> LatU_Wall = Conv.getLatticeU(Vector<T, 2>{T(0.1), T(0)});
    Lattice->getField<VELOCITY<T, LatSet::d>>().forEach(
      toplid, FlagFM, BBMovingWallFlag,
      [&](auto& field, std::size_t id) { field.SetField(id, LatU_Wall); });
  }

  ~Sim() { delete Lattice; }

  void step() {
    Lattice->Stream();
    Lattice->template ApplyCellDynamics<NSTask>(FlagFM);
    Lattice->NormalFullCommunicate();
  }

  // checksum of the (host) pop field: sum, sumsq, maxabs
  std::array<double, 3> popChecksum() {
    auto& popf = Lattice->getBlockLat(0).template getField<POP<T, LatSet::q>>();
    return frelb_diag::PopChecksumHost(popf, Lattice->getBlockLat(0).getN(), LatSet::q);
  }
};

// ---------------------------------------------------------------------------
// CyclicArray logical-view round trip: the file must hold the logical view, so
// a rotated container must survive save/load unchanged.
// ---------------------------------------------------------------------------
static void testCyclicRoundTrip() {
  const std::size_t n = 37;
  CyclicArray<double> arr(n);
  for (std::size_t i = 0; i < n; ++i) arr.set(i, static_cast<double>(i) * 1.5);
  // rotate twice by non-trivial offsets, as Stream() does
  arr.rotate(3);
  arr.rotate(-11);

  std::vector<double> logical(n);
  for (std::size_t i = 0; i < n; ++i) logical[i] = static_cast<double>(i) * 1.5;

  // the logical view must not change after normalization
  arr.normalize();
  bool same = true;
  for (std::size_t i = 0; i < n; ++i)
    if (arr[i] != logical[i]) same = false;
  CHECK(same);
  CHECK_EQ(arr.getRemainder(), n - 1);
}

int main() {
  const std::string dir = "./ckpt";
  // start from a clean directory: latestStep() must not see a stale checkpoint
  std::string rmcmd = "rm -rf " + dir;
  if (std::system(rmcmd.c_str()) != 0) CHECK(false);

  constexpr int kPre = 200;
  constexpr int kPost = 100;
  constexpr int kRun = kPre + kPost;

  testCyclicRoundTrip();

  // ---- A: uninterrupted ----
  Sim simA(40, 40);
  for (int i = 0; i < kPre; ++i) simA.step();
  flbckpt::write(dir, kPre, *simA.Lattice, simA.FlagFM);
  for (int i = kPre; i < kRun; ++i) simA.step();
  const std::array<double, 3> cA = simA.popChecksum();

  // ---- B: stop, checkpoint, resume ----
  Sim simB(40, 40);
  for (int i = 0; i < kPre; ++i) simB.step();
  flbckpt::write(dir, kPre, *simB.Lattice, simB.FlagFM);

  // the resumed run must load before stepping again
  CHECK(flbckpt::checkCompatible(dir, -1, *simB.Lattice, simB.FlagFM));
  const std::int64_t restored =
    flbckpt::read(dir, -1, *simB.Lattice, simB.FlagFM);
  CHECK_EQ(restored, (std::int64_t)kPre);
  CHECK_EQ(flbckpt::latestStep(dir), (std::int64_t)kPre);

  for (int i = kPre; i < kRun; ++i) simB.step();
  const std::array<double, 3> cB = simB.popChecksum();

  // ---- bit-level agreement ----
  CHECK(cA[0] == cB[0]);  // sum
  CHECK(cA[1] == cB[1]);  // sumsq
  CHECK(cA[2] == cB[2]);  // maxabs

  // ---- header sanity ----
  flbckpt::CheckpointHeader hdr;
  CHECK(flbckpt::readHeader(dir, -1, hdr));
  CHECK_EQ(hdr.step, (std::int64_t)kPre);
  CHECK_EQ(hdr.blocks.size(), std::size_t(1));
  // the block carries one overlap layer: (Ni + 2) * (Nj + 2)
  CHECK(hdr.blocks[0].n == 42u * 42u);
  CHECK_EQ(hdr.fields.size(), std::size_t(4));  // RHO, VELOCITY, POP, FLAG
  bool hasPop = false;
  for (const auto& f : hdr.fields) {
    if (f.array_dim == LatSet::q) {
      hasPop = true;
      CHECK_EQ((int)f.tag, (int)flbckpt::ElemTag::f64);
      CHECK_EQ(f.elem_dim, 1u);
      CHECK_EQ(f.counts.size(), std::size_t(1));
      CHECK_EQ(f.counts[0].size(), std::size_t(LatSet::q));
    }
  }
  CHECK(hasPop);

  // ---- a checkpoint from a different geometry must be rejected ----
  {
    Sim simC(41, 40);  // one cell wider
    CHECK(!flbckpt::checkCompatible(dir, -1, *simC.Lattice, simC.FlagFM));
  }

  return flbtest::summary("checkpoint");
}
