// checkpoint regression test, GPU path
//
// Same contract as tests/feature/checkpoint: a resumed run must continue
// bit-identically from an uninterrupted one.  On the GPU the interesting part
// is the device mirror:
//
//   - write() must pull the device populations back (copyToHost) before dumping
//   - read() must normalize the container rotation, re-upload the device mirror
//     object (cudev::CyclicArray caches start[] INSIDE the object, which
//     copyToDevice does NOT push) and push the scalars back
//
// The checksums are taken with PopChecksumDevice, i.e. on the device itself, so
// a stale host copy or a misaligned window cannot pass the test.
//
// Build this directory twice: default (FP32 storage) and with FP16=1
// (__half storage, ElemTag::half -- a code path the FP32 build never reaches).

#include <sys/wait.h>

#include "freelb.h"
#include "freelb.hh"
#include "../../test_util.h"
#include "utils/field_checksum.h"

using T = float;
using LatSet = D2Q9<T>;

static constexpr std::uint8_t VoidFlag = std::uint8_t(1);
static constexpr std::uint8_t AABBFlag = std::uint8_t(2);
static constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
static constexpr std::uint8_t BBMovingWallFlag = std::uint8_t(8);

struct Sim {
  AABB<T, 2> cavity;
  AABB<T, 2> toplid;
  BlockGeometry2D<T> Geo;
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM;
  BaseConverter<T> Conv;
  using FIELDS = TypePack<RHO<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>>;
  using cudevFIELDS = typename ExtractCudevFieldPack<FIELDS>::cudev_pack;
  // The POP storage strategy is a template parameter of cudev::Cell.  Under
  // FREELB_POP_FP16 the loop-based moment schemes do not compile with DirectPop
  // (moment.h computes `latset::c * cell[i]` and there is no Vector * PopRef);
  // the CSE-generated moment.ur.h hoists each cell[i] into a const T first, so
  // DirectPop+FP16 does build with -D_UNROLLFOR (benchmarks/cavity3d_cu does).
  // This target compiles the loop path, so it selects RegCell; FP32 keeps the
  // DirectPop policy used by examples/cavity2d_cu.
#if defined(FREELB_POP_FP16)
  using DevCELL = cudev::RegCell<T, LatSet, cudevFIELDS>;
#else
  using DevCELL = cudev::Cell<T, LatSet, cudevFIELDS>;
#endif
  using CELL = DevCELL;
  using BulkTask =
    tmp::Key_TypePair<AABBFlag, collision::BGK<moment::rhoU<DevCELL>, equilibrium::SecondOrder<DevCELL>>>;
  using BBTask = tmp::Key_TypePair<BouncebackFlag, collision::BounceBack<DevCELL>>;
  using BBMVTask = tmp::Key_TypePair<BBMovingWallFlag, collision::BounceBackMovingWall<DevCELL>>;
  using TaskCollection = tmp::TupleWrapper<BulkTask, BBTask, BBMVTask>;
  using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, DevCELL>;
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
    ValuePack<T, Vector<T, LatSet::d>, T> initValues(T(1), Vector<T, LatSet::d>{}, T{});
    initValues =
      ValuePack<T, Vector<T, LatSet::d>, T>(Conv.getLatRhoInit(), Vector<T, LatSet::d>{}, T{});
    Lattice = new BlockLatticeManager<T, LatSet, FIELDS>(Geo, initValues, Conv);
    Vector<T, LatSet::d> LatU_Wall = Conv.getLatticeU(Vector<T, 2>{T(0.1), T(0)});
    Lattice->getField<VELOCITY<T, LatSet::d>>().forEach(
      toplid, FlagFM, BBMovingWallFlag,
      [&](auto& field, std::size_t id) { field.SetField(id, LatU_Wall); });
    // do not forget to copy to device
    FlagFM.copyToDevice();
    Lattice->getField<POP<T, LatSet::q>>().copyToDevice();
    Lattice->getField<RHO<T>>().copyToDevice();
    Lattice->getField<VELOCITY<T, LatSet::d>>().copyToDevice();
  }

  ~Sim() { delete Lattice; }

  // one full step: collide then stream (same order as examples/cavity2d_cu)
  void step() {
    Lattice->CuDevApplyCellDynamics<NSTask, DevCELL>(FlagFM);
    Lattice->CuDevStream();
    cudaDeviceSynchronize();
  }

  // device-side pop checksum: sum, sumsq, maxabs
  std::array<double, 3> popChecksum() {
    return frelb_diag::PopChecksumDevice(Lattice->getBlockLat(0));
  }

  // host-side pop checksum (after copyToHost)
  std::array<double, 3> popChecksumHost() {
    auto& popf = Lattice->getBlockLat(0).template getField<POP<T, LatSet::q>>();
    return frelb_diag::PopChecksumHost(popf, popf.getField(0).size(), LatSet::q);
  }
};

int main(int argc, char** argv) {
  const std::string dir = "./ckpt";
  // "reject" mode: this binary is built with a different POP storage type than
  // the one that wrote ./ckpt (run `make` then `make FP16=1`); the checkpoint
  // must be rejected.  Cross-storage restore is refused rather than converted.
  const bool rejectMode = argc > 1 && std::string(argv[1]) == "reject";

  if (rejectMode) {
    Sim sim(48, 48);
    const bool accepted = flbckpt::checkCompatible("./ckpt", -1, *sim.Lattice, sim.FlagFM);
    std::cout << "[reject mode] this build stores POP as "
#if defined(FREELB_POP_FP16)
              << "half"
#else
              << "f32"
#endif
              << "; checkpoint accepted: " << (accepted ? "yes" : "no") << std::endl;
    return accepted ? 1 : 0;
  }

  {
    std::string rmcmd = "rm -rf ./ckpt";
    if (std::system(rmcmd.c_str()) != 0) CHECK(false);
  }

  constexpr int kPre = 200;
  constexpr int kPost = 100;
  constexpr int kRun = kPre + kPost;

  // ---- A: uninterrupted ----
  Sim simA(48, 48);
  for (int i = 0; i < kPre; ++i) simA.step();
  flbckpt::write(dir, kPre, *simA.Lattice, simA.FlagFM);
  for (int i = kPre; i < kRun; ++i) simA.step();
  const std::array<double, 3> cA = simA.popChecksum();

  // ---- B: stop, checkpoint, resume ----
  Sim simB(48, 48);
  for (int i = 0; i < kPre; ++i) simB.step();
  flbckpt::write(dir, kPre, *simB.Lattice, simB.FlagFM);

  CHECK(flbckpt::checkCompatible(dir, -1, *simB.Lattice, simB.FlagFM));
  const std::int64_t restored = flbckpt::read(dir, -1, *simB.Lattice, simB.FlagFM);
  CHECK_EQ(restored, (std::int64_t)kPre);
  CHECK_EQ(flbckpt::latestStep(dir), (std::int64_t)kPre);

  for (int i = kPre; i < kRun; ++i) simB.step();
  const std::array<double, 3> cB = simB.popChecksum();

  // ---- bit-level agreement ----
  // NB: PopChecksumDevice accumulates with atomicAdd, whose order is not
  // deterministic, so it can only be compared with a tolerance.  The bit-exact
  // comparison is done on the host after copyToHost (same serial order as
  // tests/feature/checkpoint).
  {
    simA.Lattice->getBlockLat(0).template getField<POP<T, LatSet::q>>().copyToHost();
    simB.Lattice->getBlockLat(0).template getField<POP<T, LatSet::q>>().copyToHost();
    const std::array<double, 3> hA = simA.popChecksumHost();
    const std::array<double, 3> hB = simB.popChecksumHost();
    CHECK(hA[0] == hB[0]);  // sum
    CHECK(hA[1] == hB[1]);  // sumsq
    CHECK(hA[2] == hB[2]);  // maxabs
  }
  // the device-side checksum must still agree to round-off
  CHECK_NEAR(cA[0], cB[0], std::abs(cA[0]) * 1e-12 + 1e-9);
  CHECK_NEAR(cA[1], cB[1], std::abs(cA[1]) * 1e-12 + 1e-9);
  CHECK(cA[2] == cB[2]);  // maxabs is order independent

  // ---- the restored device state must equal the file, not a stale host copy ----
  {
    simB.Lattice->getBlockLat(0).template getField<POP<T, LatSet::q>>().copyToHost();
    const std::array<double, 3> cHost = simB.popChecksumHost();
    // after a further 100 steps the host value must match the device value to
    // round-off (the host copy is only used for the checkpoint, not for compute)
    CHECK_NEAR(cHost[0], cB[0], std::abs(cB[0]) * 1e-12 + 1e-9);
    CHECK_NEAR(cHost[1], cB[1], std::abs(cB[1]) * 1e-12 + 1e-9);
  }

  // ---- header: element type must match this build's POP storage ----
  flbckpt::CheckpointHeader hdr;
  CHECK(flbckpt::readHeader(dir, -1, hdr));
  CHECK_EQ(hdr.step, (std::int64_t)kPre);
  for (const auto& f : hdr.fields) {
    if (f.array_dim == LatSet::q) {
#if defined(FREELB_POP_FP16)
      CHECK_EQ((int)f.tag, (int)flbckpt::ElemTag::half);
#else
      CHECK_EQ((int)f.tag, (int)flbckpt::ElemTag::f32);
#endif
      CHECK_EQ(f.elem_dim, 1u);
    }
  }

  // ---- a checkpoint from a different geometry must be rejected ----
  {
    Sim simC(48, 49);  // one cell taller
    CHECK(!flbckpt::checkCompatible(dir, -1, *simC.Lattice, simC.FlagFM));
  }

  // ---- cross-storage rejection: hand this build's checkpoint to an FP32 build ----
  // It must accept it when the parent is FP32 (same storage) and reject it when
  // the parent is FP16 (different storage): ElemTag f32 vs half.
  if (std::system("test -x ./checkpoint_f32.exe") == 0) {
    const int status = std::system("./checkpoint_f32.exe reject");
    // system() returns the wait status, not the exit code
    const int rc = WEXITSTATUS(status);
#if defined(FREELB_POP_FP16)
    CHECK(rc == 0);  // the f32 build must reject a half checkpoint
#else
    CHECK(rc == 1);  // the f32 build must accept an f32 checkpoint
#endif
  } else {
    std::cout << "[note] checkpoint_f32.exe not built; copy the FP32 binary to "
                 "checkpoint_f32.exe to enable the cross-storage check"
              << std::endl;
  }

  return flbtest::summary("checkpoint_cuda");
}