// Communication index validity: every index a comm list stores must address a
// cell of the block it belongs to.
//
// The Aver/Intp lists are built by hand in AddtoSharedAverComm() /
// AddtoSharedIntpComm() rather than through getCellIdx(), and the Intp one adds
// +1 in each direction to reach the extra interpolation source.  When the
// refined block's base is clipped by the domain edge, that +1 lands on row or
// column N of the block: the old index addressed whatever allocation followed,
// so the exchange silently corrupted the heap.
//
// This target enumerates geometries and asserts that no list stores such an
// index.  It does not run an exchange -- an out-of-range index is wrong on its
// own, whether or not the value is ever read.
//
// Layout: an Aver/Intp entry is AverSrcNum sources followed by one target, so
// the sources index the *neighbour* block and the target indexes *this* block.
// Same-level Comms are flat 1:1 pairs, even = neighbour, odd = self.
//
// Note: always look blocks up with geo.getBlock(findBlockIndex(id)), never
// geo.getBlock(id).  The latter takes a container index, and the container
// order does not match the block ids, which produces false failures.

#include <cstdio>
#include <string>
#include <vector>

#include "freelb.h"
#include "freelb.hh"
#include "../../test_util.h"

using T = FLOAT;

// One configuration's verdict.
struct Verdict {
  long pairs = 0;
  long badSrc = 0;
  long badDst = 0;
  bool bad() const { return badSrc != 0 || badDst != 0; }
};

// Same-level list: flat pairs of (neighbour cell, self cell).
template <typename BLOCK, typename GEO>
static void scanSame(const std::vector<SharedComm>& comms, const BLOCK& self, const GEO& geo,
                     Verdict& v) {
  for (const SharedComm& c : comms) {
    const auto& nbr = geo.getBlock(static_cast<int>(geo.findBlockIndex(c.SendBlockId)));
    const std::size_t n = c.SendRecvCells.size();
    for (std::size_t i = 0; i + 1 < n; i += 2) {
      ++v.pairs;
      if (c.SendRecvCells[i] >= nbr.getN()) ++v.badSrc;
      if (c.SendRecvCells[i + 1] >= self.getN()) ++v.badDst;
    }
  }
}

// Aver/Intp list: kSrc neighbour cells then one self cell.
template <typename BLOCK, typename GEO>
static void scanAverIntp(const std::vector<SharedComm>& comms, const BLOCK& self,
                         const GEO& geo, std::size_t kSrc, Verdict& v) {
  const std::size_t stride = kSrc + 1;
  for (const SharedComm& c : comms) {
    const auto& nbr = geo.getBlock(static_cast<int>(geo.findBlockIndex(c.SendBlockId)));
    const std::size_t n = c.SendRecvCells.size();
    // A list whose length is not a whole number of stencils is malformed too.
    if (n % stride != 0) {
      ++v.badSrc;
      continue;
    }
    for (std::size_t base = 0; base < n; base += stride) {
      ++v.pairs;
      for (std::size_t k = 0; k < kSrc; ++k) {
        if (c.SendRecvCells[base + k] >= nbr.getN()) {
          ++v.badSrc;
          break;
        }
      }
      if (c.SendRecvCells[base + kSrc] >= self.getN()) ++v.badDst;
    }
  }
}

template <typename BLOCK, typename GEO>
static Verdict scanBlock(const BLOCK& b, const GEO& geo, std::size_t kSrc) {
  Verdict v;
  const auto& comm = b.getCommunicator();
  scanSame(comm.Comm.Comms, b, geo, v);
  scanAverIntp(comm.Comm.AverComm, b, geo, kSrc, v);
  scanAverIntp(comm.Comm.IntpComm, b, geo, kSrc, v);
  scanAverIntp(comm.AllComm.AverComm, b, geo, kSrc, v);
  scanAverIntp(comm.AllComm.IntpComm, b, geo, kSrc, v);
  return v;
}

// MPI Intp pair invariant: for every (fine block, coarse block) pair the recv
// list on the fine side and the send list on the coarse side must describe the
// same stencil set, i.e. recv == send / kSrc, and every coarse source must
// stay inside the coarse block and every fine target inside the fine block.
// Both lists are produced by intp::BuildStencil, so a violation means the
// enumeration itself broke.
static Verdict scanMPIIntp(const BlockGeometry2D<T>& geo) {
  Verdict v;
  const std::size_t kSrc = flbcomm::AverSrcNum(2);
  for (const auto& b : geo.getBlocks()) {
    for (const DistributedComm& r : b.getCommunicator().MPIComm.IntpRecvs) {
      const auto& cblk = geo.getBlock(static_cast<int>(geo.findBlockIndex(r.TargetBlockId)));
      for (std::size_t id : r.Cells)
        if (id >= b.getN()) ++v.badDst;
      v.pairs += r.Cells.size();
    }
    for (const DistributedComm& s : b.getCommunicator().MPIComm.IntpSends) {
      const auto& fblk = geo.getBlock(static_cast<int>(geo.findBlockIndex(s.TargetBlockId)));
      // every recv entry of the fine block appears kSrc times in our send list
      std::size_t recvTotal = 0;
      for (const DistributedComm& r : fblk.getCommunicator().MPIComm.IntpRecvs)
        if (r.TargetBlockId == b.getBlockId()) recvTotal += r.Cells.size();
      if (s.Cells.size() != recvTotal * kSrc) ++v.badSrc;
      for (std::size_t id : s.Cells)
        if (id >= b.getN()) ++v.badSrc;
      v.pairs += s.Cells.size();
    }
  }
  return v;
}

static Verdict scanMPIIntp3D(const BlockGeometry3D<T>& geo) {
  Verdict v;
  const std::size_t kSrc = flbcomm::AverSrcNum(3);
  for (const auto& b : geo.getBlocks()) {
    for (const DistributedComm& r : b.getCommunicator().MPIComm.IntpRecvs) {
      const auto& cblk = geo.getBlock(static_cast<int>(geo.findBlockIndex(r.TargetBlockId)));
      for (std::size_t id : r.Cells)
        if (id >= b.getN()) ++v.badDst;
      v.pairs += r.Cells.size();
    }
    for (const DistributedComm& s : b.getCommunicator().MPIComm.IntpSends) {
      const auto& fblk = geo.getBlock(static_cast<int>(geo.findBlockIndex(s.TargetBlockId)));
      std::size_t recvTotal = 0;
      for (const DistributedComm& r : fblk.getCommunicator().MPIComm.IntpRecvs)
        if (r.TargetBlockId == b.getBlockId()) recvTotal += r.Cells.size();
      if (s.Cells.size() != recvTotal * kSrc) ++v.badSrc;
      for (std::size_t id : s.Cells)
        if (id >= b.getN()) ++v.badSrc;
      v.pairs += s.Cells.size();
    }
  }
  return v;
}

static void report(const char* name, const Verdict& v, const Verdict& m) {
  std::printf("[comm_index] %-30s pairs=%-7ld badSrc=%-5ld badDst=%-5ld mpiSrc=%-5ld mpiDst=%-5ld %s\n",
              name, v.pairs, v.badSrc, v.badDst, m.badSrc, m.badDst,
              (v.bad() || m.bad()) ? "FAIL" : "ok");
}

// ---------------------------------------------------------------------------
// 2D configurations
// ---------------------------------------------------------------------------
template <typename RefineFn>
static void run2D(const char* name, int ni, int nj, int blockLen, RefineFn&& refine) {
  AABB<T, 2> domain(Vector<T, 2>{T(0), T(0)}, Vector<T, 2>{T(ni), T(nj)});
  BlockGeometryHelper2D<T> helper(ni, nj, domain, T(1), blockLen);
  helper.forEachBlockCell([&](BasicBlock<T, 2>& b) { refine(b); });
  helper.CheckRefine();
  helper.CreateBlocks();
  helper.LoadBalancing(1);
  BlockGeometry2D<T> geo(helper);

  Verdict v;
  const std::size_t kSrc = flbcomm::AverSrcNum(2);
  for (const auto& b : geo.getBlocks()) {
    const Verdict p = scanBlock(b, geo, kSrc);
    v.pairs += p.pairs;
    v.badSrc += p.badSrc;
    v.badDst += p.badDst;
  }
  const Verdict m = scanMPIIntp(geo);
  report(name, v, m);
  CHECK_EQ(v.badSrc, 0);
  CHECK_EQ(v.badDst, 0);
  CHECK_EQ(m.badSrc, 0);
  CHECK_EQ(m.badDst, 0);
}

template <typename RefineFn>
static void run3D(const char* name, int n, int blockLen, RefineFn&& refine) {
  AABB<T, 3> domain(Vector<T, 3>{T(0), T(0), T(0)}, Vector<T, 3>{T(n), T(n), T(n)});
  BlockGeometryHelper3D<T> helper(n, n, n, domain, T(1), blockLen);
  helper.forEachBlockCell([&](BasicBlock<T, 3>& b) { refine(b); });
  helper.CheckRefine();
  helper.CreateBlocks();
  helper.LoadBalancing(1);
  BlockGeometry3D<T> geo(helper);

  Verdict v;
  const std::size_t kSrc = flbcomm::AverSrcNum(3);
  for (const auto& b : geo.getBlocks()) {
    const Verdict p = scanBlock(b, geo, kSrc);
    v.pairs += p.pairs;
    v.badSrc += p.badSrc;
    v.badDst += p.badDst;
  }
  const Verdict m = scanMPIIntp3D(geo);
  report(name, v, m);
  CHECK_EQ(v.badSrc, 0);
  CHECK_EQ(v.badDst, 0);
  CHECK_EQ(m.badSrc, 0);
  CHECK_EQ(m.badDst, 0);
}

int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);

  // --- 2D: no refinement, so no cross-level list exists at all --------------
  run2D("2D uniform", 64, 64, 16, [](BasicBlock<T, 2>&) {});

  // --- 2D: refined island, inset (the transition ring is generated) --------
  {
    AABB<T, 2> in(Vector<T, 2>{T(16), T(16)}, Vector<T, 2>{T(48), T(48)});
    run2D("2D lv1 inset", 64, 64, 16, [&](BasicBlock<T, 2>& b) {
      if (isOverlapped(b, in)) b.refine(1);
    });
  }
  // --- 2D: island touching the edge, so the refined base gets clipped ------
  {
    AABB<T, 2> in(Vector<T, 2>{T(0), T(0)}, Vector<T, 2>{T(32), T(32)});
    run2D("2D lv1 at corner", 64, 64, 16, [&](BasicBlock<T, 2>& b) {
      if (isOverlapped(b, in)) b.refine(1);
    });
  }
  {
    AABB<T, 2> in(Vector<T, 2>{T(24), T(24)}, Vector<T, 2>{T(40), T(40)});
    run2D("2D lv2 inset", 64, 64, 16, [&](BasicBlock<T, 2>& b) {
      if (isOverlapped(b, in)) b.refine(2);
    });
  }
  {
    AABB<T, 2> in(Vector<T, 2>{T(0), T(0)}, Vector<T, 2>{T(24), T(24)});
    run2D("2D lv2 at corner", 64, 64, 16, [&](BasicBlock<T, 2>& b) {
      if (isOverlapped(b, in)) b.refine(2);
    });
  }
  {
    AABB<T, 2> in(Vector<T, 2>{T(24), T(24)}, Vector<T, 2>{T(40), T(40)});
    run2D("2D lv2 centred 96", 96, 96, 16, [&](BasicBlock<T, 2>& b) {
      if (isOverlapped(b, in)) b.refine(2);
    });
  }
  {
    AABB<T, 2> in(Vector<T, 2>{T(32), T(32)}, Vector<T, 2>{T(48), T(48)});
    run2D("2D lv3 inset", 64, 64, 16, [&](BasicBlock<T, 2>& b) {
      if (isOverlapped(b, in)) b.refine(3);
    });
  }
  {
    AABB<T, 2> in(Vector<T, 2>{T(8), T(8)}, Vector<T, 2>{T(56), T(56)});
    run2D("2D ring refined", 64, 64, 16, [&](BasicBlock<T, 2>& b) {
      if (!isOverlapped(b, in)) b.refine(1);
    });
  }

  // --- 3D: same shapes, 8 sources per stencil ------------------------------
  {
    AABB<T, 3> in(Vector<T, 3>{T(24), T(24), T(24)}, Vector<T, 3>{T(40), T(40), T(40)});
    run3D("3D lv1 inset", 64, 16, [&](BasicBlock<T, 3>& b) {
      if (isOverlapped(b, in)) b.refine(1);
    });
  }
  {
    AABB<T, 3> in(Vector<T, 3>{T(24), T(24), T(24)}, Vector<T, 3>{T(40), T(40), T(40)});
    run3D("3D lv2 inset", 64, 16, [&](BasicBlock<T, 3>& b) {
      if (isOverlapped(b, in)) b.refine(2);
    });
  }
  {
    AABB<T, 3> in(Vector<T, 3>{T(0), T(0), T(0)}, Vector<T, 3>{T(24), T(24), T(24)});
    run3D("3D lv2 at corner", 64, 16, [&](BasicBlock<T, 3>& b) {
      if (isOverlapped(b, in)) b.refine(2);
    });
  }

  return flbtest::summary("block/comm_index");
}
