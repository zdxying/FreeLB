#pragma once
// ---------------------------------------------------------------------------
// POP storage strategy tags for cudev::Cell.
//
// Kept in their own header because two unrelated files need them:
//   * lbm/equilibrium.h  forward declares cudev::Cell and must supply the
//     default argument, and a default argument may only be given on the FIRST
//     declaration -- cell.h includes the generated .ur.h files, which include
//     equilibrium.h, so equilibrium.h is always seen first.
//   * data_struct/cell.h defines the class.
//
// They live outside __CUDACC__ deliberately: they are empty tag structs, so a
// host build seeing them costs nothing, and keeping the declaration and the
// default argument on the same side of the preprocessor boundary means a
// CPU-visible mention of cudev::Cell<T, LatSet, TypePack> resolves to the same
// type in a host build as in a device build.
//
// See docs/CELL_POLICY_PLAN.md.
// ---------------------------------------------------------------------------

namespace cudev {

// Direct access to the POP container; nothing is cached on the cell.
//
// block_size must stay in step with THREADS_PER_BLOCK in utils/cuda_device.h,
// which is why it is spelled out rather than referenced: that macro only exists
// on CUDA builds and this header has to be visible to the host compiler too.
struct DirectPop {
  static constexpr unsigned int block_size = 32;
};

// Keep the q populations in registers for the duration of the cell dynamics.
// Resolves the q element addresses once in the constructor and writes them back
// in flush().  The register footprint (96 regs for D3Q19) wants a multiple of the
// warp size rather than the streaming baseline's block size, hence 128.
struct RegPop {
  static constexpr unsigned int block_size = 128;
};

// ---------------------------------------------------------------------------
// A third strategy is expected to be needed, and this design is what makes it
// cheap: `LazyRegPop`, which does not move the q populations in the constructor
// but only after the flag dispatch has ruled the cell active.
//
// The current eager load is the design's main weakness.  It happens *before*
// the flag is inspected, so a void cell still pays the full q loads and q
// stores -- 152 B/cell for D3Q19 at FP32 -- for dynamics that then do nothing
// with them.  See the void-cell analysis in docs/GPU_MULTIBLOCK_PLAN.md (8.1,
// 8.2).
//
// When that gets implemented it should be a third tag plus a third PopCache
// specialisation and nothing else: csegen forwards POPPOLICY generically, so no
// new .ur.h specialisation and no csegen change are needed, and the CSE
// arithmetic is shared by all three strategies.
// ---------------------------------------------------------------------------

}  // namespace cudev
