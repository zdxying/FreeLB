# Case generator and operator DSL (tools/casegen): decision record

Status: closed (2026-11-10); this is the design record for `tools/casegen`. Phase 1–4 and the geometry patch binding (B4) are landed; usage and field tables live in `README.md` | Updated: 2026-11-10 | Target: record why the generator has this shape

> Closed: `README.md` and the code explain usage and implementation on their own, so this file
> records only the decisions and rejected alternatives behind them. It is updated only when one of
> those decisions changes. The admission rule for a new operator, field or catalogue case (it must
> have a golden or one manual verification record) is written in the README and is not repeated here.

## 1. Boundary and ownership

The generator only assembles: the field set, geometry painting, boundaries, operator selection, the
main-loop skeleton, the `.ini` and the `Makefile`. Physics kernels come from the built-in catalogue
or from `.op` operators. The output is a one-shot scaffold; the user owns it and may keep editing it
by hand.

Apart from the periodic comm storage fix recorded in §2, casegen does not touch `src/`; it does not
touch `make.mk` or `tests/run_tests.sh`; it does not do CA, free surface or AMR.

The opening is deliberately narrow: structure is compile-time type assembly and the run-time
parameters already live in the `.ini`; moving both into a new language would save a few hundred lines
of transcription at the cost of reimplementing type checking, diagnostics and a permanent
maintenance surface.

## 2. Key decisions

- **Python emits C++ directly; it does not reuse csegen.** The input basenames of `csegen` are
  hard-coded to `moment.h/equilibrium.h/force.h`; generalising it first would drag in the submodule's
  CI and its line-number contract together. An operator body is per-cell arithmetic, and the C++
  compiler already does constant folding and loop unrolling; the DSL guarantees semantics and
  `file:line` diagnostics first, and leaves CSE as a later optimisation.
- **Two syntaxes, linked by name, not merged into one file.** The structure deck is types/lists/order,
  the kernel is expressions; the `third_party/cse` front end models only 7 statement kinds, an
  unmodelled top-level token is passed through unchanged, and a single `[section]` is enough to make
  a whole block be skipped silently. The `.case` writes only `[operators] name = path`; the `.op` is
  parsed and diagnosed independently.
- **`box` is an AABB endpoint pair, not a closed range of cell indices.** The cell centre is at
  `(i+0.5)*cell_len` and `AABB::isInside` is a closed interval plus eps. `0,0,32,32` covers cell
  indices `0..31`, matching `AABB{0, Ni*Cell_Len}` in `cavity2d`; implementing it as a closed range
  of cell indices paints one extra row and the 1000-step baseline necessarily fails.
- **`boundary` only allows the full domain.** The criterion in `SetupBoundary` is "paint if any
  neighbour is outside the box", and it does not check the cell itself; a non-full domain paints
  outside the box silently.
- **One-shot scaffold; no "config as source, regenerate every time".** A single hand edit forks it,
  and it also drags the drift between the generated code and the `src/` API into the build. The
  output belongs to the user; `--check` only guards a golden directory against drift.
- **The checksum must lock the loop phase, not only aggregate quantities.** The POP sum/sumsq is
  invariant under permutation, so it alone cannot distinguish `Stream→collide` from
  `collide→Stream`; RHO sumsq is pinned as well. CUDA follows the `collide→Stream` of `cavity2d_cu`;
  against the CPU `Stream→collide`, the boundary phase difference makes `maxabs` differ by about
  `4.9e-4`, but RHO sumsq differs by only `1.5e-7`, so only RHO is pinned across back ends.
- **The `BasicComm&` of a periodic boundary must stay stable.** `LatticeComm` holds a `BasicComm&`;
  `BasicCommSet` stores in a `std::deque`, so `emplace_back` does not invalidate references to
  existing elements. This is the only `src/` change on the casegen side, and capacity patches such
  as `reserve(16)` are no longer needed.
- **Geometry is a build product, not an opaque binary in the case directory.** `[geometry.source.*]`
  in the `.case` is the only source description; casegen renders it into reviewable `.geo` text,
  `tools/geometry_build` compiles that into `.flbgeo`, and the generated Makefile rebuilds it before
  linking. `--check` compares the `.geo`, `*.flbgeo` is ignored by `.gitignore`; the prebuilt `file`
  mode is kept, for a mask produced by an external tool or not yet covered by a recipe.
- **Source order is format semantics.** `[source.<name>]` in the `.geo` is applied in file order, and
  a later solid may override an earlier domain; the parser therefore must not sort, and must not be
  swapped for an `.ini` reader that does not preserve order.
- **The patch binding belongs to `src/`; casegen only writes the `.geo`.** The `.flbgeo` persists
  role/entity/material, boundary patches and the surface band; `src/boundary/patch_binding.h` maps a
  patch kind to task flag bits and paints them into the flag field, and `tests/feature/patch_bc`
  covers the per-entity override of two patches of the same kind. `[patch.<name>]`
  (`entity`/`kind`/`box`) in the `.geo` is written into the patch table by `tools/geometry_build`;
  casegen itself does not yet emit `[patch.*]` sections.

## 3. Rejected alternatives

- **One file writing both structure and kernel**: two syntaxes, line numbers and silent skipping,
  see §2.
- **A run-time interpreter**: per-cell dispatch requires inlining, `__any__`, and no allocation, and
  a CUDA kernel requires the same source.
- **Config as source, regenerate every time**: a hand edit by the user forks it; see §2.
- **Reusing `csegen` as the operator back end**: hard-coded basenames, submodule contract; see §2.
- **A fixed driver plus generated headers**: the cavity2d/3d diff is only about 40 lines, so
  introducing a hook system for it is not worth it.
- **Committing the `.flbgeo` fixture as the source of truth**: binary, native-endian, scalar-tagged,
  unreadable in review, and a format or scalar-type upgrade rebuilds every fixture. Instead the
  `.geo` text is committed and rebuilt by the tool.
- **casegen writing the `.flbgeo` directly**: that would copy the rasterise, band and patch semantics
  to the Python side and leave two implementations against the C++ builder. `.geo` +
  `tools/geometry_build` keeps a single implementation.

## 4. Deferred and not implemented

- `open` / `anti_pressure`: an explicit error at generation time, and not listed by `--list`. They
  need inlet/outlet regions, pressure boundaries and another main loop. Exit condition: hook them up
  once there is a real open-boundary case and a fixed baseline.
- `plain` + MPI: the generated code errors out explicitly when `mpi().getSize() > 1`, until there is
  an `mpirun` checksum test for plain. Multi-layer/multi-axis `periodic` is unverified.
- FP16: unverified. Declare it only after a manual run of the CUDA/FP16 row of `AGENTS.md` §3.
- The generation-time refusals (`thread_num != 1`, `f[9]`/2D `u[2]`, `cuda = 1` + `SMAGORINSKY`) are
  guards, not gaps: do not relax them when changing the implementation.
- **Geometry consumers**: `.vti` export, the CUDA device mirror, the `examples/` migration and the
  minimum wall-thickness check are still not done; the 2D path stays as it is. casegen does not yet
  emit `[patch.*]` sections (see §2).
- Direct STL/CSG import: casegen still only generates the `.geo`, which `tools/geometry_build`
  rasterises; direct import needs a geometry kernel/boolean operations and is not done for now.

## 5. Next directions

Priority is ordered by value / maintenance cost; any new feature needs a golden or one manual
verification record, otherwise it does not go into `--list`.

| Direction | Value | Cost | Verdict |
| --- | --- | --- | --- |
| `open` / `anti_pressure` | Medium-high: the largest missing boundary family | 1–2 weeks | Do it once there is a real open-boundary case |
| geometry patch/entity → task/BC (B4) | High: let patches really take part in the boundary | — | Landed (`src/boundary/patch_binding.h`) |
| casegen emits `[patch.*]` sections | Medium: let generated cases use B4 | 1–2 days | Do it once there is a golden |
| `.vti` export | Medium: visualise/debug the geometry mask | 1–2 days | Do it once there is a consumer |
| Regional initial value `set_value` | Medium: the most useful gap after the AABB | 3–5 days | On demand |
| Catalogue extension (force variants / non-Newtonian / thermal) | Medium: 1–2 days each plus one golden | Demand-driven | Add per case |
| CA / free surface / AMR / run-time interpreter | — | Months | Not doing |

The table gives a verdict per direction; landed items are marked in the verdict column, the rest are
not landed. Acceptance for the landed directions (CI and `make test-casegen`, the DSL math functions
and `field()`, the 3D CUDA baseline, the plain+MPI trade-off, `-O3 -march=native` in the generated
Makefile, geometry file / Phase 4) lives in the golden directories and in "Exact checksums" in
`tools/casegen/README.md`; this file is updated only when a decision changes.
