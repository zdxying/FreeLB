# AGENTS.md — working on FreeLB

FreeLB is a C++17 LBM solver (single-phase, free surface, CA solidification,
runtime AMR; CPU/OpenMP/MPI and single-GPU CUDA). One author, personal project.
Pre-2026 code is the reference for *structure and terminology*; do not modernize it.

**Before touching `src/`, read `docs/ARCHITECTURE.md`** (layering, data model,
module map, build, CUDA mirror). This file holds only what the code will not
tell you: the process (§1, §5), what fails silently, what must pass.

## 1. How to work

- Author: what to build, the acceptance criterion, the trade-offs. Agent: plan,
  code and evidence (quote real command output). Before a non-trivial edit, state
  the plan in at most five lines: which files change, which §3 row you will run.
- Change what the task needs and nothing else: no drive-by refactors, no
  reformatting, no renaming sprees. A genuinely needed abstraction is fine —
  justify it in the plan. Every line of the diff should be explainable by the
  task; if you cannot say which requirement it serves, revert it.
- When asked why something fails, find the cause before naming it; a familiar
  error can have a different cause. Separate confirmed from hypothesis, and fix
  it only if asked.
- Say what you could not verify. "Not tested" is fine; "should work" is not.
- Code, comments and commit messages in English; answer the author in the
  language of the question. Documentation is written in English.
- Done means: exercised by a real run (not just compiled), the §3 row run, new
  behaviour covered by a test or checksum, comments updated, stale docs reported.
- If this file or a doc disagrees with the code, the code wins for local style;
  report the mismatch instead of choosing silently.
- No bare `MPI_Barrier` "to be safe": it synchronises without taking part in matching.
  If a call needs one, say what it matches.
- The PopLattice path is frozen (README "Umbrella headers"): do not extend it, fix it
  for a new feature, or let an active header include it.

**Style**: match the file you edit; `.clang-format` is the fallback, not a licence
to reformat. Private members `_camelCase`, device handles `dev_`, manager drivers
a PascalCase verb. New files: license block, `// filename`, `#pragma once`,
includes from `src/`. Policy template parameter over a runtime `if` in a cell loop;
kernel-callable code is `__any__` (no allocation, no `std::vector`, no per-cell
branch on `Thread_Num`). Comments explain *why*; new modules get their own namespace.

## 2. Constraints that fail silently

**Build**
- `TARGET = name.exe` must be its own literal line; `tests/run_tests.sh` parses it.
- `-D_UNROLLFOR` is probed with `findstring` in `FLAGS` only; in `CXXFLAGS` the
  codegen silently disappears (`make.mk:84`).
- `make mpi` sets both `CXXC := mpic++` and `-DMPI_ENABLED`; never set one alone.
- CUDA arch: the `-gencode` form, why `-rdc=true` is wrong and the unchecked-launch
  trap are at the top of `examples/cavity2d_cu/Makefile`.

**MPI**
- A comm plan dies when what it indexes is rebuilt: every rebuild must reset it for
  *every* depth, else stale block indices are read (silent corruption, not a crash).
  Rebuild rule: `docs/ARCHITECTURE.md` §2.
- Send buffers must outlive `MPI_Waitall`; the plan owns them for that reason.
- `MPI_ENABLED` must pair with `mpic++`, not `g++`. `MPI_RANK(x)` expands to an
  early `return` (void only); use `IF_MPI_RANK` in `main` (`utils/mpi_macros.h:46`).

**CUDA**
- `copyToDevice()` does not push container rotation state: after `normalize()`
  (checkpoint load, swap, remap) call `refreshDeviceObject()` (`io/checkpoint.hh` `loadOneBlock()`).
- `check_cuda_status()` is called around memcpy/malloc in the existing code; a launch
  path needs `cudaGetLastError()` at its end, a timed one `cudaDeviceSynchronize()` too.

**IO**
- `;` starts a comment only at line start (`io/ini_reader.h` `iniReader`); an inline `;` is
  not stripped and becomes part of the value. New keys use `getValueOr` so older `.ini` files still run.
- Write a checkpoint only at a quiescent point (comms finished, before the next
  `Stream()`); a restart assumes the same parameters and block set.

## 3. Tests and gates

`tests/<category>/<name>/` is built and run by `tests/run_tests.sh`; exit status
is the verdict. Categories: `block`, `cuda`, `feature`, `io`, `regression`.

| Change | Required |
| --- | --- |
| `docs/`, `README.md` | no regression run; make sure commands and links resolve |
| one example, test, or local fix | that target builds and runs |
| `src/` (not CUDA/MPI/IO) | `make test-cpu`, or `make test` with a GPU; quote the summary line, never a fixed count |
| `src/parallel`, MPI paths | above + `block/plain_comm_mpi`, `block/refblock`, `block/amr_comm`, `block/comm_index`, `block/checkpoint_decomp`, `feature/hidden_comm`, `block/normal_comm_correctness` |
| numerics, boundaries, equilibrium | above + `regression/cavity2d`, checksum change explained |
| `InitComm()` / `SharedComm` cell lists | above + `block/normal_comm_correctness --strict`; red by construction — record the count; defect: `docs/INTP_COMM_OVERFLOW_BUG.md` |
| CUDA, FP16, containers | `make test` all green (`--no-cuda` without a GPU) + the manual FP16 run |
| checkpoint / restart | `feature/checkpoint`, `block/checkpoint_decomp`, `cuda/checkpoint` |
| headers, includes | `tools/check_headers.sh` |
| `make.mk`, xcore (optional pool allocator), `tools/cse` | root `make`, then an example with and without xcore |

Missing tool? Run what you can, say which row was not met; never downgrade a gate
silently.

## 4. Changing numerics

- An intentional numeric change says so, updates the `tests/regression/cavity2d`
  baseline in the same commit, and justifies it.
- Iteration order and floating-point association are part of the result;
  reordering a sum in a hot loop is a numeric change, not a cleanup.
- Unit conversion stays in `lbm/unit_converter.h`. Never compare a `double`
  against a `float`-typed `T`.
- A scheme arithmetic change needs a per-configuration checksum against the
  pre-change tree, not one aggregate.
- Performance claims need the machine, the command and MLUPS (million lattice
  updates per second) before/after; neutral is a valid result.

## 5. Docs and commits

- **Architecture, not implementation**: the invariant, the boundary, the reason. Mechanism
  lives in the code. A paragraph that needs the source
  to parse is implementation detail — keep the pointer.
- **Docs are for humans**: readability is the first requirement; being parseable by an agent
  is not.
- **References go one way**: docs point to code (file, symbol, test); code comments never
  point to `docs/`. A comment carries its own reason.
- `docs/ARCHITECTURE.md` — the present: layering, data model, module map, build, CUDA
  mirror. It must stay faithful to the architecture.
- `docs/<TOPIC>_PLAN.md` — a decision: status/date/target on line 3, the choice, its
  evidence, an explicit rejected-and-why section. A defect or todo lives in the doc that
  owns the topic; there is no global TODO.
- **A code change does not edit `docs/`.** Report the stale claim; the main agent decides
  whether to land it and a doc-scoped agent writes it. Code comments are part of the code.
- A doc states the present, not its history — git carries the changelog. New material names
  its exit condition in its first sentence: the commit, test or measurement that will make
  it obsolete. Editing a doc is when you prune; past ~400 lines, prune first.
- README is user-facing, not a changelog; fix a claim your change invalidates in the same
  commit. No new `*.md` for a one-off task.

Commit: `[type] area: imperative summary` (lowercase, no trailing period), then
bullets of why/what/how, then `Verified: <command> (<counts>)`. `type` follows the
log (feat, fix, perf, refactor, test, bench, docs, build, tools); English, no AI
trailer, subject under 80 characters. Longer reasoning goes in the linked PLAN doc.

**This file's own rule**: a line earns its place only if ignoring it produces a
silent wrong result or an unmet gate; mechanism and narrative belong in the doc
it points to. Keep under 150 lines — when it grows, delete before adding.

## 6. Ask the author first

1. A new dependency, submodule or vendored code.
2. Deleting or renaming a public header, an example or a test.
3. Any `git commit` / `push` / branch / history operation.
4. Changing a default or a numerical result.
5. An on-disk format change that breaks existing files.
6. Weakening a test or skipping a gate row.

Otherwise proceed and report. Two defensible designs that differ in performance
or public API: with the author present, present the trade-off with measurements;
in an unattended run, take the one that keeps current numerics and public API,
record the alternative in one line, and continue.
