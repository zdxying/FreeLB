# casegen: generate a FreeLB case from a spec

`casegen` assembles a buildable FreeLB case from two text files:

* a `.case` file: mesh, flags, geometry, fields, tasks, schedule, checksums;
* optional `.op` files: a per-cell collision operator in a small DSL.

The output is a **one-shot scaffold**: `<name>.cpp` (or `<name>.cu`), `<name>.ini`,
`Makefile`, and one header per custom operator. The generated files belong to
you; `casegen` only rewrites them when you ask it to.

The generator never edits `src/`, `make.mk` or `tests/run_tests.sh`. It is a
pure Python 3 tool with no third-party dependencies. Design decisions and the
rejected alternatives are in [`docs/CASEGEN_PLAN.md`](../../docs/CASEGEN_PLAN.md).

## Status

Supported today:

* 2D `D2Q9`, 3D `D3Q19` / `D3Q27`;
* schedules `hidden` and `plain` (the latter with `periodic = x|y|z`);
* catalogue operators `bgk`, `bgk_force`, `mrt` (D2Q9 only), `les`;
* boundaries `normal`, `movingwall`;
* fields `RHO`, `VELOCITY`, `POP`, `FORCE`, `SMAGORINSKY`;
* writers `rho`, `velocity`, `flag`, checkpoint/restart, checksum assertions;
* custom `.op` collision operators, shared by the host and CUDA paths;
* MPI output (`mpi = 1`) with `decompose = auto` and inner-cell checksums;
* CUDA output (`cuda = 1`) for `hidden`, single-block cases.

Not supported: `open`/`anti_pressure`, CA solidification, free surface, AMR,
STL/CSG geometry.

## Invocation

```sh
# generate into a new directory (ROOT is computed from the depth)
python3 tools/casegen/casegen.py tools/casegen/cases/cavity2d.case -o examples/cavity2d_casegen

# compare a golden directory with its spec, exit 2 on drift
# (pass --selftest so the regen Makefile is part of the comparison)
python3 tools/casegen/casegen.py tools/casegen/cases/cavity2d.case -o tests/feature/casegen --check --selftest

# overwrite files whose content differs
python3 tools/casegen/casegen.py ... -o DIR --force

# emit a `regen: --check` rule into the Makefile (used by the golden tests)
python3 tools/casegen/casegen.py ... -o DIR --selftest

# print the catalogue
python3 tools/casegen/casegen.py --list
```

Exit codes: `0` success, `1` spec error or refused overwrite, `2` `--check`
found a missing or different file. `--check` renders in memory using the target
directory's `ROOT`, so it never writes a temporary tree.

`ROOT` is `../..` for `examples/<name>` and `../../..` for
`tests/<category>/<name>`. `TARGET` is always a literal `<name>.exe` line: the
test runner parses it without expanding variables.

Run the generated binary from its own directory: `<name>.ini`, `vtkoutput/` and
`ckpt/` are resolved against the current working directory. `[operators]` paths
are relative to the spec file, not to the output directory.

## `.case` format

Ordered INI dialect:

* `[section]` and `[section.name]`; region sections execute in file order;
* `key = value`, with comma-separated lists;
* a leading `#` or `;` starts a comment; inline comments are rejected;
* unknown/duplicate sections and keys are errors, so a typo never becomes a
  silently ignored option.

Minimal 2D example:

```ini
[case]
name = cavity2d
dim = 2
float = double
latset = D2Q9

[fields]
list = RHO, VELOCITY, POP

[flags]
void = 1
bulk = 2
wall = 4
moving_wall = 8

[grid]
n = 32, 32
cell_len = 1
block_cell_len = 10
thread_num = 1
decompose = single

[physics]
converter = simplified
char_u = 0.1
rt = 0.56
init_u = 0, 0

[region.cavity]
shape = aabb
box = 0, 0, 32, 32
action = paint
flag = bulk

[region.walls]
shape = aabb
box = 0, 0, 32, 32
action = boundary
flag = wall

[region.lid]
shape = aabb
box = 1, 31, 31, 32
action = paint
flag = moving_wall
only_if = wall

[region.lid_velocity]
shape = aabb
box = 1, 31, 31, 32
action = set_velocity
flag = moving_wall
value = 0.1, 0

[tasks]
bulk = bgk
wall = normal
moving_wall = movingwall

[schedule]
type = hidden
steps = 1000
output_step = 1000
tol = 0

[output]
writers = rho, velocity

[checksum]
enable = 1
expect_pop = 7728.702031
expect_rho = 1155.971535
tol = 1e-6
```

### `[case]`

| Key | Required | Value |
| --- | --- | --- |
| `name` | yes | `[A-Za-z_][A-Za-z0-9_]*`; file stem and `TARGET` |
| `dim` | yes | `2` or `3`; must match `latset` |
| `float` | yes | `double` or `float` |
| `latset` | yes | `D2Q9`, `D3Q19`, `D3Q27` |
| `cuda` | no | emit `<name>.cu` and an nvcc Makefile |
| `mpi` | no | emit an mpic++ Makefile and `-DMPI_ENABLED`; requires `decompose = auto` and cannot be combined with `cuda` |
| `selftest` | no | same as `--selftest` |

### `[fields]`, `[flags]`, `[grid]`

`[fields] list` is a non-empty, duplicate-free list of `RHO`, `VELOCITY`, `POP`,
`FORCE`, `SMAGORINSKY`. FLAG is implicit: every case has a flag manager.

`[flags]` maps names to distinct powers of two in `1..128`; `void` is mandatory
and is the initial value of every cell (it normally remains only on ghost/overlap
cells).

`[grid]`: `n = Ni, Nj[, Nk]`, `cell_len`, `block_cell_len`, `thread_num`,
`decompose = single | auto`. `thread_num` must be `1`: the OpenMP per-cell paths
are not wired, so a larger value is rejected instead of silently ignored. `auto`
uses `BlockGeometryHelper2D/3D` with a rank-independent `CreateBlocks` plus
`LoadBalancing`; `mpi = 1` requires `auto`, while `cuda = 1` requires `single`.

### `[physics]`

`converter = simplified` uses `char_u`, `rt`, `rho_init` (default 1), `init_u`.
`converter = full` uses `cell_len` (mm), `rho_ref`, `kine_visc`, `u_max`, `rt`.
Optional `force = fx, fy[, fz]` fills the `FORCE` field; `smagorinsky = value`
fills `SMAGORINSKY`. The field list must contain the corresponding field.

### Geometry: `[region.<name>]`

Only `shape = aabb`. `box` is an **AABB endpoint** list in cell units
(`x0,y0,x1,y1`, or `x0,y0,z0,x1,y1,z1`), not a closed index range. Cell `i` has
its centre at `(i+0.5)*cell_len`, so `0,0,32,32` covers indices `0..31` — the
same region as `AABB{0, Ni*Cell_Len}` in the examples.

Actions execute in section order:

* `paint`: set every cell in the box to `flag`; optional `only_if = <flag>`
  repaints only cells that already carry that flag.
* `boundary`: `SetupBoundary<LatSet>(box, flag)`. It paints every cell with a
  neighbour outside the box, **without checking that the cell itself is inside**;
  therefore `box` must be the full domain, and a `paint` must precede it. The
  generator rejects anything else.
* `set_velocity`: write `VELOCITY` on cells that carry `flag`; `value` is a
  physical velocity, converted with `BaseConv.getLatticeU`.

The standard cavity order is `paint` bulk → `boundary` walls → `paint(only_if)`
the lid.

### `[tasks]`

`<flag or flag|flag> = <operator>`. A task key may be a bitwise OR of flags; the
emitter sorts keys with more bits first so a composite key is tested before its
parts. Each flag belongs to exactly one task. For `hidden`, every painted flag
except `void` needs a task, otherwise `TaskSelector` silently skips it.

### `[schedule]`

| Key | Value |
| --- | --- |
| `type` | `hidden` or `plain` |
| `steps`, `output_step`, `tol` | required / required / optional (default 0) |
| `checkpoint_step`, `restart`, `restart_dir` | optional |
| `periodic` | optional `x`, `y` or `z`; `plain` only, serial, single block |

`hidden` applies boundary schemes as collision tasks:
`Stream(); ApplyCellDynamicsWithHiddenComm<NSTask>(FlagFM);`.
`plain` applies them through `BlockBoundaryManager`:
`RhoU → collision → (periodic) → Stream → BM.Apply(step) → Communicate`.
`periodic` requires a flag named `periodic` and a multi-cell strip; it uses
`FixedPeriodicBoundaryManager` and is rejected with MPI/CUDA.

### `[output]`, `[checksum]`

`writers = rho, velocity, flag`; `flag_writer = 0|1`. `[checksum]` enables a
post-run assertion on `POP` and/or `RHO` sumsq; a mismatch prints the actual
values and returns 1. With `[case] mpi = 1` the checksum reads only inner cells
and reduces sum/sumsq/max across ranks; pin `expect_pop_mpi`/`expect_rho_mpi`
for that build, while the non-MPI build keeps the all-cell `expect_*` path. POP
is a `CyclicArray`: its sum/sumsq are permutation-invariant, so the pinned value
is stable across the CPU/CUDA rotation, but the accumulation order differs by a
few ulps. Pin `expect_rho` for cross-backend comparisons.

## `.op` operator DSL

```ini
[operators]
my_bgk = ../operators/my_bgk.op   # key must equal `operator NAME` in the file

[tasks]
bulk = my_bgk
```

```text
operator my_bgk {
  var rho: scalar
  var u: vector
  var feq: array
  body {
    rho = moment0(f)
    u = moment1(f, rho)
    for i in q {
      feq[i] = feq2(rho, u, i)
      f[i] = omega * feq[i] + (1 - omega) * f[i]
    }
  }
}
```

Statements: `var NAME: scalar|vector|array`, assignment, `NAME[i] = expr`,
`for i in q { ... }`, leading `#` comments. Expressions: numbers, `+ - * /`,
parentheses, unary minus, indexing, calls. Variables are zero-initialised.

Scalar functions: `sqrt`, `abs`, `min`, `max`, `clamp(x, lo, hi)`; they emit
`std::sqrt/abs/min/max/clamp` and take scalar arguments. A vector variable can
be read by component with an integer literal (`u[0]`, `u[1]`, `u[2]`); components
are read-only. `field(RHO|VELOCITY|FORCE|SMAGORINSKY)` reads the corresponding
cell field (scalar or vector) and adds it to the operator's field closure;
`force()` remains an alias for `field(FORCE)`.

Builtins expand to the same arithmetic as `src/lbm`:

| Builtin | Expansion | Needs |
| --- | --- | --- |
| `moment0(f)` | `sum_i cell[i]` | POP |
| `moment1(f, rho)` | `sum_i cell[i]*c_i`, then `/= rho` | POP |
| `moment1_force(f, F)` | same, then `+= F*0.5`, then `/= rho` | POP, FORCE |
| `feq2(rho, u, i)` | `Equilibrium::Order2` (`u2` once) | POP |
| `w(i)`, `c(i)` | `latset::w/c<LatSet>(i)` | — |
| `dot(a,b)`, `norm2(a)` | `a*b`, `a.getnorm2()` | — |
| `omega`, `fomega` | `cell.getOmega()`, `cell.getfOmega()` | — |
| `force()` | `cell.get<FORCE<T, LatSet::d>>()` | FORCE |
| `field(NAME)` | `cell.get<NAME<T>>()`; VELOCITY/FORCE take `<T, d>` | NAME |
| `force_term(u, F)` | fills `std::array<T,q>` via `force::Force<CELL>::apply` | FORCE |
| `sqrt/abs/min/max/clamp` | `std::sqrt/abs/min/max/clamp` | — |
| `u[0..2]` | `(u)[k]`, read-only | — |

The emitted header is `template <typename CELL> struct PascalName { __any__
static void apply(CELL&); }`, so the same source compiles for host and device.
`c(i)` returns a vector; combine it with `dot`/`norm2`, or read a component.
Errors carry `path:line`. `my_bgk` becomes `MyBgk`.

## Exact checksums

The golden directories pin these values (relative 1e-6 unless noted):

| Case | POP sumsq | RHO sumsq |
| --- | --- | --- |
| `tests/feature/casegen` (2D BGK) | 7728.702031 | 1155.971535 |
| `tests/feature/casegen_op` (custom-only `.op` BGK) | 7728.702031 | 1155.971535 |
| `tests/feature/casegen_dsl` (2D `.op` BGK) | 7728.702031 | 1155.971535 |
| `tests/feature/casegen_3d` (16³ D3Q19) | 1327.451851 | 5831.917013 |
| `tests/feature/casegen_force` | 586.5301684 | 1155.983056 |
| `tests/feature/casegen_mrt` | 586.534128 | 1155.984267 |
| `tests/feature/casegen_les` | 586.531249 | 1155.988022 |
| `tests/feature/casegen_plain` | 287.5573693 | 1148.306514 |
| `tests/feature/casegen_mpi` (inner cells, MPI, 1e-12) | 445.51566903028549 | 1023.9878563745673 |
| `tests/feature/casegen_dsl_math` (2D `.op`, 200 steps) | 586.5034095 | 1155.987548 |
| `tests/feature/casegen_cu` (float, 1e-5, RHO only) | — | 1155.992025 |
| `tests/feature/casegen_3d_cu` (float, 1e-5, RHO only) | — | 5831.932013 |

The 3D and 200-step values are pinned to their exact grid/steps; changing either
requires re-running and updating the spec. The CUDA goldens pin the CPU(float)
reference; the GPU run may differ by a few `1e-7`.

## Tests

```sh
python3 -m unittest discover -s tools/casegen/tests   # parser, CLI, DSL
make test-casegen                                    # unit tests + every golden
make test-cpu                                        # the repository suite
```

`make test-casegen` runs `tools/casegen/check.sh`: it auto-discovers every
`tests/feature/casegen*` directory, requires a `--check --selftest` regen rule,
builds and runs each golden, and asserts its checksum. CUDA directories are
skipped when nvcc or a visible GPU is missing; MPI directories run under
`mpirun -np 2` and are skipped when `mpic++` or `mpirun` is missing. A spec
change that is not mirrored in the committed output fails the build.

## Limits

* No CA solidification, free surface, AMR, STL/CSG geometry.
* No `open`/`anti_pressure`; `periodic` is serial, single-block, one axis.
* `plain` is single-rank: the generated `main` rejects `mpi().getSize() > 1`
  until plain+MPI has a test.
* `[grid] thread_num` must be 1: the OpenMP per-cell paths are not wired, so a
  larger value is a generation error instead of a silent no-op.
* `mpi = 1` requires `decompose = auto`; it selects `mpic++` and cannot be
  combined with `cuda`. Its checksum uses inner cells and a cross-rank reduction,
  so `expect_pop_mpi`/`expect_rho_mpi` are the values to pin.
* AABB geometry only; `boundary` boxes must cover the full domain.
* Initial density is 1; initial velocity is uniform per case (plus
  `set_velocity` wall values).
* Vector components and array indices are checked against the case's `dim`/`q`
  when the operator is attached; `u[2]` in a 2D case or `f[9]` in D2Q9 is a
  generation-time error.
* `--check` compares the files it would render; it does not report extra files
  left in the output directory (for example an old `operators/*.hh`).
* CUDA is `hidden` + `single` only. It cannot carry `SMAGORINSKY`, so the `les`
  collision and `field(SMAGORINSKY)` are rejected with `cuda = 1`. The CUDA loop
  is collide→stream (matching `examples/cavity2d_cu`), while the CPU hidden loop
  is stream→collide, so the boundary phase differs by a few ulps in `maxabs`.
