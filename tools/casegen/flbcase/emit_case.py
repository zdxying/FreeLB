"""Render a CaseSpec into a buildable FreeLB case scaffold.

``render_case`` is pure: it returns ``{relative path: text}`` and never touches
the filesystem.  The caller writes the files and owns ``--check`` / ``--force``.

Only the Phase 1-3 subset is emitted: 2D D2Q9 and 3D D3Q19/D3Q27,
RHO/VELOCITY/POP/FORCE/SMAGORINSKY, the catalog BGK / BGK-force / MRT (D2Q9) /
Smagorinsky-LES collisions, normal / moving-wall bounce-back, the hidden and
plain schedules (plain with an optional serial periodic axis; it accepts a
single MPI rank only, because the plain ordering has not been validated under
MPI), single or auto
block decomposition, rho/velocity writers, checkpoint/restart and the
pop/rho checksums.  ``[case] mpi = 1`` builds with mpic++ and
``-DMPI_ENABLED``; its checksums read each block's interior and reduce across
ranks with ``MPI_SUM``/``MPI_MAX`` so overlap cells are not counted twice.

The hidden layout mirrors tests/regression/cavity2d: SimplifiedConverterFromRT
(Ni, char_u, rt), zero initial velocity, ``Stream()`` then
``ApplyCellDynamicsWithHiddenComm``, and a final ``NormalFullCommunicate()`` +
``ApplyCellDynamics<RhoUTask>`` before the checksums.  The plain schedule
follows examples/openboundary2d/openbd2d.cpp: RhoU, collision, ``Stream``,
``BlockBoundaryManager::Apply`` then ``Communicate``.  With ``periodic = x|y|z``
it also adds a ``FixedPeriodicBoundaryManager`` whose ``Apply`` runs between the
collision and the stream, exactly as examples/couette2d does.  Changing these
orders is a numerical change, not a cleanup.  The 3D layout is the same one
examples/cavity3d uses; the only dimension-dependent pieces are the vector/AABB
arity and the geometry constructor.

``[case] cuda = 1`` emits ``<name>.cu`` instead: the same host scaffold (geometry,
flags, tasks, writers, restart) plus the device mirrors, a per-step
``CuDevApplyCellDynamics`` + ``CuDevStream`` pair and a synchronous
``copyToHost`` before every host-side consumer.  The device entry points drive
block 0 only, so CUDA cases have to be single-block.
"""

from __future__ import annotations

import math
import os
import re
from typing import Dict, List, Sequence, Tuple

try:  # imported as flbcase.emit_case
    from .model import (
        BOUNDARIES,
        FIELD_TYPES,
        LATTICE_INFO,
        OPERATORS,
        CaseSpec,
        CasegenError,
        Region,
        Task,
    )
except ImportError:  # pragma: no cover - imported as a top-level module
    from model import (  # type: ignore
        BOUNDARIES,
        FIELD_TYPES,
        LATTICE_INFO,
        OPERATORS,
        CaseSpec,
        CasegenError,
        Region,
        Task,
    )

try:  # imported as flbcase.emit_case
    from .spec import PHASE
except ImportError:  # pragma: no cover - imported as a top-level module
    from spec import PHASE  # type: ignore


# Hidden communication applies boundaries as collision tasks.  The type in
# BOUNDARIES is bounceback::normal<CELL>, a per-direction boundary method for
# the plain schedule; collision.h packages the same full-way bounce-back as
# collision::BounceBack / BounceBackMovingWall.  The two forms are not
# interchangeable, which is why the task spelling is kept here.
_BOUNDARY_TASKS = {
    "normal": "collision::BounceBack<CELL>",
    "movingwall": "collision::BounceBackMovingWall<CELL>",
}


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------


def _num(value, where: str) -> str:
    """Shortest C++ literal that round-trips the value as a double."""
    if isinstance(value, bool):
        raise CasegenError(f"{where}: expected a number, got a boolean")
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise CasegenError(f"{where}: expected a number, got {value!r}")
    if not math.isfinite(number):
        raise CasegenError(f"{where}: expected a finite number, got {value!r}")
    return repr(number)


def _int(value, where: str) -> int:
    if isinstance(value, bool):
        raise CasegenError(f"{where}: expected an integer, got a boolean")
    try:
        number = int(value)
    except (TypeError, ValueError):
        raise CasegenError(f"{where}: expected an integer, got {value!r}")
    return number


def _block_count(n, block_cell_len) -> int:
    """Number of blocks along one axis for the rank-independent auto layout."""
    return max(1, int(n) // int(block_cell_len))


def _truthy(value) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def _force_expr(spec: CaseSpec) -> str:
    """FORCE initial value: the lattice-unit constant from [physics] force."""
    force = spec.physics.get("force")
    if force is None:
        return "Vector<T, LatSet::d>{}"
    comps = ", ".join("T(%s)" % _num(v, "[physics] force") for v in force)
    return "Vector<T, LatSet::d>{%s}" % comps


def _smagorinsky_expr(spec: CaseSpec) -> str:
    value = spec.physics.get("smagorinsky")
    if value is None:
        return "T{}"
    return "T(%s)" % _num(value, "[physics] smagorinsky")


def _pascal(name: str) -> str:
    parts = [p for p in re.split(r"[^0-9A-Za-z]+", name) if p]
    ident = "".join(p[:1].upper() + p[1:] for p in parts)
    if not ident:
        raise CasegenError(f"cannot derive a C++ identifier from {name!r}")
    if ident[0].isdigit():
        ident = "F" + ident
    return ident


def _require(mapping: Dict, key: str, where: str):
    if key not in mapping or mapping[key] is None:
        raise CasegenError(f"{where}: missing required key '{key}'")
    return mapping[key]


def _as_list(value) -> List:
    if value is None:
        return []
    if isinstance(value, str):
        return [part.strip() for part in value.split(",") if part.strip()]
    return list(value)


def _flag_ids(spec: CaseSpec) -> Dict[str, str]:
    ids: Dict[str, str] = {}
    used = set()
    for name in spec.flags:
        ident = _pascal(name) + "Flag"
        if ident in used:
            ident = f"{ident}_{len(ids)}"
        used.add(ident)
        ids[name] = ident
    return ids


def _region_id(region: Region, index: int) -> str:
    return f"aabb_{index}_{_pascal(region.name)}"


def _key_expr(flags: Sequence[str], flag_ids: Dict[str, str]) -> str:
    exprs = []
    for name in flags:
        if name not in flag_ids:
            raise CasegenError(f"task references unknown flag '{name}'")
        exprs.append(flag_ids[name])
    if len(exprs) == 1:
        return exprs[0]
    return "std::uint8_t(" + " | ".join(exprs) + ")"


# ---------------------------------------------------------------------------
# validation of the emitted subset
# ---------------------------------------------------------------------------


def _check_supported(spec: CaseSpec) -> None:
    if not spec.name or os.path.basename(spec.name) != spec.name or spec.name in {".", ".."}:
        raise CasegenError(f"case name {spec.name!r} is not a safe file stem")

    dim, _q = LATTICE_INFO[spec.latset]
    if spec.dim != dim:
        raise CasegenError(
            f"latset {spec.latset} needs dim = {dim}, but [case] dim = {spec.dim}"
        )
    if len(spec.n) != dim:
        raise CasegenError(
            f"[grid] n must have {dim} entries for a {dim}D case, got {spec.n}"
        )
    if any(float(v) != int(v) or int(v) <= 0 for v in spec.n):
        raise CasegenError(f"[grid] n must be {dim} positive integers, got {spec.n}")
    if float(spec.cell_len) <= 0.0:
        raise CasegenError(f"[grid] cell_len must be positive, got {spec.cell_len}")
    if int(spec.thread_num) < 1:
        raise CasegenError(f"[grid] thread_num must be at least 1, got {spec.thread_num}")
    if int(spec.block_cell_len) < 1:
        raise CasegenError(
            f"[grid] block_cell_len must be at least 1, got {spec.block_cell_len}"
        )

    schedule = spec.schedule.get("type", "hidden")
    if schedule == "open":
        raise CasegenError(
            "[schedule] type = 'open' is not implemented yet: it needs the "
            "anti_pressure boundary, which this generator does not emit"
        )
    if schedule not in {"hidden", "plain"}:
        raise CasegenError(
            f"[schedule] type = {schedule!r}: expected hidden or plain"
        )
    if spec.decompose not in {"single", "auto"}:
        raise CasegenError(
            f"[grid] decompose = {spec.decompose!r}: expected 'single' or 'auto'"
        )

    periodic_axis = spec.schedule.get("periodic")
    if periodic_axis is not None:
        if schedule != "plain":
            raise CasegenError(
                "[schedule] 'periodic' needs type = plain (it is combined with a "
                "BlockBoundaryManager)"
            )
        if periodic_axis not in "xyz"[:dim]:
            raise CasegenError(
                f"[schedule] periodic axis {periodic_axis!r} is not in 'xyz'[:{dim}]"
            )
        if int(spec.thread_num) != 1:
            raise CasegenError(
                "[schedule] periodic is serial-only: set thread_num = 1"
            )
        if spec.decompose != "single":
            raise CasegenError(
                "[schedule] periodic needs decompose = single"
            )
        if "periodic" not in spec.flags:
            raise CasegenError(
                "[schedule] periodic needs a flag named 'periodic'"
            )
        axis = "xyz".index(periodic_axis)
        for other in range(dim):
            if other != axis and int(spec.n[other]) < 3:
                raise CasegenError(
                    f"[schedule] periodic axis {periodic_axis}: grid axis {other} "
                    f"has {spec.n[other]} cells, the face strips need at least 3"
                )

    seen_bits = set()
    for name, bit in spec.flags.items():
        bit = int(bit)
        if not 0 < bit < 256 or bit & (bit - 1):
            raise CasegenError(
                f"[flags] {name} = {bit}: flags must be distinct powers of two "
                "in uint8_t range"
            )
        if bit in seen_bits:
            raise CasegenError(f"[flags] {name} = {bit}: duplicate flag bit")
        seen_bits.add(bit)
    if spec.void_name is None and 1 in seen_bits:
        raise CasegenError(
            "[flags] has no 'void' flag and bit 1 is already used; unpainted "
            "cells would collide"
        )

    unknown = [f for f in spec.fields if f not in FIELD_TYPES]
    if unknown:
        raise CasegenError(
            "unknown field(s) " + ", ".join(unknown)
            + "; supported: " + ", ".join(FIELD_TYPES)
        )
    # The rho/U task writes both fields; every schedule emits it.
    for required in ("RHO", "VELOCITY"):
        if required not in spec.fields:
            raise CasegenError(
                f"the rho/U pass needs {required} in [fields]"
            )
    if not spec.tasks:
        raise CasegenError("[tasks] must define at least one task")

    converter = spec.physics.get("converter", "simplified")
    if converter not in {"simplified", "full"}:
        raise CasegenError(
            f"[physics] converter = {converter!r}: expected 'simplified' or 'full'"
        )

    if spec.cuda and schedule != "hidden":
        raise CasegenError(
            f"[case] cuda = 1 with [schedule] type = {schedule!r}: the device "
            "path is implemented for the hidden schedule only"
        )

    if spec.cuda and spec.decompose != "single":
        # The device entry points in block_lattice.hh drive BlockLats[0] only,
        # so an auto-decomposed case would silently simulate one block.
        raise CasegenError(
            f"[grid] decompose = {spec.decompose!r} with [case] cuda = 1: the "
            "device path runs a single block, so only decompose = single is "
            "supported"
        )

    if spec.geometry is not None:
        if spec.dim != 3:
            raise CasegenError(
                f"[geometry] is only supported for dim = 3, got dim = {spec.dim}"
            )
        for key in ("domain_flag", "void_flag"):
            name = getattr(spec.geometry, key)
            if name not in spec.flags:
                raise CasegenError(
                    f"[geometry] '{key}' names undefined flag '{name}'"
                )

    if spec.cuda and "SMAGORINSKY" in spec.fields:
        # cudev mirrors the fields through cudev_FieldType; SMAGORINSKY is a
        # Data<T, ...> without one, so les and field(SMAGORINSKY) cannot be
        # instantiated on the device.
        raise CasegenError(
            "[case] cuda = 1 with SMAGORINSKY in [fields]: the CUDA device "
            "field set has no SMAGORINSKY mirror, so the les collision and "
            "field(SMAGORINSKY) cannot be built for the device"
        )

    # BaseConverter::getLatRhoInit() is 1 by construction; a different lattice
    # initial density would need a converter change, not a spec key.
    rho_init = float(spec.physics.get("rho_init", 1))
    if rho_init != 1.0:
        raise CasegenError(
            f"[physics] rho_init = {rho_init!r}: the lattice initial density is 1"
        )

    if converter == "full" and "cell_len" in spec.physics:
        if float(spec.physics["cell_len"]) != float(spec.cell_len):
            raise CasegenError(
                "[physics] cell_len is the grid cell size; it must match "
                "[grid] cell_len"
            )

    _check_regions(spec)

    for task in spec.tasks:
        _ = _task_cpp(spec, task)
        needs = _task_needs(spec, task)
        missing = needs - set(spec.fields)
        if missing:
            raise CasegenError(
                f"task '{task.op}' needs field(s) {', '.join(sorted(missing))}, "
                f"not in [fields] list"
            )

    checksum = spec.checksum
    if _truthy(checksum.get("enable", 0)):
        expect_pop = checksum.get("expect_pop")
        expect_rho = checksum.get("expect_rho")
        expect_pop_mpi = checksum.get("expect_pop_mpi")
        expect_rho_mpi = checksum.get("expect_rho_mpi")
        if expect_pop is not None and "POP" not in spec.fields:
            raise CasegenError("[checksum] expect_pop requires POP in [fields]")
        if expect_rho is not None and "RHO" not in spec.fields:
            raise CasegenError("[checksum] expect_rho requires RHO in [fields]")
        if expect_pop_mpi is not None and "POP" not in spec.fields:
            raise CasegenError("[checksum] expect_pop_mpi requires POP in [fields]")
        if expect_rho_mpi is not None and "RHO" not in spec.fields:
            raise CasegenError("[checksum] expect_rho_mpi requires RHO in [fields]")
        if (
            expect_pop is None
            and expect_rho is None
            and expect_pop_mpi is None
            and expect_rho_mpi is None
        ):
            raise CasegenError(
                "[checksum] enable = 1 needs expect_pop/expect_rho or "
                "expect_pop_mpi/expect_rho_mpi; there is nothing to assert"
            )
        tol = checksum.get("tol", 1e-6)
        if not math.isfinite(float(tol)) or float(tol) < 0.0:
            raise CasegenError(f"[checksum] tol = {tol!r} must be finite and >= 0")

    if getattr(spec, "mpi", False) and spec.cuda:
        raise CasegenError(
            "[case] mpi = 1 and cuda = 1 select different compilers "
            "(mpic++ / nvcc); generate a CPU case for MPI"
        )
    if getattr(spec, "mpi", False) and spec.decompose != "auto":
        raise CasegenError(
            "[case] mpi = 1 needs [grid] decompose = auto; a single block would "
            "be built on every rank"
        )


def _check_regions(spec: CaseSpec) -> None:
    """Reject geometry that would compile but silently misbehave at run time."""
    dim = spec.dim
    full_domain = tuple(0.0 for _ in range(dim)) + tuple(float(v) for v in spec.n)
    full_text = ", ".join(["0"] * dim + [str(v) for v in spec.n])
    void_name = spec.void_name
    for region in spec.regions:
        if region.action != "boundary":
            continue
        # SetupBoundary paints every cell that has a neighbour outside the box,
        # including cells outside the box itself; only the full domain is safe.
        box = tuple(float(v) for v in region.box)
        if box != full_domain:
            raise CasegenError(
                f"region '{region.name}': boundary boxes must cover the full "
                f"domain ({full_text}), got {region.box}"
            )

    for region in spec.regions:
        if region.flag is None:
            continue
        covering = [t for t in spec.tasks if region.flag in t.flags]
        if region.action == "paint":
            if region.flag != void_name and not covering:
                raise CasegenError(
                    f"region '{region.name}': flag '{region.flag}' has no task; "
                    "TaskSelector would silently skip those cells"
                )
        elif region.action == "boundary":
            if not any(t.op in BOUNDARIES for t in covering):
                raise CasegenError(
                    f"region '{region.name}': flag '{region.flag}' has no boundary "
                    "task"
                )
        elif region.action == "set_velocity":
            if not any(
                t.op in BOUNDARIES and BOUNDARIES[t.op]["reads_velocity"]
                for t in covering
            ):
                raise CasegenError(
                    f"region '{region.name}': flag '{region.flag}' has no boundary "
                    "task that reads VELOCITY"
                )


def _manager_schedule(spec: CaseSpec) -> bool:
    """True for schedules that run bounce-back through a BlockBoundaryManager."""
    return spec.schedule.get("type") == "plain"


def _task_needs(spec: CaseSpec, task: Task) -> set:
    if task.op in spec.operator_specs:
        return set(spec.operator_specs[task.op].needs)
    if task.op in OPERATORS:
        return set(OPERATORS[task.op]["needs"])
    if task.op in BOUNDARIES:
        return set(BOUNDARIES[task.op]["needs"])
    raise CasegenError(f"unknown operator or boundary '{task.op}'")


def _task_cpp(spec: CaseSpec, task: Task) -> str:
    if task.op in spec.operator_specs:
        return "%s<CELL>" % _pascal(spec.operator_specs[task.op].name)
    if task.op in OPERATORS:
        entry = OPERATORS[task.op]
        if entry["phase"] > PHASE:
            raise CasegenError(
                f"operator '{task.op}' is not available in phase {PHASE}"
            )
        if spec.dim not in entry.get("dims", {spec.dim}):
            raise CasegenError(
                f"operator '{task.op}' does not support dim = {spec.dim}"
            )
        return entry["cpp"]
    if task.op in BOUNDARIES:
        entry = BOUNDARIES[task.op]
        if entry["phase"] > PHASE:
            raise CasegenError(
                f"boundary '{task.op}' is not implemented in phase {PHASE}"
            )
        if task.op not in _BOUNDARY_TASKS:
            raise CasegenError(
                f"boundary '{task.op}' has no hidden-schedule collision task"
            )
        return _BOUNDARY_TASKS[task.op]
    raise CasegenError(f"unknown operator or boundary '{task.op}'")


# ---------------------------------------------------------------------------
# C++ fragments
# ---------------------------------------------------------------------------


def _box_expr(region: Region, dim: int) -> str:
    box = tuple(float(v) for v in region.box)
    if len(box) != 2 * dim:
        raise CasegenError(
            f"region '{region.name}': a {dim}D box needs {2 * dim} endpoints, "
            f"got {region.box}"
        )
    nums = [_num(v, f"region '{region.name}' box") for v in box]

    def corner(start: int) -> str:
        return ", ".join("T(%s) * Cell_Len" % v for v in nums[start:start + dim])

    return (
        "Vector<T, %d>{%s}, Vector<T, %d>{%s}"
        % (dim, corner(0), dim, corner(dim))
    )


def _aabb_lines(spec: CaseSpec, region_ids: List[str]) -> List[str]:
    lines = ["// geometry: boxes are AABB endpoints in cell units; cell i has centre"]
    lines.append("// (i + 0.5) * Cell_Len, so an endpoint of n covers cell n - 1.")
    for index, region in enumerate(spec.regions):
        if region.shape != "aabb":
            raise CasegenError(
                f"region '{region.name}': shape {region.shape!r} is not supported"
            )
        lines.append(
            "AABB<T, %d> %s(%s);"
            % (spec.dim, region_ids[index], _box_expr(region, spec.dim))
        )
    return lines


def _cpp_string(value: str) -> str:
    """A C++ string literal for a path that must be used verbatim."""
    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    return '"%s"' % escaped


def _geometry_lines(spec: CaseSpec, flag_ids: Dict[str, str]) -> List[str]:
    """Read the optional voxel geometry and apply its domain/void mask.

    ``writeVoxelFlags`` overwrites the whole flag field, so this runs after the
    region actions and before the flag-count self-check.
    """
    geometry = spec.geometry
    if geometry is None:
        return []
    if spec.dim != 3:
        raise CasegenError(
            f"[geometry] is only supported for dim = 3, got dim = {spec.dim}"
        )
    for key in ("domain_flag", "void_flag"):
        name = getattr(geometry, key)
        if name not in flag_ids:
            raise CasegenError(f"[geometry] '{key}' names undefined flag '{name}'")
    return [
        "{",
        "  auto geometry_file = flbgeo::readGeometryFile<T>(%s);"
        % _cpp_string(geometry.file),
        "  flbgeo::writeVoxelFlags(geometry_file.Geometry, FlagFM, %s, %s);"
        % (flag_ids[geometry.domain_flag], flag_ids[geometry.void_flag]),
        "}",
    ]


def _flag_action_lines(
    spec: CaseSpec, region_ids: List[str], flag_ids: Dict[str, str]
) -> List[str]:
    lines: List[str] = []
    for index, region in enumerate(spec.regions):
        if region.action == "set_velocity":
            continue
        if region.flag is None or region.flag not in flag_ids:
            raise CasegenError(
                f"region '{region.name}': unknown flag {region.flag!r}"
            )
        target = flag_ids[region.flag]
        box = region_ids[index]
        if region.action == "paint":
            if region.only_if is None:
                lines.append(
                    "FlagFM.forEach(%s, [&](FLAG& field, std::size_t id) { "
                    "field.SetField(id, %s); });" % (box, target)
                )
            else:
                if region.only_if not in flag_ids:
                    raise CasegenError(
                        f"region '{region.name}': unknown only_if flag {region.only_if!r}"
                    )
                lines.append(
                    "FlagFM.forEach(%s, [&](FLAG& field, std::size_t id) {" % box
                )
                lines.append(
                    "  if (util::isFlag(field.get(id), %s)) field.SetField(id, %s);"
                    % (flag_ids[region.only_if], target)
                )
                lines.append("});")
        elif region.action == "boundary":
            lines.append(
                "FlagFM.template SetupBoundary<LatSet>(%s, %s);" % (box, target)
            )
        else:
            raise CasegenError(
                f"region '{region.name}': action {region.action!r} is not supported"
            )
    return lines


def _velocity_lines(
    spec: CaseSpec, region_ids: List[str], flag_ids: Dict[str, str]
) -> List[str]:
    lines: List[str] = []
    for index, region in enumerate(spec.regions):
        if region.action != "set_velocity":
            continue
        if region.flag is None or region.flag not in flag_ids:
            raise CasegenError(
                f"region '{region.name}': unknown flag {region.flag!r}"
            )
        dim = spec.dim
        if region.value is None or len(tuple(region.value)) != dim:
            raise CasegenError(
                f"region '{region.name}': set_velocity needs a {dim}-vector value"
            )
        comps = ", ".join(
            "T(%s)" % _num(v, f"region '{region.name}' value") for v in region.value
        )
        name = "LatU_%d_%s" % (index, _pascal(region.name))
        lines.append(
            "Vector<T, LatSet::d> %s(BaseConv.getLatticeU("
            "Vector<T, LatSet::d>{%s}));" % (name, comps)
        )
        lines.append("NSLattice.getField<VELOCITY<T, LatSet::d>>().forEach(")
        lines.append("    %s, FlagFM, %s," % (region_ids[index], flag_ids[region.flag]))
        lines.append(
            "    [&](auto& field, std::size_t id) { field.SetField(id, %s); });"
            % name
        )
    return lines


def _count_lines(spec: CaseSpec, flag_ids: Dict[str, str]) -> List[str]:
    names = list(spec.flags)
    void_name = spec.void_name
    lines: List[str] = []
    for name in names:
        lines.append("std::size_t count_%s = 0;" % _pascal(name))
    lines.append("FlagFM.forEach([&](FLAG& field, std::size_t id) {")
    lines.append("  const std::uint8_t f = field.get(id);")
    for name in names:
        lines.append(
            "  if (util::isFlag(f, %s)) ++count_%s;"
            % (flag_ids[name], _pascal(name))
        )
    lines.append("});")
    lines.append("")
    lines.append("auto globalCount = [](std::size_t local) -> std::size_t {")
    lines.append("#ifdef MPI_ENABLED")
    lines.append("  unsigned long long in = static_cast<unsigned long long>(local);")
    lines.append("  unsigned long long total = 0;")
    lines.append(
        "  MPI_Allreduce(&in, &total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);"
    )
    lines.append("  return static_cast<std::size_t>(total);")
    lines.append("#else")
    lines.append("  return local;")
    lines.append("#endif")
    lines.append("};")
    for name in names:
        lines.append("count_%s = globalCount(count_%s);" % (_pascal(name), _pascal(name)))
    lines.append("")
    print_parts = [
        '"%s = " << count_%s << "  "' % (name, _pascal(name)) for name in names
    ]
    lines.append('std::cout << "[flags] "')
    for part in print_parts:
        lines.append("          << %s" % part)
    lines.append("          << std::endl;")
    for name in names:
        if name == void_name:
            continue
        lines.append(
            'if (count_%s == 0) { std::cerr << "[flags] flag %s has no cells\\n"; '
            "return 1; }" % (_pascal(name), name)
        )
    return lines


def _task_lines(
    spec: CaseSpec, flag_ids: Dict[str, str]
) -> Tuple[List[str], List[Tuple[int, Task]]]:
    """Collision and rho/U task aliases, plus the boundary schemes.

    The hidden schedule turns bounce-back into collision tasks; plain and
    periodic hand those schemes to a BlockBoundaryManager instead, so their
    tasks are returned separately and never enter the TaskSelector.
    """
    manager = _manager_schedule(spec)
    order = sorted(
        enumerate(spec.tasks), key=lambda item: (-len(item[1].flags), item[0])
    )
    collision: List[str] = []
    names: List[str] = []
    rho_u: List[str] = []
    rho_u_names: List[str] = []
    boundary_entries: List[Tuple[int, Task]] = []
    for index, task in order:
        if manager and task.op in BOUNDARIES:
            boundary_entries.append((index, task))
            continue
        suffix = "_".join(_pascal(f) for f in task.flags)
        name = "Task_%d_%s" % (index, suffix)
        cpp = _task_cpp(spec, task)
        key = _key_expr(task.flags, flag_ids)
        collision.append("using %s = tmp::Key_TypePair<%s, %s>;" % (name, key, cpp))
        names.append(name)
        if task.op in OPERATORS or task.op in spec.operator_specs:
            rho_name = "RhoUTask_%d_%s" % (index, suffix)
            rho_u.append(
                "using %s = tmp::Key_TypePair<%s, moment::rhoU<CELL, true>>;"
                % (rho_name, key)
            )
            rho_u_names.append(rho_name)

    if not rho_u_names:
        raise CasegenError(
            "[tasks] has no collision operator; rho/U output cannot be formed"
        )

    lines = collision + rho_u + [
        "",
        "using TaskCollection = tmp::TupleWrapper<%s>;" % ", ".join(names),
        "using NSTask = tmp::TaskSelector<TaskCollection, std::uint8_t, CELL>;",
        "using TaskCollectionRhoU = tmp::TupleWrapper<%s>;"
        % ", ".join(rho_u_names),
        "using TaskSelectorRhoU = tmp::TaskSelector<TaskCollectionRhoU, "
        "std::uint8_t, CELL>;",
    ]
    return lines, boundary_entries


def _void_flag_expr(spec: CaseSpec, flag_ids: Dict[str, str]) -> str:
    if spec.void_name is not None:
        return flag_ids[spec.void_name]
    return "std::uint8_t(1)"


def _boundary_manager_lines(
    spec: CaseSpec,
    flag_ids: Dict[str, str],
    boundary_entries: List[Tuple[int, Task]],
) -> List[str]:
    """BBLikeFixedBlockBdManager per boundary scheme + one BlockBoundaryManager."""
    lines: List[str] = []
    args: List[str] = []
    for index, task in boundary_entries:
        suffix = "_".join(_pascal(f) for f in task.flags)
        name = "NS_Bd_%d_%s" % (index, suffix)
        cpp = BOUNDARIES[task.op]["cpp"]
        key = _key_expr(task.flags, flag_ids)
        lines.append(
            "BBLikeFixedBlockBdManager<%s, "
            "BlockLatticeManager<T, LatSet, FIELDS>, "
            "BlockFieldManager<FLAG, T, LatSet::d>>" % cpp
        )
        lines.append(
            '  %s("%s", NSLattice, FlagFM, %s, %s);'
            % (name, name, key, _void_flag_expr(spec, flag_ids))
        )
        args.append(name)
    if args:
        lines.append("BlockBoundaryManager BM(%s);" % ", ".join("&" + a for a in args))
    return lines


_PERIODIC_DIRECTIONS = {"x": ("XN", "XP"), "y": ("YN", "YP"), "z": ("ZN", "ZP")}


def _periodic_face_boxes(spec: CaseSpec) -> Tuple[Tuple[float, ...], Tuple[float, ...]]:
    """Face strips in cell units: one transverse layer excluded per side."""
    axis = "xyz".index(spec.schedule["periodic"])
    n = [int(v) for v in spec.n]
    lo0 = [1] * spec.dim
    hi0 = [n[i] - 1 for i in range(spec.dim)]
    lo0[axis], hi0[axis] = 0, 1
    lo1 = [1] * spec.dim
    hi1 = [n[i] - 1 for i in range(spec.dim)]
    lo1[axis], hi1[axis] = n[axis] - 1, n[axis]
    return (
        tuple(float(v) for v in lo0 + hi0),
        tuple(float(v) for v in lo1 + hi1),
    )


def _periodic_lines(spec: CaseSpec, flag_ids: Dict[str, str]) -> List[str]:
    axis = spec.schedule["periodic"]
    neg_dir, pos_dir = _PERIODIC_DIRECTIONS[axis]
    lo_box, hi_box = _periodic_face_boxes(spec)
    dim = spec.dim

    def corner(box: Tuple[float, ...]) -> str:
        nums = [_num(v, "[schedule] periodic box") for v in box]
        return ", ".join("T(%s) * Cell_Len" % v for v in nums)

    lines = [
        "FixedPeriodicBoundaryManager<BlockLatticeManager<T, LatSet, FIELDS>, "
        "BlockFieldManager<FLAG, T, LatSet::d>>",
        '  NS_Peri("NS_Peri", NSLattice, FlagFM, %s, %s);'
        % (flag_ids["periodic"], _void_flag_expr(spec, flag_ids)),
        "AABB<T, %d> periodic_lo(Vector<T, %d>{%s}, Vector<T, %d>{%s});"
        % (dim, dim, corner(lo_box[:dim]), dim, corner(lo_box[dim:])),
        "AABB<T, %d> periodic_hi(Vector<T, %d>{%s}, Vector<T, %d>{%s});"
        % (dim, dim, corner(hi_box[:dim]), dim, corner(hi_box[dim:])),
        "NS_Peri.Setup(periodic_lo, NbrDirection::%s, periodic_hi, NbrDirection::%s);"
        % (neg_dir, pos_dir),
    ]
    return lines


def _writer_lines(spec: CaseSpec) -> Tuple[List[str], List[str]]:
    raw = _as_list(spec.output.get("writers"))
    writers: List[str] = []
    for item in raw:
        name = str(item).strip().lower()
        if not name:
            continue
        if name not in {"rho", "velocity", "flag"}:
            raise CasegenError(f"[output] writers: unknown writer '{name}'")
        if name not in writers:
            writers.append(name)
    if _truthy(spec.output.get("flag_writer", 0)) and "flag" not in writers:
        writers.append("flag")

    decls: List[str] = []
    built: List[str] = []
    if "rho" in writers:
        if "RHO" not in spec.fields:
            raise CasegenError("[output] rho writer requires RHO in [fields]")
        decls.append(
            'vtmo::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());'
        )
        built.append("RhoWriter")
    if "velocity" in writers:
        decls.append(
            'vtmo::VectorWriter VelocityWriter("Velocity", '
            "NSLattice.getField<VELOCITY<T, LatSet::d>>());"
        )
        built.append("VelocityWriter")
    if "flag" in writers:
        decls.append('vtmo::ScalarWriter FlagWriter("flag", FlagFM);')
        built.append("FlagWriter")

    if not built:
        return [], []
    use = ['vtmo::vtmWriter<T, LatSet::d> NSWriter("%s", Geo);' % spec.name]
    use.append("NSWriter.addWriterSet(%s);" % ", ".join(built))
    return decls + use, ["NSWriter.WriteBinary(MainLoopTimer());"]


def _mpi_inner_pop_lines() -> List[str]:
    """POP sum/sumsq/maxabs over every block's interior, then MPI_Allreduce.

    The inner region is ``overlap .. N - overlap`` in every direction: those
    cells tile the domain exactly once, while the halos are duplicated across
    blocks and ranks.
    """
    lines = [
        "double pop_sum = 0.0, pop_sumsq = 0.0, pop_maxabs = 0.0;",
        "for (auto& bl : NSLattice.getBlockLats()) {",
        "  const auto& popf = bl.getField<POP<T, LatSet::q>>();",
        "  const int overlap = bl.getOverlap();",
        "  if constexpr (LatSet::d == 2) {",
        "    for (int j = overlap; j < bl.getNy() - overlap; ++j) {",
        "      for (int i = overlap; i < bl.getNx() - overlap; ++i) {",
        "        const std::size_t id = std::size_t(j) * bl.getNx() + i;",
        "        for (unsigned int d = 0; d < LatSet::q; ++d) {",
        "          const double v = static_cast<double>(popf.getField(d)[id]);",
        "          pop_sum += v;",
        "          pop_sumsq += v * v;",
        "          const double a = v < 0.0 ? -v : v;",
        "          pop_maxabs = a > pop_maxabs ? a : pop_maxabs;",
        "        }",
        "      }",
        "    }",
        "  } else {",
        "    for (int k = overlap; k < bl.getNz() - overlap; ++k) {",
        "      for (int j = overlap; j < bl.getNy() - overlap; ++j) {",
        "        for (int i = overlap; i < bl.getNx() - overlap; ++i) {",
        "          const std::size_t id = std::size_t(k) * bl.getProjection()[2] +",
        "                                 std::size_t(j) * bl.getProjection()[1] + i;",
        "          for (unsigned int d = 0; d < LatSet::q; ++d) {",
        "            const double v = static_cast<double>(popf.getField(d)[id]);",
        "            pop_sum += v;",
        "            pop_sumsq += v * v;",
        "            const double a = v < 0.0 ? -v : v;",
        "            pop_maxabs = a > pop_maxabs ? a : pop_maxabs;",
        "          }",
        "        }",
        "      }",
        "    }",
        "  }",
        "}",
        "const std::array<double, 3> pop_local{pop_sum, pop_sumsq, pop_maxabs};",
        "std::array<double, 3> pop_s{0.0, 0.0, 0.0};",
        "MPI_Allreduce(pop_local.data(), pop_s.data(), 1, MPI_DOUBLE, MPI_SUM,",
        "              MPI_COMM_WORLD);",
        "MPI_Allreduce(pop_local.data() + 1, pop_s.data() + 1, 1, MPI_DOUBLE,",
        "              MPI_SUM, MPI_COMM_WORLD);",
        "MPI_Allreduce(pop_local.data() + 2, pop_s.data() + 2, 1, MPI_DOUBLE,",
        "              MPI_MAX, MPI_COMM_WORLD);",
        'frelb_diag::PrintChecksum("pop inner, MPI", pop_s);',
    ]
    return lines


def _mpi_inner_rho_lines() -> List[str]:
    """RHO sum/sumsq/maxabs over every block's interior, then MPI_Allreduce."""
    lines = [
        "double rho_sum = 0.0, rho_sumsq = 0.0, rho_maxabs = 0.0;",
        "for (auto& bl : NSLattice.getBlockLats()) {",
        "  const auto& rhof = bl.getField<RHO<T>>();",
        "  const int overlap = bl.getOverlap();",
        "  if constexpr (LatSet::d == 2) {",
        "    for (int j = overlap; j < bl.getNy() - overlap; ++j) {",
        "      for (int i = overlap; i < bl.getNx() - overlap; ++i) {",
        "        const std::size_t id = std::size_t(j) * bl.getNx() + i;",
        "        const double v = static_cast<double>(rhof.getField(0)[id]);",
        "        rho_sum += v;",
        "        rho_sumsq += v * v;",
        "        const double a = v < 0.0 ? -v : v;",
        "        rho_maxabs = a > rho_maxabs ? a : rho_maxabs;",
        "      }",
        "    }",
        "  } else {",
        "    for (int k = overlap; k < bl.getNz() - overlap; ++k) {",
        "      for (int j = overlap; j < bl.getNy() - overlap; ++j) {",
        "        for (int i = overlap; i < bl.getNx() - overlap; ++i) {",
        "          const std::size_t id = std::size_t(k) * bl.getProjection()[2] +",
        "                                 std::size_t(j) * bl.getProjection()[1] + i;",
        "          const double v = static_cast<double>(rhof.getField(0)[id]);",
        "          rho_sum += v;",
        "          rho_sumsq += v * v;",
        "          const double a = v < 0.0 ? -v : v;",
        "          rho_maxabs = a > rho_maxabs ? a : rho_maxabs;",
        "        }",
        "      }",
        "    }",
        "  }",
        "}",
        "const std::array<double, 3> rho_local{rho_sum, rho_sumsq, rho_maxabs};",
        "std::array<double, 3> rho_s{0.0, 0.0, 0.0};",
        "MPI_Allreduce(rho_local.data(), rho_s.data(), 1, MPI_DOUBLE, MPI_SUM,",
        "              MPI_COMM_WORLD);",
        "MPI_Allreduce(rho_local.data() + 1, rho_s.data() + 1, 1, MPI_DOUBLE,",
        "              MPI_SUM, MPI_COMM_WORLD);",
        "MPI_Allreduce(rho_local.data() + 2, rho_s.data() + 2, 1, MPI_DOUBLE,",
        "              MPI_MAX, MPI_COMM_WORLD);",
        'frelb_diag::PrintChecksum("rho inner, MPI", rho_s);',
    ]
    return lines


def _checksum_lines(spec: CaseSpec, cuda: bool = False) -> List[str]:
    checksum = spec.checksum
    if not _truthy(checksum.get("enable", 0)):
        return []
    tol = float(checksum.get("tol", 1e-6))
    expect_pop = checksum.get("expect_pop")
    expect_rho = checksum.get("expect_rho")

    lines: List[str] = []
    if cuda:
        lines.append("// Final fields for the checksums: one rho/U pass on the device")
        lines.append("// and a synchronous pull, so the host accumulators read exactly")
        lines.append("// the values the kernels produced.")
        lines.append(
            "NSLattice.CuDevApplyCellDynamics<TaskSelectorRhoU, CELL>(FlagFM);"
        )
        lines.append("cudaDeviceSynchronize();")
        if "POP" in spec.fields:
            lines.append(
                "NSLattice.getBlockLat(0).getField<POP<T, LatSet::q>>().copyToHost();"
            )
        if "RHO" in spec.fields:
            lines.append("NSLattice.getBlockLat(0).getField<RHO<T>>().copyToHost();")
        lines.append("{")
        lines.append("  cudaError_t err = cudaGetLastError();")
        lines.append("  if (err != cudaSuccess) {")
        lines.append(
            '    std::cerr << "[%s] final rho/U launch failed: " '
            "<< cudaGetErrorString(err) << std::endl;" % spec.name
        )
        lines.append("    return 3;")
        lines.append("  }")
        lines.append("}")
    else:
        lines.append("// Final fields for the checksums: communication first, then one")
        lines.append("// rho/U pass, exactly as in tests/regression/cavity2d.")
        lines.append("NSLattice.NormalFullCommunicate();")
        lines.append("NSLattice.ApplyCellDynamics<TaskSelectorRhoU>(FlagFM);")
    lines.append("")

    has_pop = "POP" in spec.fields
    has_rho = "RHO" in spec.fields
    mpi = bool(getattr(spec, "mpi", False))

    if has_pop:
        if mpi:
            lines.append("#ifdef MPI_ENABLED")
            lines.extend(_mpi_inner_pop_lines())
            lines.append("#else")
        lines.append("std::array<double, 3> pop_s{0.0, 0.0, 0.0};")
        lines.append("for (auto& bl : NSLattice.getBlockLats()) {")
        lines.append(
            "  const std::array<double, 3> s = frelb_diag::PopChecksumHost("
            "bl.getField<POP<T, LatSet::q>>(), bl.getN(), LatSet::q);"
        )
        lines.append("  pop_s[0] += s[0];")
        lines.append("  pop_s[1] += s[1];")
        lines.append("  pop_s[2] = s[2] > pop_s[2] ? s[2] : pop_s[2];")
        lines.append("}")
        lines.append('frelb_diag::PrintChecksum("pop", pop_s);')
        if mpi:
            lines.append("#endif")
        lines.append("")

    if has_rho:
        if mpi:
            lines.append("#ifdef MPI_ENABLED")
            lines.extend(_mpi_inner_rho_lines())
            lines.append("#else")
        lines.append("// RHO uses the same measure and the same all-block-cell order as")
        lines.append("// the pop checksum: the field manager's default forEach.")
        lines.append("double rho_sum = 0.0, rho_sumsq = 0.0, rho_maxabs = 0.0;")
        lines.append("NSLattice.getField<RHO<T>>().forEach([&](auto& field, std::size_t id) {")
        lines.append("  const double v = static_cast<double>(field.get(id));")
        lines.append("  rho_sum += v;")
        lines.append("  rho_sumsq += v * v;")
        lines.append("  const double a = v < 0.0 ? -v : v;")
        lines.append("  rho_maxabs = a > rho_maxabs ? a : rho_maxabs;")
        lines.append("});")
        lines.append(
            "const std::array<double, 3> rho_s{rho_sum, rho_sumsq, rho_maxabs};"
        )
        lines.append('frelb_diag::PrintChecksum("rho", rho_s);')
        if mpi:
            lines.append("#endif")
        lines.append("")

    if expect_pop is not None:
        lines.append(
            "constexpr double kExpectPop = %s;" % _num(expect_pop, "[checksum] expect_pop")
        )
    if expect_rho is not None:
        lines.append(
            "constexpr double kExpectRho = %s;" % _num(expect_rho, "[checksum] expect_rho")
        )
    expect_pop_mpi = checksum.get("expect_pop_mpi")
    expect_rho_mpi = checksum.get("expect_rho_mpi")
    if mpi and expect_pop_mpi is not None:
        lines.append(
            "constexpr double kExpectPopMpi = %s;"
            % _num(expect_pop_mpi, "[checksum] expect_pop_mpi")
        )
    if mpi and expect_rho_mpi is not None:
        lines.append(
            "constexpr double kExpectRhoMpi = %s;"
            % _num(expect_rho_mpi, "[checksum] expect_rho_mpi")
        )
    lines.append("constexpr double kChecksumTol = %s;" % _num(tol, "[checksum] tol"))
    lines.append("")

    finite_terms = []
    if has_pop:
        finite_terms.append("std::isfinite(pop_s[0])")
        finite_terms.append("std::isfinite(pop_s[1])")
        finite_terms.append("std::isfinite(pop_s[2])")
    if has_rho:
        finite_terms.append("std::isfinite(rho_s[0])")
        finite_terms.append("std::isfinite(rho_s[1])")
        finite_terms.append("std::isfinite(rho_s[2])")
    if finite_terms:
        lines.append("bool ok = " + " &&\n          ".join(finite_terms) + ";")
    else:  # pragma: no cover - enable without POP/RHO is rejected earlier
        lines.append("bool ok = true;")

    has_mpi_expect = expect_pop_mpi is not None or expect_rho_mpi is not None
    has_plain_expect = expect_pop is not None or expect_rho is not None
    if mpi and (has_mpi_expect or has_plain_expect):
        lines.append("#ifdef MPI_ENABLED")
    if mpi and expect_pop_mpi is not None:
        lines.append(
            "if (std::abs(pop_s[1] - kExpectPopMpi) > "
            "kChecksumTol * std::abs(kExpectPopMpi)) ok = false;"
        )
    if mpi and expect_rho_mpi is not None:
        lines.append(
            "if (std::abs(rho_s[1] - kExpectRhoMpi) > "
            "kChecksumTol * std::abs(kExpectRhoMpi)) ok = false;"
        )
    if mpi and has_plain_expect:
        lines.append("#else")
    if expect_pop is not None:
        lines.append(
            "if (std::abs(pop_s[1] - kExpectPop) > "
            "kChecksumTol * std::abs(kExpectPop)) ok = false;"
        )
    if expect_rho is not None:
        lines.append(
            "if (std::abs(rho_s[1] - kExpectRho) > "
            "kChecksumTol * std::abs(kExpectRho)) ok = false;"
        )
    if mpi and (has_mpi_expect or has_plain_expect):
        lines.append("#endif")
    lines.append("if (!ok) {")
    lines.append(
        '  std::cerr << "[checksum] mismatch (relative tolerance " '
        '<< kChecksumTol << ")" << std::endl;'
    )
    lines.append("  return 1;")
    lines.append("}")
    return lines


# ---------------------------------------------------------------------------
# C++ file
# ---------------------------------------------------------------------------


def _render_cpp(spec: CaseSpec) -> str:
    cuda = bool(spec.cuda)
    flag_ids = _flag_ids(spec)
    region_ids = [_region_id(r, i) for i, r in enumerate(spec.regions)]

    schedule = spec.schedule
    steps = _int(_require(schedule, "steps", "[schedule]"), "[schedule] steps")
    output_step = _int(
        _require(schedule, "output_step", "[schedule]"), "[schedule] output_step"
    )
    tol = _num(_require(schedule, "tol", "[schedule]"), "[schedule] tol")
    checkpoint_step = _int(schedule.get("checkpoint_step", 0), "[schedule] checkpoint_step")
    if steps <= 0:
        raise CasegenError("[schedule] steps must be positive")
    if output_step <= 0:
        raise CasegenError("[schedule] output_step must be positive")

    physics = spec.physics
    converter = physics.get("converter", "simplified")
    dim = spec.dim
    init_u = tuple(float(v) for v in physics.get("init_u", (0.0,) * dim))
    if len(init_u) != dim:
        raise CasegenError(
            f"[physics] init_u must have {dim} entries for a {dim}D case"
        )
    init_u_zero = all(v == 0.0 for v in init_u)

    lines: List[str] = []
    lines.append("// %s.%s" % (spec.name, "cu" if cuda else "cpp"))
    lines.append("// generated by casegen from %s -- do not edit"
                 % os.path.basename(spec.path or "?"))
    lines.append("// spec sha256: %s" % (spec.source_hash or "unknown"))
    lines.append("")
    lines.append('#include "freelb.h"')
    lines.append('#include "freelb.hh"')
    custom_ops = sorted(
        {task.op for task in spec.tasks if task.op in spec.operator_specs}
    )
    for name in custom_ops:
        lines.append('#include "operators/%s.hh"' % name)
    lines.append("")
    lines.append("#include <array>")
    lines.append("#include <cmath>")
    lines.append("#include <cstdint>")
    lines.append("#include <iostream>")
    lines.append("#include <string>")
    lines.append("")
    lines.append("using T = FLOAT;")
    lines.append("using LatSet = %s<T>;" % spec.latset)
    lines.append("")
    lines.append("// run-time parameters, read from %s.ini" % spec.name)
    lines.append("T RT;")
    if converter == "simplified":
        lines.append("T CharU;")
    else:
        lines.append("T rho_ref;")
        lines.append("T Kine_Visc;")
        lines.append("T U_Max;")
    lines.append("int Ni;")
    lines.append("int Nj;")
    if dim == 3:
        lines.append("int Nk;")
    lines.append("T Cell_Len;")
    lines.append("int BlockCellLen;")
    lines.append("int Thread_Num;")
    lines.append("int MaxStep;")
    lines.append("int OutputStep;")
    lines.append("T tol;")
    lines.append("int Restart;")
    lines.append("int CheckpointStep;")
    lines.append("std::string RestartDir;")
    lines.append("")
    lines.append("void readParam() {")
    lines.append('  iniReader param_reader("%s.ini");' % spec.name)
    lines.append('  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");')
    lines.append('  Ni = param_reader.getValue<int>("Mesh", "Ni");')
    lines.append('  Nj = param_reader.getValue<int>("Mesh", "Nj");')
    if dim == 3:
        lines.append('  Nk = param_reader.getValue<int>("Mesh", "Nk");')
    lines.append('  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");')
    lines.append(
        '  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");'
    )
    lines.append('  RT = param_reader.getValue<T>("LB", "RT");')
    if converter == "simplified":
        lines.append('  CharU = param_reader.getValue<T>("LB", "CharU");')
    else:
        lines.append(
            '  rho_ref = param_reader.getValue<T>("Physical_Property", "rho_ref");'
        )
        lines.append(
            '  Kine_Visc = param_reader.getValue<T>("Physical_Property", "Kine_Visc");'
        )
        lines.append('  U_Max = param_reader.getValue<T>("Init_Conditions", "U_Max");')
    lines.append(
        '  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");'
    )
    lines.append(
        '  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");'
    )
    lines.append('  tol = param_reader.getValue<T>("tolerance", "tol");')
    lines.append(
        '  Restart = param_reader.getValueOr<int>("Simulation_Settings", "Restart", 0);'
    )
    lines.append(
        '  CheckpointStep = param_reader.getValueOr<int>('
        '"Simulation_Settings", "CheckpointStep", 0);'
    )
    lines.append(
        '  RestartDir = param_reader.getValueOr<std::string>('
        '"Simulation_Settings", "RestartDir", "./ckpt");'
    )
    lines.append("")
    lines.append("  MPI_RANK(0)")
    lines.append('  std::cout << "------------Simulation Parameters:-------------\\n"')
    mesh_line = '            << "Mesh:              " << Ni << " x " << Nj'
    if dim == 3:
        mesh_line += ' << " x " << Nk'
    mesh_line += ' << "\\n"'
    lines.append(mesh_line)
    lines.append('            << "Cell_Len:          " << Cell_Len << "\\n"')
    lines.append('            << "BlockCellLen:      " << BlockCellLen << "\\n"')
    lines.append('            << "RT:                " << RT << "\\n"')
    if converter == "simplified":
        lines.append('            << "CharU:             " << CharU << "\\n"')
    lines.append('            << "TotalStep:         " << MaxStep << "\\n"')
    lines.append('            << "OutputStep:        " << OutputStep << "\\n"')
    lines.append('            << "Tolerance:         " << tol << "\\n"')
    lines.append(
        '            << "Restart:           " << (Restart ? "on" : "off") << "\\n"'
    )
    lines.append('            << "CheckpointStep:    " << CheckpointStep << "\\n"')
    lines.append('            << "RestartDir:        " << RestartDir << "\\n"')
    lines.append('            << "----------------------------------------------"')
    lines.append("            << std::endl;")
    lines.append("}")
    lines.append("")
    lines.append("int main(int argc, char* argv[]) {")

    # flag constants
    lines.append("  // flag bits from [flags]")
    for name in spec.flags:
        bit = _int(spec.flags[name], f"[flags] {name}")
        if not 0 < bit < 256:
            raise CasegenError(f"[flags] {name} = {bit} is outside uint8_t")
        lines.append(
            "  constexpr std::uint8_t %s = std::uint8_t(%d);" % (flag_ids[name], bit)
        )
    lines.append("")

    lines.append("  mpi().init(&argc, &argv);")
    if _manager_schedule(spec):
        # The plain schedule has only been validated as a serial run.  Its
        # BlockBoundaryManager and Communicate ordering under MPI is untested,
        # and a wrong ordering corrupts the result silently, so refuse ranks
        # instead of producing numbers nobody has checked.
        lines.append("  // plain + MPI has not been validated; fail before simulating")
        lines.append("  if (mpi().getSize() > 1) {")
        lines.append(
            '    std::cerr << "[%s] plain schedule with MPI is not validated; '
            'run with a single rank" << std::endl;' % spec.name
        )
        lines.append("    return 1;")
        lines.append("  }")
    lines.append("  MPI_DEBUG_WAIT;")
    lines.append('  Printer::Print_BigBanner(std::string("Initializing..."));')
    lines.append("  readParam();")
    lines.append("")
    lines.append("  BaseConverter<T> BaseConv(LatSet::cs2);")
    if converter == "simplified":
        char_u = _num(
            _require(physics, "char_u", "[physics]"), "[physics] char_u"
        )
        rt = _num(_require(physics, "rt", "[physics]"), "[physics] rt")
        lines.append(
            "  BaseConv.SimplifiedConverterFromRT(Ni, T(%s), T(%s));" % (char_u, rt)
        )
    else:
        rt = _num(_require(physics, "rt", "[physics]"), "[physics] rt")
        lines.append(
            "  BaseConv.ConvertFromRT(Cell_Len, T(%s), rho_ref, T(Ni) * Cell_Len, "
            "U_Max, Kine_Visc);" % rt
        )
    lines.append("")
    lines.append("  // ------------------ geometry ------------------")
    domain_coords = ", ".join(
        "T(%s) * Cell_Len" % name for name in ("Ni", "Nj", "Nk")[:dim]
    )
    lines.append(
        "  AABB<T, %d> domain(Vector<T, %d>{}, Vector<T, %d>{%s});"
        % (dim, dim, dim, domain_coords)
    )
    lines.extend("  " + line for line in _aabb_lines(spec, region_ids))
    lines.append("")
    # An auto layout has to be the same for every rank count: the checksums
    # sum each block's interior, and AdaptiveOptimization(np) re-splits the
    # domain per np (for 32x32/np=2 it produced two 30x15 blocks, dropping the
    # outer ring).  CreateBlocks(bx, by[, bz]) is rank independent;
    # LoadBalancing(np) then only decides which rank owns which block.
    if spec.decompose == "single":
        if dim == 3:
            lines.append("  BlockGeometry3D<T> Geo(Ni, Nj, Nk, 1, domain, Cell_Len);")
        else:
            lines.append("  BlockGeometry2D<T> Geo(Ni, Nj, 1, domain, Cell_Len);")
    elif dim == 3:
        bx = _block_count(spec.n[0], spec.block_cell_len)
        by = _block_count(spec.n[1], spec.block_cell_len)
        bz = _block_count(spec.n[2], spec.block_cell_len)
        lines.append(
            "  BlockGeometryHelper3D<T> GeoHelper(Ni, Nj, Nk, domain, Cell_Len, BlockCellLen);"
        )
        lines.append("  GeoHelper.CreateBlocks(%d, %d, %d);" % (bx, by, bz))
        lines.append("  GeoHelper.LoadBalancing(mpi().getSize());")
        lines.append("  BlockGeometry3D<T> Geo(GeoHelper);")
    else:
        bx = _block_count(spec.n[0], spec.block_cell_len)
        by = _block_count(spec.n[1], spec.block_cell_len)
        lines.append(
            "  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);"
        )
        lines.append("  GeoHelper.CreateBlocks(%d, %d);" % (bx, by))
        lines.append("  GeoHelper.LoadBalancing(mpi().getSize());")
        lines.append("  BlockGeometry2D<T> Geo(GeoHelper);")
    lines.append("")
    lines.append("  // ------------------ flag field ------------------")
    if spec.void_name:
        lines.append(
            "  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, %s);"
            % flag_ids[spec.void_name]
        )
    else:
        lines.append(
            "  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, std::uint8_t(1));"
        )
    lines.extend("  " + line for line in _flag_action_lines(spec, region_ids, flag_ids))
    geometry_lines = _geometry_lines(spec, flag_ids)
    if geometry_lines:
        lines.append("")
        lines.append("  // geometry file: the voxel record owns the domain/void split")
        lines.extend("  " + line for line in geometry_lines)
    lines.append("")
    lines.append("  // spec-vs-geometry self-check: a flag that painted nothing is a")
    lines.append("  // silently non-colliding region, so fail instead of running it.")
    lines.extend("  " + line for line in _count_lines(spec, flag_ids))
    lines.append("")
    lines.append("  // ------------------ lattice ------------------")
    fields = ", ".join(FIELD_TYPES[name] for name in spec.fields)
    lines.append("  using FIELDS = TypePack<%s>;" % fields)
    if cuda:
        # The kernels and the host share one task list built around the device
        # cell; the host cell of the same fields is only used by the writers
        # and the checksums, which go through the field managers.
        lines.append(
            "  using CELL = cudev::Cell<T, LatSet, "
            "typename ExtractCudevFieldPack<FIELDS>::cudev_pack>;"
        )
    else:
        lines.append("  using CELL = Cell<T, LatSet, FIELDS>;")
    init_values = {
        "RHO": "BaseConv.getLatRhoInit()",
        "VELOCITY": (
            "Vector<T, LatSet::d>{}"
            if init_u_zero
            else "BaseConv.getLatticeU(Vector<T, LatSet::d>{%s})"
            % ", ".join("T(%s)" % _num(v, "[physics] init_u") for v in init_u)
        ),
        "POP": "T{}",
        "FORCE": _force_expr(spec),
        "SMAGORINSKY": _smagorinsky_expr(spec),
    }
    missing_init = [name for name in spec.fields if name not in init_values]
    if missing_init:  # pragma: no cover - guarded by _check_supported
        raise CasegenError(
            "no initial value for field(s) " + ", ".join(missing_init)
        )
    lines.append(
        "  ValuePack InitValues(%s);"
        % ", ".join(init_values[name] for name in spec.fields)
    )
    lines.append(
        "  BlockLatticeManager<T, LatSet, FIELDS> NSLattice(Geo, InitValues, BaseConv);"
    )
    lines.append("")

    velocity_lines = _velocity_lines(spec, region_ids, flag_ids)
    if velocity_lines:
        lines.append("  // set_velocity actions run here: they write VELOCITY, which")
        lines.append("  // only exists after the lattice manager is built.")
        lines.extend("  " + line for line in velocity_lines)
        lines.append("")

    lines.append("  // ------------------ tasks ------------------")
    task_lines, boundary_entries = _task_lines(spec, flag_ids)
    lines.extend("  " + line for line in task_lines)
    lines.append("")
    if boundary_entries:
        lines.append("  // ------------------ boundary managers ------------------")
        lines.append("  // Plain/periodic bounce-back reads post-stream populations, so")
        lines.append("  // these schemes live outside the collision TaskSelector.")
        boundary_lines = _boundary_manager_lines(spec, flag_ids, boundary_entries)
        lines.extend("  " + line for line in boundary_lines)
        lines.append("")
    if spec.schedule.get("periodic") is not None:
        lines.append("  // ------------------ periodic boundary ------------------")
        periodic_lines = _periodic_lines(spec, flag_ids)
        lines.extend("  " + line for line in periodic_lines)
        lines.append("")

    writer_decls, writer_calls = _writer_lines(spec)
    if writer_decls:
        lines.append("  // ------------------ writers ------------------")
        lines.extend("  " + line for line in writer_decls)
        lines.append("")

    lines.append("  // ------------------ restart ------------------")
    lines.append("  // Geometry and flags are replayed above; only the fields and the")
    lines.append("  // step count are restored, so no BC object has to be rebuilt.")
    lines.append("  Timer MainLoopTimer;")
    lines.append("  if (Restart) {")
    lines.append(
        "    if (!flbckpt::checkCompatible(RestartDir, -1, NSLattice, FlagFM)) {"
    )
    lines.append(
        '      Printer::Print_BigBanner(std::string("Restart aborted: '
        'incompatible checkpoint"));'
    )
    lines.append("      return 1;")
    lines.append("    }")
    lines.append("    MainLoopTimer.set(flbckpt::read(RestartDir, -1, NSLattice, FlagFM));")
    lines.append("    IF_MPI_RANK(0) {")
    lines.append(
        '      std::cout << "[Restart] resumed at step " << MainLoopTimer() << std::endl;'
    )
    lines.append("    }")
    lines.append("  }")
    lines.append("")
    lines.append("  Timer OutputTimer;")
    if writer_calls:
        lines.extend("  " + line for line in writer_calls)
    lines.append('  Printer::Print_BigBanner(std::string("Start Calculation..."));')
    lines.append("")
    if cuda:
        # Everything that touches a field from the host (geometry replay,
        # set_velocity, a restart) has run by now, so this is the one point
        # where the device mirrors have to be refreshed.
        lines.append("  // ------------------ device mirror ------------------")
        lines.append("  // The kernels read the flag field and the population/field data")
        lines.append("  // from their device copies.")
        lines.append("  FlagFM.copyToDevice();")
        for name in spec.fields:
            lines.append(
                "  NSLattice.getField<%s>().copyToDevice();" % FIELD_TYPES[name]
            )
        lines.append("")
    lines.append("  NSLattice.EnableToleranceU();")
    lines.append("  T res = T(1);")
    lines.append("  while (MainLoopTimer() < MaxStep && res > tol) {")
    manager = _manager_schedule(spec)
    if cuda:
        lines.append("    NSLattice.CuDevApplyCellDynamics<NSTask, CELL>(FlagFM);")
        lines.append("    NSLattice.CuDevStream();")
    elif manager:
        # examples/openboundary2d/openbd2d.cpp, with the periodic copy done
        # between the collision and the stream (examples/couette2d).
        lines.append(
            "    NSLattice.ApplyCellDynamics<TaskSelectorRhoU>(MainLoopTimer(), FlagFM);"
        )
        lines.append("    NSLattice.ApplyCellDynamics<NSTask>(MainLoopTimer(), FlagFM);")
        if spec.schedule.get("periodic") is not None:
            lines.append("    NS_Peri.Apply();")
        lines.append("    NSLattice.Stream(MainLoopTimer());")
        if boundary_entries:
            lines.append("    BM.Apply(MainLoopTimer());")
        lines.append("    NSLattice.Communicate(MainLoopTimer());")
    else:
        lines.append("    NSLattice.Stream();")
        lines.append("    NSLattice.ApplyCellDynamicsWithHiddenComm<NSTask>(FlagFM);")
    lines.append("    ++MainLoopTimer;")
    lines.append("    ++OutputTimer;")
    lines.append("    if (MainLoopTimer() % OutputStep == 0) {")
    if cuda:
        lines.append(
            "      NSLattice.CuDevApplyCellDynamics<TaskSelectorRhoU, CELL>(FlagFM);"
        )
        lines.append("      cudaDeviceSynchronize();")
        # VELOCITY is required by the hidden schedule, RHO by the rho/U task;
        # both have to be on the host before the tolerance, the writers and
        # the checksums read them.
        for name in ("RHO", "VELOCITY"):
            if name in spec.fields:
                lines.append(
                    "      NSLattice.getBlockLat(0).getField<%s>().copyToHost();"
                    % FIELD_TYPES[name]
                )
    elif manager:
        lines.append("      NSLattice.ApplyCellDynamics<TaskSelectorRhoU>(FlagFM);")
    else:
        lines.append("      NSLattice.NormalFullCommunicate();")
        lines.append("      NSLattice.ApplyCellDynamics<TaskSelectorRhoU>(FlagFM);")
    lines.append("      res = NSLattice.getToleranceU(-1);")
    lines.append(
        "      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);"
    )
    lines.append("      Printer::Print_Res<T>(res);")
    lines.append("      Printer::Endl();")
    if writer_calls:
        lines.append("      NSWriter.WriteBinary(MainLoopTimer());")
    if cuda:
        lines.append("      // Quiescent snapshot point: this step's collision and streaming")
        lines.append("      // are done, the next collision has not started.")
    else:
        lines.append("      // Quiescent snapshot point: communication is done and the")
        lines.append("      // next Stream has not started.")
    lines.append(
        "      if (CheckpointStep > 0 && MainLoopTimer() % CheckpointStep == 0)"
    )
    lines.append("        flbckpt::write(RestartDir, MainLoopTimer(), NSLattice, FlagFM);")
    lines.append("    }")
    lines.append("  }")
    lines.append("")
    lines.append("  if (CheckpointStep > 0) {")
    if not cuda:
        lines.append("    NSLattice.NormalFullCommunicate();")
    lines.append("    flbckpt::write(RestartDir, MainLoopTimer(), NSLattice, FlagFM);")
    lines.append("  }")
    lines.append("")
    if cuda:
        lines.append("  // An unchecked launch is what hid the sm_89 / -rdc failure on")
        lines.append("  // this machine: the run printed the final banner while every")
        lines.append("  // kernel had failed.")
        lines.append("  {")
        lines.append("    cudaError_t err = cudaGetLastError();")
        lines.append("    if (err != cudaSuccess) {")
        lines.append(
            '      std::cerr << "[%s] kernel launch failed: " '
            "<< cudaGetErrorString(err) << std::endl;" % spec.name
        )
        lines.append("      return 3;")
        lines.append("    }")
        lines.append("  }")
        lines.append("")
    lines.append("#ifdef XCORE_ENABLED")
    lines.append("  xcore::MemoryPool::getInstance().print_status();")
    lines.append("#endif")
    lines.append('  Printer::Print_BigBanner(std::string("Calculation Complete!"));')
    lines.append("  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());")
    lines.append('  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));')
    lines.append("  Printer::Endl();")
    lines.append("")

    checksum_lines = _checksum_lines(spec, cuda)
    if checksum_lines:
        lines.append("  // ------------------ checksum ------------------")
        lines.extend("  " + line for line in checksum_lines)
        lines.append("")

    lines.append("  return 0;")
    lines.append("}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# ini and Makefile
# ---------------------------------------------------------------------------


def _render_ini(spec: CaseSpec) -> str:
    schedule = spec.schedule
    physics = spec.physics
    converter = physics.get("converter", "simplified")

    steps = _int(_require(schedule, "steps", "[schedule]"), "[schedule] steps")
    output_step = _int(
        _require(schedule, "output_step", "[schedule]"), "[schedule] output_step"
    )
    tol = _num(_require(schedule, "tol", "[schedule]"), "[schedule] tol")
    checkpoint_step = _int(schedule.get("checkpoint_step", 0), "[schedule] checkpoint_step")
    restart = 1 if _truthy(schedule.get("restart", 0)) else 0
    restart_dir = str(schedule.get("restart_dir", "./ckpt"))

    lines: List[str] = []
    lines.append("# generated by casegen from %s -- do not edit"
                 % os.path.basename(spec.path or "?"))
    lines.append("# spec sha256: %s" % (spec.source_hash or "unknown"))
    lines.append("")
    lines.append("[parallel]")
    lines.append("thread_num = %d" % _int(spec.thread_num, "[grid] thread_num"))
    lines.append("")
    lines.append("[Mesh]")
    for axis, name in enumerate(("Ni", "Nj", "Nk")[: spec.dim]):
        lines.append("%s = %d" % (name, _int(spec.n[axis], "[grid] n")))
    lines.append("Cell_Len = %s" % _num(spec.cell_len, "[grid] cell_len"))
    lines.append("BlockCellLen = %d" % _int(spec.block_cell_len, "[grid] block_cell_len"))
    lines.append("")
    if converter == "full":
        lines.append("[Physical_Property]")
        lines.append(
            "rho_ref = %s" % _num(_require(physics, "rho_ref", "[physics]"), "[physics] rho_ref")
        )
        lines.append(
            "Kine_Visc = %s"
            % _num(_require(physics, "kine_visc", "[physics]"), "[physics] kine_visc")
        )
        lines.append("")
    lines.append("[LB]")
    lines.append("RT = %s" % _num(_require(physics, "rt", "[physics]"), "[physics] rt"))
    if converter == "simplified":
        lines.append(
            "CharU = %s"
            % _num(_require(physics, "char_u", "[physics]"), "[physics] char_u")
        )
    lines.append("")
    if converter == "full":
        lines.append("[Init_Conditions]")
        lines.append(
            "U_Max = %s"
            % _num(_require(physics, "u_max", "[physics]"), "[physics] u_max")
        )
        lines.append("")
    lines.append("[Simulation_Settings]")
    lines.append("TotalStep = %d" % steps)
    lines.append("OutputStep = %d" % output_step)
    lines.append("Restart = %d" % restart)
    lines.append("CheckpointStep = %d" % checkpoint_step)
    lines.append("RestartDir = %s" % restart_dir)
    lines.append("")
    lines.append("[tolerance]")
    lines.append("tol = %s" % tol)
    lines.append("")
    return "\n".join(lines)


def _render_geo(spec: CaseSpec, out_dir: str) -> str:
    """Render the ``.geo`` recipe that tools/geometry_build turns into the mask.

    The recipe is a standalone format: source geometry is expressed in the
    same physical frame as ``domain``, so the values are written verbatim.
    STL paths are the only relative paths, rebased from the spec directory to
    the output directory so the generated file does not depend on the cwd.
    """
    geometry = spec.geometry
    if geometry is None or not geometry.sources:
        raise CasegenError("_render_geo called without geometry sources")
    origin = geometry.origin
    if len(origin) != 3:
        raise CasegenError("[geometry] origin needs 3 values")
    lo = [_num(value, "[geometry] origin") for value in origin]
    hi = [
        _num(origin[axis] + spec.n[axis] * spec.cell_len, "domain")
        for axis in range(3)
    ]
    spec_dir = os.path.dirname(os.path.abspath(spec.path))
    lines = [
        "# generated by casegen from %s -- do not edit"
        % os.path.basename(spec.path or "?"),
        "# spec sha256: %s" % (spec.source_hash or "unknown"),
        "",
        "[geometry]",
        "grid = %s"
        % ", ".join(str(_int(value, "[grid] n")) for value in spec.n),
        "domain = %s" % ", ".join(lo + hi),
        "cell_len = %s" % _num(spec.cell_len, "[grid] cell_len"),
        "band_width = %s" % _num(geometry.band_width, "[geometry] band_width"),
        "",
    ]
    for source in geometry.sources:
        lines.append("[source.%s]" % source.name)
        lines.append("shape = %s" % source.shape)
        lines.append("entity = %d" % source.entity)
        lines.append("material = %d" % source.material)
        lines.append("role = %s" % source.role)
        if source.shape == "box":
            lines.append(
                "box = %s"
                % ", ".join(
                    _num(value, f"[geometry.source.{source.name}] box")
                    for value in source.box or ()
                )
            )
        elif source.shape == "sphere":
            lines.append(
                "center = %s"
                % ", ".join(
                    _num(value, f"[geometry.source.{source.name}] center")
                    for value in source.center or ()
                )
            )
            lines.append(
                "radius = %s"
                % _num(source.radius, f"[geometry.source.{source.name}] radius")
            )
        else:
            stl_path = source.path or ""
            absolute = (
                stl_path
                if os.path.isabs(stl_path)
                else os.path.join(spec_dir, stl_path)
            )
            rel = os.path.relpath(os.path.abspath(absolute), os.path.abspath(out_dir))
            lines.append("path = %s" % rel.replace(os.sep, "/"))
            lines.append("closed = %d" % (1 if source.closed else 0))
        lines.append("")
    return "\n".join(lines)


def _root_from(out_dir: str, repo_root: str) -> str:
    rel = os.path.relpath(os.path.abspath(repo_root), os.path.abspath(out_dir))
    rel = rel.replace(os.sep, "/")
    return rel or "."


def _render_makefile(spec: CaseSpec, out_dir: str, repo_root: str) -> str:
    if spec.float_type not in {"float", "double"}:
        raise CasegenError(
            f"[case] float = {spec.float_type!r}: expected 'float' or 'double'"
        )
    root = _root_from(out_dir, repo_root)
    lines = [
        "# generated by casegen from %s -- do not edit"
        % os.path.basename(spec.path or "?"),
        "# spec sha256: %s" % (spec.source_hash or "unknown"),
        "",
        "TARGET = %s.exe" % spec.name,
        "",
        "ROOT := %s" % root,
        "",
        # The geometry recipe rule below would otherwise become the first
        # target and the default goal, so `make` would build the .flbgeo but
        # not the executable.
        ".DEFAULT_GOAL := all",
        "",
    ]
    # A golden case is built exactly as pinned, so it keeps -O2.  A normal
    # generated case is meant to be a real simulation: optimise it like the
    # hand-written examples.  CUDA already has -arch=native, so it only needs
    # -O3.
    if spec.selftest:
        cpu_opt = ["FLAGS += -O2"]
    else:
        cpu_opt = ["FLAGS += -O3 -march=native"]
    mpi = bool(getattr(spec, "mpi", False))
    if spec.cuda:
        # tests/run_tests.sh recognises a CUDA target by this literal; make.mk
        # keys the .cu rule and the -flto guard on the same variable.
        lines.extend(
            [
                "CXXC := nvcc",
                "",
                "FLAGS += -arch=native",
                "FLAGS += -O3" if not spec.selftest else "FLAGS += -O2",
                "FLAGS += -DFLOAT_TYPE=%s" % spec.float_type,
                "FLAGS += --expt-relaxed-constexpr",
            ]
        )
    elif mpi:
        # make.mk already keys the MPI build on these two: CXXC picks mpic++
        # (never g++) and -DMPI_ENABLED opens the MPI paths.  Setting only one
        # of them produces a binary that links but silently runs serially.
        lines.extend(
            [
                "CXXC := mpic++",
                "",
                "FLAGS += -DMPI_ENABLED",
            ]
        )
        lines.extend(cpu_opt)
        lines.append("FLAGS += -DFLOAT_TYPE=%s" % spec.float_type)
    else:
        lines.extend(cpu_opt)
        lines.append("FLAGS += -DFLOAT_TYPE=%s" % spec.float_type)
    lines.extend(
        [
            "",
            "OUTPUT += vtkoutput",
            "",
        ]
    )
    geometry = spec.geometry
    if geometry is not None and geometry.sources:
        geo_file = "%s.geo" % spec.name
        target_file = geometry.file
        lines.extend(
            [
                # The voxel mask is a build product of the recipe: rebuild it
                # before linking, and let `make clean` remove it with OUTPUT.
                "OUTPUT += %s" % target_file,
                "%s: %s" % (target_file, geo_file),
                "\t$(MAKE) -C $(ROOT)/tools/geometry_build",
                "\t$(ROOT)/tools/geometry_build/geometry_build.exe %s -o $@"
                % geo_file,
                # `all` rather than `$(TARGET)`: make.mk links with `$^`, so a
                # normal prerequisite would put the .flbgeo on the link line.
                "all: %s" % target_file,
                "",
            ]
        )
    if spec.selftest:
        # `all: regen` must precede make.mk's `all: $(TARGET)` so the drift
        # check runs before the build, not after it.
        here = os.path.dirname(os.path.abspath(__file__))
        casegen_py = os.path.normpath(os.path.join(here, "..", "casegen.py"))
        spec_rel = os.path.relpath(os.path.abspath(spec.path), os.path.abspath(out_dir))
        casegen_rel = os.path.relpath(casegen_py, os.path.abspath(out_dir))
        lines.extend(
            [
                "all: regen",
                "",
                ".PHONY: regen",
                "regen:",
                "\t@command -v python3 >/dev/null || "
                '{ echo "casegen: python3 required" >&2; exit 1; }',
                "\t@python3 %s %s -o . --check --selftest"
                % (
                    casegen_rel.replace(os.sep, "/"),
                    spec_rel.replace(os.sep, "/"),
                ),
                "",
            ]
        )
    lines.extend(["include $(ROOT)/make.mk", ""])
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# public entry point
# ---------------------------------------------------------------------------


def render_case(spec: CaseSpec, out_dir: str, repo_root: str) -> Dict[str, str]:
    """Render ``spec`` as ``{relative path: content}`` for ``out_dir``.

    ``out_dir`` and ``repo_root`` are only used to compute the Makefile's
    ``ROOT``, so ``--check`` against the real output directory sees the same
    text as a later write to that directory.
    """
    _check_supported(spec)
    rendered = {
        "%s.%s" % (spec.name, "cu" if spec.cuda else "cpp"): _render_cpp(spec),
        "%s.ini" % spec.name: _render_ini(spec),
        "Makefile": _render_makefile(spec, out_dir, repo_root),
    }
    if spec.geometry is not None and spec.geometry.sources:
        rendered["%s.geo" % spec.name] = _render_geo(spec, out_dir)
    used = sorted(
        {task.op for task in spec.tasks if task.op in spec.operator_specs}
    )
    if used:
        try:  # imported as flbcase.emit_case
            from .operator_dsl import render_operator
        except ImportError:  # pragma: no cover - top-level import
            from operator_dsl import render_operator  # type: ignore
        for name in used:
            rendered["operators/%s.hh" % name] = render_operator(
                spec.operator_specs[name]
            )
    return rendered
