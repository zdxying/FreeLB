"""Frozen data model shared by the parser, validator and emitters.

This module is the interface contract for tools/casegen.  It contains no I/O
and no code generation: the parser fills these dataclasses, the validator
rejects inconsistent combinations, and the emitters render C++/ini/Makefile.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Catalogue
# ---------------------------------------------------------------------------

# latset -> (dimension, q)
LATTICE_INFO = {
    "D2Q9": (2, 9),
    "D3Q19": (3, 19),
    "D3Q27": (3, 27),
}

# field name -> C++ type expression (T and LatSet are in scope in the case)
FIELD_TYPES = {
    "RHO": "RHO<T>",
    "VELOCITY": "VELOCITY<T, LatSet::d>",
    "POP": "POP<T, LatSet::q>",
    "FORCE": "FORCE<T, LatSet::d>",
    "SMAGORINSKY": "SMAGORINSKY<T>",
}

# operator name -> C++ dynamics type; needs is the field closure, bc marks a
# boundary task (applied to wall flags), dims restricts the lattice dimension.
OPERATORS = {
    "bgk": {
        "cpp": "collision::BGK<moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>",
        "needs": {"POP"},
        "bc": False,
        "dims": {2, 3},
        "phase": 1,
    },
    "bgk_force": {
        "cpp": (
            "collision::BGKForce<"
            "moment::forcerhoU<CELL, force::Force<CELL>>, "
            "equilibrium::SecondOrder<CELL>, force::Force<CELL>>"
        ),
        "needs": {"POP", "FORCE"},
        "bc": False,
        "dims": {2, 3},
        "phase": 3,
    },
    "mrt": {
        "cpp": "collision::MRT_Feq_RhoU<CELL>",
        "needs": {"POP"},
        "bc": False,
        "dims": {2},
        "phase": 3,
    },
    "les": {
        "cpp": (
            "collision::SmagorinskyBGK<"
            "moment::rhoU<CELL>, equilibrium::SecondOrder<CELL>>"
        ),
        "needs": {"POP", "RHO", "SMAGORINSKY"},
        "bc": False,
        "dims": {2, 3},
        "phase": 3,
    },
}

# boundary name -> C++ bounce-back scheme; reads_velocity matters for
# set_velocity actions, as_bc marks a scheme usable as a collision task.
BOUNDARIES = {
    "normal": {
        "cpp": "bounceback::normal<CELL>",
        "needs": {"POP"},
        "reads_velocity": False,
        "phase": 1,
    },
    "movingwall": {
        "cpp": "bounceback::movingwall<CELL>",
        "needs": {"POP", "VELOCITY"},
        "reads_velocity": True,
        "phase": 1,
    },
    "anti_pressure": {
        "cpp": "bounceback::anti_pressure<CELL>",
        "needs": {"POP", "VELOCITY"},
        "reads_velocity": True,
        # The open schedule that would drive this scheme is not implemented
        # yet; keep it out of the advertised catalogue until it is.
        "phase": 4,
    },
}

SCHEDULES = {"hidden", "plain"}
CONVERTERS = {"simplified", "full"}
DECOMPOSE = {"single", "auto"}
GEOMETRY_ACTIONS = {"paint", "boundary", "set_velocity"}
WRITERS = {"rho", "velocity", "flag"}


# ---------------------------------------------------------------------------
# Parsed case
# ---------------------------------------------------------------------------


@dataclass
class Region:
    name: str
    shape: str                       # only "aabb" in Phase 1-3
    box: Tuple[float, ...]           # AABB endpoints in cell units
    action: str                      # paint | boundary | set_velocity
    flag: Optional[str] = None
    only_if: Optional[str] = None
    value: Optional[Tuple[float, ...]] = None
    line: int = 0


@dataclass
class Task:
    flags: List[str]                 # one or more flag names, OR-ed together
    op: str                          # operator or boundary name
    line: int = 0


@dataclass
class GeometrySource:
    """One primitive in a ``[geometry.source.<name>]`` recipe."""

    name: str
    shape: str                       # box | sphere | stl
    entity: int
    material: int
    role: str                        # domain | solid
    box: Optional[Tuple[float, ...]] = None
    center: Optional[Tuple[float, ...]] = None
    radius: Optional[float] = None
    path: Optional[str] = None       # STL path as written in the spec
    closed: bool = True
    line: int = 0


@dataclass
class GeometrySpec:
    """Optional voxel-geometry file that supplies the domain/void mask."""

    file: str                        # path as written, relative to the run dir
    domain_flag: str
    void_flag: str = "void"
    band_width: float = 1.5
    origin: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    sources: List[GeometrySource] = field(default_factory=list)


@dataclass
class CaseSpec:
    name: str
    dim: int
    float_type: str
    latset: str
    fields: List[str]
    flags: Dict[str, int]
    n: Tuple[int, ...]
    cell_len: float
    block_cell_len: int
    thread_num: int
    decompose: str
    physics: Dict[str, Any]
    regions: List[Region]
    tasks: List[Task]
    schedule: Dict[str, Any]
    output: Dict[str, Any]
    checksum: Dict[str, Any]
    cuda: bool
    path: str
    source_hash: str = ""
    operators: Dict[str, str] = field(default_factory=dict)
    operator_specs: Dict[str, Any] = field(default_factory=dict)
    selftest: bool = False
    geometry: Optional[GeometrySpec] = None

    # ---- convenience -----------------------------------------------------

    @property
    def q(self) -> int:
        return LATTICE_INFO[self.latset][1]

    @property
    def void_name(self) -> Optional[str]:
        for name, bit in self.flags.items():
            if name == "void":
                return name
        return None

    @property
    def void_flag(self) -> int:
        name = self.void_name
        return self.flags[name] if name else 1

    def flag_bit(self, name: str) -> int:
        return self.flags[name]

    def task_for(self, flag_name: str) -> Optional[Task]:
        for task in self.tasks:
            if flag_name in task.flags:
                return task
        return None


@dataclass
class OperatorSpec:
    """A parsed .op file.  ``body`` is the DSL AST owned by operator_dsl."""

    name: str
    path: str
    variables: Dict[str, str]        # name -> scalar | vector | array
    body: Any
    needs: set                      # {"POP"} plus every field read via field()/force()
    source_hash: str = ""


@dataclass
class SourceSpan:
    path: str
    line: int


class CasegenError(Exception):
    """A user-facing spec error; the message already carries path:line."""

    def __init__(self, message: str, span: Optional[SourceSpan] = None):
        self.span = span
        if span is not None:
            message = f"{span.path}:{span.line}: {message}"
        super().__init__(message)
