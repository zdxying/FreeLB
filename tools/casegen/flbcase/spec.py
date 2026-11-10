"""Parse and validate a ``.case`` file into the frozen :class:`CaseSpec`.

The dialect is a flat, ordered INI variant:

* ``[section]`` headers; ``[region.<name>]`` blocks form a geometry program
  whose order is the order of the file.
* ``key = value`` lines; values are strings, integers, floats, booleans or
  comma separated lists.
* a line whose first non-blank character is ``#`` or ``;`` is a comment.
  Inline comments are *not* supported and are rejected explicitly.
* unknown or repeated sections/keys are errors: a typo that silently drops a
  setting is more expensive than a hard failure.

Every user-facing failure is a :class:`CasegenError` carrying ``path:line``.
"""

from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .model import (
    BOUNDARIES,
    CONVERTERS,
    DECOMPOSE,
    FIELD_TYPES,
    GEOMETRY_ACTIONS,
    LATTICE_INFO,
    OPERATORS,
    SCHEDULES,
    WRITERS,
    CaseSpec,
    CasegenError,
    GeometrySource,
    GeometrySpec,
    Region,
    SourceSpan,
    Task,
)

# The catalogue this validator admits.  Phase 1 ships bgk/normal/movingwall;
# phase 3 adds the operator catalogue (bgk_force, mrt, les) and the plain /
# periodic schedules.  Sche schemes still at a later phase are rejected rather
# than silently unsupported.
PHASE = 3

# Fixed sections; [region.<name>] is the only repeatable family.
_FIXED_SECTIONS = (
    "case",
    "fields",
    "flags",
    "grid",
    "physics",
    "geometry",
    "operators",
    "tasks",
    "schedule",
    "output",
    "checksum",
)
_REQUIRED_SECTIONS = (
    "case",
    "fields",
    "flags",
    "grid",
    "physics",
    "tasks",
    "schedule",
)

_NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_HEADER_RE = re.compile(
    r"^\[([A-Za-z_][A-Za-z0-9_]*)"
    r"(?:\.([A-Za-z_][A-Za-z0-9_]*))?"
    r"(?:\.([A-Za-z_][A-Za-z0-9_]*))?\]$"
)
_INT_RE = re.compile(r"^[+-]?[0-9]+$")
_FLOAT_RE = re.compile(
    r"^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"
)
_INLINE_COMMENT_RE = re.compile(r"\s[#;]")


class _Section:
    """One parsed section: its header position and raw values in file order."""

    __slots__ = ("name", "kind", "sub", "line", "values")

    def __init__(self, name: str, kind: str, sub: Optional[str], line: int):
        self.name = name
        self.kind = kind
        self.sub = sub
        self.line = line
        # key -> (raw value, line number)
        self.values: Dict[str, Tuple[str, int]] = {}

    def key_line(self, key: str) -> int:
        entry = self.values.get(key)
        return entry[1] if entry is not None else self.line


class _View:
    """Typed accessors over one section that report errors with a line span."""

    def __init__(self, path: str, section: _Section):
        self.path = path
        self.section = section

    # -- diagnostics ------------------------------------------------------

    def _span(self, key: Optional[str] = None) -> SourceSpan:
        line = self.section.line if key is None else self.section.key_line(key)
        return SourceSpan(self.path, line)

    def error(self, message: str, key: Optional[str] = None) -> None:
        raise CasegenError(message, self._span(key))

    def label(self) -> str:
        return f"[{self.section.name}]"

    # -- raw access -------------------------------------------------------

    def has(self, key: str) -> bool:
        return key in self.section.values

    def keys(self) -> Sequence[str]:
        return tuple(self.section.values)

    def reject_unknown(self, allowed: Sequence[str], hint: str = "") -> None:
        allowed_set = set(allowed)
        for key in self.section.values:
            if key not in allowed_set:
                suffix = f" ({hint})" if hint else ""
                self.error(f"unknown key '{key}' in {self.label()}{suffix}", key)

    def _value(self, key: str, required: bool, default: Any) -> Any:
        if key in self.section.values:
            return self.section.values[key][0]
        if required:
            self.error(f"missing required key '{key}' in {self.label()}", key)
        return default

    # -- typed values -----------------------------------------------------

    def get_str(
        self, key: str, required: bool = True, default: Optional[str] = None
    ) -> Optional[str]:
        raw = self._value(key, required, default)
        if key not in self.section.values:
            return default
        if not raw:
            self.error(f"empty value for '{key}' in {self.label()}", key)
        return raw

    def get_bool(
        self, key: str, required: bool = True, default: Optional[bool] = None
    ) -> Optional[bool]:
        raw = self._value(key, required, default)
        if key not in self.section.values:
            return default
        lower = raw.lower()
        if lower in ("1", "true", "yes", "on"):
            return True
        if lower in ("0", "false", "no", "off"):
            return False
        self.error(f"invalid boolean {raw!r} for '{key}' (use 0/1 or true/false)", key)
        return None  # unreachable

    def get_int(
        self,
        key: str,
        required: bool = True,
        default: Optional[int] = None,
        minimum: Optional[int] = None,
        maximum: Optional[int] = None,
    ) -> Optional[int]:
        raw = self._value(key, required, default)
        if key not in self.section.values:
            return default
        if not _INT_RE.match(raw):
            self.error(f"invalid integer {raw!r} for '{key}'", key)
        value = int(raw, 10)
        if minimum is not None and value < minimum:
            self.error(f"'{key}' must be >= {minimum}, got {value}", key)
        if maximum is not None and value > maximum:
            self.error(f"'{key}' must be <= {maximum}, got {value}", key)
        return value

    def get_float(
        self,
        key: str,
        required: bool = True,
        default: Optional[float] = None,
        positive: bool = False,
        minimum: Optional[float] = None,
        maximum: Optional[float] = None,
    ) -> Optional[float]:
        raw = self._value(key, required, default)
        if key not in self.section.values:
            return default
        if not _FLOAT_RE.match(raw):
            self.error(f"invalid number {raw!r} for '{key}'", key)
        value = float(raw)
        if not math.isfinite(value):
            self.error(f"non-finite number {raw!r} for '{key}'", key)
        if positive and value <= 0:
            self.error(f"'{key}' must be > 0, got {value}", key)
        if minimum is not None and value < minimum:
            self.error(f"'{key}' must be >= {minimum}, got {value}", key)
        if maximum is not None and value > maximum:
            self.error(f"'{key}' must be <= {maximum}, got {value}", key)
        return value

    def _list(self, key: str, required: bool, default: Any) -> Any:
        raw = self._value(key, required, default)
        if key not in self.section.values:
            return default
        if not raw.strip():
            return []
        items = [item.strip() for item in raw.split(",")]
        if any(not item for item in items):
            self.error(f"empty list item in '{key}'", key)
        return items

    def get_str_list(self, key: str, required: bool = True, default: Any = None) -> List[str]:
        return list(self._list(key, required, default if default is not None else []))

    def get_int_list(self, key: str, required: bool = True, default: Any = None) -> List[int]:
        if key not in self.section.values:
            return list(default) if default is not None else []
        items = self._list(key, required, [] if default is None else default)
        values: List[int] = []
        for item in items:
            if not _INT_RE.match(item):
                self.error(f"invalid integer {item!r} in list '{key}'", key)
            values.append(int(item, 10))
        return values

    def get_float_list(self, key: str, required: bool = True, default: Any = None) -> List[float]:
        if key not in self.section.values:
            return list(default) if default is not None else []
        items = self._list(key, required, [] if default is None else default)
        values: List[float] = []
        for item in items:
            if not _FLOAT_RE.match(item):
                self.error(f"invalid number {item!r} in list '{key}'", key)
            value = float(item)
            if not math.isfinite(value):
                self.error(f"non-finite number {item!r} in list '{key}'", key)
            values.append(value)
        return values


# ---------------------------------------------------------------------------
# Scanning
# ---------------------------------------------------------------------------


def _scan(path: str, text: str) -> Dict[str, _Section]:
    """Split the file into sections, rejecting every structural mistake."""
    sections: Dict[str, _Section] = {}
    current: Optional[_Section] = None

    for lineno, raw_line in enumerate(text.splitlines(), 1):
        line = raw_line.strip()
        if not line or line[0] in "#;":
            continue

        if line.startswith("["):
            match = _HEADER_RE.match(line)
            if match is None:
                raise CasegenError(
                    f"malformed section header {line!r} "
                    "(expected [section], [region.<name>] or "
                    "[geometry.source.<name>])",
                    SourceSpan(path, lineno),
                )
            base, sub, sub2 = match.group(1), match.group(2), match.group(3)
            if base == "region":
                if sub is None:
                    raise CasegenError(
                        "a region section needs a name: [region.<name>]",
                        SourceSpan(path, lineno),
                    )
                if sub2 is not None:
                    raise CasegenError(
                        f"[region.{sub}.{sub2}] is not a valid section; "
                        "a region name is a single identifier",
                        SourceSpan(path, lineno),
                    )
                name, kind = f"region.{sub}", "region"
            elif base == "geometry" and sub == "source":
                if sub2 is None:
                    raise CasegenError(
                        "a geometry source section needs a name: "
                        "[geometry.source.<name>]",
                        SourceSpan(path, lineno),
                    )
                name, kind = f"geometry.source.{sub2}", "geometry_source"
                sub = sub2
            else:
                if sub is not None or sub2 is not None:
                    raise CasegenError(
                        f"[{base}{'.' + sub if sub else ''}"
                        f"{'.' + sub2 if sub2 else ''}] is not a valid section",
                        SourceSpan(path, lineno),
                    )
                if base not in _FIXED_SECTIONS:
                    known = ", ".join(f"[{item}]" for item in _FIXED_SECTIONS)
                    raise CasegenError(
                        f"unknown section [{base}] (known sections: {known}, "
                        "[region.<name>], [geometry.source.<name>])",
                        SourceSpan(path, lineno),
                    )
                name, kind = base, base

            if name in sections:
                first = sections[name].line
                raise CasegenError(
                    f"duplicate section [{name}] (first defined at line {first})",
                    SourceSpan(path, lineno),
                )
            current = _Section(name, kind, sub, lineno)
            sections[name] = current
            continue

        if "=" not in line:
            raise CasegenError(
                f"expected 'key = value' or a [section] header, got {line!r}",
                SourceSpan(path, lineno),
            )
        key, _, value = line.partition("=")
        key, value = key.strip(), value.strip()
        if not key:
            raise CasegenError("empty key", SourceSpan(path, lineno))
        if any(char.isspace() for char in key):
            raise CasegenError(
                f"invalid key {key!r} (no spaces allowed)", SourceSpan(path, lineno)
            )
        if current is None:
            raise CasegenError(
                f"key '{key}' appears outside any section", SourceSpan(path, lineno)
            )
        if _INLINE_COMMENT_RE.search(value):
            raise CasegenError(
                "inline comments are not supported; put the comment on its own line",
                SourceSpan(path, lineno),
            )
        if key in current.values:
            first = current.values[key][1]
            raise CasegenError(
                f"duplicate key '{key}' in [{current.name}] (first defined at line {first})",
                SourceSpan(path, lineno),
            )
        current.values[key] = (value, lineno)

    return sections


# ---------------------------------------------------------------------------
# Building the frozen model
# ---------------------------------------------------------------------------


def _build_case(path: str, sections: Dict[str, _Section], dim: int) -> Dict[str, Any]:
    view = _View(path, sections["case"])
    view.reject_unknown(
        ("name", "dim", "float", "latset", "cuda", "selftest", "mpi")
    )

    name = view.get_str("name")
    assert name is not None
    if not _NAME_RE.match(name):
        view.error(
            f"invalid case name {name!r} (letters, digits and '_' only)", "name"
        )

    float_type = view.get_str("float").lower()
    if float_type not in ("double", "float"):
        view.error(f"unsupported float type {float_type!r} (use double or float)", "float")

    latset = view.get_str("latset").upper()
    if latset not in LATTICE_INFO:
        known = ", ".join(LATTICE_INFO)
        view.error(f"unknown lattice set {latset!r} (supported: {known})", "latset")
    if LATTICE_INFO[latset][0] != dim:
        view.error(
            f"latset {latset} is {LATTICE_INFO[latset][0]}D but dim = {dim}", "latset"
        )

    cuda = view.get_bool("cuda", required=False, default=False)
    selftest = view.get_bool("selftest", required=False, default=False)
    mpi = view.get_bool("mpi", required=False, default=False)
    return {
        "name": name,
        "float_type": float_type,
        "latset": latset,
        "cuda": bool(cuda),
        "selftest": bool(selftest),
        "mpi": bool(mpi),
    }


def _build_fields(path: str, sections: Dict[str, _Section]) -> List[str]:
    view = _View(path, sections["fields"])
    view.reject_unknown(("list",))
    raw = view.get_str_list("list")
    fields: List[str] = []
    for item in raw:
        field = item.upper()
        if field not in FIELD_TYPES:
            known = ", ".join(FIELD_TYPES)
            view.error(f"unknown field {item!r} (supported: {known})", "list")
        if field in fields:
            view.error(f"duplicate field {field!r}", "list")
        fields.append(field)
    if not fields:
        view.error("[fields] list must not be empty", "list")
    return fields


def _build_flags(path: str, sections: Dict[str, _Section]) -> Dict[str, int]:
    view = _View(path, sections["flags"])
    if not view.keys():
        view.error("[flags] must define at least the void flag")
    flags: Dict[str, int] = {}
    seen_bits: Dict[int, str] = {}
    for name in view.keys():
        if not _NAME_RE.match(name):
            view.error(f"invalid flag name {name!r}", name)
        bit = view.get_int(name, minimum=1, maximum=255)
        assert bit is not None
        if bit & (bit - 1):
            view.error(f"flag '{name}' = {bit} is not a power of two", name)
        if bit in seen_bits:
            view.error(
                f"flag '{name}' reuses bit {bit} already assigned to "
                f"'{seen_bits[bit]}'",
                name,
            )
        seen_bits[bit] = name
        flags[name] = bit
    if "void" not in flags:
        view.error("[flags] must define 'void'", None)
    return flags


def _build_geometry_source(path: str, section: _Section) -> GeometrySource:
    view = _View(path, section)
    shape = view.get_str("shape").lower()
    if shape not in ("box", "sphere", "stl"):
        view.error(
            f"unsupported source shape {shape!r} (use box, sphere or stl)", "shape"
        )
    allowed = ["shape", "entity", "material", "role"]
    if shape == "box":
        allowed.append("box")
    elif shape == "sphere":
        allowed.extend(("center", "radius"))
    else:
        allowed.extend(("path", "closed"))
    view.reject_unknown(allowed, f"shape '{shape}' does not use it")

    entity = view.get_int("entity", minimum=0)
    material = view.get_int("material", minimum=0)
    role = view.get_str("role").lower()
    if role not in ("domain", "solid"):
        view.error(f"invalid role {role!r} (use domain or solid)", "role")
    assert entity is not None and material is not None

    source = GeometrySource(
        name=section.sub or "",
        shape=shape,
        entity=entity,
        material=material,
        role=role,
        line=section.line,
    )
    if shape == "box":
        if not view.has("box"):
            view.error("missing required key 'box'", "box")
        box = view.get_float_list("box")
        if len(box) != 6:
            view.error(f"'box' needs 6 values (lo, hi per axis), got {len(box)}", "box")
        for axis in range(3):
            if not box[axis] < box[axis + 3]:
                view.error(
                    f"'box' axis {axis} must satisfy lo < hi, got "
                    f"{box[axis]} .. {box[axis + 3]}",
                    "box",
                )
        source.box = tuple(box)
    elif shape == "sphere":
        if not view.has("center"):
            view.error("missing required key 'center'", "center")
        center = view.get_float_list("center")
        if len(center) != 3:
            view.error(f"'center' needs 3 values, got {len(center)}", "center")
        source.center = tuple(center)
        source.radius = view.get_float("radius", positive=True)
    else:
        source.path = view.get_str("path")
        assert source.path is not None
        source.closed = bool(view.get_bool("closed", required=False, default=True))
    return source


def _build_geometry_sources(
    path: str, sections: Dict[str, _Section]
) -> List[GeometrySource]:
    return [
        _build_geometry_source(path, section)
        for section in sections.values()
        if section.kind == "geometry_source"
    ]


def _build_geometry(
    path: str, sections: Dict[str, _Section], dim: int, flags: Dict[str, int]
) -> Optional[GeometrySpec]:
    source_sections = [
        section for section in sections.values() if section.kind == "geometry_source"
    ]
    if "geometry" not in sections:
        if source_sections:
            first = source_sections[0]
            raise CasegenError(
                f"[geometry.source.{first.sub}] needs a [geometry] section "
                "with 'file'",
                SourceSpan(path, first.line),
            )
        return None
    view = _View(path, sections["geometry"])
    view.reject_unknown(("file", "domain_flag", "void_flag", "band_width", "origin"))
    if dim != 3:
        view.error(
            f"[geometry] is only supported for dim = 3, got dim = {dim}", None
        )
    file = view.get_str("file")
    assert file is not None
    domain_flag = view.get_str("domain_flag")
    assert domain_flag is not None
    void_flag = view.get_str("void_flag", required=False, default="void")
    assert void_flag is not None
    for key, name in (("domain_flag", domain_flag), ("void_flag", void_flag)):
        if name not in flags:
            known = ", ".join(flags)
            view.error(
                f"[geometry] '{key}' names undefined flag {name!r} "
                f"(defined: {known})",
                key,
            )
    band_width = view.get_float(
        "band_width", required=False, default=1.5, positive=True
    )
    origin = view.get_float_list(
        "origin", required=False, default=[0.0, 0.0, 0.0]
    )
    if len(origin) != 3:
        view.error(f"'origin' needs 3 values, got {len(origin)}", "origin")
    return GeometrySpec(
        file=file,
        domain_flag=domain_flag,
        void_flag=void_flag,
        band_width=float(band_width),
        origin=tuple(float(value) for value in origin),
        sources=_build_geometry_sources(path, sections),
    )


def _build_grid(path: str, sections: Dict[str, _Section], dim: int) -> Dict[str, Any]:
    view = _View(path, sections["grid"])
    view.reject_unknown(("n", "cell_len", "block_cell_len", "thread_num", "decompose"))
    n = view.get_int_list("n")
    if len(n) != dim:
        view.error(f"'n' needs {dim} values for dim = {dim}, got {len(n)}", "n")
    if any(value < 1 for value in n):
        view.error(f"every 'n' entry must be >= 1, got {n}", "n")
    cell_len = view.get_float("cell_len", positive=True)
    block_cell_len = view.get_int("block_cell_len", minimum=1)
    thread_num = view.get_int("thread_num", required=False, default=1, minimum=1)
    decompose = view.get_str("decompose", required=False, default="single").lower()
    if decompose not in DECOMPOSE:
        view.error(f"unknown decompose {decompose!r} (use single or auto)", "decompose")
    return {
        "n": tuple(n),
        "cell_len": cell_len,
        "block_cell_len": block_cell_len,
        "thread_num": thread_num,
        "decompose": decompose,
    }


def _build_physics(path: str, sections: Dict[str, _Section], dim: int) -> Dict[str, Any]:
    view = _View(path, sections["physics"])
    converter = view.get_str("converter").lower()
    if converter not in CONVERTERS:
        view.error(f"unknown converter {converter!r} (use simplified or full)", "converter")

    if converter == "simplified":
        allowed = ("converter", "char_u", "rt", "rho_init", "init_u", "force", "smagorinsky")
        view.reject_unknown(allowed, "converter 'simplified' does not use it")
        physics: Dict[str, Any] = {
            "converter": converter,
            "char_u": view.get_float("char_u", positive=True),
            "rt": view.get_float("rt", positive=True),
        }
    else:
        allowed = (
            "converter",
            "cell_len",
            "rho_ref",
            "kine_visc",
            "u_max",
            "rt",
            "rho_init",
            "init_u",
            "force",
            "smagorinsky",
        )
        view.reject_unknown(allowed, "converter 'full' does not use it")
        physics = {
            "converter": converter,
            "cell_len": view.get_float("cell_len", positive=True),
            "rho_ref": view.get_float("rho_ref", positive=True),
            "kine_visc": view.get_float("kine_visc", positive=True),
            "u_max": view.get_float("u_max", positive=True),
            "rt": view.get_float("rt", positive=True),
        }

    physics["rho_init"] = view.get_float("rho_init", required=False, default=1.0, positive=True)
    init_u = view.get_float_list("init_u", required=False, default=[0.0] * dim)
    if len(init_u) != dim:
        view.error(f"'init_u' needs {dim} values, got {len(init_u)}", "init_u")
    physics["init_u"] = tuple(init_u)

    # `force` is a constant lattice-unit force per unit volume written into
    # the FORCE field; no conversion exists, matching how examples/dambreak2d
    # seeds CONSTFORCE.  Absent means a zero force.
    force = view.get_float_list("force", required=False, default=None)
    if force:
        if len(force) != dim:
            view.error(f"'force' needs {dim} values, got {len(force)}", "force")
        physics["force"] = tuple(force)
    else:
        physics["force"] = None

    smagorinsky = view.get_float(
        "smagorinsky", required=False, default=None, minimum=0.0
    )
    physics["smagorinsky"] = smagorinsky
    return physics


def _build_regions(
    path: str, sections: Dict[str, _Section], dim: int
) -> List[Region]:
    regions: List[Region] = []
    for section in sections.values():
        if section.kind != "region":
            continue
        view = _View(path, section)
        view.reject_unknown(("shape", "box", "action", "flag", "only_if", "value"))

        shape = view.get_str("shape").lower()
        if shape != "aabb":
            view.error(f"unsupported shape {shape!r} (only aabb is implemented)", "shape")

        box = view.get_float_list("box")
        if len(box) != 2 * dim:
            view.error(
                f"'box' needs {2 * dim} values (lo, hi per axis), got {len(box)}", "box"
            )

        action = view.get_str("action").lower()
        if action not in GEOMETRY_ACTIONS:
            known = ", ".join(sorted(GEOMETRY_ACTIONS))
            view.error(f"unknown action {action!r} (supported: {known})", "action")

        flag = view.get_str("flag")
        assert flag is not None
        only_if = view.get_str("only_if", required=False)
        if only_if is not None and not _NAME_RE.match(only_if):
            view.error(f"invalid flag name {only_if!r}", "only_if")

        value = tuple(view.get_float_list("value")) if view.has("value") else None

        if action != "paint" and only_if is not None:
            view.error("'only_if' is only valid for action = paint", "only_if")
        if action == "set_velocity":
            if value is None:
                view.error("action = set_velocity requires 'value'", "value")
            if len(value) != dim:
                view.error(f"'value' needs {dim} components, got {len(value)}", "value")
        elif value is not None:
            view.error(f"'value' is not valid for action = {action}", "value")

        regions.append(
            Region(
                name=section.sub,
                shape=shape,
                box=tuple(box),
                action=action,
                flag=flag,
                only_if=only_if,
                value=value,
                line=section.line,
            )
        )
    return regions


def _build_tasks(path: str, sections: Dict[str, _Section]) -> List[Task]:
    view = _View(path, sections["tasks"])
    tasks: List[Task] = []
    for expr in view.keys():
        components = [item.strip() for item in expr.split("|")]
        if any(not item for item in components):
            view.error(f"invalid flag expression {expr!r}", expr)
        if any(not _NAME_RE.match(item) for item in components):
            view.error(f"invalid flag name in {expr!r}", expr)
        if len(set(components)) != len(components):
            view.error(f"flag '{components[0]}' listed twice in {expr!r}", expr)
        op = view.get_str(expr)
        assert op is not None
        tasks.append(
            Task(flags=components, op=op.lower(), line=view.section.key_line(expr))
        )
    return tasks


def _build_schedule(path: str, sections: Dict[str, _Section]) -> Dict[str, Any]:
    view = _View(path, sections["schedule"])
    view.reject_unknown(
        ("type", "steps", "output_step", "tol", "checkpoint_step", "restart",
         "restart_dir", "periodic")
    )
    schedule_type = view.get_str("type").lower()
    if schedule_type == "open":
        view.error(
            "schedule 'open' is not implemented yet: it needs the anti_pressure "
            "boundary, which the generator does not emit",
            "type",
        )
    if schedule_type not in SCHEDULES:
        view.error(
            f"unknown schedule type {schedule_type!r} "
            f"(supported: {', '.join(sorted(SCHEDULES))})",
            "type",
        )
    periodic_axis = view.get_str("periodic", required=False)
    if periodic_axis is not None:
        periodic_axis = periodic_axis.lower()
        if periodic_axis not in ("x", "y", "z"):
            view.error(
                f"invalid periodic axis {periodic_axis!r} (use x, y or z)", "periodic"
            )
    steps = view.get_int("steps", minimum=1)
    output_step = view.get_int("output_step", minimum=1)
    tol = view.get_float("tol", required=False, default=0.0, minimum=0.0)
    checkpoint_step = view.get_int("checkpoint_step", required=False, default=0, minimum=0)
    restart = view.get_bool("restart", required=False, default=False)
    restart_dir = view.get_str("restart_dir", required=False, default="./ckpt")
    return {
        "type": schedule_type,
        "steps": steps,
        "output_step": output_step,
        "tol": tol,
        "checkpoint_step": checkpoint_step,
        "restart": bool(restart),
        "restart_dir": restart_dir,
        "periodic": periodic_axis,
    }


def _build_output(path: str, sections: Dict[str, _Section]) -> Dict[str, Any]:
    if "output" not in sections:
        return {"writers": [], "flag_writer": 0}
    view = _View(path, sections["output"])
    view.reject_unknown(("writers", "flag_writer"))
    writers: List[str] = []
    for item in view.get_str_list("writers", required=False, default=[]):
        writer = item.lower()
        if writer not in WRITERS:
            known = ", ".join(sorted(WRITERS))
            view.error(f"unknown writer {item!r} (supported: {known})", "writers")
        if writer in writers:
            view.error(f"duplicate writer {writer!r}", "writers")
        writers.append(writer)
    flag_writer = view.get_int("flag_writer", required=False, default=0, minimum=0, maximum=255)
    return {"writers": writers, "flag_writer": flag_writer}


def _build_checksum(path: str, sections: Dict[str, _Section]) -> Dict[str, Any]:
    if "checksum" not in sections:
        return {
            "enable": False,
            "expect_pop": None,
            "expect_rho": None,
            "expect_pop_mpi": None,
            "expect_rho_mpi": None,
            "tol": 1e-6,
        }
    view = _View(path, sections["checksum"])
    view.reject_unknown(
        ("enable", "expect_pop", "expect_rho", "expect_pop_mpi", "expect_rho_mpi", "tol")
    )
    enable = view.get_bool("enable", required=False, default=False)
    expect_pop = view.get_float("expect_pop", required=False)
    expect_rho = view.get_float("expect_rho", required=False)
    expect_pop_mpi = view.get_float("expect_pop_mpi", required=False)
    expect_rho_mpi = view.get_float("expect_rho_mpi", required=False)
    tol = view.get_float("tol", required=False, default=1e-6, positive=True)
    if (
        enable
        and expect_pop is None
        and expect_rho is None
        and expect_pop_mpi is None
        and expect_rho_mpi is None
    ):
        view.error(
            "checksum is enabled but neither expect_pop nor expect_rho is given "
            "(and no expect_pop_mpi/expect_rho_mpi)",
            "enable",
        )
    return {
        "enable": bool(enable),
        "expect_pop": expect_pop,
        "expect_rho": expect_rho,
        "expect_pop_mpi": expect_pop_mpi,
        "expect_rho_mpi": expect_rho_mpi,
        "tol": tol,
    }


def _build_operators(path: str, sections: Dict[str, _Section]) -> Dict[str, str]:
    """Map operator name -> .op path (relative to the spec file)."""
    if "operators" not in sections:
        return {}
    view = _View(path, sections["operators"])
    operators: Dict[str, str] = {}
    for name in view.keys():
        if not _NAME_RE.match(name):
            view.error(f"invalid operator name {name!r}", name)
        rel = view.get_str(name)
        if not rel:
            view.error("operator path must not be empty", name)
        operators[name] = rel
    return operators


def _build_spec(path: str, sections: Dict[str, _Section]) -> CaseSpec:
    for name in _REQUIRED_SECTIONS:
        if name not in sections:
            raise CasegenError(
                f"missing required section [{name}]", SourceSpan(path, 1)
            )

    case = _View(path, sections["case"])
    dim = case.get_int("dim")
    if dim not in (2, 3):
        case.error(f"'dim' must be 2 or 3, got {dim}", "dim")
    assert dim is not None

    case_values = _build_case(path, sections, dim)
    fields = _build_fields(path, sections)
    flags = _build_flags(path, sections)
    geometry = _build_geometry(path, sections, dim, flags)
    grid = _build_grid(path, sections, dim)
    physics = _build_physics(path, sections, dim)
    regions = _build_regions(path, sections, dim)
    tasks = _build_tasks(path, sections)
    operators = _build_operators(path, sections)
    schedule = _build_schedule(path, sections)
    output = _build_output(path, sections)
    checksum = _build_checksum(path, sections)

    spec = CaseSpec(
        name=case_values["name"],
        dim=dim,
        float_type=case_values["float_type"],
        latset=case_values["latset"],
        fields=fields,
        flags=flags,
        n=grid["n"],
        cell_len=grid["cell_len"],
        block_cell_len=grid["block_cell_len"],
        thread_num=grid["thread_num"],
        decompose=grid["decompose"],
        physics=physics,
        regions=regions,
        tasks=tasks,
        schedule=schedule,
        output=output,
        checksum=checksum,
        cuda=case_values["cuda"],
        path=path,
        operators=operators,
        selftest=case_values["selftest"],
        geometry=geometry,
    )
    # [case] mpi is a case-file switch, not part of the frozen block model;
    # CaseSpec is a plain dataclass, so it travels on the instance.
    spec.mpi = bool(case_values["mpi"])
    return spec


# ---------------------------------------------------------------------------
# Semantic validation (generation time, before any file is written)
# ---------------------------------------------------------------------------


def parse_case(path: Any) -> CaseSpec:
    """Read, parse and validate ``path``.  Raises :class:`CasegenError`."""
    spec_path = str(path)
    try:
        data = Path(path).read_bytes()
    except OSError as exc:
        raise CasegenError(f"cannot read spec '{spec_path}': {exc}") from exc
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        line = data[: exc.start].count(b"\n") + 1
        raise CasegenError(
            f"spec is not valid UTF-8: {exc.reason}", SourceSpan(spec_path, line)
        ) from exc

    sections = _scan(spec_path, text)
    spec = _build_spec(spec_path, sections)
    spec.source_hash = hashlib.sha256(data).hexdigest()
    validate_case(spec, sections)
    return spec


def _span_for(
    path: str, sections: Optional[Dict[str, _Section]], section: str, key: Optional[str] = None
) -> SourceSpan:
    line = 0
    if sections is not None and section in sections:
        sec = sections[section]
        line = sec.line if key is None else sec.key_line(key)
    return SourceSpan(path, line)


def _task_bits(spec: CaseSpec, task: Task) -> int:
    bits = 0
    for name in task.flags:
        if name not in spec.flags:
            raise CasegenError(
                f"task references undefined flag '{name}'",
                SourceSpan(spec.path, task.line),
            )
        bits |= spec.flags[name]
    return bits


def _check_box(spec: CaseSpec, region: Region, span: SourceSpan) -> None:
    """Cell-centre containment, matching AABB::isInside's closed interval."""
    eps = 1e-9
    for axis in range(spec.dim):
        lo, hi = region.box[axis], region.box[axis + spec.dim]
        size = spec.n[axis]
        if not (0.0 <= lo < hi <= float(size)):
            raise CasegenError(
                f"region '{region.name}' box axis {axis} must satisfy "
                f"0 <= {lo} < {hi} <= {size}",
                span,
            )
        first = max(0, math.ceil(lo - 0.5 - eps))
        last = min(size - 1, math.floor(hi - 0.5 + eps))
        if first > last:
            raise CasegenError(
                f"region '{region.name}' box axis {axis} ({lo}..{hi}) "
                "contains no cell centre",
                span,
            )


def _check_full_domain(spec: CaseSpec, region: Region, span: SourceSpan) -> None:
    for axis in range(spec.dim):
        lo, hi = region.box[axis], region.box[axis + spec.dim]
        if abs(lo) > 1e-9 or abs(hi - spec.n[axis]) > 1e-9:
            raise CasegenError(
                f"boundary region '{region.name}' must cover the whole domain "
                f"(axis {axis}: {lo}..{hi}, expected 0..{spec.n[axis]})",
                span,
            )


def validate_case(spec: CaseSpec, sections: Optional[Dict[str, _Section]] = None) -> None:
    """Check the generation-time validation rules against ``spec``.

    ``sections`` is the raw scan; it only supplies line numbers.  Task keys are
    sorted here (more bits first), which is the order the emitter must use.
    """
    path = spec.path
    fields = set(spec.fields)
    dim = spec.dim

    # 1. lattice/dimension, flag bits ------------------------------------
    if spec.latset not in LATTICE_INFO:
        raise CasegenError(
            f"unknown lattice set {spec.latset!r}", _span_for(path, sections, "case", "latset")
        )
    if LATTICE_INFO[spec.latset][0] != dim:
        raise CasegenError(
            f"latset {spec.latset} is {LATTICE_INFO[spec.latset][0]}D but dim = {dim}",
            _span_for(path, sections, "case", "latset"),
        )
    if not set(spec.fields) <= set(FIELD_TYPES):
        unknown = sorted(set(spec.fields) - set(FIELD_TYPES))
        raise CasegenError(
            f"unknown field(s): {', '.join(unknown)}",
            _span_for(path, sections, "fields", "list"),
        )
    if spec.thread_num != 1:
        message = "OpenMP is not wired; set thread_num = 1"
        if spec.schedule.get("periodic") is not None:
            message += "; the periodic boundary is serial-only"
        raise CasegenError(
            message,
            _span_for(path, sections, "grid", "thread_num"),
        )
    if getattr(spec, "mpi", False) and spec.cuda:
        raise CasegenError(
            "mpi = 1 and cuda = 1 select different compilers (mpic++ / nvcc); "
            "generate a CPU case for MPI",
            _span_for(path, sections, "case", "mpi"),
        )
    seen_bits: Dict[int, str] = {}
    for name, bit in spec.flags.items():
        if not (isinstance(bit, int) and 0 < bit <= 255 and not (bit & (bit - 1))):
            raise CasegenError(
                f"flag '{name}' = {bit} must be a power of two in 1..128",
                _span_for(path, sections, "flags", name),
            )
        if bit in seen_bits:
            raise CasegenError(
                f"flag '{name}' reuses bit {bit} already assigned to '{seen_bits[bit]}'",
                _span_for(path, sections, "flags", name),
            )
        seen_bits[bit] = name
    if "void" not in spec.flags:
        raise CasegenError(
            "the [flags] section must define 'void'", _span_for(path, sections, "flags")
        )
    if spec.geometry is not None:
        if dim != 3:
            raise CasegenError(
                f"[geometry] is only supported for dim = 3, got dim = {dim}",
                _span_for(path, sections, "geometry"),
            )
        for key in ("domain_flag", "void_flag"):
            name = getattr(spec.geometry, key)
            if name not in spec.flags:
                raise CasegenError(
                    f"[geometry] '{key}' names undefined flag '{name}'",
                    _span_for(path, sections, "geometry", key),
                )

    # 2-3. catalogue, field closure, task-key overlap ---------------------
    catalogue: Dict[str, Tuple[str, Dict[str, Any]]] = {}
    for name, meta in OPERATORS.items():
        if meta["phase"] <= PHASE:
            catalogue[name] = ("collision", meta)
    for name, meta in BOUNDARIES.items():
        if meta["phase"] <= PHASE:
            catalogue[name] = ("boundary", meta)

    if not spec.tasks:
        raise CasegenError(
            "the [tasks] section must map at least one flag to an operator",
            _span_for(path, sections, "tasks"),
        )

    seen_keys: List[Tuple[int, Task]] = []
    for task in spec.tasks:
        span = SourceSpan(path, task.line)
        bits = _task_bits(spec, task)
        if bits & spec.flags.get("void", 0):
            raise CasegenError(
                "a task cannot include the void flag (void cells do not collide)", span
            )
        if task.op == "anti_pressure":
            raise CasegenError(
                "boundary operator 'anti_pressure' is not implemented yet; it is "
                "planned together with the 'open' schedule",
                span,
            )
        if task.op in spec.operators:
            # Custom .op operators: the field closure lives in the parsed
            # operator (op.needs), which casegen.py attaches after validation.
            # Checked before the catalogue so a custom name that also exists
            # there is not rejected by its dims/needs metadata.
            pass
        elif task.op in catalogue:
            kind, meta = catalogue[task.op]
            missing = sorted(meta["needs"] - fields)
            if missing:
                raise CasegenError(
                    f"operator '{task.op}' needs field(s) {', '.join(missing)}, "
                    "which [fields] does not list",
                    span,
                )
            if kind == "collision" and dim not in meta.get("dims", {dim}):
                raise CasegenError(
                    f"operator '{task.op}' does not support dim = {dim}", span
                )
        else:
            allowed = ", ".join(sorted(catalogue))
            extra = (
                "; custom: " + ", ".join(sorted(spec.operators))
                if spec.operators
                else ""
            )
            raise CasegenError(
                f"unknown operator '{task.op}'; phase {PHASE} catalogue: "
                f"{allowed}{extra}",
                span,
            )
        for other_bits, other in seen_keys:
            if bits & other_bits:
                raise CasegenError(
                    f"task '{'|'.join(task.flags)}' overlaps task "
                    f"'{'|'.join(other.flags)}'; a flag may belong to one task only",
                    span,
                )
        seen_keys.append((bits, task))

    # More bits first: a composite key must be tested before its parts.
    spec.tasks.sort(key=lambda item: bin(_task_bits(spec, item)).count("1"), reverse=True)

    # Region references and ordering --------------------------------------
    assigned: set = set()
    for region in spec.regions:
        span = _span_for(path, sections, f"region.{region.name}")
        if region.flag not in spec.flags:
            raise CasegenError(
                f"region '{region.name}' uses undefined flag '{region.flag}'", span
            )
        if region.only_if is not None:
            if region.only_if not in spec.flags:
                raise CasegenError(
                    f"region '{region.name}' only_if uses undefined flag "
                    f"'{region.only_if}'",
                    span,
                )
            if region.only_if not in assigned:
                raise CasegenError(
                    f"region '{region.name}' only_if flag '{region.only_if}' has not "
                    "been written by an earlier region",
                    span,
                )
        if region.action == "boundary" and not assigned:
            raise CasegenError(
                f"boundary region '{region.name}' must follow at least one paint region",
                span,
            )
        _check_box(spec, region, _span_for(path, sections, f"region.{region.name}", "box"))
        if region.action == "boundary":
            _check_full_domain(
                spec, region, _span_for(path, sections, f"region.{region.name}", "box")
            )
        if region.action in ("paint", "boundary"):
            assigned.add(region.flag)

    # 2. action -> task mapping, 4. schedule/boundary interaction ---------
    schedule_type = spec.schedule.get("type")
    if schedule_type not in SCHEDULES:
        raise CasegenError(
            f"unknown schedule type {schedule_type!r}",
            _span_for(path, sections, "schedule", "type"),
        )
    if schedule_type == "open":
        raise CasegenError(
            "schedule 'open' is not implemented yet: it needs the anti_pressure "
            "boundary, which the generator does not emit",
            _span_for(path, sections, "schedule", "type"),
        )

    periodic_axis = spec.schedule.get("periodic")
    if periodic_axis is not None:
        if schedule_type != "plain":
            raise CasegenError(
                "'periodic' is only valid with [schedule] type = plain; the "
                "hidden schedule has no BlockBoundaryManager to combine it with",
                _span_for(path, sections, "schedule", "periodic"),
            )
        axis = "xyz".index(periodic_axis)
        if axis >= dim:
            raise CasegenError(
                f"periodic axis {periodic_axis!r} does not exist in a {dim}D case",
                _span_for(path, sections, "schedule", "periodic"),
            )
        if spec.thread_num != 1:
            raise CasegenError(
                "the periodic boundary is serial-only: set [grid] thread_num = 1",
                _span_for(path, sections, "grid", "thread_num"),
            )
        if spec.decompose != "single":
            raise CasegenError(
                "the periodic boundary does not support block decomposition: set "
                "[grid] decompose = single",
                _span_for(path, sections, "grid", "decompose"),
            )
        if "periodic" not in spec.flags:
            raise CasegenError(
                "the periodic boundary needs a flag named 'periodic' for the faces",
                _span_for(path, sections, "flags"),
            )
        for other in range(dim):
            if other != axis and spec.n[other] < 3:
                raise CasegenError(
                    f"periodic axis {periodic_axis}: grid axis {other} has "
                    f"{spec.n[other]} cells; the two face strips need at least 3",
                    _span_for(path, sections, "grid", "n"),
                )

    if spec.physics.get("force") is not None and "FORCE" not in fields:
        raise CasegenError(
            "[physics] force writes the FORCE field; add FORCE to [fields]",
            _span_for(path, sections, "physics", "force"),
        )
    if (
        spec.physics.get("smagorinsky") is not None
        and "SMAGORINSKY" not in fields
    ):
        raise CasegenError(
            "[physics] smagorinsky writes the SMAGORINSKY field; add SMAGORINSKY "
            "to [fields]",
            _span_for(path, sections, "physics", "smagorinsky"),
        )

    for region in spec.regions:
        span = _span_for(path, sections, f"region.{region.name}")
        task = spec.task_for(region.flag)
        if region.action == "paint":
            if region.flag != spec.void_name and task is None:
                raise CasegenError(
                    f"paint region '{region.name}' uses flag '{region.flag}' with no "
                    "task; TaskSelector would silently skip those cells",
                    span,
                )
        elif region.action == "boundary":
            if task is None:
                raise CasegenError(
                    f"boundary region '{region.name}' flag '{region.flag}' has no "
                    "boundary scheme; the bounce-back would be silently skipped",
                    span,
                )
            if task.op not in BOUNDARIES:
                raise CasegenError(
                    f"boundary region '{region.name}' flag '{region.flag}' maps to "
                    f"'{task.op}', which is not a boundary operator",
                    span,
                )
        elif region.action == "set_velocity":
            if "VELOCITY" not in fields:
                raise CasegenError(
                    f"set_velocity region '{region.name}' needs VELOCITY in [fields]",
                    span,
                )
            if task is None:
                raise CasegenError(
                    f"set_velocity region '{region.name}' flag '{region.flag}' has no "
                    "task",
                    span,
                )
            meta = BOUNDARIES.get(task.op)
            if meta is None or not meta.get("reads_velocity"):
                raise CasegenError(
                    f"set_velocity region '{region.name}' flag '{region.flag}' maps to "
                    f"'{task.op}', which does not read VELOCITY",
                    span,
                )

    # Output writers need their fields ------------------------------------
    writers = set(spec.output.get("writers", []))
    if "rho" in writers and "RHO" not in fields:
        raise CasegenError(
            "writer 'rho' needs RHO in [fields]",
            _span_for(path, sections, "output", "writers"),
        )
    if "velocity" in writers and "VELOCITY" not in fields:
        raise CasegenError(
            "writer 'velocity' needs VELOCITY in [fields]",
            _span_for(path, sections, "output", "writers"),
        )

    # Checksum thresholds --------------------------------------------------
    if spec.checksum.get("enable"):
        expectations = (
            spec.checksum.get("expect_pop"),
            spec.checksum.get("expect_rho"),
            spec.checksum.get("expect_pop_mpi"),
            spec.checksum.get("expect_rho_mpi"),
        )
        if all(value is None for value in expectations):
            raise CasegenError(
                "checksum is enabled but neither expect_pop nor expect_rho is "
                "given (and no expect_pop_mpi/expect_rho_mpi)",
                _span_for(path, sections, "checksum", "enable"),
            )
