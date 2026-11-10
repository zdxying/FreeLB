#!/usr/bin/env python3
"""Command line front end for the FreeLB case generator.

Usage:
    casegen.py <spec.case> -o <dir> [--check] [--force]
    casegen.py --list

``--check`` renders the case in memory and compares every generated file with
the copy in ``<dir>``; it exits 2 when anything is missing or stale, 0 when the
directory is up to date.  Without ``--check`` the generator refuses to
overwrite files whose content differs unless ``--force`` is given.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Iterable, List, Mapping, Optional, Sequence, Tuple

_HERE = Path(__file__).resolve().parent
_REPO_ROOT = _HERE.parent.parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from flbcase.model import (  # noqa: E402  (path bootstrap above)
    BOUNDARIES,
    CONVERTERS,
    DECOMPOSE,
    FIELD_TYPES,
    GEOMETRY_ACTIONS,
    LATTICE_INFO,
    OPERATORS,
    SCHEDULES,
    WRITERS,
    CasegenError,
    SourceSpan,
)
from flbcase.spec import PHASE, parse_case  # noqa: E402

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_STALE = 2


# ---------------------------------------------------------------------------
# Rendering boundary
# ---------------------------------------------------------------------------


def _load_renderer() -> Any:
    """Import the emitter lazily: ``--list`` and bad specs stay import-light."""
    from flbcase.emit_case import render_case

    return render_case


def _render(
    render_case: Any, spec: Any, out_dir: Path, repo_root: Path
) -> List[Tuple[str, Any]]:
    rendered = render_case(spec, out_dir, repo_root)
    if rendered is None:
        raise CasegenError("render_case returned no files")
    if isinstance(rendered, Mapping):
        items: Iterable[Any] = rendered.items()
    else:
        items = rendered
    result: List[Tuple[str, Any]] = []
    seen = set()
    for item in items:
        try:
            rel, content = item
        except (TypeError, ValueError):
            raise CasegenError(
                "render_case must return a mapping or (path, content) pairs"
            )
        rel = str(rel)
        if rel in seen:
            raise CasegenError(f"render_case returned '{rel}' twice")
        if not isinstance(content, (str, bytes)):
            raise CasegenError(f"render_case content for '{rel}' is not text or bytes")
        seen.add(rel)
        result.append((rel, content))
    return result


def _target_path(out_dir: Path, rel: str) -> Path:
    candidate = Path(rel)
    if candidate.is_absolute():
        try:
            candidate.relative_to(out_dir)
        except ValueError as exc:
            raise CasegenError(
                f"rendered path '{rel}' is outside the output directory"
            ) from exc
        return candidate
    if ".." in candidate.parts:
        raise CasegenError(f"rendered path '{rel}' escapes the output directory")
    return out_dir / candidate


def _as_bytes(content: Any) -> bytes:
    return content.encode("utf-8") if isinstance(content, str) else content


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="casegen.py",
        description="Generate a FreeLB case from a .case spec.",
    )
    parser.add_argument("spec", nargs="?", help="input .case file")
    parser.add_argument("-o", "--output", metavar="DIR", help="output directory")
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare rendered files with DIR and exit 2 when they differ",
    )
    parser.add_argument(
        "--force", action="store_true", help="overwrite files whose content differs"
    )
    parser.add_argument(
        "-l", "--list", action="store_true", dest="catalog", help="print the catalogue"
    )
    parser.add_argument(
        "--selftest",
        action="store_true",
        help="emit a 'regen --check' rule into the generated Makefile",
    )
    return parser


def _attach_operators(spec: Any) -> None:
    """Parse the .op files referenced by [operators] and attach their ASTs."""
    if not spec.operators:
        return
    from flbcase.operator_dsl import parse_operator

    base = Path(spec.path).resolve().parent
    for name, rel in spec.operators.items():
        op_path = (base / rel).resolve()
        if not op_path.is_file():
            raise CasegenError(
                f"[operators] {name}: cannot read '{op_path}'",
                SourceSpan(spec.path, 1),
            )
        op = parse_operator(str(op_path), q=spec.q, dim=spec.dim)
        if op.name != name:
            raise CasegenError(
                f"[operators] {name} points at a file declaring operator "
                f"'{op.name}'",
                SourceSpan(spec.path, 1),
            )
        spec.operator_specs[name] = op


def _print_catalogue() -> None:
    print("lattices:")
    for name, (dim, q) in LATTICE_INFO.items():
        print(f"  {name}  dim={dim} q={q}")
    print("fields:")
    for name, cpp in FIELD_TYPES.items():
        print(f"  {name}  {cpp}")
    print(f"collision operators (phase {PHASE}):")
    for name, meta in OPERATORS.items():
        if meta["phase"] <= PHASE:
            needs = ", ".join(sorted(meta["needs"]))
            dims = ", ".join(str(item) for item in sorted(meta["dims"]))
            print(f"  {name}  needs={needs} dims={dims}")
    print(f"boundary operators (phase {PHASE}):")
    for name, meta in BOUNDARIES.items():
        if meta["phase"] <= PHASE:
            needs = ", ".join(sorted(meta["needs"]))
            reads = "yes" if meta["reads_velocity"] else "no"
            print(f"  {name}  needs={needs} reads_velocity={reads}")
    print("geometry actions:")
    print("  " + ", ".join(sorted(GEOMETRY_ACTIONS)))
    print("schedules:")
    print("  " + ", ".join(sorted(SCHEDULES)))
    print("decompose:")
    print("  " + ", ".join(sorted(DECOMPOSE)))
    print("converters:")
    print("  " + ", ".join(sorted(CONVERTERS)))
    print("writers:")
    print("  " + ", ".join(sorted(WRITERS)))


def _check(
    spec: Any, render_case: Any, out_dir: Path, repo_root: Path
) -> int:
    rendered = _render(render_case, spec, out_dir, repo_root)
    stale: List[Tuple[str, str]] = []
    for rel, content in rendered:
        target = _target_path(out_dir, rel)
        if not target.is_file():
            stale.append((rel, "missing"))
            continue
        try:
            current = target.read_bytes()
        except OSError as exc:
            stale.append((rel, f"unreadable ({exc})"))
            continue
        if current != _as_bytes(content):
            stale.append((rel, "differs"))

    if not stale:
        print(f"casegen: {out_dir} is up to date ({len(rendered)} file(s))")
        return EXIT_OK
    for rel, reason in stale:
        print(f"casegen: {rel}: {reason}", file=sys.stderr)
    print(
        f"casegen: {len(stale)} of {len(rendered)} file(s) out of date in {out_dir}",
        file=sys.stderr,
    )
    return EXIT_STALE


def _write(
    spec: Any, render_case: Any, out_dir: Path, repo_root: Path, force: bool
) -> int:
    rendered = _render(render_case, spec, out_dir, repo_root)
    conflicts: List[str] = []
    for rel, content in rendered:
        target = _target_path(out_dir, rel)
        if target.is_file():
            try:
                current = target.read_bytes()
            except OSError:
                conflicts.append(rel)
                continue
            if current != _as_bytes(content):
                conflicts.append(rel)
    if conflicts and not force:
        for rel in conflicts:
            print(f"casegen: {rel}: exists with different content", file=sys.stderr)
        print(
            "casegen: refusing to overwrite "
            f"{len(conflicts)} file(s); use --force to replace them",
            file=sys.stderr,
        )
        return EXIT_ERROR

    written = 0
    for rel, content in rendered:
        target = _target_path(out_dir, rel)
        data = _as_bytes(content)
        if target.is_file():
            try:
                if target.read_bytes() == data:
                    continue
            except OSError:
                pass
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        written += 1
    print(f"casegen: wrote {written} file(s) to {out_dir}")
    return EXIT_OK


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.catalog:
        _print_catalogue()
        return EXIT_OK
    if args.spec is None:
        parser.error("a spec file is required unless --list is given")
    if args.output is None:
        parser.error("-o/--output is required unless --list is given")
    if args.check and args.force:
        parser.error("--check and --force are mutually exclusive")

    out_dir = Path(args.output).expanduser().resolve()
    try:
        spec = parse_case(args.spec)
        if args.selftest:
            spec.selftest = True
        _attach_operators(spec)
        render_case = _load_renderer()
        if args.check:
            return _check(spec, render_case, out_dir, _REPO_ROOT)
        return _write(spec, render_case, out_dir, _REPO_ROOT, args.force)
    except CasegenError as exc:
        print(f"casegen: {exc}", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":
    sys.exit(main())
