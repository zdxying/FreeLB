"""Tests for the .case parser/validator and the casegen CLI.

Run with:
    python3 -m unittest discover -s tools/casegen/tests -p 'test_spec.py'
"""

from __future__ import annotations

import io
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

_HERE = Path(__file__).resolve().parent
_TOOL = _HERE.parent
if str(_TOOL) not in sys.path:
    sys.path.insert(0, str(_TOOL))

import casegen  # noqa: E402
from flbcase.model import CasegenError  # noqa: E402
from flbcase.spec import parse_case  # noqa: E402

BASE = """\
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
rho_init = 1
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
checkpoint_step = 0
restart = 0
restart_dir = ./ckpt

[output]
writers = rho, velocity
flag_writer = 0

[checksum]
enable = 1
expect_pop = 7728.702031
expect_rho = 1155.971535
tol = 1e-6
"""


class CaseFileTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = Path(self._tmp.name)

    def parse_text(self, text: str) -> object:
        path = self.dir / "test.case"
        path.write_text(text, encoding="utf-8")
        return parse_case(path)

    def mutate(self, old: str, new: str) -> str:
        self.assertIn(old, BASE)
        return BASE.replace(old, new, 1)

    def assert_error(self, text: str, needle: str, line: int | None = None) -> None:
        with self.assertRaises(CasegenError) as caught:
            self.parse_text(text)
        message = str(caught.exception)
        self.assertIn(needle, message)
        self.assertIn(str(self.dir / "test.case"), message)
        if line is not None:
            self.assertIn(f":{line}:", message)


class ValidSpecTest(CaseFileTest):
    def test_parses_the_cavity_spec(self) -> None:
        spec = self.parse_text(BASE)
        self.assertEqual(spec.name, "cavity2d")
        self.assertEqual((spec.dim, spec.float_type, spec.latset), (2, "double", "D2Q9"))
        self.assertEqual(spec.fields, ["RHO", "VELOCITY", "POP"])
        self.assertEqual(
            spec.flags, {"void": 1, "bulk": 2, "wall": 4, "moving_wall": 8}
        )
        self.assertEqual(spec.n, (32, 32))
        self.assertEqual(spec.cell_len, 1.0)
        self.assertEqual(spec.block_cell_len, 10)
        self.assertEqual(spec.thread_num, 1)
        self.assertEqual(spec.decompose, "single")
        self.assertEqual(spec.physics["converter"], "simplified")
        self.assertAlmostEqual(spec.physics["rt"], 0.56)
        self.assertEqual(spec.physics["init_u"], (0.0, 0.0))
        self.assertEqual([region.name for region in spec.regions],
                         ["cavity", "walls", "lid", "lid_velocity"])
        self.assertEqual(spec.regions[2].only_if, "wall")
        self.assertEqual(spec.regions[3].value, (0.1, 0.0))
        self.assertEqual([task.op for task in spec.tasks],
                         ["bgk", "normal", "movingwall"])
        self.assertEqual(spec.schedule["type"], "hidden")
        self.assertEqual(spec.schedule["steps"], 1000)
        self.assertFalse(spec.schedule["restart"])
        self.assertEqual(spec.output, {"writers": ["rho", "velocity"], "flag_writer": 0})
        self.assertTrue(spec.checksum["enable"])
        self.assertFalse(spec.cuda)
        self.assertEqual(len(spec.source_hash), 64)
        self.assertEqual(spec.path, str(self.dir / "test.case"))

    def test_comments_and_blank_lines_are_ignored(self) -> None:
        text = "# leading comment\n; another\n\n" + BASE.replace(
            "thread_num = 1", "  ; indent comment\nthread_num = 1"
        )
        spec = self.parse_text(text)
        self.assertEqual(spec.thread_num, 1)

    def test_tasks_with_more_flag_bits_come_first(self) -> None:
        text = base_replace_tasks(
            BASE,
            "wall|moving_wall = movingwall\nbulk = bgk",
        )
        spec = self.parse_text(text)
        self.assertEqual([task.flags for task in spec.tasks],
                         [["wall", "moving_wall"], ["bulk"]])

    def test_full_converter(self) -> None:
        text = BASE.replace(
            """converter = simplified
char_u = 0.1
rt = 0.56
rho_init = 1
init_u = 0, 0""",
            """converter = full
cell_len = 0.1
rho_ref = 0.001
kine_visc = 0.02
u_max = 0.1
rt = 0.56
rho_init = 1
init_u = 0, 0""",
        )
        spec = self.parse_text(text)
        self.assertEqual(spec.physics["converter"], "full")
        self.assertAlmostEqual(spec.physics["kine_visc"], 0.02)

    def test_output_and_checksum_are_optional(self) -> None:
        text = BASE.split("[output]")[0]
        spec = self.parse_text(text)
        self.assertEqual(spec.output["writers"], [])
        self.assertFalse(spec.checksum["enable"])


class InvalidSpecTest(CaseFileTest):
    def test_unknown_section(self) -> None:
        self.assert_error(BASE + "\n[bogus]\nx = 1\n", "unknown section [bogus]")

    def test_malformed_header(self) -> None:
        self.assert_error(self.mutate("[grid]", "[grid"), "malformed section header")

    def test_duplicate_section(self) -> None:
        self.assert_error(
            self.mutate("[grid]", "[grid]\n[grid]"), "duplicate section [grid]"
        )

    def test_duplicate_region_section(self) -> None:
        self.assert_error(BASE + "\n[region.cavity]\n", "duplicate section [region.cavity]")

    def test_unknown_key(self) -> None:
        self.assert_error(
            self.mutate("thread_num = 1", "thread_num = 1\nbogus = 2"),
            "unknown key 'bogus'",
        )

    def test_duplicate_key(self) -> None:
        self.assert_error(
            self.mutate("thread_num = 1", "thread_num = 1\nthread_num = 2"),
            "duplicate key 'thread_num'",
        )

    def test_invalid_integer(self) -> None:
        self.assert_error(
            self.mutate("thread_num = 1", "thread_num = two"), "invalid integer"
        )

    def test_invalid_float(self) -> None:
        self.assert_error(self.mutate("rt = 0.56", "rt = abc"), "invalid number")

    def test_inline_comment_is_rejected(self) -> None:
        self.assert_error(
            self.mutate("thread_num = 1", "thread_num = 1 ; one thread"),
            "inline comments are not supported",
        )

    def test_missing_required_key(self) -> None:
        self.assert_error(self.mutate("latset = D2Q9\n", ""), "missing required key 'latset'")

    def test_latset_dimension_mismatch(self) -> None:
        self.assert_error(
            self.mutate("latset = D2Q9", "latset = D3Q19"), "latset D3Q19 is 3D"
        )

    def test_unknown_latset(self) -> None:
        self.assert_error(self.mutate("latset = D2Q9", "latset = D2Q5"), "unknown lattice")

    def test_unknown_field(self) -> None:
        self.assert_error(
            self.mutate("list = RHO, VELOCITY, POP", "list = RHO, VELOCITY, POP, BOGUS"),
            "unknown field",
        )

    def test_unknown_writer(self) -> None:
        self.assert_error(
            self.mutate("writers = rho, velocity", "writers = rho, pressure"),
            "unknown writer",
        )

    def test_flag_not_power_of_two(self) -> None:
        self.assert_error(self.mutate("wall = 4", "wall = 6"), "not a power of two")

    def test_duplicate_flag_bit(self) -> None:
        self.assert_error(
            self.mutate("moving_wall = 8", "moving_wall = 4"), "reuses bit 4"
        )

    def test_missing_void_flag(self) -> None:
        self.assert_error(self.mutate("void = 1\n", ""), "must define 'void'")

    def test_paint_flag_without_task(self) -> None:
        self.assert_error(self.mutate("bulk = bgk\n", ""), "with no task")

    def test_boundary_box_must_cover_domain(self) -> None:
        text = self.mutate(
            "[region.walls]\nshape = aabb\nbox = 0, 0, 32, 32",
            "[region.walls]\nshape = aabb\nbox = 0, 0, 32, 31",
        )
        self.assert_error(text, "must cover the whole domain")

    def test_box_outside_domain(self) -> None:
        text = self.mutate(
            "[region.cavity]\nshape = aabb\nbox = 0, 0, 32, 32",
            "[region.cavity]\nshape = aabb\nbox = 0, 0, 33, 32",
        )
        self.assert_error(text, "must satisfy")

    def test_box_without_cell_centre(self) -> None:
        text = self.mutate("box = 1, 31, 31, 32", "box = 1.1, 31, 1.4, 32")
        self.assert_error(text, "contains no cell centre")

    def test_boundary_region_before_paint(self) -> None:
        text = self.mutate("action = paint\nflag = bulk", "action = boundary\nflag = bulk")
        self.assert_error(text, "must follow at least one paint region")

    def test_only_if_uses_undefined_flag(self) -> None:
        self.assert_error(
            self.mutate("only_if = wall", "only_if = ghost"), "undefined flag 'ghost'"
        )

    def test_only_if_flag_not_written_yet(self) -> None:
        self.assert_error(
            self.mutate("only_if = wall", "only_if = moving_wall"),
            "has not been written by an earlier region",
        )

    def test_only_if_on_set_velocity(self) -> None:
        text = self.mutate(
            "action = set_velocity\nflag = moving_wall\nvalue = 0.1, 0",
            "action = set_velocity\nflag = moving_wall\nonly_if = wall\nvalue = 0.1, 0",
        )
        self.assert_error(text, "only valid for action = paint")

    def test_set_velocity_requires_velocity_field(self) -> None:
        text = self.mutate("list = RHO, VELOCITY, POP", "list = RHO, POP")
        text = text.replace("moving_wall = movingwall", "moving_wall = normal")
        self.assert_error(text, "needs VELOCITY")

    def test_set_velocity_requires_velocity_reading_boundary(self) -> None:
        self.assert_error(
            self.mutate("moving_wall = movingwall", "moving_wall = normal"),
            "does not read VELOCITY",
        )

    def test_hidden_boundary_flag_needs_boundary_task(self) -> None:
        self.assert_error(
            self.mutate("wall = normal", "wall = bgk"), "not a boundary operator"
        )

    def test_plain_schedule_accepts_boundary_scheme(self) -> None:
        spec = self.parse_text(self.mutate("type = hidden", "type = plain"))
        self.assertEqual(spec.schedule["type"], "plain")
        # The boundary flag keeps its scheme; the emitter turns it into a
        # BlockBoundaryManager instead of a collision task.
        self.assertEqual(spec.task_for("wall").op, "normal")

    def test_task_references_undefined_flag(self) -> None:
        self.assert_error(
            self.mutate("bulk = bgk", "bogus = bgk"), "undefined flag 'bogus'"
        )

    def test_overlapping_task_keys(self) -> None:
        text = self.mutate("wall = normal\n", "bulk|wall = bgk\n")
        self.assert_error(text, "overlaps task")

    def test_task_cannot_include_void(self) -> None:
        self.assert_error(
            self.mutate("bulk = bgk", "bulk|void = bgk"), "cannot include the void flag"
        )

    def test_unknown_operator(self) -> None:
        self.assert_error(
            self.mutate("bulk = bgk", "bulk = bkg_typo"), "unknown operator 'bkg_typo'"
        )

    def test_les_needs_smagorinsky_field(self) -> None:
        self.assert_error(
            self.mutate("bulk = bgk", "bulk = les"), "needs field(s) SMAGORINSKY"
        )

    def test_mrt_is_d2q9_only(self) -> None:
        text = self.mutate("dim = 2", "dim = 3").replace(
            "latset = D2Q9", "latset = D3Q19"
        )
        text = text.replace("n = 32, 32", "n = 8, 8, 8")
        text = text.replace("init_u = 0, 0", "init_u = 0, 0, 0")
        text = text.replace("box = 0, 0, 32, 32", "box = 0, 0, 0, 8, 8, 8")
        text = text.replace("box = 1, 31, 31, 32", "box = 1, 1, 7, 7, 7, 8")
        text = text.replace("value = 0.1, 0", "value = 0.1, 0, 0")
        text = text.replace("bulk = bgk", "bulk = mrt")
        self.assert_error(text, "does not support dim = 3")

    def test_force_requires_force_field(self) -> None:
        self.assert_error(
            self.mutate("init_u = 0, 0", "init_u = 0, 0\nforce = 0, -1e-5"),
            "add FORCE to [fields]",
        )

    def test_smagorinsky_requires_its_field(self) -> None:
        self.assert_error(
            self.mutate("init_u = 0, 0", "init_u = 0, 0\nsmagorinsky = 0.1"),
            "add SMAGORINSKY to [fields]",
        )

    def test_open_schedule_is_rejected(self) -> None:
        self.assert_error(
            self.mutate("type = hidden", "type = open"),
            "schedule 'open' is not implemented",
        )

    def test_anti_pressure_is_rejected(self) -> None:
        self.assert_error(
            self.mutate("wall = normal", "wall = anti_pressure"),
            "anti_pressure' is not implemented",
        )

    def test_periodic_axis_parses_with_plain_schedule(self) -> None:
        text = base_replace_tasks(
            BASE, "bulk|periodic = bgk\nwall = normal\nmoving_wall = movingwall"
        )
        text = text.replace("moving_wall = 8", "moving_wall = 8\nperiodic = 16")
        text = text.replace("type = hidden", "type = plain\nperiodic = x")
        text = text.replace(
            "[region.lid]\nshape = aabb\nbox = 1, 31, 31, 32\naction = paint\n"
            "flag = moving_wall\nonly_if = wall",
            "[region.lid]\nshape = aabb\nbox = 1, 31, 31, 32\naction = paint\n"
            "flag = moving_wall\nonly_if = wall\n\n"
            "[region.sides]\nshape = aabb\nbox = 0, 1, 1, 31\n"
            "action = paint\nflag = periodic\nonly_if = wall\n\n"
            "[region.sides_right]\nshape = aabb\nbox = 31, 1, 32, 31\n"
            "action = paint\nflag = periodic\nonly_if = wall",
        )
        spec = self.parse_text(text)
        self.assertEqual(spec.schedule["type"], "plain")
        self.assertEqual(spec.schedule["periodic"], "x")

    def test_periodic_needs_the_plain_schedule(self) -> None:
        self.assert_error(
            self.mutate("type = hidden", "type = hidden\nperiodic = x"),
            "only valid with [schedule] type = plain",
        )

    def test_periodic_axis_must_match_the_grid(self) -> None:
        self.assert_error(
            self.mutate("type = hidden", "type = plain\nperiodic = z"),
            "does not exist in a 2D case",
        )

    def test_periodic_boundary_is_serial(self) -> None:
        text = self.mutate("type = hidden", "type = plain\nperiodic = x")
        text = text.replace("thread_num = 1", "thread_num = 2")
        self.assert_error(text, "serial-only")

    def test_task_needs_missing_field(self) -> None:
        text = self.mutate("list = RHO, VELOCITY, POP", "list = RHO, VELOCITY")
        self.assert_error(text, "needs field(s) POP")

    def test_physics_converter_key_mismatch(self) -> None:
        self.assert_error(
            self.mutate("converter = simplified", "converter = full"),
            "does not use it",
        )

    def test_missing_required_section(self) -> None:
        block = BASE[BASE.index("[grid]"):BASE.index("[physics]")]
        self.assert_error(BASE.replace(block, "", 1), "missing required section [grid]")

    def test_checksum_without_expectation(self) -> None:
        text = self.mutate(
            "expect_pop = 7728.702031\nexpect_rho = 1155.971535\n", ""
        )
        self.assert_error(text, "neither expect_pop nor expect_rho")

    def test_cuda_flag_is_parsed(self) -> None:
        spec = self.parse_text(self.mutate("latset = D2Q9", "latset = D2Q9\ncuda = 1"))
        self.assertTrue(spec.cuda)

    def test_grid_n_length(self) -> None:
        self.assert_error(self.mutate("n = 32, 32", "n = 32"), "needs 2 values")


def base_replace_tasks(text: str, tasks: str) -> str:
    start = text.index("[tasks]")
    end = text.index("[schedule]")
    return text[:start] + "[tasks]\n" + tasks + "\n\n" + text[end:]


class CliTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = Path(self._tmp.name)
        self.spec = self.dir / "case.case"
        self.spec.write_text(BASE, encoding="utf-8")
        self.out = self.dir / "out"

    def run_cli(self, argv: list[str]):
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            code = casegen.main(argv)
        return code, stdout.getvalue(), stderr.getvalue()

    def fake_renderer(self, files: dict[str, str]):
        calls = []

        def render(spec, out_dir, repo_root):
            calls.append((spec, out_dir, repo_root))
            return files

        return render, calls

    def test_list_prints_catalogue_without_spec(self) -> None:
        code, out, err = self.run_cli(["--list"])
        self.assertEqual(code, 0)
        self.assertEqual(err, "")
        self.assertIn("D2Q9", out)
        self.assertIn("movingwall", out)
        self.assertIn("bgk", out)
        # Phase 3 operators are part of the catalogue now; the unimplemented
        # anti_pressure boundary must stay out of it.
        for name in ("bgk_force", "mrt", "les"):
            self.assertIn(name, out)
        self.assertNotIn("anti_pressure", out)

    def test_write_then_refuse_then_force(self) -> None:
        render, calls = self.fake_renderer({"case.cpp": "A\n", "operators/x.hh": "X\n"})
        with mock.patch.object(casegen, "_load_renderer", return_value=render):
            code, out, err = self.run_cli([str(self.spec), "-o", str(self.out)])
            self.assertEqual((code, err), (0, ""))
            self.assertEqual((self.out / "case.cpp").read_text(), "A\n")
            self.assertEqual((self.out / "operators/x.hh").read_text(), "X\n")
            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0][1], self.out.resolve())
            self.assertEqual(calls[0][2], Path(casegen.__file__).resolve().parent.parent.parent)

            code, out, err = self.run_cli([str(self.spec), "-o", str(self.out)])
            self.assertEqual((code, out), (0, f"casegen: wrote 0 file(s) to {self.out}\n"))

            render, _ = self.fake_renderer({"case.cpp": "B\n", "operators/x.hh": "X\n"})
            with mock.patch.object(casegen, "_load_renderer", return_value=render):
                code, out, err = self.run_cli([str(self.spec), "-o", str(self.out)])
                self.assertEqual(code, 1)
                self.assertIn("refusing to overwrite", err)
                self.assertEqual((self.out / "case.cpp").read_text(), "A\n")

                code, out, err = self.run_cli(
                    [str(self.spec), "-o", str(self.out), "--force"]
                )
                self.assertEqual((code, err), (0, ""))
                self.assertEqual((self.out / "case.cpp").read_text(), "B\n")

    def test_check_reports_stale_and_missing_with_exit_2(self) -> None:
        render, _ = self.fake_renderer({"case.cpp": "A\n", "Makefile": "M\n"})
        with mock.patch.object(casegen, "_load_renderer", return_value=render):
            code, out, err = self.run_cli(
                [str(self.spec), "-o", str(self.out), "--check"]
            )
            self.assertEqual(code, 2)
            self.assertIn("missing", err)
            self.assertFalse(self.out.exists())

            (self.out).mkdir()
            (self.out / "case.cpp").write_text("A\n")
            (self.out / "Makefile").write_text("OLD\n")
            code, out, err = self.run_cli(
                [str(self.spec), "-o", str(self.out), "--check"]
            )
            self.assertEqual(code, 2)
            self.assertIn("Makefile: differs", err)

            (self.out / "Makefile").write_text("M\n")
            code, out, err = self.run_cli(
                [str(self.spec), "-o", str(self.out), "--check"]
            )
            self.assertEqual(code, 0)
            self.assertIn("up to date", out)

    def test_check_accepts_list_of_pairs(self) -> None:
        def render(spec, out_dir, repo_root):
            return [("a.txt", "a"), ("b/c.txt", b"c")]

        with mock.patch.object(casegen, "_load_renderer", return_value=render):
            code, _, err = self.run_cli([str(self.spec), "-o", str(self.out)])
            self.assertEqual((code, err), (0, ""))
            self.assertEqual((self.out / "b/c.txt").read_bytes(), b"c")

    def test_check_and_force_conflict(self) -> None:
        with self.assertRaises(SystemExit) as caught:
            self.run_cli([str(self.spec), "-o", str(self.out), "--check", "--force"])
        self.assertEqual(caught.exception.code, 2)

    def test_invalid_spec_does_not_load_the_emitter(self) -> None:
        self.spec.write_text(self.spec.read_text().replace("dim = 2", "dim = two"))
        with mock.patch.object(casegen, "_load_renderer") as loader:
            code, out, err = self.run_cli([str(self.spec), "-o", str(self.out)])
        self.assertEqual(code, 1)
        self.assertIn("invalid integer", err)
        self.assertIn(f"{self.spec}:3:", err)
        loader.assert_not_called()

    def test_emitter_error_is_reported(self) -> None:
        def render(spec, out_dir, repo_root):
            raise CasegenError("unsupported combination", None)

        with mock.patch.object(casegen, "_load_renderer", return_value=render):
            code, out, err = self.run_cli([str(self.spec), "-o", str(self.out)])
        self.assertEqual(code, 1)
        self.assertIn("unsupported combination", err)

    def test_rendered_path_escape_is_rejected(self) -> None:
        def render(spec, out_dir, repo_root):
            return {"../escape.txt": "x"}

        with mock.patch.object(casegen, "_load_renderer", return_value=render):
            code, out, err = self.run_cli([str(self.spec), "-o", str(self.out)])
        self.assertEqual(code, 1)
        self.assertIn("escapes the output directory", err)
        self.assertFalse((self.dir / "escape.txt").exists())


if __name__ == "__main__":
    unittest.main()
