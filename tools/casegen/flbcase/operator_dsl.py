"""Per-cell operator DSL: parse a ``.op`` file and render a C++ header.

The language is deliberately small.  Every builtin expands to the arithmetic of
the reference scheme in ``src/lbm`` -- same operands, same association -- so a
generated operator can be bit-identical to the hand-written one it mirrors.
Parsing, type checking and emission all live here; ``model.OperatorSpec`` only
carries the result.

Grammar (comments start a line with ``#``)::

    operator NAME {
      var n: scalar | vector | array
      ...
      body {
        n = expr
        n[k] = expr
        for i in q { ... }
      }
    }

``f`` is the population alias of the cell and ``q`` the direction count.
Vectors combine with ``+ -`` and with scalars; arrays are element-wise only.
``u[k]`` reads a vector component (integer literal only, ``0 <= k < dim``).
``field(NAME)`` reads RHO, VELOCITY, FORCE or SMAGORINSKY from the cell;
``sqrt``, ``abs``, ``min``, ``max`` and ``clamp`` are scalar math builtins.
``moment0``, ``moment1``, ``moment1_force`` and ``force_term`` emit statements,
so they must be the whole right-hand side of an assignment.

Every diagnostic is a ``model.CasegenError`` carrying ``<path>:<line>``.
"""

from __future__ import annotations

import hashlib
import os
import re
from dataclasses import dataclass, field
from difflib import get_close_matches
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Union

from . import model


# ---------------------------------------------------------------------------
# Names the generated code owns
# ---------------------------------------------------------------------------

_CPP_KEYWORDS = frozenset(
    """
    alignas alignof and asm auto bool break case catch char class const
    constexpr continue default delete do double else enum explicit export
    extern false float for friend goto if inline int long mutable namespace
    new noexcept not nullptr operator or private protected public register
    return short signed sizeof static struct switch template this throw true
    try typedef typename union unsigned using virtual void volatile while xor
    """.split()
)

# `f` is the population alias, `q` the loop bound, `i` the loop index used by
# the builtin reductions; the rest are type names or namespaces the generated
# code spells explicitly.  `field` is call syntax, never a variable.
_RESERVED_NAMES = frozenset(
    {"f", "q", "omega", "fomega", "cell", "T", "LatSet", "CELL",
     "std", "force", "field", "latset", "Vector"}
)

# Locals declared by the emitted reductions and by feq2's hoisted terms.
_RESERVED_PATTERN = re.compile(r"^(?:u2|uc|rho_value)(?:_\d+)?$")

_REDUCERS = frozenset({"moment0", "moment1", "moment1_force", "force_term"})

# Scalar builtins: name -> arity.  They emit the std:: function of the same
# name, so every argument and the result stay T.
_MATH_SCALAR_FUNCS = {"sqrt": 1, "abs": 1, "min": 2, "max": 2, "clamp": 3}

# Fields readable through field(NAME).  POP is deliberately absent: the
# populations are addressed through `f`, moment0() and moment1().
_FIELD_KINDS = {
    "RHO": "scalar",
    "VELOCITY": "vector",
    "FORCE": "vector",
    "SMAGORINSKY": "scalar",
}

_CALLS = (
    frozenset({"feq2", "w", "c", "dot", "norm2", "force", "field"})
    | frozenset(_MATH_SCALAR_FUNCS)
)
_BUILTINS = _REDUCERS | _CALLS

_VAR_TYPES = ("scalar", "vector", "array")


def _error(path: str, line: int, message: str) -> "model.CasegenError":
    return model.CasegenError(message, model.SourceSpan(path, line))


def pascal_name(name: str) -> str:
    """C++ struct name of an operator: ``my_bgk`` -> ``MyBgk``.

    The emitter needs the same name when it wires the task into the generated
    case, so this is public rather than a renderer detail.
    """
    return "".join(part[:1].upper() + part[1:] for part in name.split("_") if part)


def _kind_label(kind: str) -> str:
    return "loop index" if kind == "index" else kind


# ---------------------------------------------------------------------------
# AST
# ---------------------------------------------------------------------------


@dataclass
class _Literal:
    text: str
    is_int: bool
    line: int


@dataclass
class _Name:
    name: str
    line: int


@dataclass
class _Index:
    base: str
    index: "_Expr"
    line: int


@dataclass
class _Unary:
    op: str
    operand: "_Expr"
    line: int


@dataclass
class _Binary:
    op: str
    left: "_Expr"
    right: "_Expr"
    line: int


@dataclass
class _Call:
    name: str
    args: List["_Expr"]
    line: int


@dataclass
class _Assign:
    target: Union[_Name, _Index]
    expr: "_Expr"
    line: int


@dataclass
class _For:
    index: str
    body: List[Union[_Assign, "_For"]] = field(default_factory=list)
    line: int = 0


_Expr = Union[_Literal, _Name, _Index, _Unary, _Binary, _Call]
_Stmt = Union[_Assign, _For]


@dataclass
class _Parsed:
    name: str
    variables: Dict[str, str]
    body: List[_Stmt]


# ---------------------------------------------------------------------------
# Lexer
# ---------------------------------------------------------------------------


@dataclass
class _Token:
    kind: str          # ident | int | float | sym | eof
    text: str
    line: int


_SYMBOLS = frozenset("{}[]():;,+-*/=")


def _lex(text: str, path: str) -> List[_Token]:
    tokens: List[_Token] = []
    i = 0
    line = 1
    line_has_token = False
    n = len(text)
    while i < n:
        ch = text[i]
        if ch == "\n":
            line += 1
            line_has_token = False
            i += 1
            continue
        if ch in " \t\r":
            i += 1
            continue
        if ch == "#":
            if line_has_token:
                raise _error(path, line, "a comment must start a line: put '#' first")
            while i < n and text[i] != "\n":
                i += 1
            continue
        line_has_token = True
        if ch.isalpha() or ch == "_":
            j = i
            while j < n and (text[j].isalnum() or text[j] == "_"):
                j += 1
            tokens.append(_Token("ident", text[i:j], line))
            i = j
            continue
        if ch in "0123456789" or (ch == "." and i + 1 < n and text[i + 1] in "0123456789"):
            j = i
            is_float = False
            while j < n and text[j] in "0123456789":
                j += 1
            if j < n and text[j] == ".":
                is_float = True
                j += 1
                while j < n and text[j] in "0123456789":
                    j += 1
            if j < n and text[j] in "eE":
                k = j + 1
                if k < n and text[k] in "+-":
                    k += 1
                if k >= n or text[k] not in "0123456789":
                    raise _error(path, line, f"malformed exponent in number '{text[i:j + 2]}'")
                is_float = True
                j = k
                while j < n and text[j] in "0123456789":
                    j += 1
            tokens.append(_Token("float" if is_float else "int", text[i:j], line))
            i = j
            continue
        if ch in _SYMBOLS:
            tokens.append(_Token("sym", ch, line))
            i += 1
            continue
        raise _error(path, line, f"unexpected character {ch!r}")
    tokens.append(_Token("eof", "", line))
    return tokens


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


class _Parser:
    def __init__(self, tokens: List[_Token], path: str) -> None:
        self._tokens = tokens
        self._pos = 0
        self._path = path

    # -- token helpers ----------------------------------------------------

    def _peek(self, ahead: int = 0) -> _Token:
        index = min(self._pos + ahead, len(self._tokens) - 1)
        return self._tokens[index]

    def _take(self) -> _Token:
        token = self._tokens[self._pos]
        self._pos += 1
        return token

    def _at_ident(self, text: str) -> bool:
        token = self._peek()
        return token.kind == "ident" and token.text == text

    def _at_sym(self, text: str) -> bool:
        token = self._peek()
        return token.kind == "sym" and token.text == text

    def _error(self, token: _Token, message: str) -> "model.CasegenError":
        return _error(self._path, token.line, message)

    def _expect_sym(self, text: str) -> _Token:
        token = self._peek()
        if token.kind != "sym" or token.text != text:
            found = "end of file" if token.kind == "eof" else f"'{token.text}'"
            raise self._error(token, f"expected '{text}', found {found}")
        return self._take()

    def _expect_ident(self, what: str) -> _Token:
        token = self._peek()
        if token.kind != "ident":
            found = "end of file" if token.kind == "eof" else f"'{token.text}'"
            raise self._error(token, f"expected {what}, found {found}")
        return self._take()

    # -- grammar ----------------------------------------------------------

    def parse(self) -> _Parsed:
        if not self._at_ident("operator"):
            raise self._error(self._peek(), "expected 'operator NAME {'")
        self._take()
        name_token = self._expect_ident("the operator name")
        name = name_token.text
        if name in _CPP_KEYWORDS:
            raise self._error(name_token, f"'{name}' is a C++ keyword")
        self._expect_sym("{")

        variables: Dict[str, str] = {}
        while self._at_ident("var"):
            self._take()
            var_token = self._expect_ident("a variable name")
            varname = var_token.text
            self._check_decl_name(varname, var_token)
            if varname in variables:
                raise self._error(var_token, f"variable '{varname}' is declared twice")
            self._expect_sym(":")
            kind_token = self._expect_ident("a type (scalar, vector or array)")
            if kind_token.text not in _VAR_TYPES:
                raise self._error(
                    kind_token,
                    f"unknown type '{kind_token.text}'; expected scalar, vector or array",
                )
            variables[varname] = kind_token.text

        if not self._at_ident("body"):
            token = self._peek()
            if token.kind == "eof":
                raise self._error(token, "expected 'body {'")
            raise self._error(token, "expected 'body {' after the variable declarations")
        self._take()
        self._expect_sym("{")
        body = self._parse_statements()
        self._expect_sym("}")
        self._expect_sym("}")
        token = self._peek()
        if token.kind != "eof":
            raise self._error(token, "unexpected text after the operator definition")
        return _Parsed(name=name, variables=variables, body=body)

    def _check_decl_name(self, name: str, token: _Token, loop_index: bool = False) -> None:
        if name in _CPP_KEYWORDS:
            raise self._error(token, f"'{name}' is a C++ keyword")
        # `i` names the loop of the emitted reductions; a direction loop may
        # still use it (reductions are rejected inside one).
        if name in _RESERVED_NAMES or _RESERVED_PATTERN.match(name) or (name == "i" and not loop_index):
            raise self._error(token, f"'{name}' is reserved by the generated code")

    def _parse_statements(self) -> List[_Stmt]:
        statements: List[_Stmt] = []
        while True:
            token = self._peek()
            if token.kind == "eof":
                raise self._error(token, "expected '}' to close the block")
            if token.kind == "sym" and token.text == "}":
                return statements
            if token.kind == "sym" and token.text == ";":
                raise self._error(token, "statements are separated by lines, not ';'")
            statements.append(self._parse_statement())

    def _parse_statement(self) -> _Stmt:
        token = self._peek()
        if token.kind == "ident" and token.text == "for":
            return self._parse_for()
        if token.kind == "ident" and token.text in ("var", "body", "operator"):
            raise self._error(token, f"unexpected '{token.text}' inside a body")
        if token.kind != "ident":
            found = "end of file" if token.kind == "eof" else f"'{token.text}'"
            raise self._error(token, f"expected an assignment or 'for i in q {{', found {found}")

        self._take()
        target: Union[_Name, _Index] = _Name(token.text, token.line)
        if self._at_sym("["):
            self._take()
            index = self._parse_expr()
            self._expect_sym("]")
            target = _Index(token.text, index, token.line)
        self._expect_sym("=")
        expr = self._parse_expr()
        return _Assign(target, expr, token.line)

    def _parse_for(self) -> _For:
        for_token = self._take()
        index_token = self._expect_ident("a loop index")
        self._check_decl_name(index_token.text, index_token, loop_index=True)
        if not self._at_ident("in"):
            raise self._error(self._peek(), "expected 'in q' after the loop index")
        self._take()
        if not self._at_ident("q"):
            raise self._error(self._peek(), "expected 'q' (the direction count) as the loop bound")
        self._take()
        self._expect_sym("{")
        body = self._parse_statements()
        self._expect_sym("}")
        return _For(index=index_token.text, body=body, line=for_token.line)

    # -- expressions ------------------------------------------------------

    def _parse_expr(self) -> _Expr:
        left = self._parse_term()
        while self._peek().kind == "sym" and self._peek().text in "+-":
            op_token = self._take()
            right = self._parse_term()
            left = _Binary(op_token.text, left, right, op_token.line)
        return left

    def _parse_term(self) -> _Expr:
        left = self._parse_unary()
        while self._peek().kind == "sym" and self._peek().text in "*/":
            op_token = self._take()
            right = self._parse_unary()
            left = _Binary(op_token.text, left, right, op_token.line)
        return left

    def _parse_unary(self) -> _Expr:
        if self._at_sym("-"):
            op_token = self._take()
            return _Unary("-", self._parse_unary(), op_token.line)
        return self._parse_postfix()

    def _parse_postfix(self) -> _Expr:
        node = self._parse_primary()
        while self._at_sym("["):
            bracket = self._take()
            if isinstance(node, _Index):
                raise self._error(bracket, "nested subscripts are not supported")
            if not isinstance(node, _Name):
                raise self._error(bracket, "only a variable can be subscripted")
            index = self._parse_expr()
            self._expect_sym("]")
            node = _Index(node.name, index, node.line)
        return node

    def _parse_primary(self) -> _Expr:
        token = self._peek()
        if token.kind in ("int", "float"):
            self._take()
            return _Literal(token.text, token.kind == "int", token.line)
        if token.kind == "ident":
            self._take()
            if self._at_sym("("):
                self._take()
                args: List[_Expr] = []
                if not self._at_sym(")"):
                    args.append(self._parse_expr())
                    while self._at_sym(","):
                        self._take()
                        args.append(self._parse_expr())
                self._expect_sym(")")
                return _Call(token.text, args, token.line)
            return _Name(token.text, token.line)
        if token.kind == "sym" and token.text == "(":
            self._take()
            expr = self._parse_expr()
            self._expect_sym(")")
            return expr
        found = "end of file" if token.kind == "eof" else f"'{token.text}'"
        raise self._error(token, f"expected an expression, found {found}")


# ---------------------------------------------------------------------------
# Type checker
# ---------------------------------------------------------------------------


def _walk_expr(node: _Expr):
    yield node
    if isinstance(node, _Unary):
        yield from _walk_expr(node.operand)
    elif isinstance(node, _Binary):
        yield from _walk_expr(node.left)
        yield from _walk_expr(node.right)
    elif isinstance(node, _Index):
        yield from _walk_expr(node.index)
    elif isinstance(node, _Call):
        for arg in node.args:
            yield from _walk_expr(arg)


def _called(node: _Expr, name: str) -> List[_Call]:
    return [n for n in _walk_expr(node) if isinstance(n, _Call) and n.name == name]


def _used_names(statements: Sequence[_Stmt]) -> Set[str]:
    """Every identifier referenced by the body (used to hoist omega/fomega)."""
    names: Set[str] = set()
    for statement in statements:
        if isinstance(statement, _For):
            names |= _used_names(statement.body)
            continue
        names |= {n.name for n in _walk_expr(statement.expr) if isinstance(n, _Name)}
        if isinstance(statement.target, _Index):
            names |= {n.name for n in _walk_expr(statement.target.index) if isinstance(n, _Name)}
    return names


def _used_calls(statements: Sequence[_Stmt]) -> Set[str]:
    """Names of the builtins called by the body (selects the std headers)."""
    names: Set[str] = set()
    for statement in statements:
        if isinstance(statement, _For):
            names |= _used_calls(statement.body)
            continue
        names |= {n.name for n in _walk_expr(statement.expr) if isinstance(n, _Call)}
    return names


class _Checker:
    def __init__(
        self,
        path: str,
        variables: Dict[str, str],
        q: Optional[int] = None,
        dim: Optional[int] = None,
    ) -> None:
        self._path = path
        self._vars = dict(variables)
        self._q = q
        self._dim = dim
        self._defined: Set[str] = set()
        self._loop: Optional[str] = None
        self._loop_defined: Set[str] = set()
        self.needs: Set[str] = {"POP"}

    # -- helpers ----------------------------------------------------------

    def _error(self, line: int, message: str) -> "model.CasegenError":
        return _error(self._path, line, message)

    def _kind(self, name: str) -> Optional[str]:
        if name == "f":
            return "array"
        return self._vars.get(name)

    def _unknown(self, node: _Name) -> "model.CasegenError":
        name = node.name
        if name in model.FIELD_TYPES:
            if name == "POP":
                return self._error(
                    node.line,
                    "'POP' is the population alias f; use moment0(f) or f[i]",
                )
            return self._error(
                node.line,
                f"field '{name}' cannot be read as a variable; use field({name})",
            )
        if name in _BUILTINS:
            return self._error(node.line, f"'{name}' is a builtin function; call it as {name}(...)")
        if name in ("var", "body", "for", "in", "operator"):
            return self._error(node.line, f"unexpected keyword '{name}'")
        candidates = get_close_matches(name, list(self._vars) + ["f", "omega", "fomega"], n=1)
        hint = f"; did you mean '{candidates[0]}'?" if candidates else ""
        return self._error(node.line, f"unknown variable '{name}'{hint}")

    def _require_defined(self, node: _Name) -> None:
        if node.name not in self._defined:
            raise self._error(
                node.line,
                f"variable '{node.name}' is read before it is assigned",
            )

    # -- entry ------------------------------------------------------------

    def check(self, body: Sequence[_Stmt]) -> Set[str]:
        for statement in body:
            self._check_statement(statement)
        return set(self.needs)

    def _check_statement(self, statement: _Stmt) -> None:
        if isinstance(statement, _For):
            self._check_for(statement)
            return
        self._check_assign(statement)

    def _check_for(self, statement: _For) -> None:
        if self._loop is not None:
            raise self._error(
                statement.line,
                f"loops cannot be nested: 'for {self._loop} in q' already covers every direction",
            )
        if statement.index in self._vars:
            raise self._error(
                statement.line,
                f"loop index '{statement.index}' shadows a declared variable",
            )
        outer_loop, outer_defined = self._loop, self._loop_defined
        self._loop = statement.index
        self._loop_defined = set(self._defined)
        try:
            for inner in statement.body:
                self._check_statement(inner)
        finally:
            self._loop, self._loop_defined = outer_loop, outer_defined

    def _check_assign(self, statement: _Assign) -> None:
        target = statement.target
        if isinstance(target, _Name):
            self._check_name_target(statement, target)
        else:
            self._check_index_target(statement, target)

    def _check_name_target(self, statement: _Assign, target: _Name) -> None:
        name = target.name
        if name in ("omega", "fomega", "q"):
            raise self._error(statement.line, f"'{name}' cannot be assigned")
        kind = self._kind(name)
        if kind is None:
            raise self._unknown(target)
        if kind == "array":
            self._check_force_term_target(statement, name)
            return
        expr = statement.expr
        if isinstance(expr, _Call) and expr.name in _REDUCERS:
            got = self._check_reducer(expr, target=name)
        else:
            got = self._check_value(expr)
        if got != kind:
            raise self._error(
                statement.line,
                f"cannot assign {got} to {kind} variable '{name}'",
            )
        self._defined.add(name)

    def _check_index_target(self, statement: _Assign, target: _Index) -> None:
        name = target.base
        kind = self._kind(name)
        if kind is None:
            raise self._unknown(_Name(name, target.line))
        if kind == "vector":
            raise self._error(
                target.line,
                f"vector '{name}' cannot be subscripted on the left-hand side; "
                f"components u[0]..u[2] are read-only, so update the whole "
                f"vector with +, - or a scalar factor",
            )
        if kind != "array":
            raise self._error(target.line, f"'{name}' is a {_kind_label(kind)} and cannot be subscripted")
        self._check_index(target.index, "the subscript", self._q, "q")
        expr = statement.expr
        if isinstance(expr, _Call) and expr.name in _REDUCERS:
            raise self._error(
                statement.line,
                f"{expr.name}() produces a whole-cell value; assign it to a variable, "
                f"not to an array element",
            )
        got = self._check_value(expr)
        if got != "scalar":
            raise self._error(statement.line, f"cannot assign {got} to an array element")
        self._defined.add(name)

    def _check_force_term_target(self, statement: _Assign, name: str) -> None:
        expr = statement.expr
        if name == "f":
            raise self._error(
                statement.line,
                "f is the population alias, not an array object; "
                "force_term() must fill a declared array variable",
            )
        if not (isinstance(expr, _Call) and expr.name == "force_term"):
            raise self._error(
                statement.line,
                f"array '{name}' cannot be assigned as a whole; assign elements "
                f"({name}[i] = ...) or use {name} = force_term(u, F)",
            )
        if self._loop is not None:
            raise self._error(
                statement.line,
                "force_term() fills the whole array and cannot be used inside a direction loop",
            )
        self._check_arity(expr, 2, "force_term(u, F)")
        u_kind = self._check_value(expr.args[0])
        self._expect_vector(u_kind, expr.args[0], "the first argument of force_term()")
        f_kind = self._check_value(expr.args[1])
        self._expect_vector(f_kind, expr.args[1], "the second argument of force_term()")
        if any(isinstance(n, _Name) and n.name == name for n in _walk_expr(expr)):
            raise self._error(
                statement.line,
                f"the target '{name}' cannot also be an argument of force_term()",
            )
        self.needs.add("FORCE")
        self._defined.add(name)

    # -- expressions ------------------------------------------------------

    def _check_value(self, node: _Expr, expected: Optional[str] = None) -> str:
        kind = self._check_expr(node)
        if kind == "array":
            raise self._error(
                node.line,
                "an array can only be used element-wise, e.g. a[i]",
            )
        if kind == "index":
            raise self._error(node.line, "a loop index cannot be used as a value")
        if expected is not None and kind != expected:
            raise self._error(node.line, f"expected {expected}, found {kind}")
        return kind

    def _expect_vector(self, kind: str, node: _Expr, what: str) -> None:
        if kind != "vector":
            raise self._error(node.line, f"{what} must be a vector, found {kind}")

    def _check_expr(self, node: _Expr) -> str:
        if isinstance(node, _Literal):
            return "scalar"
        if isinstance(node, _Name):
            return self._check_name(node)
        if isinstance(node, _Index):
            base = _Name(node.base, node.line)
            kind = self._check_name(base)
            if kind == "vector":
                self._check_component(node)
                return "scalar"
            if kind == "array":
                if node.base != "f":
                    self._require_defined(base)
            else:
                raise self._error(node.line, f"'{node.base}' is a {_kind_label(kind)} and cannot be subscripted")
            self._check_index(node.index, "the subscript", self._q, "q")
            return "scalar"
        if isinstance(node, _Unary):
            kind = self._check_value(node.operand)
            if kind == "vector":
                raise self._error(
                    node.line,
                    "unary minus is not defined for vectors; write (T{0} - a) instead",
                )
            return "scalar"
        if isinstance(node, _Binary):
            return self._check_binary(node)
        if isinstance(node, _Call):
            return self._check_call(node)
        raise self._error(0, "internal: unknown expression node")  # pragma: no cover

    def _check_name(self, node: _Name) -> str:
        name = node.name
        if name in ("omega", "fomega"):
            return "scalar"
        if name == self._loop:
            return "index"
        kind = self._kind(name)
        if kind is None:
            raise self._unknown(node)
        if kind != "array":
            self._require_defined(node)
        return kind

    def _check_component(self, node: _Index) -> None:
        """A vector component is a read-only integer literal ``u[0..dim-1]``."""
        index = node.index
        if not isinstance(index, _Literal) or not index.is_int:
            raise self._error(
                index.line,
                f"a component of vector '{node.base}' must be an integer "
                f"literal 0, 1 or 2",
            )
        value = int(index.text)
        limit = self._dim if self._dim is not None else 3
        if not 0 <= value < limit:
            allowed = (
                f"use 0, 1 or 2"
                if self._dim is None
                else f"use 0..{limit - 1} (dim = {limit})"
            )
            raise self._error(
                index.line,
                f"component index {value} of vector '{node.base}' is out of "
                f"range; {allowed}",
            )

    def _check_binary(self, node: _Binary) -> str:
        left = self._check_value(node.left)
        right = self._check_value(node.right)
        if left == "vector" and right == "vector":
            if node.op in "+-":
                return "vector"
            if node.op == "*":
                raise self._error(
                    node.line,
                    "vector * vector is not defined; use dot(a, b)",
                )
            raise self._error(
                node.line,
                "vector / vector is not defined; use dot(a, b) and a scalar divisor",
            )
        if left == "vector" or right == "vector":
            if left == "scalar" and right == "vector" and node.op == "/":
                raise self._error(
                    node.line,
                    "scalar / vector is not defined; multiply by a scalar factor instead",
                )
            return "vector"
        return "scalar"

    def _check_call(self, node: _Call) -> str:
        name = node.name
        if name in _REDUCERS:
            raise self._error(
                node.line,
                f"{name}() expands to statements and can only be the whole "
                f"right-hand side of an assignment",
            )
        if name == "field":
            return self._check_field(node)
        if name in _MATH_SCALAR_FUNCS:
            return self._check_scalar_call(node)
        if name == "feq2":
            self._check_arity(node, 3, "feq2(rho, u, i)")
            self._check_value(node.args[0], "scalar")
            velocity = self._check_value(node.args[1])
            self._expect_vector(velocity, node.args[1], "the second argument of feq2()")
            self._check_index(
                node.args[2], "the third argument of feq2()", self._q, "q"
            )
            if self._loop is not None:
                for inner in _walk_expr(node.args[1]):
                    if isinstance(inner, _Name) and inner.name in self._vars:
                        if inner.name not in self._loop_defined:
                            raise self._error(
                                inner.line,
                                f"feq2()'s velocity uses '{inner.name}', which is assigned "
                                f"inside the loop; u2 is hoisted before it",
                            )
            if _called(node.args[1], "feq2"):
                raise self._error(node.line, "feq2() cannot be nested in its own velocity")
            return "scalar"
        if name in ("w", "c"):
            self._check_arity(node, 1, f"{name}(i)")
            self._check_index(node.args[0], f"the argument of {name}()", self._q, "q")
            return "scalar" if name == "w" else "vector"
        if name == "dot":
            self._check_arity(node, 2, "dot(a, b)")
            self._expect_vector(self._check_value(node.args[0]), node.args[0], "dot()'s first argument")
            self._expect_vector(self._check_value(node.args[1]), node.args[1], "dot()'s second argument")
            return "scalar"
        if name == "norm2":
            self._check_arity(node, 1, "norm2(a)")
            self._expect_vector(self._check_value(node.args[0]), node.args[0], "norm2()'s argument")
            return "scalar"
        if name == "force":
            self._check_arity(node, 0, "force()")
            self.needs.add("FORCE")
            return "vector"
        if name == "force_term":
            raise self._error(
                node.line,
                "force_term(u, F) fills an array and can only be the right-hand "
                "side of `array = force_term(u, F)`",
            )
        candidates = get_close_matches(name, sorted(_BUILTINS), n=1)
        hint = f"; did you mean '{candidates[0]}'?" if candidates else ""
        raise self._error(node.line, f"unknown function '{name}'{hint}")

    def _check_field(self, node: _Call) -> str:
        self._check_arity(node, 1, "field(RHO|VELOCITY|FORCE|SMAGORINSKY)")
        arg = node.args[0]
        if not isinstance(arg, _Name):
            raise self._error(
                arg.line,
                "field() takes a field name: RHO, VELOCITY, FORCE or SMAGORINSKY",
            )
        name = arg.name
        if name not in _FIELD_KINDS:
            candidates = get_close_matches(name, sorted(_FIELD_KINDS), n=1)
            hint = f"; did you mean '{candidates[0]}'?" if candidates else ""
            raise self._error(
                arg.line,
                f"unknown field '{name}'; field() supports RHO, VELOCITY, "
                f"FORCE and SMAGORINSKY{hint}",
            )
        self.needs.add(name)
        return _FIELD_KINDS[name]

    def _check_scalar_call(self, node: _Call) -> str:
        name = node.name
        arity = _MATH_SCALAR_FUNCS[name]
        signatures = {
            "sqrt": "sqrt(x)",
            "abs": "abs(x)",
            "min": "min(a, b)",
            "max": "max(a, b)",
            "clamp": "clamp(x, lo, hi)",
        }
        self._check_arity(node, arity, signatures[name])
        for arg in node.args:
            self._check_value(arg, "scalar")
        return "scalar"

    def _check_reducer(self, node: _Call, target: str) -> str:
        name = node.name
        if name == "force_term":
            raise self._error(
                node.line,
                "force_term(u, F) fills an array and can only be the right-hand "
                "side of `array = force_term(u, F)`",
            )
        if self._loop is not None:
            raise self._error(
                node.line,
                f"{name}() reduces over all directions and cannot be used inside "
                f"a direction loop",
            )
        if name == "moment0":
            self._check_arity(node, 1, "moment0(f)")
            self._check_populations(node, node.args[0], name)
            return "scalar"
        if name == "moment1":
            self._check_arity(node, 2, "moment1(f, rho)")
            self._check_populations(node, node.args[0], name)
            self._check_value(node.args[1], "scalar")
            self._reject_feq2_args(node)
            return "vector"
        if name == "moment1_force":
            self._check_arity(node, 2, "moment1_force(f, F)")
            self._check_populations(node, node.args[0], name)
            self._expect_vector(
                self._check_value(node.args[1]), node.args[1],
                "the second argument of moment1_force()",
            )
            self._reject_feq2_args(node)
            if any(isinstance(n, _Name) and n.name == target for n in _walk_expr(node.args[1])):
                raise self._error(
                    node.line,
                    f"the target '{target}' cannot also be its force argument: "
                    f"the target is cleared before the force is read",
                )
            return "vector"
        # force_term never reaches here: _check_force_term_target handles it.
        raise self._error(node.line, f"{name}() cannot be assigned element-wise")  # pragma: no cover

    def _check_populations(self, node: _Call, arg: _Expr, name: str) -> None:
        if not (isinstance(arg, _Name) and arg.name == "f"):
            raise self._error(
                arg.line,
                f"the first argument of {name}() must be f, the populations",
            )

    def _reject_feq2_args(self, node: _Call) -> None:
        for arg in node.args[1:]:
            if _called(arg, "feq2"):
                raise self._error(
                    node.line,
                    f"feq2() cannot be used inside a {node.name}() argument",
                )

    def _check_arity(self, node: _Call, count: int, signature: str) -> None:
        if len(node.args) != count:
            raise self._error(
                node.line,
                f"{node.name}() takes {count} argument{'s' if count != 1 else ''} "
                f"({signature}), got {len(node.args)}",
            )

    def _check_index(
        self,
        node: _Expr,
        what: str,
        limit: Optional[int] = None,
        limit_name: str = "q",
    ) -> None:
        if isinstance(node, _Literal):
            if not node.is_int:
                raise self._error(
                    node.line,
                    f"{what} must be the loop index or an integer literal, "
                    f"found the float literal {node.text}",
                )
            value = int(node.text)
            if limit is not None and not 0 <= value < limit:
                raise self._error(
                    node.line,
                    f"index {value} is out of range for {what}; valid indices "
                    f"are 0..{limit - 1} ({limit_name} = {limit})",
                )
            return
        if isinstance(node, _Name):
            if node.name == self._loop and self._loop is not None:
                return
            raise self._error(
                node.line,
                f"{what} must be the loop index or an integer literal; "
                f"'{node.name}' is neither",
            )
        raise self._error(
            node.line,
            f"{what} must be the loop index or an integer literal",
        )


# ---------------------------------------------------------------------------
# C++ emission
# ---------------------------------------------------------------------------


class _Feq2Ctx:
    """Allocates feq2's hoisted `u2` and per-call `uc` locals."""

    def __init__(self, renderer: "_Renderer", allow_new_u2: bool) -> None:
        self._renderer = renderer
        self._allow_new = allow_new_u2
        self._u2: Dict[str, str] = {} if allow_new_u2 else dict(renderer._loop_u2 or {})
        self._uc_seq = 0
        self.decls: List[str] = []

    def u2_name(self, velocity: str) -> str:
        name = self._u2.get(velocity)
        if name is not None:
            return name
        if not self._allow_new:
            # The loop pre-scan walks the same AST, so this is an internal bug.
            raise RuntimeError(f"internal: feq2 velocity '{velocity}' was not hoisted")
        name = self._renderer.new_u2_name()
        self._u2[velocity] = name
        self.decls.append(f"const T {name} = ({velocity}).getnorm2();")
        return name

    def uc_name(self, velocity: str, index: str) -> str:
        name = "uc" if self._uc_seq == 0 else f"uc_{self._uc_seq}"
        self._uc_seq += 1
        self.decls.append(f"const T {name} = {velocity} * latset::c<LatSet>({index});")
        return name


class _Renderer:
    def __init__(self, spec: model.OperatorSpec) -> None:
        self._spec = spec
        self._vars = spec.variables
        self._u2_seq = 0
        self._loop: Optional[str] = None
        self._loop_u2: Optional[Dict[str, str]] = None

    # -- names ------------------------------------------------------------

    def new_u2_name(self) -> str:
        name = "u2" if self._u2_seq == 0 else f"u2_{self._u2_seq}"
        self._u2_seq += 1
        return name

    @staticmethod
    def _declare(name: str, kind: str) -> str:
        if kind == "scalar":
            return f"T {name}{{}};"
        if kind == "vector":
            return f"Vector<T, LatSet::d> {name}{{}};"
        return f"std::array<T, LatSet::q> {name}{{}};"

    # -- data -------------------------------------------------------------

    def render(self) -> str:
        spec = self._spec
        struct = pascal_name(spec.name)
        lines: List[str] = [
            f"// generated by casegen from {os.path.basename(spec.path)} -- do not edit",
            f"// operator {spec.name}: per-cell dynamics expanded from the .op DSL.",
            "// Include this header after freelb.h.",
            "",
            "#pragma once",
            "",
        ]
        if _used_calls(spec.body) & set(_MATH_SCALAR_FUNCS):
            lines += [
                "#include <algorithm>",
                "#include <array>",
                "#include <cmath>",
            ]
        else:
            lines.append("#include <array>")
        lines += [
            "",
            '#include "lbm/lattice_set.h"',
        ]
        if "FORCE" in spec.needs:
            lines.append('#include "lbm/force.h"')
        lines += [
            "",
            "template <typename CELL>",
            f"struct {struct} {{",
            "  using T = typename CELL::FloatType;",
            "  using LatSet = typename CELL::LatticeSet;",
            "",
            "  __any__ static void apply(CELL& cell) {",
        ]
        for name, kind in spec.variables.items():
            lines.append("    " + self._declare(name, kind))
        used = _used_names(spec.body)
        if "omega" in used:
            lines.append("    const T omega = cell.getOmega();")
        if "fomega" in used:
            lines.append("    const T fomega = cell.getfOmega();")
        if spec.variables or used:
            lines.append("")
        for statement in spec.body:
            lines.extend(self._statement(statement, 2))
        lines += [
            "  }",
            "};",
            "",
        ]
        return "\n".join(lines)

    # -- statements -------------------------------------------------------

    def _statement(self, statement: _Stmt, level: int) -> List[str]:
        if isinstance(statement, _For):
            return self._for_lines(statement, level)
        return self._assign_lines(statement, level)

    def _for_lines(self, statement: _For, level: int) -> List[str]:
        base = "  " * level
        velocities: List[str] = []
        for inner in statement.body:
            for call in _called(inner.expr, "feq2"):
                velocity = self._expr(call.args[1])
                if velocity not in velocities:
                    velocities.append(velocity)

        lines: List[str] = []
        loop_u2: Dict[str, str] = {}
        for velocity in velocities:
            name = self.new_u2_name()
            loop_u2[velocity] = name
            lines.append(f"{base}const T {name} = ({velocity}).getnorm2();")

        index = statement.index
        lines.append(
            f"{base}for (unsigned int {index} = 0; {index} < LatSet::q; ++{index}) {{"
        )
        previous_loop, previous_u2 = self._loop, self._loop_u2
        self._loop, self._loop_u2 = index, loop_u2
        try:
            for inner in statement.body:
                lines.extend(self._statement(inner, level + 1))
        finally:
            self._loop, self._loop_u2 = previous_loop, previous_u2
        lines.append(f"{base}}}")
        return lines

    def _assign_lines(self, statement: _Assign, level: int) -> List[str]:
        base = "  " * level
        target = statement.target
        expr = statement.expr

        if isinstance(expr, _Call):
            if expr.name == "force_term":
                # Only reachable for `array = force_term(u, F)`; the checker
                # rejects every other use.
                name = target.name if isinstance(target, _Name) else ""
                return [
                    f"{base}force::Force<CELL>::apply("
                    f"{self._expr(expr.args[0])}, {self._expr(expr.args[1])}, {name});"
                ]
            if expr.name in _REDUCERS:
                return self._reducer_lines(statement, level)

        target_code = self._target_code(target)
        if not _called(expr, "feq2"):
            return [f"{base}{target_code} = {self._expr(expr)};"]

        ctx = _Feq2Ctx(self, allow_new_u2=self._loop is None)
        value = self._expr(expr, ctx)
        lines = [f"{base}{{"]
        lines += [f"{base}  {decl}" for decl in ctx.decls]
        lines.append(f"{base}  {target_code} = {value};")
        lines.append(f"{base}}}")
        return lines

    def _target_code(self, target: Union[_Name, _Index]) -> str:
        if isinstance(target, _Name):
            return target.name
        # f is the population alias, not an object: cell[i] is the lvalue.
        if target.base == "f":
            return f"cell[{self._index_code(target.index)}]"
        return f"{target.base}[{self._index_code(target.index)}]"

    def _reducer_lines(self, statement: _Assign, level: int) -> List[str]:
        base = "  " * level
        call = statement.expr
        assert isinstance(call, _Call)
        target = self._target_code(statement.target)

        if call.name == "moment0":
            return [
                f"{base}{target} = T{{}};",
                f"{base}for (unsigned int i = 0; i < LatSet::q; ++i) {{",
                f"{base}  {target} += cell[i];",
                f"{base}}}",
            ]

        if call.name == "moment1":
            rho = self._expr(call.args[1])
            return [
                f"{base}{target}.clear();",
                f"{base}for (unsigned int i = 0; i < LatSet::q; ++i) {{",
                f"{base}  {target} += latset::c<LatSet>(i) * cell[i];",
                f"{base}}}",
                f"{base}{target} /= {rho};",
            ]

        # moment1_force: the internal rho is a local of a block, so it can
        # never clobber a variable of the same name from the .op.
        force = self._expr(call.args[1])
        return [
            f"{base}{{",
            f"{base}  T rho_value{{}};",
            f"{base}  {target}.clear();",
            f"{base}  for (unsigned int i = 0; i < LatSet::q; ++i) {{",
            f"{base}    rho_value += cell[i];",
            f"{base}    {target} += latset::c<LatSet>(i) * cell[i];",
            f"{base}  }}",
            f"{base}  {target} += {force} * T{{0.5}};",
            f"{base}  {target} /= rho_value;",
            f"{base}}}",
        ]

    # -- expressions ------------------------------------------------------

    @staticmethod
    def _index_code(node: _Expr) -> str:
        if isinstance(node, _Literal):
            return node.text
        assert isinstance(node, _Name)
        return node.name

    def _expr(self, node: _Expr, ctx: Optional[_Feq2Ctx] = None) -> str:
        if isinstance(node, _Literal):
            return f"T{{{node.text}}}"
        if isinstance(node, _Name):
            if node.name == "omega":
                return "omega"
            if node.name == "fomega":
                return "fomega"
            return node.name
        if isinstance(node, _Index):
            if node.base == "f":
                return f"cell[{self._index_code(node.index)}]"
            return f"{node.base}[{self._index_code(node.index)}]"
        if isinstance(node, _Unary):
            return f"-({self._expr(node.operand, ctx)})"
        if isinstance(node, _Binary):
            return f"({self._expr(node.left, ctx)} {node.op} {self._expr(node.right, ctx)})"
        if isinstance(node, _Call):
            return self._call_expr(node, ctx)
        raise RuntimeError("internal: unknown expression node")  # pragma: no cover

    def _call_expr(self, node: _Call, ctx: Optional[_Feq2Ctx]) -> str:
        name = node.name
        if name == "feq2":
            if ctx is None:
                raise RuntimeError("internal: feq2() needs an emission context")
            rho = self._expr(node.args[0], ctx)
            velocity = self._expr(node.args[1], ctx)
            index = self._index_code(node.args[2])
            u2 = ctx.u2_name(velocity)
            uc = ctx.uc_name(velocity, index)
            # Same association as equilibrium::SecondOrder::apply / Order2.
            return (
                f"(latset::w<LatSet>({index}) * {rho} * "
                f"(T{{1}} + LatSet::InvCs2 * {uc} + {uc} * {uc} * T{{0.5}} * "
                f"LatSet::InvCs4 - LatSet::InvCs2 * {u2} * T{{0.5}}))"
            )
        if name == "w":
            return f"latset::w<LatSet>({self._index_code(node.args[0])})"
        if name == "c":
            return f"latset::c<LatSet>({self._index_code(node.args[0])})"
        if name == "dot":
            return f"({self._expr(node.args[0], ctx)} * {self._expr(node.args[1], ctx)})"
        if name == "norm2":
            return f"({self._expr(node.args[0], ctx)}).getnorm2()"
        if name == "field":
            return self._field_expr(node.args[0].name)
        if name in _MATH_SCALAR_FUNCS:
            args = ", ".join(self._expr(arg, ctx) for arg in node.args)
            return f"std::{name}({args})"
        if name == "force":
            return self._field_expr("FORCE")
        raise RuntimeError(f"internal: cannot emit builtin '{name}'")  # pragma: no cover

    @staticmethod
    def _field_expr(name: str) -> str:
        return f"cell.template get<{model.FIELD_TYPES[name]}>()"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def parse_operator(
    path, q: Optional[int] = None, dim: Optional[int] = None
) -> model.OperatorSpec:
    """Parse and type-check a ``.op`` file.

    Returns a ``model.OperatorSpec`` whose ``body`` is the checked AST and
    whose ``needs`` is the field closure (``POP`` plus ``FORCE`` when the
    operator reads the force field).  Errors are ``model.CasegenError`` with
    ``<path>:<line>`` in the message.  When ``q``/``dim`` are given, literal
    lattice and component indices are range-checked against them.
    """
    path = os.fspath(path)
    try:
        raw = Path(path).read_bytes()
    except OSError as exc:
        raise model.CasegenError(f"cannot read operator file '{path}': {exc.strerror or exc}")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise model.CasegenError(f"{path}: not valid UTF-8: {exc}")

    parsed = _Parser(_lex(text, path), path).parse()
    needs = _Checker(path, parsed.variables, q=q, dim=dim).check(parsed.body)
    return model.OperatorSpec(
        name=parsed.name,
        path=path,
        variables=dict(parsed.variables),
        body=parsed.body,
        needs=needs,
        source_hash=hashlib.sha256(raw).hexdigest(),
    )


def render_operator(op: model.OperatorSpec) -> str:
    """Render a parsed operator as a C++ header.

    The result is a self-contained header defining ``struct PascalName`` with
    ``__any__ static void apply(CELL& cell)``.  The spec is re-checked so a
    hand-built AST cannot emit unchecked code; ``op.needs`` is refreshed from
    that pass.
    """
    if not isinstance(op, model.OperatorSpec):
        raise TypeError("render_operator() expects a model.OperatorSpec")
    op.needs = _Checker(op.path, op.variables).check(op.body)
    return _Renderer(op).render()
