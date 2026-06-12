"""Safe user-defined formula engine (deterministic core).

Compiles a user formula once and evaluates it per DataFrame row with a
positive AST allowlist, lazy if()/switch(), name-bound column aliases, and
correct NaN propagation. No arbitrary code execution.
"""
from __future__ import annotations

import ast
import io
import math
import re
import tokenize
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from simpleeval import SimpleEval
import simpleeval as _se

MAX_FORMULA_LEN = 500
MAX_STRING_LEN = 1000
ALIAS_PREFIX = "__c"  # reserved; raw formulas containing it are rejected


class FormulaError(Exception):
    """Raised on compile/validation failure (user-facing message)."""


class _PropagateNaN(Exception):
    """Internal: an unguarded referenced value was missing this row."""


ALLOWED_CONSTANTS = {"pi": math.pi, "e": math.e}


def _is_missing(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))


def _require(x):
    """Unguarded column reference: raise if missing so the row -> NaN."""
    if _is_missing(x):
        raise _PropagateNaN()
    return x


def _coalesce(x, fallback):
    return fallback if _is_missing(x) else x


def _isempty(x) -> bool:
    return _is_missing(x)


def _log(x, base=math.e):
    return math.log(x, base)


def _substring(s, start, end=None):
    out = str(s)[start:end]
    if len(out) > MAX_STRING_LEN:
        raise FormulaError("string result too long")
    return out


def _concat(*parts):
    out = "".join(str(p) for p in parts)
    if len(out) > MAX_STRING_LEN:
        raise FormulaError("string result too long")
    return out


ALLOWED_FUNCTIONS = {
    "abs": abs, "round": round, "ceil": math.ceil, "floor": math.floor,
    "sqrt": math.sqrt, "exp": math.exp, "ln": math.log, "log10": math.log10,
    "log": _log, "pow": pow, "min": min, "max": max,
    "mod": (lambda a, b: a % b),
    "isempty": _isempty, "coalesce": _coalesce,
    "concat": _concat, "contains": (lambda s, sub: str(sub) in str(s)),
    "lower": (lambda s: str(s).lower()), "upper": (lambda s: str(s).upper()),
    "len": (lambda s: len(str(s))), "substring": _substring,
    "not_": (lambda x: not x),
    # internal helpers (not user-callable by name; injected by the transform)
    "_require": _require,
}

# Function names a user may type (excludes internal helpers).
# __if__ and __switch__ are the pre-parse placeholders for the DSL keywords;
# they are listed here so _validate_ast accepts them before the conditional
# transformer rewrites them to IfExp nodes.
USER_FUNCTIONS = frozenset(ALLOWED_FUNCTIONS) - {"_require"} | {"__if__", "__switch__"}


# ---------------------------------------------------------------------------
# Task 2: AST allowlist validator
# ---------------------------------------------------------------------------

_ALLOWED_NODES = (
    ast.Expression, ast.Constant, ast.Name, ast.Load,
    ast.BinOp, ast.UnaryOp, ast.BoolOp, ast.Compare, ast.IfExp, ast.Call,
    ast.Add, ast.Sub, ast.Mult, ast.Div, ast.FloorDiv, ast.Mod, ast.Pow,
    ast.USub, ast.UAdd, ast.Not, ast.And, ast.Or,
    ast.Eq, ast.NotEq, ast.Lt, ast.LtE, ast.Gt, ast.GtE,
)


_DSL_KEYWORD_MAP = {"if": "__if__", "switch": "__switch__"}


def _substitute_dsl_keywords(src: str) -> str:
    """Replace DSL `if(` and `switch(` with safe placeholder names using tokenize.

    Only NAME tokens immediately followed by an OP `(` token are replaced.
    STRING tokens are never altered, so literals like "if(x>0)" are preserved.
    """
    tokens = []
    try:
        token_list = list(tokenize.generate_tokens(io.StringIO(src).readline))
    except tokenize.TokenError:
        # If tokenization fails (e.g. unclosed string), fall through to ast.parse
        # which will produce a clearer SyntaxError.
        return src

    for i, tok in enumerate(token_list):
        tok_type, tok_string, tok_start, tok_end, tok_line = tok
        if (
            tok_type == tokenize.NAME
            and tok_string in _DSL_KEYWORD_MAP
            and i + 1 < len(token_list)
            and token_list[i + 1][0] == tokenize.OP
            and token_list[i + 1][1] == "("
        ):
            tokens.append(tokenize.TokenInfo(
                type=tok_type,
                string=_DSL_KEYWORD_MAP[tok_string],
                start=tok_start,
                end=tok_end,
                line=tok_line,
            ))
        else:
            tokens.append(tok)
    return tokenize.untokenize(tokens)


def _parse(src: str) -> ast.Expression:
    if len(src) > MAX_FORMULA_LEN:
        raise FormulaError(f"Formula too long (max {MAX_FORMULA_LEN} chars).")
    # Pre-substitute DSL keywords that are reserved in Python so ast.parse succeeds.
    # Uses tokenize to avoid corrupting string literals that contain "if(" or "switch(".
    src = _substitute_dsl_keywords(src)
    try:
        return ast.parse(src, mode="eval")
    except SyntaxError as exc:
        raise FormulaError(f"Syntax error: {exc.msg}")


def _validate_ast(tree: ast.AST, allowed_names: Optional[set] = None) -> None:
    """Reject any node/operator/call outside the positive allowlist."""
    for node in ast.walk(tree):
        if not isinstance(node, _ALLOWED_NODES):
            raise FormulaError(
                f"Disallowed expression element: {type(node).__name__}")
        if isinstance(node, ast.Call):
            if node.keywords or any(isinstance(a, ast.Starred) for a in node.args):
                raise FormulaError("Function keyword/starred args are not allowed.")
            if not isinstance(node.func, ast.Name):
                raise FormulaError("Only direct function calls are allowed.")
            if node.func.id not in USER_FUNCTIONS and node.func.id != "_require":
                raise FormulaError(f"Unknown function: {node.func.id}()")
        if isinstance(node, ast.Name) and allowed_names is not None:
            # callable names are validated via Call above
            is_call_target = False
            if node.id in USER_FUNCTIONS or node.id == "_require":
                is_call_target = True
            if not is_call_target and node.id not in allowed_names:
                raise FormulaError(f"Unknown name: {node.id}")


# ---------------------------------------------------------------------------
# Task 3: Bracket → name-bound alias transform
# ---------------------------------------------------------------------------

_BRACKET_RE = re.compile(r"\[([^\[\]]+)\]")


def _resolve_brackets(src: str, columns: list[str]) -> tuple[str, dict[str, str]]:
    """Replace [Column Name] with a safe alias; return (rewritten, alias->name)."""
    if ALIAS_PREFIX in src:
        raise FormulaError(f"Formula may not contain the reserved token {ALIAS_PREFIX!r}.")
    if len(columns) != len(set(columns)):
        raise FormulaError("Duplicate column names — rename them before using formulas.")
    colset = set(columns)
    name_to_alias: dict[str, str] = {}
    alias_map: dict[str, str] = {}

    def repl(m: re.Match) -> str:
        name = m.group(1).strip()
        if name not in colset:
            raise FormulaError(f"Unknown column: [{name}]")
        if name not in name_to_alias:
            alias = f"{ALIAS_PREFIX}{len(name_to_alias)}__"
            name_to_alias[name] = alias
            alias_map[alias] = name
        return name_to_alias[name]

    rewritten = _BRACKET_RE.sub(repl, src)
    return rewritten, alias_map


# ---------------------------------------------------------------------------
# Task 4: if()/switch() → lazy IfExp transform + caret-to-pow
# ---------------------------------------------------------------------------

class _ConditionalTransformer(ast.NodeTransformer):
    def visit_Call(self, node: ast.Call) -> ast.AST:
        self.generic_visit(node)
        # Match the placeholder names that _parse substituted for the DSL keywords.
        if isinstance(node.func, ast.Name) and node.func.id == "__if__":
            if len(node.args) != 3:
                raise FormulaError("if() takes exactly 3 arguments: if(cond, then, else)")
            cond, then, els = node.args
            return ast.copy_location(ast.IfExp(test=cond, body=then, orelse=els), node)
        if isinstance(node.func, ast.Name) and node.func.id == "__switch__":
            args = node.args
            if len(args) < 4 or (len(args) % 2) != 0:
                raise FormulaError(
                    "switch() takes value, then case/result pairs, then a default: "
                    "switch(v, c1, r1, ..., default)")
            value, *rest = args
            default = rest[-1]
            pairs = rest[:-1]
            expr: ast.AST = default
            for i in range(len(pairs) - 2, -1, -2):
                case, result = pairs[i], pairs[i + 1]
                test = ast.Compare(left=value, ops=[ast.Eq()], comparators=[case])
                expr = ast.IfExp(test=test, body=result, orelse=expr)
            return ast.copy_location(ast.fix_missing_locations(expr), node)
        return node


def _transform_conditionals(tree: ast.Expression) -> ast.Expression:
    tree = _ConditionalTransformer().visit(tree)
    ast.fix_missing_locations(tree)
    return tree


class _CaretTransformer(ast.NodeTransformer):
    """In our DSL, `^` means power. Python parses it as BitXor — rewrite to Pow."""
    def visit_BinOp(self, node: ast.BinOp) -> ast.AST:
        self.generic_visit(node)
        if isinstance(node.op, ast.BitXor):
            node.op = ast.Pow()
        return node


def _caret_to_pow(tree: ast.Expression) -> ast.Expression:
    tree = _CaretTransformer().visit(tree)
    ast.fix_missing_locations(tree)
    return tree


# ---------------------------------------------------------------------------
# Task 5: NaN-guard transform — wrap unguarded refs with _require
# ---------------------------------------------------------------------------

_GUARDING_FUNCS = {"isempty", "coalesce"}


class _NaNGuardTransformer(ast.NodeTransformer):
    def __init__(self, aliases: set[str]):
        self._aliases = aliases
        self._guard_depth = 0

    def visit_Call(self, node: ast.Call) -> ast.AST:
        guarding = isinstance(node.func, ast.Name) and node.func.id in _GUARDING_FUNCS
        node.func = self.visit(node.func)
        new_args = []
        if guarding:
            # Increment depth around ALL arguments so that any alias nested
            # inside a guarding call (e.g. coalesce(sqrt([A]), 0)) is left raw.
            self._guard_depth += 1
            for arg in node.args:
                new_args.append(self.visit(arg))
            self._guard_depth -= 1
        else:
            for arg in node.args:
                new_args.append(self.visit(arg))
        node.args = new_args
        return node

    def visit_Name(self, node: ast.Name) -> ast.AST:
        if node.id in self._aliases and self._guard_depth == 0:
            return ast.copy_location(
                ast.Call(func=ast.Name(id="_require", ctx=ast.Load()),
                         args=[node], keywords=[]),
                node)
        return node


def _wrap_nan_guards(tree: ast.Expression, aliases: set[str]) -> ast.Expression:
    tree = _NaNGuardTransformer(aliases).visit(tree)
    ast.fix_missing_locations(tree)
    return tree


# ---------------------------------------------------------------------------
# Task 6: CompiledFormula + compile()
# ---------------------------------------------------------------------------

@dataclass
class CompiledFormula:
    source: str
    code: str = ""                       # final, unparsed, validated expression
    alias_map: dict = field(default_factory=dict)   # alias -> real column name
    referenced_cols: set = field(default_factory=set)
    expected_columns: list = field(default_factory=list)
    error: Optional[str] = None


_RESERVED_TOKENS = ("__if__", "__switch__", "_require", ALIAS_PREFIX)


def compile(formula: str, columns: list[str]) -> CompiledFormula:
    """Single compiler path: bracket-resolve -> validate -> transform -> unparse.

    Never raises FormulaError; returns it on `.error` so the UI can show it.
    """
    cf = CompiledFormula(source=formula, expected_columns=list(columns))
    try:
        # Reject any raw formula that contains engine-internal reserved tokens.
        for reserved in _RESERVED_TOKENS:
            if reserved in formula:
                cf.error = (
                    f"Formula contains reserved token {reserved!r} — "
                    "use the public DSL functions instead."
                )
                return cf
        rewritten, alias_map = _resolve_brackets(formula, columns)
        aliases = set(alias_map)
        allowed_names = aliases | set(ALLOWED_CONSTANTS)
        tree = _parse(rewritten)
        tree = _caret_to_pow(tree)                        # `^` means power in the DSL
        _validate_ast(tree, allowed_names)               # user-facing allowlist
        tree = _transform_conditionals(tree)             # if()/switch() -> IfExp
        tree = _wrap_nan_guards(tree, aliases)           # unguarded refs -> _require
        ast.fix_missing_locations(tree)
        cf.alias_map = alias_map
        cf.referenced_cols = set(alias_map.values())
        cf.code = ast.unparse(tree)
    except FormulaError as exc:
        cf.error = str(exc)
    return cf


# ---------------------------------------------------------------------------
# Task 7: evaluate()
# ---------------------------------------------------------------------------

_se.MAX_POWER = 10 ** 6  # DoS guard for the 1 GB instance


def _new_evaluator() -> SimpleEval:
    s = SimpleEval(functions=dict(ALLOWED_FUNCTIONS), names={})
    return s


def _scalar(v):
    """Normalize any missing form (np.nan/pd.NA/None) to float nan."""
    if pd.isna(v):
        return float("nan")
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    return v


def evaluate(compiled: CompiledFormula, df: pd.DataFrame) -> tuple[pd.Series, int]:
    """Evaluate a compiled formula per row. Returns (result Series, failed rows)."""
    if compiled.error:
        raise FormulaError(compiled.error)
    missing = [c for c in compiled.referenced_cols if c not in df.columns]
    if missing:
        raise FormulaError(
            f"Column(s) no longer present: {', '.join(missing)} — recompile.")
    if len(df) == 0:
        return pd.Series([], dtype="float64"), 0

    s = _new_evaluator()
    # alias -> the actual column Series, in row order
    alias_cols = {alias: df[name].to_numpy(dtype=object)
                  for alias, name in compiled.alias_map.items()}
    results = []
    failed = 0
    n = len(df)
    for i in range(n):
        names = {alias: _scalar(col[i]) for alias, col in alias_cols.items()}
        names.update(ALLOWED_CONSTANTS)
        s.names = names
        try:
            # compiled.code is the engine-generated, AST-validated string — not raw user input.
            results.append(s.eval(compiled.code))
        except _PropagateNaN:
            results.append(float("nan"))
            failed += 1
        except Exception:
            results.append(float("nan"))
            failed += 1
    return pd.Series(results, index=df.index), failed
