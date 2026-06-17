"""Tests for formula_engine — Tasks 1 through 8."""
import ast

import numpy as np
import pandas as pd
import pytest

from molecular_calculator.services import formula_engine as fe


# ---------------------------------------------------------------------------
# Task 1: whitelist presence
# ---------------------------------------------------------------------------

def test_whitelists_exist():
    assert "pi" in fe.ALLOWED_CONSTANTS and "e" in fe.ALLOWED_CONSTANTS
    for name in ["abs", "round", "ceil", "floor", "sqrt", "exp", "ln", "log10",
                 "log", "pow", "min", "max", "mod", "isempty", "coalesce",
                 "concat", "contains", "lower", "upper", "len", "substring"]:
        assert name in fe.ALLOWED_FUNCTIONS, name


# ---------------------------------------------------------------------------
# Task 2: AST allowlist validator
# ---------------------------------------------------------------------------

def test_validator_allows_safe_nodes():
    fe._validate_ast(fe._parse("1 + 2 * 3 - (4 / 2) ** 2 % 3"))
    fe._validate_ast(fe._parse("a > 5 and (b <= 3 or not c)"))
    fe._validate_ast(fe._parse("sqrt(a) + log10(b)"))


@pytest.mark.parametrize("src", [
    "__import__('os')",
    "a.__class__",
    "().__class__.__bases__",
    "data[0]",
    "[x for x in range(3)]",
    "{1: 2}",
    "(lambda: 1)()",
    "f'{a}'",
    "foo(**a)",
    "min(*a)",
    "(x := 5)",
    "evil()",            # call to non-whitelisted name
])
def test_validator_rejects_unsafe(src):
    with pytest.raises(fe.FormulaError):
        fe._validate_ast(fe._parse(src))


# ---------------------------------------------------------------------------
# Task 3: bracket resolution
# ---------------------------------------------------------------------------

def test_bracket_resolution_and_aliases():
    src = "[Molecular Weight] / 1000 + [LogP]"
    rewritten, alias_map = fe._resolve_brackets(src, ["Molecular Weight", "LogP", "TPSA"])
    # aliases are safe identifiers mapped to real names
    assert set(alias_map.values()) == {"Molecular Weight", "LogP"}
    assert "[" not in rewritten and "]" not in rewritten
    for alias in alias_map:
        assert alias.startswith(fe.ALIAS_PREFIX)


def test_bracket_unknown_column():
    with pytest.raises(fe.FormulaError):
        fe._resolve_brackets("[Nope]", ["LogP"])


def test_reserved_prefix_rejected():
    with pytest.raises(fe.FormulaError):
        fe._resolve_brackets("__c0__ + 1", ["LogP"])


def test_duplicate_columns_rejected():
    with pytest.raises(fe.FormulaError):
        fe._resolve_brackets("[LogP]", ["LogP", "LogP"])


# ---------------------------------------------------------------------------
# Task 4: conditionals + caret transforms
# ---------------------------------------------------------------------------

def _roundtrip(src):
    tree = fe._parse(src)
    tree = fe._transform_conditionals(tree)
    return ast.unparse(tree)


def test_if_becomes_ternary():
    out = _roundtrip("if(a > 1, 10, 20)")
    # ast.unparse renders IfExp as "body if test else orelse" (no parens needed).
    # Plan had "(" in out but ast.unparse 3.8+ omits unnecessary parens — fixed.
    assert "if" in out and "else" in out  # IfExp unparsed


def test_switch_nested_ternary():
    out = _roundtrip("switch(a, 1, 100, 2, 200, 999)")
    # nested IfExp: a==1 ? 100 : (a==2 ? 200 : 999)
    assert out.count("if") >= 2 and "999" in out


def test_if_wrong_arity():
    with pytest.raises(fe.FormulaError):
        _roundtrip("if(a, 1)")


def test_caret_is_power():
    tree = fe._caret_to_pow(fe._parse("a ^ 0.3"))
    assert "**" in ast.unparse(tree)  # ^ rewritten to power, not bitwise-xor


# ---------------------------------------------------------------------------
# Task 5: NaN-guard transform
# ---------------------------------------------------------------------------

def test_guard_wraps_unguarded_only():
    tree = fe._parse("__c0__ + coalesce(__c1__, 0) + isempty(__c2__)")
    tree = fe._wrap_nan_guards(tree, {"__c0__", "__c1__", "__c2__"})
    out = ast.unparse(tree)
    assert "_require(__c0__)" in out          # unguarded -> wrapped
    assert "_require(__c1__)" not in out       # inside coalesce -> raw
    assert "_require(__c2__)" not in out       # inside isempty -> raw


# ---------------------------------------------------------------------------
# Task 6: compile()
# ---------------------------------------------------------------------------

def test_compile_success_and_metadata():
    cf = fe.compile("[A] + coalesce([B], 0)", ["A", "B", "C"])
    assert cf.error is None
    assert cf.expected_columns == ["A", "B", "C"]
    assert cf.referenced_cols == {"A", "B"}


def test_compile_reports_error_not_raises():
    cf = fe.compile("[A] + nope()", ["A"])
    assert cf.error is not None and "nope" in cf.error


# ---------------------------------------------------------------------------
# Task 7: evaluate()
# ---------------------------------------------------------------------------

def _ev(formula, df):
    cf = fe.compile(formula, list(df.columns))
    assert cf.error is None, cf.error
    return fe.evaluate(cf, df)


def test_evaluate_basic_math():
    df = pd.DataFrame({"MW": [100.0, 200.0]})
    out, report = _ev("[MW] / 1000", df)
    assert list(out) == [0.1, 0.2] and report.errored == 0 and report.blank == 0


def test_nan_propagates_through_comparison():
    df = pd.DataFrame({"MW": [600.0, np.nan]})
    out, report = _ev("if([MW] > 500, 1, 0)", df)
    assert out.iloc[0] == 1
    assert pd.isna(out.iloc[1]) and report.blank == 1 and report.errored == 0


def test_if_laziness_untaken_branch_safe():
    df = pd.DataFrame({"x": [0.0], "y": [5.0]})
    out, report = _ev("if([x] != 0, [y] / [x], 0)", df)
    assert out.iloc[0] == 0 and report.errored == 0  # no ZeroDivision, branch not taken


def test_coalesce_escape_hatch():
    df = pd.DataFrame({"p": [np.nan]})
    out, report = _ev("coalesce([p], 0) + 1", df)
    assert out.iloc[0] == 1 and report.errored == 0 and report.blank == 0


def test_pd_na_sanitized():
    df = pd.DataFrame({"v": pd.array([1, pd.NA], dtype="Int64")})
    out, report = _ev("[v] * 2", df)
    assert out.iloc[0] == 2 and pd.isna(out.iloc[1])


def test_empty_df():
    df = pd.DataFrame({"A": []})
    out, report = _ev("[A] + 1", df)
    assert len(out) == 0 and report.errored == 0 and report.blank == 0


def test_evaluate_missing_column_errors():
    cf = fe.compile("[A] + 1", ["A"])
    with pytest.raises(fe.FormulaError):
        fe.evaluate(cf, pd.DataFrame({"B": [1]}))


# ---------------------------------------------------------------------------
# Task 8: security-escape regression suite
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("formula", [
    "__import__('os').system('echo hi')",
    "[A].__class__",
    "().__class__.__bases__[0].__subclasses__()",
    "[A][0]",
    "[x for x in [1]]",
    "(lambda: 1)()",
    "f'{[A]}'",
    "open('x')",
    "eval('1')",
    "[A] ** 999999999",   # power DoS guarded by MAX_POWER at eval
])
def test_security_blocked(formula):
    cf = fe.compile(formula, ["A"])
    if cf.error is None:
        # if it compiled, evaluation must not escape — should fail-safe to NaN
        out, report = fe.evaluate(cf, pd.DataFrame({"A": [2.0]}))
        assert pd.isna(out.iloc[0])
    else:
        assert cf.error  # rejected at compile


# ---------------------------------------------------------------------------
# New tests for review fixes
# ---------------------------------------------------------------------------

def test_string_literal_with_if_preserved():
    """HIGH 1: string literals containing 'if(' must not be mangled."""
    df = pd.DataFrame({"Name": ["zif(x)q"]})
    # contains([Name], "if(x)") should return True because "if(x)" is in "zif(x)q"
    cf = fe.compile('contains([Name], "if(x)")', list(df.columns))
    assert cf.error is None, cf.error
    out, report = fe.evaluate(cf, df)
    assert out.iloc[0] == True and report.errored == 0


def test_reserved_namespace_rejected():
    """MED 2: formulas containing engine-internal tokens must be rejected at compile."""
    cf = fe.compile("__if__(1, 2, 3)", ["A"])
    assert cf.error is not None


def test_coalesce_with_subexpr_nan_returns_fallback():
    """MED 3: coalesce(sqrt([A]), 0) with A=NaN must return 0, not propagate NaN."""
    df = pd.DataFrame({"A": [np.nan]})
    cf = fe.compile("coalesce(sqrt([A]), 0)", list(df.columns))
    assert cf.error is None, cf.error
    out, report = fe.evaluate(cf, df)
    assert out.iloc[0] == 0.0 and report.errored == 0 and report.blank == 0


# ---------------------------------------------------------------------------
# Task 1A: Exponent guard + result bounding (#1, #2)
# ---------------------------------------------------------------------------

def _df(n=3):
    return pd.DataFrame({"A": [1.0] * n})


def test_huge_exponent_does_not_crash_series():
    cf = fe.compile("10 ** 999999", ["A"])
    series, report = fe.evaluate(cf, _df())
    # No OverflowError; offending rows are NaN and counted as errored.
    assert series.isna().all()
    assert report.errored == 3


def test_large_operand_power_is_rejected_fast():
    # 999999 exceeds MAX_POWER, so simpleeval refuses it (no multi-second compute).
    cf = fe.compile("999999 ** 999999", ["A"])
    series, report = fe.evaluate(cf, _df())
    assert series.isna().all()
    assert report.errored == 3


def test_normal_power_still_works():
    cf = fe.compile("[A] ** 2", ["A"])
    series, report = fe.evaluate(cf, pd.DataFrame({"A": [2.0, 3.0]}))
    assert series.tolist() == [4.0, 9.0]
    assert report.errored == 0


# ---------------------------------------------------------------------------
# Task 1B: Whitespace column disambiguation (#5) + reserved-token fix (#6)
# ---------------------------------------------------------------------------

def test_whitespace_columns_are_disambiguated():
    df = pd.DataFrame({"A": [10.0], "A ": [20.0]})
    # exact matches resolve to the right column
    s_a, _ = fe.evaluate(fe.compile("[A]", ["A", "A "]), df)
    s_b, _ = fe.evaluate(fe.compile("[A ]", ["A", "A "]), df)
    assert s_a.tolist() == [10.0]
    assert s_b.tolist() == [20.0]


def test_ambiguous_stripped_reference_errors():
    cf = fe.compile("[ A ]", ["A", "A "])  # neither exact; strips to 'A' -> ambiguous
    assert cf.error is not None and "mbiguous" in cf.error


def test_column_named_with_double_underscore_c_is_allowed():
    df = pd.DataFrame({"logP__corrected": [1.5]})
    cf = fe.compile("[logP__corrected] * 2", ["logP__corrected"])
    assert cf.error is None
    series, _ = fe.evaluate(cf, df)
    assert series.tolist() == [3.0]


# ---------------------------------------------------------------------------
# Task 1C: Separate "blank" vs "errored" rows (#7)
# ---------------------------------------------------------------------------

def test_eval_report_distinguishes_blank_from_error():
    df = pd.DataFrame({"A": [1.0, float("nan"), 4.0]})
    # division by zero is a real error; the NaN row is a blank input
    cf = fe.compile("1 / ([A] - [A])", ["A"])  # always 1/0 where A present
    series, report = fe.evaluate(cf, df)
    assert report.blank == 1     # the NaN-input row
    assert report.errored == 2   # the two 1/0 rows
    assert report.first_error is not None
