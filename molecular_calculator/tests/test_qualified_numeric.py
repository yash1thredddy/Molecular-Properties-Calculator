"""Tests for detecting and (optionally) stripping qualified numeric columns.

Assay exports (StarDrop, etc.) often store censored measurements as strings
with a comparison operator: ">12", "<0.5", ">=100". pandas then keeps the whole
column as text, so it is excluded from numeric analysis. Removing the operator
changes the meaning of the data, so the app only DETECTS such columns
automatically and leaves the actual stripping to an explicit user choice.

Detection and stripping are column-agnostic: they apply to any column, not just
a specific one like MIC.
"""

import numpy as np
import pandas as pd

from molecular_calculator.ui.pages.gmm_analysis import (
    find_qualified_columns,
    strip_qualifiers,
    _strip_qualifier,
)


def test_strip_qualifier_variants():
    assert _strip_qualifier(">12") == 12.0
    assert _strip_qualifier("< 0.5") == 0.5
    assert _strip_qualifier(">=100") == 100.0
    assert _strip_qualifier("<= 3.2") == 3.2
    assert _strip_qualifier("12") == 12.0
    assert _strip_qualifier("1.5e-3") == 0.0015
    assert _strip_qualifier(7) == 7.0
    assert _strip_qualifier(7.5) == 7.5
    assert _strip_qualifier("abc") is None
    assert _strip_qualifier("100-200") is None  # range, not a single number
    assert _strip_qualifier(None) is None
    assert _strip_qualifier(True) is None  # booleans are not measurements


def test_find_qualified_columns_detects_without_mutating():
    df = pd.DataFrame({
        "SMILES": ["CCO", "c1ccccc1", "CCN"],
        "MIC min (µM)": [">12", "5", ">8"],
    })
    report = find_qualified_columns(df)

    assert report == {"MIC min (µM)": 2}
    # detection must NOT change the data
    assert df["MIC min (µM)"].tolist() == [">12", "5", ">8"]


def test_find_qualified_columns_is_column_agnostic():
    """Any column with comparison qualifiers is flagged, not just MIC."""
    df = pd.DataFrame({
        "IC50 (nM)": [">100", "12", "<5"],
        "Solubility": ["<0.1", "2.3", "1.0"],
        "Name": ["a", "b", "c"],          # text -> ignored
        "MW": [120.0, 240.0, 88.0],        # already numeric -> ignored
    })
    report = find_qualified_columns(df)
    assert set(report) == {"IC50 (nM)", "Solubility"}
    assert report["IC50 (nM)"] == 2
    assert report["Solubility"] == 1


def test_clean_numeric_string_column_not_flagged():
    """No comparison symbols => nothing to ask the user about."""
    df = pd.DataFrame({"vals": ["1", "2", "3"]})
    assert find_qualified_columns(df) == {}


def test_text_column_not_flagged():
    df = pd.DataFrame({"Status": ["active", "inactive", "active"]})
    assert find_qualified_columns(df) == {}


def test_strip_qualifiers_converts_selected_columns():
    df = pd.DataFrame({
        "MIC min (µM)": [">10", ">12", "3.5", "20", "<0.5"],
        "Name": ["a", "b", "c", "d", "e"],
    })
    out = strip_qualifiers(df, ["MIC min (µM)"])

    assert out["MIC min (µM)"].tolist() == [10.0, 12.0, 3.5, 20.0, 0.5]
    assert np.issubdtype(out["MIC min (µM)"].dtype, np.number)
    # untouched column stays as-is
    assert out["Name"].tolist() == ["a", "b", "c", "d", "e"]
    # original DataFrame not mutated
    assert df["MIC min (µM)"].tolist() == [">10", ">12", "3.5", "20", "<0.5"]


def test_strip_qualifiers_ignores_unknown_columns():
    df = pd.DataFrame({"MIC": [">12", "5"]})
    out = strip_qualifiers(df, ["does_not_exist"])
    assert out["MIC"].tolist() == [">12", "5"]
