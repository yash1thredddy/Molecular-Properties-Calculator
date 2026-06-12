import numpy as np
import pandas as pd
from molecular_calculator.services import formula_presets as fp
from molecular_calculator.services import formula_engine as fe


def test_presets_have_required_fields():
    for p in fp.PRESETS:
        assert p.name and p.expression and p.group and p.citation is not None


def test_percent_inhibition_alias_detection():
    df = pd.DataFrame({"%inhibition": [98.0], "Molecular_Weight": [180.0]})
    assert fp.detect_activity_column(df, "percent_inhibition") == "%inhibition"


def test_resolve_placeholder_to_column():
    expr = fp.resolve_placeholders(
        "([activity:percent_inhibition]/100) / ([Molecular_Weight]/1000)",
        {"percent_inhibition": "%inhibition"})
    assert "[%inhibition]" in expr and "[activity:" not in expr


def test_pei_value_correct():
    expr = fp.resolve_placeholders(
        "([activity:percent_inhibition]/100) / ([Molecular_Weight]/1000)",
        {"percent_inhibition": "Percent_Inhibition"})
    df = pd.DataFrame({"Percent_Inhibition": [98.0], "Molecular_Weight": [180.0]})
    cf = fe.compile(expr, list(df.columns))
    assert cf.error is None, cf.error
    out, _ = fe.evaluate(cf, df)
    # (0.98) / (0.18) = 5.444...
    assert abs(out.iloc[0] - 5.4444) < 0.01


def test_lipinski_violation_count():
    df = pd.DataFrame({"Molecular_Weight": [600.0], "LogP": [6.0],
                       "HB_Donors": [2], "HB_Acceptors": [3]})
    p = next(x for x in fp.PRESETS if x.name == "Lipinski_Violation_Count")
    cf = fe.compile(p.expression, list(df.columns))
    out, _ = fe.evaluate(cf, df)
    assert out.iloc[0] == 2  # MW>500 and LogP>5


def test_raw_concentration_wrapped_to_pscale():
    # mapping a raw IC50 (nM) column to a pIC50 placeholder must inject -log10(x/1e9)
    expr = fp.resolve_placeholders(
        "[activity:pIC50] - [LogP]",
        {"pIC50": "IC50_nM"},
        raw_concentration={"pIC50": True})
    assert "log10" in expr and "1000000000" in expr
