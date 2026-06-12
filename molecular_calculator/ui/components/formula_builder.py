"""Shared Streamlit formula-builder component.

Renders an expander that lets users define custom DataFrame columns via a
safe formula DSL, apply presets, preview results live, and export/import
formula sets as JSON.

Usage::

    from molecular_calculator.ui.components.formula_builder import render_formula_builder

    # Call at any point where `df` is the working DataFrame.
    df = render_formula_builder(df, page_key="batch")
    # `df` now includes any previously applied formula columns.

Public exports:
    render_formula_builder(df, page_key) -> pd.DataFrame
    _schema_hash(df) -> str   # used by page integrations for stale-state detection
"""
from __future__ import annotations

import hashlib
import json

import pandas as pd
import streamlit as st

from molecular_calculator.services import formula_engine as fe
from molecular_calculator.services import formula_presets as fp


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def _schema_hash(df: pd.DataFrame) -> str:
    """Return an MD5 hex digest of the column names (order-sensitive).

    Used by page integrations to detect when the working DataFrame schema
    has changed so they can invalidate stale downstream state.
    """
    return hashlib.md5("|".join(map(str, df.columns)).encode()).hexdigest()


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _apply_saved_formulas(df: pd.DataFrame, defs: list[dict]) -> pd.DataFrame:
    """Re-evaluate all saved formula definitions against `df`.

    Formulas that fail to compile or evaluate (e.g. referenced column was
    removed) are silently skipped so the page remains usable.
    """
    out = df.copy()
    for d in defs:
        cf = fe.compile(d["expression"], list(out.columns))
        if cf.error:
            continue
        series, _ = fe.evaluate(cf, out)
        out[d["name"]] = series
    return out


# ---------------------------------------------------------------------------
# Main component
# ---------------------------------------------------------------------------

def render_formula_builder(df: pd.DataFrame, page_key: str) -> pd.DataFrame:
    """Render the formula-builder expander and return the augmented DataFrame.

    Args:
        df: The current working DataFrame.
        page_key: A short, page-unique string used to namespace all Streamlit
            widget keys and session-state entries (e.g. ``"batch"``, ``"gmm"``).

    Returns:
        A copy of ``df`` with any applied formula columns appended.
    """
    defs_key = f"{page_key}_formula_defs"
    if defs_key not in st.session_state:
        st.session_state[defs_key] = []

    with st.expander("➕ Custom formula columns", expanded=False):
        cols = list(df.columns)

        # ------------------------------------------------------------------
        # Layout: formula input (left) | column picker + function ref (right)
        # ------------------------------------------------------------------
        c1, c2 = st.columns([2, 1])

        with c2:
            st.caption("Insert a column reference")
            picked = st.selectbox(
                "Columns",
                cols,
                key=f"{page_key}_fcol",
                label_visibility="collapsed",
            )
            st.code(f"[{picked}]", language=None)
            st.caption(
                "**Functions:** abs round ceil floor sqrt exp ln log10 log "
                "pow min max mod if() switch() isempty coalesce concat "
                "contains lower upper len substring"
            )

        with c1:
            # --------------------------------------------------------------
            # Preset gallery
            # --------------------------------------------------------------
            preset_names = ["(none)"] + [p.name for p in fp.PRESETS]
            chosen = st.selectbox(
                "Preset (optional)",
                preset_names,
                key=f"{page_key}_preset",
            )

            default_expr: str = ""
            default_name: str = ""

            if chosen != "(none)":
                preset = next(p for p in fp.PRESETS if p.name == chosen)
                expr = preset.expression
                ok = True
                mapping: dict[str, str] = {}
                raw_flags: dict[str, bool] = {}

                for key in preset.needs_activity:
                    detected = fp.detect_activity_column(df, key)
                    sel = st.selectbox(
                        f"Map activity:{key}",
                        ["(none)"] + cols,
                        index=(cols.index(detected) + 1) if (detected and detected in cols) else 0,
                        key=f"{page_key}_act_{key}",
                    )
                    if sel == "(none)":
                        ok = False
                    else:
                        mapping[key] = sel

                        # Task 16, Step 5 — unit radio for p-scale activity keys
                        if key in ("pIC50", "pKi"):
                            unit = st.radio(
                                f"'{sel}' is:",
                                ["already log-scale (pIC50/pKi)", "raw IC50/Ki in nM"],
                                key=f"{page_key}_unit_{key}",
                                horizontal=True,
                            )
                            raw_flags[key] = unit.startswith("raw")

                if ok:
                    try:
                        # Pass raw_concentration so nM columns get -log10 injection
                        expr = fp.resolve_placeholders(
                            expr, mapping, raw_concentration=raw_flags
                        )
                        default_expr, default_name = expr, preset.name
                    except KeyError:
                        st.warning("Map all activity columns to use this preset.")
                else:
                    st.info(
                        f"{preset.description}  \nNeeds: "
                        f"{', '.join(preset.needs_activity)}"
                    )

            new_name = st.text_input(
                "New column name",
                value=default_name,
                key=f"{page_key}_fname",
            )
            formula = st.text_area(
                "Formula",
                value=default_expr,
                key=f"{page_key}_fexpr",
                height=80,
            )

        # ------------------------------------------------------------------
        # Live validation + preview
        # ------------------------------------------------------------------
        if formula.strip():
            cf = fe.compile(formula, cols)
            if cf.error:
                st.error(f"⚠️ {cf.error}")
            else:
                preview_df = df.head(10)
                series, failed = fe.evaluate(cf, preview_df)
                st.success("✅ Valid")
                prev = preview_df.copy()
                prev[new_name or "result"] = series.values
                st.dataframe(
                    prev[[*cols[:2], new_name or "result"]].head(10),
                    use_container_width=True,
                )
                if failed:
                    st.caption(
                        f"{failed} of {len(preview_df)} preview rows blank "
                        "(missing values)."
                    )

                if new_name in cols:
                    st.caption(
                        "⚠️ Name already exists — choose another or it will overwrite."
                    )

                # Task 17 — large-dataset performance warning
                if len(df) > 25000:
                    st.caption(
                        f"Note: applying to {len(df):,} rows is single-threaded "
                        "and may take a few seconds."
                    )

                if st.button(
                    "Apply column",
                    key=f"{page_key}_apply",
                    disabled=(not new_name or len(df) == 0),
                ):
                    st.session_state[defs_key].append(
                        {"name": new_name, "expression": formula}
                    )
                    st.rerun()

        # ------------------------------------------------------------------
        # Applied formula list
        # ------------------------------------------------------------------
        if st.session_state[defs_key]:
            st.caption("Applied formulas")
            for idx, d in enumerate(list(st.session_state[defs_key])):
                fc1, fc2 = st.columns([5, 1])
                fc1.code(f"{d['name']} = {d['expression']}", language=None)
                if fc2.button("✕", key=f"{page_key}_del_{idx}"):
                    st.session_state[defs_key].pop(idx)
                    st.rerun()

            st.download_button(
                "⬇ Export formulas (JSON)",
                json.dumps(
                    {"version": 1, "formulas": st.session_state[defs_key]},
                    indent=2,
                ),
                file_name="formulas.json",
                key=f"{page_key}_dl",
            )

        # ------------------------------------------------------------------
        # JSON import
        # ------------------------------------------------------------------
        uploaded = st.file_uploader(
            "Import formulas (JSON)",
            type=["json"],
            key=f"{page_key}_ul",
        )
        if uploaded is not None:
            try:
                data = json.load(uploaded)
                st.session_state[defs_key] = data.get("formulas", [])
                st.rerun()
            except Exception:
                st.error("Invalid formula JSON.")

    # Always (re)apply saved formulas to the working df before returning
    return _apply_saved_formulas(df, st.session_state[defs_key])
