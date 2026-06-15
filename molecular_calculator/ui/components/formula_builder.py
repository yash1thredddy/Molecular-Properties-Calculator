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
import html
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


def _formula_input(default_expr: str, cols: list[str], page_key: str, gen: int = 0) -> str:
    """Render the formula input and return its current text.

    Prefers the custom Monaco ``formula_editor`` component (column + function
    autocomplete); falls back to ``st.text_area`` if the component is not
    built/importable, so the page never breaks.

    The widget key is suffixed with a hash of ``default_expr`` so that applying
    a preset (which changes ``default_expr``) remounts the input seeded with the
    preset text — Monaco's ``defaultValue`` does not re-sync after mount, and a
    keyed ``text_area`` is likewise sticky, so a stable key would ignore the new
    seed. While the user types (no preset change) the seed is constant, so the
    key is stable and edits persist across reruns. ``gen`` is a form-generation
    counter: bumping it (after Apply) changes the key so the editor remounts
    empty, clearing the formula for the next column.
    """
    seed = hashlib.md5(default_expr.encode()).hexdigest()[:8]
    st.markdown("**Formula**")
    try:
        from molecular_calculator.components.formula_editor import formula_editor

        st.caption(
            "Start typing a column name or function — autocomplete suggests "
            "`[Column]` references and `func()` calls. Pick from the list rather "
            "than typing the `[` yourself (the suggestion inserts the full "
            "`[Column]`)."
        )
        return formula_editor(
            value=default_expr,
            columns=cols,
            functions=fe.editor_function_names(),
            height=120,
            key=f"{page_key}_editor_{gen}_{seed}",
        )
    except Exception:
        # Component not built or failed to load — graceful text_area fallback.
        return st.text_area(
            "Formula",
            value=default_expr,
            key=f"{page_key}_fexpr_{gen}_{seed}",
            height=80,
            label_visibility="collapsed",
        )


def _chip(label: str, highlight: bool) -> str:
    """Return an HTML chip span for a column name (label is HTML-escaped)."""
    safe = html.escape(str(label))
    style = (
        "background:#FF6B6B;color:#fff;border:1px solid #FF6B6B;"
        if highlight
        else "background:#F0F2F6;color:#262730;border:1px solid #E0E0E0;"
    )
    return (
        f"<span style='{style}padding:1px 8px;border-radius:10px;"
        f"font-size:0.78em;font-family:monospace'>{safe}</span>"
    )


def _render_column_chips(cols: list[str], working_name: str) -> None:
    """Show every column available to a formula as a chip, highlighting the
    column currently being created (``working_name``).

    Replaces the old single-column reference picker: the chips give an
    at-a-glance list of valid ``[Column]`` names while the editor's autocomplete
    handles actual insertion. The new column isn't in ``cols`` yet, so it's shown
    as a pending highlighted chip (or, if the name matches an existing column,
    that existing column is highlighted to flag the overwrite).
    """
    if not cols and not working_name:
        return
    chips = [_chip(c, bool(working_name) and c == working_name) for c in cols]
    if working_name and working_name not in cols:
        chips.append(_chip(f"{working_name} (new)", True))
    st.caption("Columns available to reference:")
    st.markdown(
        "<div style='display:flex;flex-wrap:wrap;gap:4px;margin-bottom:6px'>"
        + "".join(chips)
        + "</div>",
        unsafe_allow_html=True,
    )


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

    # Form-generation counter: bumped after each Apply so the name + formula
    # inputs remount empty for the next column. Without this, a just-applied
    # column lingers in the form and spuriously trips the "name already exists"
    # check (the column now exists because Apply added it).
    gen_key = f"{page_key}_formgen"
    gen = st.session_state.get(gen_key, 0)

    # Apply already-saved formulas up front so their columns are part of `cols`
    # below. That lets the editor's autocomplete, the live preview, and any
    # further formulas reference previously-created columns (formula chaining).
    df = _apply_saved_formulas(df, st.session_state[defs_key])

    with st.expander("➕ Custom calculations", expanded=False):
        cols = list(df.columns)

        # Preset gallery. The Monaco editor's autocomplete replaces the old
        # column-reference picker and static function list.
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
            key=f"{page_key}_fname_{gen}",
        )

        # Reference: every column available to the formula, with the column being
        # created highlighted. Replaces the old single-column reference picker.
        _render_column_chips(cols, new_name)

        formula = _formula_input(default_expr, cols, page_key, gen)

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
                st.caption(
                    "Preview only — this is a draft of the first rows. The column "
                    "is added to your data when you click **Apply column** below."
                )
                prev = preview_df.copy()
                prev[new_name or "result"] = series.values
                st.dataframe(
                    prev[[*cols[:2], new_name or "result"]].head(10),
                    width='stretch',
                )
                if failed:
                    st.caption(
                        f"{failed} of {len(preview_df)} preview rows blank "
                        "(missing values)."
                    )

                if new_name in cols:
                    st.caption(
                        f"ℹ️ A column named **{new_name}** already exists — "
                        "applying will overwrite it."
                    )

                # Task 17 — large-dataset performance warning
                if len(df) > 25000:
                    st.caption(
                        f"Note: applying to {len(df):,} rows is single-threaded "
                        "and may take a few seconds."
                    )

                # Highlight the button (primary colour) only once it's ready to
                # apply — a valid formula (we're in the no-error branch) plus a
                # column name and non-empty data.
                ready = bool(new_name) and len(df) > 0
                if st.button(
                    "Apply column",
                    key=f"{page_key}_apply",
                    type="primary" if ready else "secondary",
                    disabled=not ready,
                ):
                    overwrote = new_name in cols
                    new_def = {"name": new_name, "expression": formula}
                    # Replace an existing definition of the same name (so the
                    # applied-formulas list doesn't accumulate duplicates), else
                    # append a new one.
                    defs = st.session_state[defs_key]
                    for i, d in enumerate(defs):
                        if d["name"] == new_name:
                            defs[i] = new_def
                            break
                    else:
                        defs.append(new_def)
                    # Reset the form (name + formula) so the next column starts
                    # fresh and the just-applied name doesn't trip the warning.
                    st.session_state[gen_key] = gen + 1
                    st.session_state.pop(f"{page_key}_preset", None)
                    st.toast(
                        f"{'Updated' if overwrote else 'Added'} column "
                        f"'{new_name}'",
                        icon="✅",
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

    # `df` already has the saved formulas applied (done up front for chaining).
    return df
