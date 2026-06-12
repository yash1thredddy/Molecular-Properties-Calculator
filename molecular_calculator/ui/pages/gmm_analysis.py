# molecular_calculator/ui/pages/gmm_analysis.py
"""GMM Analysis page.

Fits a Gaussian Mixture Model to user-selected numeric columns of an uploaded
CSV/XLSX and presents the resulting groups in plain language. Two modes:
single-property (density overlay) and multi-property (clustering scatter).
"""

import numpy as np
import pandas as pd
import streamlit as st

from molecular_calculator.models.gmm import (
    DEFAULT_RANDOM_STATE, MIN_COMPONENTS, MAX_COMPONENTS,
    prepare_numeric_data, best_fit_k, bic_aic_sweep, gmm_sentinel_check, GMMAnalysis,
)
from molecular_calculator.core.molecular_calculator import MolecularCalculator
from molecular_calculator.ui.components.gmm_charts import (
    create_density_overlay, create_cluster_scatter, create_bic_aic_plot,
)
from molecular_calculator.ui.components.file_uploader import (
    render_file_uploader, validate_uploaded_file, read_uploaded_file,
    detect_smiles_column, create_download_button,
)
from molecular_calculator.utils.validators import DataFrameValidator
from molecular_calculator.utils.session_state import SessionState


def render_gmm_page() -> None:
    """Render the GMM Analysis page."""
    st.header("GMM Analysis")
    st.write(
        "Upload a CSV or Excel file and discover the natural **groups** hidden in "
        "your data using a Gaussian Mixture Model. Pick one property to see how it "
        "splits, or several properties to find clusters across them."
    )

    df = _load_dataframe()
    if df is None:
        st.info("👆 Upload a file to get started")
        return

    st.success(f"✅ Loaded {len(df):,} rows, {len(df.columns)} columns")
    with st.expander("📋 Data Preview (first 50 rows)", expanded=True):
        st.dataframe(df.head(50), width="stretch")

    # Optional: compute molecular properties first when a SMILES column exists.
    df = _maybe_calculate_properties(df)
    SessionState.set("gmm_working_df", df)

    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if len(numeric_cols) < 1:
        st.warning("No numeric columns found. Upload data with numeric values, or "
                   "use 'Calculate properties first' above.")
        return

    st.markdown("---")
    st.subheader("1. Choose analysis mode")
    mode = st.radio(
        "How many properties do you want to analyze?",
        ["Single property", "Multiple properties"],
        horizontal=True, key="gmm_mode",
        help="Single property shows a distribution split into groups. "
             "Multiple properties finds clusters across several columns.",
    )

    if mode == "Multiple properties" and len(numeric_cols) < 2:
        st.warning("Multi-property mode needs at least 2 numeric columns. "
                   "Switch to Single property, or add more numeric columns.")
        return

    st.subheader("2. Select columns")
    selected_cols, x_name, y_name = _select_columns(df, numeric_cols, mode)
    if not selected_cols:
        return

    # Controls + run are added in the next task.
    SessionState.set("gmm_selected_cols", selected_cols)
    SessionState.set("gmm_axes", (x_name, y_name))

    st.markdown("---")
    st.subheader("3. Settings")
    is_single = mode == "Single property"

    if "gmm_seed" not in st.session_state:
        st.session_state["gmm_seed"] = DEFAULT_RANDOM_STATE

    # Drop any previously-chosen log-transform columns that are no longer in the
    # current selection. Otherwise the keyed multiselect below (and this
    # prepare_numeric_data call) would receive stale values absent from its
    # options, which Streamlit rejects with an exception.
    if "gmm_log_cols" in st.session_state:
        st.session_state["gmm_log_cols"] = [
            c for c in st.session_state["gmm_log_cols"] if c in selected_cols
        ]

    prepared = prepare_numeric_data(
        df, selected_cols,
        log_transform_cols=st.session_state.get("gmm_log_cols", []),
    )
    if prepared.n_dropped:
        st.caption(f"ℹ️ {prepared.n_dropped} row(s) with missing/invalid values were skipped.")

    ca, cb = st.columns(2)
    with ca:
        standardize = st.checkbox(
            "Standardize columns (recommended)", value=not is_single,
            key="gmm_standardize",
            help="Puts every property on the same scale so none dominates. "
                 "Off for single-property so the chart stays in real units.",
        )
        st.multiselect(
            "Log-transform these columns (for skewed data)",
            selected_cols, default=[], key="gmm_log_cols",
        )
    with cb:
        def _new_seed():
            st.session_state["gmm_seed"] = int(np.random.randint(0, 2**31))

        def _reset_seed():
            st.session_state["gmm_seed"] = DEFAULT_RANDOM_STATE

        st.number_input(
            "Random seed", min_value=0, max_value=2_147_483_647, step=1,
            key="gmm_seed", help="Fixing the seed makes results reproducible.",
        )
        sc1, sc2 = st.columns(2)
        sc1.button("🎲 Refit (new seed)", key="gmm_new_seed", on_click=_new_seed)
        sc2.button("↺ Reset seed", key="gmm_reset_seed", on_click=_reset_seed)

    seed = int(st.session_state["gmm_seed"])
    standardize = bool(standardize)

    ok, sentinel = gmm_sentinel_check(prepared.values, MIN_COMPONENTS)
    if not ok:
        st.warning(sentinel)
        return

    suggested = _cached_best_fit_k(prepared.values, seed, standardize)
    st.caption(f"💡 BIC suggests **{suggested}** groups.")

    if "gmm_n_slider" not in st.session_state:
        st.session_state["gmm_n_slider"] = suggested

    def _auto_k():
        # Captures `prepared.values` and `standardize` from THIS render. Safe because
        # Streamlit serializes interactions — only one widget fires per rerun, so the
        # captured values always reflect the last settled UI state when this fires.
        st.session_state["gmm_n_slider"] = _cached_best_fit_k(
            prepared.values, int(st.session_state["gmm_seed"]), standardize,
        )

    cc, cd = st.columns([3, 1])
    with cc:
        n_components = st.slider(
            "Number of groups", min_value=MIN_COMPONENTS, max_value=MAX_COMPONENTS,
            key="gmm_n_slider",
        )
    with cd:
        st.write("")
        st.write("")
        st.button("Auto (BIC)", key="gmm_auto_k", on_click=_auto_k)

    run = st.button("▶ Run GMM analysis", type="primary", key="gmm_run")
    if run:
        _run_and_store(prepared, df, n_components, seed, standardize)

    _render_results_if_available(df)


def _load_dataframe():
    """Upload, validate, and read a file into a DataFrame (or None)."""
    uploaded = render_file_uploader(key="gmm_uploader")
    if uploaded is None:
        return None
    ok, error, warnings = validate_uploaded_file(uploaded)
    if not ok:
        st.error(f"❌ {error}")
        return None
    for w in warnings:
        st.warning(w)
    if SessionState.file_changed(uploaded):
        for k in ("gmm_analysis", "gmm_kept_index", "gmm_prepared", "gmm_props_df",
                  "gmm_run_params", "gmm_n_slider", "gmm_log_cols"):
            SessionState.clear(k)
    df = read_uploaded_file(uploaded)
    if df is None:
        return None
    result = DataFrameValidator.validate(df)
    if not result.is_valid:
        st.error(f"❌ {result.errors[0]}")
        return None
    for w in result.warnings:
        st.warning(w)
    return df


def _maybe_calculate_properties(df):
    """If a SMILES column is detected, offer to compute molecular properties and
    append them as new numeric columns. Reuses the batch pipeline. Returns the
    (possibly augmented) DataFrame."""
    smiles_col = detect_smiles_column(df)
    if smiles_col is None:
        return df

    with st.container(border=True):
        st.markdown("#### 🧪 Calculate molecular properties from SMILES")
        st.caption(
            f"Detected a SMILES column (**{smiles_col}**). Optionally calculate "
            "molecular properties (MW, LogP, TPSA, QED, and more) and add them to "
            "the columns available for GMM — handy when your file has structures "
            "but not numeric data yet."
        )
        if st.button("Calculate properties", key="gmm_calc_props"):
            props = ["Molecular_Weight", "LogP", "TPSA", "HB_Donors", "HB_Acceptors",
                     "Rotatable_Bonds", "Aromatic_Rings", "QED"]
            # Reuse the core facade's public, UI-free batch API rather than the batch
            # page's private helper. Parallel for larger files, sequential below.
            batch_fn = (
                MolecularCalculator.process_batch_parallel
                if len(df) > 50 else MolecularCalculator.process_batch
            )
            with st.spinner("Calculating properties..."):
                results = batch_fn(
                    df, smiles_col, selected_properties=set(props),
                    enable_online_lookup=True,
                )
            if results is not None and not results.empty:
                SessionState.set("gmm_props_df", results)
                st.success("Properties calculated and added to the column list.")
    augmented = SessionState.get("gmm_props_df")
    return augmented if augmented is not None else df


def _select_columns(df, numeric_cols, mode):
    """Return (selected_cols, x_name, y_name). In single-property mode y_name is None
    and x_name is the chosen column; in multi mode x_name/y_name are the scatter axes."""
    if mode == "Single property":
        col = st.selectbox("Property to analyze", numeric_cols, key="gmm_single_col")
        return [col], col, None

    selected = st.multiselect(
        "Properties to cluster on (pick 2 or more)",
        numeric_cols, default=numeric_cols[:2], key="gmm_multi_cols",
    )
    if len(selected) < 2:
        st.warning("Select at least 2 properties for multi-property clustering.")
        return [], None, None

    c1, c2 = st.columns(2)
    with c1:
        x_name = st.selectbox("Scatter X axis", selected, index=0, key="gmm_x_axis")
    with c2:
        y_choices = [c for c in selected if c != x_name]
        y_name = st.selectbox("Scatter Y axis", y_choices, index=0, key="gmm_y_axis")
    return selected, x_name, y_name


@st.cache_data(show_spinner=False)
def _cached_best_fit_k(values, random_state, standardize):
    """Cached BIC-based K suggestion so widget reruns don't recompute the fit."""
    return best_fit_k(values, random_state=random_state, standardize=standardize)


@st.cache_data(show_spinner=False)
def _cached_gmm_fit(values, column_names, n_components, random_state, standardize):
    """Cached GMM fit keyed on data + params; avoids refitting on unrelated reruns."""
    return GMMAnalysis(values, list(column_names), n_components=n_components,
                       random_state=random_state, standardize=standardize)


def _run_and_store(prepared, df, n_components, seed, standardize):
    ok, sentinel = gmm_sentinel_check(prepared.values, n_components)
    if not ok:
        st.warning(sentinel)
        return
    analysis = _cached_gmm_fit(prepared.values, tuple(prepared.column_names),
                               n_components, seed, standardize)
    SessionState.set("gmm_analysis", analysis)
    SessionState.set("gmm_kept_index", prepared.kept_index)
    SessionState.set("gmm_prepared", prepared)
    SessionState.set("gmm_run_params", {"seed": seed, "standardize": standardize})
    if not analysis.converged:
        st.info("The model did not fully converge. Try fewer groups or different columns.")


def _render_results_if_available(df):
    analysis = SessionState.get("gmm_analysis")
    prepared = SessionState.get("gmm_prepared")
    kept_index = SessionState.get("gmm_kept_index")
    if analysis is None or prepared is None:
        return

    selected_cols = analysis.column_names
    x_name, y_name = st.session_state.get("gmm_axes", (selected_cols[0], None))

    st.markdown("---")
    st.subheader("Results")

    # 1. Plain-language headline.
    st.markdown(
        f"### Your data splits into **{analysis.n_components} natural groups**."
    )
    with st.expander("What does this mean?"):
        st.markdown(
            "- A **group** is a cluster of rows whose values look similar.\n"
            "- **Group size %** is the share of your data that falls in each group.\n"
            "- **Confidence** is how sure the model is about each row's group "
            "(100% = a clear fit; near 50% = sits between groups).\n"
            "- We used a Gaussian Mixture Model (scikit-learn) with a fixed random "
            "seed so the result is reproducible. The number of groups was chosen "
            "with the BIC score (lower = better fit without over-complicating)."
        )

    # 2. Per-group summary table — averages in REAL units from the original rows.
    summary = _summary_table(df, kept_index, analysis, selected_cols)
    st.markdown("**Group summary**")
    st.dataframe(summary, width="stretch", hide_index=True)

    # 3. Headline chart.
    if analysis.n_features == 1:
        # Label the density axis log(...) when the column was log-transformed, so the
        # chart's x-axis/legend stays honest and consistent with the multi-property
        # scatter (the real-unit averages live in the summary table above).
        density_label = (
            f"log({selected_cols[0]})"
            if selected_cols[0] in prepared.logged_columns
            else selected_cols[0]
        )
        fig = create_density_overlay(prepared.values[:, 0], analysis, density_label)
        st.plotly_chart(fig, width="stretch")
    elif x_name in selected_cols and y_name in selected_cols:
        xi = selected_cols.index(x_name)
        yi = selected_cols.index(y_name)
        xlabel = f"log({x_name})" if x_name in prepared.logged_columns else x_name
        ylabel = f"log({y_name})" if y_name in prepared.logged_columns else y_name
        fig = create_cluster_scatter(analysis, xi, yi, xlabel, ylabel)
        st.plotly_chart(fig, width="stretch")
    else:
        st.info("Column selection changed since the last run — click ▶ Run GMM "
                "analysis to refresh the chart.")

    # 4. Model-quality plot. Use the params the model was ACTUALLY fit with
    # (stored at run time) so the curve and its recommended-K marker stay
    # consistent even if the user nudges the live widgets afterward.
    run_params = SessionState.get("gmm_run_params", {})
    seed = int(run_params.get("seed", DEFAULT_RANDOM_STATE))
    standardize = bool(run_params.get("standardize", analysis.scaler is not None))
    # Sweep a WIDER range than the fit slider (MIN..MAX) on purpose: the
    # model-quality curve is for inspection, so the user can see BIC/AIC behavior
    # beyond the selectable range and judge where the real minimum sits.
    sweep = _cached_bic_aic_sweep(prepared.values, 1, 10, seed, standardize)
    if not sweep.empty:
        st.markdown("**Model quality (how many groups?)**")
        st.plotly_chart(
            create_bic_aic_plot(sweep, recommended_k=analysis.n_components),
            width="stretch",
        )

    # 5. Confidence & outliers.
    outliers = analysis.outlier_mask(percentile=5.0)
    n_low_conf = int((analysis.confidence < 0.6).sum())
    st.markdown(
        f"**Confidence:** {n_low_conf} row(s) sit between groups (confidence < 60%). "
        f"**Unusual rows:** {int(outliers.sum())} row(s) fit no group well and may be "
        "worth a closer look."
    )

    # 6. Download labeled CSV.
    labeled = analysis.labeled_dataframe(df, kept_index)
    create_download_button(
        labeled, filename="gmm_results.csv",
        label="⬇ Download data with group labels", key="gmm_download",
    )


@st.cache_data(show_spinner=False)
def _cached_bic_aic_sweep(values, k_min, k_max, random_state, standardize):
    """Cached BIC/AIC sweep so the model-quality plot doesn't recompute on reruns."""
    return bic_aic_sweep(values, k_min=k_min, k_max=k_max,
                         random_state=random_state, standardize=standardize)


def _summary_table(df, kept_index, analysis, selected_cols):
    """Per-group table: size, % and the real-unit average of each selected column,
    computed from the ORIGINAL (pre-transform) kept rows so values are honest."""
    kept = df.loc[kept_index, selected_cols].apply(pd.to_numeric, errors="coerce")
    kept = kept.copy()
    kept["__group__"] = analysis.labels + 1
    n = len(kept)
    rows = []
    for g in range(1, analysis.n_components + 1):
        members = kept[kept["__group__"] == g]
        row = {"Group": g, "Count": len(members),
               "% of data": round(100.0 * len(members) / n, 1)}
        for col in selected_cols:
            row[f"Avg {col}"] = round(float(members[col].mean()), 3)
        rows.append(row)
    return pd.DataFrame(rows)
