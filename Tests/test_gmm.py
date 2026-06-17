# Tests/test_gmm.py
"""Unit tests for the GMM analysis engine (Streamlit-free)."""
import numpy as np
import pandas as pd
import pytest

from molecular_calculator.models.gmm import (
    MIN_COMPONENTS, MAX_COMPONENTS, DEFAULT_COMPONENTS, DEFAULT_RANDOM_STATE, N_INIT,
    gmm_sentinel_check,
)


def test_constants_match_impulator_reference():
    assert (MIN_COMPONENTS, MAX_COMPONENTS, DEFAULT_COMPONENTS) == (2, 6, 3)
    assert DEFAULT_RANDOM_STATE == 42
    assert N_INIT >= 5


def test_sentinel_too_few_samples():
    values = np.array([[1.0], [2.0]])  # 2 rows
    ok, msg = gmm_sentinel_check(values, n_components=3)
    assert ok is False
    assert "least 4 rows" in msg  # n_components + 1


def test_sentinel_zero_variation():
    values = np.array([[5.0], [5.0], [5.0], [5.0]])
    ok, msg = gmm_sentinel_check(values, n_components=2)
    assert ok is False
    assert "No variation" in msg


def test_sentinel_too_few_unique():
    values = np.array([[1.0], [1.0], [2.0], [2.0]])  # 2 unique rows
    ok, msg = gmm_sentinel_check(values, n_components=3)
    assert ok is False
    assert "distinct" in msg


def test_sentinel_passes_on_good_data():
    rng = np.random.RandomState(0)
    values = rng.normal(size=(50, 1))
    ok, msg = gmm_sentinel_check(values, n_components=3)
    assert ok is True
    assert msg is None


def test_sentinel_rejects_below_min_components():
    rng = np.random.RandomState(0)
    values = rng.normal(size=(50, 1))
    ok, msg = gmm_sentinel_check(values, n_components=1)
    assert ok is False
    assert "At least 2 groups" in msg


def test_sentinel_row_count_boundary():
    rng = np.random.RandomState(0)
    # exactly n_components rows -> rejected; n_components+1 rows -> passes row check
    at_boundary = rng.normal(size=(3, 1))
    ok_at, _ = gmm_sentinel_check(at_boundary, n_components=3)
    assert ok_at is False
    above = rng.normal(size=(4, 1))
    ok_above, msg_above = gmm_sentinel_check(above, n_components=3)
    assert ok_above is True
    assert msg_above is None


def test_sentinel_passes_on_multidimensional_data():
    rng = np.random.RandomState(1)
    values = rng.normal(size=(60, 3))  # genuine N-D (3 columns)
    ok, msg = gmm_sentinel_check(values, n_components=2)
    assert ok is True
    assert msg is None


from molecular_calculator.models.gmm import prepare_numeric_data, PreparedData


def test_prepare_drops_nan_and_inf_rows():
    df = pd.DataFrame({
        "A": [1.0, 2.0, np.nan, 4.0, np.inf],
        "B": [10.0, 20.0, 30.0, 40.0, 50.0],
        "name": ["a", "b", "c", "d", "e"],
    })
    prepared = prepare_numeric_data(df, ["A", "B"])
    assert isinstance(prepared, PreparedData)
    assert prepared.values.shape == (3, 2)        # rows 0,1,3 kept
    assert prepared.n_dropped == 2
    assert list(prepared.kept_index) == [0, 1, 3]  # rows 2 (NaN) and 4 (inf) dropped
    assert prepared.column_names == ["A", "B"]


def test_prepare_coerces_numeric_strings():
    df = pd.DataFrame({"A": ["1", "2", "3"], "B": [4.0, 5.0, 6.0]})
    prepared = prepare_numeric_data(df, ["A", "B"])
    assert prepared.values.shape == (3, 2)
    np.testing.assert_allclose(prepared.values[:, 0], [1.0, 2.0, 3.0])


def test_prepare_applies_log1p_only_to_requested_columns():
    df = pd.DataFrame({"A": [0.0, 9.0], "B": [0.0, 9.0]})
    prepared = prepare_numeric_data(df, ["A", "B"], log_transform_cols=["A"])
    np.testing.assert_allclose(prepared.values[:, 0], np.log1p([0.0, 9.0]))
    np.testing.assert_allclose(prepared.values[:, 1], [0.0, 9.0])
    assert prepared.logged_columns == ["A"]


def test_prepare_skips_log_on_negative_column():
    df = pd.DataFrame({"A": [-1.0, 2.0, 3.0]})
    prepared = prepare_numeric_data(df, ["A"], log_transform_cols=["A"])
    np.testing.assert_allclose(prepared.values[:, 0], [-1.0, 2.0, 3.0])  # unchanged
    assert prepared.logged_columns == []  # not applied


def test_prepare_raises_on_empty_columns():
    df = pd.DataFrame({"A": [1.0, 2.0]})
    with pytest.raises(ValueError):
        prepare_numeric_data(df, [])


def test_prepare_all_nan_column_drops_everything():
    df = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [np.nan, np.nan, np.nan]})
    prepared = prepare_numeric_data(df, ["A", "B"])
    assert prepared.values.shape == (0, 2)
    assert prepared.n_dropped == 3


from molecular_calculator.models.gmm import best_fit_k, bic_aic_sweep


def _two_cluster_data(seed=0, n=120):
    rng = np.random.RandomState(seed)
    a = rng.normal(loc=0.0, scale=0.5, size=(n // 2, 2))
    b = rng.normal(loc=8.0, scale=0.5, size=(n // 2, 2))
    return np.vstack([a, b])


def test_best_fit_k_recovers_two_clusters():
    values = _two_cluster_data()
    assert best_fit_k(values, k_min=2, k_max=6) == 2


def test_best_fit_k_is_reproducible():
    values = _two_cluster_data()
    assert best_fit_k(values) == best_fit_k(values)


def test_best_fit_k_never_exceeds_sample_count():
    # With only 3 samples, best_fit_k must stay below n_samples to match
    # gmm_sentinel_check (which requires n_samples > n_components strictly),
    # so "Auto (BIC)" never suggests a K the UI would then reject.
    values = np.array([[1.0], [2.0], [100.0]])  # 3 rows
    assert best_fit_k(values, k_min=2, k_max=6) < values.shape[0]


def test_facade_batch_augments_with_gmm_property_set():
    # The GMM page relies on the public MolecularCalculator facade (not the batch
    # page's private _process_batch) to add numeric property columns.
    import pandas as pd
    from molecular_calculator.core.molecular_calculator import MolecularCalculator
    df = pd.DataFrame({"SMILES": ["CCO", "c1ccccc1", "CC(=O)O"]})
    out = MolecularCalculator.process_batch(
        df, "SMILES", selected_properties={"Molecular_Weight", "LogP", "TPSA"})
    for prop in ("Molecular_Weight", "LogP", "TPSA"):
        assert prop in out.columns
    assert len(out) == 3


def test_bic_aic_sweep_shape_and_columns():
    values = _two_cluster_data()
    sweep = bic_aic_sweep(values, k_min=1, k_max=5)
    assert list(sweep.columns) == ["n_groups", "bic", "aic"]
    assert list(sweep["n_groups"]) == [1, 2, 3, 4, 5]
    # BIC for the correct K (2) should beat K=1
    bic_by_k = dict(zip(sweep["n_groups"], sweep["bic"]))
    assert bic_by_k[2] < bic_by_k[1]


from molecular_calculator.models.gmm import GMMAnalysis


def test_gmm_analysis_recovers_two_groups():
    values = _two_cluster_data()
    a = GMMAnalysis(values, ["X", "Y"], n_components=2)
    assert a.n_components == 2
    assert a.n_features == 2
    assert a.converged is True
    assert set(np.unique(a.labels)) == {0, 1}
    # confidence is a valid posterior max in [0, 1]
    assert a.confidence.min() >= 0.0 and a.confidence.max() <= 1.0
    assert a.confidence.shape == (values.shape[0],)


def test_gmm_group1_is_lowest_mean_on_first_column():
    values = _two_cluster_data()
    a = GMMAnalysis(values, ["X", "Y"], n_components=2)
    centers = a.means_real_units  # ascending-mean order on column 0
    assert centers[0, 0] < centers[1, 0]


def test_gmm_means_real_units_are_inverse_scaled():
    values = _two_cluster_data()
    a = GMMAnalysis(values, ["X", "Y"], n_components=2, standardize=True)
    centers = a.means_real_units
    # group centers should land near the true cluster means (0 and 8), not z-scores
    lows, highs = sorted(centers[:, 0])
    assert -1.5 < lows < 1.5
    assert 6.5 < highs < 9.5


def test_gmm_reproducible_labels():
    values = _two_cluster_data()
    a1 = GMMAnalysis(values, ["X", "Y"], n_components=2)
    a2 = GMMAnalysis(values, ["X", "Y"], n_components=2)
    np.testing.assert_array_equal(a1.labels, a2.labels)


def test_gmm_rejects_out_of_range_components():
    values = _two_cluster_data()
    with pytest.raises(ValueError):
        GMMAnalysis(values, ["X", "Y"], n_components=99)


def test_gmm_single_feature_smoke():
    rng = np.random.RandomState(3)
    values = np.concatenate([rng.normal(0, 0.5, 60), rng.normal(8, 0.5, 60)]).reshape(-1, 1)
    a = GMMAnalysis(values, ["score"], n_components=2, standardize=False)
    assert a.n_features == 1
    assert a.means_real_units.shape == (2, 1)
    assert set(np.unique(a.labels)) == {0, 1}


def test_gmm_means_real_units_without_standardize():
    values = _two_cluster_data()
    a = GMMAnalysis(values, ["X", "Y"], n_components=2, standardize=False)
    assert a.scaler is None  # means returned directly in original units
    centers = a.means_real_units
    lows, highs = sorted(centers[:, 0])
    assert -1.5 < lows < 1.5
    assert 6.5 < highs < 9.5


def test_gmm_rejects_column_name_count_mismatch():
    values = _two_cluster_data()  # 2 features
    with pytest.raises(ValueError):
        GMMAnalysis(values, ["only_one"], n_components=2)


def test_component_curves_1d_shapes():
    rng = np.random.RandomState(1)
    values = np.concatenate([rng.normal(0, 1, 60), rng.normal(10, 1, 60)]).reshape(-1, 1)
    a = GMMAnalysis(values, ["score"], n_components=2, standardize=False)
    grid = np.linspace(values.min(), values.max(), 200)
    means, weights, sigmas, pdfs = a.component_curves(grid)
    assert means.shape == (2,)
    assert pdfs.shape == (2, 200)
    assert means[0] < means[1]  # ascending


def test_component_curves_rejects_multifeature():
    values = _two_cluster_data()
    a = GMMAnalysis(values, ["X", "Y"], n_components=2)
    with pytest.raises(ValueError):
        a.component_curves(np.linspace(0, 1, 10))


def test_outlier_mask_flags_bottom_fraction():
    values = _two_cluster_data()
    a = GMMAnalysis(values, ["X", "Y"], n_components=2)
    mask = a.outlier_mask(percentile=5.0)
    assert mask.dtype == bool
    assert mask.shape == (values.shape[0],)
    assert 0 < mask.sum() <= max(1, int(0.10 * values.shape[0]))


def test_labeled_dataframe_aligns_and_marks_dropped_rows():
    df = pd.DataFrame({"A": [0.0, 0.1, 8.0, 8.1, np.nan], "B": [0.0, 0.2, 8.0, 8.2, 1.0]})
    prepared = prepare_numeric_data(df, ["A", "B"])
    a = GMMAnalysis(prepared.values, prepared.column_names, n_components=2)
    out = a.labeled_dataframe(df, prepared.kept_index)
    assert "GMM_Group" in out.columns and "GMM_Confidence_%" in out.columns
    assert len(out) == len(df)
    assert pd.isna(out.loc[4, "GMM_Group"])          # dropped row stays NA
    assert set(out.loc[[0, 1, 2, 3], "GMM_Group"]) == {1, 2}  # 1-based labels
    kept_conf = out.loc[[0, 1, 2, 3], "GMM_Confidence_%"]
    assert kept_conf.notna().all()
    assert ((kept_conf > 0) & (kept_conf <= 100)).all()


from molecular_calculator.ui.components.gmm_charts import (
    create_density_overlay, create_cluster_scatter, create_bic_aic_plot,
)


def test_density_overlay_has_histogram_and_component_traces():
    rng = np.random.RandomState(2)
    values = np.concatenate([rng.normal(0, 1, 60), rng.normal(10, 1, 60)]).reshape(-1, 1)
    a = GMMAnalysis(values, ["score"], n_components=2, standardize=False)
    fig = create_density_overlay(values.ravel(), a, "score")
    # 1 histogram + 2 component curves
    assert len(fig.data) == 3
    import plotly.graph_objects as go
    assert isinstance(fig.data[0], go.Histogram)
    assert all(isinstance(t, go.Scatter) for t in fig.data[1:])


def test_cluster_scatter_has_group_and_center_traces():
    values = _two_cluster_data()
    a = GMMAnalysis(values, ["X", "Y"], n_components=2)
    fig = create_cluster_scatter(a, x_idx=0, y_idx=1, x_name="X", y_name="Y")
    # 2 group traces + 1 centers trace
    assert len(fig.data) == 3


def test_bic_aic_plot_has_two_lines():
    values = _two_cluster_data()
    sweep = bic_aic_sweep(values, k_min=1, k_max=5)
    fig = create_bic_aic_plot(sweep, recommended_k=2)
    names = {t.name for t in fig.data}
    assert {"BIC", "AIC"}.issubset(names)
    assert len(fig.layout.shapes) == 1
    assert fig.layout.shapes[0].x0 == 2  # recommended_k


from molecular_calculator.models.gmm import gmm_sentinel_check


def test_sentinel_rejects_nan_inf():
    ok, msg = gmm_sentinel_check([[1.0], [2.0], [np.nan], [4.0], [5.0]], 2)
    assert ok is False
    assert "missing" in msg.lower() or "non-numeric" in msg.lower()


from molecular_calculator.models.gmm import prepare_numeric_data, GMMAnalysis


def test_labeled_dataframe_handles_duplicate_index():
    df = pd.DataFrame({"x": [1.0, 1.1, 5.0, 5.1, 9.0, 9.2]}, index=[0, 0, 1, 1, 2, 2])
    prep = prepare_numeric_data(df, ["x"])
    analysis = GMMAnalysis(prep.values, prep.column_names, n_components=2, random_state=42)
    out = analysis.labeled_dataframe(df, prep.kept_positions)
    assert out["GMM_Group"].notna().sum() == 6
