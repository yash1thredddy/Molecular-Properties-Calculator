# molecular_calculator/models/gmm.py
"""Gaussian Mixture Model analysis engine.

Streamlit-free, fully unit-testable. Ported and generalized from the
Impulator 3 ``imp_gmm.py`` (1-D IMP-score clustering) to support N-D
clustering over arbitrary numeric columns. Uses scikit-learn's
``GaussianMixture`` (covariance_type='full') with a fixed random seed for
reproducibility, BIC for model selection, and ``StandardScaler`` so that
columns on different scales contribute comparably.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)

MIN_COMPONENTS: int = 2
MAX_COMPONENTS: int = 6
DEFAULT_COMPONENTS: int = 3
DEFAULT_RANDOM_STATE: int = 42
N_INIT: int = 5


@dataclass
class PreparedData:
    """Cleaned, model-ready numeric data plus provenance for the UI."""
    values: np.ndarray
    kept_index: pd.Index
    n_dropped: int
    column_names: List[str]
    logged_columns: List[str] = field(default_factory=list)
    kept_positions: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))


def prepare_numeric_data(
    df: pd.DataFrame,
    columns: List[str],
    log_transform_cols: Optional[List[str]] = None,
) -> PreparedData:
    """Select ``columns``, coerce to numeric, drop NaN/inf rows, optionally log1p.

    Log-transform is applied per requested column ONLY when that column has no
    negative values (log of a negative is undefined). Skipped columns are
    reported in ``logged_columns`` so the UI can label axes honestly.
    """
    columns = list(columns)
    if not columns:
        raise ValueError("columns must not be empty")
    requested_log = list(log_transform_cols or [])

    sub = df[columns].apply(pd.to_numeric, errors="coerce")
    sub = sub.replace([np.inf, -np.inf], np.nan)
    mask = sub.notna().all(axis=1)
    kept_positions = np.flatnonzero(mask.to_numpy())
    kept = sub[mask]
    n_dropped = int((~mask).sum())
    if kept.empty:
        logger.warning(
            "prepare_numeric_data: all %d row(s) dropped for columns %s "
            "(check for an all-missing or non-numeric column)", len(df), columns,
        )
    values = kept.to_numpy(dtype=float).copy()

    applied_log: List[str] = []
    for i, col in enumerate(columns):
        if col in requested_log:
            col_vals = values[:, i]
            if col_vals.size and np.any(col_vals < 0):
                logger.warning("Skipping log-transform for %r: contains negatives", col)
                continue
            values[:, i] = np.log1p(col_vals)
            applied_log.append(col)

    return PreparedData(
        values=values,
        kept_index=kept.index,
        n_dropped=n_dropped,
        column_names=columns,
        logged_columns=applied_log,
        kept_positions=kept_positions,
    )


def gmm_sentinel_check(values, n_components: int) -> Tuple[bool, Optional[str]]:
    """Return (ok, message). ``ok`` is False when the data cannot support a fit.

    Generalizes Impulator's three sentinel conditions to N-D:
      * invalid request (n_components < MIN_COMPONENTS)
      * too few rows (n_samples <= n_components)
      * no variation (all rows identical)
      * too few distinct rows (unique rows < n_components)
    Messages are written in plain language for a non-technical user.

    Note: uniqueness is checked with exact float equality (``np.unique(axis=0)``),
    so this is a coarse guard — near-constant columns differing only by
    floating-point noise are treated as distinct.
    """
    if n_components < MIN_COMPONENTS:
        return False, (
            f"At least {MIN_COMPONENTS} groups are required; got {n_components}."
        )
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if arr.size and not np.isfinite(arr).all():
        return False, (
            "Some selected values are missing or non-numeric. Remove or fill those "
            "rows/columns before grouping."
        )
    n = arr.shape[0]
    # GMM requires n_samples > n_components (strictly), so n == n_components is rejected.
    if n <= n_components:
        return False, (
            f"Not enough data — fitting {n_components} groups needs at least "
            f"{n_components + 1} rows, but only {n} usable rows were found. "
            f"Try fewer groups."
        )
    unique_rows = np.unique(arr, axis=0)
    if unique_rows.shape[0] == 1:
        return False, (
            "No variation — every row has identical values for the selected "
            "columns. GMM needs variation to find groups."
        )
    if unique_rows.shape[0] < n_components:
        return False, (
            f"Too few distinct values — only {unique_rows.shape[0]} unique rows "
            f"but {n_components} groups requested. Try fewer groups."
        )
    return True, None


def _scale(values: np.ndarray, standardize: bool):
    """Return (X, scaler). When standardize is False, scaler is None."""
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if standardize:
        scaler = StandardScaler().fit(arr)
        return scaler.transform(arr), scaler
    return arr, None


def best_fit_k(values, *, k_min: int = MIN_COMPONENTS, k_max: int = MAX_COMPONENTS,
               random_state: int = DEFAULT_RANDOM_STATE, standardize: bool = True) -> int:
    """Return the K in [k_min, k_max] minimizing BIC. Falls back to DEFAULT_COMPONENTS clamped into [k_min, k_max].

    Ported from Impulator's ``best_fit_k`` and generalized to N-D. Skips
    non-converged fits. All candidates share the same seed and covariance type
    so their BIC values are directly comparable.
    """
    X, _ = _scale(values, standardize)
    best_k: Optional[int] = None
    best_bic = float("inf")
    for k in range(int(k_min), int(k_max) + 1):
        # GMM needs n_samples > k (strictly), matching gmm_sentinel_check. Skipping
        # k == n_samples avoids suggesting a K the UI would then reject.
        if X.shape[0] <= k:
            break
        try:
            cand = GaussianMixture(
                n_components=k, covariance_type="full",
                n_init=N_INIT, random_state=int(random_state),
            ).fit(X)
            if not cand.converged_:
                continue
            bic = float(cand.bic(X))
            if bic < best_bic:
                best_bic, best_k = bic, k
        except Exception as exc:  # singular covariance etc. — skip this K
            logger.debug("best_fit_k: K=%d skipped (%s)", k, exc)
            continue
    return best_k if best_k is not None else min(max(DEFAULT_COMPONENTS, int(k_min)), int(k_max))


def bic_aic_sweep(values, *, k_min: int = 1, k_max: int = 10,
                  random_state: int = DEFAULT_RANDOM_STATE, standardize: bool = True) -> pd.DataFrame:
    """Return a DataFrame with columns n_groups, bic, aic across the K range.

    Used to draw the model-quality plot. Wider default range (1-10) than the
    fit slider (2-6) so the user can see the full curve and its minimum.

    Non-converged fits are included here (BIC is still finite and useful for
    visualizing the curve shape); ``best_fit_k`` excludes them when selecting
    the winner.
    """
    X, _ = _scale(values, standardize)
    rows = []
    for k in range(int(k_min), int(k_max) + 1):
        if X.shape[0] < k:
            break
        try:
            m = GaussianMixture(
                n_components=k, covariance_type="full",
                n_init=N_INIT, random_state=int(random_state),
            ).fit(X)
            rows.append({"n_groups": k, "bic": float(m.bic(X)), "aic": float(m.aic(X))})
        except Exception as exc:
            logger.debug("bic_aic_sweep: K=%d skipped (%s)", k, exc)
            continue
    return pd.DataFrame(rows, columns=["n_groups", "bic", "aic"])


class GMMAnalysis:
    """Fit and interpret a Gaussian Mixture Model over one or more numeric columns.

    Group labels are remapped so that **group 0 is the lowest-mean group on the
    first selected column** (stable, human-friendly ordering; mirrors the
    Impulator ascending-mean contract). ``means_real_units`` and ``weights`` are
    returned in this same order.

    The fittable range is ``[MIN_COMPONENTS, MAX_COMPONENTS]`` (2-6), intentionally
    narrower than ``bic_aic_sweep``'s wider visualization range — the sweep curve is
    for inspection, the fit is the engine.
    """

    def __init__(self, values, column_names, *, n_components: int,
                 random_state: int = DEFAULT_RANDOM_STATE, standardize: bool = True):
        if not (MIN_COMPONENTS <= int(n_components) <= MAX_COMPONENTS):
            raise ValueError(
                f"n_components must be in [{MIN_COMPONENTS}, {MAX_COMPONENTS}]; "
                f"got {n_components}"
            )
        self.raw_values = np.asarray(values, dtype=float)
        if self.raw_values.ndim == 1:
            self.raw_values = self.raw_values.reshape(-1, 1)
        if self.raw_values.ndim != 2:
            raise ValueError("values must be 2-D: (n_samples, n_features)")

        self.column_names = list(column_names)
        self.n_components = int(n_components)
        self.n_features = self.raw_values.shape[1]
        if len(self.column_names) != self.n_features:
            raise ValueError(
                f"column_names has {len(self.column_names)} entries but values has "
                f"{self.n_features} feature column(s)."
            )

        X, self.scaler = _scale(self.raw_values, standardize)
        self._X = X
        self.model = GaussianMixture(
            n_components=self.n_components, covariance_type="full",
            n_init=N_INIT, random_state=int(random_state),
        ).fit(X)
        self.converged = bool(self.model.converged_)
        if not self.converged:
            logger.warning(
                "GMM did not converge (n_components=%s, n_samples=%s)",
                self.n_components, X.shape[0],
            )

        # Ascending-mean ordering on the FIRST feature -> stable group ids.
        self._order = np.argsort(self.model.means_[:, 0])
        remap = {int(old): int(new) for new, old in enumerate(self._order)}
        raw_labels = self.model.predict(X)
        self.labels = np.array([remap[int(lab)] for lab in raw_labels])

        proba = self.model.predict_proba(X)
        # proba columns are in the model's internal component order (NOT the
        # remapped ascending-mean order). .max() is order-independent, so this is
        # safe — but never use proba.argmax() for group ids; use self.labels.
        self.confidence = proba.max(axis=1)
        self.log_likelihood = self.model.score_samples(X)

    @property
    def weights(self) -> np.ndarray:
        return self.model.weights_[self._order]

    @property
    def means_real_units(self) -> np.ndarray:
        means = self.model.means_[self._order]
        if self.scaler is not None:
            means = self.scaler.inverse_transform(means)
        return means

    def component_curves(self, x_grid):
        """Per-component weighted Gaussian PDFs on a 1-D grid, ascending-mean order.

        Single-property only. Returns (means, weights, sigmas, pdfs) where
        ``pdfs[k] = weights[k] * norm.pdf(x_grid, means[k], sigmas[k])``. Ported
        directly from Impulator's ``component_curves``. ``x_grid`` is in the
        caller's real (un-standardized) units, so when this analysis was fit with
        ``standardize=True`` the model's z-space means/sigmas are converted back
        to real units before building the PDFs.
        """
        if self.n_features != 1:
            raise ValueError("component_curves is only defined for single-property (1-D) analyses")
        order = self._order
        means = self.model.means_.flatten()[order]
        weights = self.model.weights_[order]
        # covariances_ is (K, 1, 1) for covariance_type='full' with n_features==1;
        # [:, 0, 0] is self-documenting and fails loudly if a multi-feature array
        # is ever passed here by mistake.
        sigmas = np.sqrt(self.model.covariances_[:, 0, 0])[order]
        if self.scaler is not None:
            # Model was fit in standardized space; map means/sigmas back to the
            # real units of x_grid so the overlaid Gaussians line up with the data.
            means = self.scaler.inverse_transform(means.reshape(-1, 1)).flatten()
            sigmas = sigmas * float(self.scaler.scale_[0])
        x = np.asarray(x_grid, dtype=float)
        pdfs = np.vstack([
            weights[k] * norm.pdf(x, means[k], sigmas[k]) for k in range(len(means))
        ])
        return means, weights, sigmas, pdfs

    def outlier_mask(self, percentile: float = 5.0) -> np.ndarray:
        """Boolean mask of 'unusual' rows: bottom ``percentile`` by log-likelihood.

        Low ``score_samples`` means the row fits no group well — a candidate
        outlier worth a second look.
        """
        threshold = np.percentile(self.log_likelihood, percentile)
        mask = self.log_likelihood <= threshold
        if mask.all():
            # Degenerate case: all rows share (nearly) the same likelihood, so the
            # threshold flags everything. No row is meaningfully an outlier.
            logger.warning("outlier_mask: all log-likelihoods tied; flagging none")
            return np.zeros(len(self.log_likelihood), dtype=bool)
        return mask

    def labeled_dataframe(self, original_df: pd.DataFrame, kept_positions) -> pd.DataFrame:
        """Return a copy of ``original_df`` with GMM_Group (1-based) and
        GMM_Confidence_% columns. Dropped rows stay NA. ``kept_positions`` are the
        integer row positions that survived prepare_numeric_data — positional
        assignment is robust to non-unique/duplicate indexes."""
        out = original_df.copy()
        out["GMM_Group"] = pd.NA
        out["GMM_Confidence_%"] = pd.NA
        gcol = out.columns.get_loc("GMM_Group")
        ccol = out.columns.get_loc("GMM_Confidence_%")
        out.iloc[kept_positions, gcol] = (self.labels + 1).astype(int)
        out.iloc[kept_positions, ccol] = np.round(self.confidence * 100, 1)
        return out


__all__ = [
    "MIN_COMPONENTS", "MAX_COMPONENTS", "DEFAULT_COMPONENTS",
    "DEFAULT_RANDOM_STATE", "N_INIT",
    "PreparedData", "prepare_numeric_data",
    "best_fit_k", "bic_aic_sweep", "gmm_sentinel_check", "GMMAnalysis",
]
