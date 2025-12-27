"""
3D OLS Regression Analysis Module

This module provides comprehensive 3D Ordinary Least Squares regression
with detailed statistical output similar to statsmodels.

Developed by: Yashwanth Reddy for ITR-UIC
Part of: Chemo-Informatics Toolkit
"""

import pandas as pd
import numpy as np
from datetime import datetime
from typing import Dict, Any, Optional, Tuple
from scipy import stats

# Import ThreeDOLSRegression from the legacy molecular_calculator.py file
# Using importlib to avoid naming conflict with the molecular_calculator package
import importlib.util
import os as _os

_module_path = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'molecular_calculator.py')
_spec = importlib.util.spec_from_file_location("molecular_calculator_legacy", _module_path)
if _spec is None or _spec.loader is None:
    raise ImportError("Could not load molecular_calculator.py")
_legacy_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_legacy_module)
ThreeDOLSRegression = _legacy_module.ThreeDOLSRegression


class RegressionSummary:
    """
    Generate comprehensive regression summary similar to statsmodels output
    """

    def __init__(self, model: ThreeDOLSRegression, x_name: str, y_name: str, z_name: str):
        """
        Initialize regression summary

        Args:
            model: Fitted ThreeDOLSRegression model
            x_name: Name of X variable
            y_name: Name of Y variable
            z_name: Name of Z variable (dependent)
        """
        self.model = model
        self.x_name = x_name
        self.y_name = y_name
        self.z_name = z_name

        # Calculate additional statistics
        self._calculate_statistics()

    def _calculate_statistics(self):
        """Calculate comprehensive regression statistics"""
        n = self.model.n
        k = 2  # Number of independent variables (X and Y)

        # Degrees of freedom
        self.df_model = k
        self.df_residuals = n - k - 1
        self.df_total = n - 1

        # Sum of squares
        z_mean = np.mean(self.model.z)
        z_predicted = self.model.predict(self.model.x, self.model.y)

        self.ss_total = np.sum((self.model.z - z_mean) ** 2)
        self.ss_residual = np.sum(self.model.residuals ** 2)
        self.ss_model = self.ss_total - self.ss_residual

        # Mean squared errors
        self.mse_model = self.ss_model / self.df_model if self.df_model > 0 else 0
        self.mse_residual = self.ss_residual / self.df_residuals if self.df_residuals > 0 else 0

        # R-squared and Adjusted R-squared
        self.r_squared = self.model.r_squared
        if self.df_total > 0:
            self.adj_r_squared = 1 - (1 - self.r_squared) * (self.df_total / self.df_residuals)
        else:
            self.adj_r_squared = self.r_squared

        # F-statistic
        if self.mse_residual > 0:
            self.f_statistic = self.mse_model / self.mse_residual
            self.f_pvalue = 1 - stats.f.cdf(self.f_statistic, self.df_model, self.df_residuals)
        else:
            self.f_statistic = np.inf
            self.f_pvalue = 0.0

        # Log-Likelihood (assuming normal errors) — use sigma^2_ML = RSS / n for statsmodels parity
        if self.ss_residual > 0:
            sigma2_mle = self.ss_residual / n
            self.log_likelihood = -0.5 * n * (np.log(2 * np.pi) + np.log(sigma2_mle) + 1)
        else:
            self.log_likelihood = np.inf

        # AIC and BIC (p = number of parameters = k predictors + intercept)
        p = k + 1
        self.aic = -2 * self.log_likelihood + 2 * p
        self.bic = -2 * self.log_likelihood + np.log(n) * p

        # Standard errors for coefficients
        self._calculate_coefficient_statistics()

        # Omnibus test for normality of residuals
        self._calculate_diagnostic_tests()

    def _calculate_coefficient_statistics(self):
        """Calculate standard errors, t-statistics, and p-values for coefficients"""
        n = self.model.n

        # Design matrix X with intercept
        X = np.column_stack([np.ones(n), self.model.x, self.model.y])

        # Variance-covariance matrix: σ² * (X'X)^-1
        XTX = X.T @ X
        try:
            XTX_inv = np.linalg.inv(XTX)
            var_covar = self.mse_residual * XTX_inv

            # Standard errors
            self.se_b0 = np.sqrt(var_covar[0, 0])
            self.se_b1 = np.sqrt(var_covar[1, 1])
            self.se_b2 = np.sqrt(var_covar[2, 2])

            # t-statistics
            self.t_b0 = self.model.b0 / self.se_b0 if self.se_b0 > 0 else 0
            self.t_b1 = self.model.b1 / self.se_b1 if self.se_b1 > 0 else 0
            self.t_b2 = self.model.b2 / self.se_b2 if self.se_b2 > 0 else 0

            # p-values (two-tailed)
            self.p_b0 = 2 * (1 - stats.t.cdf(abs(self.t_b0), self.df_residuals))
            self.p_b1 = 2 * (1 - stats.t.cdf(abs(self.t_b1), self.df_residuals))
            self.p_b2 = 2 * (1 - stats.t.cdf(abs(self.t_b2), self.df_residuals))

            # Confidence intervals (95%)
            t_critical = stats.t.ppf(0.975, self.df_residuals)

            self.ci_b0 = (self.model.b0 - t_critical * self.se_b0, self.model.b0 + t_critical * self.se_b0)
            self.ci_b1 = (self.model.b1 - t_critical * self.se_b1, self.model.b1 + t_critical * self.se_b1)
            self.ci_b2 = (self.model.b2 - t_critical * self.se_b2, self.model.b2 + t_critical * self.se_b2)

        except np.linalg.LinAlgError:
            # Singular matrix - set default values
            self.se_b0 = self.se_b1 = self.se_b2 = np.nan
            self.t_b0 = self.t_b1 = self.t_b2 = np.nan
            self.p_b0 = self.p_b1 = self.p_b2 = np.nan
            self.ci_b0 = self.ci_b1 = self.ci_b2 = (np.nan, np.nan)

    def _calculate_diagnostic_tests(self):
        """Calculate diagnostic tests for residuals"""
        residuals = self.model.residuals

        # Jarque-Bera and Omnibus (D’Agostino-Pearson) normality tests
        try:
            from scipy.stats import jarque_bera, normaltest
            self.jb_statistic, self.jb_pvalue = jarque_bera(residuals)
            self.omni_statistic, self.omni_pvalue = normaltest(residuals)
        except Exception:
            self.jb_statistic = np.nan
            self.jb_pvalue = np.nan
            self.omni_statistic = np.nan
            self.omni_pvalue = np.nan

        # Durbin-Watson statistic (test for autocorrelation)
        diff_residuals = np.diff(residuals)
        self.durbin_watson = np.sum(diff_residuals ** 2) / np.sum(residuals ** 2)

        # Condition Number
        n = self.model.n
        X = np.column_stack([np.ones(n), self.model.x, self.model.y])
        try:
            eigenvalues = np.linalg.eigvalsh(X.T @ X)
            self.condition_number = np.sqrt(eigenvalues.max() / eigenvalues.min())
        except:
            self.condition_number = np.nan

    def get_summary_dataframe(self) -> pd.DataFrame:
        """
        Get coefficient summary as a pandas DataFrame

        Returns:
            DataFrame with coefficient statistics
        """
        data = {
            'Variable': ['const', self.x_name, self.y_name],
            'Coefficient': [self.model.b0, self.model.b1, self.model.b2],
            'Std Error': [self.se_b0, self.se_b1, self.se_b2],
            't-statistic': [self.t_b0, self.t_b1, self.t_b2],
            'P>|t|': [self.p_b0, self.p_b1, self.p_b2],
            '[0.025': [self.ci_b0[0], self.ci_b1[0], self.ci_b2[0]],
            '0.975]': [self.ci_b0[1], self.ci_b1[1], self.ci_b2[1]]
        }

        return pd.DataFrame(data)

    def get_summary_text(self) -> str:
        """
        Generate statsmodels-style text summary

        Returns:
            Formatted text summary
        """
        summary = []
        summary.append("=" * 78)
        summary.append("                          OLS Regression Results")
        summary.append("=" * 78)

        # Header information
        summary.append(f"Dep. Variable:        {self.z_name:>20s}   R-squared:           {self.r_squared:>10.3f}")
        summary.append(f"Model:                             OLS   Adj. R-squared:      {self.adj_r_squared:>10.3f}")
        summary.append(f"Method:                 Least Squares   F-statistic:         {self.f_statistic:>10.3f}")
        summary.append(f"Date:                {datetime.now().strftime('%a, %d %b %Y'):>20s}   Prob (F-statistic):  {self.f_pvalue:>10.3e}")
        summary.append(f"Time:                    {datetime.now().strftime('%H:%M:%S'):>16s}   Log-Likelihood:      {self.log_likelihood:>10.3f}")
        summary.append(f"No. Observations:        {self.model.n:>16d}   AIC:                 {self.aic:>10.3f}")
        summary.append(f"Df Residuals:            {self.df_residuals:>16d}   BIC:                 {self.bic:>10.3f}")
        summary.append(f"Df Model:                {self.df_model:>16d}")
        summary.append(f"Covariance Type:            nonrobust")
        summary.append("=" * 78)

        # Coefficient table
        summary.append(f"{'Variable':>14s} {'Coefficient':>12s} {'Std Error':>12s} {'t':>8s} {'P>|t|':>8s} {'[0.025':>10s} {'0.975]':>10s}")
        summary.append("-" * 78)

        # Format coefficients
        coef_data = [
            ('const', self.model.b0, self.se_b0, self.t_b0, self.p_b0, self.ci_b0),
            (self.x_name[:14], self.model.b1, self.se_b1, self.t_b1, self.p_b1, self.ci_b1),
            (self.y_name[:14], self.model.b2, self.se_b2, self.t_b2, self.p_b2, self.ci_b2)
        ]

        for name, coef, se, t, p, ci in coef_data:
            summary.append(f"{name:>14s} {coef:>12.4f} {se:>12.3f} {t:>8.3f} {p:>8.3f} {ci[0]:>10.3f} {ci[1]:>10.3f}")

        summary.append("=" * 78)

        # Diagnostic statistics
        summary.append(f"Omnibus:              {getattr(self, 'omni_statistic', np.nan):>16.3f}   Durbin-Watson:       {self.durbin_watson:>10.3f}")
        summary.append(f"Prob(Omnibus):        {getattr(self, 'omni_pvalue', np.nan):>16.3f}   Jarque-Bera (JB):    {self.jb_statistic:>10.3f}")
        summary.append(f"Skew:                 {stats.skew(self.model.residuals):>16.3f}   Prob(JB):            {self.jb_pvalue:>10.3f}")
        summary.append(f"Kurtosis:             {stats.kurtosis(self.model.residuals, fisher=False):>16.3f}   Cond. No.            {self.condition_number:>10.1f}")
        summary.append("=" * 78)

        # Notes
        summary.append("\nNotes:")
        summary.append("[1] Standard Errors assume that the covariance matrix of the errors is correctly specified.")

        return "\n".join(summary)

    def get_statistics_dict(self) -> Dict[str, Any]:
        """
        Get all statistics as a dictionary

        Returns:
            Dictionary containing all regression statistics
        """
        return {
            # Model statistics
            'r_squared': self.r_squared,
            'adj_r_squared': self.adj_r_squared,
            'f_statistic': self.f_statistic,
            'f_pvalue': self.f_pvalue,
            'log_likelihood': self.log_likelihood,
            'aic': self.aic,
            'bic': self.bic,

            # Coefficients
            'b0': self.model.b0,
            'b1': self.model.b1,
            'b2': self.model.b2,

            # Standard errors
            'se_b0': self.se_b0,
            'se_b1': self.se_b1,
            'se_b2': self.se_b2,

            # t-statistics
            't_b0': self.t_b0,
            't_b1': self.t_b1,
            't_b2': self.t_b2,

            # p-values
            'p_b0': self.p_b0,
            'p_b1': self.p_b1,
            'p_b2': self.p_b2,

            # Sample statistics
            'n_observations': self.model.n,
            'df_model': self.df_model,
            'df_residuals': self.df_residuals,

            # Diagnostic tests
            'durbin_watson': self.durbin_watson,
            'jarque_bera': self.jb_statistic,
            'jarque_bera_pvalue': self.jb_pvalue,
            'condition_number': self.condition_number,

            # RMSE
            'rmse': self.model.rmse
        }


def perform_3d_regression(df: pd.DataFrame, x_col: str, y_col: str, z_col: str) -> Tuple[ThreeDOLSRegression, RegressionSummary]:
    """
    Perform 3D OLS regression and generate comprehensive summary

    Args:
        df: DataFrame containing the data
        x_col: Name of X variable column
        y_col: Name of Y variable column
        z_col: Name of Z variable (dependent) column

    Returns:
        Tuple of (model, summary)
    """
    # Extract data
    x = df[x_col].values
    y = df[y_col].values
    z = df[z_col].values

    # Fit model
    model = ThreeDOLSRegression(x, y, z)

    # Generate summary
    summary = RegressionSummary(model, x_col, y_col, z_col)

    return model, summary


def _build_design_matrix(x: np.ndarray, y: np.ndarray, include_interaction: bool, include_quadratic: bool) -> np.ndarray:
    """Build design matrix with optional interaction and quadratic terms."""
    terms = [np.ones_like(x), x, y]
    if include_interaction:
        terms.append(x * y)
    if include_quadratic:
        terms.append(x * x)
        terms.append(y * y)
    return np.column_stack(terms)


def _evaluate_3d_pair(
    df: pd.DataFrame,
    z_col: str,
    x_col: str,
    y_col: str,
    *,
    include_interaction: bool = False,
    include_quadratic: bool = False,
    cv_folds: int = 5,
) -> Dict[str, Any]:
    """
    Evaluate a 3D OLS model Z ~ 1 + X + Y and return key metrics.

    Uses the same conventions as RegressionSummary for llf/AIC/BIC.
    """
    # Drop rows with non-finite values in any of the 3 columns
    sub = df[[z_col, x_col, y_col]].copy()
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna()
    if len(sub) < 3:
        return {
            'x': x_col, 'y': y_col, 'n': len(sub), 'ok': False,
            'reason': 'insufficient_rows'
        }

    z = sub[z_col].to_numpy(dtype=float)
    x = sub[x_col].to_numpy(dtype=float)
    y = sub[y_col].to_numpy(dtype=float)

    n = len(sub)
    # Design matrix with intercept and optional terms
    X = _build_design_matrix(x, y, include_interaction, include_quadratic)
    p = X.shape[1]                      # number of parameters (incl. intercept)
    k = p - 1                           # number of predictors

    # Solve normal equations
    XT = X.T
    XTX = XT @ X
    try:
        XTX_inv = np.linalg.inv(XTX)
    except np.linalg.LinAlgError:
        return {
            'x': x_col, 'y': y_col, 'n': n, 'ok': False,
            'reason': 'singular_matrix'
        }
    beta = XTX_inv @ (XT @ z)
    b0, b1, b2 = beta.tolist()

    # Fitted and residuals
    zhat = X @ beta
    resid = z - zhat

    # Sums of squares
    zbar = np.mean(z)
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((z - zbar) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    df_res = n - p
    mse_res = ss_res / df_res if df_res > 0 else np.nan

    # Coefficient SEs, t, p
    se = np.sqrt(np.diag(mse_res * XTX_inv)) if np.isfinite(mse_res) else np.full(p, np.nan)
    with np.errstate(divide='ignore', invalid='ignore'):
        tvals = beta / se
    try:
        pvals = 2 * (1 - stats.t.cdf(np.abs(tvals), df_res))
    except Exception:
        pvals = np.array([np.nan, np.nan, np.nan])

    # F-statistic
    ss_model = ss_tot - ss_res
    mse_model = ss_model / k if k > 0 else np.nan
    f_stat = (mse_model / mse_res) if (mse_res and mse_res > 0) else np.nan
    f_pvalue = 1 - stats.f.cdf(f_stat, k, df_res) if np.isfinite(f_stat) else np.nan

    # Log-likelihood (use MLE sigma^2 = RSS / n)
    sigma2_mle = ss_res / n
    llf = -0.5 * n * (np.log(2 * np.pi) + np.log(sigma2_mle) + 1) if sigma2_mle > 0 else -np.inf
    aic = -2 * llf + 2 * p
    bic = -2 * llf + np.log(n) * p

    # Simple collinearity signal: VIF for X and Y using r_xy
    r_xy = np.corrcoef(x, y)[0, 1]
    r_xy2 = r_xy * r_xy
    with np.errstate(divide='ignore', invalid='ignore'):
        vif_x = float(1.0 / (1.0 - r_xy2)) if r_xy2 < 0.999999 else np.inf
        vif_y = vif_x

    # K-fold CV RMSE (deterministic split by index to avoid RNG dependency)
    rmse_cv = np.nan
    if cv_folds and cv_folds >= 3 and n > cv_folds:
        fold_sizes = np.full(cv_folds, n // cv_folds, dtype=int)
        fold_sizes[: n % cv_folds] += 1
        indices = np.arange(n)
        current = 0
        se_sum = 0.0
        m_sum = 0
        for fold_size in fold_sizes:
            start, stop = current, current + fold_size
            test_idx = indices[start:stop]
            train_idx = np.concatenate([indices[:start], indices[stop:]])
            current = stop

            X_tr, z_tr = X[train_idx], z[train_idx]
            X_te, z_te = X[test_idx], z[test_idx]
            XT_tr = X_tr.T
            XTX_tr = XT_tr @ X_tr
            try:
                beta_tr = np.linalg.inv(XTX_tr) @ (XT_tr @ z_tr)
            except np.linalg.LinAlgError:
                beta_tr = None
            if beta_tr is None:
                continue
            z_pred = X_te @ beta_tr
            err = z_te - z_pred
            se_sum += float(np.sum(err ** 2))
            m_sum += len(z_te)
        if m_sum > 0:
            rmse_cv = float(np.sqrt(se_sum / m_sum))

    model_spec = 'linear'
    if include_interaction and include_quadratic:
        model_spec = 'linear + xy + x2+y2'
    elif include_interaction:
        model_spec = 'linear + xy'
    elif include_quadratic:
        model_spec = 'linear + x2+y2'

    return {
        'x': x_col,
        'y': y_col,
        'n': n,
        'ok': True,
        # Return first three coefficients if present for compatibility
        'b0': float(beta[0]) if p >= 1 else np.nan,
        'b1': float(beta[1]) if p >= 2 else np.nan,
        'b2': float(beta[2]) if p >= 3 else np.nan,
        'se_b0': float(se[0]) if p >= 1 else np.nan,
        'se_b1': float(se[1]) if p >= 2 else np.nan,
        'se_b2': float(se[2]) if p >= 3 else np.nan,
        't_b0': float(tvals[0]) if p >= 1 else np.nan,
        't_b1': float(tvals[1]) if p >= 2 else np.nan,
        't_b2': float(tvals[2]) if p >= 3 else np.nan,
        'p_b0': float(pvals[0]) if p >= 1 else np.nan,
        'p_b1': float(pvals[1]) if p >= 2 else np.nan,
        'p_b2': float(pvals[2]) if p >= 3 else np.nan,
        'r2': float(r2),
        'adj_r2': float(1 - (1 - r2) * ((n - 1) / df_res)) if df_res > 0 else float(r2),
        'f_stat': float(f_stat) if np.isfinite(f_stat) else np.nan,
        'f_pvalue': float(f_pvalue) if np.isfinite(f_pvalue) else np.nan,
        'llf': float(llf), 'aic': float(aic), 'bic': float(bic),
        'rmse_cv': rmse_cv,
        'vif_x': float(vif_x), 'vif_y': float(vif_y),
        'r_xy': float(r_xy),
        'model_spec': model_spec,
        'p_params': int(p)
    }


def suggest_best_3d_pairs(
    df: pd.DataFrame,
    z_col: str,
    candidate_cols: Optional[list] = None,
    top_n: int = 5,
    *,
    include_interaction: bool = False,
    include_quadratic: bool = False,
    cv_folds: int = 5,
) -> Tuple[pd.DataFrame, Optional[Dict[str, Any]]]:
    """
    Suggest the best (X, Y) pairs for a given dependent variable Z.

    - Considers numeric columns (or `candidate_cols` if provided), excluding `z_col`.
    - Scores pairs primarily by highest adjusted R^2, then lowest BIC.
    - Returns a DataFrame of top pairs (proof) and the top suggestion dict.
    """
    # Identify numeric candidate columns
    if candidate_cols is None:
        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        candidates = [c for c in numeric_cols if c != z_col]
    else:
        candidates = [c for c in candidate_cols if c != z_col]

    # Evaluate all pairs
    results = []
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            x_col, y_col = candidates[i], candidates[j]
            res = _evaluate_3d_pair(
                df, z_col, x_col, y_col,
                include_interaction=include_interaction,
                include_quadratic=include_quadratic,
                cv_folds=cv_folds,
            )
            if res.get('ok'):
                results.append(res)

    if not results:
        return pd.DataFrame(), None

    res_df = pd.DataFrame(results)

    # Rank: highest adj_r2, then lowest bic, then lowest rmse_cv
    res_df = res_df.sort_values(
        by=['adj_r2', 'bic', 'rmse_cv'], ascending=[False, True, True]
    ).reset_index(drop=True)
    top = res_df.iloc[0].to_dict() if len(res_df) > 0 else None

    # Keep only essential proof columns for display
    proof_cols = [
        'x', 'y', 'n', 'model_spec', 'p_params', 'adj_r2', 'r2', 'rmse_cv', 'aic', 'bic',
        'f_stat', 'f_pvalue', 'p_b1', 'p_b2', 'vif_x', 'vif_y', 'r_xy'
    ]
    proof_df = res_df[proof_cols].copy()
    return proof_df.head(top_n), top
