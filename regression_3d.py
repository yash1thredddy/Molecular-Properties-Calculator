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
from molecular_calculator import ThreeDOLSRegression


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

        # Log-Likelihood (assuming normal errors)
        if self.mse_residual > 0:
            self.log_likelihood = -0.5 * n * (np.log(2 * np.pi) + np.log(self.mse_residual) + 1)
        else:
            self.log_likelihood = np.inf

        # AIC and BIC
        self.aic = -2 * self.log_likelihood + 2 * (k + 1)
        self.bic = -2 * self.log_likelihood + np.log(n) * (k + 1)

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

        # Jarque-Bera test for normality
        try:
            from scipy.stats import jarque_bera
            self.jb_statistic, self.jb_pvalue = jarque_bera(residuals)
        except:
            self.jb_statistic = np.nan
            self.jb_pvalue = np.nan

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
        summary.append(f"Omnibus:              {self.jb_statistic:>16.3f}   Durbin-Watson:       {self.durbin_watson:>10.3f}")
        summary.append(f"Prob(Omnibus):        {self.jb_pvalue:>16.3f}   Jarque-Bera (JB):    {self.jb_statistic:>10.3f}")
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
