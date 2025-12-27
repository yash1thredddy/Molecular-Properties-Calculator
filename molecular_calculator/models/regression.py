"""3D OLS Regression model.

Provides the ThreeDOLSRegression class for fitting planes to 3D data.
"""

from typing import Dict, Union, Tuple, Optional
import numpy as np


class ThreeDOLSRegression:
    """
    3D Ordinary Least Squares (OLS) Regression Calculator.

    Fits a plane to 3D data: Z = b0 + b1*X + b2*Y
    Based on analytical formulas from partial differentiation of the cost function.

    Attributes:
        x: X values array
        y: Y values array
        z: Z values array (dependent variable)
        b0: Intercept coefficient
        b1: X coefficient
        b2: Y coefficient
        r_squared: Coefficient of determination
        rmse: Root mean squared error
        residuals: Array of residuals (observed - predicted)
        n: Number of data points

    Example:
        >>> x = np.array([1, 2, 3, 4, 5])
        >>> y = np.array([2, 3, 4, 5, 6])
        >>> z = np.array([3, 5, 7, 9, 11])
        >>> model = ThreeDOLSRegression(x, y, z)
        >>> print(f"R² = {model.r_squared:.4f}")
        >>> z_pred = model.predict(2.5, 3.5)
    """

    def __init__(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        original_indices: Optional[np.ndarray] = None
    ):
        """
        Initialize 3D OLS regression with data.

        Args:
            x: 1D array of X values
            y: 1D array of Y values
            z: 1D array of Z values (dependent variable)
            original_indices: Optional array of original indices (for tracking valid data points)

        Raises:
            ValueError: If less than 3 valid data points
            ValueError: If X and Y are perfectly collinear
        """
        # Convert to numpy arrays and remove any NaN values
        self.x = np.array(x, dtype=float)
        self.y = np.array(y, dtype=float)
        self.z = np.array(z, dtype=float)

        # Store or create original indices
        if original_indices is not None:
            self.original_indices = np.array(original_indices)
        else:
            self.original_indices = np.arange(len(self.x))

        # Create mask for valid data (no NaN or Inf)
        valid_mask = (
            np.isfinite(self.x) &
            np.isfinite(self.y) &
            np.isfinite(self.z)
        )

        self.x = self.x[valid_mask]
        self.y = self.y[valid_mask]
        self.z = self.z[valid_mask]
        self.valid_indices = self.original_indices[valid_mask]

        self.n = len(self.x)

        if self.n < 3:
            raise ValueError("At least 3 valid data points are required for 3D OLS regression")

        # Calculate coefficients
        self.b0, self.b1, self.b2 = self._calculate_coefficients()

        # Calculate statistics
        self.r_squared = self._calculate_r_squared()
        self.residuals = self._calculate_residuals()
        self.rmse = self._calculate_rmse()

    def _calculate_coefficients(self) -> Tuple[float, float, float]:
        """
        Calculate OLS coefficients using analytical formulas.

        Based on formulas from: https://sapiencespace.com/breaking-down-3d-linear-regression/

        Returns:
            Tuple of (b0, b1, b2) coefficients

        Raises:
            ValueError: If X and Y are perfectly collinear
        """
        # Calculate means
        x_mean = np.mean(self.x)
        y_mean = np.mean(self.y)
        z_mean = np.mean(self.z)

        # Calculate deviations from mean
        x_dev = self.x - x_mean
        y_dev = self.y - y_mean
        z_dev = self.z - z_mean

        # Calculate intermediate sums (Sa, Sb, Sc, Sd, Se)
        Sa = np.sum(x_dev ** 2)  # Σ(x - x_mean)²
        Sb = np.sum(y_dev ** 2)  # Σ(y - y_mean)²
        Sc = np.sum(x_dev * y_dev)  # Σ(x - x_mean)(y - y_mean)
        Sd = np.sum(x_dev * z_dev)  # Σ(x - x_mean)(z - z_mean)
        Se = np.sum(y_dev * z_dev)  # Σ(y - y_mean)(z - z_mean)

        # Calculate denominator for slope coefficients
        denominator = Sa * Sb - Sc ** 2

        if abs(denominator) < 1e-10:
            raise ValueError("X and Y variables are perfectly collinear")

        # Calculate slope coefficients
        b1 = (Sd * Sb - Sc * Se) / denominator
        b2 = (Sa * Se - Sc * Sd) / denominator

        # Calculate intercept
        b0 = z_mean - b1 * x_mean - b2 * y_mean

        return b0, b1, b2

    def _calculate_residuals(self) -> np.ndarray:
        """Calculate residuals (observed - predicted)."""
        z_predicted = self.predict(self.x, self.y)
        return self.z - z_predicted

    def _calculate_r_squared(self) -> float:
        """
        Calculate R² (coefficient of determination).

        R² = 1 - (SS_res / SS_tot)
        where SS_res = Σ(observed - predicted)²
              SS_tot = Σ(observed - mean)²

        Returns:
            R² value between -∞ and 1 (closer to 1 is better)
        """
        z_predicted = self.predict(self.x, self.y)
        z_mean = np.mean(self.z)

        ss_res = np.sum((self.z - z_predicted) ** 2)
        ss_tot = np.sum((self.z - z_mean) ** 2)

        if ss_tot < 1e-10:
            return 1.0 if ss_res < 1e-10 else 0.0

        return 1 - (ss_res / ss_tot)

    def _calculate_rmse(self) -> float:
        """
        Calculate Root Mean Squared Error.

        RMSE = sqrt(Σ(observed - predicted)² / n)

        Returns:
            RMSE value
        """
        return np.sqrt(np.mean(self.residuals ** 2))

    def predict(
        self,
        x: Union[float, np.ndarray],
        y: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """
        Predict Z values using the fitted plane equation.

        Args:
            x: X value(s) for prediction
            y: Y value(s) for prediction

        Returns:
            Predicted Z value(s)
        """
        return self.b0 + self.b1 * np.array(x) + self.b2 * np.array(y)

    def get_equation_string(self, decimals: int = 4) -> str:
        """
        Get the plane equation as a formatted string.

        Args:
            decimals: Number of decimal places for coefficients

        Returns:
            Equation string in format: Z = b0 + b1*X + b2*Y
        """
        b0_str = f"{self.b0:.{decimals}f}"
        b1_str = f"{abs(self.b1):.{decimals}f}"
        b2_str = f"{abs(self.b2):.{decimals}f}"

        b1_sign = "+" if self.b1 >= 0 else "-"
        b2_sign = "+" if self.b2 >= 0 else "-"

        return f"Z = {b0_str} {b1_sign} {b1_str}·X {b2_sign} {b2_str}·Y"

    def get_statistics(self) -> Dict[str, float]:
        """
        Get regression statistics.

        Returns:
            Dictionary containing regression statistics
        """
        return {
            'b0': self.b0,
            'b1': self.b1,
            'b2': self.b2,
            'R²': self.r_squared,
            'RMSE': self.rmse,
            'n': self.n
        }

    def get_plane_mesh(self, num_points: int = 20) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate mesh grid for plotting the fitted plane.

        Args:
            num_points: Number of points in each dimension for the mesh

        Returns:
            Tuple of (X_mesh, Y_mesh, Z_mesh) for 3D plotting
        """
        x_min, x_max = self.x.min(), self.x.max()
        y_min, y_max = self.y.min(), self.y.max()

        x_range = x_max - x_min
        y_range = y_max - y_min

        # Add 10% padding
        x_min -= 0.1 * x_range
        x_max += 0.1 * x_range
        y_min -= 0.1 * y_range
        y_max += 0.1 * y_range

        X_mesh = np.linspace(x_min, x_max, num_points)
        Y_mesh = np.linspace(y_min, y_max, num_points)
        X_mesh, Y_mesh = np.meshgrid(X_mesh, Y_mesh)

        Z_mesh = self.predict(X_mesh, Y_mesh)

        return X_mesh, Y_mesh, Z_mesh
