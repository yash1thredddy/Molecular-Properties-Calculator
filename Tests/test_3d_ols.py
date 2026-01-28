"""
Test script for 3D OLS Regression implementation

This script tests the ThreeDOLSRegression class with sample data.
All tests use proper assertions that will fail if results are incorrect.
"""

import unittest
import numpy as np
from molecular_calculator import ThreeDOLSRegression


class TestThreeDOLSBasic(unittest.TestCase):
    """Test basic 3D OLS regression functionality with synthetic data."""

    def setUp(self):
        """Set up test data with known true model: Z = 2 + 3*X + 4*Y + noise."""
        np.random.seed(42)
        self.n_points = 30

        self.X = np.random.uniform(0, 10, self.n_points)
        self.Y = np.random.uniform(0, 10, self.n_points)

        # True model coefficients
        self.true_b0 = 2.0
        self.true_b1 = 3.0
        self.true_b2 = 4.0

        # Add controlled noise
        noise = np.random.normal(0, 2, self.n_points)
        self.Z = self.true_b0 + self.true_b1 * self.X + self.true_b2 * self.Y + noise

    def test_model_fits_successfully(self):
        """Test that the model fits without errors."""
        model = ThreeDOLSRegression(self.X, self.Y, self.Z)
        self.assertIsNotNone(model)

    def test_coefficients_near_true_values(self):
        """Test that fitted coefficients are reasonably close to true values.

        With random noise (std=2) on 30 points, we allow some tolerance.
        Expected values come from the known generating model: Z = 2 + 3*X + 4*Y + noise.
        """
        model = ThreeDOLSRegression(self.X, self.Y, self.Z)
        stats = model.get_statistics()

        # Coefficients should be within 1.5 of true values given noise level
        # These tolerances account for random sampling variation
        self.assertAlmostEqual(
            stats['b0'], self.true_b0, delta=2.0,
            msg=f"Intercept {stats['b0']:.4f} too far from true value {self.true_b0}"
        )
        self.assertAlmostEqual(
            stats['b1'], self.true_b1, delta=1.0,
            msg=f"X coefficient {stats['b1']:.4f} too far from true value {self.true_b1}"
        )
        self.assertAlmostEqual(
            stats['b2'], self.true_b2, delta=1.0,
            msg=f"Y coefficient {stats['b2']:.4f} too far from true value {self.true_b2}"
        )

    def test_r_squared_is_high(self):
        """Test that R-squared indicates good model fit.

        For synthetic data with moderate noise, R-squared should be > 0.9.
        """
        model = ThreeDOLSRegression(self.X, self.Y, self.Z)
        stats = model.get_statistics()

        self.assertGreaterEqual(
            stats['R\u00b2'], 0.90,
            msg=f"R-squared {stats['R\u00b2']:.4f} indicates poor fit (expected >= 0.90)"
        )
        self.assertLessEqual(
            stats['R\u00b2'], 1.0,
            msg=f"R-squared {stats['R\u00b2']:.4f} exceeds maximum of 1.0"
        )

    def test_rmse_is_reasonable(self):
        """Test that RMSE is in expected range given noise level."""
        model = ThreeDOLSRegression(self.X, self.Y, self.Z)
        stats = model.get_statistics()

        # With noise std=2, RMSE should be around 2 (slightly less due to fitting)
        self.assertGreater(
            stats['RMSE'], 0,
            msg="RMSE should be positive"
        )
        self.assertLess(
            stats['RMSE'], 5,
            msg=f"RMSE {stats['RMSE']:.4f} is unexpectedly high"
        )

    def test_n_observations_correct(self):
        """Test that number of observations is recorded correctly."""
        model = ThreeDOLSRegression(self.X, self.Y, self.Z)
        stats = model.get_statistics()

        self.assertEqual(
            stats['n'], self.n_points,
            msg=f"Expected {self.n_points} observations, got {stats['n']}"
        )

    def test_prediction_accuracy(self):
        """Test that predictions are close to expected values."""
        model = ThreeDOLSRegression(self.X, self.Y, self.Z)

        # Test prediction at a specific point
        test_x, test_y = 5.0, 7.0
        predicted_z = model.predict(test_x, test_y)

        # Expected value from the true model (without noise)
        expected_z = self.true_b0 + self.true_b1 * test_x + self.true_b2 * test_y

        # Prediction should be reasonably close (within 5 units given noise)
        self.assertAlmostEqual(
            predicted_z, expected_z, delta=5.0,
            msg=f"Prediction {predicted_z:.4f} too far from expected {expected_z:.4f}"
        )

    def test_equation_string_format(self):
        """Test that equation string is properly formatted."""
        model = ThreeDOLSRegression(self.X, self.Y, self.Z)
        equation = model.get_equation_string(decimals=4)

        self.assertIsInstance(equation, str)
        self.assertIn('X', equation, msg="Equation should contain X variable")
        self.assertIn('Y', equation, msg="Equation should contain Y variable")
        # Should contain numeric coefficients
        self.assertRegex(equation, r'\d+\.\d+', msg="Equation should contain decimal numbers")


class TestThreeDOLSSARExample(unittest.TestCase):
    """Test 3D OLS with SAR-like data (Structure-Activity Relationship)."""

    def setUp(self):
        """Set up SAR-like test data.

        Simulates a typical SAR relationship where biological activity (pKi)
        depends on molecular properties (LogP and TPSA).
        Model: pKi = 6.0 + 0.5*LogP - 0.02*TPSA + noise

        This is a common pattern in drug discovery where:
        - Activity increases with lipophilicity (positive LogP coefficient)
        - Activity decreases with polar surface area (negative TPSA coefficient)
        """
        np.random.seed(123)
        self.n_compounds = 30

        # Molecular properties in typical drug-like ranges
        self.LogP = np.random.uniform(0, 5, self.n_compounds)  # Lipophilicity
        self.TPSA = np.random.uniform(20, 140, self.n_compounds)  # Polar Surface Area

        # True model coefficients for SAR
        self.true_intercept = 6.0
        self.true_logp_coef = 0.5
        self.true_tpsa_coef = -0.02

        # Biological activity with noise
        noise = np.random.normal(0, 0.3, self.n_compounds)
        self.pKi = (self.true_intercept +
                    self.true_logp_coef * self.LogP +
                    self.true_tpsa_coef * self.TPSA +
                    noise)

    def test_sar_model_fits(self):
        """Test that SAR model fits without errors."""
        model = ThreeDOLSRegression(self.LogP, self.TPSA, self.pKi)
        self.assertIsNotNone(model)

    def test_sar_coefficient_signs(self):
        """Test that coefficient signs match expected SAR relationship.

        In drug discovery, we typically expect:
        - Positive LogP coefficient (more lipophilic = higher activity)
        - Negative TPSA coefficient (more polar = lower activity)
        """
        model = ThreeDOLSRegression(self.LogP, self.TPSA, self.pKi)
        stats = model.get_statistics()

        self.assertGreater(
            stats['b1'], 0,
            msg=f"LogP coefficient should be positive, got {stats['b1']:.4f}"
        )
        self.assertLess(
            stats['b2'], 0,
            msg=f"TPSA coefficient should be negative, got {stats['b2']:.4f}"
        )

    def test_sar_coefficient_magnitudes(self):
        """Test that coefficient magnitudes are close to true values.

        True coefficients: LogP=0.5, TPSA=-0.02
        """
        model = ThreeDOLSRegression(self.LogP, self.TPSA, self.pKi)
        stats = model.get_statistics()

        self.assertAlmostEqual(
            stats['b1'], self.true_logp_coef, delta=0.2,
            msg=f"LogP coefficient {stats['b1']:.4f} too far from true value {self.true_logp_coef}"
        )
        self.assertAlmostEqual(
            stats['b2'], self.true_tpsa_coef, delta=0.01,
            msg=f"TPSA coefficient {stats['b2']:.4f} too far from true value {self.true_tpsa_coef}"
        )

    def test_sar_r_squared_excellent(self):
        """Test that SAR model has excellent fit (low noise data)."""
        model = ThreeDOLSRegression(self.LogP, self.TPSA, self.pKi)
        stats = model.get_statistics()

        self.assertGreaterEqual(
            stats['R\u00b2'], 0.85,
            msg=f"R-squared {stats['R\u00b2']:.4f} indicates suboptimal fit for low-noise SAR data"
        )

    def test_sar_prediction_drug_like_compound(self):
        """Test prediction for a hypothetical optimized compound.

        Predicts activity for a compound with:
        - LogP = 3.0 (good balance of lipophilicity)
        - TPSA = 60 (within drug-like range)
        """
        model = ThreeDOLSRegression(self.LogP, self.TPSA, self.pKi)

        optimal_logp = 3.0
        optimal_tpsa = 60.0
        predicted_pki = model.predict(optimal_logp, optimal_tpsa)

        # Expected from true model: 6.0 + 0.5*3.0 - 0.02*60.0 = 6.3
        expected_pki = (self.true_intercept +
                        self.true_logp_coef * optimal_logp +
                        self.true_tpsa_coef * optimal_tpsa)

        self.assertAlmostEqual(
            predicted_pki, expected_pki, delta=0.5,
            msg=f"Predicted pKi {predicted_pki:.2f} too far from expected {expected_pki:.2f}"
        )

    def test_sar_rmse_low(self):
        """Test that RMSE is low for SAR data with controlled noise."""
        model = ThreeDOLSRegression(self.LogP, self.TPSA, self.pKi)
        stats = model.get_statistics()

        # With noise std=0.3, RMSE should be around 0.3
        self.assertLess(
            stats['RMSE'], 0.5,
            msg=f"RMSE {stats['RMSE']:.4f} is too high for low-noise SAR data"
        )


class TestThreeDOLSEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions for 3D OLS regression."""

    def test_perfect_fit_no_noise(self):
        """Test with perfectly linear data (no noise)."""
        X = np.array([1, 2, 3, 4, 5])
        Y = np.array([2, 4, 5, 4, 6])
        Z = 1 + 2*X + 3*Y  # Perfect linear relationship

        model = ThreeDOLSRegression(X, Y, Z)
        stats = model.get_statistics()

        # R-squared should be exactly 1.0 (or very close)
        self.assertAlmostEqual(
            stats['R\u00b2'], 1.0, places=10,
            msg="Perfect linear data should have R-squared = 1.0"
        )

        # RMSE should be 0 (or very close)
        self.assertAlmostEqual(
            stats['RMSE'], 0.0, places=10,
            msg="Perfect linear data should have RMSE = 0.0"
        )

    def test_minimum_data_points(self):
        """Test with minimum required data points (4 for 3D regression)."""
        # X and Y must NOT be collinear for valid 3D regression
        X = np.array([1, 2, 3, 4])
        Y = np.array([4, 2, 3, 1])  # Different from X to avoid collinearity
        Z = np.array([2, 4, 6, 8])

        # Should not raise an error
        model = ThreeDOLSRegression(X, Y, Z)
        stats = model.get_statistics()

        self.assertEqual(stats['n'], 4)

    def test_large_values(self):
        """Test with large numerical values."""
        np.random.seed(42)
        n = 20
        X = np.random.uniform(1e6, 1e7, n)
        Y = np.random.uniform(1e6, 1e7, n)
        Z = 1e6 + 0.1*X + 0.2*Y + np.random.normal(0, 1e4, n)

        model = ThreeDOLSRegression(X, Y, Z)
        stats = model.get_statistics()

        # Should still compute valid statistics
        self.assertGreater(stats['R\u00b2'], 0)
        self.assertTrue(np.isfinite(stats['RMSE']))

    def test_small_values(self):
        """Test with small numerical values."""
        np.random.seed(42)
        n = 20
        # Use values in 0.01-0.1 range to avoid numerical precision issues
        # while still testing "small" values relative to typical data
        X = np.random.uniform(0.01, 0.1, n)
        np.random.seed(43)  # Different seed for Y
        Y = np.random.uniform(0.01, 0.1, n)
        Z = 0.05 + 2*X + 3*Y + np.random.normal(0, 0.001, n)

        model = ThreeDOLSRegression(X, Y, Z)
        stats = model.get_statistics()

        # Should still compute valid statistics
        self.assertTrue(np.isfinite(stats['R²']))
        self.assertTrue(np.isfinite(stats['RMSE']))

    def test_negative_values(self):
        """Test with negative values in data."""
        np.random.seed(42)
        n = 20
        X = np.random.uniform(-10, 10, n)
        Y = np.random.uniform(-10, 10, n)
        Z = 5 + 2*X - 3*Y + np.random.normal(0, 1, n)

        model = ThreeDOLSRegression(X, Y, Z)
        stats = model.get_statistics()

        # Should handle negative values correctly
        self.assertGreater(stats['R\u00b2'], 0.8)
        self.assertTrue(np.isfinite(stats['b0']))
        self.assertTrue(np.isfinite(stats['b1']))
        self.assertTrue(np.isfinite(stats['b2']))


class TestThreeDOLSStatistics(unittest.TestCase):
    """Test statistical properties of 3D OLS regression."""

    def test_statistics_keys_present(self):
        """Test that all expected statistics keys are returned."""
        X = np.array([1, 2, 3, 4, 5])
        Y = np.array([2, 4, 5, 4, 6])
        Z = 1 + 2*X + 3*Y

        model = ThreeDOLSRegression(X, Y, Z)
        stats = model.get_statistics()

        required_keys = ['b0', 'b1', 'b2', 'R\u00b2', 'RMSE', 'n']
        for key in required_keys:
            self.assertIn(
                key, stats,
                msg=f"Missing required statistic: {key}"
            )

    def test_r_squared_bounds(self):
        """Test that R-squared is always between 0 and 1."""
        # Run multiple tests with different random data
        for seed in range(5):
            np.random.seed(seed)
            n = 30
            X = np.random.uniform(0, 10, n)
            Y = np.random.uniform(0, 10, n)
            Z = np.random.uniform(0, 10, n)  # Random Z (low correlation expected)

            model = ThreeDOLSRegression(X, Y, Z)
            stats = model.get_statistics()

            self.assertGreaterEqual(
                stats['R\u00b2'], 0,
                msg=f"R-squared should be >= 0, got {stats['R\u00b2']}"
            )
            self.assertLessEqual(
                stats['R\u00b2'], 1,
                msg=f"R-squared should be <= 1, got {stats['R\u00b2']}"
            )

    def test_rmse_non_negative(self):
        """Test that RMSE is always non-negative."""
        np.random.seed(42)
        n = 30
        X = np.random.uniform(0, 10, n)
        Y = np.random.uniform(0, 10, n)
        Z = np.random.uniform(0, 10, n)

        model = ThreeDOLSRegression(X, Y, Z)
        stats = model.get_statistics()

        self.assertGreaterEqual(
            stats['RMSE'], 0,
            msg="RMSE should be non-negative"
        )


if __name__ == '__main__':
    unittest.main(verbosity=2)
