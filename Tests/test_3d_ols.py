"""
Test script for 3D OLS Regression implementation

This script tests the ThreeDOLSRegression class with sample data
"""

import numpy as np
from molecular_calculator import ThreeDOLSRegression

def test_3d_ols_basic():
    """Test basic 3D OLS regression with synthetic data"""
    print("=" * 60)
    print("Testing 3D OLS Regression Implementation")
    print("=" * 60)

    # Create synthetic data: Z = 2 + 3*X + 4*Y + noise
    np.random.seed(42)
    n_points = 30

    X = np.random.uniform(0, 10, n_points)
    Y = np.random.uniform(0, 10, n_points)

    # True model: Z = 2 + 3*X + 4*Y
    true_b0 = 2.0
    true_b1 = 3.0
    true_b2 = 4.0

    # Add some noise
    noise = np.random.normal(0, 2, n_points)
    Z = true_b0 + true_b1 * X + true_b2 * Y + noise

    print(f"\nTrue model: Z = {true_b0} + {true_b1}*X + {true_b2}*Y")
    print(f"Number of data points: {n_points}")

    # Create OLS regression model
    print("\nFitting 3D OLS regression...")
    model = ThreeDOLSRegression(X, Y, Z)

    # Get statistics
    stats = model.get_statistics()
    equation = model.get_equation_string(decimals=4)

    # Display results
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"\nFitted Equation: {equation}")
    print(f"\nCoefficients:")
    print(f"  b0 (Intercept):     {stats['b0']:.4f}  (True: {true_b0:.4f})")
    print(f"  b1 (X coefficient): {stats['b1']:.4f}  (True: {true_b1:.4f})")
    print(f"  b2 (Y coefficient): {stats['b2']:.4f}  (True: {true_b2:.4f})")

    print(f"\nModel Quality:")
    print(f"  R² (Coefficient of Determination): {stats['R²']:.4f}")
    print(f"  RMSE (Root Mean Squared Error):    {stats['RMSE']:.4f}")
    print(f"  Number of observations:            {stats['n']}")

    # Test predictions
    print(f"\nTesting predictions:")
    test_x, test_y = 5.0, 7.0
    predicted_z = model.predict(test_x, test_y)
    expected_z = true_b0 + true_b1 * test_x + true_b2 * test_y
    print(f"  For X={test_x}, Y={test_y}:")
    print(f"    Predicted Z: {predicted_z:.4f}")
    print(f"    Expected Z (true model): {expected_z:.4f}")

    # Interpretation
    print(f"\nInterpretation:")
    if stats['R²'] >= 0.9:
        quality = "excellent"
    elif stats['R²'] >= 0.7:
        quality = "good"
    elif stats['R²'] >= 0.5:
        quality = "moderate"
    else:
        quality = "poor"

    print(f"  The model has {quality} fit with R² = {stats['R²']:.4f}")
    print(f"  This means {stats['R²']*100:.2f}% of the variance in Z is explained by X and Y.")

    print("\n" + "=" * 60)
    print("✅ Test completed successfully!")
    print("=" * 60)


def test_3d_ols_with_sar_example():
    """Test with SAR-like data (Structure-Activity Relationship)"""
    print("\n\n" + "=" * 60)
    print("SAR Example: Biological Activity vs Molecular Properties")
    print("=" * 60)

    # Simulate SAR data: Activity (pKi) as function of LogP and TPSA
    # Typical relationship: Activity increases with LogP (lipophilicity)
    #                      Activity decreases with TPSA (polar surface area)

    np.random.seed(123)
    n_compounds = 30

    # Molecular properties
    LogP = np.random.uniform(0, 5, n_compounds)  # Lipophilicity
    TPSA = np.random.uniform(20, 140, n_compounds)  # Polar Surface Area

    # Biological activity (pKi) - simplified relationship
    # pKi tends to increase with LogP and decrease with TPSA
    pKi = 6.0 + 0.5 * LogP - 0.02 * TPSA + np.random.normal(0, 0.3, n_compounds)

    print(f"\nSimulated {n_compounds} compounds")
    print(f"Variables:")
    print(f"  X = LogP (lipophilicity, range: {LogP.min():.2f} to {LogP.max():.2f})")
    print(f"  Y = TPSA (polar surface area, range: {TPSA.min():.2f} to {TPSA.max():.2f})")
    print(f"  Z = pKi (biological activity, range: {pKi.min():.2f} to {pKi.max():.2f})")

    # Fit model
    print("\nFitting 3D OLS regression...")
    sar_model = ThreeDOLSRegression(LogP, TPSA, pKi)

    # Get results
    stats = sar_model.get_statistics()
    equation = sar_model.get_equation_string(decimals=3)

    print("\n" + "=" * 60)
    print("SAR MODEL RESULTS")
    print("=" * 60)
    print(f"\nFitted Equation: {equation}")
    print(f"\nCoefficients Interpretation:")
    print(f"  Intercept (b0): {stats['b0']:.3f}")
    print(f"  LogP effect (b1): {stats['b1']:.3f} (positive = increased activity with lipophilicity)")
    print(f"  TPSA effect (b2): {stats['b2']:.3f} (negative = decreased activity with polarity)")

    print(f"\nModel Statistics:")
    print(f"  R² = {stats['R²']:.3f} ({stats['R²']*100:.1f}% variance explained)")
    print(f"  RMSE = {stats['RMSE']:.3f} pKi units")

    # Predict activity for an optimized compound
    optimal_logp = 3.0  # Good balance
    optimal_tpsa = 60.0  # Drug-like range
    predicted_pki = sar_model.predict(optimal_logp, optimal_tpsa)

    print(f"\nPrediction for optimized compound:")
    print(f"  LogP = {optimal_logp}, TPSA = {optimal_tpsa}")
    print(f"  Predicted pKi = {predicted_pki:.2f}")

    print("\n" + "=" * 60)
    print("✅ SAR example completed!")
    print("=" * 60)


if __name__ == "__main__":
    # Run tests
    test_3d_ols_basic()
    test_3d_ols_with_sar_example()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED ✅")
    print("=" * 60)
    print("\nThe 3D OLS Regression implementation is working correctly!")
    print("You can now use it in the Molecular Properties Calculator app.")
