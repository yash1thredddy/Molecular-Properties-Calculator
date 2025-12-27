"""Data models for molecular calculator.

This module exports all data classes and enums for molecular structures,
properties, and calculation results.
"""

from .molecule import (
    InputFormat,
    PropertyGroup,
    MoleculeInput,
    MolecularProperties,
    CalculationResult,
    ConversionResult,
    LigandEfficiencyIndices,
    PROPERTY_GROUPS,
    LEI_PROPERTY_GROUP,
)
from .regression import ThreeDOLSRegression
from .explanations import PropertyExplanations
from .regression_3d import RegressionSummary, perform_3d_regression, suggest_best_3d_pairs

__all__ = [
    # Enums
    "InputFormat",
    "PropertyGroup",
    # Data classes
    "MoleculeInput",
    "MolecularProperties",
    "CalculationResult",
    "ConversionResult",
    "LigandEfficiencyIndices",
    # Regression
    "ThreeDOLSRegression",
    "RegressionSummary",
    "perform_3d_regression",
    "suggest_best_3d_pairs",
    # Documentation
    "PropertyExplanations",
    # Constants
    "PROPERTY_GROUPS",
    "LEI_PROPERTY_GROUP",
]
