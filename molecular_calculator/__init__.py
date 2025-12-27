"""Molecular Properties Calculator Package.

A production-ready application for calculating and visualizing
molecular properties from SMILES, InChI, and InChI Key formats.

Example:
    >>> from molecular_calculator import MolecularCalculator
    >>> calc = MolecularCalculator()
    >>> props = calc.calculate_molecular_properties("CCO")
    >>> print(props['Molecular_Weight'])
    46.069
"""

__version__ = "2.0.0"
__author__ = "Development Team"

# Core exports
from molecular_calculator.core import MolecularCalculator

# Model exports
from molecular_calculator.models import (
    InputFormat,
    MolecularProperties,
    CalculationResult,
    ConversionResult,
    PROPERTY_GROUPS,
    ThreeDOLSRegression,
    PropertyExplanations,
)

# Service exports
from molecular_calculator.services import (
    ConversionService,
    PropertyCalculator,
    ChemicalAPIClient,
)

# Utility exports
from molecular_calculator.utils import (
    InputValidator,
    FileValidator,
    DataFrameValidator,
    SessionState,
    InvalidSMILESError,
    PropertyCalculationError,
    ConversionError,
)

__all__ = [
    # Version info
    "__version__",
    "__author__",
    # Core
    "MolecularCalculator",
    # Models
    "InputFormat",
    "MolecularProperties",
    "CalculationResult",
    "ConversionResult",
    "PROPERTY_GROUPS",
    "ThreeDOLSRegression",
    "PropertyExplanations",
    # Services
    "ConversionService",
    "PropertyCalculator",
    "ChemicalAPIClient",
    # Utils
    "InputValidator",
    "FileValidator",
    "DataFrameValidator",
    "SessionState",
    "InvalidSMILESError",
    "PropertyCalculationError",
    "ConversionError",
]
