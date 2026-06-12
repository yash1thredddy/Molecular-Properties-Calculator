"""Services layer for molecular calculator.

This module exports service classes for API access, format conversion,
and property calculation.
"""

from .api_client import (
    APIResponse,
    ChemicalAPIClient,
    get_api_client,
)
from .conversion_service import (
    ConversionService,
    get_conversion_service,
)
from .property_service import (
    PropertyCalculator,
    get_property_calculator,
)
from .ligand_efficiency import (
    DependencyChecker,
    LigandEfficiencyCalculator,
    get_lei_descriptions,
)
from .assay_interference import (
    InterferenceFlags,
    calculate_interference_flags,
    calculate_batch_interference_flags,
    get_interference_flags_from_smiles,
    get_interference_summary,
    FLAG_DESCRIPTIONS,
)

__all__ = [
    # API Client
    "APIResponse",
    "ChemicalAPIClient",
    "get_api_client",
    # Conversion Service
    "ConversionService",
    "get_conversion_service",
    # Property Service
    "PropertyCalculator",
    "get_property_calculator",
    # Ligand Efficiency
    "DependencyChecker",
    "LigandEfficiencyCalculator",
    "get_lei_descriptions",
    # Assay Interference
    "InterferenceFlags",
    "calculate_interference_flags",
    "calculate_batch_interference_flags",
    "get_interference_flags_from_smiles",
    "get_interference_summary",
    "FLAG_DESCRIPTIONS",
]
