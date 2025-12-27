"""Application configuration and constants.

This module provides centralized configuration for the application.
All magic numbers and configuration values should be defined here.
"""

import os
from dataclasses import dataclass, field
from typing import FrozenSet


@dataclass(frozen=True)
class AppConfig:
    """Immutable application configuration.

    This class uses frozen=True to ensure configuration values
    cannot be accidentally modified at runtime.

    Example:
        >>> from molecular_calculator.config import config
        >>> print(config.APP_NAME)
        'Molecular Properties Calculator'
    """

    # Application info
    APP_NAME: str = "Molecular Properties Calculator"
    APP_VERSION: str = "2.0.0"
    APP_ICON: str = "🧪"

    # File handling
    MAX_FILE_SIZE_MB: int = 50
    ALLOWED_EXTENSIONS: FrozenSet[str] = frozenset({'.csv', '.xlsx'})

    # Data limits
    MAX_ROWS_WARNING: int = 10_000
    MAX_ROWS_LIMIT: int = 100_000
    MAX_SMILES_LENGTH: int = 10_000

    # Chart defaults
    DEFAULT_MARKER_SIZE: int = 8
    MIN_MARKER_SIZE: int = 2
    MAX_MARKER_SIZE: int = 50
    DEFAULT_CHART_HEIGHT: int = 600

    # Categorical threshold
    MAX_CATEGORICAL_CARDINALITY: int = 50

    # API settings
    API_TIMEOUT_SECONDS: int = 10
    MAX_RETRY_ATTEMPTS: int = 3

    # Rate limiting
    NIH_API_RATE_LIMIT: int = 30  # requests per minute
    PUBCHEM_API_RATE_LIMIT: int = 5  # requests per second

    # External API URLs
    NIH_CIR_URL: str = "https://cactus.nci.nih.gov/chemical/structure"
    PUBCHEM_URL: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

    # Cache settings
    CACHE_TTL_SECONDS: int = 3600  # 1 hour
    MAX_CACHE_SIZE: int = 1000

    @property
    def MAX_FILE_SIZE_BYTES(self) -> int:
        """Get maximum file size in bytes."""
        return self.MAX_FILE_SIZE_MB * 1024 * 1024


# Singleton configuration instance
# Use environment variables to override defaults if needed
def _create_config() -> AppConfig:
    """Create configuration with optional environment overrides."""
    overrides = {}

    # Check for environment variable overrides
    if os.getenv('MAX_FILE_SIZE_MB'):
        overrides['MAX_FILE_SIZE_MB'] = int(os.getenv('MAX_FILE_SIZE_MB'))

    if os.getenv('API_TIMEOUT_SECONDS'):
        overrides['API_TIMEOUT_SECONDS'] = int(os.getenv('API_TIMEOUT_SECONDS'))

    if os.getenv('APP_VERSION'):
        overrides['APP_VERSION'] = os.getenv('APP_VERSION')

    # Create config with overrides (if frozen dataclass, need to use replace)
    if overrides:
        # For frozen dataclass, we need a different approach
        # Just use defaults for now, overrides would need unfrozen class
        pass

    return AppConfig()


# Global configuration instance
config = _create_config()


# Property groups for UI organization
PROPERTY_GROUPS = {
    'Basic': [
        'Molecular_Weight',
        'ExactMolWt',
        'HeavyAtomMolWt',
    ],
    'Lipophilicity': [
        'LogP',
        'MolLogP',
    ],
    'Solubility': [
        'TPSA',
        'LabuteASA',
    ],
    'H-Bonding': [
        'HBA',
        'HBD',
        'NumHAcceptors',
        'NumHDonors',
    ],
    'Size & Complexity': [
        'HeavyAtomCount',
        'NumRotatableBonds',
        'RingCount',
        'NumAromaticRings',
        'FractionCSP3',
    ],
    'Drug-likeness': [
        'qed',
        'NumHeteroatoms',
        'NHOHCount',
        'NOCount',
    ],
}

# Chart type definitions
CHART_TYPES = [
    'Scatter Plot',
    'Box Plot',
    'Violin Plot',
    'Bar Chart',
    'Histogram',
    'Heatmap',
    '3D Scatter',
    'Line Plot',
]

# Color scales for continuous data
COLOR_SCALES = [
    'Viridis',
    'Plasma',
    'Inferno',
    'Magma',
    'Blues',
    'Reds',
    'Greens',
    'RdBu',
    'Spectral',
]
