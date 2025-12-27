"""Application configuration and constants.

This module provides centralized configuration for the application.
All magic numbers and configuration values should be defined here.

Environment variables can override defaults (read at module import time):
- MAX_FILE_SIZE_MB: Override max file size
- API_TIMEOUT_SECONDS: Override API timeout
- APP_VERSION: Override version string
- MAX_ROWS_LIMIT: Override max rows limit
- CACHE_TTL_SECONDS: Override cache TTL

Note:
    Environment variables are read when this module is first imported,
    not when AppConfig instances are created. To change configuration
    at runtime, set environment variables before importing this module.
"""

import os
from dataclasses import dataclass, field
from typing import FrozenSet


def _get_int_env(name: str, default: int) -> int:
    """Get integer environment variable or return default."""
    val = os.getenv(name)
    if val is not None:
        try:
            return int(val)
        except ValueError:
            pass
    return default


def _get_str_env(name: str, default: str) -> str:
    """Get string environment variable or return default."""
    return os.getenv(name, default)


@dataclass(frozen=True)
class AppConfig:
    """Immutable application configuration.

    This class uses frozen=True to ensure configuration values
    cannot be accidentally modified at runtime.

    Important:
        Environment variables are read at module import time (when
        default_factory lambdas are evaluated), not at instance creation.
        Set environment variables before importing this module.

    Example:
        >>> import os
        >>> os.environ['MAX_FILE_SIZE_MB'] = '100'  # Before import!
        >>> from molecular_calculator.config import config
        >>> print(config.MAX_FILE_SIZE_MB)
        100
    """

    # Application info
    APP_NAME: str = "Molecular Properties Calculator"
    APP_VERSION: str = field(
        default_factory=lambda: _get_str_env('APP_VERSION', "2.0.0")
    )
    APP_ICON: str = "🧪"

    # File handling
    MAX_FILE_SIZE_MB: int = field(
        default_factory=lambda: _get_int_env('MAX_FILE_SIZE_MB', 50)
    )
    ALLOWED_EXTENSIONS: FrozenSet[str] = frozenset({'.csv', '.xlsx'})

    # Data limits
    MAX_ROWS_WARNING: int = 10_000
    MAX_ROWS_LIMIT: int = field(
        default_factory=lambda: _get_int_env('MAX_ROWS_LIMIT', 100_000)
    )
    MAX_SMILES_LENGTH: int = 10_000

    # Chart defaults
    DEFAULT_MARKER_SIZE: int = 8
    MIN_MARKER_SIZE: int = 2
    MAX_MARKER_SIZE: int = 50
    DEFAULT_CHART_HEIGHT: int = 600

    # Categorical threshold
    MAX_CATEGORICAL_CARDINALITY: int = 50

    # API settings
    API_TIMEOUT_SECONDS: int = field(
        default_factory=lambda: _get_int_env('API_TIMEOUT_SECONDS', 10)
    )
    MAX_RETRY_ATTEMPTS: int = 3

    # Rate limiting
    NIH_API_RATE_LIMIT: int = 30  # requests per minute
    PUBCHEM_API_RATE_LIMIT: int = 5  # requests per second

    # External API URLs
    NIH_CIR_URL: str = "https://cactus.nci.nih.gov/chemical/structure"
    PUBCHEM_URL: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

    # Cache settings
    CACHE_TTL_SECONDS: int = field(
        default_factory=lambda: _get_int_env('CACHE_TTL_SECONDS', 3600)
    )
    MAX_CACHE_SIZE: int = 1000

    @property
    def MAX_FILE_SIZE_BYTES(self) -> int:
        """Get maximum file size in bytes."""
        return self.MAX_FILE_SIZE_MB * 1024 * 1024


def _create_config() -> AppConfig:
    """Create configuration instance (environment overrides applied at construction)."""
    return AppConfig()


# Global configuration instance
config = _create_config()


# NOTE: PROPERTY_GROUPS has been moved to molecular_calculator.models.molecule
# Import from there if needed: from molecular_calculator.models import PROPERTY_GROUPS

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
