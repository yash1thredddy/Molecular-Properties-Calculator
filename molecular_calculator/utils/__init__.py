"""Utilities package."""

from .exceptions import (
    MolecularCalculatorError,
    InvalidSMILESError,
    PropertyCalculationError,
    FileValidationError,
    APIError,
    RateLimitError,
    ConversionError,
)
from .validators import InputValidator, FileValidator, DataFrameValidator, ValidationResult
from .session_state import SessionState
from .export import (
    export_chart,
    export_chart_to_bytes,
    get_download_data,
    create_export_buttons,
    EXPORT_FORMATS,
)
from .suggestions import (
    detect_smiles_column,
    detect_id_column,
    detect_column_type,
    suggest_correlated_pairs,
    suggest_regression_variables,
    recommend_visualization_variables,
    get_column_stats,
)
from .error_handler import (
    handle_errors,
    handle_calculation_errors,
    handle_api_errors,
    handle_file_errors,
    error_boundary,
    suppress_warnings,
    get_error_details,
)
from .rate_limiter import (
    RateLimiter,
    nih_limiter,
    pubchem_limiter,
    rate_limited,
    get_limiter_status,
    reset_all_limiters,
)
from .cache import (
    TTLCache,
    cached,
    cache_molecule,
    cache_conversion,
    molecule_cache,
    conversion_cache,
    api_cache,
    get_all_cache_stats,
    clear_all_caches,
)
from .sanitizer import (
    sanitize_smiles,
    sanitize_inchi,
    sanitize_inchi_key,
    sanitize_column_name,
    sanitize_filename,
    sanitize_html,
    sanitize_numeric,
    sanitize_integer,
    is_valid_smiles_format,
    is_valid_inchi_key,
    detect_identifier_type,
)
from .monitoring import (
    PerformanceMetrics,
    app_metrics,
    check_rdkit_available,
    check_dependencies,
    get_system_info,
    get_health_status,
    track_performance,
    setup_logging,
    log_operation,
)

__all__ = [
    # Exceptions
    "MolecularCalculatorError",
    "InvalidSMILESError",
    "PropertyCalculationError",
    "FileValidationError",
    "APIError",
    "RateLimitError",
    "ConversionError",
    # Validators
    "InputValidator",
    "FileValidator",
    "DataFrameValidator",
    "ValidationResult",
    # Session State
    "SessionState",
    # Export utilities
    "export_chart",
    "export_chart_to_bytes",
    "get_download_data",
    "create_export_buttons",
    "EXPORT_FORMATS",
    # Suggestions
    "detect_smiles_column",
    "detect_id_column",
    "detect_column_type",
    "suggest_correlated_pairs",
    "suggest_regression_variables",
    "recommend_visualization_variables",
    "get_column_stats",
    # Error handling
    "handle_errors",
    "handle_calculation_errors",
    "handle_api_errors",
    "handle_file_errors",
    "error_boundary",
    "suppress_warnings",
    "get_error_details",
    # Rate limiting
    "RateLimiter",
    "nih_limiter",
    "pubchem_limiter",
    "rate_limited",
    "get_limiter_status",
    "reset_all_limiters",
    # Caching
    "TTLCache",
    "cached",
    "cache_molecule",
    "cache_conversion",
    "molecule_cache",
    "conversion_cache",
    "api_cache",
    "get_all_cache_stats",
    "clear_all_caches",
    # Sanitization
    "sanitize_smiles",
    "sanitize_inchi",
    "sanitize_inchi_key",
    "sanitize_column_name",
    "sanitize_filename",
    "sanitize_html",
    "sanitize_numeric",
    "sanitize_integer",
    "is_valid_smiles_format",
    "is_valid_inchi_key",
    "detect_identifier_type",
    # Monitoring
    "PerformanceMetrics",
    "app_metrics",
    "check_rdkit_available",
    "check_dependencies",
    "get_system_info",
    "get_health_status",
    "track_performance",
    "setup_logging",
    "log_operation",
]
