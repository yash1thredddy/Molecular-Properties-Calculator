# ==============================================================================
# Monitoring and Health Check Utilities
# ==============================================================================
"""
Application monitoring, health checks, and metrics utilities.
Provides system status and performance tracking.
"""

import time
import logging
import platform
import sys
from threading import Lock
from typing import Dict, Any, Optional, Deque
from datetime import datetime
from functools import wraps
from collections import deque

logger = logging.getLogger(__name__)


# ==============================================================================
# Performance Metrics
# ==============================================================================

class PerformanceMetrics:
    """
    Track performance metrics for operations.

    Thread-safe metrics tracking with locking.

    Example:
        metrics = PerformanceMetrics()

        with metrics.track("calculation"):
            result = expensive_calculation()

        print(metrics.get_stats("calculation"))
    """

    def __init__(self):
        self._metrics: Dict[str, Deque[float]] = {}
        self._counts: Dict[str, int] = {}
        self._errors: Dict[str, int] = {}
        self._lock = Lock()

    def record(self, operation: str, duration: float, success: bool = True) -> None:
        """Record a metric for an operation (thread-safe)."""
        with self._lock:
            if operation not in self._metrics:
                # Use deque with maxlen for automatic size limiting
                self._metrics[operation] = deque(maxlen=1000)
                self._counts[operation] = 0
                self._errors[operation] = 0

            self._metrics[operation].append(duration)
            self._counts[operation] += 1

            if not success:
                self._errors[operation] += 1

    def track(self, operation: str):
        """Context manager to track operation duration."""
        return _MetricsContext(self, operation)

    def get_stats(self, operation: str) -> Dict[str, Any]:
        """Get statistics for an operation (thread-safe)."""
        with self._lock:
            if operation not in self._metrics:
                return {}

            durations = list(self._metrics[operation])  # Copy for thread safety
            count = self._counts[operation]
            errors = self._errors[operation]

        if not durations:
            return {}

        return {
            "count": count,
            "errors": errors,
            "error_rate": errors / count if count > 0 else 0,
            "min_ms": min(durations) * 1000,
            "max_ms": max(durations) * 1000,
            "avg_ms": (sum(durations) / len(durations)) * 1000,
            "total_ms": sum(durations) * 1000,
        }

    def get_all_stats(self) -> Dict[str, Dict[str, Any]]:
        """Get statistics for all operations (thread-safe)."""
        with self._lock:
            operations = list(self._metrics.keys())
        return {op: self.get_stats(op) for op in operations}

    def reset(self) -> None:
        """Reset all metrics (thread-safe)."""
        with self._lock:
            self._metrics.clear()
            self._counts.clear()
            self._errors.clear()


class _MetricsContext:
    """Context manager for tracking metrics."""

    def __init__(self, metrics: PerformanceMetrics, operation: str):
        self.metrics = metrics
        self.operation = operation
        self.start_time: float = 0
        self.success = True

    def __enter__(self):
        self.start_time = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        duration = time.perf_counter() - self.start_time
        self.success = exc_type is None
        self.metrics.record(self.operation, duration, self.success)
        return False  # Don't suppress exceptions


# Global metrics instance
app_metrics = PerformanceMetrics()


# ==============================================================================
# Health Checks
# ==============================================================================

def check_rdkit_available() -> Dict[str, Any]:
    """Check if RDKit is available and working."""
    try:
        from rdkit import Chem
        from rdkit import __version__ as rdkit_version

        # Quick test
        mol = Chem.MolFromSmiles("CCO")
        if mol is None:
            return {"status": "error", "message": "RDKit failed to parse test molecule"}

        return {
            "status": "ok",
            "version": rdkit_version,
            "message": "RDKit is working"
        }
    except ImportError:
        return {"status": "error", "message": "RDKit not installed"}
    except Exception as e:
        return {"status": "error", "message": str(e)}


def check_dependencies() -> Dict[str, Dict[str, Any]]:
    """Check status of all required dependencies."""
    checks = {}

    # RDKit
    checks["rdkit"] = check_rdkit_available()

    # Pandas
    try:
        import pandas as pd
        checks["pandas"] = {"status": "ok", "version": pd.__version__}
    except ImportError:
        checks["pandas"] = {"status": "error", "message": "Not installed"}

    # NumPy
    try:
        import numpy as np
        checks["numpy"] = {"status": "ok", "version": np.__version__}
    except ImportError:
        checks["numpy"] = {"status": "error", "message": "Not installed"}

    # Plotly
    try:
        import plotly
        checks["plotly"] = {"status": "ok", "version": plotly.__version__}
    except ImportError:
        checks["plotly"] = {"status": "error", "message": "Not installed"}

    # Streamlit
    try:
        import streamlit as st
        checks["streamlit"] = {"status": "ok", "version": st.__version__}
    except ImportError:
        checks["streamlit"] = {"status": "error", "message": "Not installed"}

    # Kaleido (for chart export)
    try:
        import kaleido
        version = getattr(kaleido, '__version__', 'unknown')
        checks["kaleido"] = {"status": "ok", "version": version}
    except ImportError:
        checks["kaleido"] = {"status": "warning", "message": "Not installed (chart export disabled)"}
    except Exception:
        checks["kaleido"] = {"status": "warning", "message": "Version check failed"}

    # scikit-learn
    try:
        import sklearn
        checks["scikit-learn"] = {"status": "ok", "version": sklearn.__version__}
    except ImportError:
        checks["scikit-learn"] = {"status": "error", "message": "Not installed"}

    # statsmodels
    try:
        import statsmodels
        checks["statsmodels"] = {"status": "ok", "version": statsmodels.__version__}
    except ImportError:
        checks["statsmodels"] = {"status": "error", "message": "Not installed"}

    return checks


def get_system_info() -> Dict[str, Any]:
    """Get system information."""
    return {
        "python_version": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "machine": platform.machine(),
    }


def get_health_status() -> Dict[str, Any]:
    """
    Get overall application health status.

    Returns:
        Dictionary with health information
    """
    dependencies = check_dependencies()

    # Calculate overall status
    all_ok = all(dep["status"] == "ok" for dep in dependencies.values())
    has_errors = any(dep["status"] == "error" for dep in dependencies.values())

    if all_ok:
        overall_status = "healthy"
    elif has_errors:
        overall_status = "unhealthy"
    else:
        overall_status = "degraded"

    return {
        "status": overall_status,
        "timestamp": datetime.now().isoformat(),
        "system": get_system_info(),
        "dependencies": dependencies,
    }


# ==============================================================================
# Decorators
# ==============================================================================

def track_performance(operation_name: Optional[str] = None):
    """
    Decorator to track function performance.

    Args:
        operation_name: Name for the operation (default: function name)

    Example:
        @track_performance("molecule_calculation")
        def calculate_properties(smiles):
            ...
    """
    def decorator(func):
        name = operation_name or func.__name__

        @wraps(func)
        def wrapper(*args, **kwargs):
            with app_metrics.track(name):
                return func(*args, **kwargs)

        return wrapper
    return decorator


# ==============================================================================
# Logging Utilities
# ==============================================================================

_logging_configured = False
_logging_lock = Lock()


def setup_logging(
    level: int = logging.INFO,
    format_string: Optional[str] = None,
) -> None:
    """
    Setup application logging.

    Thread-safe and idempotent - safe to call multiple times.
    Delegates to the centralized config.logging_config.setup_logging.

    Args:
        level: Logging level (int or string like "INFO")
        format_string: Custom format string
    """
    global _logging_configured

    with _logging_lock:
        # Prevent duplicate handlers when called multiple times
        if _logging_configured:
            return

        # Convert int level to string for the config version
        level_name = logging.getLevelName(level) if isinstance(level, int) else level

        # Delegate to centralized logging config
        from molecular_calculator.config.logging_config import setup_logging as config_setup_logging
        config_setup_logging(level=level_name, log_format=format_string)

        _logging_configured = True


def log_operation(
    operation: str,
    success: bool = True,
    details: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Log an operation with structured data.

    Args:
        operation: Operation name
        success: Whether operation succeeded
        details: Additional details to log
    """
    log_data = {
        "operation": operation,
        "success": success,
        "timestamp": datetime.now().isoformat(),
    }

    if details:
        log_data.update(details)

    if success:
        logger.info(f"Operation completed: {log_data}")
    else:
        logger.error(f"Operation failed: {log_data}")


# ==============================================================================
# Streamlit Integration
# ==============================================================================

def render_health_dashboard():
    """Render a health status dashboard in Streamlit."""
    import streamlit as st

    health = get_health_status()

    # Overall status
    status = health["status"]
    if status == "healthy":
        st.success(f"System Status: {status.upper()}")
    elif status == "degraded":
        st.warning(f"System Status: {status.upper()}")
    else:
        st.error(f"System Status: {status.upper()}")

    # Dependencies
    st.subheader("Dependencies")
    deps = health["dependencies"]

    for dep_name, dep_info in deps.items():
        dep_status = dep_info["status"]
        version = dep_info.get("version", "N/A")
        message = dep_info.get("message", "")

        if dep_status == "ok":
            st.markdown(f"- **{dep_name}**: {version}")
        elif dep_status == "warning":
            st.markdown(f"- **{dep_name}**: {message}")
        else:
            st.markdown(f"- **{dep_name}**: {message}")

    # Performance metrics
    st.subheader("Performance Metrics")
    metrics = app_metrics.get_all_stats()

    if metrics:
        for op, stats in metrics.items():
            st.markdown(f"**{op}**: {stats['count']} calls, avg {stats['avg_ms']:.2f}ms")
    else:
        st.info("No performance data collected yet")
