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
from typing import Dict, Any, Optional, List
from datetime import datetime
from functools import wraps

logger = logging.getLogger(__name__)


# ==============================================================================
# Performance Metrics
# ==============================================================================

class PerformanceMetrics:
    """
    Track performance metrics for operations.

    Example:
        metrics = PerformanceMetrics()

        with metrics.track("calculation"):
            result = expensive_calculation()

        print(metrics.get_stats("calculation"))
    """

    def __init__(self):
        self._metrics: Dict[str, List[float]] = {}
        self._counts: Dict[str, int] = {}
        self._errors: Dict[str, int] = {}

    def record(self, operation: str, duration: float, success: bool = True) -> None:
        """Record a metric for an operation."""
        if operation not in self._metrics:
            self._metrics[operation] = []
            self._counts[operation] = 0
            self._errors[operation] = 0

        self._metrics[operation].append(duration)
        self._counts[operation] += 1

        if not success:
            self._errors[operation] += 1

        # Keep only last 1000 measurements
        if len(self._metrics[operation]) > 1000:
            self._metrics[operation] = self._metrics[operation][-1000:]

    def track(self, operation: str):
        """Context manager to track operation duration."""
        return _MetricsContext(self, operation)

    def get_stats(self, operation: str) -> Dict[str, Any]:
        """Get statistics for an operation."""
        if operation not in self._metrics:
            return {}

        durations = self._metrics[operation]
        if not durations:
            return {}

        return {
            "count": self._counts[operation],
            "errors": self._errors[operation],
            "error_rate": self._errors[operation] / self._counts[operation] if self._counts[operation] > 0 else 0,
            "min_ms": min(durations) * 1000,
            "max_ms": max(durations) * 1000,
            "avg_ms": (sum(durations) / len(durations)) * 1000,
            "total_ms": sum(durations) * 1000,
        }

    def get_all_stats(self) -> Dict[str, Dict[str, Any]]:
        """Get statistics for all operations."""
        return {op: self.get_stats(op) for op in self._metrics}

    def reset(self) -> None:
        """Reset all metrics."""
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
        checks["kaleido"] = {"status": "ok", "version": kaleido.__version__}
    except ImportError:
        checks["kaleido"] = {"status": "warning", "message": "Not installed (chart export disabled)"}

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

def setup_logging(
    level: int = logging.INFO,
    format_string: Optional[str] = None,
) -> None:
    """
    Setup application logging.

    Args:
        level: Logging level
        format_string: Custom format string
    """
    if format_string is None:
        format_string = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    logging.basicConfig(
        level=level,
        format=format_string,
        handlers=[logging.StreamHandler()]
    )

    # Reduce noise from external libraries
    logging.getLogger('urllib3').setLevel(logging.WARNING)
    logging.getLogger('requests').setLevel(logging.WARNING)


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
