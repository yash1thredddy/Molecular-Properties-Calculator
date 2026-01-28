"""Tests for utility modules.

Tests for:
- sanitizer.py
- cache.py
- rate_limiter.py
- monitoring.py
- error_handler.py
"""

import pytest
import time
from unittest.mock import Mock, patch

# ==============================================================================
# Sanitizer Tests
# ==============================================================================

class TestSanitizer:
    """Tests for sanitizer module."""

    def test_sanitize_smiles_valid(self):
        """Test sanitizing valid SMILES."""
        from molecular_calculator.utils.sanitizer import sanitize_smiles

        assert sanitize_smiles("CCO") == "CCO"
        assert sanitize_smiles("  CCO  ") == "CCO"  # Strips whitespace
        assert sanitize_smiles("c1ccccc1") == "c1ccccc1"  # Aromatic

    def test_sanitize_smiles_invalid(self):
        """Test sanitizing invalid SMILES."""
        from molecular_calculator.utils.sanitizer import sanitize_smiles

        assert sanitize_smiles("") is None
        assert sanitize_smiles(None) is None
        assert sanitize_smiles(123) is None

    def test_sanitize_smiles_removes_invalid_chars(self):
        """Test that invalid characters are removed."""
        from molecular_calculator.utils.sanitizer import sanitize_smiles

        # Invalid chars should be stripped
        result = sanitize_smiles("CCO<script>")
        assert "<" not in result
        assert ">" not in result

    def test_is_valid_smiles_format(self):
        """Test SMILES format validation."""
        from molecular_calculator.utils.sanitizer import is_valid_smiles_format

        assert is_valid_smiles_format("CCO") is True
        assert is_valid_smiles_format("c1ccccc1") is True
        assert is_valid_smiles_format("") is False
        assert is_valid_smiles_format("((()") is False  # Unbalanced

    def test_sanitize_inchi_key(self):
        """Test InChI Key sanitization."""
        from molecular_calculator.utils.sanitizer import sanitize_inchi_key

        # Valid InChI Key
        assert sanitize_inchi_key("LFQSCWFLJHTTHZ-UHFFFAOYSA-N") == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

        # Should convert to uppercase
        assert sanitize_inchi_key("lfqscwfljhtthz-uhfffaoysa-n") == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

        # Invalid
        assert sanitize_inchi_key("invalid") is None
        assert sanitize_inchi_key("") is None

    def test_sanitize_column_name(self):
        """Test column name sanitization."""
        from molecular_calculator.utils.sanitizer import sanitize_column_name

        assert sanitize_column_name("valid_name") == "valid_name"
        assert sanitize_column_name("Column Name") == "Column_Name"
        assert sanitize_column_name("") == "unnamed"
        assert sanitize_column_name("  ") == "unnamed"

    def test_sanitize_filename(self):
        """Test filename sanitization."""
        from molecular_calculator.utils.sanitizer import sanitize_filename

        assert sanitize_filename("test.csv") == "test.csv"
        # Path traversal attacks should be sanitized - keeps last part
        result = sanitize_filename("../../../etc/passwd")
        assert ".." not in result
        assert "/" not in result
        assert sanitize_filename("") == "file"

    def test_sanitize_html(self):
        """Test HTML sanitization."""
        from molecular_calculator.utils.sanitizer import sanitize_html

        assert "&lt;" in sanitize_html("<script>")
        assert "&gt;" in sanitize_html("</script>")
        assert sanitize_html("") == ""

    def test_sanitize_numeric(self):
        """Test numeric sanitization."""
        from molecular_calculator.utils.sanitizer import sanitize_numeric

        assert sanitize_numeric(42) == 42.0
        assert sanitize_numeric("42.5") == 42.5
        assert sanitize_numeric("invalid") == 0.0  # Default
        assert sanitize_numeric("invalid", default=99.0) == 99.0
        assert sanitize_numeric(150, max_val=100) == 100.0
        assert sanitize_numeric(-50, min_val=0) == 0.0

    def test_detect_identifier_type(self):
        """Test identifier type detection."""
        from molecular_calculator.utils.sanitizer import detect_identifier_type

        assert detect_identifier_type("CCO") == "smiles"
        assert detect_identifier_type("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3") == "inchi"
        assert detect_identifier_type("LFQSCWFLJHTTHZ-UHFFFAOYSA-N") == "inchi_key"
        assert detect_identifier_type("") == "unknown"


# ==============================================================================
# Cache Tests
# ==============================================================================

class TestCache:
    """Tests for cache module."""

    def test_ttl_cache_set_get(self):
        """Test basic cache set/get."""
        from molecular_calculator.utils.cache import TTLCache

        cache = TTLCache(maxsize=100, ttl=3600)
        cache.set("key1", "value1")

        assert cache.get("key1") == "value1"
        assert cache.get("nonexistent") is None

    def test_ttl_cache_expiration(self):
        """Test cache expiration using time mocking for reliability.

        This test uses mocking instead of time.sleep() to avoid flaky behavior
        due to timing variations across different systems and CI environments.
        """
        from molecular_calculator.utils.cache import TTLCache

        # Capture real time BEFORE any operations
        base_time = time.time()

        cache = TTLCache(maxsize=100, ttl=60)  # 60 second TTL
        cache.set("key1", "value1")

        # Verify value is accessible
        assert cache.get("key1") == "value1"

        # Mock time.time() to simulate expiration
        # The cache stores entry time and checks if (current_time - entry_time) > ttl
        with patch('time.time') as mock_time:
            # Set time to 61 seconds after cache was set (beyond TTL)
            mock_time.return_value = base_time + 61
            assert cache.get("key1") is None, "Cache entry should expire after TTL"

    def test_ttl_cache_expiration_boundary(self):
        """Test cache expiration at boundary conditions."""
        from molecular_calculator.utils.cache import TTLCache

        # Capture real time BEFORE any operations
        base_time = time.time()

        cache = TTLCache(maxsize=100, ttl=10)  # 10 second TTL
        cache.set("key1", "value1")

        # Before expiration
        with patch('time.time') as mock_time:
            mock_time.return_value = base_time + 9  # Just before TTL
            assert cache.get("key1") == "value1", "Cache should still be valid before TTL"

        # Right at expiration
        with patch('time.time') as mock_time:
            mock_time.return_value = base_time + 11  # Just after TTL
            assert cache.get("key1") is None, "Cache should expire after TTL"

    def test_ttl_cache_maxsize(self):
        """Test cache eviction at maxsize."""
        from molecular_calculator.utils.cache import TTLCache

        cache = TTLCache(maxsize=2, ttl=3600)
        cache.set("key1", "value1")
        cache.set("key2", "value2")
        cache.set("key3", "value3")  # Should evict key1

        assert cache.get("key1") is None  # Evicted
        assert cache.get("key2") == "value2"
        assert cache.get("key3") == "value3"

    def test_ttl_cache_stats(self):
        """Test cache statistics."""
        from molecular_calculator.utils.cache import TTLCache

        cache = TTLCache(maxsize=100, ttl=3600)
        cache.set("key1", "value1")
        cache.get("key1")  # Hit
        cache.get("key1")  # Hit
        cache.get("nonexistent")  # Miss

        stats = cache.stats()
        assert stats["hits"] == 2
        assert stats["misses"] == 1
        assert stats["hit_rate"] == 2/3

    def test_ttl_cache_delete(self):
        """Test cache deletion."""
        from molecular_calculator.utils.cache import TTLCache

        cache = TTLCache(maxsize=100, ttl=3600)
        cache.set("key1", "value1")

        assert cache.delete("key1") is True
        assert cache.get("key1") is None
        assert cache.delete("nonexistent") is False

    def test_ttl_cache_clear(self):
        """Test cache clearing."""
        from molecular_calculator.utils.cache import TTLCache

        cache = TTLCache(maxsize=100, ttl=3600)
        cache.set("key1", "value1")
        cache.set("key2", "value2")
        cache.clear()

        assert cache.get("key1") is None
        assert cache.get("key2") is None
        assert len(cache) == 0


# ==============================================================================
# Rate Limiter Tests
# ==============================================================================

class TestRateLimiter:
    """Tests for rate limiter module."""

    def test_rate_limiter_acquire(self):
        """Test basic rate limiter acquisition."""
        from molecular_calculator.utils.rate_limiter import RateLimiter
        from molecular_calculator.utils.exceptions import RateLimitError

        limiter = RateLimiter(max_requests=5, window_seconds=60, name="test")

        # Should allow 5 requests
        for _ in range(5):
            assert limiter.acquire() is True

        # 6th should raise RateLimitError when blocking=False
        try:
            limiter.acquire(blocking=False)
            assert False, "Should have raised RateLimitError"
        except RateLimitError:
            pass  # Expected

    def test_rate_limiter_window_reset(self):
        """Test that rate limiter resets after window."""
        from molecular_calculator.utils.rate_limiter import RateLimiter

        limiter = RateLimiter(max_requests=2, window_seconds=0.1, name="test")

        assert limiter.acquire() is True
        assert limiter.acquire() is True
        # After 2 requests, should not be allowed
        assert limiter.is_allowed() is False

        time.sleep(0.15)  # Wait for window reset
        assert limiter.is_allowed() is True  # Should work again

    def test_rate_limiter_remaining(self):
        """Test remaining count."""
        from molecular_calculator.utils.rate_limiter import RateLimiter

        limiter = RateLimiter(max_requests=5, window_seconds=60, name="test")

        assert limiter.remaining_requests() == 5
        limiter.acquire()
        assert limiter.remaining_requests() == 4

    def test_rate_limiter_repr(self):
        """Test limiter string representation."""
        from molecular_calculator.utils.rate_limiter import RateLimiter

        limiter = RateLimiter(max_requests=5, window_seconds=60, name="test")
        limiter.acquire()
        limiter.acquire()

        repr_str = repr(limiter)
        assert "test" in repr_str
        assert "max_requests=5" in repr_str


# ==============================================================================
# Monitoring Tests
# ==============================================================================

class TestMonitoring:
    """Tests for monitoring module."""

    def test_performance_metrics_record(self):
        """Test recording performance metrics."""
        from molecular_calculator.utils.monitoring import PerformanceMetrics

        metrics = PerformanceMetrics()
        metrics.record("test_op", 0.5, success=True)
        metrics.record("test_op", 0.3, success=True)
        metrics.record("test_op", 0.4, success=False)

        stats = metrics.get_stats("test_op")
        assert stats["count"] == 3
        assert stats["errors"] == 1
        assert stats["error_rate"] == 1/3

    def test_performance_metrics_track_context(self):
        """Test context manager for tracking."""
        from molecular_calculator.utils.monitoring import PerformanceMetrics

        metrics = PerformanceMetrics()

        with metrics.track("test_op"):
            time.sleep(0.01)  # Simulate work

        stats = metrics.get_stats("test_op")
        assert stats["count"] == 1
        assert stats["avg_ms"] >= 10  # At least 10ms

    def test_check_rdkit_available(self):
        """Test RDKit availability check."""
        from molecular_calculator.utils.monitoring import check_rdkit_available

        result = check_rdkit_available()
        assert result["status"] in ["ok", "error"]
        assert "message" in result or "version" in result

    def test_get_system_info(self):
        """Test system info retrieval."""
        from molecular_calculator.utils.monitoring import get_system_info

        info = get_system_info()
        assert "python_version" in info
        assert "platform" in info

    def test_get_health_status(self):
        """Test health status check."""
        from molecular_calculator.utils.monitoring import get_health_status

        health = get_health_status()
        assert health["status"] in ["healthy", "degraded", "unhealthy"]
        assert "timestamp" in health
        assert "dependencies" in health


# ==============================================================================
# Error Handler Tests
# ==============================================================================

class TestErrorHandler:
    """Tests for error handler module."""

    def test_error_boundary_success(self):
        """Test error boundary with successful operation."""
        from molecular_calculator.utils.error_handler import error_boundary

        with error_boundary("test"):
            result = 1 + 1

        # Should complete without error
        assert result == 2

    def test_error_boundary_exception(self):
        """Test error boundary catching exception."""
        from molecular_calculator.utils.error_handler import error_boundary

        with error_boundary("test", show_error=False):
            raise ValueError("Test error")

        # Should not propagate exception

    def test_suppress_warnings(self):
        """Test warning suppression."""
        from molecular_calculator.utils.error_handler import suppress_warnings
        import warnings

        with suppress_warnings():
            warnings.warn("This should be suppressed")

        # Should complete without showing warning

    def test_get_error_details(self):
        """Test error details extraction."""
        from molecular_calculator.utils.error_handler import get_error_details

        try:
            raise ValueError("Test error message")
        except ValueError as e:
            details = get_error_details(e)

        assert details["type"] == "ValueError"
        assert "Test error message" in details["message"]


# ==============================================================================
# Suggestions Tests
# ==============================================================================

class TestSuggestions:
    """Tests for suggestions module."""

    def test_detect_smiles_column(self):
        """Test SMILES column detection."""
        import pandas as pd
        from molecular_calculator.utils.suggestions import detect_smiles_column

        df = pd.DataFrame({
            "SMILES": ["CCO", "CCC", "CCCC"],
            "Name": ["Ethanol", "Propane", "Butane"]
        })

        assert detect_smiles_column(df) == "SMILES"

    def test_detect_id_column(self):
        """Test ID column detection."""
        import pandas as pd
        from molecular_calculator.utils.suggestions import detect_id_column

        df = pd.DataFrame({
            "ID": [1, 2, 3],
            "Compound_Name": ["A", "B", "C"],
            "Value": [1.0, 2.0, 3.0]
        })

        result = detect_id_column(df)
        assert result in ["ID", "Compound_Name"]

    def test_detect_column_type(self):
        """Test column type detection by name patterns."""
        from molecular_calculator.utils.suggestions import detect_column_type

        # detect_column_type only takes column name, not DataFrame
        # Patterns are checked in order: smiles -> id -> activity -> property
        # Test values must only match their intended pattern category

        # "SMILES" matches "smiles" in SMILES_COLUMN_PATTERNS
        assert detect_column_type("SMILES") == "smiles"
        # "chembl_id" matches "chembl_id" in ID_COLUMN_PATTERNS (no SMILES overlap)
        assert detect_column_type("chembl_id") == "id"
        # "IC50" matches "ic50" in ACTIVITY_COLUMN_PATTERNS
        assert detect_column_type("IC50") == "activity"
        # "tpsa" matches "tpsa" in PROPERTY_COLUMN_PATTERNS
        # (avoid "logp" since "log" matches activity patterns first)
        assert detect_column_type("tpsa") == "property"
        # "unknown_col" matches no patterns
        assert detect_column_type("unknown_col") == "unknown"

    def test_get_column_stats(self):
        """Test column statistics."""
        import pandas as pd
        from molecular_calculator.utils.suggestions import get_column_stats

        df = pd.DataFrame({
            "values": [1.0, 2.0, 3.0, 4.0, 5.0]
        })

        stats = get_column_stats(df, "values")
        assert stats["mean"] == 3.0
        assert stats["min"] == 1.0
        assert stats["max"] == 5.0


# ==============================================================================
# Export Tests
# ==============================================================================

class TestExport:
    """Tests for export module."""

    def test_export_formats_defined(self):
        """Test that export formats are defined."""
        from molecular_calculator.utils.export import EXPORT_FORMATS

        assert "png" in EXPORT_FORMATS
        assert "svg" in EXPORT_FORMATS
        assert "html" in EXPORT_FORMATS

    def test_get_download_data_html(self):
        """Test HTML export data generation."""
        import plotly.graph_objects as go
        from molecular_calculator.utils.export import get_download_data

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=[1, 2, 3], y=[1, 2, 3]))

        # get_download_data returns (bytes, mime_type, file_extension)
        data, mime, ext = get_download_data(fig, "html")

        assert data is not None
        assert mime == "text/html"
        assert ext == ".html"
