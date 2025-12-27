# ==============================================================================
# Rate Limiter Utilities
# ==============================================================================
"""
Rate limiting utilities for API calls.
Prevents exceeding external API rate limits.
"""

import time
import logging
from collections import deque
from threading import Lock
from typing import Optional
from functools import wraps

from .exceptions import RateLimitError

logger = logging.getLogger(__name__)


class RateLimiter:
    """
    Token bucket rate limiter for API calls.

    Implements a sliding window rate limiter that tracks requests
    within a time window and enforces rate limits.

    Attributes:
        max_requests: Maximum requests allowed in the window
        window_seconds: Time window in seconds
        name: Optional name for logging

    Example:
        limiter = RateLimiter(max_requests=30, window_seconds=60)

        if limiter.is_allowed():
            make_api_call()
        else:
            wait_time = limiter.time_until_allowed()
            time.sleep(wait_time)
    """

    def __init__(
        self,
        max_requests: int = 30,
        window_seconds: float = 60.0,
        name: Optional[str] = None,
    ):
        """
        Initialize the rate limiter.

        Args:
            max_requests: Maximum requests allowed in window
            window_seconds: Size of sliding window in seconds
            name: Optional name for identification
        """
        self.max_requests = max_requests
        self.window_seconds = window_seconds
        self.name = name or "RateLimiter"
        self._timestamps: deque = deque()
        self._lock = Lock()

    def _cleanup_old_requests(self) -> None:
        """Remove timestamps outside the current window."""
        current_time = time.time()
        cutoff = current_time - self.window_seconds

        while self._timestamps and self._timestamps[0] < cutoff:
            self._timestamps.popleft()

    def is_allowed(self) -> bool:
        """
        Check if a request is allowed under the rate limit.

        Returns:
            True if request is allowed, False otherwise
        """
        with self._lock:
            self._cleanup_old_requests()
            return len(self._timestamps) < self.max_requests

    def record_request(self) -> None:
        """Record a request timestamp."""
        with self._lock:
            self._timestamps.append(time.time())

    def acquire(self, blocking: bool = True, timeout: Optional[float] = None) -> bool:
        """
        Acquire permission to make a request.

        Args:
            blocking: If True, wait until allowed; if False, return immediately
            timeout: Maximum time to wait in seconds (None = no limit)

        Returns:
            True if acquired, False if timed out or not blocking

        Raises:
            RateLimitError: If blocking is False and rate limit exceeded
        """
        start_time = time.time()

        while True:
            if self.is_allowed():
                self.record_request()
                return True

            if not blocking:
                raise RateLimitError(
                    f"Rate limit exceeded for {self.name}. "
                    f"Max {self.max_requests} requests per {self.window_seconds}s"
                )

            # Check timeout
            if timeout is not None:
                elapsed = time.time() - start_time
                if elapsed >= timeout:
                    return False

            # Wait a bit before retrying
            wait_time = min(0.1, self.time_until_allowed())
            time.sleep(wait_time)

    def time_until_allowed(self) -> float:
        """
        Get time in seconds until next request is allowed.

        Returns:
            Time in seconds (0 if allowed now)
        """
        with self._lock:
            self._cleanup_old_requests()

            if len(self._timestamps) < self.max_requests:
                return 0.0

            # Time until oldest request expires
            oldest = self._timestamps[0]
            wait_time = (oldest + self.window_seconds) - time.time()
            return max(0.0, wait_time)

    def remaining_requests(self) -> int:
        """
        Get number of requests remaining in current window.

        Returns:
            Number of available requests
        """
        with self._lock:
            self._cleanup_old_requests()
            return max(0, self.max_requests - len(self._timestamps))

    def reset(self) -> None:
        """Clear all recorded timestamps."""
        with self._lock:
            self._timestamps.clear()

    def __repr__(self) -> str:
        return (
            f"RateLimiter(name={self.name!r}, "
            f"max_requests={self.max_requests}, "
            f"window_seconds={self.window_seconds}, "
            f"remaining={self.remaining_requests()})"
        )


# ==============================================================================
# Pre-configured Rate Limiters
# ==============================================================================

# NIH CIR API: ~30 requests per minute recommended
nih_limiter = RateLimiter(
    max_requests=30,
    window_seconds=60.0,
    name="NIH_CIR"
)

# PubChem API: ~5 requests per second recommended
pubchem_limiter = RateLimiter(
    max_requests=5,
    window_seconds=1.0,
    name="PubChem"
)


# ==============================================================================
# Decorator
# ==============================================================================

def rate_limited(limiter: RateLimiter, blocking: bool = True):
    """
    Decorator to apply rate limiting to a function.

    Args:
        limiter: RateLimiter instance to use
        blocking: Whether to block until allowed

    Example:
        @rate_limited(nih_limiter)
        def call_nih_api(inchi_key):
            ...
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            limiter.acquire(blocking=blocking)
            return func(*args, **kwargs)
        return wrapper
    return decorator


# ==============================================================================
# Utility Functions
# ==============================================================================

def get_limiter_status() -> dict:
    """
    Get status of all pre-configured rate limiters.

    Returns:
        Dictionary with limiter status information
    """
    return {
        "nih_cir": {
            "remaining": nih_limiter.remaining_requests(),
            "max": nih_limiter.max_requests,
            "window_seconds": nih_limiter.window_seconds,
        },
        "pubchem": {
            "remaining": pubchem_limiter.remaining_requests(),
            "max": pubchem_limiter.max_requests,
            "window_seconds": pubchem_limiter.window_seconds,
        },
    }


def reset_all_limiters() -> None:
    """Reset all pre-configured rate limiters."""
    nih_limiter.reset()
    pubchem_limiter.reset()
    logger.info("All rate limiters reset")
