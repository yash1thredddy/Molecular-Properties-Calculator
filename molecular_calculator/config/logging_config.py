"""Logging configuration for the application.

This module provides centralized logging setup with consistent
formatting and configurable log levels.
"""

import logging
import sys
from typing import Optional
from pathlib import Path


def setup_logging(
    level: str = "INFO",
    log_file: Optional[str] = None,
    log_format: Optional[str] = None
) -> logging.Logger:
    """Configure application logging.

    Sets up logging with console output and optional file output.
    Configures consistent formatting and suppresses noisy third-party loggers.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional file path for log output
        log_format: Optional custom format string

    Returns:
        Root logger instance

    Example:
        >>> from molecular_calculator.config import setup_logging
        >>> setup_logging(level="DEBUG")
        >>> import logging
        >>> logger = logging.getLogger(__name__)
        >>> logger.info("Application started")
    """
    # Default format
    if log_format is None:
        log_format = (
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

    # Create formatter
    formatter = logging.Formatter(
        log_format,
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Get root logger
    root_logger = logging.getLogger()

    # Clear existing handlers
    root_logger.handlers.clear()

    # Set level
    log_level = getattr(logging, level.upper(), logging.INFO)
    root_logger.setLevel(log_level)

    # Create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    console_handler.setLevel(log_level)
    root_logger.addHandler(console_handler)

    # Add file handler if specified
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setFormatter(formatter)
        file_handler.setLevel(log_level)
        root_logger.addHandler(file_handler)

    # Suppress noisy third-party loggers
    _configure_third_party_loggers()

    return root_logger


def _configure_third_party_loggers() -> None:
    """Configure third-party library loggers to reduce noise."""
    noisy_loggers = [
        'urllib3',
        'matplotlib',
        'PIL',
        'streamlit',
        'watchdog',
        'fsevents',
        'asyncio',
        'rdkit',
    ]

    for logger_name in noisy_loggers:
        logging.getLogger(logger_name).setLevel(logging.WARNING)


def get_logger(name: str) -> logging.Logger:
    """Get a logger with the given name.

    This is a convenience function that ensures the logger
    is properly configured.

    Args:
        name: Logger name (typically __name__)

    Returns:
        Configured logger instance

    Example:
        >>> from molecular_calculator.config.logging_config import get_logger
        >>> logger = get_logger(__name__)
        >>> logger.info("Processing started")
    """
    return logging.getLogger(name)


class LoggerAdapter(logging.LoggerAdapter):
    """Custom logger adapter for adding context to log messages.

    Example:
        >>> logger = get_logger(__name__)
        >>> adapted = LoggerAdapter(logger, {'user': 'john'})
        >>> adapted.info("Processing file")
        # Output: 2024-01-01 12:00:00 - module - INFO - [user=john] Processing file
    """

    def process(self, msg, kwargs):
        """Add context to log message."""
        context = ' '.join(f'[{k}={v}]' for k, v in self.extra.items())
        return f'{context} {msg}' if context else msg, kwargs
