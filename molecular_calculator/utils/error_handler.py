# ==============================================================================
# Centralized Error Handling
# ==============================================================================
"""
Error handling utilities including decorators and context managers.
Provides consistent error handling patterns across the application.
"""

import functools
import logging
import traceback
from contextlib import contextmanager
from typing import Callable, Optional, Type, Tuple, Any, TypeVar

import streamlit as st

from .exceptions import (
    MolecularCalculatorError,
    InvalidSMILESError,
    PropertyCalculationError,
    FileValidationError,
    APIError,
    RateLimitError,
    ConversionError,
)

logger = logging.getLogger(__name__)

# Type variable for generic return types
T = TypeVar('T')


# ==============================================================================
# Error Decorators
# ==============================================================================

def handle_errors(
    default_return: Any = None,
    show_user_error: bool = True,
    log_level: str = "error",
    reraise: bool = False,
    catch: Tuple[Type[Exception], ...] = (Exception,),
) -> Callable:
    """
    Decorator for consistent error handling.

    Args:
        default_return: Value to return on error (if not reraising)
        show_user_error: Whether to display error in Streamlit UI
        log_level: Logging level ('debug', 'info', 'warning', 'error')
        reraise: Whether to re-raise the exception after handling
        catch: Tuple of exception types to catch

    Returns:
        Decorator function

    Example:
        @handle_errors(default_return=None, show_user_error=True)
        def risky_operation():
            ...
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except catch as e:
                # Log the error
                log_func = getattr(logger, log_level, logger.error)
                log_func(f"Error in {func.__name__}: {e}", exc_info=True)

                # Show user-friendly error if requested
                if show_user_error:
                    _show_user_error(e)

                # Re-raise or return default
                if reraise:
                    raise
                return default_return

        return wrapper
    return decorator


def handle_calculation_errors(func: Callable) -> Callable:
    """
    Decorator specifically for molecular calculation functions.

    Catches calculation-specific exceptions and provides appropriate feedback.

    Example:
        @handle_calculation_errors
        def calculate_properties(smiles):
            ...
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except InvalidSMILESError as e:
            logger.warning(f"Invalid SMILES in {func.__name__}: {e}")
            st.error(f"Invalid molecular structure: {e}")
            return None
        except PropertyCalculationError as e:
            logger.error(f"Calculation error in {func.__name__}: {e}")
            st.error(f"Could not calculate properties: {e}")
            return None
        except MolecularCalculatorError as e:
            logger.error(f"Calculator error in {func.__name__}: {e}")
            st.error(f"Calculation failed: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error in {func.__name__}: {e}", exc_info=True)
            st.error(f"An unexpected error occurred: {str(e)}")
            return None

    return wrapper


def handle_api_errors(
    max_retries: int = 3,
    default_return: Any = None,
) -> Callable:
    """
    Decorator for API call error handling with retry logic.

    Args:
        max_retries: Maximum number of retry attempts
        default_return: Value to return if all retries fail

    Example:
        @handle_api_errors(max_retries=3)
        def call_external_api():
            ...
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            last_error = None

            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except RateLimitError as e:
                    logger.warning(f"Rate limit hit in {func.__name__}, attempt {attempt + 1}")
                    last_error = e
                    if attempt < max_retries - 1:
                        import time
                        time.sleep(2 ** attempt)  # Exponential backoff
                except APIError as e:
                    logger.error(f"API error in {func.__name__}: {e}")
                    last_error = e
                    break  # Don't retry on general API errors
                except Exception as e:
                    logger.error(f"Unexpected API error in {func.__name__}: {e}")
                    last_error = e
                    break

            # All retries exhausted
            if last_error:
                logger.error(f"All retries failed for {func.__name__}: {last_error}")
            return default_return

        return wrapper
    return decorator


def handle_file_errors(func: Callable) -> Callable:
    """
    Decorator for file operation error handling.

    Example:
        @handle_file_errors
        def read_uploaded_file(file):
            ...
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except FileValidationError as e:
            logger.warning(f"File validation error in {func.__name__}: {e}")
            st.error(f"Invalid file: {e}")
            return None
        except FileNotFoundError as e:
            logger.error(f"File not found in {func.__name__}: {e}")
            st.error("The requested file could not be found.")
            return None
        except PermissionError as e:
            logger.error(f"Permission error in {func.__name__}: {e}")
            st.error("Permission denied when accessing the file.")
            return None
        except Exception as e:
            logger.error(f"File error in {func.__name__}: {e}", exc_info=True)
            st.error(f"Error processing file: {str(e)}")
            return None

    return wrapper


# ==============================================================================
# Context Managers
# ==============================================================================

@contextmanager
def error_boundary(
    operation_name: str = "operation",
    show_error: bool = True,
    default_return: Any = None,
):
    """
    Context manager for error boundaries.

    Args:
        operation_name: Name of the operation for error messages
        show_error: Whether to show error in Streamlit UI
        default_return: Value to yield on error

    Example:
        with error_boundary("data processing"):
            result = process_data(df)
    """
    try:
        yield
    except MolecularCalculatorError as e:
        logger.error(f"Error during {operation_name}: {e}")
        if show_error:
            st.error(f"Error during {operation_name}: {e}")
    except Exception as e:
        logger.error(f"Unexpected error during {operation_name}: {e}", exc_info=True)
        if show_error:
            st.error(f"An error occurred during {operation_name}. Please try again.")


@contextmanager
def suppress_warnings():
    """
    Context manager to suppress RDKit and numpy warnings temporarily.

    Example:
        with suppress_warnings():
            mol = Chem.MolFromSmiles(smiles)
    """
    import warnings
    import numpy as np

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        old_settings = np.seterr(all='ignore')
        try:
            yield
        finally:
            np.seterr(**old_settings)


# ==============================================================================
# Helper Functions
# ==============================================================================

def _show_user_error(error: Exception) -> None:
    """
    Display a user-friendly error message in Streamlit.

    Args:
        error: The exception to display
    """
    error_messages = {
        InvalidSMILESError: "The molecular structure could not be parsed. Please check the input format.",
        PropertyCalculationError: "Could not calculate one or more molecular properties.",
        FileValidationError: "The uploaded file is invalid or in an unsupported format.",
        APIError: "Could not connect to the external service. Please try again later.",
        RateLimitError: "Too many requests. Please wait a moment and try again.",
        ConversionError: "Could not convert the molecular identifier.",
    }

    error_type = type(error)
    message = error_messages.get(error_type, str(error))

    st.error(message)


def get_error_details(error: Exception) -> dict:
    """
    Get detailed error information for logging/debugging.

    Args:
        error: The exception to analyze

    Returns:
        Dictionary with error details
    """
    return {
        "type": type(error).__name__,
        "message": str(error),
        "traceback": traceback.format_exc(),
        "args": error.args,
    }


def log_and_raise(
    error: Exception,
    message: Optional[str] = None,
    log_level: str = "error",
) -> None:
    """
    Log an error and re-raise it.

    Args:
        error: The exception to log and raise
        message: Optional additional message
        log_level: Logging level

    Raises:
        The original exception
    """
    log_func = getattr(logger, log_level, logger.error)
    full_message = f"{message}: {error}" if message else str(error)
    log_func(full_message, exc_info=True)
    raise error
