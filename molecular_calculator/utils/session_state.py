"""Session state management utilities.

This module provides centralized session state management for Streamlit,
with features like:
- Default value initialization
- Type-safe access
- File change detection
- Mode-specific state clearing
"""

import hashlib
import logging
from typing import Any, Optional, Dict, List

logger = logging.getLogger(__name__)


class SessionState:
    """Centralized session state management.

    This class provides a clean interface for managing Streamlit session state.
    It handles initialization, access, and cleanup of session state values.

    Note: This class uses static methods because Streamlit session state
    is global and doesn't benefit from instance-based access.

    Example:
        >>> from molecular_calculator.utils import SessionState
        >>> SessionState.init_defaults()
        >>> SessionState.set('my_key', 'my_value')
        >>> value = SessionState.get('my_key')
    """

    # Default values for session state keys
    DEFAULTS: Dict[str, Any] = {
        # Data storage
        'batch_results_df': None,
        'current_smiles_col': None,
        'uploaded_file_hash': None,
        'selected_properties': [],

        # UI state
        'mode': 'single',
        'chart_type': 'Scatter Plot',
        'show_structure': True,

        # Processing state
        'is_processing': False,
        'last_error': None,
    }

    # Keys associated with each mode
    MODE_KEYS: Dict[str, List[str]] = {
        'batch': [
            'batch_results_df',
            'current_smiles_col',
            'uploaded_file_hash',
            'selected_properties',
        ],
        'single': [
            'single_molecule_result',
            'single_smiles',
        ],
        'visualization': [
            'viz_data',
            'viz_config',
        ],
    }

    @classmethod
    def _get_session_state(cls):
        """Get Streamlit session state (lazy import for testing)."""
        try:
            import streamlit as st
            return st.session_state
        except ImportError:
            # Fallback for testing without Streamlit
            if not hasattr(cls, '_mock_state'):
                cls._mock_state = {}
            return cls._mock_state

    @classmethod
    def init_defaults(cls) -> None:
        """Initialize default session state values.

        Call this at the start of your Streamlit app to ensure
        all expected keys exist with sensible defaults.

        Example:
            >>> # At the top of your Streamlit app
            >>> SessionState.init_defaults()
        """
        session_state = cls._get_session_state()

        for key, default in cls.DEFAULTS.items():
            if key not in session_state:
                session_state[key] = default
                logger.debug(f"Initialized session state key: {key}")

    @classmethod
    def get(cls, key: str, default: Any = None) -> Any:
        """Get a value from session state.

        Args:
            key: Session state key
            default: Default value if key not found

        Returns:
            Value from session state or default

        Example:
            >>> value = SessionState.get('my_key', 'default_value')
        """
        session_state = cls._get_session_state()
        return session_state.get(key, default)

    @classmethod
    def set(cls, key: str, value: Any) -> None:
        """Set a value in session state.

        Args:
            key: Session state key
            value: Value to set

        Example:
            >>> SessionState.set('my_key', 'my_value')
        """
        session_state = cls._get_session_state()
        session_state[key] = value
        logger.debug(f"Set session state: {key} = {type(value).__name__}")

    @classmethod
    def clear(cls, key: str) -> None:
        """Clear a session state key.

        Removes the key from session state entirely.

        Args:
            key: Session state key to clear

        Example:
            >>> SessionState.clear('temporary_data')
        """
        session_state = cls._get_session_state()
        if key in session_state:
            del session_state[key]
            logger.debug(f"Cleared session state key: {key}")

    @classmethod
    def clear_mode(cls, mode: str) -> None:
        """Clear all state associated with a specific mode.

        This is useful when switching between modes to ensure
        stale data doesn't persist.

        Args:
            mode: Mode to clear ('batch', 'single', 'visualization')

        Example:
            >>> SessionState.clear_mode('batch')
        """
        keys = cls.MODE_KEYS.get(mode, [])
        for key in keys:
            cls.clear(key)
        logger.debug(f"Cleared session state for mode: {mode}")

    @classmethod
    def clear_all(cls) -> None:
        """Clear all session state.

        Use with caution - this removes all stored data.
        """
        session_state = cls._get_session_state()
        keys = list(session_state.keys())
        for key in keys:
            try:
                del session_state[key]
            except KeyError:
                pass
        logger.debug("Cleared all session state")

    @classmethod
    def file_changed(cls, file) -> bool:
        """Check if an uploaded file has changed.

        Uses MD5 hash to detect when a file with the same name
        has different content.

        Args:
            file: Uploaded file object (must have getvalue() method)

        Returns:
            True if file is new or changed, False otherwise

        Example:
            >>> if SessionState.file_changed(uploaded_file):
            ...     # Reprocess the file
            ...     process_file(uploaded_file)
        """
        if file is None:
            return False

        try:
            # Calculate file hash
            file_content = file.getvalue()
            current_hash = hashlib.md5(file_content).hexdigest()

            # Compare with stored hash
            stored_hash = cls.get('uploaded_file_hash')

            if current_hash != stored_hash:
                cls.set('uploaded_file_hash', current_hash)
                logger.debug("File change detected")
                return True

            return False

        except Exception as e:
            logger.warning(f"Error checking file change: {e}")
            return True  # Assume changed if we can't check

    @classmethod
    def has(cls, key: str) -> bool:
        """Check if a key exists in session state.

        Args:
            key: Session state key

        Returns:
            True if key exists
        """
        session_state = cls._get_session_state()
        return key in session_state

    @classmethod
    def get_or_set(cls, key: str, default: Any) -> Any:
        """Get value if exists, otherwise set and return default.

        This is useful for lazy initialization.

        Args:
            key: Session state key
            default: Default value to set if key doesn't exist

        Returns:
            Existing value or newly set default

        Example:
            >>> # Initialize on first access
            >>> cache = SessionState.get_or_set('cache', {})
        """
        if not cls.has(key):
            cls.set(key, default)
        return cls.get(key)

    @classmethod
    def update(cls, key: str, updates: Dict[str, Any]) -> None:
        """Update a dictionary value in session state.

        Args:
            key: Session state key (must contain a dict)
            updates: Dictionary of updates to apply

        Example:
            >>> SessionState.set('config', {'a': 1})
            >>> SessionState.update('config', {'b': 2})
            >>> SessionState.get('config')
            {'a': 1, 'b': 2}
        """
        current = cls.get(key, {})
        if isinstance(current, dict):
            current.update(updates)
            cls.set(key, current)
        else:
            logger.warning(f"Cannot update non-dict value at key: {key}")
