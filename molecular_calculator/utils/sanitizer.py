# ==============================================================================
# Input Sanitization Utilities
# ==============================================================================
"""
Input sanitization utilities for molecular data.

This module focuses on CLEANING and TRANSFORMING data to be safe.
For VALIDATION (checking if data is valid), use validators.py instead.

The key difference:
- sanitize_smiles() -> cleans/transforms, returns cleaned string or None
- InputValidator.validate_smiles() -> validates, returns ValidationResult

Use sanitizer for cleaning user input before processing.
Use validators for checking if data meets requirements.
"""

import re
import logging
from typing import Optional, List, Any

# Import shared patterns from validators to avoid duplication
from .validators import InputValidator

logger = logging.getLogger(__name__)


# ==============================================================================
# SMILES Sanitization
# ==============================================================================

# Valid SMILES characters (atoms, bonds, rings, branches, stereochemistry)
VALID_SMILES_CHARS = set(
    'BCNOPSFIHcnopsbrCl' +  # Atoms
    '0123456789' +          # Ring numbers
    '()[]' +                # Branches and brackets
    '=#-+\\/@.' +           # Bonds and stereochemistry
    '%'                     # Extended ring numbers
)

# Maximum reasonable SMILES length
MAX_SMILES_LENGTH = 5000

# Re-use patterns from validators
INCHI_KEY_PATTERN = InputValidator.INCHI_KEY_PATTERN
INCHI_PATTERN = InputValidator.INCHI_PATTERN


def sanitize_smiles(smiles: str) -> Optional[str]:
    """
    Sanitize a SMILES string by removing invalid characters.

    Args:
        smiles: Raw SMILES string

    Returns:
        Sanitized SMILES or None if invalid
    """
    if not smiles or not isinstance(smiles, str):
        return None

    # Strip whitespace
    smiles = smiles.strip()

    if not smiles:
        return None

    # Check length
    if len(smiles) > MAX_SMILES_LENGTH:
        logger.warning(f"SMILES too long: {len(smiles)} chars")
        return None

    # Remove any non-printable characters
    smiles = ''.join(c for c in smiles if c.isprintable())

    # Validate characters
    invalid_chars = set(smiles) - VALID_SMILES_CHARS
    if invalid_chars:
        logger.warning(f"Invalid SMILES characters removed: {invalid_chars}")
        smiles = ''.join(c for c in smiles if c in VALID_SMILES_CHARS)

    return smiles if smiles else None


def is_valid_smiles_format(smiles: str) -> bool:
    """
    Quick format check for SMILES (not chemical validity).

    Args:
        smiles: SMILES string to check

    Returns:
        True if format appears valid
    """
    if not smiles or not isinstance(smiles, str):
        return False

    smiles = smiles.strip()

    # Basic checks
    if len(smiles) == 0 or len(smiles) > MAX_SMILES_LENGTH:
        return False

    # Must contain at least one atom
    if not re.search(r'[BCNOPSFIHcnopsbrCl]', smiles):
        return False

    # Balanced parentheses
    if smiles.count('(') != smiles.count(')'):
        return False

    # Balanced brackets
    if smiles.count('[') != smiles.count(']'):
        return False

    return True


# ==============================================================================
# InChI Sanitization
# ==============================================================================
# Note: INCHI_PATTERN and INCHI_KEY_PATTERN are imported from validators above


def sanitize_inchi(inchi: str) -> Optional[str]:
    """
    Sanitize an InChI string.

    Args:
        inchi: Raw InChI string

    Returns:
        Sanitized InChI or None if invalid
    """
    if not inchi or not isinstance(inchi, str):
        return None

    inchi = inchi.strip()

    # Must start with InChI=
    if not inchi.startswith('InChI='):
        return None

    # Basic format validation
    if INCHI_PATTERN.match(inchi):
        return inchi

    return None


def sanitize_inchi_key(inchi_key: str) -> Optional[str]:
    """
    Sanitize an InChI Key.

    Args:
        inchi_key: Raw InChI Key string

    Returns:
        Sanitized InChI Key or None if invalid
    """
    if not inchi_key or not isinstance(inchi_key, str):
        return None

    # Convert to uppercase and strip
    inchi_key = inchi_key.strip().upper()

    # Validate format: XXXXXXXXXXXXXX-XXXXXXXXXX-X
    if INCHI_KEY_PATTERN.match(inchi_key):
        return inchi_key

    return None


def is_valid_inchi_key(inchi_key: str) -> bool:
    """
    Check if string is a valid InChI Key format.

    Args:
        inchi_key: String to check

    Returns:
        True if valid InChI Key format
    """
    return sanitize_inchi_key(inchi_key) is not None


# ==============================================================================
# General Text Sanitization
# ==============================================================================

def sanitize_column_name(name: str) -> str:
    """
    Sanitize a column name for safe use.

    Args:
        name: Raw column name

    Returns:
        Sanitized column name
    """
    if not name or not isinstance(name, str):
        return "unnamed"

    # Strip and limit length
    name = name.strip()[:100]

    # Replace problematic characters
    name = re.sub(r'[^\w\s\-_]', '', name)

    # Replace spaces with underscores
    name = name.replace(' ', '_')

    # Remove consecutive underscores
    name = re.sub(r'_+', '_', name)

    # Strip leading/trailing underscores
    name = name.strip('_')

    return name or "unnamed"


def sanitize_filename(filename: str) -> str:
    """
    Sanitize a filename for safe file operations.

    Args:
        filename: Raw filename

    Returns:
        Sanitized filename
    """
    if not filename or not isinstance(filename, str):
        return "file"

    # Get just the filename, not path
    filename = filename.replace('\\', '/').split('/')[-1]

    # Remove null bytes and other dangerous characters
    filename = re.sub(r'[\x00-\x1f\x7f<>:"/\\|?*]', '', filename)

    # Limit length
    filename = filename[:200]

    # Don't allow hidden files or parent directory references
    filename = filename.lstrip('.')

    return filename or "file"


def sanitize_html(text: str) -> str:
    """
    Escape HTML special characters.

    Args:
        text: Raw text

    Returns:
        HTML-safe text
    """
    if not text or not isinstance(text, str):
        return ""

    replacements = {
        '&': '&amp;',
        '<': '&lt;',
        '>': '&gt;',
        '"': '&quot;',
        "'": '&#x27;',
    }

    for char, replacement in replacements.items():
        text = text.replace(char, replacement)

    return text


# ==============================================================================
# Numeric Sanitization
# ==============================================================================

def sanitize_numeric(
    value: Any,
    min_val: Optional[float] = None,
    max_val: Optional[float] = None,
    default: float = 0.0
) -> float:
    """
    Sanitize a numeric value.

    Args:
        value: Raw value
        min_val: Minimum allowed value
        max_val: Maximum allowed value
        default: Default if conversion fails

    Returns:
        Sanitized float value
    """
    try:
        result = float(value)

        # Check for infinity/NaN
        if not (-float('inf') < result < float('inf')):
            return default

        # Apply bounds
        if min_val is not None:
            result = max(min_val, result)
        if max_val is not None:
            result = min(max_val, result)

        return result

    except (TypeError, ValueError):
        return default


def sanitize_integer(
    value: Any,
    min_val: Optional[int] = None,
    max_val: Optional[int] = None,
    default: int = 0
) -> int:
    """
    Sanitize an integer value.

    Args:
        value: Raw value
        min_val: Minimum allowed value
        max_val: Maximum allowed value
        default: Default if conversion fails

    Returns:
        Sanitized integer value
    """
    try:
        result = int(float(value))

        if min_val is not None:
            result = max(min_val, result)
        if max_val is not None:
            result = min(max_val, result)

        return result

    except (TypeError, ValueError):
        return default


# ==============================================================================
# Batch Sanitization
# ==============================================================================

def sanitize_smiles_list(smiles_list: List[str]) -> List[Optional[str]]:
    """
    Sanitize a list of SMILES strings.

    Args:
        smiles_list: List of raw SMILES

    Returns:
        List of sanitized SMILES (None for invalid)
    """
    return [sanitize_smiles(s) for s in smiles_list]


def detect_identifier_type(identifier: str) -> str:
    """
    Detect the type of molecular identifier.

    This is an alias for InputValidator.detect_format() for backwards compatibility.

    Args:
        identifier: Molecular identifier string

    Returns:
        Type: 'smiles', 'inchi', 'inchi_key', or 'unknown'
    """
    return InputValidator.detect_format(identifier)
