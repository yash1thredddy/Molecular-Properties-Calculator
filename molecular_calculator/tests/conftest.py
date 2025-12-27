"""Pytest configuration and fixtures.

This module provides shared fixtures for all tests in the package.
Fixtures are automatically discovered by pytest.
"""

import pytest
import pandas as pd
from io import BytesIO
from typing import List
from unittest.mock import MagicMock


# ============================================================================
# SMILES Fixtures
# ============================================================================

@pytest.fixture
def sample_smiles() -> List[str]:
    """Sample valid SMILES strings for testing.

    Returns a list of common molecules that should always parse correctly.
    """
    return [
        'C',              # Methane
        'CC',             # Ethane
        'CCO',            # Ethanol
        'c1ccccc1',       # Benzene
        'CC(=O)O',        # Acetic acid
        'CCN(CC)CC',      # Triethylamine
        'CC(C)C',         # Isobutane
        'c1ccc2ccccc2c1', # Naphthalene
    ]


@pytest.fixture
def sample_smiles_with_names() -> List[tuple]:
    """Sample SMILES with their common names."""
    return [
        ('C', 'Methane'),
        ('CC', 'Ethane'),
        ('CCO', 'Ethanol'),
        ('c1ccccc1', 'Benzene'),
        ('CC(=O)O', 'Acetic acid'),
        ('CC(=O)Oc1ccccc1C(=O)O', 'Aspirin'),
    ]


@pytest.fixture
def invalid_smiles() -> List[str]:
    """Invalid SMILES strings for testing error handling."""
    return [
        '',              # Empty string
        '   ',           # Whitespace only
        'invalid',       # Not valid SMILES syntax
        'X1234',         # Invalid element
        '(((',           # Unbalanced parentheses
        'C1CC',          # Unclosed ring
    ]


@pytest.fixture
def edge_case_smiles() -> List[str]:
    """Edge case SMILES for testing robustness."""
    return [
        '[H][H]',        # Hydrogen molecule
        '[Na+]',         # Sodium ion
        '[O-]',          # Oxide ion
        'C#N',           # Hydrogen cyanide
        'C=C',           # Ethylene
        '[C@@H]',        # Chiral center
    ]


# ============================================================================
# InChI Fixtures
# ============================================================================

@pytest.fixture
def sample_inchi() -> str:
    """Sample valid InChI string (ethanol)."""
    return 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'


@pytest.fixture
def sample_inchi_key() -> str:
    """Sample valid InChI Key (ethanol)."""
    return 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'


@pytest.fixture
def invalid_inchi_keys() -> List[str]:
    """Invalid InChI Keys for testing.

    Note: Lowercase keys are NOT invalid - validator normalizes to uppercase.
    See test_validate_inchi_key_lowercase for lowercase handling.
    """
    return [
        '',
        'INVALID-KEY',
        'LFQSCWFLJHTTHZ-UHFFFAOYSA',      # Missing last character (26 chars)
        'LFQSCWFLJHTTHZ-UHFFFAOYSA-NN',   # Extra character (28 chars)
        'LFQSCWFLJHTTHZ-UHFFFAOYSA-1',    # Invalid last character (number)
    ]


# ============================================================================
# DataFrame Fixtures
# ============================================================================

@pytest.fixture
def sample_dataframe(sample_smiles) -> pd.DataFrame:
    """Sample DataFrame with SMILES and properties."""
    return pd.DataFrame({
        'SMILES': sample_smiles[:5],
        'Name': ['Methane', 'Ethane', 'Ethanol', 'Benzene', 'Acetic acid'],
        'Category': ['Alkane', 'Alkane', 'Alcohol', 'Aromatic', 'Acid'],
        'MW': [16.04, 30.07, 46.07, 78.11, 60.05],
        'LogP': [-0.64, 1.02, -0.14, 1.99, -0.17],
    })


@pytest.fixture
def large_dataframe() -> pd.DataFrame:
    """Large DataFrame for performance testing."""
    n = 15000  # Above warning threshold
    return pd.DataFrame({
        'SMILES': ['CCO'] * n,
        'Index': range(n),
    })


@pytest.fixture
def empty_dataframe() -> pd.DataFrame:
    """Empty DataFrame for testing edge cases."""
    return pd.DataFrame()


@pytest.fixture
def dataframe_with_missing() -> pd.DataFrame:
    """DataFrame with missing values."""
    return pd.DataFrame({
        'SMILES': ['CCO', None, 'c1ccccc1', '', 'CC'],
        'Name': ['Ethanol', 'Unknown', None, 'Empty', 'Ethane'],
        'MW': [46.07, None, 78.11, None, 30.07],
    })


# ============================================================================
# File Fixtures
# ============================================================================

@pytest.fixture
def mock_csv_file() -> MagicMock:
    """Mock uploaded CSV file."""
    content = b'SMILES,Name\nCCO,Ethanol\nCC,Ethane\n'
    file = MagicMock()
    file.name = 'test.csv'
    file.size = len(content)
    file.getvalue.return_value = content
    file.read.return_value = content
    return file


@pytest.fixture
def mock_xlsx_file() -> MagicMock:
    """Mock uploaded XLSX file."""
    # Create minimal XLSX content
    df = pd.DataFrame({'SMILES': ['CCO'], 'Name': ['Ethanol']})
    buffer = BytesIO()
    df.to_excel(buffer, index=False, engine='openpyxl')
    content = buffer.getvalue()

    file = MagicMock()
    file.name = 'test.xlsx'
    file.size = len(content)
    file.getvalue.return_value = content
    return file


@pytest.fixture
def mock_large_file() -> MagicMock:
    """Mock file that exceeds size limit."""
    file = MagicMock()
    file.name = 'large.csv'
    file.size = 100 * 1024 * 1024  # 100 MB
    return file


@pytest.fixture
def mock_invalid_file() -> MagicMock:
    """Mock file with invalid extension."""
    file = MagicMock()
    file.name = 'data.txt'
    file.size = 1000
    return file


# ============================================================================
# Security Test Fixtures
# ============================================================================

@pytest.fixture
def xss_payloads() -> List[str]:
    """XSS attack payloads for security testing."""
    return [
        '<script>alert("xss")</script>',
        '"><script>alert(1)</script>',
        "javascript:alert('xss')",
        '<img src=x onerror=alert(1)>',
        '<svg onload=alert(1)>',
    ]


@pytest.fixture
def sql_injection_payloads() -> List[str]:
    """SQL injection payloads for security testing."""
    return [
        "'; DROP TABLE users; --",
        "1' OR '1'='1",
        "1; DELETE FROM molecules",
        "' UNION SELECT * FROM secrets --",
    ]


@pytest.fixture
def path_traversal_payloads() -> List[str]:
    """Path traversal payloads for security testing."""
    return [
        '../../../etc/passwd',
        '..\\..\\..\\windows\\system32',
        '/etc/passwd',
        'C:\\Windows\\System32',
    ]


# ============================================================================
# Calculator Fixtures (for Phase 2)
# ============================================================================

@pytest.fixture
def property_names() -> List[str]:
    """List of common property names."""
    return [
        'Molecular_Weight',
        'LogP',
        'TPSA',
        'HBA',
        'HBD',
        'RotatableBonds',
        'RingCount',
    ]


# ============================================================================
# Helper Functions
# ============================================================================

def create_mock_file(
    filename: str,
    content: bytes,
    size: int = None
) -> MagicMock:
    """Create a mock file object.

    Args:
        filename: Name of the file
        content: File content as bytes
        size: Override size (defaults to len(content))

    Returns:
        Mock file object
    """
    file = MagicMock()
    file.name = filename
    file.size = size if size is not None else len(content)
    file.getvalue.return_value = content
    file.read.return_value = content
    return file
