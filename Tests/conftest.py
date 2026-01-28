"""Pytest configuration and shared fixtures for Tests/ directory.

This file provides shared fixtures that can be reused across all test modules.
Fixtures are automatically discovered by pytest.

The fixtures are organized into categories:
- SMILES fixtures: Valid and invalid SMILES strings
- InChI fixtures: InChI and InChI Key test data
- DataFrame fixtures: Sample DataFrames for batch processing tests
- File fixtures: Mock file objects for upload testing
- Security fixtures: Payloads for security testing
- Molecule fixtures: Well-characterized molecules with known properties
"""

import sys
from pathlib import Path
from typing import List, Dict, Any, Tuple
from io import BytesIO
from unittest.mock import MagicMock

import pytest
import pandas as pd

# Add project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


# =============================================================================
# SMILES Fixtures
# =============================================================================

@pytest.fixture
def sample_smiles() -> List[str]:
    """Sample valid SMILES strings for testing.

    Returns a list of common molecules that should always parse correctly.
    Molecules selected for structural diversity.
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
def sample_smiles_with_names() -> List[Tuple[str, str]]:
    """Sample SMILES with their common names.

    Useful for tests that need molecule identification.
    """
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
        'ZZZZZ',         # Not valid elements
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
        '[CH3]',         # Methyl radical (if supported)
    ]


@pytest.fixture
def stereochemistry_smiles() -> List[Tuple[str, str]]:
    """SMILES with stereochemistry for testing chiral handling."""
    return [
        ('C[C@H](O)CC', 'R-2-butanol'),
        ('C[C@@H](O)CC', 'S-2-butanol'),
        ('C/C=C/C', 'E-2-butene'),
        ('C/C=C\\C', 'Z-2-butene'),
        ('C[C@H](O)[C@@H](C)O', 'Threitol-like'),
    ]


# =============================================================================
# InChI Fixtures
# =============================================================================

@pytest.fixture
def sample_inchi() -> str:
    """Sample valid InChI string (ethanol).

    Source: PubChem CID 702
    """
    return 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'


@pytest.fixture
def sample_inchi_key() -> str:
    """Sample valid InChI Key (ethanol).

    Source: PubChem CID 702
    """
    return 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'


@pytest.fixture
def invalid_inchi_keys() -> List[str]:
    """Invalid InChI Keys for testing error handling.

    Note: Lowercase keys are NOT invalid - validator normalizes to uppercase.
    """
    return [
        '',
        'INVALID-KEY',
        'LFQSCWFLJHTTHZ-UHFFFAOYSA',      # Missing last character (26 chars)
        'LFQSCWFLJHTTHZ-UHFFFAOYSA-NN',   # Extra character (28 chars)
        'LFQSCWFLJHTTHZ-UHFFFAOYSA-1',    # Invalid last character (number)
        'ABC',                              # Too short
    ]


@pytest.fixture
def molecule_identifiers() -> Dict[str, Dict[str, str]]:
    """Complete identifier sets for test molecules.

    Each molecule has SMILES, InChI, and InChI Key for cross-conversion testing.
    Sources: PubChem
    """
    return {
        'ethanol': {
            'smiles': 'CCO',
            'inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
            'inchi_key': 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N',
        },
        'aspirin': {
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
            'inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
        },
        'caffeine': {
            'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'inchi': 'InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3',
            'inchi_key': 'RYYVLZVUVIJVGH-UHFFFAOYSA-N',
        },
    }


# =============================================================================
# DataFrame Fixtures
# =============================================================================

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
def dataframe_with_activity() -> pd.DataFrame:
    """DataFrame with activity data for LEI testing."""
    return pd.DataFrame({
        'SMILES': ['CCO', 'c1ccccc1', 'CC(=O)O', 'CCN'],
        'Name': ['Ethanol', 'Benzene', 'Acetic acid', 'Ethylamine'],
        'pKi': [5.0, 6.5, 4.2, 5.8],
        'MW': [46.07, 78.11, 60.05, 45.08],
        'TPSA': [20.23, 0.0, 37.3, 26.02],
    })


@pytest.fixture
def large_dataframe() -> pd.DataFrame:
    """Large DataFrame for performance testing."""
    n = 15000  # Above typical warning threshold
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
    """DataFrame with missing values for testing NaN handling."""
    return pd.DataFrame({
        'SMILES': ['CCO', None, 'c1ccccc1', '', 'CC'],
        'Name': ['Ethanol', 'Unknown', None, 'Empty', 'Ethane'],
        'MW': [46.07, None, 78.11, None, 30.07],
    })


@pytest.fixture
def dataframe_various_column_names() -> pd.DataFrame:
    """DataFrame with various column naming conventions."""
    return pd.DataFrame({
        'canonical_smiles': ['CCO', 'CC'],
        'compound_id': ['ETH001', 'ETH002'],
        'Molecular_Weight': [46.07, 30.07],
        'pIC50': [5.5, 4.2],
    })


# =============================================================================
# File Fixtures
# =============================================================================

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


# =============================================================================
# Security Test Fixtures
# =============================================================================

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


# =============================================================================
# Well-Characterized Molecule Fixtures
# =============================================================================

@pytest.fixture
def aspirin_data() -> Dict[str, Any]:
    """Well-characterized aspirin data.

    Reference: PubChem CID 2244, DrugBank DB00945
    """
    return {
        'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
        'inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
        'name': 'Aspirin',
        'expected_mw': 180.158,
        'expected_logp_range': (1.0, 2.0),
        'expected_hbd': 1,
        'expected_hba': 3,
        'pubchem_cid': 2244,
    }


@pytest.fixture
def caffeine_data() -> Dict[str, Any]:
    """Well-characterized caffeine data.

    Reference: PubChem CID 2519, DrugBank DB00201
    """
    return {
        'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
        'inchi': 'InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3',
        'inchi_key': 'RYYVLZVUVIJVGH-UHFFFAOYSA-N',
        'name': 'Caffeine',
        'expected_mw': 194.191,
        'expected_logp_range': (-2.0, 0.5),
        'expected_hbd': 0,
        'expected_hba': 6,
        'pubchem_cid': 2519,
    }


@pytest.fixture
def ethanol_data() -> Dict[str, Any]:
    """Well-characterized ethanol data.

    Reference: PubChem CID 702
    """
    return {
        'smiles': 'CCO',
        'inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
        'inchi_key': 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N',
        'name': 'Ethanol',
        'expected_mw': 46.069,
        'expected_logp_range': (-0.5, 0.0),
        'expected_hbd': 1,
        'expected_hba': 1,
        'pubchem_cid': 702,
    }


@pytest.fixture
def ibuprofen_data() -> Dict[str, Any]:
    """Well-characterized ibuprofen data (clean drug, no interference).

    Reference: PubChem CID 3672
    Useful for testing clean molecules that shouldn't trigger interference flags.
    """
    return {
        'smiles': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
        'inchi_key': 'HEFNNWSXXWATRW-UHFFFAOYSA-N',
        'name': 'Ibuprofen',
        'expected_mw': 206.285,
        'expected_logp_range': (3.0, 4.5),
        'pubchem_cid': 3672,
        'is_clean_drug': True,  # No PAINS/interference flags expected
    }


# =============================================================================
# Property Groups Fixture
# =============================================================================

@pytest.fixture
def property_names() -> List[str]:
    """List of common property names calculated by MolecularCalculator."""
    return [
        'Molecular_Weight',
        'LogP',
        'TPSA',
        'HB_Donors',
        'HB_Acceptors',
        'Rotatable_Bonds',
        'Aromatic_Rings',
        'Heavy_Atom_Count',
        'QED',
        'Lipinski_Violations',
        'Veber_Violations',
    ]


# =============================================================================
# Helper Functions
# =============================================================================

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


def create_test_dataframe(
    smiles_list: List[str],
    names: List[str] = None,
    activities: List[float] = None
) -> pd.DataFrame:
    """Create a test DataFrame with optional activity data.

    Args:
        smiles_list: List of SMILES strings
        names: Optional list of molecule names
        activities: Optional list of pKi values

    Returns:
        DataFrame suitable for testing
    """
    data = {'SMILES': smiles_list}
    if names:
        data['Name'] = names
    if activities:
        data['pKi'] = activities
    return pd.DataFrame(data)
