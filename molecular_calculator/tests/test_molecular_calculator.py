"""Tests for the refactored MolecularCalculator.

This module tests the core molecular calculator functionality including:
- Property calculation
- Format conversion
- Batch processing
- Singleton pattern
- Parallel processing
- Input sanitization
"""

import pytest
import pandas as pd

from molecular_calculator.core import MolecularCalculator, get_calculator
from molecular_calculator.models import (
    InputFormat,
    MolecularProperties,
    CalculationResult,
    ConversionResult,
    PROPERTY_GROUPS,
)
from molecular_calculator.services import (
    ConversionService,
    PropertyCalculator,
)
from molecular_calculator.utils.sanitizer import (
    sanitize_smiles,
    sanitize_inchi,
    sanitize_inchi_key,
    sanitize_column_name,
    sanitize_filename,
    sanitize_html,
)


class TestMolecularCalculator:
    """Tests for MolecularCalculator class."""

    @pytest.fixture
    def calculator(self):
        """Create a MolecularCalculator instance."""
        return MolecularCalculator()

    # ========================================================================
    # Property Calculation Tests
    # ========================================================================

    def test_calculate_molecular_properties_ethanol(self, calculator):
        """Test property calculation for ethanol."""
        props = calculator.calculate_molecular_properties("CCO")

        assert 'Molecular_Weight' in props
        assert abs(props['Molecular_Weight'] - 46.069) < 0.01
        assert props['Heavy_Atom_Count'] == 3
        assert props['HB_Donors'] == 1
        assert props['HB_Acceptors'] == 1

    def test_calculate_molecular_properties_benzene(self, calculator):
        """Test property calculation for benzene."""
        props = calculator.calculate_molecular_properties("c1ccccc1")

        assert abs(props['Molecular_Weight'] - 78.114) < 0.01
        assert props['Aromatic_Rings'] == 1
        assert props['Heavy_Atom_Count'] == 6

    def test_calculate_molecular_properties_aspirin(self, calculator):
        """Test property calculation for aspirin."""
        props = calculator.calculate_molecular_properties("CC(=O)Oc1ccccc1C(=O)O")

        assert 'Molecular_Weight' in props
        assert abs(props['Molecular_Weight'] - 180.159) < 0.1
        assert props['Aromatic_Rings'] == 1

    def test_calculate_molecular_properties_empty(self, calculator):
        """Test property calculation with empty SMILES."""
        props = calculator.calculate_molecular_properties("")
        assert props == {}

    def test_calculate_molecular_properties_none(self, calculator):
        """Test property calculation with None."""
        props = calculator.calculate_molecular_properties(None)
        assert props == {}

    def test_calculate_molecular_properties_invalid(self, calculator):
        """Test property calculation with invalid SMILES."""
        props = calculator.calculate_molecular_properties("invalid_smiles")
        assert props == {}

    def test_calculate_returns_result_object(self, calculator):
        """Test that calculate method returns CalculationResult."""
        result = calculator.calculate("CCO")

        assert isinstance(result, CalculationResult)
        assert result.success
        assert result.smiles == "CCO"
        assert result.properties is not None
        assert isinstance(result.properties, MolecularProperties)

    def test_calculate_invalid_returns_error(self, calculator):
        """Test that calculate returns error for invalid SMILES."""
        result = calculator.calculate("not_valid")

        assert isinstance(result, CalculationResult)
        assert not result.success
        assert result.error is not None

    # ========================================================================
    # Lipinski Rule Tests
    # ========================================================================

    def test_lipinski_compliant_molecule(self, calculator):
        """Test Lipinski rule compliance for a drug-like molecule."""
        props = calculator.calculate_molecular_properties("CCO")  # Ethanol

        assert props['Lipinski_Violations'] == 0

    def test_lipinski_violation_detection(self, calculator):
        """Test Lipinski violation detection."""
        # Large molecule that violates MW > 500
        large_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        props = calculator.calculate_molecular_properties(large_smiles)

        # Should have at least one violation
        assert 'Lipinski_Violations' in props

    # ========================================================================
    # Format Detection Tests
    # ========================================================================

    def test_detect_input_format_smiles(self):
        """Test format detection for SMILES."""
        assert MolecularCalculator.detect_input_format("CCO") == "smiles"
        assert MolecularCalculator.detect_input_format("c1ccccc1") == "smiles"

    def test_detect_input_format_inchi(self):
        """Test format detection for InChI."""
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        assert MolecularCalculator.detect_input_format(inchi) == "inchi"

    def test_detect_input_format_inchi_key(self):
        """Test format detection for InChI Key."""
        inchi_key = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        assert MolecularCalculator.detect_input_format(inchi_key) == "inchi_key"

    # ========================================================================
    # Format Conversion Tests
    # ========================================================================

    def test_convert_smiles_passthrough(self, calculator):
        """Test that SMILES passes through unchanged."""
        result = MolecularCalculator.convert_to_smiles("CCO", "smiles")
        assert result == "CCO"

    def test_convert_inchi_to_smiles(self, calculator):
        """Test InChI to SMILES conversion."""
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = MolecularCalculator.convert_to_smiles(inchi, "inchi")

        assert result is not None
        # Verify it's valid SMILES by calculating properties
        props = calculator.calculate_molecular_properties(result)
        assert 'Molecular_Weight' in props

    def test_convert_instance_method(self, calculator):
        """Test convert instance method."""
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = calculator.convert(inchi, InputFormat.INCHI)

        assert isinstance(result, ConversionResult)
        assert result.success
        assert result.smiles is not None

    def test_calculate_from_any_format(self, calculator):
        """Test calculating properties from InChI directly."""
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = calculator.calculate_from_any_format(inchi)

        assert result.success
        assert result.properties is not None
        assert abs(result.properties.molecular_weight - 46.069) < 0.01

    # ========================================================================
    # Column Detection Tests
    # ========================================================================

    def test_detect_smiles_column_uppercase(self):
        """Test SMILES column detection with uppercase name."""
        df = pd.DataFrame({'SMILES': ['CCO'], 'Name': ['Ethanol']})
        assert MolecularCalculator.detect_smiles_column(df) == 'SMILES'

    def test_detect_smiles_column_lowercase(self):
        """Test SMILES column detection with lowercase name."""
        df = pd.DataFrame({'smiles': ['CCO'], 'name': ['Ethanol']})
        assert MolecularCalculator.detect_smiles_column(df) == 'smiles'

    def test_detect_smiles_column_mixed_case(self):
        """Test SMILES column detection with mixed case."""
        df = pd.DataFrame({'Smiles': ['CCO'], 'Name': ['Ethanol']})
        assert MolecularCalculator.detect_smiles_column(df) == 'Smiles'

    def test_detect_smiles_column_smi(self):
        """Test SMILES column detection with 'smi' name."""
        df = pd.DataFrame({'smi': ['CCO'], 'name': ['Ethanol']})
        assert MolecularCalculator.detect_smiles_column(df) == 'smi'

    def test_detect_smiles_column_not_found(self):
        """Test SMILES column detection when not present."""
        df = pd.DataFrame({'Molecule': ['CCO'], 'Name': ['Ethanol']})
        assert MolecularCalculator.detect_smiles_column(df) is None

    # ========================================================================
    # Property Groups Tests
    # ========================================================================

    def test_get_property_groups_default(self):
        """Test getting default property groups."""
        groups = MolecularCalculator.get_property_groups()

        assert 'Basic Properties' in groups
        assert 'Lipinski Properties' in groups
        assert 'Drug-likeness' in groups
        assert 'Molecular_Weight' in groups['Basic Properties']

    def test_get_property_groups_with_lei(self):
        """Test getting property groups with LEI included."""
        groups = MolecularCalculator.get_property_groups(include_lei=True)

        assert 'Ligand Efficiency Indices' in groups
        assert 'NSEI' in groups['Ligand Efficiency Indices']

    # ========================================================================
    # Batch Processing Tests
    # ========================================================================

    def test_process_batch_basic(self):
        """Test basic batch processing."""
        df = pd.DataFrame({
            'SMILES': ['C', 'CC', 'CCO'],
            'Name': ['Methane', 'Ethane', 'Ethanol']
        })

        result = MolecularCalculator.process_batch(df, 'SMILES')

        assert 'Molecular_Weight' in result.columns
        assert len(result) == 3
        # Check ethanol MW
        assert abs(result.loc[2, 'Molecular_Weight'] - 46.069) < 0.01

    def test_process_batch_with_selection(self):
        """Test batch processing with selected properties."""
        df = pd.DataFrame({
            'SMILES': ['CCO'],
            'Name': ['Ethanol']
        })

        result = MolecularCalculator.process_batch(
            df,
            'SMILES',
            selected_properties={'Molecular_Weight', 'LogP'}
        )

        assert 'Molecular_Weight' in result.columns
        assert 'LogP' in result.columns
        # Should not have unselected properties
        assert 'TPSA' not in result.columns or pd.isna(result.loc[0, 'TPSA'])

    def test_process_batch_with_invalid_smiles(self):
        """Test batch processing handles invalid SMILES."""
        df = pd.DataFrame({
            'SMILES': ['CCO', 'invalid', 'c1ccccc1'],
            'Name': ['Ethanol', 'Invalid', 'Benzene']
        })

        result = MolecularCalculator.process_batch(df, 'SMILES')

        assert len(result) == 3
        # Valid molecules should have MW
        assert pd.notna(result.loc[0, 'Molecular_Weight'])
        assert pd.notna(result.loc[2, 'Molecular_Weight'])

    def test_process_batch_with_nan(self):
        """Test batch processing handles NaN values."""
        df = pd.DataFrame({
            'SMILES': ['CCO', None, 'c1ccccc1'],
            'Name': ['Ethanol', 'Missing', 'Benzene']
        })

        result = MolecularCalculator.process_batch(df, 'SMILES')

        assert len(result) == 3
        # Valid molecules should have MW
        assert pd.notna(result.loc[0, 'Molecular_Weight'])
        assert pd.notna(result.loc[2, 'Molecular_Weight'])


class TestPropertyCalculator:
    """Tests for PropertyCalculator service."""

    @pytest.fixture
    def property_calc(self):
        """Create a PropertyCalculator instance."""
        return PropertyCalculator()

    def test_calculate_returns_result(self, property_calc):
        """Test calculate returns CalculationResult."""
        result = property_calc.calculate("CCO")

        assert isinstance(result, CalculationResult)
        assert result.success
        assert result.properties is not None

    def test_calculate_as_dict(self, property_calc):
        """Test calculate_as_dict returns dictionary."""
        result = property_calc.calculate_as_dict("CCO")

        assert isinstance(result, dict)
        assert 'Molecular_Weight' in result

    def test_calculate_qed(self, property_calc):
        """Test QED calculation."""
        result = property_calc.calculate("CCO")

        assert result.properties.qed is not None
        assert 0 <= result.properties.qed <= 1

    def test_calculate_all_properties_present(self, property_calc):
        """Test that all expected properties are calculated."""
        result = property_calc.calculate("c1ccccc1")
        props = result.properties

        # Check basic properties
        assert props.molecular_weight is not None
        assert props.heavy_atom_count is not None
        assert props.atom_count is not None

        # Check Lipinski properties
        assert props.logp is not None
        assert props.hb_donors is not None
        assert props.hb_acceptors is not None
        assert props.tpsa is not None

        # Check ring properties
        assert props.aromatic_rings is not None

    def test_get_all_property_names(self, property_calc):
        """Test getting all property names."""
        names = PropertyCalculator.get_all_property_names()

        assert isinstance(names, list)
        assert 'Molecular_Weight' in names
        assert 'LogP' in names


class TestConversionService:
    """Tests for ConversionService."""

    @pytest.fixture
    def service(self):
        """Create a ConversionService instance."""
        return ConversionService()

    def test_detect_format_smiles(self, service):
        """Test detecting SMILES format."""
        assert service.detect_format("CCO") == InputFormat.SMILES
        assert service.detect_format("c1ccccc1") == InputFormat.SMILES

    def test_detect_format_inchi(self, service):
        """Test detecting InChI format."""
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        assert service.detect_format(inchi) == InputFormat.INCHI

    def test_detect_format_inchi_key(self, service):
        """Test detecting InChI Key format."""
        key = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        assert service.detect_format(key) == InputFormat.INCHI_KEY

    def test_detect_format_unknown(self, service):
        """Test detecting unknown format."""
        assert service.detect_format("") == InputFormat.UNKNOWN
        assert service.detect_format(None) == InputFormat.UNKNOWN

    def test_to_smiles_from_smiles(self, service):
        """Test SMILES passthrough."""
        result = service.to_smiles("CCO")

        assert result.success
        assert result.smiles == "CCO"
        assert result.source_format == InputFormat.SMILES

    def test_to_smiles_from_inchi(self, service):
        """Test InChI to SMILES conversion."""
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = service.to_smiles(inchi)

        assert result.success
        assert result.smiles is not None
        assert result.source_format == InputFormat.INCHI

    def test_to_smiles_invalid_smiles(self, service):
        """Test handling invalid SMILES."""
        result = service.to_smiles("not_valid_smiles", InputFormat.SMILES)

        assert not result.success
        assert result.error is not None

    def test_validate_smiles_valid(self, service):
        """Test SMILES validation with valid input."""
        assert service.validate_smiles("CCO")
        assert service.validate_smiles("c1ccccc1")

    def test_validate_smiles_invalid(self, service):
        """Test SMILES validation with invalid input."""
        assert not service.validate_smiles("invalid")
        assert not service.validate_smiles("")
        assert not service.validate_smiles(None)


class TestMolecularProperties:
    """Tests for MolecularProperties dataclass."""

    def test_to_dict(self):
        """Test converting properties to dictionary."""
        props = MolecularProperties(
            molecular_weight=46.069,
            heavy_atom_count=3,
            logp=-0.14
        )

        d = props.to_dict()

        assert d['Molecular_Weight'] == 46.069
        assert d['Heavy_Atom_Count'] == 3
        assert d['LogP'] == -0.14

    def test_to_dict_filters_none(self):
        """Test that None values are filtered from dict."""
        props = MolecularProperties(
            molecular_weight=46.069,
            heavy_atom_count=None
        )

        d = props.to_dict()

        assert 'Molecular_Weight' in d
        assert 'Heavy_Atom_Count' not in d

    def test_from_dict(self):
        """Test creating properties from dictionary."""
        d = {
            'Molecular_Weight': 46.069,
            'Heavy_Atom_Count': 3,
            'LogP': -0.14
        }

        props = MolecularProperties.from_dict(d)

        assert props.molecular_weight == 46.069
        assert props.heavy_atom_count == 3
        assert props.logp == -0.14


# ============================================================================
# Singleton Pattern Tests
# ============================================================================

class TestGetCalculator:
    """Tests for get_calculator() singleton function."""

    def test_get_calculator_returns_instance(self):
        """Test that get_calculator returns a MolecularCalculator instance."""
        calc = get_calculator()
        assert isinstance(calc, MolecularCalculator)

    def test_get_calculator_returns_same_instance(self):
        """Test that get_calculator returns the same instance each time."""
        calc1 = get_calculator()
        calc2 = get_calculator()
        assert calc1 is calc2

    def test_get_calculator_is_functional(self):
        """Test that the singleton can calculate properties."""
        calc = get_calculator()
        result = calc.calculate("CCO")
        assert result.success
        assert result.properties.molecular_weight is not None


# ============================================================================
# Parallel Processing Tests
# ============================================================================

class TestParallelBatchProcessing:
    """Tests for process_batch_parallel() method."""

    def test_parallel_batch_basic(self):
        """Test basic parallel batch processing."""
        df = pd.DataFrame({
            'SMILES': ['C', 'CC', 'CCO', 'c1ccccc1'],
            'Name': ['Methane', 'Ethane', 'Ethanol', 'Benzene']
        })

        result = MolecularCalculator.process_batch_parallel(df, 'SMILES')

        assert 'Molecular_Weight' in result.columns
        assert len(result) == 4
        # Check ethanol MW
        assert abs(result.loc[2, 'Molecular_Weight'] - 46.069) < 0.01

    def test_parallel_batch_with_selection(self):
        """Test parallel batch with property selection."""
        df = pd.DataFrame({
            'SMILES': ['CCO', 'c1ccccc1'],
            'Name': ['Ethanol', 'Benzene']
        })

        result = MolecularCalculator.process_batch_parallel(
            df,
            'SMILES',
            selected_properties={'Molecular_Weight', 'LogP'}
        )

        assert 'Molecular_Weight' in result.columns
        assert 'LogP' in result.columns

    def test_parallel_batch_handles_invalid(self):
        """Test parallel batch handles invalid SMILES."""
        df = pd.DataFrame({
            'SMILES': ['CCO', 'invalid_smiles', 'c1ccccc1'],
            'Name': ['Ethanol', 'Invalid', 'Benzene']
        })

        result = MolecularCalculator.process_batch_parallel(df, 'SMILES')

        assert len(result) == 3
        # Valid molecules should have MW
        assert pd.notna(result.loc[0, 'Molecular_Weight'])
        assert pd.notna(result.loc[2, 'Molecular_Weight'])

    def test_parallel_batch_handles_nan(self):
        """Test parallel batch handles NaN values."""
        df = pd.DataFrame({
            'SMILES': ['CCO', None, 'c1ccccc1'],
            'Name': ['Ethanol', 'Missing', 'Benzene']
        })

        result = MolecularCalculator.process_batch_parallel(df, 'SMILES')

        assert len(result) == 3
        assert pd.notna(result.loc[0, 'Molecular_Weight'])
        assert pd.notna(result.loc[2, 'Molecular_Weight'])

    def test_parallel_batch_progress_callback(self):
        """Test progress callback is called."""
        df = pd.DataFrame({
            'SMILES': ['C', 'CC', 'CCO'],
            'Name': ['Methane', 'Ethane', 'Ethanol']
        })

        progress_calls = []

        def callback(completed, total):
            progress_calls.append((completed, total))

        MolecularCalculator.process_batch_parallel(
            df,
            'SMILES',
            progress_callback=callback
        )

        assert len(progress_calls) == 3
        assert progress_calls[-1] == (3, 3)

    def test_parallel_batch_preserves_order(self):
        """Test parallel processing preserves row order."""
        df = pd.DataFrame({
            'SMILES': ['C', 'CC', 'CCC', 'CCCC', 'CCCCC'],
            'ID': [1, 2, 3, 4, 5]
        })

        result = MolecularCalculator.process_batch_parallel(df, 'SMILES')

        # Check order is preserved
        assert list(result['ID']) == [1, 2, 3, 4, 5]
        # Check MW increases with chain length
        mws = result['Molecular_Weight'].tolist()
        assert mws == sorted(mws)


# ============================================================================
# Input Sanitization Tests
# ============================================================================

class TestSanitization:
    """Tests for input sanitization functions."""

    # SMILES Sanitization
    def test_sanitize_smiles_valid(self):
        """Test sanitizing valid SMILES."""
        assert sanitize_smiles("CCO") == "CCO"
        assert sanitize_smiles("c1ccccc1") == "c1ccccc1"
        assert sanitize_smiles("  CCO  ") == "CCO"

    def test_sanitize_smiles_rejects_invalid_chars(self):
        """Test that SMILES with invalid characters are rejected."""
        # SMILES with invalid characters should be rejected (return None)
        # to prevent corrupted data from being processed
        result = sanitize_smiles("CCO!")
        assert result is None  # Invalid char '!' causes rejection

    def test_sanitize_smiles_empty(self):
        """Test sanitizing empty SMILES."""
        assert sanitize_smiles("") is None
        assert sanitize_smiles("   ") is None
        assert sanitize_smiles(None) is None

    def test_sanitize_smiles_too_long(self):
        """Test sanitizing overly long SMILES."""
        long_smiles = "C" * 6000
        assert sanitize_smiles(long_smiles) is None

    # InChI Sanitization
    def test_sanitize_inchi_valid(self):
        """Test sanitizing valid InChI."""
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        assert sanitize_inchi(inchi) == inchi

    def test_sanitize_inchi_with_whitespace(self):
        """Test sanitizing InChI with whitespace."""
        inchi = "  InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3  "
        assert sanitize_inchi(inchi) == "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"

    def test_sanitize_inchi_invalid_prefix(self):
        """Test sanitizing InChI without proper prefix."""
        assert sanitize_inchi("1S/C2H6O/c1-2-3") is None
        assert sanitize_inchi("INCHI=1S/C2H6O/c1-2-3") is None

    def test_sanitize_inchi_empty(self):
        """Test sanitizing empty InChI."""
        assert sanitize_inchi("") is None
        assert sanitize_inchi(None) is None

    # InChI Key Sanitization
    def test_sanitize_inchi_key_valid(self):
        """Test sanitizing valid InChI Key."""
        key = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        assert sanitize_inchi_key(key) == key

    def test_sanitize_inchi_key_lowercase(self):
        """Test InChI Key is converted to uppercase."""
        key = "lfqscwfljhtthz-uhfffaoysa-n"
        assert sanitize_inchi_key(key) == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

    def test_sanitize_inchi_key_with_whitespace(self):
        """Test sanitizing InChI Key with whitespace."""
        key = "  LFQSCWFLJHTTHZ-UHFFFAOYSA-N  "
        assert sanitize_inchi_key(key) == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

    def test_sanitize_inchi_key_invalid_format(self):
        """Test sanitizing invalid InChI Key formats."""
        assert sanitize_inchi_key("INVALID-KEY") is None
        assert sanitize_inchi_key("LFQSCWFLJHTTHZ-UHFFFAOYSA") is None  # Too short
        assert sanitize_inchi_key("") is None
        assert sanitize_inchi_key(None) is None

    # Column Name Sanitization
    def test_sanitize_column_name_valid(self):
        """Test sanitizing valid column names."""
        assert sanitize_column_name("SMILES") == "SMILES"
        assert sanitize_column_name("Molecular Weight") == "Molecular_Weight"

    def test_sanitize_column_name_special_chars(self):
        """Test sanitizing column names with special characters."""
        assert sanitize_column_name("Col<>Name") == "ColName"
        assert sanitize_column_name("Test@Column!") == "TestColumn"

    def test_sanitize_column_name_empty(self):
        """Test sanitizing empty column names."""
        assert sanitize_column_name("") == "unnamed"
        assert sanitize_column_name(None) == "unnamed"

    # Filename Sanitization
    def test_sanitize_filename_valid(self):
        """Test sanitizing valid filenames."""
        assert sanitize_filename("data.csv") == "data.csv"
        assert sanitize_filename("results_2024.xlsx") == "results_2024.xlsx"

    def test_sanitize_filename_path_traversal(self):
        """Test that path traversal is prevented."""
        # sanitize_filename extracts just the filename, stripping all path components
        result = sanitize_filename("../../../etc/passwd")
        assert ".." not in result
        assert "/" not in result
        # Should get just the filename portion
        assert result == "passwd"

        result2 = sanitize_filename("..\\windows\\system32")
        assert ".." not in result2
        assert result2 == "system32"

    def test_sanitize_filename_special_chars(self):
        """Test sanitizing filenames with special characters."""
        result = sanitize_filename("file<>:name.csv")
        assert "<" not in result
        assert ">" not in result
        assert ":" not in result

    def test_sanitize_filename_hidden(self):
        """Test that hidden files are not allowed."""
        assert sanitize_filename(".hidden") == "hidden"

    # HTML Sanitization
    def test_sanitize_html_valid(self):
        """Test sanitizing HTML content."""
        assert sanitize_html("Hello World") == "Hello World"
        assert sanitize_html("<script>alert(1)</script>") == "&lt;script&gt;alert(1)&lt;/script&gt;"

    def test_sanitize_html_special_chars(self):
        """Test HTML special characters are escaped."""
        assert "&" in sanitize_html("A & B")
        assert "&lt;" in sanitize_html("a < b")
        assert "&gt;" in sanitize_html("a > b")
        assert "&quot;" in sanitize_html('say "hi"')

    def test_sanitize_html_empty(self):
        """Test sanitizing empty HTML."""
        assert sanitize_html("") == ""
        assert sanitize_html(None) == ""


# ============================================================================
# Logging Setup Tests
# ============================================================================

class TestLoggingSetup:
    """Tests for logging configuration."""

    def test_setup_logging_import(self):
        """Test that setup_logging can be imported."""
        from molecular_calculator.config import setup_logging
        assert callable(setup_logging)

    def test_setup_logging_from_monitoring(self):
        """Test that monitoring.setup_logging works."""
        from molecular_calculator.utils.monitoring import setup_logging
        assert callable(setup_logging)

    def test_setup_logging_is_idempotent(self):
        """Test that setup_logging can be called multiple times."""
        from molecular_calculator.utils.monitoring import setup_logging
        import logging

        # Should not raise
        setup_logging(level=logging.INFO)
        setup_logging(level=logging.DEBUG)

        # Should still work
        logger = logging.getLogger(__name__)
        logger.info("Test message")
