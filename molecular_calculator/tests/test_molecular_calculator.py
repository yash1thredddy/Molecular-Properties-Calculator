"""Tests for the refactored MolecularCalculator.

This module tests the core molecular calculator functionality including:
- Property calculation
- Format conversion
- Batch processing
"""

import pytest
import pandas as pd

from molecular_calculator.core import MolecularCalculator
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
