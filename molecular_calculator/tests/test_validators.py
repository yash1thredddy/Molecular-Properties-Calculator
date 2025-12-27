"""Tests for validation utilities.

This module tests the validators in molecular_calculator.utils.validators.
"""

import pytest
import pandas as pd

from molecular_calculator.utils.validators import (
    InputValidator,
    FileValidator,
    DataFrameValidator,
    ValidationResult,
)


class TestValidationResult:
    """Tests for ValidationResult class."""

    def test_valid_result(self):
        """Test creating a valid result."""
        result = ValidationResult(is_valid=True)
        assert result.is_valid
        assert len(result.errors) == 0
        assert len(result.warnings) == 0

    def test_invalid_result(self):
        """Test creating an invalid result."""
        result = ValidationResult(
            is_valid=False,
            errors=["Error 1", "Error 2"]
        )
        assert not result.is_valid
        assert len(result.errors) == 2

    def test_bool_conversion(self):
        """Test boolean conversion."""
        valid = ValidationResult(is_valid=True)
        invalid = ValidationResult(is_valid=False)

        assert bool(valid) is True
        assert bool(invalid) is False

        # Can use in if statements
        if valid:
            passed = True
        else:
            passed = False
        assert passed

    def test_add_error(self):
        """Test adding errors."""
        result = ValidationResult(is_valid=True)
        assert result.is_valid

        result.add_error("Something went wrong")
        assert not result.is_valid
        assert "Something went wrong" in result.errors

    def test_add_warning(self):
        """Test adding warnings doesn't affect validity."""
        result = ValidationResult(is_valid=True)

        result.add_warning("This might be slow")
        assert result.is_valid  # Still valid
        assert "This might be slow" in result.warnings


class TestInputValidator:
    """Tests for InputValidator class."""

    # ========================================================================
    # SMILES Validation Tests
    # ========================================================================

    def test_validate_smiles_valid(self, sample_smiles):
        """Test validation of valid SMILES strings."""
        for smiles in sample_smiles:
            result = InputValidator.validate_smiles(smiles)
            assert result.is_valid, f"Expected valid: {smiles}"
            assert len(result.errors) == 0

    def test_validate_smiles_empty(self):
        """Test validation of empty SMILES."""
        result = InputValidator.validate_smiles('')
        assert not result.is_valid
        assert 'empty' in result.errors[0].lower()

    def test_validate_smiles_none(self):
        """Test validation of None SMILES."""
        result = InputValidator.validate_smiles(None)
        assert not result.is_valid

    def test_validate_smiles_whitespace(self):
        """Test validation of whitespace-only SMILES."""
        result = InputValidator.validate_smiles('   ')
        assert not result.is_valid

    def test_validate_smiles_too_long(self):
        """Test validation of excessively long SMILES."""
        long_smiles = 'C' * 11000
        result = InputValidator.validate_smiles(long_smiles)
        assert not result.is_valid
        assert 'too long' in result.errors[0].lower()

    def test_validate_smiles_with_spaces(self):
        """Test that leading/trailing spaces are handled."""
        result = InputValidator.validate_smiles('  CCO  ')
        assert result.is_valid

    def test_validate_smiles_invalid_characters(self):
        """Test validation catches invalid characters."""
        result = InputValidator.validate_smiles('CCO<script>')
        assert not result.is_valid
        assert 'invalid characters' in result.errors[0].lower()

    # ========================================================================
    # InChI Validation Tests
    # ========================================================================

    def test_validate_inchi_valid(self, sample_inchi):
        """Test validation of valid InChI."""
        result = InputValidator.validate_inchi(sample_inchi)
        assert result.is_valid

    def test_validate_inchi_empty(self):
        """Test validation of empty InChI."""
        result = InputValidator.validate_inchi('')
        assert not result.is_valid

    def test_validate_inchi_missing_prefix(self):
        """Test validation requires InChI= prefix."""
        result = InputValidator.validate_inchi('1S/C2H6O/c1-2-3/h3H,2H2,1H3')
        assert not result.is_valid
        assert "InChI=" in result.errors[0]

    # ========================================================================
    # InChI Key Validation Tests
    # ========================================================================

    def test_validate_inchi_key_valid(self, sample_inchi_key):
        """Test validation of valid InChI Key."""
        result = InputValidator.validate_inchi_key(sample_inchi_key)
        assert result.is_valid

    def test_validate_inchi_key_invalid(self, invalid_inchi_keys):
        """Test validation of invalid InChI Keys."""
        for key in invalid_inchi_keys:
            if key:  # Skip empty string test (separate test)
                result = InputValidator.validate_inchi_key(key)
                assert not result.is_valid, f"Expected invalid: {key}"

    def test_validate_inchi_key_empty(self):
        """Test validation of empty InChI Key."""
        result = InputValidator.validate_inchi_key('')
        assert not result.is_valid

    def test_validate_inchi_key_lowercase(self):
        """Test that lowercase InChI Keys are validated (case normalized)."""
        result = InputValidator.validate_inchi_key(
            'lfqscwfljhtthz-uhfffaoysa-n'
        )
        # Should be valid because we normalize to uppercase
        assert result.is_valid

    # ========================================================================
    # Format Detection Tests
    # ========================================================================

    def test_detect_format_smiles(self):
        """Test format detection for SMILES."""
        assert InputValidator.detect_format('CCO') == 'smiles'
        assert InputValidator.detect_format('c1ccccc1') == 'smiles'

    def test_detect_format_inchi(self, sample_inchi):
        """Test format detection for InChI."""
        assert InputValidator.detect_format(sample_inchi) == 'inchi'

    def test_detect_format_inchi_key(self, sample_inchi_key):
        """Test format detection for InChI Key."""
        assert InputValidator.detect_format(sample_inchi_key) == 'inchi_key'

    def test_detect_format_unknown(self):
        """Test format detection for unknown formats."""
        assert InputValidator.detect_format('') == 'unknown'
        assert InputValidator.detect_format(None) == 'unknown'
        assert InputValidator.detect_format('random text here') == 'unknown'

    def test_detect_format_edge_cases(self):
        """Test format detection edge cases."""
        # InChI Key takes precedence (more specific pattern)
        assert InputValidator.detect_format(
            'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
        ) == 'inchi_key'

    # ========================================================================
    # Sanitization Tests
    # ========================================================================

    def test_sanitize_html_script(self, xss_payloads):
        """Test HTML sanitization removes script tags."""
        for payload in xss_payloads:
            sanitized = InputValidator.sanitize_html(payload)
            assert '<script>' not in sanitized
            assert 'javascript:' not in sanitized.lower() or '&' in sanitized

    def test_sanitize_html_preserves_safe_text(self):
        """Test sanitization preserves normal text."""
        text = "Hello, world! This is a test."
        sanitized = InputValidator.sanitize_html(text)
        assert sanitized == text

    def test_sanitize_html_escapes_angles(self):
        """Test angle brackets are escaped."""
        text = "<div>content</div>"
        sanitized = InputValidator.sanitize_html(text)
        assert '&lt;' in sanitized
        assert '&gt;' in sanitized

    def test_sanitize_html_handles_non_string(self):
        """Test sanitization handles non-string input."""
        assert InputValidator.sanitize_html(123) == '123'
        assert InputValidator.sanitize_html(None) == 'None'

    # ========================================================================
    # Safety Check Tests
    # ========================================================================

    def test_is_safe_input_normal(self):
        """Test safe input detection for normal text."""
        assert InputValidator.is_safe_input('CCO')
        assert InputValidator.is_safe_input('Normal text')
        assert InputValidator.is_safe_input('file.csv')

    def test_is_safe_input_sql_injection(self, sql_injection_payloads):
        """Test unsafe input detection for SQL injection."""
        for payload in sql_injection_payloads:
            assert not InputValidator.is_safe_input(payload), \
                f"Should detect: {payload}"

    def test_is_safe_input_xss(self):
        """Test unsafe input detection for XSS."""
        assert not InputValidator.is_safe_input('<script>alert(1)</script>')
        assert not InputValidator.is_safe_input('javascript:alert(1)')


class TestFileValidator:
    """Tests for FileValidator class."""

    def test_validate_upload_valid_csv(self, mock_csv_file):
        """Test validation of valid CSV file."""
        result = FileValidator.validate_upload(mock_csv_file)
        assert result.is_valid

    def test_validate_upload_valid_xlsx(self, mock_xlsx_file):
        """Test validation of valid XLSX file."""
        result = FileValidator.validate_upload(mock_xlsx_file)
        assert result.is_valid

    def test_validate_upload_none(self):
        """Test validation of None file."""
        result = FileValidator.validate_upload(None)
        assert not result.is_valid
        assert 'No file uploaded' in result.errors[0]

    def test_validate_upload_too_large(self, mock_large_file):
        """Test validation of oversized file."""
        result = FileValidator.validate_upload(mock_large_file)
        assert not result.is_valid
        assert 'too large' in result.errors[0].lower()

    def test_validate_upload_invalid_type(self, mock_invalid_file):
        """Test validation of invalid file type."""
        result = FileValidator.validate_upload(mock_invalid_file)
        assert not result.is_valid
        assert 'Invalid file type' in result.errors[0]

    def test_sanitize_filename_normal(self):
        """Test filename sanitization for normal names."""
        assert FileValidator.sanitize_filename('data.csv') == 'data.csv'
        assert FileValidator.sanitize_filename('my_file.xlsx') == 'my_file.xlsx'

    def test_sanitize_filename_path_traversal(self, path_traversal_payloads):
        """Test filename sanitization removes path traversal."""
        for payload in path_traversal_payloads:
            sanitized = FileValidator.sanitize_filename(payload)
            assert '..' not in sanitized
            assert '/' not in sanitized
            assert '\\' not in sanitized

    def test_sanitize_filename_special_chars(self):
        """Test filename sanitization removes special characters."""
        sanitized = FileValidator.sanitize_filename('file<>:"|?*.csv')
        assert '<' not in sanitized
        assert '>' not in sanitized
        assert ':' not in sanitized

    def test_sanitize_filename_empty(self):
        """Test filename sanitization handles empty string."""
        assert FileValidator.sanitize_filename('') == 'unnamed'
        assert FileValidator.sanitize_filename(None) == 'unnamed'


class TestDataFrameValidator:
    """Tests for DataFrameValidator class."""

    def test_validate_valid_df(self, sample_dataframe):
        """Test validation of valid DataFrame."""
        result = DataFrameValidator.validate(sample_dataframe)
        assert result.is_valid

    def test_validate_empty_df(self, empty_dataframe):
        """Test validation of empty DataFrame."""
        result = DataFrameValidator.validate(empty_dataframe)
        assert not result.is_valid
        assert 'empty' in result.errors[0].lower()

    def test_validate_none_df(self):
        """Test validation of None DataFrame."""
        result = DataFrameValidator.validate(None)
        assert not result.is_valid

    def test_validate_large_df_warning(self, large_dataframe):
        """Test validation adds warning for large DataFrames."""
        result = DataFrameValidator.validate(large_dataframe)
        # May or may not be valid depending on limit
        # But should have warning
        assert len(result.warnings) > 0 or not result.is_valid

    def test_validate_required_columns_present(self, sample_dataframe):
        """Test validation with required columns present."""
        result = DataFrameValidator.validate(
            sample_dataframe,
            required_columns=['SMILES', 'Name']
        )
        assert result.is_valid

    def test_validate_required_columns_missing(self, sample_dataframe):
        """Test validation with required columns missing."""
        result = DataFrameValidator.validate(
            sample_dataframe,
            required_columns=['SMILES', 'NonExistent']
        )
        assert not result.is_valid
        assert 'NonExistent' in result.errors[0]

    def test_column_exists(self, sample_dataframe):
        """Test column existence check."""
        assert DataFrameValidator.column_exists(sample_dataframe, 'SMILES')
        assert DataFrameValidator.column_exists(sample_dataframe, 'Name')
        assert not DataFrameValidator.column_exists(sample_dataframe, 'Missing')

    def test_column_exists_none_df(self):
        """Test column existence with None DataFrame."""
        assert not DataFrameValidator.column_exists(None, 'SMILES')

    def test_get_numeric_columns(self, sample_dataframe):
        """Test getting numeric columns."""
        numeric = DataFrameValidator.get_numeric_columns(sample_dataframe)
        assert 'MW' in numeric
        assert 'LogP' in numeric
        assert 'SMILES' not in numeric
        assert 'Name' not in numeric

    def test_get_categorical_columns(self, sample_dataframe):
        """Test getting categorical columns."""
        categorical = DataFrameValidator.get_categorical_columns(sample_dataframe)
        assert 'Category' in categorical
        assert 'Name' in categorical
        assert 'MW' not in categorical
