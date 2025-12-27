"""Refactored Molecular Calculator - Main Calculator Class.

This module provides the main MolecularCalculator class that orchestrates
property calculation, format conversion, and batch processing.

This is a facade that delegates to specialized services while maintaining
backwards compatibility with the original API.
"""

import logging
from typing import Dict, Any, Optional, List, Set

import pandas as pd

from molecular_calculator.models import (
    InputFormat,
    MolecularProperties,
    CalculationResult,
    ConversionResult,
    PROPERTY_GROUPS,
    LEI_PROPERTY_GROUP,
)
from molecular_calculator.services import (
    ConversionService,
    PropertyCalculator,
    get_conversion_service,
    get_property_calculator,
)
from molecular_calculator.utils.exceptions import InvalidSMILESError

logger = logging.getLogger(__name__)


class MolecularCalculator:
    """Main class for molecular property calculations.

    This class provides a unified interface for:
    - Converting between molecular formats (SMILES, InChI, InChI Key)
    - Calculating molecular properties
    - Processing batch data

    The class delegates to specialized services internally but maintains
    a simple API for users.

    Example:
        >>> calc = MolecularCalculator()
        >>> props = calc.calculate_molecular_properties("CCO")
        >>> print(props['Molecular_Weight'])
        46.069
    """

    def __init__(
        self,
        conversion_service: ConversionService = None,
        property_calculator: PropertyCalculator = None,
        suppress_warnings: bool = True
    ):
        """Initialize the molecular calculator.

        Args:
            conversion_service: Optional custom conversion service
            property_calculator: Optional custom property calculator
            suppress_warnings: Whether to suppress RDKit warnings
        """
        self._conversion_service = conversion_service
        self._property_calculator = property_calculator
        self._suppress_warnings = suppress_warnings

    @property
    def conversion_service(self) -> ConversionService:
        """Get the conversion service (lazy initialization)."""
        if self._conversion_service is None:
            self._conversion_service = get_conversion_service()
        return self._conversion_service

    @property
    def property_calculator(self) -> PropertyCalculator:
        """Get the property calculator (lazy initialization)."""
        if self._property_calculator is None:
            self._property_calculator = get_property_calculator()
        return self._property_calculator

    # =========================================================================
    # Static methods for backwards compatibility
    # =========================================================================

    @staticmethod
    def suppress_rdkit_warnings(suppress: bool = True) -> None:
        """Suppress or enable RDKit warning messages.

        Args:
            suppress: If True, suppress warnings; if False, enable warnings
        """
        if suppress:
            PropertyCalculator.disable_rdkit_warnings()
        else:
            PropertyCalculator.enable_rdkit_warnings()

    @staticmethod
    def detect_input_format(input_text: str) -> str:
        """Auto-detect input format.

        Args:
            input_text: Input molecular structure string

        Returns:
            Detected format ('inchi', 'inchi_key', or 'smiles')
        """
        service = get_conversion_service()
        format_enum = service.detect_format(input_text)
        return format_enum.value

    @staticmethod
    def convert_inchi_key_to_smiles(
        inchi_key: str,
        timeout: int = 10
    ) -> Optional[str]:
        """Convert InChI Key to SMILES using Chemical Identifier Resolver.

        Args:
            inchi_key: InChI Key string
            timeout: Request timeout in seconds

        Returns:
            SMILES string or None if conversion fails
        """
        from molecular_calculator.services.api_client import ChemicalAPIClient

        client = ChemicalAPIClient(timeout=timeout)
        response = client.inchi_key_to_smiles(inchi_key)

        return response.data if response.success else None

    @staticmethod
    def convert_to_smiles(
        input_text: str,
        input_type: str,
        enable_online_lookup: bool = True
    ) -> Optional[str]:
        """Convert InChI or InChI Key to SMILES.

        Args:
            input_text: Input molecular structure string
            input_type: Type of input ('smiles', 'inchi', 'inchi_key')
            enable_online_lookup: Allow online database lookup for InChI Keys

        Returns:
            SMILES string or None if conversion fails
        """
        # Map string type to enum
        format_map = {
            'smiles': InputFormat.SMILES,
            'inchi': InputFormat.INCHI,
            'inchi_key': InputFormat.INCHI_KEY,
        }
        input_format = format_map.get(input_type.lower(), InputFormat.UNKNOWN)

        service = get_conversion_service()
        result = service.to_smiles(
            input_text,
            input_format=input_format,
            enable_online_lookup=enable_online_lookup
        )

        return result.smiles if result.success else None

    @staticmethod
    def calculate_molecular_properties(smiles: str) -> Dict[str, Any]:
        """Calculate various molecular properties from SMILES.

        Args:
            smiles: SMILES string

        Returns:
            Dictionary containing calculated properties
        """
        if not smiles:
            return {}

        calculator = get_property_calculator()
        return calculator.calculate_as_dict(smiles)

    @staticmethod
    def detect_smiles_column(df: pd.DataFrame) -> Optional[str]:
        """Detect SMILES column with case-insensitive matching.

        Args:
            df: Pandas DataFrame

        Returns:
            Column name containing SMILES or None if not found
        """
        possible_names = [
            'smiles', 'SMILES', 'Smiles', 'smi', 'SMI',
            'canonical_smiles', 'CANONICAL_SMILES'
        ]

        for col in df.columns:
            if col in possible_names:
                return col
            if col.lower() in [name.lower() for name in possible_names]:
                return col

        return None

    @staticmethod
    def get_property_groups(include_lei: bool = False) -> Dict[str, List[str]]:
        """Get organized property groups for UI display.

        Args:
            include_lei: Include Ligand Efficiency Indices (requires pKi)

        Returns:
            Dictionary with property groups and their properties
        """
        groups = dict(PROPERTY_GROUPS)

        if include_lei:
            groups.update(LEI_PROPERTY_GROUP)

        return groups

    # =========================================================================
    # Instance methods
    # =========================================================================

    def calculate(self, smiles: str) -> CalculationResult:
        """Calculate properties for a single SMILES.

        Args:
            smiles: SMILES string

        Returns:
            CalculationResult with properties
        """
        return self.property_calculator.calculate(smiles)

    def convert(
        self,
        input_text: str,
        input_format: InputFormat = None,
        enable_online_lookup: bool = True
    ) -> ConversionResult:
        """Convert input to SMILES format.

        Args:
            input_text: Input molecular structure string
            input_format: Optional explicit format (auto-detected if None)
            enable_online_lookup: Allow API calls for InChI Key conversion

        Returns:
            ConversionResult with SMILES if successful
        """
        return self.conversion_service.to_smiles(
            input_text,
            input_format=input_format,
            enable_online_lookup=enable_online_lookup
        )

    def calculate_from_any_format(
        self,
        input_text: str,
        enable_online_lookup: bool = True
    ) -> CalculationResult:
        """Calculate properties from any supported format.

        Automatically detects and converts the input format.

        Args:
            input_text: Input in SMILES, InChI, or InChI Key format
            enable_online_lookup: Allow API calls for InChI Key conversion

        Returns:
            CalculationResult with properties
        """
        # Convert to SMILES first
        conversion_result = self.convert(
            input_text,
            enable_online_lookup=enable_online_lookup
        )

        if not conversion_result.success:
            return CalculationResult(
                success=False,
                error=f"Conversion failed: {conversion_result.error}"
            )

        # Calculate properties
        return self.calculate(conversion_result.smiles)

    @classmethod
    def process_batch(
        cls,
        df: pd.DataFrame,
        smiles_col: str,
        selected_properties: Set[str] = None,
        enable_online_lookup: bool = True
    ) -> pd.DataFrame:
        """Process batch of molecules and calculate properties.

        Args:
            df: DataFrame containing molecules
            smiles_col: Name of column containing SMILES
            selected_properties: Set of properties to calculate (None for all)
            enable_online_lookup: Allow API calls for InChI Key conversion

        Returns:
            DataFrame with calculated properties
        """
        calculator = cls()
        results = []

        for idx, row in df.iterrows():
            smiles = row[smiles_col]

            # Skip empty values
            if pd.isna(smiles):
                results.append({})
                continue

            smiles_str = str(smiles)

            # Auto-detect and convert if needed
            input_format = calculator.conversion_service.detect_format(smiles_str)

            if input_format != InputFormat.SMILES:
                conversion_result = calculator.convert(
                    smiles_str,
                    input_format=input_format,
                    enable_online_lookup=enable_online_lookup
                )
                if conversion_result.success:
                    smiles_str = conversion_result.smiles
                else:
                    results.append({})
                    continue

            # Calculate properties
            properties = calculator.property_calculator.calculate_as_dict(
                smiles_str,
                selected_properties
            )
            results.append(properties)

        # Create results DataFrame
        results_df = pd.DataFrame(results)

        # Combine with original DataFrame
        final_df = pd.concat([df.reset_index(drop=True), results_df], axis=1)

        return final_df
