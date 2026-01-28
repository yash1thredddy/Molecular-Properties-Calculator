"""Refactored Molecular Calculator - Main Calculator Class.

This module provides the main MolecularCalculator class that orchestrates
property calculation, format conversion, and batch processing.

This is a facade that delegates to specialized services while maintaining
backwards compatibility with the original API.

Usage:
    # Preferred: Use the default singleton instance
    >>> from molecular_calculator import get_calculator
    >>> calc = get_calculator()
    >>> result = calc.calculate("CCO")

    # Alternative: Create your own instance
    >>> calc = MolecularCalculator()
    >>> result = calc.calculate("CCO")

    # Legacy: Static methods (kept for backwards compatibility)
    >>> props = MolecularCalculator.calculate_molecular_properties("CCO")
"""

import logging
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Any, Optional, List, Set, Callable

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


def _process_single_molecule(
    idx_and_smiles: tuple,
    conversion_service: ConversionService,
    property_calculator: PropertyCalculator,
    selected_properties: Optional[Set[str]],
    enable_online_lookup: bool
) -> tuple:
    """Process a single molecule - designed for parallel execution.

    This helper function is extracted from process_batch_parallel to improve
    code organization and testability.

    Args:
        idx_and_smiles: Tuple of (original_index, smiles_value)
        conversion_service: Service for format detection and conversion
        property_calculator: Service for property calculation
        selected_properties: Optional set of properties to calculate
        enable_online_lookup: Whether to enable online InChI Key lookup

    Returns:
        Tuple of (original_index, properties_dict or error_dict)
    """
    idx, smiles_value = idx_and_smiles

    if pd.isna(smiles_value):
        return idx, {'_error': 'Empty or NaN value'}

    smiles_str = str(smiles_value)

    # Auto-detect and convert if needed
    input_format = conversion_service.detect_format(smiles_str)

    if input_format != InputFormat.SMILES:
        conversion_result = conversion_service.to_smiles(
            smiles_str,
            input_format=input_format,
            enable_online_lookup=enable_online_lookup
        )
        if conversion_result.success:
            smiles_str = conversion_result.smiles
        else:
            return idx, {'_error': conversion_result.error or 'Conversion failed'}

    # Calculate properties
    properties = property_calculator.calculate_as_dict(
        smiles_str,
        selected_properties
    )
    return idx, properties


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

        # Apply suppress_warnings setting
        if suppress_warnings:
            self.suppress_rdkit_warnings(True)

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
            DataFrame with calculated properties.
            Returns empty DataFrame if input is None or empty.

        Raises:
            ValueError: If column not found in DataFrame
        """
        # Validate DataFrame is not empty
        if df is None or df.empty:
            logger.warning("Empty DataFrame provided to process_batch")
            return df if df is not None else pd.DataFrame()

        # Validate column exists
        if smiles_col not in df.columns:
            raise ValueError(f"Column '{smiles_col}' not found in DataFrame")

        calculator = cls()
        results = []

        for idx, row in df.iterrows():
            smiles = row[smiles_col]

            # Skip empty values
            if pd.isna(smiles):
                logger.warning(f"Row {idx}: Empty or NaN value in column '{smiles_col}'")
                results.append({'_error': 'Empty or NaN value'})
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
                    error_msg = conversion_result.error or 'Conversion failed'
                    logger.warning(f"Row {idx}: Conversion failed - {error_msg}")
                    results.append({'_error': error_msg})
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

    @classmethod
    def process_batch_parallel(
        cls,
        df: pd.DataFrame,
        smiles_col: str,
        selected_properties: Set[str] = None,
        enable_online_lookup: bool = True,
        max_workers: int = None,
        progress_callback: Callable[[int, int], None] = None
    ) -> pd.DataFrame:
        """Process batch of molecules with parallel execution.

        This is an optimized version of process_batch that uses ThreadPoolExecutor
        for parallel processing. Typically 2-4x faster for large datasets.

        Args:
            df: DataFrame containing molecules
            smiles_col: Name of column containing SMILES
            selected_properties: Set of properties to calculate (None for all)
            enable_online_lookup: Allow API calls for InChI Key conversion
            max_workers: Maximum parallel workers (default: min(32, cpu_count + 4))
            progress_callback: Optional callback(completed, total) for progress

        Returns:
            DataFrame with calculated properties

        Raises:
            ValueError: If DataFrame is empty or column not found
        """
        import os

        # Validate DataFrame is not empty
        if df is None or df.empty:
            logger.warning("Empty DataFrame provided to process_batch_parallel")
            return df if df is not None else pd.DataFrame()

        # Validate column exists
        if smiles_col not in df.columns:
            raise ValueError(f"Column '{smiles_col}' not found in DataFrame")

        if max_workers is None:
            max_workers = min(32, (os.cpu_count() or 1) + 4)

        # Get services once for reuse
        conversion_service = get_conversion_service()
        property_calculator = get_property_calculator()

        # Create a wrapper that captures the closure variables for the extracted function
        def process_single(idx_and_smiles: tuple) -> tuple:
            """Wrapper for _process_single_molecule with captured context."""
            return _process_single_molecule(
                idx_and_smiles,
                conversion_service,
                property_calculator,
                selected_properties,
                enable_online_lookup
            )

        # Prepare work items with position index for O(1) lookup
        work_items = []
        for pos, (idx, row) in enumerate(df.iterrows()):
            work_items.append((pos, idx, row[smiles_col]))

        total = len(work_items)
        results = [None] * total  # Pre-allocate for ordered results

        # Process in parallel
        completed = 0
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Map future to position index for O(1) lookup
            futures = {
                executor.submit(process_single, (item[1], item[2])): item[0]
                for item in work_items
            }

            for future in as_completed(futures):
                # Get position directly from futures mapping - O(1) lookup
                position = futures[future]

                try:
                    idx, props = future.result()
                    results[position] = props
                except Exception as e:
                    # Log error and continue processing other molecules
                    logger.warning(f"Error processing molecule at position {position}: {e}")
                    results[position] = {'_error': str(e)}

                completed += 1

                if progress_callback:
                    progress_callback(completed, total)

        # Create results DataFrame
        results_df = pd.DataFrame(results)

        # Combine with original DataFrame
        final_df = pd.concat([df.reset_index(drop=True), results_df], axis=1)

        return final_df


# =============================================================================
# Module-level singleton for convenience
# =============================================================================

_default_calculator: Optional[MolecularCalculator] = None
_calculator_lock = threading.Lock()


def get_calculator() -> MolecularCalculator:
    """Get the default MolecularCalculator singleton instance.

    This is the recommended way to use the calculator for most use cases.
    The instance is lazily created on first call and reused thereafter.

    Thread-safe using double-checked locking pattern.

    Returns:
        MolecularCalculator: The default calculator instance

    Example:
        >>> calc = get_calculator()
        >>> result = calc.calculate("CCO")
        >>> print(result.properties.molecular_weight)
    """
    global _default_calculator
    if _default_calculator is None:
        with _calculator_lock:
            # Double-check after acquiring lock
            if _default_calculator is None:
                _default_calculator = MolecularCalculator()
    return _default_calculator
