"""Chemical structure format conversion service.

This module provides services for converting between different
molecular representation formats (SMILES, InChI, InChI Key).
"""

import re
import logging
from typing import Optional

from rdkit import Chem

from molecular_calculator.models import InputFormat, ConversionResult
from molecular_calculator.services.api_client import get_api_client, ChemicalAPIClient
from molecular_calculator.utils.exceptions import ConversionError, InvalidSMILESError

logger = logging.getLogger(__name__)


class ConversionService:
    """Service for converting between molecular formats.

    Supports conversion from InChI and InChI Key formats to SMILES.
    InChI Key conversion requires external API lookup.

    Example:
        >>> service = ConversionService()
        >>> result = service.to_smiles("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
        >>> if result.success:
        ...     print(f"SMILES: {result.smiles}")
    """

    # InChI Key pattern: exactly 27 characters in format XXXXXXXXXXXXXX-YYYYYYYYYY-Z
    INCHI_KEY_PATTERN = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')

    def __init__(self, api_client: ChemicalAPIClient = None):
        """Initialize the conversion service.

        Args:
            api_client: Optional API client for InChI Key lookups
        """
        self._api_client = api_client

    @property
    def api_client(self) -> ChemicalAPIClient:
        """Get the API client (lazy initialization)."""
        if self._api_client is None:
            self._api_client = get_api_client()
        return self._api_client

    def detect_format(self, input_text: str) -> InputFormat:
        """Auto-detect the input format of a molecular structure string.

        Args:
            input_text: Input string to analyze

        Returns:
            Detected InputFormat enum value

        Example:
            >>> service = ConversionService()
            >>> service.detect_format("CCO")
            InputFormat.SMILES
            >>> service.detect_format("InChI=1S/C2H6O/c1-2-3")
            InputFormat.INCHI
        """
        if not input_text or not isinstance(input_text, str):
            return InputFormat.UNKNOWN

        input_text = input_text.strip()

        # Check InChI Key first (most specific format)
        if self.INCHI_KEY_PATTERN.match(input_text.upper()):
            return InputFormat.INCHI_KEY

        # Check InChI
        if input_text.startswith('InChI='):
            return InputFormat.INCHI

        # Assume SMILES for anything else that looks valid
        # A more robust check would use RDKit, but we want to be fast here
        if input_text and not input_text.startswith('<'):
            return InputFormat.SMILES

        return InputFormat.UNKNOWN

    def _convert_inchi_to_smiles(self, inchi: str) -> ConversionResult:
        """Convert InChI string to SMILES.

        Args:
            inchi: InChI string

        Returns:
            ConversionResult with SMILES if successful
        """
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                return ConversionResult(
                    success=False,
                    source_format=InputFormat.INCHI,
                    error="RDKit could not parse InChI"
                )

            # Sanitize molecule to resolve stereochemical conflicts
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                # If sanitization fails, try without stereo
                try:
                    Chem.RemoveStereochemistry(mol)
                    Chem.SanitizeMol(mol)
                except Exception as e:
                    return ConversionResult(
                        success=False,
                        source_format=InputFormat.INCHI,
                        error=f"Molecule sanitization failed: {e}"
                    )

            smiles = Chem.MolToSmiles(mol)
            logger.debug(f"Converted InChI to SMILES: {smiles}")

            return ConversionResult(
                success=True,
                smiles=smiles,
                source_format=InputFormat.INCHI
            )

        except Exception as e:
            logger.error(f"InChI conversion failed: {e}")
            return ConversionResult(
                success=False,
                source_format=InputFormat.INCHI,
                error=str(e)
            )

    def _convert_inchi_key_to_smiles(
        self,
        inchi_key: str,
        enable_online_lookup: bool = True
    ) -> ConversionResult:
        """Convert InChI Key to SMILES via API lookup.

        Args:
            inchi_key: InChI Key string
            enable_online_lookup: Whether to allow API calls

        Returns:
            ConversionResult with SMILES if successful
        """
        if not enable_online_lookup:
            return ConversionResult(
                success=False,
                source_format=InputFormat.INCHI_KEY,
                error="Online lookup disabled - cannot convert InChI Key"
            )

        try:
            response = self.api_client.inchi_key_to_smiles(inchi_key)

            if response.success:
                logger.debug(f"Converted InChI Key via {response.source}: {response.data}")
                return ConversionResult(
                    success=True,
                    smiles=response.data,
                    source_format=InputFormat.INCHI_KEY
                )
            else:
                return ConversionResult(
                    success=False,
                    source_format=InputFormat.INCHI_KEY,
                    error=response.error
                )

        except Exception as e:
            logger.error(f"InChI Key conversion failed: {e}")
            return ConversionResult(
                success=False,
                source_format=InputFormat.INCHI_KEY,
                error=str(e)
            )

    def to_smiles(
        self,
        input_text: str,
        input_format: InputFormat = None,
        enable_online_lookup: bool = True
    ) -> ConversionResult:
        """Convert input to SMILES format.

        Automatically detects input format if not specified.

        Args:
            input_text: Input molecular structure string
            input_format: Optional explicit format (auto-detected if None)
            enable_online_lookup: Allow API calls for InChI Key conversion

        Returns:
            ConversionResult with SMILES if successful

        Example:
            >>> service = ConversionService()
            >>> result = service.to_smiles("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
            >>> print(result.smiles)  # "CCO"
        """
        if not input_text or not isinstance(input_text, str):
            return ConversionResult(
                success=False,
                error="Input text is empty or invalid"
            )

        input_text = input_text.strip()

        # Auto-detect format if not specified
        if input_format is None:
            input_format = self.detect_format(input_text)

        # Handle each format
        if input_format == InputFormat.SMILES:
            # Validate SMILES
            mol = Chem.MolFromSmiles(input_text)
            if mol is not None:
                return ConversionResult(
                    success=True,
                    smiles=input_text,
                    source_format=InputFormat.SMILES
                )
            else:
                return ConversionResult(
                    success=False,
                    source_format=InputFormat.SMILES,
                    error="Invalid SMILES string"
                )

        elif input_format == InputFormat.INCHI:
            return self._convert_inchi_to_smiles(input_text)

        elif input_format == InputFormat.INCHI_KEY:
            return self._convert_inchi_key_to_smiles(
                input_text,
                enable_online_lookup=enable_online_lookup
            )

        else:
            return ConversionResult(
                success=False,
                source_format=InputFormat.UNKNOWN,
                error="Unknown input format"
            )

    def validate_smiles(self, smiles: str) -> bool:
        """Validate that a SMILES string is parseable.

        Args:
            smiles: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        if not smiles or not isinstance(smiles, str):
            return False

        try:
            mol = Chem.MolFromSmiles(smiles.strip())
            return mol is not None
        except Exception:
            return False


# Singleton instance for convenience
_conversion_service: Optional[ConversionService] = None


def get_conversion_service() -> ConversionService:
    """Get the shared conversion service instance.

    Returns:
        ConversionService singleton instance
    """
    global _conversion_service
    if _conversion_service is None:
        _conversion_service = ConversionService()
    return _conversion_service
