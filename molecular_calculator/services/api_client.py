"""External API client for chemical database lookups.

This module provides a service for interacting with external chemical
databases like NIH CIR and PubChem for InChI Key to SMILES conversion.
"""

import logging
from typing import Optional
from dataclasses import dataclass

import requests

from molecular_calculator.config.settings import config
from molecular_calculator.utils.exceptions import APIError, RateLimitError

logger = logging.getLogger(__name__)


@dataclass
class APIResponse:
    """Standardized response from API calls.

    Attributes:
        success: Whether the API call succeeded
        data: Response data (typically SMILES string)
        source: Which API returned the result
        error: Error message if call failed
    """
    success: bool
    data: Optional[str] = None
    source: Optional[str] = None
    error: Optional[str] = None


class ChemicalAPIClient:
    """Client for external chemical database APIs.

    Provides methods for looking up chemical structures from InChI Keys
    using NIH CIR and PubChem APIs with automatic fallback.

    Example:
        >>> client = ChemicalAPIClient()
        >>> response = client.inchi_key_to_smiles("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
        >>> if response.success:
        ...     print(f"SMILES: {response.data}")
    """

    # API endpoints
    NIH_CIR_URL = "https://cactus.nci.nih.gov/chemical/structure/{}/smiles"
    PUBCHEM_URL = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}"
        "/property/IsomericSMILES/JSON"
    )

    def __init__(self, timeout: int = None):
        """Initialize the API client.

        Args:
            timeout: Request timeout in seconds (uses config default if None)
        """
        self.timeout = timeout or config.API_TIMEOUT

    def _validate_smiles(self, smiles: str) -> bool:
        """Validate that a SMILES string is parseable.

        Uses RDKit to verify the SMILES can be parsed into a molecule.

        Args:
            smiles: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception:
            return False

    def _query_nih_cir(self, inchi_key: str) -> APIResponse:
        """Query NIH Chemical Identifier Resolver.

        Args:
            inchi_key: InChI Key to look up

        Returns:
            APIResponse with result
        """
        url = self.NIH_CIR_URL.format(inchi_key)
        logger.debug(f"Querying NIH CIR for: {inchi_key}")

        try:
            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 200:
                smiles = response.text.strip()
                # Validate response is actually SMILES
                if (smiles and
                    not smiles.startswith("Error") and
                    not smiles.startswith("<")):
                    if self._validate_smiles(smiles):
                        logger.debug(f"NIH CIR success: {smiles}")
                        return APIResponse(
                            success=True,
                            data=smiles,
                            source="NIH CIR"
                        )
                    else:
                        logger.warning(f"NIH CIR returned invalid SMILES: {smiles}")

            elif response.status_code == 429:
                raise RateLimitError(
                    "NIH CIR rate limit exceeded",
                    api_name="NIH CIR"
                )

            return APIResponse(
                success=False,
                error=f"NIH CIR returned status {response.status_code}"
            )

        except RateLimitError:
            raise
        except requests.exceptions.Timeout:
            logger.warning("NIH CIR request timed out")
            return APIResponse(success=False, error="Request timed out")
        except requests.exceptions.RequestException as e:
            logger.warning(f"NIH CIR request failed: {e}")
            return APIResponse(success=False, error=str(e))

    def _query_pubchem(self, inchi_key: str) -> APIResponse:
        """Query PubChem API.

        Args:
            inchi_key: InChI Key to look up

        Returns:
            APIResponse with result
        """
        url = self.PUBCHEM_URL.format(inchi_key)
        logger.debug(f"Querying PubChem for: {inchi_key}")

        try:
            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    properties = data['PropertyTable']['Properties']
                    if properties and 'IsomericSMILES' in properties[0]:
                        smiles = properties[0]['IsomericSMILES']
                        if self._validate_smiles(smiles):
                            logger.debug(f"PubChem success: {smiles}")
                            return APIResponse(
                                success=True,
                                data=smiles,
                                source="PubChem"
                            )

            elif response.status_code == 429:
                raise RateLimitError(
                    "PubChem rate limit exceeded",
                    api_name="PubChem"
                )

            return APIResponse(
                success=False,
                error=f"PubChem returned status {response.status_code}"
            )

        except RateLimitError:
            raise
        except requests.exceptions.Timeout:
            logger.warning("PubChem request timed out")
            return APIResponse(success=False, error="Request timed out")
        except requests.exceptions.RequestException as e:
            logger.warning(f"PubChem request failed: {e}")
            return APIResponse(success=False, error=str(e))
        except (ValueError, KeyError) as e:
            logger.warning(f"PubChem response parsing failed: {e}")
            return APIResponse(success=False, error=f"Invalid response: {e}")

    def inchi_key_to_smiles(
        self,
        inchi_key: str,
        use_fallback: bool = True
    ) -> APIResponse:
        """Convert InChI Key to SMILES using external APIs.

        Tries NIH CIR first, then falls back to PubChem if enabled.

        Args:
            inchi_key: InChI Key to convert
            use_fallback: Whether to try PubChem if NIH CIR fails

        Returns:
            APIResponse with SMILES if successful

        Raises:
            RateLimitError: If rate limit is exceeded on all APIs
        """
        if not inchi_key or not isinstance(inchi_key, str):
            return APIResponse(success=False, error="Invalid InChI Key")

        inchi_key = inchi_key.strip().upper()

        # Try NIH CIR first
        response = self._query_nih_cir(inchi_key)
        if response.success:
            return response

        # Try PubChem as fallback
        if use_fallback:
            logger.debug("Falling back to PubChem")
            response = self._query_pubchem(inchi_key)
            if response.success:
                return response

        return APIResponse(
            success=False,
            error="Could not convert InChI Key using available APIs"
        )


# Singleton instance for convenience
_api_client: Optional[ChemicalAPIClient] = None


def get_api_client() -> ChemicalAPIClient:
    """Get the shared API client instance.

    Returns:
        ChemicalAPIClient singleton instance
    """
    global _api_client
    if _api_client is None:
        _api_client = ChemicalAPIClient()
    return _api_client
