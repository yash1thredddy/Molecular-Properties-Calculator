"""External API client for chemical database lookups.

This module provides a service for interacting with external chemical
databases like NIH CIR and PubChem for InChI Key to SMILES conversion.
"""

import logging
import time
from typing import Optional, Callable
from dataclasses import dataclass
from urllib.parse import quote

import requests

from molecular_calculator.config.settings import config
from molecular_calculator.utils.exceptions import APIError, RateLimitError
from molecular_calculator.utils.sanitizer import sanitize_inchi_key

# Import rate limiter with fallback
try:
    from molecular_calculator.utils.rate_limiter import nih_limiter, pubchem_limiter
    RATE_LIMITER_AVAILABLE = True
except ImportError:
    RATE_LIMITER_AVAILABLE = False

# Import cache with fallback
try:
    from molecular_calculator.utils.cache import api_cache
    CACHE_AVAILABLE = True
except ImportError:
    CACHE_AVAILABLE = False

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
                     Bounded between 1 and 120 seconds.
        """
        default_timeout = config.API_TIMEOUT_SECONDS
        if timeout is None:
            self.timeout = default_timeout
        else:
            # Bound timeout between 1 and 120 seconds for security
            self.timeout = min(max(timeout, 1), 120)

    def _retry_with_backoff(
        self,
        query_func: Callable,
        *args,
        max_retries: int = None,
        **kwargs
    ) -> APIResponse:
        """Retry an API query with exponential backoff.

        Args:
            query_func: The query function to retry
            *args: Positional arguments for the query function
            max_retries: Maximum number of retry attempts (uses config default if None)
            **kwargs: Keyword arguments for the query function

        Returns:
            APIResponse from the query function
        """
        max_retries = max_retries or config.MAX_RETRY_ATTEMPTS

        for attempt in range(max_retries):
            try:
                response = query_func(*args, **kwargs)

                # Return on success or non-retryable failures
                if response.success or attempt == max_retries - 1:
                    return response

                # Retry on transient failures (timeouts, connection errors)
                if response.error and any(x in response.error.lower() for x in ['timeout', 'connection', 'network']):
                    wait_time = 2 ** attempt  # Exponential backoff: 1s, 2s, 4s
                    logger.info(f"Retry attempt {attempt + 1}/{max_retries} after {wait_time}s")
                    time.sleep(wait_time)
                    continue

                # Don't retry other errors
                return response

            except RateLimitError:
                # Don't retry rate limit errors
                raise
            except requests.exceptions.Timeout:
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    logger.info(f"Request timed out, retry {attempt + 1}/{max_retries} after {wait_time}s")
                    time.sleep(wait_time)
                else:
                    return APIResponse(success=False, error="Request timed out after retries")
            except requests.exceptions.RequestException as e:
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    logger.info(f"Request failed: {e}, retry {attempt + 1}/{max_retries} after {wait_time}s")
                    time.sleep(wait_time)
                else:
                    return APIResponse(success=False, error=f"Request failed after retries: {e}")

        return APIResponse(success=False, error="Max retries exceeded")

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
        # Validate InChI Key format before URL construction (defense in depth)
        if not inchi_key or sanitize_inchi_key(inchi_key) is None:
            return APIResponse(success=False, error="Invalid InChI Key format")

        # URL-encode to prevent SSRF attacks
        url = self.NIH_CIR_URL.format(quote(inchi_key, safe=''))
        logger.debug(f"Querying NIH CIR for: {inchi_key}")

        # Apply rate limiting if available
        if RATE_LIMITER_AVAILABLE:
            if not nih_limiter.acquire():
                raise RateLimitError(
                    "NIH CIR rate limit reached, please wait",
                    api_name="NIH CIR"
                )

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
        # Validate InChI Key format before URL construction (defense in depth)
        if not inchi_key or sanitize_inchi_key(inchi_key) is None:
            return APIResponse(success=False, error="Invalid InChI Key format")

        # URL-encode to prevent SSRF attacks
        url = self.PUBCHEM_URL.format(quote(inchi_key, safe=''))
        logger.debug(f"Querying PubChem for: {inchi_key}")

        # Apply rate limiting if available
        if RATE_LIMITER_AVAILABLE:
            if not pubchem_limiter.acquire():
                raise RateLimitError(
                    "PubChem rate limit reached, please wait",
                    api_name="PubChem"
                )

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
        Results are cached to avoid redundant API calls.

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

        # Sanitize and validate InChI Key before API calls
        sanitized_key = sanitize_inchi_key(inchi_key)
        if sanitized_key is None:
            return APIResponse(
                success=False,
                error="Invalid InChI Key format. Expected format: XXXXXXXXXXXXXX-XXXXXXXXXX-X"
            )

        inchi_key = sanitized_key

        # Check cache first
        if CACHE_AVAILABLE:
            cached = api_cache.get(f"inchi2smiles:{inchi_key}")
            if cached is not None:
                logger.debug(f"Cache hit for InChI Key: {inchi_key}")
                return APIResponse(success=True, data=cached, source="cache")

        # Try NIH CIR first with retry logic
        response = self._retry_with_backoff(self._query_nih_cir, inchi_key)
        if response.success:
            # Cache successful result
            if CACHE_AVAILABLE:
                api_cache.set(f"inchi2smiles:{inchi_key}", response.data)
            return response

        # Try PubChem as fallback with retry logic
        if use_fallback:
            logger.debug("Falling back to PubChem")
            response = self._retry_with_backoff(self._query_pubchem, inchi_key)
            if response.success:
                # Cache successful result
                if CACHE_AVAILABLE:
                    api_cache.set(f"inchi2smiles:{inchi_key}", response.data)
                return response

        return APIResponse(
            success=False,
            error="Could not convert InChI Key using available APIs"
        )


# Singleton instance for convenience
_api_client: Optional[ChemicalAPIClient] = None
_api_client_lock = __import__('threading').Lock()


def get_api_client() -> ChemicalAPIClient:
    """Get the shared API client instance.

    Thread-safe singleton accessor using double-checked locking.

    Returns:
        ChemicalAPIClient singleton instance
    """
    global _api_client
    if _api_client is None:
        with _api_client_lock:
            # Double-check after acquiring lock
            if _api_client is None:
                _api_client = ChemicalAPIClient()
    return _api_client
