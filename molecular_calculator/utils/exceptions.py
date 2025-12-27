"""Custom exceptions for the application.

This module defines a hierarchy of exceptions for better error handling
and more informative error messages.

Exception Hierarchy:
    MolecularCalculatorError (base)
    ├── InvalidSMILESError
    ├── PropertyCalculationError
    ├── FileValidationError
    ├── ConversionError
    └── APIError
        └── RateLimitError
"""


class MolecularCalculatorError(Exception):
    """Base exception for molecular calculator.

    All custom exceptions in this application should inherit from this class.
    This allows catching all application-specific errors with a single except clause.

    Example:
        >>> try:
        ...     risky_operation()
        ... except MolecularCalculatorError as e:
        ...     handle_error(e)
    """

    def __init__(self, message: str = "An error occurred in the molecular calculator"):
        self.message = message
        super().__init__(self.message)


class InvalidSMILESError(MolecularCalculatorError):
    """Raised when a SMILES string is invalid.

    This exception includes the invalid SMILES string for debugging purposes.

    Attributes:
        smiles: The invalid SMILES string
        message: Human-readable error message

    Example:
        >>> raise InvalidSMILESError("invalid_smiles", "Could not parse SMILES")
    """

    def __init__(self, smiles: str, message: str = None):
        self.smiles = smiles
        self.message = message or f"Invalid SMILES: {smiles}"
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return f"InvalidSMILESError(smiles={self.smiles!r}, message={self.message!r})"


class PropertyCalculationError(MolecularCalculatorError):
    """Raised when a molecular property calculation fails.

    This can occur when RDKit encounters an error calculating a specific
    descriptor for a given molecule.

    Attributes:
        property_name: Name of the property that failed (optional)
        smiles: SMILES string of the molecule (optional)

    Example:
        >>> raise PropertyCalculationError(
        ...     "Failed to calculate LogP",
        ...     property_name="LogP",
        ...     smiles="CCO"
        ... )
    """

    def __init__(
        self,
        message: str = "Property calculation failed",
        property_name: str = None,
        smiles: str = None
    ):
        self.property_name = property_name
        self.smiles = smiles
        super().__init__(message)


class FileValidationError(MolecularCalculatorError):
    """Raised when file validation fails.

    This can occur for various reasons:
    - File is too large
    - File type is not allowed
    - File is corrupted or unreadable
    - File is empty

    Attributes:
        filename: Name of the file that failed validation (optional)

    Example:
        >>> raise FileValidationError("File exceeds maximum size of 50MB")
    """

    def __init__(self, message: str = "File validation failed", filename: str = None):
        self.filename = filename
        super().__init__(message)


class ConversionError(MolecularCalculatorError):
    """Raised when chemical format conversion fails.

    This can occur when converting between formats like:
    - InChI to SMILES
    - InChI Key to SMILES
    - Name to SMILES

    Attributes:
        input_value: The value that failed to convert
        source_format: The source format (e.g., 'inchi_key')
        target_format: The target format (e.g., 'smiles')

    Example:
        >>> raise ConversionError(
        ...     "Could not convert InChI Key",
        ...     input_value="INVALID-KEY",
        ...     source_format="inchi_key",
        ...     target_format="smiles"
        ... )
    """

    def __init__(
        self,
        message: str = "Format conversion failed",
        input_value: str = None,
        source_format: str = None,
        target_format: str = None
    ):
        self.input_value = input_value
        self.source_format = source_format
        self.target_format = target_format
        super().__init__(message)


class APIError(MolecularCalculatorError):
    """Raised when an external API call fails.

    This is the base class for API-related errors.

    Attributes:
        api_name: Name of the API that failed (e.g., 'NIH CIR', 'PubChem')
        status_code: HTTP status code (if applicable)

    Example:
        >>> raise APIError("API request failed", api_name="PubChem", status_code=500)
    """

    def __init__(
        self,
        message: str = "API call failed",
        api_name: str = None,
        status_code: int = None
    ):
        self.api_name = api_name
        self.status_code = status_code
        super().__init__(message)


class RateLimitError(APIError):
    """Raised when API rate limit is exceeded.

    This exception indicates that too many requests have been made
    to an external API in a short time period.

    Attributes:
        retry_after: Seconds to wait before retrying (if known)

    Example:
        >>> raise RateLimitError(
        ...     "Rate limit exceeded for NIH CIR API",
        ...     api_name="NIH CIR",
        ...     retry_after=60
        ... )
    """

    def __init__(
        self,
        message: str = "Rate limit exceeded",
        api_name: str = None,
        retry_after: int = None
    ):
        self.retry_after = retry_after
        super().__init__(message, api_name=api_name)
