"""Custom exceptions for bias analysis."""


class BiasAnalysisError(Exception):
    """Base exception for bias analysis errors."""
    pass


class InvalidParameterError(BiasAnalysisError):
    """Raised when a parameter value is invalid."""
    pass


class ConvergenceError(BiasAnalysisError):
    """Raised when an iterative method fails to converge."""
    pass


class DataValidationError(BiasAnalysisError):
    """Raised when input data fails validation."""
    pass


class InsufficientDataError(BiasAnalysisError):
    """Raised when insufficient data is provided."""
    pass
