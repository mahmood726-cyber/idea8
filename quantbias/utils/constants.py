"""Constants and enumerations for bias analysis."""

from enum import Enum


class EffectMeasure(Enum):
    """Supported effect measures."""
    RISK_RATIO = "rr"
    ODDS_RATIO = "or"
    HAZARD_RATIO = "hr"
    RISK_DIFFERENCE = "rd"
    STANDARDIZED_MEAN_DIFFERENCE = "smd"


class BiasType(Enum):
    """Types of bias."""
    SELECTION = "selection"
    INFORMATION = "information"
    CONFOUNDING = "confounding"
    UNMEASURED_CONFOUNDING = "unmeasured_confounding"
    MEASUREMENT_ERROR = "measurement_error"


class CIMethod(Enum):
    """Confidence interval calculation methods."""
    NORMAL = "normal"
    BOOTSTRAP = "bootstrap"
    PERCENTILE = "percentile"


# Default values
DEFAULT_ALPHA = 0.05
DEFAULT_N_BOOTSTRAP = 10000
DEFAULT_RANDOM_SEED = 42

# Numerical constants
EPSILON = 1e-10  # Small value for numerical stability
MAX_ITERATIONS = 1000  # Maximum iterations for iterative methods
CONVERGENCE_TOL = 1e-6  # Convergence tolerance
