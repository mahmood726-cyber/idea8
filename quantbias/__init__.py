"""
QuantBias: Quantitative Bias Analysis for Observational Studies

A comprehensive package for conducting quantitative bias analysis in observational
research and meta-analysis.

All formulas validated against published examples (VanderWeele & Ding 2017,
Lash et al. 2021, Greenland & Kleinbaum 1983).
"""

__version__ = "0.2.0"
__author__ = "Research Team"
__license__ = "MIT"

# Core E-values
from quantbias.evalues import EValue, calculate_evalue

# Extended E-values (NEW in v0.2.0)
from quantbias.evalues_extended import (
    EValueExtended,
    calculate_evalue_rd,
    calculate_evalue_smd,
)

# Sensitivity analysis
from quantbias.sensitivity import UnmeasuredConfounding, SensitivityAnalysis

# Bias parameters (CORRECTED in v0.2.0)
from quantbias.bias_parameters import (
    BiasParameter,
    SelectionBias,
    MeasurementError,
    Confounding,
)

# Monte Carlo analysis
from quantbias.monte_carlo import MonteCarloAnalysis, MultipleBiasModeling

# Bias correction
from quantbias.bias_correction import BiasCorrection

# Statistical inference (NEW in v0.2.0)
from quantbias.inference import BiasInference, CorrectedEstimate

__all__ = [
    # E-values
    "EValue",
    "calculate_evalue",
    "EValueExtended",
    "calculate_evalue_rd",
    "calculate_evalue_smd",
    # Sensitivity analysis
    "UnmeasuredConfounding",
    "SensitivityAnalysis",
    # Bias parameters
    "BiasParameter",
    "SelectionBias",
    "MeasurementError",
    "Confounding",
    # Monte Carlo
    "MonteCarloAnalysis",
    "MultipleBiasModeling",
    # Bias correction
    "BiasCorrection",
    # Statistical inference
    "BiasInference",
    "CorrectedEstimate",
]
