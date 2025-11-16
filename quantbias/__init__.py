"""
QuantBias: Quantitative Bias Analysis for Observational Studies

A comprehensive package for conducting quantitative bias analysis in observational
research and meta-analysis.
"""

__version__ = "0.1.0"
__author__ = "Research Team"
__license__ = "MIT"

from quantbias.evalues import EValue, calculate_evalue
from quantbias.sensitivity import UnmeasuredConfounding, SensitivityAnalysis
from quantbias.bias_parameters import (
    BiasParameter,
    SelectionBias,
    MeasurementError,
    Confounding,
)
from quantbias.monte_carlo import MonteCarloAnalysis
from quantbias.bias_correction import BiasCorrection

__all__ = [
    "EValue",
    "calculate_evalue",
    "UnmeasuredConfounding",
    "SensitivityAnalysis",
    "BiasParameter",
    "SelectionBias",
    "MeasurementError",
    "Confounding",
    "MonteCarloAnalysis",
    "BiasCorrection",
]
