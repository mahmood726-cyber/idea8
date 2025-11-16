"""
Comprehensive bias correction framework.

Integrates all bias analysis methods for complete bias correction workflow.
"""

import numpy as np
import pandas as pd
from typing import List, Optional, Tuple, Dict
from dataclasses import dataclass

from quantbias.bias_parameters import BiasParameter
from quantbias.monte_carlo import MonteCarloAnalysis, MultipleBiasModeling
from quantbias.sensitivity import UnmeasuredConfounding
from quantbias.evalues import EValue
from quantbias.utils.helpers import validate_rr, validate_ci


@dataclass
class CorrectedEffect:
    """
    Container for bias-corrected effect estimates.

    Attributes
    ----------
    observed_rr : float
        Original observed RR
    corrected_rr : float
        Bias-corrected RR
    corrected_ci : tuple of float
        Corrected confidence interval
    bias_factors : dict
        Individual bias factors applied
    total_bias_factor : float
        Combined bias factor
    evalue : float
        E-value for corrected estimate
    """
    observed_rr: float
    corrected_rr: float
    corrected_ci: Tuple[float, float]
    bias_factors: Dict[str, float]
    total_bias_factor: float
    evalue: float

    def __str__(self) -> str:
        return (
            f"Bias Correction Results:\n"
            f"  Observed RR: {self.observed_rr:.3f}\n"
            f"  Corrected RR: {self.corrected_rr:.3f}\n"
            f"  Corrected 95% CI: ({self.corrected_ci[0]:.3f}, {self.corrected_ci[1]:.3f})\n"
            f"  Total bias factor: {self.total_bias_factor:.3f}\n"
            f"  E-value (corrected): {self.evalue:.3f}\n"
        )


class BiasCorrection:
    """
    Comprehensive bias correction workflow.

    Provides unified interface for:
    - Deterministic bias correction
    - Probabilistic bias analysis
    - Sensitivity analysis
    - E-value calculation
    """

    def __init__(
        self,
        observed_rr: float,
        observed_ci: Optional[Tuple[float, float]] = None,
        biases: Optional[List[BiasParameter]] = None
    ):
        """
        Initialize bias correction.

        Parameters
        ----------
        observed_rr : float
            Observed risk ratio
        observed_ci : tuple of float, optional
            Observed confidence interval
        biases : list of BiasParameter, optional
            List of bias parameters to correct for
        """
        validate_rr(observed_rr)
        self.observed_rr = observed_rr
        self.observed_ci = observed_ci

        if observed_ci is not None:
            validate_ci(observed_ci, observed_rr)

        self.biases = biases or []

    def correct(
        self,
        method: str = "deterministic"
    ) -> CorrectedEffect:
        """
        Apply bias correction.

        Parameters
        ----------
        method : str
            Correction method: 'deterministic' or 'probabilistic'

        Returns
        -------
        CorrectedEffect
            Bias-corrected results
        """
        if method == "deterministic":
            return self._deterministic_correction()
        elif method == "probabilistic":
            return self._probabilistic_correction()
        else:
            raise ValueError(f"Unknown method: {method}")

    def _deterministic_correction(self) -> CorrectedEffect:
        """
        Apply deterministic bias correction.

        Uses point estimates of bias parameters.
        """
        corrected_rr = self.observed_rr
        bias_factors = {}

        # Apply each bias correction sequentially
        for bias in self.biases:
            previous_rr = corrected_rr
            corrected_rr = bias.apply_bias(corrected_rr)
            bias_factor = previous_rr / corrected_rr
            bias_factors[bias.name] = bias_factor

        total_bias_factor = self.observed_rr / corrected_rr

        # Correct confidence interval
        if self.observed_ci is not None:
            corrected_ci_lower = self.observed_ci[0] / total_bias_factor
            corrected_ci_upper = self.observed_ci[1] / total_bias_factor
            corrected_ci = (corrected_ci_lower, corrected_ci_upper)
        else:
            corrected_ci = (np.nan, np.nan)

        # Calculate E-value for corrected estimate
        evalue = EValue.calculate_evalue_rr(corrected_rr)

        return CorrectedEffect(
            observed_rr=self.observed_rr,
            corrected_rr=corrected_rr,
            corrected_ci=corrected_ci,
            bias_factors=bias_factors,
            total_bias_factor=total_bias_factor,
            evalue=evalue
        )

    def _probabilistic_correction(self, n_iterations: int = 10000) -> CorrectedEffect:
        """
        Apply probabilistic bias correction.

        Uses Monte Carlo simulation with parameter uncertainty.
        """
        mc = MonteCarloAnalysis(
            observed_rr=self.observed_rr,
            bias_parameters=self.biases,
            n_iterations=n_iterations
        )

        results = mc.run()

        # Use median as point estimate
        corrected_rr = results.median
        corrected_ci = (results.ci_lower, results.ci_upper)

        # Calculate bias factors (approximate)
        bias_factors = {
            bias.name: self.observed_rr / corrected_rr
            for bias in self.biases
        }

        total_bias_factor = self.observed_rr / corrected_rr
        evalue = EValue.calculate_evalue_rr(corrected_rr)

        return CorrectedEffect(
            observed_rr=self.observed_rr,
            corrected_rr=corrected_rr,
            corrected_ci=corrected_ci,
            bias_factors=bias_factors,
            total_bias_factor=total_bias_factor,
            evalue=evalue
        )

    def sensitivity_analysis(
        self,
        parameter_ranges: Optional[Dict[str, Tuple[float, float]]] = None
    ) -> pd.DataFrame:
        """
        Perform sensitivity analysis across parameter ranges.

        Parameters
        ----------
        parameter_ranges : dict, optional
            Ranges for each bias parameter

        Returns
        -------
        pd.DataFrame
            Sensitivity analysis results
        """
        if not self.biases:
            raise ValueError("No bias parameters specified")

        # Default ranges if not provided
        if parameter_ranges is None:
            parameter_ranges = {
                'rr_confounder': (1.0, 5.0),
                'prevalence_confounder': (0.1, 0.9),
            }

        # Create UnmeasuredConfounding analysis
        uc = UnmeasuredConfounding(self.observed_rr, self.observed_ci)

        # Generate sensitivity grid
        rr_range = parameter_ranges.get('rr_confounder', (1.0, 5.0))
        prev_range = parameter_ranges.get('prevalence_confounder', (0.1, 0.9))

        rr_values = np.linspace(rr_range[0], rr_range[1], 10)
        prev_values = np.linspace(prev_range[0], prev_range[1], 10)

        return uc.sensitivity_grid(
            rr_confounder_outcome=rr_values,
            prevalence_exposed=prev_values
        )

    def evalue_analysis(self) -> Dict[str, float]:
        """
        Calculate E-values for observed and corrected estimates.

        Returns
        -------
        dict
            E-values for various estimates
        """
        evalue_observed = EValue.calculate(
            observed_rr=self.observed_rr,
            confidence_interval=self.observed_ci
        )

        corrected = self.correct()

        evalue_corrected = EValue.calculate_evalue_rr(corrected.corrected_rr)

        return {
            'evalue_observed_point': evalue_observed.point_estimate,
            'evalue_observed_ci': evalue_observed.ci_lower,
            'evalue_corrected': evalue_corrected,
            'interpretation': (
                f"Original estimate requires E-value of {evalue_observed.point_estimate:.2f}. "
                f"After bias correction, requires E-value of {evalue_corrected:.2f}."
            )
        }

    def full_report(self) -> str:
        """
        Generate comprehensive bias analysis report.

        Returns
        -------
        str
            Full analysis report
        """
        report = []
        report.append("=" * 60)
        report.append("COMPREHENSIVE BIAS ANALYSIS REPORT")
        report.append("=" * 60)
        report.append("")

        # Original estimate
        report.append("1. OBSERVED EFFECT")
        report.append(f"   Risk Ratio: {self.observed_rr:.3f}")
        if self.observed_ci:
            report.append(f"   95% CI: ({self.observed_ci[0]:.3f}, {self.observed_ci[1]:.3f})")
        report.append("")

        # E-values
        report.append("2. E-VALUE ANALYSIS")
        evalue_results = self.evalue_analysis()
        report.append(f"   E-value (point): {evalue_results['evalue_observed_point']:.3f}")
        report.append(f"   E-value (CI): {evalue_results['evalue_observed_ci']:.3f}")
        report.append("")

        # Bias correction
        if self.biases:
            report.append("3. BIAS CORRECTION")
            report.append(f"   Number of bias sources: {len(self.biases)}")
            for bias in self.biases:
                report.append(f"   - {bias.name}")

            corrected = self.correct()
            report.append("")
            report.append(f"   Corrected RR: {corrected.corrected_rr:.3f}")
            report.append(f"   Corrected CI: ({corrected.corrected_ci[0]:.3f}, {corrected.corrected_ci[1]:.3f})")
            report.append(f"   Total bias factor: {corrected.total_bias_factor:.3f}")
            report.append("")

            report.append("4. BIAS DECOMPOSITION")
            for bias_name, factor in corrected.bias_factors.items():
                report.append(f"   {bias_name}: {factor:.3f}")
            report.append("")

        # Interpretation
        report.append("5. INTERPRETATION")
        report.append(evalue_results['interpretation'])
        report.append("")

        report.append("=" * 60)

        return "\n".join(report)
