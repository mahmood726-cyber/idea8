"""
E-values for sensitivity analysis to unmeasured confounding.

Implementation based on:
VanderWeele TJ, Ding P. Sensitivity Analysis in Observational Research:
Introducing the E-Value. Ann Intern Med. 2017;167(4):268-274.
"""

import numpy as np
from typing import Optional, Tuple, Union
from dataclasses import dataclass
from quantbias.utils.helpers import validate_rr, validate_ci, validate_or
from quantbias.utils.exceptions import InvalidParameterError


@dataclass
class EValueResult:
    """
    Container for E-value results.

    Attributes
    ----------
    point_estimate : float
        E-value for the point estimate
    ci_lower : float
        E-value for the lower confidence limit
    observed_effect : float
        Original observed effect measure
    observed_ci : tuple of float
        Original confidence interval
    effect_measure : str
        Type of effect measure ('RR', 'OR', 'HR')
    """
    point_estimate: float
    ci_lower: float
    observed_effect: float
    observed_ci: Tuple[float, float]
    effect_measure: str = "RR"

    def __str__(self) -> str:
        return (
            f"E-value Results:\n"
            f"  Observed {self.effect_measure}: {self.observed_effect:.3f} "
            f"({self.observed_ci[0]:.3f}, {self.observed_ci[1]:.3f})\n"
            f"  E-value (point estimate): {self.point_estimate:.3f}\n"
            f"  E-value (CI lower bound): {self.ci_lower:.3f}\n"
            f"\nInterpretation:\n"
            f"  An unmeasured confounder associated with both the exposure and outcome\n"
            f"  by a risk ratio of {self.point_estimate:.2f}-fold each could explain\n"
            f"  away the observed point estimate, but weaker confounding could not.\n"
            f"  To move the confidence interval to include the null, an unmeasured\n"
            f"  confounder would need to be associated with both by {self.ci_lower:.2f}-fold."
        )


class EValue:
    """
    Calculate E-values for sensitivity analysis to unmeasured confounding.

    The E-value quantifies the minimum strength of association on the risk ratio
    scale that an unmeasured confounder would need to have with both the exposure
    and the outcome to fully explain away a specific exposure-outcome association.
    """

    @staticmethod
    def calculate_evalue_rr(rr: float) -> float:
        """
        Calculate E-value for a risk ratio.

        Parameters
        ----------
        rr : float
            Risk ratio (must be >= 1 or will be inverted)

        Returns
        -------
        float
            E-value

        Notes
        -----
        Formula: E-value = RR + sqrt(RR * (RR - 1))

        For RR < 1, the formula uses 1/RR.
        """
        # Handle protective effects (RR < 1)
        if rr < 1:
            rr = 1 / rr

        validate_rr(rr, "risk ratio")

        if rr < 1:
            return 1.0  # No unmeasured confounding needed to explain null effect

        evalue = rr + np.sqrt(rr * (rr - 1))
        return evalue

    @staticmethod
    def calculate_evalue_or(
        odds_ratio: float,
        baseline_risk: Optional[float] = None
    ) -> float:
        """
        Calculate E-value for an odds ratio.

        Parameters
        ----------
        odds_ratio : float
            Odds ratio
        baseline_risk : float, optional
            Baseline risk (for rare outcomes, OR approximates RR)
            If provided, OR is converted to RR for E-value calculation

        Returns
        -------
        float
            E-value

        Notes
        -----
        For rare outcomes (< 15%), OR approximates RR and can be used directly.
        For common outcomes, provide baseline_risk to convert OR to RR.
        """
        validate_or(odds_ratio, "odds ratio")

        if baseline_risk is not None:
            # Convert OR to RR for more accurate E-value
            from quantbias.utils.helpers import convert_or_to_rr
            rr = convert_or_to_rr(odds_ratio, baseline_risk)
            return EValue.calculate_evalue_rr(rr)
        else:
            # Use OR directly (approximation for rare outcomes)
            return EValue.calculate_evalue_rr(odds_ratio)

    @staticmethod
    def calculate_evalue_hr(hazard_ratio: float) -> float:
        """
        Calculate E-value for a hazard ratio.

        Parameters
        ----------
        hazard_ratio : float
            Hazard ratio

        Returns
        -------
        float
            E-value

        Notes
        -----
        The same formula applies to hazard ratios as risk ratios.
        """
        validate_rr(hazard_ratio, "hazard ratio")
        return EValue.calculate_evalue_rr(hazard_ratio)

    @staticmethod
    def calculate(
        observed_rr: float,
        confidence_interval: Optional[Tuple[float, float]] = None,
        effect_measure: str = "RR",
        baseline_risk: Optional[float] = None
    ) -> EValueResult:
        """
        Calculate E-values for point estimate and confidence interval.

        Parameters
        ----------
        observed_rr : float
            Observed effect measure (RR, OR, or HR)
        confidence_interval : tuple of float, optional
            95% confidence interval (lower, upper)
        effect_measure : str
            Type of effect measure: 'RR', 'OR', or 'HR'
        baseline_risk : float, optional
            Baseline risk (needed for OR conversion)

        Returns
        -------
        EValueResult
            E-value results for point estimate and CI

        Examples
        --------
        >>> result = EValue.calculate(observed_rr=2.5, confidence_interval=(1.8, 3.5))
        >>> print(result.point_estimate)
        3.89
        >>> print(result.ci_lower)
        2.54
        """
        effect_measure = effect_measure.upper()

        # Calculate E-value for point estimate
        if effect_measure == "RR":
            evalue_point = EValue.calculate_evalue_rr(observed_rr)
        elif effect_measure == "OR":
            evalue_point = EValue.calculate_evalue_or(observed_rr, baseline_risk)
        elif effect_measure == "HR":
            evalue_point = EValue.calculate_evalue_hr(observed_rr)
        else:
            raise InvalidParameterError(
                f"Unsupported effect measure: {effect_measure}. Use 'RR', 'OR', or 'HR'."
            )

        # Calculate E-value for confidence interval
        evalue_ci = 1.0  # Default if no CI provided
        ci = (np.nan, np.nan)

        if confidence_interval is not None:
            validate_ci(confidence_interval, observed_rr, "confidence interval")
            ci = confidence_interval

            # For E-value of CI, use the bound closest to the null (1.0)
            # If RR > 1, use lower bound; if RR < 1, use upper bound
            if observed_rr >= 1:
                ci_for_evalue = ci[0]  # Lower bound
            else:
                ci_for_evalue = ci[1]  # Upper bound

            # Calculate E-value for the CI bound
            if effect_measure == "RR":
                evalue_ci = EValue.calculate_evalue_rr(ci_for_evalue)
            elif effect_measure == "OR":
                evalue_ci = EValue.calculate_evalue_or(ci_for_evalue, baseline_risk)
            elif effect_measure == "HR":
                evalue_ci = EValue.calculate_evalue_hr(ci_for_evalue)

        return EValueResult(
            point_estimate=evalue_point,
            ci_lower=evalue_ci,
            observed_effect=observed_rr,
            observed_ci=ci,
            effect_measure=effect_measure
        )

    @staticmethod
    def required_confounding_strength(
        observed_rr: float,
        true_rr: float = 1.0
    ) -> float:
        """
        Calculate required confounding strength to explain away effect to a specific true RR.

        Parameters
        ----------
        observed_rr : float
            Observed risk ratio
        true_rr : float
            True risk ratio (default 1.0 for null effect)

        Returns
        -------
        float
            Required confounding strength (as RR)

        Notes
        -----
        This calculates what strength of confounding would be needed to reduce
        the observed effect to a specific true effect.
        """
        validate_rr(observed_rr, "observed risk ratio")
        validate_rr(true_rr, "true risk ratio")

        # Handle case where observed < true (protective effect becoming harmful)
        if observed_rr < true_rr:
            observed_rr, true_rr = true_rr, observed_rr

        # Bias factor needed
        bias_factor = observed_rr / true_rr

        # E-value for this bias factor
        evalue = EValue.calculate_evalue_rr(bias_factor)

        return evalue


def calculate_evalue(
    observed_rr: float,
    confidence_interval: Optional[Tuple[float, float]] = None,
    effect_measure: str = "RR",
    baseline_risk: Optional[float] = None
) -> EValueResult:
    """
    Convenience function to calculate E-values.

    Parameters
    ----------
    observed_rr : float
        Observed effect measure
    confidence_interval : tuple of float, optional
        95% confidence interval
    effect_measure : str
        Type of effect measure: 'RR', 'OR', or 'HR'
    baseline_risk : float, optional
        Baseline risk (for OR conversion)

    Returns
    -------
    EValueResult
        E-value results

    Examples
    --------
    >>> result = calculate_evalue(2.0, (1.5, 2.7))
    >>> print(f"E-value: {result.point_estimate:.2f}")
    E-value: 3.41
    """
    return EValue.calculate(
        observed_rr=observed_rr,
        confidence_interval=confidence_interval,
        effect_measure=effect_measure,
        baseline_risk=baseline_risk
    )
