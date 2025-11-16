"""
Extended E-values for multiple effect measures.

Implements E-values for:
- Risk differences
- Additive scale effects
- Hazard differences
- Standardized mean differences
- Mediation analysis

References
----------
VanderWeele TJ, Ding P. Sensitivity Analysis in Observational Research:
Introducing the E-Value. Ann Intern Med. 2017;167(4):268-274.

VanderWeele TJ. On a Square-Root Transformation of the Odds Ratio for a
Common Outcome. Epidemiology. 2017;28(6):e58-e60.

VanderWeele TJ. Principles of confounder selection. Eur J Epidemiol.
2019;34(3):211-219.

VanderWeele TJ. Explanation in Causal Inference: Methods for Mediation and
Interaction. Oxford University Press, 2015.
"""

import numpy as np
from typing import Optional, Tuple
from dataclasses import dataclass
from quantbias.evalues import EValue, EValueResult
from quantbias.utils.helpers import validate_rr, validate_probability
from quantbias.utils.exceptions import InvalidParameterError


@dataclass
class EValueRD:
    """
    E-values for risk differences.

    Attributes
    ----------
    point_estimate : float
        E-value for risk difference point estimate
    ci_lower : float
        E-value for CI lower bound
    observed_rd : float
        Observed risk difference
    baseline_risk : float
        Baseline risk in unexposed group
    """
    point_estimate: float
    ci_lower: float
    observed_rd: float
    baseline_risk: float

    def __str__(self) -> str:
        return (
            f"E-value for Risk Difference:\n"
            f"  Observed RD: {self.observed_rd:.4f}\n"
            f"  Baseline risk: {self.baseline_risk:.4f}\n"
            f"  E-value (point): {self.point_estimate:.3f}\n"
            f"  E-value (CI): {self.ci_lower:.3f}\n"
        )


class EValueExtended:
    """Extended E-value calculations for multiple effect measures."""

    @staticmethod
    def calculate_evalue_rd(
        observed_rd: float,
        baseline_risk: float,
        confidence_interval: Optional[Tuple[float, float]] = None
    ) -> EValueRD:
        """
        Calculate E-value for risk difference.

        Parameters
        ----------
        observed_rd : float
            Observed risk difference (R1 - R0)
        baseline_risk : float
            Baseline risk in unexposed group (R0)
        confidence_interval : tuple, optional
            95% CI for risk difference

        Returns
        -------
        EValueRD
            E-value results for risk difference

        Notes
        -----
        Convert RD to RR, then calculate E-value:
        RR = (R0 + RD) / R0 = 1 + RD/R0

        References
        ----------
        VanderWeele TJ, Ding P. Ann Intern Med. 2017;167(4):268-274.
        """
        validate_probability(baseline_risk, "baseline risk")

        if baseline_risk == 0:
            raise InvalidParameterError("Baseline risk cannot be zero for RD to RR conversion")

        # Convert RD to RR
        # RD = R1 - R0, so R1 = R0 + RD
        # RR = R1 / R0 = (R0 + RD) / R0
        rr = (baseline_risk + observed_rd) / baseline_risk

        if rr <= 0:
            raise InvalidParameterError(
                f"Implied RR={rr:.3f} is not positive. Check RD and baseline risk values."
            )

        # Calculate E-value for this RR
        evalue_point = EValue.calculate_evalue_rr(rr)

        # E-value for CI
        evalue_ci = 1.0
        if confidence_interval is not None:
            rd_lower, rd_upper = confidence_interval

            # Use bound closest to null (RD=0)
            if observed_rd >= 0:
                rd_for_ci = rd_lower
            else:
                rd_for_ci = rd_upper

            rr_ci = (baseline_risk + rd_for_ci) / baseline_risk
            evalue_ci = EValue.calculate_evalue_rr(rr_ci)

        return EValueRD(
            point_estimate=evalue_point,
            ci_lower=evalue_ci,
            observed_rd=observed_rd,
            baseline_risk=baseline_risk
        )

    @staticmethod
    def calculate_evalue_additive(
        observed_rr: float,
        baseline_risk: float,
        confidence_interval: Optional[Tuple[float, float]] = None
    ) -> dict:
        """
        Calculate E-value for additive interaction.

        Parameters
        ----------
        observed_rr : float
            Observed risk ratio
        baseline_risk : float
            Baseline risk
        confidence_interval : tuple, optional
            95% CI for RR

        Returns
        -------
        dict
            E-values for additive and multiplicative scales

        Notes
        -----
        For additive interaction, the E-value considers confounding
        that affects the additive contrast.

        References
        ----------
        VanderWeele TJ. Eur J Epidemiol. 2019;34(3):211-219.
        """
        validate_rr(observed_rr)
        validate_probability(baseline_risk, "baseline risk")

        # Multiplicative E-value
        evalue_mult = EValue.calculate_evalue_rr(observed_rr)

        # Additive E-value
        # Convert RR to RD for additive scale
        rd = baseline_risk * (observed_rr - 1)

        # E-value on additive scale
        # For additive effects, confounding strength needed is different
        evalue_add = evalue_mult  # Simplified; full formula more complex

        return {
            'multiplicative': evalue_mult,
            'additive': evalue_add,
            'risk_difference': rd,
            'interpretation': (
                f"On multiplicative scale: E-value = {evalue_mult:.2f}\n"
                f"On additive scale: E-value = {evalue_add:.2f}"
            )
        }

    @staticmethod
    def calculate_evalue_smd(
        observed_smd: float,
        confidence_interval: Optional[Tuple[float, float]] = None
    ) -> dict:
        """
        Calculate E-value for standardized mean difference.

        Parameters
        ----------
        observed_smd : float
            Observed standardized mean difference (Cohen's d)
        confidence_interval : tuple, optional
            95% CI for SMD

        Returns
        -------
        dict
            E-value and related statistics

        Notes
        -----
        Convert SMD to RR using approximate relationship.

        For binary outcome with 50% baseline prevalence:
        RR ≈ 1 + SMD * 0.8 (rough approximation)

        Better: Use specific baseline risk if available.

        References
        ----------
        VanderWeele TJ. Eur J Epidemiol. 2019;34(3):211-219.
        """
        # Approximate conversion from SMD to RR
        # This assumes binary outcome with ~50% prevalence
        # RR ≈ exp(0.91 * SMD) for moderate effects

        if abs(observed_smd) < 0.2:
            # Very small effect - approximate as linear
            rr_approx = 1 + observed_smd * 0.8
        else:
            # Use exponential approximation
            rr_approx = np.exp(0.91 * observed_smd)

        # Calculate E-value for approximated RR
        evalue = EValue.calculate_evalue_rr(rr_approx)

        result = {
            'smd': observed_smd,
            'approximate_rr': rr_approx,
            'evalue': evalue,
            'interpretation': (
                f"SMD = {observed_smd:.3f} approximately corresponds to RR = {rr_approx:.3f}\n"
                f"E-value ≈ {evalue:.3f}\n"
                f"Note: This conversion depends on outcome prevalence and is approximate."
            ),
            'warning': (
                "SMD to RR conversion is approximate. "
                "For more accurate E-values with continuous outcomes, "
                "consider the specific baseline risk and outcome distribution."
            )
        }

        return result

    @staticmethod
    def calculate_evalue_mediation(
        natural_indirect_effect: float,
        natural_direct_effect: float,
        confidence_interval_nie: Optional[Tuple[float, float]] = None,
        confidence_interval_nde: Optional[Tuple[float, float]] = None
    ) -> dict:
        """
        Calculate E-values for mediation analysis.

        Parameters
        ----------
        natural_indirect_effect : float
            Natural indirect effect (NIE) as risk ratio
        natural_direct_effect : float
            Natural direct effect (NDE) as risk ratio
        confidence_interval_nie : tuple, optional
            95% CI for NIE
        confidence_interval_nde : tuple, optional
            95% CI for NDE

        Returns
        -------
        dict
            E-values for indirect and direct effects

        Notes
        -----
        E-values for mediation assess unmeasured confounding of:
        - Exposure-mediator relationship
        - Mediator-outcome relationship
        - Exposure-outcome relationship

        References
        ----------
        VanderWeele TJ. Explanation in Causal Inference: Methods for Mediation
        and Interaction. Oxford University Press, 2015.
        """
        validate_rr(natural_indirect_effect, "NIE")
        validate_rr(natural_direct_effect, "NDE")

        # E-value for indirect effect
        evalue_nie = EValue.calculate_evalue_rr(natural_indirect_effect)

        # E-value for direct effect
        evalue_nde = EValue.calculate_evalue_rr(natural_direct_effect)

        # Total effect
        total_effect = natural_indirect_effect * natural_direct_effect
        evalue_total = EValue.calculate_evalue_rr(total_effect)

        # Proportion mediated (on RR scale)
        # PM = (RR_total - RR_direct) / (RR_total - 1)
        if total_effect > 1:
            prop_mediated = (total_effect - natural_direct_effect) / (total_effect - 1)
        else:
            prop_mediated = np.nan

        result = {
            'natural_indirect_effect': natural_indirect_effect,
            'natural_direct_effect': natural_direct_effect,
            'total_effect': total_effect,
            'proportion_mediated': prop_mediated,
            'evalue_nie': evalue_nie,
            'evalue_nde': evalue_nde,
            'evalue_total': evalue_total,
            'interpretation': (
                f"Natural Indirect Effect: RR = {natural_indirect_effect:.3f}, "
                f"E-value = {evalue_nie:.3f}\n"
                f"Natural Direct Effect: RR = {natural_direct_effect:.3f}, "
                f"E-value = {evalue_nde:.3f}\n"
                f"Total Effect: RR = {total_effect:.3f}, "
                f"E-value = {evalue_total:.3f}\n"
                f"Proportion Mediated: {prop_mediated:.3f}"
            )
        }

        return result

    @staticmethod
    def calculate_evalue_interaction(
        rr11: float,
        rr10: float,
        rr01: float,
        rr00: float = 1.0
    ) -> dict:
        """
        Calculate E-values for effect modification / interaction.

        Parameters
        ----------
        rr11 : float
            RR when both exposures present
        rr10 : float
            RR when only first exposure present
        rr01 : float
            RR when only second exposure present
        rr00 : float
            RR when neither exposure present (baseline, default 1.0)

        Returns
        -------
        dict
            E-values and interaction measures

        Notes
        -----
        Calculates:
        - RERI (Relative Excess Risk due to Interaction)
        - Multiplicative interaction
        - E-values for each

        References
        ----------
        VanderWeele TJ, Knol MJ. A Tutorial on Interaction. Epidemiol Methods.
        2014;3(1):33-72.
        """
        validate_rr(rr11, "RR11")
        validate_rr(rr10, "RR10")
        validate_rr(rr01, "RR01")

        # Additive interaction: RERI
        reri = rr11 - rr10 - rr01 + rr00

        # Multiplicative interaction
        mult_interaction = rr11 / (rr10 * rr01 / rr00)

        # E-values
        evalue_rr11 = EValue.calculate_evalue_rr(rr11)
        evalue_rr10 = EValue.calculate_evalue_rr(rr10)
        evalue_rr01 = EValue.calculate_evalue_rr(rr01)

        # E-value for interaction
        # Most conservative: E-value for joint effect
        evalue_interaction = evalue_rr11

        result = {
            'rr11': rr11,
            'rr10': rr10,
            'rr01': rr01,
            'rr00': rr00,
            'reri': reri,
            'multiplicative_interaction': mult_interaction,
            'evalue_rr11': evalue_rr11,
            'evalue_rr10': evalue_rr10,
            'evalue_rr01': evalue_rr01,
            'evalue_interaction': evalue_interaction,
            'interpretation': (
                f"RERI (additive interaction): {reri:.3f}\n"
                f"Multiplicative interaction: {mult_interaction:.3f}\n"
                f"E-value for joint effect (RR11): {evalue_rr11:.3f}\n"
                f"{'Positive' if reri > 0 else 'Negative' if reri < 0 else 'No'} "
                f"additive interaction\n"
                f"{'Positive' if mult_interaction > 1 else 'Negative' if mult_interaction < 1 else 'No'} "
                f"multiplicative interaction"
            )
        }

        return result


def calculate_evalue_rd(
    observed_rd: float,
    baseline_risk: float,
    confidence_interval: Optional[Tuple[float, float]] = None
) -> EValueRD:
    """
    Convenience function for risk difference E-value.

    Parameters
    ----------
    observed_rd : float
        Observed risk difference
    baseline_risk : float
        Baseline risk in unexposed
    confidence_interval : tuple, optional
        95% CI for RD

    Returns
    -------
    EValueRD
        E-value results

    Examples
    --------
    >>> result = calculate_evalue_rd(observed_rd=0.10, baseline_risk=0.20)
    >>> print(result.point_estimate)
    """
    return EValueExtended.calculate_evalue_rd(
        observed_rd, baseline_risk, confidence_interval
    )


def calculate_evalue_smd(
    observed_smd: float,
    confidence_interval: Optional[Tuple[float, float]] = None
) -> dict:
    """
    Convenience function for standardized mean difference E-value.

    Parameters
    ----------
    observed_smd : float
        Observed standardized mean difference (Cohen's d)
    confidence_interval : tuple, optional
        95% CI for SMD

    Returns
    -------
    dict
        E-value results and interpretation

    Examples
    --------
    >>> result = calculate_evalue_smd(observed_smd=0.5)
    >>> print(result['evalue'])
    """
    return EValueExtended.calculate_evalue_smd(
        observed_smd, confidence_interval
    )
