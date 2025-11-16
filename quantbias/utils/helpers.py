"""Helper functions for bias analysis."""

import numpy as np
from typing import Tuple, Union, Optional
from quantbias.utils.exceptions import InvalidParameterError
from quantbias.utils.constants import EPSILON


def validate_rr(rr: float, parameter_name: str = "risk ratio") -> None:
    """
    Validate risk ratio value.

    Parameters
    ----------
    rr : float
        Risk ratio value
    parameter_name : str
        Name of parameter for error message

    Raises
    ------
    InvalidParameterError
        If RR is not positive
    """
    if rr <= 0:
        raise InvalidParameterError(f"{parameter_name} must be positive, got {rr}")


def validate_or(odds_ratio: float, parameter_name: str = "odds ratio") -> None:
    """
    Validate odds ratio value.

    Parameters
    ----------
    odds_ratio : float
        Odds ratio value
    parameter_name : str
        Name of parameter for error message

    Raises
    ------
    InvalidParameterError
        If OR is not positive
    """
    if odds_ratio <= 0:
        raise InvalidParameterError(f"{parameter_name} must be positive, got {odds_ratio}")


def validate_probability(p: float, parameter_name: str = "probability") -> None:
    """
    Validate probability value.

    Parameters
    ----------
    p : float
        Probability value
    parameter_name : str
        Name of parameter for error message

    Raises
    ------
    InvalidParameterError
        If probability is not in [0, 1]
    """
    if not 0 <= p <= 1:
        raise InvalidParameterError(
            f"{parameter_name} must be between 0 and 1, got {p}"
        )


def validate_ci(
    ci: Tuple[float, float],
    point_estimate: Optional[float] = None,
    parameter_name: str = "confidence interval"
) -> None:
    """
    Validate confidence interval.

    Parameters
    ----------
    ci : tuple of float
        Confidence interval (lower, upper)
    point_estimate : float, optional
        Point estimate (should be within CI)
    parameter_name : str
        Name of parameter for error message

    Raises
    ------
    InvalidParameterError
        If CI is invalid
    """
    lower, upper = ci

    if lower >= upper:
        raise InvalidParameterError(
            f"{parameter_name} lower bound must be less than upper bound"
        )

    if lower <= 0 or upper <= 0:
        raise InvalidParameterError(
            f"{parameter_name} bounds must be positive"
        )

    if point_estimate is not None:
        if not (lower <= point_estimate <= upper):
            raise InvalidParameterError(
                f"Point estimate {point_estimate} not within {parameter_name} ({lower}, {upper})"
            )


def convert_or_to_rr(
    odds_ratio: float,
    baseline_risk: float
) -> float:
    """
    Convert odds ratio to risk ratio.

    Uses the relationship: RR = OR / (1 - p0 + p0 * OR)
    where p0 is the baseline risk.

    Parameters
    ----------
    odds_ratio : float
        Odds ratio
    baseline_risk : float
        Baseline risk (probability in unexposed group)

    Returns
    -------
    float
        Risk ratio

    References
    ----------
    Zhang J, Yu KF. What's the relative risk? A method of correcting the odds ratio
    in cohort studies of common outcomes. JAMA. 1998;280(19):1690-1691.
    """
    validate_or(odds_ratio)
    validate_probability(baseline_risk, "baseline risk")

    rr = odds_ratio / (1 - baseline_risk + baseline_risk * odds_ratio)
    return rr


def convert_rr_to_or(
    risk_ratio: float,
    baseline_risk: float
) -> float:
    """
    Convert risk ratio to odds ratio.

    Parameters
    ----------
    risk_ratio : float
        Risk ratio
    baseline_risk : float
        Baseline risk (probability in unexposed group)

    Returns
    -------
    float
        Odds ratio
    """
    validate_rr(risk_ratio)
    validate_probability(baseline_risk, "baseline risk")

    # Calculate risk in exposed: R1 = RR * R0
    risk_exposed = risk_ratio * baseline_risk

    if risk_exposed >= 1:
        raise InvalidParameterError(
            f"Exposed risk ({risk_exposed}) cannot be >= 1"
        )

    # OR = (R1 / (1 - R1)) / (R0 / (1 - R0))
    odds_exposed = risk_exposed / (1 - risk_exposed + EPSILON)
    odds_unexposed = baseline_risk / (1 - baseline_risk + EPSILON)

    odds_ratio = odds_exposed / (odds_unexposed + EPSILON)
    return odds_ratio


def compute_standard_error(
    lower: float,
    upper: float,
    alpha: float = 0.05
) -> float:
    """
    Compute standard error from confidence interval.

    Assumes log-normal distribution for ratio measures.

    Parameters
    ----------
    lower : float
        Lower confidence bound
    upper : float
        Upper confidence bound
    alpha : float
        Significance level (default 0.05 for 95% CI)

    Returns
    -------
    float
        Standard error on log scale
    """
    from scipy import stats

    z_score = stats.norm.ppf(1 - alpha / 2)
    se = (np.log(upper) - np.log(lower)) / (2 * z_score)
    return se


def compute_ci_from_se(
    estimate: float,
    se: float,
    alpha: float = 0.05
) -> Tuple[float, float]:
    """
    Compute confidence interval from standard error.

    Assumes log-normal distribution for ratio measures.

    Parameters
    ----------
    estimate : float
        Point estimate
    se : float
        Standard error on log scale
    alpha : float
        Significance level (default 0.05 for 95% CI)

    Returns
    -------
    tuple of float
        Lower and upper confidence bounds
    """
    from scipy import stats

    z_score = stats.norm.ppf(1 - alpha / 2)
    log_estimate = np.log(estimate)

    lower = np.exp(log_estimate - z_score * se)
    upper = np.exp(log_estimate + z_score * se)

    return (lower, upper)


def format_ci(
    estimate: float,
    ci: Tuple[float, float],
    decimals: int = 2
) -> str:
    """
    Format estimate and confidence interval for display.

    Parameters
    ----------
    estimate : float
        Point estimate
    ci : tuple of float
        Confidence interval (lower, upper)
    decimals : int
        Number of decimal places

    Returns
    -------
    str
        Formatted string
    """
    return f"{estimate:.{decimals}f} ({ci[0]:.{decimals}f}, {ci[1]:.{decimals}f})"


def apply_delta_method(
    func,
    estimate: float,
    se: float,
    epsilon: float = 1e-6
) -> Tuple[float, float]:
    """
    Apply delta method for variance of transformed estimates.

    Parameters
    ----------
    func : callable
        Transformation function
    estimate : float
        Original estimate
    se : float
        Standard error of original estimate
    epsilon : float
        Small value for numerical derivative

    Returns
    -------
    tuple of float
        Transformed estimate and standard error
    """
    # Transform estimate
    transformed_estimate = func(estimate)

    # Approximate derivative
    derivative = (func(estimate + epsilon) - func(estimate - epsilon)) / (2 * epsilon)

    # Transform SE using delta method
    transformed_se = abs(derivative) * se

    return transformed_estimate, transformed_se
