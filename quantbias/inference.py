"""
Statistical inference with proper CI propagation for bias-corrected estimates.

Implements:
- Bootstrap CI for bias-corrected estimates
- Delta method for variance propagation
- Accounting for additional uncertainty from bias parameters

References
----------
Fox MP, Lash TL. On the Need for Quantitative Bias Analysis in the Peer-Review Process.
Am J Epidemiol. 2017;185(10):865-868.

Greenland S. Interval estimation by simulation as an alternative to and extension of
confidence intervals. Int J Epidemiol. 2004;33(6):1389-1397.

Lash TL, Fox MP, MacLehose RF. Applying Quantitative Bias Analysis to Epidemiologic Data.
2nd ed. Springer, 2021. Chapter 11: Probabilistic Bias Analysis.
"""

import numpy as np
from typing import Optional, Tuple, Callable, List
from dataclasses import dataclass
from scipy import stats
from quantbias.utils.helpers import validate_rr, validate_ci
from quantbias.utils.exceptions import InvalidParameterError


@dataclass
class CorrectedEstimate:
    """
    Container for bias-corrected estimate with proper uncertainty.

    Attributes
    ----------
    point_estimate : float
        Bias-corrected point estimate
    ci_lower : float
        Lower confidence bound (accounting for bias uncertainty)
    ci_upper : float
        Upper confidence bound (accounting for bias uncertainty)
    se_conventional : float
        Conventional standard error (from original study)
    se_bias : float
        Additional standard error from bias correction
    se_total : float
        Total standard error (conventional + bias)
    method : str
        Method used for CI ('bootstrap', 'delta', 'simulation')
    """
    point_estimate: float
    ci_lower: float
    ci_upper: float
    se_conventional: float
    se_bias: float
    se_total: float
    method: str

    def __str__(self) -> str:
        return (
            f"Bias-Corrected Estimate:\n"
            f"  Point estimate: {self.point_estimate:.3f}\n"
            f"  95% CI: ({self.ci_lower:.3f}, {self.ci_upper:.3f})\n"
            f"  SE (conventional): {self.se_conventional:.4f}\n"
            f"  SE (bias): {self.se_bias:.4f}\n"
            f"  SE (total): {self.se_total:.4f}\n"
            f"  Method: {self.method}\n"
        )


class BiasInference:
    """
    Statistical inference for bias-corrected estimates.

    Provides proper confidence intervals that account for:
    1. Random error from original study
    2. Uncertainty in bias parameters
    3. Combined uncertainty
    """

    @staticmethod
    def bootstrap_ci(
        observed_effect: float,
        observed_ci: Tuple[float, float],
        bias_correction_func: Callable,
        n_bootstrap: int = 10000,
        alpha: float = 0.05,
        random_state: Optional[int] = None
    ) -> CorrectedEstimate:
        """
        Bootstrap confidence interval for bias-corrected estimate.

        Parameters
        ----------
        observed_effect : float
            Observed effect measure (e.g., RR)
        observed_ci : tuple
            Original confidence interval
        bias_correction_func : callable
            Function that takes observed effect and returns corrected effect
            Should sample bias parameters internally
        n_bootstrap : int
            Number of bootstrap iterations
        alpha : float
            Significance level (default 0.05 for 95% CI)
        random_state : int, optional
            Random seed

        Returns
        -------
        CorrectedEstimate
            Bias-corrected estimate with bootstrap CI

        Notes
        -----
        Bootstrap procedure:
        1. Sample from distribution of observed effect (using SE from CI)
        2. For each sample, apply bias correction with sampled bias parameters
        3. Calculate percentile CI from bootstrap distribution

        This properly accounts for both random error and bias parameter uncertainty.

        References
        ----------
        Lash et al. 2021, Chapter 11.
        Greenland S. Int J Epidemiol. 2004;33(6):1389-1397.
        """
        validate_rr(observed_effect)
        validate_ci(observed_ci, observed_effect)

        rng = np.random.default_rng(random_state)

        # Calculate SE from observed CI (log scale for ratio measures)
        from quantbias.utils.helpers import compute_standard_error
        se_log = compute_standard_error(observed_ci[0], observed_ci[1], alpha)

        # Bootstrap distribution
        bootstrap_estimates = np.zeros(n_bootstrap)

        for i in range(n_bootstrap):
            # Sample observed effect from its distribution
            sampled_log_effect = rng.normal(np.log(observed_effect), se_log)
            sampled_effect = np.exp(sampled_log_effect)

            # Apply bias correction (function should handle parameter sampling)
            corrected_effect = bias_correction_func(sampled_effect)

            bootstrap_estimates[i] = corrected_effect

        # Point estimate (median of bootstrap distribution)
        corrected_point = np.median(bootstrap_estimates)

        # Confidence interval (percentile method)
        ci_lower = np.percentile(bootstrap_estimates, 100 * alpha / 2)
        ci_upper = np.percentile(bootstrap_estimates, 100 * (1 - alpha / 2))

        # Standard errors
        se_conventional = se_log  # From original study
        se_total_log = np.std(np.log(bootstrap_estimates))
        se_bias = np.sqrt(max(0, se_total_log**2 - se_conventional**2))

        return CorrectedEstimate(
            point_estimate=corrected_point,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            se_conventional=se_conventional,
            se_bias=se_bias,
            se_total=se_total_log,
            method='bootstrap'
        )

    @staticmethod
    def delta_method_ci(
        observed_effect: float,
        observed_ci: Tuple[float, float],
        bias_factor: float,
        bias_factor_se: float = 0.0,
        alpha: float = 0.05
    ) -> CorrectedEstimate:
        """
        Delta method confidence interval for bias-corrected estimate.

        Parameters
        ----------
        observed_effect : float
            Observed effect measure
        observed_ci : tuple
            Original confidence interval
        bias_factor : float
            Bias factor (RR_observed / RR_true)
        bias_factor_se : float
            Standard error of bias factor (on log scale)
        alpha : float
            Significance level

        Returns
        -------
        CorrectedEstimate
            Bias-corrected estimate with delta method CI

        Notes
        -----
        Delta method approximation:
        Var(RR_corrected) ≈ Var(RR_obs) + Var(BF)

        On log scale:
        SE(log RR_corrected) = sqrt(SE(log RR_obs)^2 + SE(log BF)^2)

        This assumes independence between random error and bias factor.

        References
        ----------
        Greenland S. Bayesian perspectives for epidemiological research: I.
        Foundations and basic methods. Int J Epidemiol. 2006;35(3):765-775.
        """
        validate_rr(observed_effect)
        validate_ci(observed_ci, observed_effect)
        validate_rr(bias_factor, "bias factor")

        # Correct point estimate
        corrected_point = observed_effect / bias_factor

        # Standard errors (log scale)
        from quantbias.utils.helpers import compute_standard_error, compute_ci_from_se

        se_obs_log = compute_standard_error(observed_ci[0], observed_ci[1], alpha)

        # Total SE combines random error and bias uncertainty
        se_total_log = np.sqrt(se_obs_log**2 + bias_factor_se**2)

        # Construct CI
        ci_lower, ci_upper = compute_ci_from_se(corrected_point, se_total_log, alpha)

        return CorrectedEstimate(
            point_estimate=corrected_point,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            se_conventional=se_obs_log,
            se_bias=bias_factor_se,
            se_total=se_total_log,
            method='delta'
        )

    @staticmethod
    def simulation_ci(
        observed_effect: float,
        observed_ci: Tuple[float, float],
        bias_parameters: List[dict],
        n_simulations: int = 10000,
        alpha: float = 0.05,
        random_state: Optional[int] = None
    ) -> CorrectedEstimate:
        """
        Simulation-based CI for bias correction with multiple parameters.

        Parameters
        ----------
        observed_effect : float
            Observed effect measure
        observed_ci : tuple
            Original confidence interval
        bias_parameters : list of dict
            List of bias parameter specifications with distributions
        n_simulations : int
            Number of simulation iterations
        alpha : float
            Significance level
        random_state : int, optional
            Random seed

        Returns
        -------
        CorrectedEstimate
            Bias-corrected estimate with simulation-based CI

        Notes
        -----
        Simulation procedure (Lash et al. 2021, Algorithm 11.1):
        1. For each iteration:
           a. Sample observed effect from its distribution
           b. Sample each bias parameter from its distribution
           c. Calculate bias-corrected effect
        2. Calculate median and percentiles from distribution

        This is the most flexible approach for complex bias models.

        References
        ----------
        Lash TL et al. 2021, Chapter 11, Algorithm 11.1.
        Greenland S. Int J Epidemiol. 2004;33(6):1389-1397.
        """
        validate_rr(observed_effect)
        validate_ci(observed_ci, observed_effect)

        rng = np.random.default_rng(random_state)

        # Calculate SE from observed CI
        from quantbias.utils.helpers import compute_standard_error
        se_log = compute_standard_error(observed_ci[0], observed_ci[1], alpha)

        # Simulation
        simulated_estimates = np.zeros(n_simulations)

        for i in range(n_simulations):
            # Sample observed effect
            sampled_log_effect = rng.normal(np.log(observed_effect), se_log)
            sampled_effect = np.exp(sampled_log_effect)

            # Apply bias corrections
            corrected = sampled_effect

            for bias_param in bias_parameters:
                # Sample bias parameter
                if 'distribution' in bias_param:
                    dist_type = bias_param['distribution']['type']
                    params = bias_param['distribution']['params']

                    if dist_type == 'lognormal':
                        bias_sample = rng.lognormal(params['mean'], params['sigma'])
                    elif dist_type == 'uniform':
                        bias_sample = rng.uniform(params['low'], params['high'])
                    elif dist_type == 'beta':
                        bias_sample = rng.beta(params['a'], params['b'])
                    else:
                        bias_sample = bias_param['value']
                else:
                    bias_sample = bias_param['value']

                # Apply correction
                if bias_param['type'] == 'multiplicative':
                    corrected = corrected / bias_sample
                elif bias_param['type'] == 'additive':
                    corrected = corrected - bias_sample

            simulated_estimates[i] = corrected

        # Results
        corrected_point = np.median(simulated_estimates)
        ci_lower = np.percentile(simulated_estimates, 100 * alpha / 2)
        ci_upper = np.percentile(simulated_estimates, 100 * (1 - alpha / 2))

        # Standard errors
        se_conventional = se_log
        se_total_log = np.std(np.log(simulated_estimates))
        se_bias = np.sqrt(max(0, se_total_log**2 - se_conventional**2))

        return CorrectedEstimate(
            point_estimate=corrected_point,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            se_conventional=se_conventional,
            se_bias=se_bias,
            se_total=se_total_log,
            method='simulation'
        )

    @staticmethod
    def variance_inflation_factor(
        se_conventional: float,
        se_bias: float
    ) -> float:
        """
        Calculate variance inflation factor from bias correction.

        Parameters
        ----------
        se_conventional : float
            Standard error from original study
        se_bias : float
            Additional standard error from bias correction

        Returns
        -------
        float
            Variance inflation factor (VIF)

        Notes
        -----
        VIF = Var(total) / Var(conventional)
            = (SE_conventional^2 + SE_bias^2) / SE_conventional^2
            = 1 + (SE_bias / SE_conventional)^2

        VIF > 1.5 indicates substantial additional uncertainty from bias.

        References
        ----------
        Fox MP, Lash TL. Am J Epidemiol. 2017;185(10):865-868.
        """
        if se_conventional <= 0:
            raise InvalidParameterError("Conventional SE must be positive")

        vif = 1 + (se_bias / se_conventional)**2

        return vif

    @staticmethod
    def compare_methods(
        observed_effect: float,
        observed_ci: Tuple[float, float],
        bias_factor: float,
        bias_factor_se: float = 0.1
    ) -> dict:
        """
        Compare different CI methods for same bias correction.

        Parameters
        ----------
        observed_effect : float
            Observed effect measure
        observed_ci : tuple
            Original confidence interval
        bias_factor : float
            Bias factor
        bias_factor_se : float
            SE of bias factor

        Returns
        -------
        dict
            Comparison of different methods
        """
        # Naive (incorrect) method
        naive_point = observed_effect / bias_factor
        naive_ci_lower = observed_ci[0] / bias_factor
        naive_ci_upper = observed_ci[1] / bias_factor

        # Delta method
        delta_result = BiasInference.delta_method_ci(
            observed_effect, observed_ci, bias_factor, bias_factor_se
        )

        # Calculate widths
        naive_width = naive_ci_upper - naive_ci_lower
        delta_width = delta_result.ci_upper - delta_result.ci_lower

        # Variance inflation
        vif = BiasInference.variance_inflation_factor(
            delta_result.se_conventional,
            delta_result.se_bias
        )

        return {
            'naive': {
                'point': naive_point,
                'ci': (naive_ci_lower, naive_ci_upper),
                'width': naive_width
            },
            'delta': {
                'point': delta_result.point_estimate,
                'ci': (delta_result.ci_lower, delta_result.ci_upper),
                'width': delta_width
            },
            'width_inflation': delta_width / naive_width,
            'variance_inflation_factor': vif,
            'interpretation': (
                f"Naive CI width: {naive_width:.3f}\n"
                f"Correct CI width: {delta_width:.3f}\n"
                f"Width inflation: {delta_width/naive_width:.2f}x\n"
                f"Variance inflation factor: {vif:.2f}\n"
                f"{'WARNING: Substantial uncertainty from bias correction!' if vif > 1.5 else 'Moderate additional uncertainty'}"
            )
        }
