"""
Bias parameter specification for quantitative bias analysis.

Implements bias parameters for:
- Selection bias
- Information bias (measurement error)
- Confounding bias
"""

import numpy as np
from typing import Optional, Tuple, Union, Dict, Any
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
from quantbias.utils.helpers import validate_probability, validate_rr
from quantbias.utils.exceptions import InvalidParameterError
from quantbias.utils.constants import BiasType


@dataclass
class BiasParameter(ABC):
    """
    Abstract base class for bias parameters.

    All bias parameters should inherit from this class and implement
    the apply_bias method.
    """
    name: str = ""
    bias_type: BiasType = field(default=BiasType.CONFOUNDING)

    @abstractmethod
    def apply_bias(self, observed_effect: float) -> float:
        """
        Apply bias correction to observed effect.

        Parameters
        ----------
        observed_effect : float
            Observed effect measure

        Returns
        -------
        float
            Bias-corrected effect measure
        """
        pass

    @abstractmethod
    def validate(self) -> None:
        """Validate bias parameter values."""
        pass


@dataclass
class SelectionBias(BiasParameter):
    """
    Selection bias parameters.

    Models selection bias using sensitivity and specificity of selection
    into the study.

    Attributes
    ----------
    sensitivity : float
        Probability of selection given exposed and diseased
    specificity : float
        Probability of selection given unexposed and not diseased
    selection_probability : float
        Overall probability of selection into study
    sensitivity_distribution : tuple, optional
        Distribution for Monte Carlo (dist_name, params)
    specificity_distribution : tuple, optional
        Distribution for Monte Carlo (dist_name, params)

    References
    ----------
    Lash TL, Fox MP, Fink AK. Applying Quantitative Bias Analysis to
    Epidemiologic Data. Springer, 2009.
    """
    sensitivity: float = 0.8
    specificity: float = 0.8
    selection_probability: float = 0.5
    sensitivity_distribution: Optional[Tuple[str, Dict[str, Any]]] = None
    specificity_distribution: Optional[Tuple[str, Dict[str, Any]]] = None

    def __post_init__(self):
        self.name = "Selection Bias"
        self.bias_type = BiasType.SELECTION
        self.validate()

    def validate(self) -> None:
        """Validate selection bias parameters."""
        validate_probability(self.sensitivity, "sensitivity")
        validate_probability(self.specificity, "specificity")
        validate_probability(self.selection_probability, "selection probability")

    def apply_bias(self, observed_rr: float) -> float:
        """
        Apply selection bias correction.

        Parameters
        ----------
        observed_rr : float
            Observed risk ratio

        Returns
        -------
        float
            Bias-corrected risk ratio

        Notes
        -----
        The correction formula accounts for differential selection based on
        exposure and outcome status.
        """
        validate_rr(observed_rr)

        # Selection bias correction
        # Simplified model: assumes selection depends on exposure-outcome combination
        se = self.sensitivity
        sp = self.specificity

        # Bias factor
        bias_factor = (se * sp) / ((1 - se) * (1 - sp))

        corrected_rr = observed_rr / bias_factor
        return max(corrected_rr, 0.001)  # Ensure positive

    def sample_parameters(self, n: int = 1, random_state: Optional[int] = None) -> Dict[str, np.ndarray]:
        """
        Sample parameters from distributions for Monte Carlo analysis.

        Parameters
        ----------
        n : int
            Number of samples
        random_state : int, optional
            Random seed

        Returns
        -------
        dict
            Sampled parameter values
        """
        rng = np.random.default_rng(random_state)

        params = {}

        if self.sensitivity_distribution:
            dist_name, dist_params = self.sensitivity_distribution
            if dist_name == 'beta':
                params['sensitivity'] = rng.beta(
                    dist_params['a'],
                    dist_params['b'],
                    size=n
                )
            elif dist_name == 'uniform':
                params['sensitivity'] = rng.uniform(
                    dist_params['low'],
                    dist_params['high'],
                    size=n
                )
        else:
            params['sensitivity'] = np.full(n, self.sensitivity)

        if self.specificity_distribution:
            dist_name, dist_params = self.specificity_distribution
            if dist_name == 'beta':
                params['specificity'] = rng.beta(
                    dist_params['a'],
                    dist_params['b'],
                    size=n
                )
            elif dist_name == 'uniform':
                params['specificity'] = rng.uniform(
                    dist_params['low'],
                    dist_params['high'],
                    size=n
                )
        else:
            params['specificity'] = np.full(n, self.specificity)

        return params


@dataclass
class MeasurementError(BiasParameter):
    """
    Information bias / measurement error parameters.

    Models misclassification of exposure and/or outcome.

    Attributes
    ----------
    sensitivity_exposure : float
        Sensitivity of exposure measurement
    specificity_exposure : float
        Specificity of exposure measurement
    sensitivity_outcome : float
        Sensitivity of outcome measurement
    specificity_outcome : float
        Specificity of outcome measurement
    nondifferential : bool
        Whether misclassification is non-differential (default True)

    References
    ----------
    Greenland S. The effect of misclassification in the presence of covariates.
    Am J Epidemiol. 1980;112(4):564-569.
    """
    sensitivity_exposure: float = 0.9
    specificity_exposure: float = 0.9
    sensitivity_outcome: float = 0.95
    specificity_outcome: float = 0.95
    nondifferential: bool = True
    sensitivity_exposure_dist: Optional[Tuple[str, Dict[str, Any]]] = None
    specificity_exposure_dist: Optional[Tuple[str, Dict[str, Any]]] = None
    sensitivity_outcome_dist: Optional[Tuple[str, Dict[str, Any]]] = None
    specificity_outcome_dist: Optional[Tuple[str, Dict[str, Any]]] = None

    def __post_init__(self):
        self.name = "Measurement Error"
        self.bias_type = BiasType.MEASUREMENT_ERROR
        self.validate()

    def validate(self) -> None:
        """Validate measurement error parameters."""
        validate_probability(self.sensitivity_exposure, "exposure sensitivity")
        validate_probability(self.specificity_exposure, "exposure specificity")
        validate_probability(self.sensitivity_outcome, "outcome sensitivity")
        validate_probability(self.specificity_outcome, "outcome specificity")

    def apply_bias(self, observed_rr: float) -> float:
        """
        Apply measurement error correction.

        Parameters
        ----------
        observed_rr : float
            Observed risk ratio

        Returns
        -------
        float
            Bias-corrected risk ratio

        Notes
        -----
        For non-differential misclassification, the bias is typically toward
        the null. This method applies the correction formula from Greenland (1980).
        """
        validate_rr(observed_rr)

        se_exp = self.sensitivity_exposure
        sp_exp = self.specificity_exposure
        se_out = self.sensitivity_outcome
        sp_out = self.specificity_outcome

        if self.nondifferential:
            # Non-differential misclassification correction
            # Exposure misclassification
            exp_bias_factor = (se_exp + sp_exp - 1)

            # Outcome misclassification
            out_bias_factor = (se_out + sp_out - 1)

            # Combined bias factor
            bias_factor = exp_bias_factor * out_bias_factor

            # Correct the log-RR
            log_observed_rr = np.log(observed_rr)
            log_corrected_rr = log_observed_rr / bias_factor

            corrected_rr = np.exp(log_corrected_rr)
        else:
            # Differential misclassification (more complex)
            # Simplified correction
            bias_factor = (se_exp * sp_out) / ((1 - se_exp) * (1 - sp_out))
            corrected_rr = observed_rr / bias_factor

        return max(corrected_rr, 0.001)

    def sample_parameters(self, n: int = 1, random_state: Optional[int] = None) -> Dict[str, np.ndarray]:
        """Sample parameters for Monte Carlo analysis."""
        rng = np.random.default_rng(random_state)
        params = {}

        # Helper function to sample from distribution
        def sample_dist(dist_spec, default_value, param_name):
            if dist_spec:
                dist_name, dist_params = dist_spec
                if dist_name == 'beta':
                    return rng.beta(dist_params['a'], dist_params['b'], size=n)
                elif dist_name == 'uniform':
                    return rng.uniform(dist_params['low'], dist_params['high'], size=n)
            return np.full(n, default_value)

        params['sensitivity_exposure'] = sample_dist(
            self.sensitivity_exposure_dist, self.sensitivity_exposure, 'sensitivity_exposure'
        )
        params['specificity_exposure'] = sample_dist(
            self.specificity_exposure_dist, self.specificity_exposure, 'specificity_exposure'
        )
        params['sensitivity_outcome'] = sample_dist(
            self.sensitivity_outcome_dist, self.sensitivity_outcome, 'sensitivity_outcome'
        )
        params['specificity_outcome'] = sample_dist(
            self.specificity_outcome_dist, self.specificity_outcome, 'specificity_outcome'
        )

        return params


@dataclass
class Confounding(BiasParameter):
    """
    Unmeasured confounding parameters.

    Models unmeasured confounding using the relationship between confounder,
    exposure, and outcome.

    Attributes
    ----------
    rr_confounder_exposure : float
        Risk ratio for confounder-exposure association
    rr_confounder_outcome : float
        Risk ratio for confounder-outcome association
    prevalence_confounder_unexposed : float
        Prevalence of confounder among unexposed
    prevalence_confounder_exposed : float, optional
        Prevalence of confounder among exposed (calculated if not provided)

    References
    ----------
    Greenland S, Lash TL. Bias Analysis. In: Modern Epidemiology, 3rd ed.
    Lippincott Williams & Wilkins, 2008.
    """
    rr_confounder_exposure: float = 2.0
    rr_confounder_outcome: float = 2.0
    prevalence_confounder_unexposed: float = 0.3
    prevalence_confounder_exposed: Optional[float] = None
    rr_conf_exp_dist: Optional[Tuple[str, Dict[str, Any]]] = None
    rr_conf_out_dist: Optional[Tuple[str, Dict[str, Any]]] = None
    prev_conf_unexp_dist: Optional[Tuple[str, Dict[str, Any]]] = None

    def __post_init__(self):
        self.name = "Unmeasured Confounding"
        self.bias_type = BiasType.CONFOUNDING
        self.validate()

        # Calculate prevalence in exposed if not provided
        if self.prevalence_confounder_exposed is None:
            # Use RR to calculate prevalence in exposed
            p0 = self.prevalence_confounder_unexposed
            rr_ce = self.rr_confounder_exposure
            # p1 = p0 * RR_CE / (1 + p0 * (RR_CE - 1))
            self.prevalence_confounder_exposed = min(
                p0 * rr_ce / (1 + p0 * (rr_ce - 1)),
                0.99
            )

    def validate(self) -> None:
        """Validate confounding parameters."""
        validate_rr(self.rr_confounder_exposure, "RR confounder-exposure")
        validate_rr(self.rr_confounder_outcome, "RR confounder-outcome")
        validate_probability(
            self.prevalence_confounder_unexposed,
            "prevalence confounder unexposed"
        )

    def apply_bias(self, observed_rr: float) -> float:
        """
        Apply confounding bias correction.

        Parameters
        ----------
        observed_rr : float
            Observed risk ratio

        Returns
        -------
        float
            Bias-corrected risk ratio

        Notes
        -----
        Uses the formula:
        RR_true = RR_obs / RR_bias
        where RR_bias is the bias factor from unmeasured confounding.
        """
        validate_rr(observed_rr)

        p0 = self.prevalence_confounder_unexposed
        p1 = self.prevalence_confounder_exposed
        rr_co = self.rr_confounder_outcome

        # Bias factor formula (Greenland & Lash)
        bias_factor = (p1 * rr_co + (1 - p1)) / (p0 * rr_co + (1 - p0))

        corrected_rr = observed_rr / bias_factor
        return max(corrected_rr, 0.001)

    def sample_parameters(self, n: int = 1, random_state: Optional[int] = None) -> Dict[str, np.ndarray]:
        """Sample parameters for Monte Carlo analysis."""
        rng = np.random.default_rng(random_state)
        params = {}

        def sample_dist(dist_spec, default_value):
            if dist_spec:
                dist_name, dist_params = dist_spec
                if dist_name == 'lognormal':
                    return rng.lognormal(dist_params['mean'], dist_params['sigma'], size=n)
                elif dist_name == 'uniform':
                    return rng.uniform(dist_params['low'], dist_params['high'], size=n)
                elif dist_name == 'beta':
                    return rng.beta(dist_params['a'], dist_params['b'], size=n)
            return np.full(n, default_value)

        params['rr_confounder_exposure'] = sample_dist(
            self.rr_conf_exp_dist, self.rr_confounder_exposure
        )
        params['rr_confounder_outcome'] = sample_dist(
            self.rr_conf_out_dist, self.rr_confounder_outcome
        )
        params['prevalence_confounder_unexposed'] = sample_dist(
            self.prev_conf_unexp_dist, self.prevalence_confounder_unexposed
        )

        return params
