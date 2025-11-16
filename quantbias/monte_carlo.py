"""
Monte Carlo sensitivity analysis for probabilistic bias analysis.

Implements:
- Probabilistic bias analysis
- Multiple bias modeling
- Uncertainty propagation
"""

import numpy as np
import pandas as pd
from typing import List, Optional, Tuple, Dict, Union
from dataclasses import dataclass
from scipy import stats
from tqdm import tqdm

from quantbias.bias_parameters import BiasParameter, SelectionBias, MeasurementError, Confounding
from quantbias.utils.helpers import validate_rr, compute_ci_from_se
from quantbias.utils.constants import DEFAULT_RANDOM_SEED


@dataclass
class MonteCarloResult:
    """
    Container for Monte Carlo analysis results.

    Attributes
    ----------
    observed_rr : float
        Original observed RR
    adjusted_rr_samples : np.ndarray
        Distribution of bias-adjusted RRs
    median : float
        Median adjusted RR
    mean : float
        Mean adjusted RR
    ci_lower : float
        Lower confidence bound (2.5th percentile)
    ci_upper : float
        Upper confidence bound (97.5th percentile)
    percentiles : dict
        Additional percentiles
    """
    observed_rr: float
    adjusted_rr_samples: np.ndarray
    median: float
    mean: float
    ci_lower: float
    ci_upper: float
    percentiles: Dict[float, float]

    def summary(self) -> str:
        """Return summary string."""
        return (
            f"Monte Carlo Bias Analysis Results:\n"
            f"  Observed RR: {self.observed_rr:.3f}\n"
            f"  Adjusted RR (median): {self.median:.3f}\n"
            f"  Adjusted RR (mean): {self.mean:.3f}\n"
            f"  95% Simulation Interval: ({self.ci_lower:.3f}, {self.ci_upper:.3f})\n"
            f"  Interquartile Range: ({self.percentiles[25]:.3f}, {self.percentiles[75]:.3f})\n"
        )


class MonteCarloAnalysis:
    """
    Monte Carlo sensitivity analysis.

    Performs probabilistic bias analysis by:
    1. Specifying uncertainty distributions for bias parameters
    2. Sampling from these distributions
    3. Applying bias corrections
    4. Analyzing distribution of bias-adjusted estimates
    """

    def __init__(
        self,
        observed_rr: float,
        bias_parameters: List[BiasParameter],
        n_iterations: int = 10000,
        random_state: Optional[int] = None
    ):
        """
        Initialize Monte Carlo analysis.

        Parameters
        ----------
        observed_rr : float
            Observed risk ratio
        bias_parameters : list of BiasParameter
            List of bias parameters to simulate
        n_iterations : int
            Number of Monte Carlo iterations
        random_state : int, optional
            Random seed for reproducibility
        """
        validate_rr(observed_rr)
        self.observed_rr = observed_rr
        self.bias_parameters = bias_parameters
        self.n_iterations = n_iterations
        self.random_state = random_state or DEFAULT_RANDOM_SEED
        self.rng = np.random.default_rng(self.random_state)

        # Results storage
        self.results: Optional[MonteCarloResult] = None

    def run(self, show_progress: bool = True) -> MonteCarloResult:
        """
        Run Monte Carlo simulation.

        Parameters
        ----------
        show_progress : bool
            Show progress bar

        Returns
        -------
        MonteCarloResult
            Results of Monte Carlo analysis
        """
        adjusted_rrs = np.zeros(self.n_iterations)

        iterator = tqdm(range(self.n_iterations), desc="Monte Carlo simulation") \
            if show_progress else range(self.n_iterations)

        for i in iterator:
            # Start with observed RR
            current_rr = self.observed_rr

            # Apply each bias correction
            for bias_param in self.bias_parameters:
                # Sample bias parameter values
                if isinstance(bias_param, SelectionBias):
                    current_rr = self._apply_selection_bias_sample(current_rr, bias_param, i)
                elif isinstance(bias_param, MeasurementError):
                    current_rr = self._apply_measurement_error_sample(current_rr, bias_param, i)
                elif isinstance(bias_param, Confounding):
                    current_rr = self._apply_confounding_sample(current_rr, bias_param, i)
                else:
                    # Generic bias parameter
                    current_rr = bias_param.apply_bias(current_rr)

            adjusted_rrs[i] = current_rr

        # Compute summary statistics
        median = np.median(adjusted_rrs)
        mean = np.mean(adjusted_rrs)
        ci_lower = np.percentile(adjusted_rrs, 2.5)
        ci_upper = np.percentile(adjusted_rrs, 97.5)

        percentiles = {
            5: np.percentile(adjusted_rrs, 5),
            10: np.percentile(adjusted_rrs, 10),
            25: np.percentile(adjusted_rrs, 25),
            50: median,
            75: np.percentile(adjusted_rrs, 75),
            90: np.percentile(adjusted_rrs, 90),
            95: np.percentile(adjusted_rrs, 95),
        }

        self.results = MonteCarloResult(
            observed_rr=self.observed_rr,
            adjusted_rr_samples=adjusted_rrs,
            median=median,
            mean=mean,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            percentiles=percentiles
        )

        return self.results

    def _apply_selection_bias_sample(
        self,
        rr: float,
        bias_param: SelectionBias,
        iteration: int
    ) -> float:
        """Apply sampled selection bias correction."""
        # Sample parameters
        params = bias_param.sample_parameters(n=1, random_state=self.random_state + iteration)

        # Create temporary bias parameter with sampled values
        temp_bias = SelectionBias(
            s11=params['s11'][0],
            s10=params['s10'][0],
            s01=params['s01'][0],
            s00=params['s00'][0]
        )

        return temp_bias.apply_bias(rr)

    def _apply_measurement_error_sample(
        self,
        rr: float,
        bias_param: MeasurementError,
        iteration: int
    ) -> float:
        """Apply sampled measurement error correction."""
        params = bias_param.sample_parameters(n=1, random_state=self.random_state + iteration)

        temp_bias = MeasurementError(
            sensitivity_exposure=params['sensitivity_exposure'][0],
            specificity_exposure=params['specificity_exposure'][0],
            sensitivity_outcome=params['sensitivity_outcome'][0],
            specificity_outcome=params['specificity_outcome'][0],
            nondifferential=bias_param.nondifferential
        )

        return temp_bias.apply_bias(rr)

    def _apply_confounding_sample(
        self,
        rr: float,
        bias_param: Confounding,
        iteration: int
    ) -> float:
        """Apply sampled confounding correction."""
        params = bias_param.sample_parameters(n=1, random_state=self.random_state + iteration)

        # Get sampled parameters
        p0 = params['prevalence_confounder_unexposed'][0]
        p1 = params['prevalence_confounder_exposed'][0]
        rr_cd0 = params['rr_confounder_outcome_unexposed'][0]
        rr_cd1 = params['rr_confounder_outcome_exposed'][0]

        temp_bias = Confounding(
            rr_confounder_outcome_unexposed=rr_cd0,
            rr_confounder_outcome_exposed=rr_cd1,
            prevalence_confounder_unexposed=p0,
            prevalence_confounder_exposed=p1,
            rr_confounder_exposure=bias_param.rr_confounder_exposure
        )

        return temp_bias.apply_bias(rr)

    def get_distribution_summary(self) -> Dict[str, float]:
        """
        Get summary statistics of adjusted RR distribution.

        Returns
        -------
        dict
            Distribution summary with keys: mean, median, std, ci_lower, ci_upper
        """
        if self.results is None:
            raise ValueError("Must run simulation first using .run()")

        summary_stats = {
            'mean': self.results.mean,
            'median': self.results.median,
            'std': np.std(self.results.adjusted_rr_samples),
            'ci_lower': self.results.ci_lower,
            'ci_upper': self.results.ci_upper
        }

        return summary_stats

    def probability_below_threshold(self, threshold: float = 1.0) -> float:
        """
        Calculate probability that adjusted RR is below threshold.

        Parameters
        ----------
        threshold : float
            Threshold value (default 1.0 for null)

        Returns
        -------
        float
            Probability
        """
        if self.results is None:
            raise ValueError("Must run simulation first using .run()")

        return np.mean(self.results.adjusted_rr_samples < threshold)

    def probability_above_threshold(self, threshold: float = 1.0) -> float:
        """
        Calculate probability that adjusted RR is above threshold.

        Parameters
        ----------
        threshold : float
            Threshold value

        Returns
        -------
        float
            Probability
        """
        return 1 - self.probability_below_threshold(threshold)

    def probability_rr_greater_than(self, threshold: float) -> float:
        """
        Calculate probability that adjusted RR is greater than threshold.

        Parameters
        ----------
        threshold : float
            Threshold value

        Returns
        -------
        float
            Probability
        """
        return self.probability_above_threshold(threshold)

    def probability_rr_less_than(self, threshold: float) -> float:
        """
        Calculate probability that adjusted RR is less than threshold.

        Parameters
        ----------
        threshold : float
            Threshold value

        Returns
        -------
        float
            Probability
        """
        return self.probability_below_threshold(threshold)

    def get_samples(self) -> np.ndarray:
        """
        Get array of bias-adjusted RR samples.

        Returns
        -------
        np.ndarray
            Array of adjusted RR values
        """
        if self.results is None:
            raise ValueError("Must run simulation first using .run()")

        return self.results.adjusted_rr_samples


class MultipleBiasModeling:
    """
    Framework for modeling multiple biases simultaneously.

    Allows for:
    - Sequential bias correction
    - Correlated bias parameters
    - Joint uncertainty assessment
    """

    def __init__(
        self,
        observed_rr: float,
        observed_ci: Optional[Tuple[float, float]] = None
    ):
        """
        Initialize multiple bias modeling.

        Parameters
        ----------
        observed_rr : float
            Observed risk ratio
        observed_ci : tuple of float, optional
            Confidence interval
        """
        validate_rr(observed_rr)
        self.observed_rr = observed_rr
        self.observed_ci = observed_ci
        self.bias_parameters: List[BiasParameter] = []

    def add_bias(self, bias_parameter: BiasParameter) -> None:
        """
        Add a bias parameter to the model.

        Parameters
        ----------
        bias_parameter : BiasParameter
            Bias parameter to add
        """
        self.bias_parameters.append(bias_parameter)

    def deterministic_correction(self) -> float:
        """
        Apply deterministic bias correction using point estimates.

        Returns
        -------
        float
            Bias-corrected RR
        """
        corrected_rr = self.observed_rr

        for bias_param in self.bias_parameters:
            corrected_rr = bias_param.apply_bias(corrected_rr)

        return corrected_rr

    def probabilistic_correction(
        self,
        n_iterations: int = 10000,
        random_state: Optional[int] = None,
        show_progress: bool = True
    ) -> MonteCarloResult:
        """
        Apply probabilistic bias correction using Monte Carlo.

        Parameters
        ----------
        n_iterations : int
            Number of iterations
        random_state : int, optional
            Random seed
        show_progress : bool
            Show progress bar

        Returns
        -------
        MonteCarloResult
            Monte Carlo results
        """
        mc = MonteCarloAnalysis(
            observed_rr=self.observed_rr,
            bias_parameters=self.bias_parameters,
            n_iterations=n_iterations,
            random_state=random_state
        )

        return mc.run(show_progress=show_progress)

    def bias_decomposition(self) -> pd.DataFrame:
        """
        Decompose total bias into components.

        Returns
        -------
        pd.DataFrame
            Bias decomposition showing effect of each bias source
        """
        results = []

        # Original effect
        current_rr = self.observed_rr
        results.append({
            'Step': 'Observed',
            'Bias_source': 'None',
            'RR': current_rr,
            'Bias_factor': 1.0
        })

        # Apply each bias sequentially and track
        for i, bias_param in enumerate(self.bias_parameters):
            previous_rr = current_rr
            current_rr = bias_param.apply_bias(current_rr)
            bias_factor = previous_rr / current_rr

            results.append({
                'Step': f'After {bias_param.name}',
                'Bias_source': bias_param.name,
                'RR': current_rr,
                'Bias_factor': bias_factor
            })

        return pd.DataFrame(results)
