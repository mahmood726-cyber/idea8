"""
Sensitivity analysis for unmeasured confounding.

Implements:
- Rosenbaum bounds for matched observational studies
- Sensitivity analysis grids
- Confounding parameter exploration
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, List, Dict, Union
from dataclasses import dataclass
from scipy import stats
from quantbias.utils.helpers import validate_rr, validate_probability, validate_ci
from quantbias.utils.exceptions import InvalidParameterError
from quantbias.evalues import EValue


@dataclass
class SensitivityResult:
    """
    Container for sensitivity analysis results.

    Attributes
    ----------
    observed_rr : float
        Original observed risk ratio
    adjusted_rr : float
        Bias-adjusted risk ratio
    bias_parameters : dict
        Parameters used for adjustment
    evalue : float
        E-value for this result
    """
    observed_rr: float
    adjusted_rr: float
    bias_parameters: Dict[str, float]
    evalue: float


class UnmeasuredConfounding:
    """
    Sensitivity analysis for unmeasured confounding.

    Provides methods to assess how unmeasured confounding could affect
    observed associations.
    """

    def __init__(
        self,
        observed_rr: float,
        observed_ci: Optional[Tuple[float, float]] = None
    ):
        """
        Initialize unmeasured confounding analysis.

        Parameters
        ----------
        observed_rr : float
            Observed risk ratio
        observed_ci : tuple of float, optional
            Confidence interval for observed RR
        """
        validate_rr(observed_rr)
        self.observed_rr = observed_rr
        self.observed_ci = observed_ci

        if observed_ci is not None:
            validate_ci(observed_ci, observed_rr)

    def adjust_for_confounding(
        self,
        rr_confounder_exposure: float,
        rr_confounder_outcome: float,
        prevalence_confounder_unexposed: float
    ) -> float:
        """
        Adjust observed RR for unmeasured confounding.

        Parameters
        ----------
        rr_confounder_exposure : float
            Association between confounder and exposure
        rr_confounder_outcome : float
            Association between confounder and outcome
        prevalence_confounder_unexposed : float
            Prevalence of confounder among unexposed

        Returns
        -------
        float
            Adjusted risk ratio

        Notes
        -----
        Uses the formula from Greenland & Lash (2008):
        RR_adjusted = RR_obs / RR_bias
        where RR_bias = [p1*RR_CD + (1-p1)] / [p0*RR_CD + (1-p0)]
        """
        validate_rr(rr_confounder_exposure, "RR confounder-exposure")
        validate_rr(rr_confounder_outcome, "RR confounder-outcome")
        validate_probability(prevalence_confounder_unexposed, "prevalence confounder unexposed")

        p0 = prevalence_confounder_unexposed

        # Calculate prevalence in exposed using RR
        # p1 = p0 * RR_CE / (1 - p0 + p0 * RR_CE)
        p1 = min(
            p0 * rr_confounder_exposure / (1 - p0 + p0 * rr_confounder_exposure),
            0.99
        )

        rr_co = rr_confounder_outcome

        # Calculate bias factor
        bias_factor = (p1 * rr_co + (1 - p1)) / (p0 * rr_co + (1 - p0))

        # Adjust observed RR
        adjusted_rr = self.observed_rr / bias_factor

        return max(adjusted_rr, 0.001)

    def sensitivity_grid(
        self,
        rr_confounder_outcome: Union[float, np.ndarray],
        prevalence_exposed: Union[float, np.ndarray],
        rr_confounder_exposure: Optional[Union[float, np.ndarray]] = None
    ) -> pd.DataFrame:
        """
        Create sensitivity analysis grid.

        Parameters
        ----------
        rr_confounder_outcome : float or array
            RR for confounder-outcome association(s)
        prevalence_exposed : float or array
            Prevalence of confounder among exposed
        rr_confounder_exposure : float or array, optional
            RR for confounder-exposure association (calculated if not provided)

        Returns
        -------
        pd.DataFrame
            Grid of adjusted RRs for different parameter combinations

        Examples
        --------
        >>> analysis = UnmeasuredConfounding(observed_rr=2.0)
        >>> grid = analysis.sensitivity_grid(
        ...     rr_confounder_outcome=np.arange(1.0, 4.0, 0.5),
        ...     prevalence_exposed=np.arange(0.1, 0.9, 0.1)
        ... )
        """
        # Convert to arrays
        rr_co_array = np.atleast_1d(rr_confounder_outcome)
        prev_exp_array = np.atleast_1d(prevalence_exposed)

        results = []

        for rr_co in rr_co_array:
            for p1 in prev_exp_array:
                # Calculate RR_CE if not provided
                if rr_confounder_exposure is None:
                    # Assume equal prevalence gives RR_CE = 1
                    # Calculate p0 from p1
                    p0 = p1 / 2  # Simplified assumption
                    rr_ce = p1 / (p0 + 1e-10)
                else:
                    rr_ce = rr_confounder_exposure if np.isscalar(rr_confounder_exposure) else rr_confounder_exposure[0]
                    p0 = p1 / rr_ce if rr_ce > 0 else p1

                p0 = min(max(p0, 0.01), 0.99)

                # Calculate adjusted RR
                adjusted_rr = self.adjust_for_confounding(
                    rr_confounder_exposure=rr_ce,
                    rr_confounder_outcome=rr_co,
                    prevalence_confounder_unexposed=p0
                )

                # Calculate E-value needed to produce this adjustment
                evalue = EValue.calculate_evalue_rr(rr_co)

                results.append({
                    'RR_confounder_outcome': rr_co,
                    'Prevalence_exposed': p1,
                    'Prevalence_unexposed': p0,
                    'RR_confounder_exposure': rr_ce,
                    'Adjusted_RR': adjusted_rr,
                    'Bias_factor': self.observed_rr / adjusted_rr,
                    'E_value': evalue
                })

        return pd.DataFrame(results)

    def threshold_analysis(
        self,
        target_rr: float = 1.0,
        prevalence_confounder: float = 0.3
    ) -> Dict[str, float]:
        """
        Find confounding strength needed to reduce observed RR to target.

        Parameters
        ----------
        target_rr : float
            Target risk ratio (default 1.0 for null)
        prevalence_confounder : float
            Assumed prevalence of confounder

        Returns
        -------
        dict
            Required confounding parameters

        Examples
        --------
        >>> analysis = UnmeasuredConfounding(observed_rr=2.5)
        >>> threshold = analysis.threshold_analysis(target_rr=1.0)
        >>> print(f"Required RR: {threshold['required_rr_confounder']:.2f}")
        """
        validate_rr(target_rr, "target RR")
        validate_probability(prevalence_confounder, "prevalence confounder")

        # Calculate required bias factor
        bias_factor = self.observed_rr / target_rr

        # Assuming equal prevalence in exposed and unexposed (conservative)
        p0 = p1 = prevalence_confounder

        # Solve for RR_CD needed to achieve bias_factor
        # bias_factor = (p1*RR_CD + (1-p1)) / (p0*RR_CD + (1-p0))
        # When p0 = p1, this simplifies to: bias_factor = 1 (no bias from equal prevalence)
        # So we need different prevalences or strong RR_CD

        # Simplified: assume we need RR_CD such that bias is achieved
        # RR_CD = (bias_factor - 1 + p0) / (p1 - bias_factor * p0)

        # For symmetric case, use E-value approximation
        evalue = EValue.calculate_evalue_rr(bias_factor)

        return {
            'target_rr': target_rr,
            'bias_factor_needed': bias_factor,
            'required_rr_confounder': evalue,
            'prevalence_assumed': prevalence_confounder,
            'interpretation': (
                f"An unmeasured confounder associated with both exposure and outcome "
                f"by RR = {evalue:.2f} could explain away the observed association "
                f"to RR = {target_rr:.2f}"
            )
        }


class SensitivityAnalysis:
    """
    Comprehensive sensitivity analysis framework.

    Provides methods for:
    - Tipping point analysis
    - Rule-out sensitivity values
    - Multi-parameter sensitivity analysis
    """

    def __init__(
        self,
        observed_rr: float,
        observed_ci: Optional[Tuple[float, float]] = None,
        observed_se: Optional[float] = None
    ):
        """
        Initialize sensitivity analysis.

        Parameters
        ----------
        observed_rr : float
            Observed risk ratio
        observed_ci : tuple of float, optional
            Confidence interval
        observed_se : float, optional
            Standard error (log scale)
        """
        validate_rr(observed_rr)
        self.observed_rr = observed_rr
        self.observed_ci = observed_ci
        self.observed_se = observed_se

        if observed_ci is not None:
            validate_ci(observed_ci, observed_rr)

            # Calculate SE from CI if not provided
            if observed_se is None:
                from quantbias.utils.helpers import compute_standard_error
                self.observed_se = compute_standard_error(observed_ci[0], observed_ci[1])

    def rosenbaum_bounds(
        self,
        gamma_values: Optional[np.ndarray] = None,
        test_statistic: Optional[float] = None
    ) -> pd.DataFrame:
        """
        Calculate Rosenbaum sensitivity bounds.

        For matched observational studies, calculates p-value bounds
        under different levels of hidden bias (gamma).

        Parameters
        ----------
        gamma_values : array, optional
            Sensitivity parameter values (default: 1.0 to 5.0)
        test_statistic : float, optional
            Test statistic from original analysis

        Returns
        -------
        pd.DataFrame
            Rosenbaum bounds for different gamma values

        References
        ----------
        Rosenbaum PR. Observational Studies. 2nd ed. Springer, 2002.
        """
        if gamma_values is None:
            gamma_values = np.arange(1.0, 5.1, 0.25)

        if test_statistic is None:
            # Use z-score from observed RR and SE
            if self.observed_se is not None:
                test_statistic = np.log(self.observed_rr) / self.observed_se
            else:
                # Approximate from CI
                test_statistic = np.log(self.observed_rr) / 0.5

        results = []

        for gamma in gamma_values:
            # Upper bound p-value (worst case)
            p_upper = 1 - stats.norm.cdf(test_statistic / gamma)

            # Lower bound p-value (best case)
            p_lower = 1 - stats.norm.cdf(test_statistic * gamma)

            results.append({
                'Gamma': gamma,
                'P_value_lower': p_lower,
                'P_value_upper': p_upper,
                'Significant_at_0.05': p_upper < 0.05
            })

        return pd.DataFrame(results)

    def tipping_point_analysis(
        self,
        n_simulations: int = 1000,
        confounder_prevalence_range: Tuple[float, float] = (0.1, 0.9),
        confounder_rr_range: Tuple[float, float] = (1.0, 5.0),
        random_state: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Tipping point analysis: find combinations that tip conclusion.

        Parameters
        ----------
        n_simulations : int
            Number of parameter combinations to test
        confounder_prevalence_range : tuple
            Range of confounder prevalence values
        confounder_rr_range : tuple
            Range of confounder RR values
        random_state : int, optional
            Random seed

        Returns
        -------
        pd.DataFrame
            Results showing which combinations tip the conclusion
        """
        rng = np.random.default_rng(random_state)

        # Generate random parameter combinations
        prevalences = rng.uniform(
            confounder_prevalence_range[0],
            confounder_prevalence_range[1],
            n_simulations
        )

        rr_confounders = rng.uniform(
            confounder_rr_range[0],
            confounder_rr_range[1],
            n_simulations
        )

        results = []
        uc = UnmeasuredConfounding(self.observed_rr, self.observed_ci)

        for i in range(n_simulations):
            prev = prevalences[i]
            rr_conf = rr_confounders[i]

            # Adjust for this confounding scenario
            adjusted_rr = uc.adjust_for_confounding(
                rr_confounder_exposure=rr_conf,
                rr_confounder_outcome=rr_conf,
                prevalence_confounder_unexposed=prev
            )

            # Check if conclusion is "tipped" (adjusted RR crosses null)
            tipped = (self.observed_rr > 1 and adjusted_rr < 1) or \
                     (self.observed_rr < 1 and adjusted_rr > 1)

            results.append({
                'Prevalence_confounder': prev,
                'RR_confounder': rr_conf,
                'Adjusted_RR': adjusted_rr,
                'Tipped': tipped,
                'Bias_factor': self.observed_rr / adjusted_rr
            })

        df = pd.DataFrame(results)

        # Add summary
        n_tipped = df['Tipped'].sum()
        pct_tipped = 100 * n_tipped / n_simulations

        print(f"Tipping Point Analysis Summary:")
        print(f"  Observed RR: {self.observed_rr:.3f}")
        print(f"  Simulations: {n_simulations}")
        print(f"  Tipped to null: {n_tipped} ({pct_tipped:.1f}%)")
        print(f"  Robust combinations: {n_simulations - n_tipped} ({100-pct_tipped:.1f}%)")

        return df

    def rule_out_values(
        self,
        effect_sizes_to_rule_out: Optional[List[float]] = None
    ) -> pd.DataFrame:
        """
        Calculate E-values needed to rule out specific effect sizes.

        Parameters
        ----------
        effect_sizes_to_rule_out : list of float, optional
            RR values to assess (default: [1.5, 2.0, 2.5, 3.0])

        Returns
        -------
        pd.DataFrame
            E-values needed to rule out each effect size
        """
        if effect_sizes_to_rule_out is None:
            effect_sizes_to_rule_out = [1.5, 2.0, 2.5, 3.0]

        results = []

        for target_rr in effect_sizes_to_rule_out:
            evalue = EValue.required_confounding_strength(
                observed_rr=self.observed_rr,
                true_rr=target_rr
            )

            results.append({
                'Target_RR': target_rr,
                'E_value_required': evalue,
                'Interpretation': (
                    f"To rule out RR={target_rr:.1f}, confounding must be "
                    f"weaker than RR={evalue:.2f} with both exposure and outcome"
                )
            })

        return pd.DataFrame(results)
