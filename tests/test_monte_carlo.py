"""Tests for Monte Carlo sensitivity analysis."""

import pytest
import numpy as np
from quantbias import MonteCarloAnalysis, SelectionBias, Confounding


class TestMonteCarloAnalysis:
    """Test Monte Carlo sensitivity analysis."""

    def test_monte_carlo_initialization(self):
        """Test initialization of Monte Carlo analysis."""
        bias = SelectionBias(
            s11=0.8,
            s10=0.9,
            s01=0.7,
            s00=0.6
        )

        mc = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=100,
            random_state=42
        )

        assert mc.observed_rr == 2.5
        assert mc.n_iterations == 100
        assert len(mc.bias_parameters) == 1

    def test_monte_carlo_run(self):
        """Test running Monte Carlo simulation."""
        bias = SelectionBias(
            s11=0.8,
            s10=0.9,
            s01=0.7,
            s00=0.6,
            s11_distribution=('beta', {'a': 8, 'b': 2}),
        )

        mc = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=1000,
            random_state=42
        )

        results = mc.run(show_progress=False)

        assert results is not None
        assert results.observed_rr == 2.5
        assert len(results.adjusted_rr_samples) == 1000
        assert results.median > 0
        assert results.mean > 0
        assert results.ci_lower < results.median
        assert results.ci_upper > results.median

    def test_monte_carlo_multiple_biases(self):
        """Test Monte Carlo with multiple bias sources."""
        biases = [
            SelectionBias(
                s11=0.85,
                s10=0.90,
                s01=0.80,
                s00=0.60
            ),
            Confounding(
                rr_confounder_outcome_unexposed=2.0,
                prevalence_confounder_unexposed=0.30,
                rr_confounder_exposure=1.5
            )
        ]

        mc = MonteCarloAnalysis(
            observed_rr=2.8,
            bias_parameters=biases,
            n_iterations=500,
            random_state=42
        )

        results = mc.run(show_progress=False)

        # Should have corrected for both biases
        # Note: Without distributions, all samples are identical (point estimates)
        # The correction direction depends on the specific bias parameters
        assert results.median != 2.8  # Some correction applied
        assert results.median > 0

    def test_monte_carlo_distribution_summary(self):
        """Test distribution summary."""
        bias = Confounding(
            rr_confounder_outcome_unexposed=2.0,
            prevalence_confounder_unexposed=0.3,
            rr_confounder_exposure=2.0
        )

        mc = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=500,
            random_state=42
        )

        results = mc.run(show_progress=False)
        summary = mc.get_distribution_summary()

        assert 'mean' in summary
        assert 'median' in summary
        assert 'ci_lower' in summary
        assert 'ci_upper' in summary
        assert 'std' in summary

    def test_probability_calculations(self):
        """Test probability calculations."""
        bias = Confounding(
            rr_confounder_outcome_unexposed=2.0,
            prevalence_confounder_unexposed=0.3,
            rr_confounder_exposure=2.0
        )

        mc = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=1000,
            random_state=42
        )

        results = mc.run(show_progress=False)

        # Test probability that true RR > threshold
        prob_gt_2 = mc.probability_rr_greater_than(2.0)
        assert 0 <= prob_gt_2 <= 1

        # Test probability that true RR < threshold
        prob_lt_3 = mc.probability_rr_less_than(3.0)
        assert 0 <= prob_lt_3 <= 1

    def test_get_samples(self):
        """Test getting samples from analysis."""
        bias = SelectionBias(s11=0.8, s10=0.9, s01=0.7, s00=0.6)

        mc = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=500,
            random_state=42
        )

        results = mc.run(show_progress=False)
        samples = mc.get_samples()

        assert len(samples) == 500
        assert all(s > 0 for s in samples)
        assert isinstance(samples, np.ndarray)

    def test_reproducibility(self):
        """Test that results are reproducible with same random state."""
        bias = SelectionBias(
            s11=0.8,
            s10=0.9,
            s01=0.7,
            s00=0.6,
            s11_distribution=('beta', {'a': 8, 'b': 2})
        )

        mc1 = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=500,
            random_state=42
        )

        mc2 = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=500,
            random_state=42
        )

        results1 = mc1.run(show_progress=False)
        results2 = mc2.run(show_progress=False)

        assert np.allclose(results1.adjusted_rr_samples, results2.adjusted_rr_samples)


class TestMonteCarloResult:
    """Test MonteCarloResult class."""

    def test_monte_carlo_result_summary(self):
        """Test Monte Carlo result summary."""
        bias = SelectionBias(s11=0.8, s10=0.9, s01=0.7, s00=0.6)

        mc = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=1000,
            random_state=42
        )

        results = mc.run(show_progress=False)

        # Test that result object has expected attributes
        assert hasattr(results, 'observed_rr')
        assert hasattr(results, 'mean')
        assert hasattr(results, 'median')
        assert hasattr(results, 'ci_lower')
        assert hasattr(results, 'ci_upper')
        assert hasattr(results, 'adjusted_rr_samples')

        # Test values are reasonable
        assert results.observed_rr == 2.5
        assert results.mean > 0
        assert results.median > 0
        # Note: Without distributions, all samples are identical so CI bounds = median
        assert results.ci_lower <= results.median <= results.ci_upper


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
