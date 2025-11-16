"""Tests for Monte Carlo sensitivity analysis."""

import pytest
import numpy as np
from quantbias import MonteCarloAnalysis, SelectionBias, Confounding


class TestMonteCarloAnalysis:
    """Test Monte Carlo sensitivity analysis."""

    def test_monte_carlo_initialization(self):
        """Test initialization of Monte Carlo analysis."""
        bias = SelectionBias(
            sensitivity=0.8,
            specificity=0.9,
            selection_probability=0.6
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
            sensitivity=0.8,
            specificity=0.9,
            selection_probability=0.6,
            sensitivity_distribution=('beta', {'a': 8, 'b': 2}),
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
                sensitivity=0.85,
                specificity=0.90,
                selection_probability=0.60
            ),
            Confounding(
                rr_confounder_exposure=1.5,
                rr_confounder_outcome=2.0,
                prevalence_confounder_unexposed=0.30
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
        assert results.median < 2.8  # Some correction expected
        assert results.median > 0

    def test_monte_carlo_distribution_summary(self):
        """Test distribution summary."""
        bias = Confounding(
            rr_confounder_exposure=2.0,
            rr_confounder_outcome=2.0,
            prevalence_confounder_unexposed=0.3
        )

        mc = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=500,
            random_state=42
        )

        results = mc.run(show_progress=False)
        summary = mc.get_distribution_summary()

        assert 'Statistic' in summary.columns
        assert 'Value' in summary.columns
        assert len(summary) > 0

    def test_probability_calculations(self):
        """Test probability calculations."""
        bias = Confounding(
            rr_confounder_exposure=3.0,
            rr_confounder_outcome=3.0,
            prevalence_confounder_unexposed=0.4
        )

        mc = MonteCarloAnalysis(
            observed_rr=2.0,
            bias_parameters=[bias],
            n_iterations=1000,
            random_state=42
        )

        mc.run(show_progress=False)

        prob_below_null = mc.probability_below_threshold(1.0)
        prob_above_null = mc.probability_above_threshold(1.0)

        assert 0 <= prob_below_null <= 1
        assert 0 <= prob_above_null <= 1
        assert np.isclose(prob_below_null + prob_above_null, 1.0, atol=0.01)

    def test_get_samples(self):
        """Test getting samples array."""
        bias = SelectionBias(sensitivity=0.8, specificity=0.9, selection_probability=0.6)

        mc = MonteCarloAnalysis(
            observed_rr=2.0,
            bias_parameters=[bias],
            n_iterations=100,
            random_state=42
        )

        mc.run(show_progress=False)
        samples = mc.get_samples()

        assert len(samples) == 100
        assert all(s > 0 for s in samples)

    def test_reproducibility(self):
        """Test that results are reproducible with same seed."""
        bias = SelectionBias(
            sensitivity=0.8,
            specificity=0.9,
            selection_probability=0.6,
            sensitivity_distribution=('beta', {'a': 8, 'b': 2}),
        )

        mc1 = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=100,
            random_state=42
        )

        mc2 = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=100,
            random_state=42
        )

        results1 = mc1.run(show_progress=False)
        results2 = mc2.run(show_progress=False)

        assert np.allclose(results1.adjusted_rr_samples, results2.adjusted_rr_samples)


class TestMonteCarloResult:
    """Test Monte Carlo result object."""

    def test_monte_carlo_result_summary(self):
        """Test summary string generation."""
        bias = SelectionBias(sensitivity=0.8, specificity=0.9, selection_probability=0.6)

        mc = MonteCarloAnalysis(
            observed_rr=2.0,
            bias_parameters=[bias],
            n_iterations=100,
            random_state=42
        )

        results = mc.run(show_progress=False)
        summary = results.summary()

        assert isinstance(summary, str)
        assert "Observed RR" in summary
        assert "Adjusted RR" in summary


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
