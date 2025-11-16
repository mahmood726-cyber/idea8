"""Tests for sensitivity analysis."""

import pytest
import numpy as np
import pandas as pd
from quantbias import UnmeasuredConfounding, SensitivityAnalysis


class TestUnmeasuredConfounding:
    """Test unmeasured confounding analysis."""

    def test_initialization(self):
        """Test initialization."""
        uc = UnmeasuredConfounding(observed_rr=2.0)
        assert uc.observed_rr == 2.0

    def test_adjust_for_confounding(self):
        """Test confounding adjustment."""
        uc = UnmeasuredConfounding(observed_rr=2.5)

        adjusted = uc.adjust_for_confounding(
            rr_confounder_exposure=2.0,
            rr_confounder_outcome=2.0,
            prevalence_confounder_unexposed=0.3
        )

        assert adjusted > 0
        assert adjusted != 2.5  # Should be different from observed

    def test_sensitivity_grid(self):
        """Test sensitivity grid generation."""
        uc = UnmeasuredConfounding(observed_rr=2.0)

        grid = uc.sensitivity_grid(
            rr_confounder_outcome=np.arange(1.0, 3.1, 0.5),
            prevalence_exposed=np.arange(0.2, 0.7, 0.1)
        )

        assert isinstance(grid, pd.DataFrame)
        assert 'Adjusted_RR' in grid.columns
        assert 'E_value' in grid.columns
        assert len(grid) > 0

    def test_threshold_analysis(self):
        """Test threshold analysis."""
        uc = UnmeasuredConfounding(observed_rr=2.5)

        threshold = uc.threshold_analysis(target_rr=1.0, prevalence_confounder=0.3)

        assert 'required_rr_confounder' in threshold
        assert 'bias_factor_needed' in threshold
        assert threshold['bias_factor_needed'] == pytest.approx(2.5, rel=0.01)


class TestSensitivityAnalysis:
    """Test comprehensive sensitivity analysis."""

    def test_initialization(self):
        """Test initialization."""
        sa = SensitivityAnalysis(observed_rr=2.0, observed_ci=(1.5, 2.7))
        assert sa.observed_rr == 2.0
        assert sa.observed_ci == (1.5, 2.7)

    def test_rosenbaum_bounds(self):
        """Test Rosenbaum bounds calculation."""
        sa = SensitivityAnalysis(observed_rr=2.0, observed_ci=(1.5, 2.7))

        bounds = sa.rosenbaum_bounds(gamma_values=np.arange(1.0, 3.1, 0.5))

        assert isinstance(bounds, pd.DataFrame)
        assert 'Gamma' in bounds.columns
        assert 'P_value_lower' in bounds.columns
        assert 'P_value_upper' in bounds.columns
        assert len(bounds) > 0

    def test_tipping_point_analysis(self):
        """Test tipping point analysis."""
        sa = SensitivityAnalysis(observed_rr=2.0)

        tipping = sa.tipping_point_analysis(
            n_simulations=100,
            random_state=42
        )

        assert isinstance(tipping, pd.DataFrame)
        assert 'Tipped' in tipping.columns
        assert 'Adjusted_RR' in tipping.columns
        assert len(tipping) == 100

    def test_rule_out_values(self):
        """Test rule-out analysis."""
        sa = SensitivityAnalysis(observed_rr=2.8)

        rule_out = sa.rule_out_values(effect_sizes_to_rule_out=[1.5, 2.0, 2.5])

        assert isinstance(rule_out, pd.DataFrame)
        assert 'Target_RR' in rule_out.columns
        assert 'E_value_required' in rule_out.columns
        assert len(rule_out) == 3


class TestSensitivityEdgeCases:
    """Test edge cases in sensitivity analysis."""

    def test_strong_confounding(self):
        """Test with very strong confounding."""
        uc = UnmeasuredConfounding(observed_rr=2.0)

        adjusted = uc.adjust_for_confounding(
            rr_confounder_exposure=5.0,
            rr_confounder_outcome=5.0,
            prevalence_confounder_unexposed=0.5
        )

        # Strong confounding should substantially reduce estimate
        assert adjusted < 2.0

    def test_weak_confounding(self):
        """Test with very weak confounding."""
        uc = UnmeasuredConfounding(observed_rr=2.0)

        adjusted = uc.adjust_for_confounding(
            rr_confounder_exposure=1.1,
            rr_confounder_outcome=1.1,
            prevalence_confounder_unexposed=0.3
        )

        # Weak confounding should barely change estimate
        assert np.isclose(adjusted, 2.0, rtol=0.1)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
