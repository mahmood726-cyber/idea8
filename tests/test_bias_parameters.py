"""Tests for bias parameters."""

import pytest
import numpy as np
from quantbias.bias_parameters import (
    SelectionBias,
    MeasurementError,
    Confounding,
)
from quantbias.utils.exceptions import InvalidParameterError


class TestSelectionBias:
    """Test selection bias parameters."""

    def test_selection_bias_initialization(self):
        """Test initialization of selection bias."""
        bias = SelectionBias(
            s11=0.8,
            s10=0.9,
            s01=0.7,
            s00=0.6
        )
        assert bias.s11 == 0.8
        assert bias.s10 == 0.9
        assert bias.s01 == 0.7
        assert bias.s00 == 0.6

    def test_selection_bias_apply(self):
        """Test applying selection bias correction."""
        bias = SelectionBias(
            s11=0.8,
            s10=0.8,
            s01=0.7,
            s00=0.5
        )
        observed_rr = 2.0
        corrected_rr = bias.apply_bias(observed_rr)

        # Corrected RR should differ from observed
        assert corrected_rr != observed_rr
        assert corrected_rr > 0

    def test_selection_bias_validation(self):
        """Test validation of selection bias parameters."""
        with pytest.raises(InvalidParameterError):
            SelectionBias(
                s11=1.5,  # Invalid: > 1
                s10=0.9,
                s01=0.7,
                s00=0.6
            )

        with pytest.raises(InvalidParameterError):
            SelectionBias(
                s11=0.8,
                s10=-0.1,  # Invalid: < 0
                s01=0.7,
                s00=0.6
            )

    def test_selection_bias_sampling(self):
        """Test parameter sampling for Monte Carlo."""
        bias = SelectionBias(
            s11=0.8,
            s10=0.9,
            s01=0.7,
            s00=0.6,
            s11_distribution=('beta', {'a': 8, 'b': 2}),
            s10_distribution=('beta', {'a': 9, 'b': 1})
        )

        samples = bias.sample_parameters(n=100, random_state=42)

        assert 's11' in samples
        assert 's10' in samples
        assert len(samples['s11']) == 100
        assert all(0 <= s <= 1 for s in samples['s11'])


class TestMeasurementError:
    """Test measurement error parameters."""

    def test_measurement_error_initialization(self):
        """Test initialization of measurement error."""
        bias = MeasurementError(
            sensitivity_exposure=0.95,
            specificity_exposure=0.90,
            sensitivity_outcome=0.98,
            specificity_outcome=0.95
        )
        assert bias.sensitivity_exposure == 0.95

    def test_measurement_error_apply(self):
        """Test applying measurement error correction."""
        bias = MeasurementError(
            sensitivity_exposure=0.90,
            specificity_exposure=0.90,
            sensitivity_outcome=0.95,
            specificity_outcome=0.95,
            nondifferential=True
        )
        observed_rr = 2.0
        corrected_rr = bias.apply_bias(observed_rr)

        # Non-differential misclassification biases toward null
        # So corrected should be further from null than observed
        assert corrected_rr > observed_rr

    def test_measurement_error_validation(self):
        """Test validation of measurement error parameters."""
        with pytest.raises(InvalidParameterError):
            MeasurementError(
                sensitivity_exposure=1.2,  # Invalid
                specificity_exposure=0.90,
                sensitivity_outcome=0.98,
                specificity_outcome=0.95
            )

    def test_measurement_error_sampling(self):
        """Test parameter sampling."""
        bias = MeasurementError(
            sensitivity_exposure=0.95,
            specificity_exposure=0.90,
            sensitivity_outcome=0.98,
            specificity_outcome=0.95,
            sensitivity_exposure_dist=('beta', {'a': 19, 'b': 1}),
        )

        samples = bias.sample_parameters(n=50, random_state=42)

        assert 'sensitivity_exposure' in samples
        assert len(samples['sensitivity_exposure']) == 50


class TestConfounding:
    """Test confounding parameters."""

    def test_confounding_initialization(self):
        """Test initialization of confounding bias."""
        bias = Confounding(
            rr_confounder_outcome_unexposed=2.0,
            prevalence_confounder_unexposed=0.3,
            rr_confounder_exposure=2.0
        )
        assert bias.rr_confounder_exposure == 2.0
        assert bias.rr_confounder_outcome_unexposed == 2.0
        assert bias.prevalence_confounder_unexposed == 0.3

    def test_confounding_apply(self):
        """Test applying confounding correction."""
        bias = Confounding(
            rr_confounder_outcome_unexposed=2.0,
            prevalence_confounder_unexposed=0.3,
            rr_confounder_exposure=2.0
        )
        observed_rr = 2.5
        corrected_rr = bias.apply_bias(observed_rr)

        # Confounding should reduce the observed effect
        assert corrected_rr < observed_rr
        assert corrected_rr > 0

    def test_confounding_validation(self):
        """Test validation of confounding parameters."""
        with pytest.raises(InvalidParameterError):
            Confounding(
                rr_confounder_outcome_unexposed=2.0,
                rr_confounder_exposure=-1.0,  # Invalid
                prevalence_confounder_unexposed=0.3
            )

        with pytest.raises(InvalidParameterError):
            Confounding(
                rr_confounder_outcome_unexposed=2.0,
                rr_confounder_exposure=2.0,
                prevalence_confounder_unexposed=1.5  # Invalid: > 1
            )

    def test_confounding_sampling(self):
        """Test parameter sampling."""
        bias = Confounding(
            rr_confounder_outcome_unexposed=2.0,
            prevalence_confounder_unexposed=0.3,
            rr_confounder_exposure=2.0,
            rr_conf_out_unexp_dist=('lognormal', {'mean': np.log(2.0), 'sigma': 0.2}),
        )

        samples = bias.sample_parameters(n=100, random_state=42)

        assert 'rr_confounder_outcome_unexposed' in samples
        assert len(samples['rr_confounder_outcome_unexposed']) == 100
        assert all(rr > 0 for rr in samples['rr_confounder_outcome_unexposed'])

    def test_confounding_prevalence_calculation(self):
        """Test automatic prevalence calculation."""
        bias = Confounding(
            rr_confounder_outcome_unexposed=2.0,
            prevalence_confounder_unexposed=0.3,
            rr_confounder_exposure=2.0,
            prevalence_confounder_exposed=None  # Should be calculated
        )

        # Should have calculated prevalence in exposed
        assert bias.prevalence_confounder_exposed is not None
        assert 0 < bias.prevalence_confounder_exposed < 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
