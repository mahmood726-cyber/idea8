"""Tests for E-value calculations."""

import pytest
import numpy as np
from quantbias import EValue, calculate_evalue
from quantbias.utils.exceptions import InvalidParameterError


class TestEValue:
    """Test E-value calculations."""

    def test_evalue_rr_basic(self):
        """Test basic E-value calculation for RR."""
        # RR = 2.0 should give E-value = 2 + sqrt(2*1) = 3.414
        evalue = EValue.calculate_evalue_rr(2.0)
        assert np.isclose(evalue, 3.414, atol=0.01)

    def test_evalue_rr_null(self):
        """Test E-value for null effect (RR = 1.0)."""
        evalue = EValue.calculate_evalue_rr(1.0)
        assert evalue == 1.0

    def test_evalue_protective(self):
        """Test E-value for protective effect (RR < 1)."""
        # Should invert RR < 1
        evalue = EValue.calculate_evalue_rr(0.5)
        expected_evalue = EValue.calculate_evalue_rr(2.0)
        assert np.isclose(evalue, expected_evalue, atol=0.001)

    def test_evalue_with_ci(self):
        """Test E-value calculation with confidence interval."""
        result = calculate_evalue(
            observed_rr=2.5,
            confidence_interval=(1.8, 3.5)
        )

        assert result.observed_effect == 2.5
        assert result.point_estimate > 1.0
        assert result.ci_lower > 1.0
        assert result.ci_lower < result.point_estimate

    def test_evalue_odds_ratio(self):
        """Test E-value for odds ratio."""
        # For rare outcomes, OR approximates RR
        evalue_or = EValue.calculate_evalue_or(2.0)
        evalue_rr = EValue.calculate_evalue_rr(2.0)
        assert np.isclose(evalue_or, evalue_rr, atol=0.01)

    def test_evalue_hazard_ratio(self):
        """Test E-value for hazard ratio."""
        evalue_hr = EValue.calculate_evalue_hr(2.0)
        evalue_rr = EValue.calculate_evalue_rr(2.0)
        assert np.isclose(evalue_hr, evalue_rr, atol=0.01)

    def test_invalid_rr(self):
        """Test that invalid RR raises error."""
        with pytest.raises(InvalidParameterError):
            EValue.calculate_evalue_rr(-1.0)

        with pytest.raises(InvalidParameterError):
            EValue.calculate_evalue_rr(0.0)

    def test_required_confounding_strength(self):
        """Test calculation of required confounding strength."""
        # To reduce RR=2.0 to null (1.0)
        required = EValue.required_confounding_strength(
            observed_rr=2.0,
            true_rr=1.0
        )
        assert required > 1.0

    def test_evalue_result_string(self):
        """Test string representation of EValue result."""
        result = calculate_evalue(
            observed_rr=2.0,
            confidence_interval=(1.5, 2.7)
        )
        string_repr = str(result)
        assert "E-value" in string_repr
        assert "2.0" in string_repr


class TestEValueEdgeCases:
    """Test edge cases for E-value calculations."""

    def test_very_large_rr(self):
        """Test E-value for very large RR."""
        evalue = EValue.calculate_evalue_rr(10.0)
        assert evalue > 10.0

    def test_very_small_protective_rr(self):
        """Test E-value for very small protective RR."""
        evalue = EValue.calculate_evalue_rr(0.1)
        # Should be equivalent to RR=10
        expected = EValue.calculate_evalue_rr(10.0)
        assert np.isclose(evalue, expected, atol=0.01)

    def test_rr_just_above_null(self):
        """Test E-value for RR just above null."""
        evalue = EValue.calculate_evalue_rr(1.01)
        assert evalue > 1.0
        assert evalue < 2.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
