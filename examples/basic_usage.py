"""
Basic usage examples for quantbias package.

This script demonstrates fundamental bias analysis workflows.
"""

import numpy as np
from quantbias import (
    EValue,
    calculate_evalue,
    UnmeasuredConfounding,
    BiasCorrection,
    SelectionBias,
    MeasurementError,
    Confounding,
)


def example_1_evalue_calculation():
    """Example 1: Calculate E-values for an observed association."""
    print("=" * 70)
    print("EXAMPLE 1: E-value Calculation")
    print("=" * 70)

    # Observed effect from a study
    observed_rr = 2.5
    confidence_interval = (1.8, 3.5)

    # Calculate E-value
    result = calculate_evalue(
        observed_rr=observed_rr,
        confidence_interval=confidence_interval,
        effect_measure="RR"
    )

    print(result)
    print()


def example_2_sensitivity_analysis():
    """Example 2: Sensitivity analysis for unmeasured confounding."""
    print("=" * 70)
    print("EXAMPLE 2: Sensitivity Analysis for Unmeasured Confounding")
    print("=" * 70)

    # Observed association
    observed_rr = 2.0

    # Create sensitivity analysis
    analysis = UnmeasuredConfounding(observed_rr=observed_rr)

    # Test a specific confounding scenario
    adjusted_rr = analysis.adjust_for_confounding(
        rr_confounder_exposure=2.0,
        rr_confounder_outcome=2.0,
        prevalence_confounder_unexposed=0.3
    )

    print(f"Observed RR: {observed_rr:.3f}")
    print(f"Adjusted RR (assuming RR_conf = 2.0): {adjusted_rr:.3f}")
    print()

    # Threshold analysis
    threshold = analysis.threshold_analysis(target_rr=1.0, prevalence_confounder=0.3)
    print("Threshold Analysis:")
    print(f"  {threshold['interpretation']}")
    print()


def example_3_simple_bias_correction():
    """Example 3: Simple bias correction for selection bias."""
    print("=" * 70)
    print("EXAMPLE 3: Selection Bias Correction")
    print("=" * 70)

    # Define selection bias
    selection_bias = SelectionBias(
        sensitivity=0.80,  # 80% of true exposed+diseased are selected
        specificity=0.85,  # 85% of true unexposed+non-diseased are selected
        selection_probability=0.60
    )

    # Apply correction
    observed_rr = 2.5
    corrected_rr = selection_bias.apply_bias(observed_rr)

    print(f"Observed RR: {observed_rr:.3f}")
    print(f"Selection-corrected RR: {corrected_rr:.3f}")
    print(f"Bias factor: {observed_rr / corrected_rr:.3f}")
    print()


def example_4_multiple_bias_correction():
    """Example 4: Correcting for multiple bias sources."""
    print("=" * 70)
    print("EXAMPLE 4: Multiple Bias Correction")
    print("=" * 70)

    # Define multiple bias sources
    biases = [
        SelectionBias(
            sensitivity=0.85,
            specificity=0.90,
            selection_probability=0.60
        ),
        MeasurementError(
            sensitivity_exposure=0.95,
            specificity_exposure=0.90,
            sensitivity_outcome=0.98,
            specificity_outcome=0.95,
            nondifferential=True
        ),
        Confounding(
            rr_confounder_exposure=1.5,
            rr_confounder_outcome=2.0,
            prevalence_confounder_unexposed=0.30
        )
    ]

    # Apply comprehensive bias correction
    correction = BiasCorrection(
        observed_rr=2.5,
        observed_ci=(1.8, 3.5),
        biases=biases
    )

    # Deterministic correction
    result = correction.correct(method="deterministic")
    print(result)
    print()

    # Bias decomposition
    print("Bias Decomposition:")
    for bias_name, factor in result.bias_factors.items():
        print(f"  {bias_name}: {factor:.3f}")
    print()


def example_5_sensitivity_grid():
    """Example 5: Create sensitivity analysis grid."""
    print("=" * 70)
    print("EXAMPLE 5: Sensitivity Analysis Grid")
    print("=" * 70)

    analysis = UnmeasuredConfounding(observed_rr=2.0)

    # Create grid of confounding scenarios
    grid = analysis.sensitivity_grid(
        rr_confounder_outcome=np.arange(1.0, 4.1, 0.5),
        prevalence_exposed=np.arange(0.2, 0.8, 0.1)
    )

    print("Sensitivity Grid (first 10 rows):")
    print(grid[['RR_confounder_outcome', 'Prevalence_exposed', 'Adjusted_RR', 'E_value']].head(10))
    print()

    # Find scenarios that eliminate the association
    null_scenarios = grid[grid['Adjusted_RR'] < 1.1]
    print(f"Number of scenarios reducing RR to near-null: {len(null_scenarios)}/{len(grid)}")
    print()


def example_6_evalue_interpretation():
    """Example 6: E-value interpretation for different scenarios."""
    print("=" * 70)
    print("EXAMPLE 6: E-value Interpretation")
    print("=" * 70)

    scenarios = [
        {"name": "Strong effect", "rr": 3.0, "ci": (2.5, 3.8)},
        {"name": "Moderate effect", "rr": 1.8, "ci": (1.5, 2.2)},
        {"name": "Weak effect", "rr": 1.3, "ci": (1.1, 1.6)},
    ]

    for scenario in scenarios:
        result = calculate_evalue(
            observed_rr=scenario["rr"],
            confidence_interval=scenario["ci"]
        )

        print(f"{scenario['name']}:")
        print(f"  RR = {scenario['rr']:.2f} ({scenario['ci'][0]:.2f}, {scenario['ci'][1]:.2f})")
        print(f"  E-value (point): {result.point_estimate:.2f}")
        print(f"  E-value (CI): {result.ci_lower:.2f}")
        print()


def example_7_odds_ratio_conversion():
    """Example 7: E-values for odds ratios."""
    print("=" * 70)
    print("EXAMPLE 7: E-values for Odds Ratios")
    print("=" * 70)

    # For common outcomes, provide baseline risk
    observed_or = 2.5
    baseline_risk = 0.30  # 30% baseline risk

    # Calculate E-value with OR-to-RR conversion
    result = calculate_evalue(
        observed_rr=observed_or,
        effect_measure="OR",
        baseline_risk=baseline_risk
    )

    print(f"Observed OR: {observed_or:.2f}")
    print(f"Baseline risk: {baseline_risk:.2f}")
    print(f"Approximate RR: {observed_or:.2f} (using OR approximation)")
    print(f"E-value: {result.point_estimate:.2f}")
    print()


if __name__ == "__main__":
    # Run all examples
    examples = [
        example_1_evalue_calculation,
        example_2_sensitivity_analysis,
        example_3_simple_bias_correction,
        example_4_multiple_bias_correction,
        example_5_sensitivity_grid,
        example_6_evalue_interpretation,
        example_7_odds_ratio_conversion,
    ]

    for example_func in examples:
        example_func()
        print("\n")
