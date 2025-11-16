"""
Advanced bias analysis examples.

Demonstrates sophisticated analytical techniques for bias analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from quantbias import (
    SensitivityAnalysis,
    UnmeasuredConfounding,
    BiasCorrection,
    Confounding,
    SelectionBias,
)
from quantbias.visualization import (
    plot_sensitivity_contour,
    plot_tipping_point,
    plot_rosenbaum_bounds,
)


def example_1_tipping_point_analysis():
    """Example 1: Tipping point analysis."""
    print("=" * 70)
    print("EXAMPLE 1: Tipping Point Analysis")
    print("=" * 70)

    # Create sensitivity analysis
    analysis = SensitivityAnalysis(
        observed_rr=2.2,
        observed_ci=(1.6, 3.0)
    )

    # Run tipping point analysis
    tipping_results = analysis.tipping_point_analysis(
        n_simulations=1000,
        confounder_prevalence_range=(0.1, 0.9),
        confounder_rr_range=(1.0, 5.0),
        random_state=42
    )

    print("Tipping point analysis complete!")
    print()

    # Analyze results
    tipped = tipping_results[tipping_results['Tipped']]
    print(f"Scenarios that tipped the conclusion: {len(tipped)}/{len(tipping_results)}")
    print()

    # Find minimum confounding needed to tip
    if len(tipped) > 0:
        min_rr_to_tip = tipped['RR_confounder'].min()
        print(f"Minimum RR of confounder to tip: {min_rr_to_tip:.2f}")
        print()


def example_2_rosenbaum_bounds():
    """Example 2: Rosenbaum sensitivity bounds."""
    print("=" * 70)
    print("EXAMPLE 2: Rosenbaum Sensitivity Bounds")
    print("=" * 70)

    analysis = SensitivityAnalysis(
        observed_rr=1.8,
        observed_ci=(1.3, 2.5)
    )

    # Calculate Rosenbaum bounds
    bounds = analysis.rosenbaum_bounds(
        gamma_values=np.arange(1.0, 4.1, 0.2)
    )

    print("Rosenbaum Bounds (first 10 rows):")
    print(bounds.head(10))
    print()

    # Find tipping point
    tipping_gamma = bounds[~bounds['Significant_at_0.05']].iloc[0]['Gamma'] \
        if any(~bounds['Significant_at_0.05']) else None

    if tipping_gamma:
        print(f"Tipping point (Gamma where significance is lost): {tipping_gamma:.2f}")
        print(f"Interpretation: Hidden bias stronger than Gamma={tipping_gamma:.2f} could")
        print(f"  explain away the observed association.")
    else:
        print("Association remains significant across all tested Gamma values.")
    print()


def example_3_rule_out_analysis():
    """Example 3: Rule-out sensitivity values."""
    print("=" * 70)
    print("EXAMPLE 3: Rule-Out Analysis")
    print("=" * 70)

    analysis = SensitivityAnalysis(
        observed_rr=2.8,
        observed_ci=(2.1, 3.7)
    )

    # Calculate E-values to rule out specific effect sizes
    rule_out = analysis.rule_out_values(
        effect_sizes_to_rule_out=[1.5, 2.0, 2.5, 3.0]
    )

    print("E-values Required to Rule Out Different Effect Sizes:")
    print(rule_out[['Target_RR', 'E_value_required']])
    print()


def example_4_comprehensive_sensitivity_grid():
    """Example 4: Comprehensive sensitivity grid analysis."""
    print("=" * 70)
    print("EXAMPLE 4: Comprehensive Sensitivity Grid")
    print("=" * 70)

    analysis = UnmeasuredConfounding(
        observed_rr=2.5,
        observed_ci=(1.9, 3.3)
    )

    # Create comprehensive grid
    grid = analysis.sensitivity_grid(
        rr_confounder_outcome=np.linspace(1.0, 5.0, 20),
        prevalence_exposed=np.linspace(0.1, 0.9, 20)
    )

    print(f"Created sensitivity grid with {len(grid)} scenarios")
    print()

    # Analyze grid
    scenarios_near_null = grid[grid['Adjusted_RR'].between(0.8, 1.2)]
    scenarios_strong_effect = grid[grid['Adjusted_RR'] > 2.0]

    print(f"Scenarios producing near-null effect (RR 0.8-1.2): {len(scenarios_near_null)}")
    print(f"Scenarios maintaining strong effect (RR > 2.0): {len(scenarios_strong_effect)}")
    print()

    # Find typical scenario that eliminates effect
    if len(scenarios_near_null) > 0:
        typical_null = scenarios_near_null.iloc[len(scenarios_near_null) // 2]
        print("Typical confounding scenario that eliminates effect:")
        print(f"  RR confounder-outcome: {typical_null['RR_confounder_outcome']:.2f}")
        print(f"  Prevalence in exposed: {typical_null['Prevalence_exposed']:.2f}")
        print(f"  E-value: {typical_null['E_value']:.2f}")
    print()


def example_5_full_bias_report():
    """Example 5: Generate full bias analysis report."""
    print("=" * 70)
    print("EXAMPLE 5: Full Bias Analysis Report")
    print("=" * 70)

    # Define comprehensive bias scenario
    biases = [
        SelectionBias(
            sensitivity=0.85,
            specificity=0.88,
            selection_probability=0.65
        ),
        Confounding(
            rr_confounder_exposure=1.8,
            rr_confounder_outcome=2.0,
            prevalence_confounder_unexposed=0.25
        )
    ]

    # Create bias correction with full analysis
    correction = BiasCorrection(
        observed_rr=2.8,
        observed_ci=(2.1, 3.7),
        biases=biases
    )

    # Generate full report
    report = correction.full_report()
    print(report)
    print()


def example_6_sensitivity_visualization():
    """Example 6: Create sensitivity analysis visualizations."""
    print("=" * 70)
    print("EXAMPLE 6: Sensitivity Analysis Visualization")
    print("=" * 70)

    analysis = UnmeasuredConfounding(observed_rr=2.5)

    # Create sensitivity grid
    grid = analysis.sensitivity_grid(
        rr_confounder_outcome=np.linspace(1.0, 5.0, 30),
        prevalence_exposed=np.linspace(0.1, 0.9, 30)
    )

    # Create contour plot
    fig = plot_sensitivity_contour(grid, figsize=(10, 8))
    plt.savefig('/tmp/sensitivity_contour.png', dpi=300, bbox_inches='tight')
    print("Sensitivity contour plot saved to: /tmp/sensitivity_contour.png")
    print()


def example_7_threshold_analysis_detailed():
    """Example 7: Detailed threshold analysis."""
    print("=" * 70)
    print("EXAMPLE 7: Detailed Threshold Analysis")
    print("=" * 70)

    analysis = UnmeasuredConfounding(observed_rr=3.0)

    # Test multiple target RR values
    targets = [0.8, 1.0, 1.5, 2.0, 2.5]

    print("Confounding Strength Required for Different Target RRs:\n")
    print(f"{'Target RR':<12} {'Bias Factor':<15} {'Required RR Confounder':<25}")
    print("-" * 70)

    for target in targets:
        result = analysis.threshold_analysis(
            target_rr=target,
            prevalence_confounder=0.30
        )

        print(f"{target:<12.2f} {result['bias_factor_needed']:<15.3f} "
              f"{result['required_rr_confounder']:<25.3f}")

    print()


def example_8_comparative_scenarios():
    """Example 8: Compare different bias scenarios."""
    print("=" * 70)
    print("EXAMPLE 8: Comparative Bias Scenarios")
    print("=" * 70)

    observed_rr = 2.5

    scenarios = [
        {
            'name': 'Weak confounding',
            'bias': Confounding(
                rr_confounder_exposure=1.3,
                rr_confounder_outcome=1.5,
                prevalence_confounder_unexposed=0.20
            )
        },
        {
            'name': 'Moderate confounding',
            'bias': Confounding(
                rr_confounder_exposure=2.0,
                rr_confounder_outcome=2.0,
                prevalence_confounder_unexposed=0.30
            )
        },
        {
            'name': 'Strong confounding',
            'bias': Confounding(
                rr_confounder_exposure=3.0,
                rr_confounder_outcome=3.0,
                prevalence_confounder_unexposed=0.40
            )
        },
    ]

    print(f"Observed RR: {observed_rr:.3f}\n")
    print(f"{'Scenario':<25} {'Corrected RR':<15} {'Bias Factor':<15}")
    print("-" * 70)

    for scenario in scenarios:
        corrected_rr = scenario['bias'].apply_bias(observed_rr)
        bias_factor = observed_rr / corrected_rr

        print(f"{scenario['name']:<25} {corrected_rr:<15.3f} {bias_factor:<15.3f}")

    print()


if __name__ == "__main__":
    examples = [
        example_1_tipping_point_analysis,
        example_2_rosenbaum_bounds,
        example_3_rule_out_analysis,
        example_4_comprehensive_sensitivity_grid,
        example_5_full_bias_report,
        example_6_sensitivity_visualization,
        example_7_threshold_analysis_detailed,
        example_8_comparative_scenarios,
    ]

    for example_func in examples:
        try:
            example_func()
            print("\n")
        except Exception as e:
            print(f"Error in {example_func.__name__}: {e}")
            import traceback
            traceback.print_exc()
            print("\n")
