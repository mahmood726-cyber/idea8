"""
Monte Carlo sensitivity analysis examples.

Demonstrates probabilistic bias analysis with uncertainty propagation.
"""

import numpy as np
import matplotlib.pyplot as plt
from quantbias import (
    MonteCarloAnalysis,
    SelectionBias,
    MeasurementError,
    Confounding,
)
from quantbias.visualization import plot_monte_carlo_distribution


def example_1_basic_monte_carlo():
    """Example 1: Basic Monte Carlo sensitivity analysis."""
    print("=" * 70)
    print("EXAMPLE 1: Basic Monte Carlo Analysis")
    print("=" * 70)

    # Define bias with uncertainty distributions
    selection_bias = SelectionBias(
        sensitivity=0.80,
        specificity=0.85,
        selection_probability=0.60,
        # Add uncertainty distributions
        sensitivity_distribution=('beta', {'a': 8, 'b': 2}),  # Beta(8, 2) ~ 80%
        specificity_distribution=('beta', {'a': 8.5, 'b': 1.5})  # Beta(8.5, 1.5) ~ 85%
    )

    # Run Monte Carlo
    mc = MonteCarloAnalysis(
        observed_rr=2.5,
        bias_parameters=[selection_bias],
        n_iterations=10000,
        random_state=42
    )

    results = mc.run(show_progress=True)
    print(results.summary())
    print()

    # Calculate probabilities
    prob_null = mc.probability_below_threshold(1.0)
    prob_strong = mc.probability_above_threshold(2.0)

    print(f"Probability corrected RR < 1.0: {prob_null:.3f}")
    print(f"Probability corrected RR > 2.0: {prob_strong:.3f}")
    print()


def example_2_multiple_biases_monte_carlo():
    """Example 2: Monte Carlo with multiple bias sources."""
    print("=" * 70)
    print("EXAMPLE 2: Multiple Biases Monte Carlo")
    print("=" * 70)

    # Define multiple biases with uncertainty
    biases = [
        SelectionBias(
            sensitivity=0.85,
            specificity=0.90,
            selection_probability=0.60,
            sensitivity_distribution=('beta', {'a': 17, 'b': 3}),
            specificity_distribution=('beta', {'a': 18, 'b': 2})
        ),
        MeasurementError(
            sensitivity_exposure=0.95,
            specificity_exposure=0.90,
            sensitivity_outcome=0.98,
            specificity_outcome=0.95,
            nondifferential=True,
            sensitivity_exposure_dist=('beta', {'a': 19, 'b': 1}),
            specificity_exposure_dist=('beta', {'a': 18, 'b': 2}),
            sensitivity_outcome_dist=('beta', {'a': 49, 'b': 1}),
            specificity_outcome_dist=('beta', {'a': 19, 'b': 1})
        ),
        Confounding(
            rr_confounder_exposure=1.8,
            rr_confounder_outcome=2.2,
            prevalence_confounder_unexposed=0.25,
            # Uncertainty in confounding strength
            rr_conf_exp_dist=('lognormal', {'mean': np.log(1.8), 'sigma': 0.2}),
            rr_conf_out_dist=('lognormal', {'mean': np.log(2.2), 'sigma': 0.2}),
            prev_conf_unexp_dist=('beta', {'a': 5, 'b': 15})
        )
    ]

    mc = MonteCarloAnalysis(
        observed_rr=2.8,
        bias_parameters=biases,
        n_iterations=10000,
        random_state=42
    )

    results = mc.run(show_progress=True)
    print(results.summary())
    print()

    # Distribution summary
    summary = mc.get_distribution_summary()
    print("Distribution Summary:")
    print(summary)
    print()


def example_3_visualization():
    """Example 3: Visualize Monte Carlo results."""
    print("=" * 70)
    print("EXAMPLE 3: Monte Carlo Visualization")
    print("=" * 70)

    # Simple bias scenario
    bias = Confounding(
        rr_confounder_exposure=2.0,
        rr_confounder_outcome=2.5,
        prevalence_confounder_unexposed=0.30,
        rr_conf_exp_dist=('lognormal', {'mean': np.log(2.0), 'sigma': 0.3}),
        rr_conf_out_dist=('lognormal', {'mean': np.log(2.5), 'sigma': 0.3}),
        prev_conf_unexp_dist=('beta', {'a': 6, 'b': 14})
    )

    mc = MonteCarloAnalysis(
        observed_rr=3.0,
        bias_parameters=[bias],
        n_iterations=10000,
        random_state=42
    )

    results = mc.run()

    # Create visualization
    fig = plot_monte_carlo_distribution(results)
    plt.savefig('/tmp/monte_carlo_example.png', dpi=300, bbox_inches='tight')
    print("Plot saved to: /tmp/monte_carlo_example.png")
    print()


def example_4_sensitivity_to_parameters():
    """Example 4: Sensitivity to different parameter assumptions."""
    print("=" * 70)
    print("EXAMPLE 4: Sensitivity to Parameter Assumptions")
    print("=" * 70)

    # Test different confounding strengths
    rr_confounder_values = [1.5, 2.0, 2.5, 3.0]

    print("Comparison of different confounding strengths:\n")
    print(f"{'RR Conf':<10} {'Median':<10} {'95% CI':<25} {'P(RR<1.0)':<12}")
    print("-" * 70)

    for rr_conf in rr_confounder_values:
        bias = Confounding(
            rr_confounder_exposure=rr_conf,
            rr_confounder_outcome=rr_conf,
            prevalence_confounder_unexposed=0.30,
            rr_conf_exp_dist=('lognormal', {'mean': np.log(rr_conf), 'sigma': 0.2}),
            rr_conf_out_dist=('lognormal', {'mean': np.log(rr_conf), 'sigma': 0.2}),
        )

        mc = MonteCarloAnalysis(
            observed_rr=2.5,
            bias_parameters=[bias],
            n_iterations=5000,
            random_state=42
        )

        results = mc.run(show_progress=False)
        prob_null = mc.probability_below_threshold(1.0)

        print(f"{rr_conf:<10.1f} {results.median:<10.3f} "
              f"({results.ci_lower:.3f}, {results.ci_upper:.3f}){'':<8} {prob_null:<12.3f}")

    print()


def example_5_uncertainty_quantification():
    """Example 5: Quantifying total uncertainty."""
    print("=" * 70)
    print("EXAMPLE 5: Uncertainty Quantification")
    print("=" * 70)

    bias = SelectionBias(
        sensitivity=0.80,
        specificity=0.85,
        selection_probability=0.60,
        sensitivity_distribution=('beta', {'a': 8, 'b': 2}),
        specificity_distribution=('beta', {'a': 8.5, 'b': 1.5})
    )

    mc = MonteCarloAnalysis(
        observed_rr=2.0,
        bias_parameters=[bias],
        n_iterations=10000,
        random_state=42
    )

    results = mc.run(show_progress=False)

    # Analyze uncertainty
    samples = results.adjusted_rr_samples

    print("Uncertainty Analysis:")
    print(f"  Observed RR: {results.observed_rr:.3f}")
    print(f"  Median corrected RR: {results.median:.3f}")
    print(f"  Mean corrected RR: {results.mean:.3f}")
    print(f"  Standard deviation: {np.std(samples):.3f}")
    print(f"  Coefficient of variation: {np.std(samples) / results.mean:.3f}")
    print()

    # Percentiles
    print("Percentiles of corrected RR:")
    for p in [5, 25, 50, 75, 95]:
        value = np.percentile(samples, p)
        print(f"  {p}th percentile: {value:.3f}")
    print()

    # Range analysis
    range_width = results.ci_upper - results.ci_lower
    relative_range = range_width / results.median

    print(f"95% Interval width: {range_width:.3f}")
    print(f"Relative interval width: {relative_range:.3f}")
    print()


if __name__ == "__main__":
    examples = [
        example_1_basic_monte_carlo,
        example_2_multiple_biases_monte_carlo,
        example_3_visualization,
        example_4_sensitivity_to_parameters,
        example_5_uncertainty_quantification,
    ]

    for example_func in examples:
        try:
            example_func()
            print("\n")
        except Exception as e:
            print(f"Error in {example_func.__name__}: {e}")
            print("\n")
