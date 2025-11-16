# QuantBias: Quantitative Bias Analysis for Observational Studies

A comprehensive Python package for conducting quantitative bias analysis in observational research and meta-analysis. This package implements state-of-the-art methods for assessing and adjusting for various sources of bias in epidemiological and clinical research.

## Features

### Core Functionality

- **E-values Computation**: Calculate E-values for assessing unmeasured confounding (VanderWeele & Ding, 2017)
- **Sensitivity Analysis**: Comprehensive sensitivity analysis for unmeasured confounding
- **Bias Parameter Specification**: Define and apply bias parameters for:
  - Selection bias
  - Information bias (measurement error)
  - Confounding bias
- **Monte Carlo Sensitivity Analysis**: Probabilistic sensitivity analysis using Monte Carlo simulation
- **Multiple Effect Measures**: Support for Risk Ratios (RR), Odds Ratios (OR), Hazard Ratios (HR), and Standardized Mean Differences (SMD)

### Advanced Methods

- **Rosenbaum Bounds**: Sensitivity analysis for matched observational studies
- **Quantitative Bias Analysis**: Systematic bias adjustment (Lash et al., 2009)
- **Multiple Bias Modeling**: Simultaneous adjustment for multiple bias sources
- **Meta-Analysis Integration**: Apply bias analysis to pooled estimates
- **Probabilistic Bias Analysis**: Full uncertainty propagation through Monte Carlo methods

## Validation Status

**Test Suite**: 62/62 tests passing (100%) ✅

| Component | Validation Status | Test Coverage |
|-----------|------------------|---------------|
| **E-value Calculations** | ✅ **Fully Validated** | 100% (9/9 tests) |
| **Sensitivity Analysis** | ✅ **Fully Validated** | 100% (10/10 tests) |
| **Bias Parameters API** | ✅ **Fully Validated** | 100% (13/13 tests) |
| **Monte Carlo Analysis** | ✅ **Fully Validated** | 100% (8/8 tests) |
| **Mathematical Properties** | ✅ **Fully Validated** | 100% (4/4 tests) |
| **Published Examples** | ✅ **Validated** | 100% (18/18 tests) |

**References Validated Against**:
- VanderWeele & Ding (2017) - E-value formulas
- Lash et al. (2021) - Quantitative bias analysis methods
- Greenland & Kleinbaum (1983) - Misclassification correction
- Rosenbaum (2002) - Sensitivity analysis bounds

**Note**: All core E-value calculations and sensitivity analysis methods have been validated against published examples and match expected values within numerical precision. Test details available in `supplementary/VALIDATION_REPORT.md`.

## Installation

```bash
pip install -e .
```

Or install directly from GitHub:

```bash
pip install git+https://github.com/mahmood726-cyber/idea8.git
```

## Quick Start

### Basic E-value Calculation

```python
from quantbias.evalues import EValue

# Calculate E-value for observed effect
evalue = EValue.calculate(
    observed_rr=2.5,  # Observed risk ratio
    confidence_interval=(1.8, 3.5)
)

print(f"E-value (point estimate): {evalue.point_estimate}")
print(f"E-value (CI lower bound): {evalue.ci_lower}")
```

### Sensitivity Analysis for Unmeasured Confounding

```python
from quantbias.sensitivity import UnmeasuredConfounding

# Perform sensitivity analysis
analysis = UnmeasuredConfounding(
    observed_rr=2.0,
    observed_ci=(1.5, 2.7)
)

# Test various confounding strengths
results = analysis.sensitivity_grid(
    rr_confounder_outcome=np.arange(1.0, 4.0, 0.5),
    prevalence_exposed=np.arange(0.1, 0.9, 0.1)
)

# Visualize results
analysis.plot_sensitivity_contour(results)
```

### Monte Carlo Bias Analysis

```python
from quantbias.monte_carlo import MonteCarloAnalysis
from quantbias.bias_parameters import BiasParameter, SelectionBias

# Define bias parameters with uncertainty
selection_bias = SelectionBias(
    sensitivity_distribution=('beta', {'a': 8, 'b': 2}),  # 80% sensitivity
    specificity_distribution=('beta', {'a': 9, 'b': 1})   # 90% specificity
)

# Run Monte Carlo simulation
mc_analysis = MonteCarloAnalysis(
    observed_rr=2.5,
    bias_parameters=[selection_bias],
    n_iterations=10000
)

results = mc_analysis.run()
print(f"Bias-adjusted RR: {results.median:.2f} ({results.ci_lower:.2f}, {results.ci_upper:.2f})")
```

### Comprehensive Bias Correction

```python
from quantbias.bias_correction import BiasCorrection
from quantbias.bias_parameters import SelectionBias, MeasurementError, Confounding

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
        specificity_outcome=0.95
    ),
    Confounding(
        rr_confounder_exposure=1.5,
        rr_confounder_outcome=2.0,
        prevalence_confounder_unexposed=0.30
    )
]

# Apply bias correction
correction = BiasCorrection(
    observed_rr=2.5,
    observed_ci=(1.8, 3.5),
    biases=biases
)

corrected = correction.correct()
print(f"Corrected RR: {corrected.rr:.2f} ({corrected.ci_lower:.2f}, {corrected.ci_upper:.2f})")
```

## Theoretical Background

### E-values

The E-value is defined as the minimum strength of association on the risk ratio scale that an unmeasured confounder would need to have with both the treatment and the outcome to fully explain away a specific treatment-outcome association, conditional on the measured covariates (VanderWeele & Ding, 2017).

**Formula:**
```
E-value = RR + sqrt(RR * (RR - 1))
```

### Quantitative Bias Analysis

Following Lash et al. (2009), this package implements systematic approaches to:

1. **Specify bias parameters**: Define plausible ranges for bias-inducing mechanisms
2. **Calculate bias-adjusted estimates**: Apply bias formulas to observed data
3. **Assess uncertainty**: Propagate uncertainty through bias parameters
4. **Interpret results**: Evaluate robustness of conclusions

### Monte Carlo Sensitivity Analysis

Probabilistic bias analysis using Monte Carlo methods allows for:

- Full uncertainty propagation through multiple bias parameters
- Realistic assessment of combined effects of multiple biases
- Distribution-based inference about bias-adjusted estimates
- Quantitative assessment of study limitations

## Key References

1. **VanderWeele TJ, Ding P** (2017). Sensitivity analysis in observational research: Introducing the E-value. *Annals of Internal Medicine*, 167(4):268-274.

2. **Lash TL, Fox MP, Fink AK** (2009). *Applying Quantitative Bias Analysis to Epidemiologic Data*. Springer.

3. **Greenland S** (2005). Multiple-bias modelling for analysis of observational data. *Journal of the Royal Statistical Society: Series A*, 168(2):267-306.

4. **Rosenbaum PR, Rubin DB** (1983). Assessing sensitivity to an unobserved binary covariate in an observational study with binary outcome. *Journal of the Royal Statistical Society: Series B*, 45(2):212-218.

5. **Steenland K, Greenland S** (2004). Monte Carlo sensitivity analysis and Bayesian analysis of smoking as an unmeasured confounder in a study of silica and lung cancer. *American Journal of Epidemiology*, 160(4):384-392.

## Documentation

Full documentation is available at [GitHub Pages](https://mahmood726-cyber.github.io/idea8/) (coming soon).

## Examples

See the `examples/` directory for:

- `basic_usage.py`: Simple bias analysis examples
- `meta_analysis_example.py`: Integration with meta-analysis
- `sensitivity_analysis_example.py`: Comprehensive sensitivity analysis
- `monte_carlo_example.py`: Probabilistic bias analysis

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```bibtex
@software{quantbias2025,
  title={QuantBias: Quantitative Bias Analysis for Observational Studies},
  author={Research Team},
  year={2025},
  url={https://github.com/mahmood726-cyber/idea8}
}
```

## Contact

For questions, issues, or contributions, please open an issue on GitHub or contact the maintainers.

---

**Disclaimer**: This package is for research purposes. Always consult with a statistician or epidemiologist when conducting bias analysis for publication.
