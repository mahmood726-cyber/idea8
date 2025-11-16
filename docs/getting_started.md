# Getting Started with QuantBias

This guide will help you get started with the QuantBias package for quantitative bias analysis.

## Installation

```bash
pip install -e .
```

## Quick Start

### 1. E-value Calculation

The E-value quantifies the minimum strength of unmeasured confounding needed to explain away an observed association.

```python
from quantbias import calculate_evalue

# Calculate E-value for an observed risk ratio
result = calculate_evalue(
    observed_rr=2.5,
    confidence_interval=(1.8, 3.5)
)

print(result)
# Output shows E-value for point estimate and CI
```

### 2. Sensitivity Analysis

Assess how unmeasured confounding could affect your results:

```python
from quantbias import UnmeasuredConfounding
import numpy as np

# Create analysis
analysis = UnmeasuredConfounding(observed_rr=2.0)

# Test specific confounding scenario
adjusted_rr = analysis.adjust_for_confounding(
    rr_confounder_exposure=2.0,
    rr_confounder_outcome=2.0,
    prevalence_confounder_unexposed=0.3
)

print(f"Adjusted RR: {adjusted_rr:.3f}")
```

### 3. Bias Correction

Correct for multiple known biases:

```python
from quantbias import BiasCorrection, SelectionBias, MeasurementError, Confounding

# Define bias sources
biases = [
    SelectionBias(sensitivity=0.85, specificity=0.90),
    MeasurementError(sensitivity_exposure=0.95, specificity_exposure=0.90),
    Confounding(rr_confounder_exposure=1.5, rr_confounder_outcome=2.0)
]

# Apply correction
correction = BiasCorrection(
    observed_rr=2.5,
    observed_ci=(1.8, 3.5),
    biases=biases
)

result = correction.correct()
print(result)
```

### 4. Monte Carlo Sensitivity Analysis

Propagate uncertainty through bias parameters:

```python
from quantbias import MonteCarloAnalysis, SelectionBias

# Define bias with uncertainty
bias = SelectionBias(
    sensitivity=0.80,
    specificity=0.85,
    sensitivity_distribution=('beta', {'a': 8, 'b': 2}),
    specificity_distribution=('beta', {'a': 8.5, 'b': 1.5})
)

# Run Monte Carlo simulation
mc = MonteCarloAnalysis(
    observed_rr=2.5,
    bias_parameters=[bias],
    n_iterations=10000
)

results = mc.run()
print(results.summary())
```

## Key Concepts

### E-values

The E-value is the minimum strength of association that an unmeasured confounder would need to have with both the exposure and outcome to explain away the observed association.

- **Point estimate E-value**: Confounding needed to reduce the point estimate to null
- **CI E-value**: Confounding needed to shift the CI to include null

### Bias Parameters

Three main types of bias:

1. **Selection Bias**: Differential selection into the study
2. **Information Bias (Measurement Error)**: Misclassification of exposure or outcome
3. **Confounding**: Unmeasured confounders

### Sensitivity Analysis

Techniques to assess robustness:

- **Sensitivity grids**: Test multiple parameter combinations
- **Tipping point analysis**: Find values that change conclusions
- **Rosenbaum bounds**: Formal sensitivity analysis for matched studies

### Monte Carlo Analysis

Probabilistic bias analysis:

- Specify uncertainty distributions for bias parameters
- Sample from distributions
- Generate distribution of bias-corrected estimates

## Next Steps

- See `examples/` directory for detailed examples
- Read methodology documentation in `docs/methodology.md`
- Check API reference in `docs/api_reference.md`

## References

1. VanderWeele TJ, Ding P. Sensitivity Analysis in Observational Research: Introducing the E-Value. Ann Intern Med. 2017;167(4):268-274.

2. Lash TL, Fox MP, Fink AK. Applying Quantitative Bias Analysis to Epidemiologic Data. Springer, 2009.

3. Greenland S. Multiple-bias modelling for analysis of observational data. J R Stat Soc Ser A. 2005;168(2):267-306.
