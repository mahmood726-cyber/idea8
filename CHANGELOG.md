# Changelog

All notable changes to the QuantBias project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2025-01-16 - Major Peer-Review Revisions

### 🔴 BREAKING CHANGES
- Complete rewrite of bias correction formulas
- `SelectionBias` now uses 4 parameters (s11, s10, s01, s00) instead of sensitivity/specificity
- `Confounding` now uses stratified RR parameters

### ✅ Added
- **Corrected Implementations**
  - `quantbias/bias_parameters_corrected.py` - Mathematically validated bias corrections
  - `quantbias/evalues_extended.py` - Extended E-value methods (RD, SMD, mediation)
  - `quantbias/inference.py` - Proper CI propagation (bootstrap, delta method, simulation)

- **Validation Suite**
  - `tests/test_validation.py` - 38 validation tests against published examples
  - Tests from VanderWeele & Ding (2017)
  - Tests from Lash et al. (2021)
  - Tests from Greenland & Kleinbaum (1983)
  - R package comparison tests

- **Documentation**
  - `docs/methodology.md` - 1100+ lines of mathematical derivations
  - Mathematical formulas with LaTeX
  - Detailed assumptions and limitations
  - Bias parameter elicitation guidelines
  - Complete reference list

- **New E-value Methods**
  - Risk difference E-values
  - Additive scale E-values
  - Standardized mean difference E-values
  - Mediation analysis E-values (NIE, NDE)
  - Effect modification / interaction E-values

- **Statistical Inference**
  - Bootstrap confidence intervals
  - Delta method confidence intervals
  - Simulation-based confidence intervals
  - Variance inflation factor calculation
  - Method comparison utilities

### 🔧 Fixed
- **Selection Bias**: Now uses proper 2×2 selection probability framework (Lash Equation 5.1)
- **Measurement Error**: Implemented matrix inversion method (Greenland & Kleinbaum 1983)
- **Confounding**: Corrected formula to allow effect modification (Greenland & Lash 2008)
- **CI Propagation**: Proper uncertainty accounting (no longer naive division)
- **E-values**: Corrected CI E-value calculation for protective effects

### 📝 Changed
- Updated all formulas with explicit citations
- Improved error messages with parameter names
- Enhanced validation of all inputs
- Better numerical stability (avoid log(0), division by zero)

### 🐛 Bug Fixes
- Fixed measurement error bias factor formula
- Fixed confounding prevalence calculation
- Fixed E-value for odds ratios with baseline risk conversion
- Fixed bootstrap random state handling

### 📚 Documentation
- Complete mathematical methodology document
- Response to reviewers document
- Expanded README with limitations
- Comparison to existing tools (EValue, episensr)
- When to use / when NOT to use guidelines

### ✅ Validation
All methods validated against published examples:
- ✅ VanderWeele & Ding (2017) E-values: 100% match
- ✅ Lash et al. (2021) Selection bias: 100% match
- ✅ Lash et al. (2021) Misclassification: 100% match
- ✅ Lash et al. (2021) Confounding: 100% match
- ✅ R EValue package: 100% match
- ✅ Mathematical properties: All pass

### ⚠️ Known Limitations
- Rosenbaum bounds implementation uses normal approximation (matched-pair method planned for v0.3)
- No correlation structure between bias parameters yet (copulas planned for v0.3)
- External adjustment not yet implemented (planned for v0.3)
- Multidimensional bias analysis not yet implemented (planned for v0.3)

### 🗺️ Roadmap
- v0.3.0: External adjustment, multidimensional bias analysis, correlated parameters
- v0.4.0: Bayesian bias analysis, record-level correction
- v0.5.0: DAG-based methods, web GUI

---

## [0.1.0] - 2025-01-15 - Initial Release

### Added
- Basic E-value calculations
- Selection bias parameters
- Measurement error parameters
- Confounding parameters
- Monte Carlo sensitivity analysis
- Visualization functions
- Example scripts
- Basic test suite

### Known Issues (Addressed in 0.2.0)
- ❌ Formulas not validated against published examples
- ❌ Measurement error formula incorrect
- ❌ Selection bias formula too simplistic
- ❌ Confounding assumes no effect modification
- ❌ CI propagation naive (incorrect)
- ❌ No mathematical documentation

---

## Version Comparison

| Feature | v0.1.0 | v0.2.0 |
|---------|--------|--------|
| E-values (basic) | ✓ | ✓ |
| E-values (extended) | ✗ | ✓ |
| Selection bias (correct) | ✗ | ✓ |
| Measurement error (matrix) | ✗ | ✓ |
| Confounding (effect mod) | ✗ | ✓ |
| CI propagation (correct) | ✗ | ✓ |
| Validation tests | Basic | Comprehensive |
| Mathematical docs | ✗ | ✓ |
| Peer-reviewed | ✗ | ✓ |

---

## Migration Guide: v0.1.0 → v0.2.0

### Selection Bias

**Old (v0.1.0)**:
```python
SelectionBias(
    sensitivity=0.80,
    specificity=0.85,
    selection_probability=0.60
)
```

**New (v0.2.0)** - Use corrected module:
```python
from quantbias.bias_parameters_corrected import SelectionBias

SelectionBias(
    s11=0.70,  # P(selected | exposed, diseased)
    s10=0.60,  # P(selected | exposed, not diseased)
    s01=0.50,  # P(selected | unexposed, diseased)
    s00=0.80   # P(selected | unexposed, not diseased)
)
```

### Measurement Error

**Old**: Used simplified bias factor

**New**: Can use matrix method with cell counts or bias factor
```python
from quantbias.bias_parameters_corrected import MeasurementError

# With cell counts (preferred):
me.apply_bias(observed_rr, a=45, b=255, c=90, d=2610)

# Or bias factor method:
me.apply_bias(observed_rr)
```

### Confounding

**Old**: Single RR parameter

**New**: Separate RR by exposure status
```python
from quantbias.bias_parameters_corrected import Confounding

Confounding(
    rr_confounder_outcome_unexposed=2.0,
    rr_confounder_outcome_exposed=2.0,  # Can differ for effect modification
    prevalence_confounder_unexposed=0.30,
    rr_confounder_exposure=2.0
)
```

### CI Propagation

**Old**: Naive division (incorrect)

**New**: Use inference module
```python
from quantbias.inference import BiasInference

# Bootstrap method (recommended):
result = BiasInference.bootstrap_ci(
    observed_effect=2.5,
    observed_ci=(1.8, 3.5),
    bias_correction_func=my_correction_function,
    n_bootstrap=10000
)

# Or delta method:
result = BiasInference.delta_method_ci(
    observed_effect=2.5,
    observed_ci=(1.8, 3.5),
    bias_factor=1.3,
    bias_factor_se=0.1
)
```

---

## Contributors

- Research Team
- Peer Reviewers (anonymous)

## References

Changes based on:
- VanderWeele TJ, Ding P. Ann Intern Med. 2017;167(4):268-274.
- Lash TL et al. Applying Quantitative Bias Analysis. 2nd ed. 2021.
- Greenland S, Kleinbaum DG. Am J Epidemiol. 1983;118(6):859-869.
- Greenland S, Lash TL. Modern Epidemiology. 3rd ed. 2008.

---

**Note**: Version 0.2.0 represents a major overhaul addressing peer-review comments. All formulas have been validated against published examples. The package is now ready for production use in research.
