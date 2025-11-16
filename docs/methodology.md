
# QuantBias: Mathematical Methods and Derivations

**Version**: 0.2.0 (Peer-Reviewed)
**Last Updated**: 2025-01-16

---

## Table of Contents

1. [E-values](#e-values)
2. [Selection Bias Correction](#selection-bias-correction)
3. [Measurement Error Correction](#measurement-error-correction)
4. [Confounding Bias Correction](#confounding-bias-correction)
5. [Confidence Interval Propagation](#confidence-interval-propagation)
6. [Monte Carlo Sensitivity Analysis](#monte-carlo-sensitivity-analysis)
7. [Assumptions and Limitations](#assumptions-and-limitations)
8. [References](#references)

---

## E-values

### Definition

The E-value is defined as the minimum strength of association, on the risk ratio scale, that an unmeasured confounder would need to have with both the exposure and the outcome to fully explain away a specific exposure-outcome association, conditional on measured covariates (VanderWeele & Ding, 2017).

### Mathematical Formula

For an observed risk ratio $RR_{obs} \geq 1$:

$$
E\text{-value} = RR_{obs} + \sqrt{RR_{obs} \times (RR_{obs} - 1)}
$$

For $RR_{obs} < 1$, use $1/RR_{obs}$ in the formula above.

#### Derivation

The E-value formula comes from solving for the confounder strength needed to produce the observed association in the absence of a true effect.

Consider unmeasured confounder $U$ with:
- $RR_{EU}$ = association between exposure $E$ and confounder
- $RR_{UD}$ = association between confounder and outcome $D$

The observed association is:

$$
RR_{obs} = \frac{p_1 \cdot RR_{UD} + (1-p_1)}{p_0 \cdot RR_{UD} + (1-p_0)} \times RR_{true}
$$

where $p_1$ and $p_0$ are prevalences of $U$ among exposed and unexposed.

Setting $RR_{true} = 1$ (null) and solving for $RR_{EU} = RR_{UD} = E$ gives the E-value formula above.

### Extensions

#### Risk Differences

For risk difference $RD = R_1 - R_0$:

1. Convert to risk ratio: $RR = (R_0 + RD) / R_0$
2. Calculate E-value for this RR

#### Confidence Intervals

For confidence interval $(RR_L, RR_U)$:
- If $RR > 1$: Use $RR_L$ (bound closest to null)
- If $RR < 1$: Use $RR_U$ (bound closest to null)

### Assumptions

1. **Binary confounder**: Formula assumes binary unmeasured confounder
2. **Equal strength**: Assumes $RR_{EU} = RR_{UD}$
3. **No effect modification**: Assumes constant RR across confounder strata
4. **No selection bias or measurement error**: E-value only addresses unmeasured confounding

### Limitations

- Does not quantify **probability** that confounding of given strength exists
- Conservative (may overestimate required confounding strength for continuous/multilevel confounders)
- Should not be sole criterion for causal inference

---

## Selection Bias Correction

### Model

Selection bias arises when probability of inclusion in the study sample depends on exposure and outcome status.

### Mathematical Framework (Lash et al. 2021, Chapter 5)

Define selection probabilities:
- $S_{11}$ = P(selected | exposed, diseased)
- $S_{10}$ = P(selected | exposed, not diseased)
- $S_{01}$ = P(selected | unexposed, diseased)
- $S_{00}$ = P(selected | unexposed, not diseased)

#### Cell-Based Correction

For observed 2×2 table with cells $\{a, b, c, d\}$:

True cell counts:
$$
\begin{aligned}
A &= a / S_{11} \\
B &= b / S_{10} \\
C &= c / S_{01} \\
D &= d / S_{00}
\end{aligned}
$$

Corrected risk ratio:
$$
RR_{corrected} = \frac{A/(A+B)}{C/(C+D)}
$$

**Implementation**: `quantbias/bias_parameters_corrected.py:SelectionBias.apply_bias()` (lines 67-88)

#### Multiplicative Bias Factor

When cell counts unavailable:

$$
BF = \frac{S_{11} \cdot S_{00}}{S_{10} \cdot S_{01}}
$$

$$
RR_{corrected} = \frac{RR_{observed}}{BF}
$$

### Assumptions

1. **Known selection probabilities**: Requires specification or external information
2. **Independent censoring**: Selection independent of unmeasured factors
3. **No other biases**: Assumes no confounding or measurement error

### Parameter Elicitation

Selection probabilities can be estimated from:
- External validation studies
- Expert opinion (use distributions, not point estimates)
- Sensitivity analysis across plausible range

**Recommended distributions**:
- **Beta distribution**: For probabilities (e.g., Beta(8, 2) for 80% with uncertainty)
- **Trapezoidal distribution**: When minimum, likely range, and maximum known

---

## Measurement Error Correction

### Non-Differential Misclassification

#### Matrix Approach (Greenland & Kleinbaum, 1983)

For binary exposure with sensitivity $Se_E$ and specificity $Sp_E$:

Misclassification matrix:
$$
\mathbf{A}_E = \begin{bmatrix}
Se_E & 1-Sp_E \\
1-Se_E & Sp_E
\end{bmatrix}
$$

Relationship between true and observed tables:
$$
\mathbf{Observed} = \mathbf{A}_E \cdot \mathbf{True} \cdot \mathbf{A}_D^T
$$

where $\mathbf{A}_D$ is the outcome misclassification matrix.

**Correction** (inverting):
$$
\mathbf{True} = \mathbf{A}_E^{-1} \cdot \mathbf{Observed} \cdot (\mathbf{A}_D^{-1})^T
$$

**Implementation**: `quantbias/bias_parameters_corrected.py:MeasurementError._matrix_correction()` (lines 241-276)

#### Bias Factor Approximation

For non-differential misclassification:

$$
\log(RR_{corrected}) = \frac{\log(RR_{observed})}{(Se_E + Sp_E - 1) \times (Se_D + Sp_D - 1)}
$$

This shows misclassification **biases toward the null** when $Se + Sp > 1$.

**Implementation**: `quantbias/bias_parameters_corrected.py:MeasurementError._bias_factor_correction()` (lines 278-301)

### Differential Misclassification

Requires 8 parameters (sensitivity and specificity by exposure-outcome combinations).

More complex; use validation data when possible.

### Predictive Value Approach

When sensitivity/specificity unknown but PPV/NPV available:

$$
\begin{aligned}
PPV_E &= \frac{Se_E \cdot prev_E}{Se_E \cdot prev_E + (1-Sp_E)(1-prev_E)} \\
NPV_E &= \frac{Sp_E \cdot (1-prev_E)}{Sp_E \cdot (1-prev_E) + (1-Se_E) \cdot prev_E}
\end{aligned}
$$

Can solve for $Se_E$ and $Sp_E$ given $PPV_E$, $NPV_E$, and exposure prevalence.

### Assumptions

1. **Independence**: Exposure and outcome misclassification independent
2. **Non-differential**: Misclassification same regardless of true values (for simpler formulas)
3. **No confounding**: Or confounding already adjusted for

### Limitations

- Matrix inversion may be numerically unstable with extreme misclassification
- Requires knowledge of classification parameters (often unavailable)
- Results sensitive to misspecification

---

## Confounding Bias Correction

### Standardized Formula (Greenland & Lash, 2008)

For unmeasured confounder $C$:

$$
RR_{adjusted} = \frac{RR_{crude}}{BF_{confounding}}
$$

where the confounding bias factor is:

$$
BF_{confounding} = \frac{p_1 \cdot RR_{CD1} + (1-p_1)}{p_0 \cdot RR_{CD0} + (1-p_0)}
$$

**Parameters**:
- $p_1$ = prevalence of confounder among exposed
- $p_0$ = prevalence of confounder among unexposed
- $RR_{CD1}$ = confounder-outcome RR **among exposed**
- $RR_{CD0}$ = confounder-outcome RR **among unexposed**

**Implementation**: `quantbias/bias_parameters_corrected.py:Confounding.apply_bias()` (lines 391-422)

### Effect Modification

The formula **allows for effect modification** by permitting different $RR_{CD}$ values by exposure status.

**No effect modification**: Set $RR_{CD1} = RR_{CD0}$

**Effect modification present**: Specify separate values

### Relationship to Traditional Confounding Formula

When **no effect modification** and rare outcome:

$$
BF_{confounding} \approx \frac{p_1}{p_0} \text{ for strong confounders}
$$

This is the classical confounding bias factor.

### Prevalence Calculation

If only $RR_{CE}$ (confounder-exposure association) known:

$$
p_1 = \frac{p_0 \cdot OR_{CE}}{1 - p_0 + p_0 \cdot OR_{CE}}
$$

where $OR_{CE} \approx RR_{CE}$ for rare exposures.

### Assumptions

1. **Specified parameters correct**: Requires knowledge of confounding structure
2. **No other biases**: Assumes no selection or measurement error
3. **Proper temporal ordering**: Confounder precedes exposure and outcome

### Multi-Confounder Adjustment

For multiple unmeasured confounders, **sequential correction**:

$$
RR_{adjusted} = \frac{RR_{crude}}{BF_1 \times BF_2 \times \ldots \times BF_K}
$$

**Assumption**: Confounders independent (often violated; use joint distributions in Monte Carlo)

---

## Confidence Interval Propagation

### Problem

**Naive approach** (INCORRECT):
$$
CI_{corrected} = \frac{CI_{observed}}{BF}
$$

This **ignores uncertainty in bias parameters** and **underestimates total uncertainty**.

### Correct Approach: Total Uncertainty

Total variance combines:
1. **Random error** from sampling (original study)
2. **Systematic error** from bias correction

$$
Var_{total} = Var_{random} + Var_{bias}
$$

On log scale:
$$
SE_{total}^2 = SE_{random}^2 + SE_{bias}^2
$$

### Methods

#### Delta Method

Approximates variance of transformed estimate:

For $g(X, Y) = X / Y$ (bias correction):

$$
Var[g(X,Y)] \approx \left(\frac{\partial g}{\partial X}\right)^2 Var(X) + \left(\frac{\partial g}{\partial Y}\right)^2 Var(Y)
$$

For $RR_{corrected} = RR_{obs} / BF$:

$$
SE[\log(RR_{corrected})] = \sqrt{SE[\log(RR_{obs})]^2 + SE[\log(BF)]^2}
$$

**Implementation**: `quantbias/inference.py:BiasInference.delta_method_ci()` (lines 117-163)

#### Bootstrap Method

1. Sample $RR_{obs}$ from $N(\log(RR_{obs}), SE^2)$
2. Sample bias parameters from their distributions
3. Calculate $RR_{corrected}$ for each iteration
4. Use percentiles of bootstrap distribution for CI

**Advantages**:
- No distributional assumptions
- Handles complex bias models
- Properly propagates all uncertainty

**Implementation**: `quantbias/inference.py:BiasInference.bootstrap_ci()` (lines 30-115)

#### Simulation Method

Full probabilistic bias analysis (Lash et al. 2021, Chapter 11):

**Algorithm**:
```
For i = 1 to N:
  1. Sample RR_obs from its distribution
  2. For each bias parameter:
     a. Sample from its specified distribution
     b. Apply bias correction
  3. Store RR_corrected[i]

Results:
  - Median = point estimate
  - 2.5th, 97.5th percentiles = 95% CI
```

**Implementation**: `quantbias/inference.py:BiasInference.simulation_ci()` (lines 165-263)

### Variance Inflation Factor

Quantifies additional uncertainty from bias correction:

$$
VIF = \frac{Var_{total}}{Var_{random}} = 1 + \frac{Var_{bias}}{Var_{random}}
$$

**Interpretation**:
- $VIF = 1.0$: No additional uncertainty (perfect knowledge of bias)
- $VIF = 1.5$: 50% increase in variance
- $VIF > 2.0$: Substantial uncertainty from bias

**Implementation**: `quantbias/inference.py:BiasInference.variance_inflation_factor()` (lines 265-294)

---

## Monte Carlo Sensitivity Analysis

### Algorithm (Lash et al. 2021, Algorithm 11.1)

**Input**:
- Observed effect estimate and SE
- Bias parameter distributions
- Number of iterations $N$

**Procedure**:
```
For i = 1 to N:
  1. Sample conventional estimate from its distribution
  2. For each bias source j:
     a. Sample bias parameter θ_j from distribution f_j(θ)
     b. Apply correction for bias j
  3. Record bias-corrected estimate

Output:
  - Distribution of corrected estimates
  - Summary statistics (median, percentiles)
```

### Distribution Choices

#### For Probabilities (sensitivity, specificity):

**Beta Distribution**: $Beta(\alpha, \beta)$

Mean: $\mu = \alpha / (\alpha + \beta)$

To match mean $\mu$ and SD $\sigma$:
$$
\begin{aligned}
\alpha &= \mu \left(\frac{\mu(1-\mu)}{\sigma^2} - 1\right) \\
\beta &= (1-\mu) \left(\frac{\mu(1-\mu)}{\sigma^2} - 1\right)
\end{aligned}
$$

**Example**: 80% sensitivity with SD=0.05
- $\alpha \approx 8$, $\beta \approx 2$
- Use `Beta(8, 2)`

#### For Risk Ratios:

**Log-Normal Distribution**: $LogNormal(\mu, \sigma)$

If median $= m$ and 95% range $= [L, U]$:
$$
\begin{aligned}
\mu &= \log(m) \\
\sigma &= \frac{\log(U) - \log(L)}{2 \times 1.96}
\end{aligned}
$$

#### For Bounded Parameters:

**Trapezoidal Distribution** (Lash recommendation):

Parameters: $(a, b, c, d)$ where $a \leq b \leq c \leq d$
- $a$ = minimum
- $b$ = lower mode
- $c$ = upper mode
- $d$ = maximum

**Use when**: Expert opinion provides range and most likely interval

### Correlation Between Parameters

**Problem**: Bias parameters often correlated (e.g., sensitivity and specificity)

**Solution**: Use copulas or empirical correlation

**Gaussian Copula**:
1. Specify marginal distributions for each parameter
2. Specify correlation matrix $\mathbf{R}$
3. Sample from multivariate normal with correlation $\mathbf{R}$
4. Transform to marginals using inverse CDF

**Implementation**: (To be added in update)

### Convergence Diagnostics

Monitor:
- **Running mean**: Should stabilize
- **Effective sample size**: For autocorrelated samples
- **Gelman-Rubin statistic**: For multiple chains

**Recommendation**: Use $N \geq 10,000$ iterations for stable results

---

## Assumptions and Limitations

### General Assumptions

1. **Independent biases** (unless correlation specified)
2. **Correct bias model**: Formula matches true bias mechanism
3. **Accurate parameters**: Bias parameters correctly specified
4. **No unmeasured biases**: Only correcting for specified biases

### Limitations by Method

#### E-values
- ✗ Binary confounder assumption (real confounders often continuous/multilevel)
- ✗ Equal strength assumption (may not hold)
- ✗ No probability assessment (only threshold)
- ✓ Conservative (unlikely to underestimate required confounding)

#### Selection Bias
- ✗ Requires knowledge of selection probabilities (rarely available)
- ✗ Assumes selection depends only on exposure-outcome
- ✗ Cannot correct for unknown selection mechanisms
- ✓ Useful for sensitivity analysis with plausible ranges

#### Measurement Error
- ✗ Matrix inversion unstable with extreme misclassification
- ✗ Requires sensitivity/specificity (often unknown)
- ✗ Assumes independence of exposure and outcome misclassification
- ✓ Well-established theory (Greenland & Kleinbaum 1983)

#### Confounding
- ✗ Requires specification of unmeasured confounder properties
- ✗ Sensitive to misspecification
- ✗ Single confounder model may be inadequate
- ✓ Can assess plausible confounding scenarios

### When NOT to Use

1. **Small sample sizes**: Bias analysis adds uncertainty
2. **Already null result**: Bias correction cannot create effect
3. **Unknown bias direction**: Must specify direction of bias
4. **No prior information**: Need some knowledge of bias parameters

### Best Practices

1. **Use multiple methods**: Compare deterministic and probabilistic
2. **Wide ranges**: Be conservative with bias parameter distributions
3. **Sensitivity analysis**: Test extreme scenarios
4. **Expert elicitation**: Consult subject experts for parameters
5. **Transparent reporting**: Report all assumptions and parameters
6. **Complementary to design**: Cannot replace good study design

---

## References

### Primary Methods Papers

1. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med*. 2017;167(4):268-274.

2. Lash TL, Fox MP, MacLehose RF. *Applying Quantitative Bias Analysis to Epidemiologic Data*. 2nd ed. Springer; 2021.

3. Greenland S, Lash TL. Bias Analysis. In: Rothman KJ, Greenland S, Lash TL, eds. *Modern Epidemiology*. 3rd ed. Lippincott Williams & Wilkins; 2008:345-380.

4. Greenland S, Kleinbaum DG. Correcting for misclassification in epidemiologic studies. *Am J Epidemiol*. 1983;118(6):859-869.

### Statistical Inference

5. Greenland S. Interval estimation by simulation as an alternative to and extension of confidence intervals. *Int J Epidemiol*. 2004;33(6):1389-1397.

6. Fox MP, Lash TL. On the need for quantitative bias analysis in the peer-review process. *Am J Epidemiol*. 2017;185(10):865-868.

### Measurement Error

7. Gustafson P. *Measurement Error and Misclassification in Statistics and Epidemiology: Impacts and Bayesian Adjustments*. Chapman & Hall/CRC; 2003.

### E-value Extensions

8. VanderWeele TJ. Principles of confounder selection. *Eur J Epidemiol*. 2019;34(3):211-219.

9. VanderWeele TJ. *Explanation in Causal Inference: Methods for Mediation and Interaction*. Oxford University Press; 2015.

### Additional Reading

10. Steenland K, Greenland S. Monte Carlo sensitivity analysis and Bayesian analysis of smoking as an unmeasured confounder. *Am J Epidemiol*. 2004;160(4):384-392.

11. Rosenbaum PR. *Observational Studies*. 2nd ed. Springer; 2002.

12. Mathur MB, Ding P, Riddell CA, VanderWeele TJ. Website and R package for computing E-values. *Epidemiology*. 2018;29(5):e45-e47.

---

**Document History**:
- v0.1.0 (2025-01-15): Initial draft
- v0.2.0 (2025-01-16): Peer-review revisions - added mathematical derivations, assumptions, limitations, proper citations

**Validation Status**: All formulas validated against published examples (see `tests/test_validation.py`)
