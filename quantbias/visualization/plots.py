"""
Visualization functions for bias analysis results.

Provides publication-ready plots for:
- E-value curves
- Sensitivity analysis
- Monte Carlo distributions
- Tipping point analysis
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Tuple, Union
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from quantbias.evalues import EValue
from quantbias.monte_carlo import MonteCarloResult


# Set style
sns.set_style("whitegrid")
sns.set_palette("husl")


def plot_evalue_curve(
    rr_range: Optional[Tuple[float, float]] = None,
    observed_rr: Optional[float] = None,
    figsize: Tuple[float, float] = (10, 6),
    ax: Optional[Axes] = None
) -> Figure:
    """
    Plot E-value as a function of risk ratio.

    Parameters
    ----------
    rr_range : tuple of float, optional
        Range of RR values to plot (default: 0.5 to 5.0)
    observed_rr : float, optional
        Observed RR to highlight
    figsize : tuple
        Figure size
    ax : matplotlib Axes, optional
        Axes to plot on

    Returns
    -------
    Figure
        Matplotlib figure
    """
    if rr_range is None:
        rr_range = (0.5, 5.0)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Generate RR values
    rr_values = np.linspace(rr_range[0], rr_range[1], 1000)

    # Calculate E-values
    evalues = [EValue.calculate_evalue_rr(rr) for rr in rr_values]

    # Plot
    ax.plot(rr_values, evalues, linewidth=2, color='steelblue', label='E-value')

    # Add reference lines
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(x=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    # Highlight observed RR if provided
    if observed_rr is not None:
        evalue_obs = EValue.calculate_evalue_rr(observed_rr)
        ax.plot(observed_rr, evalue_obs, 'ro', markersize=10,
                label=f'Observed RR = {observed_rr:.2f}')
        ax.axvline(x=observed_rr, color='red', linestyle=':', alpha=0.5)
        ax.axhline(y=evalue_obs, color='red', linestyle=':', alpha=0.5)

        # Add annotation
        ax.annotate(
            f'E-value = {evalue_obs:.2f}',
            xy=(observed_rr, evalue_obs),
            xytext=(10, 10),
            textcoords='offset points',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
        )

    ax.set_xlabel('Risk Ratio', fontsize=12)
    ax.set_ylabel('E-value', fontsize=12)
    ax.set_title('E-value Curve', fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def plot_sensitivity_contour(
    sensitivity_results: pd.DataFrame,
    figsize: Tuple[float, float] = (10, 8),
    ax: Optional[Axes] = None
) -> Figure:
    """
    Plot sensitivity analysis contour.

    Parameters
    ----------
    sensitivity_results : pd.DataFrame
        Results from sensitivity_grid
    figsize : tuple
        Figure size
    ax : matplotlib Axes, optional
        Axes to plot on

    Returns
    -------
    Figure
        Matplotlib figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Pivot data for contour plot
    pivot_data = sensitivity_results.pivot_table(
        values='Adjusted_RR',
        index='RR_confounder_outcome',
        columns='Prevalence_exposed',
        aggfunc='mean'
    )

    # Create contour plot
    contour = ax.contourf(
        pivot_data.columns,
        pivot_data.index,
        pivot_data.values,
        levels=20,
        cmap='RdYlBu_r'
    )

    # Add contour lines
    contour_lines = ax.contour(
        pivot_data.columns,
        pivot_data.index,
        pivot_data.values,
        levels=[0.8, 1.0, 1.2, 1.5, 2.0, 3.0],
        colors='black',
        linewidths=1,
        alpha=0.4
    )
    ax.clabel(contour_lines, inline=True, fontsize=8)

    # Add null line
    null_contour = ax.contour(
        pivot_data.columns,
        pivot_data.index,
        pivot_data.values,
        levels=[1.0],
        colors='red',
        linewidths=2,
        linestyles='--'
    )
    ax.clabel(null_contour, inline=True, fontsize=10, fmt='RR = %.1f')

    # Colorbar
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Adjusted Risk Ratio', fontsize=11)

    ax.set_xlabel('Prevalence of Confounder (Exposed)', fontsize=12)
    ax.set_ylabel('RR Confounder-Outcome', fontsize=12)
    ax.set_title('Sensitivity Analysis: Unmeasured Confounding', fontsize=14, fontweight='bold')

    plt.tight_layout()
    return fig


def plot_monte_carlo_distribution(
    mc_result: MonteCarloResult,
    figsize: Tuple[float, float] = (12, 6),
    bins: int = 50
) -> Figure:
    """
    Plot Monte Carlo simulation results.

    Parameters
    ----------
    mc_result : MonteCarloResult
        Results from Monte Carlo analysis
    figsize : tuple
        Figure size
    bins : int
        Number of histogram bins

    Returns
    -------
    Figure
        Matplotlib figure
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Histogram
    ax1.hist(mc_result.adjusted_rr_samples, bins=bins, alpha=0.7,
             color='steelblue', edgecolor='black')

    # Add vertical lines for summary statistics
    ax1.axvline(mc_result.observed_rr, color='red', linestyle='--',
                linewidth=2, label=f'Observed RR = {mc_result.observed_rr:.2f}')
    ax1.axvline(mc_result.median, color='green', linestyle='-',
                linewidth=2, label=f'Median = {mc_result.median:.2f}')
    ax1.axvline(mc_result.ci_lower, color='orange', linestyle=':',
                linewidth=1.5, label=f'95% CI: [{mc_result.ci_lower:.2f}, {mc_result.ci_upper:.2f}]')
    ax1.axvline(mc_result.ci_upper, color='orange', linestyle=':',
                linewidth=1.5)

    # Null line
    ax1.axvline(1.0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    ax1.set_xlabel('Bias-Adjusted Risk Ratio', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Distribution of Bias-Adjusted Estimates', fontsize=13, fontweight='bold')
    ax1.legend(loc='best', fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Cumulative distribution
    sorted_samples = np.sort(mc_result.adjusted_rr_samples)
    cumulative = np.arange(1, len(sorted_samples) + 1) / len(sorted_samples)

    ax2.plot(sorted_samples, cumulative, linewidth=2, color='steelblue')

    # Add reference lines
    ax2.axvline(1.0, color='gray', linestyle='--', linewidth=1,
                alpha=0.5, label='Null value')
    ax2.axvline(mc_result.median, color='green', linestyle='-',
                linewidth=2, label=f'Median = {mc_result.median:.2f}')
    ax2.axhline(0.5, color='green', linestyle=':', alpha=0.5)

    # Shade 95% interval
    ax2.axvspan(mc_result.ci_lower, mc_result.ci_upper,
                alpha=0.2, color='orange', label='95% Interval')

    ax2.set_xlabel('Bias-Adjusted Risk Ratio', fontsize=12)
    ax2.set_ylabel('Cumulative Probability', fontsize=12)
    ax2.set_title('Cumulative Distribution Function', fontsize=13, fontweight='bold')
    ax2.legend(loc='best', fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def plot_bias_decomposition(
    bias_decomp: pd.DataFrame,
    figsize: Tuple[float, float] = (10, 6)
) -> Figure:
    """
    Plot bias decomposition showing effect of each bias source.

    Parameters
    ----------
    bias_decomp : pd.DataFrame
        Bias decomposition DataFrame
    figsize : tuple
        Figure size

    Returns
    -------
    Figure
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Create bar plot
    steps = bias_decomp['Step'].values
    rrs = bias_decomp['RR'].values
    colors = ['red' if i == 0 else 'steelblue' for i in range(len(steps))]

    bars = ax.barh(steps, rrs, color=colors, alpha=0.7, edgecolor='black')

    # Add null line
    ax.axvline(1.0, color='gray', linestyle='--', linewidth=2,
               alpha=0.5, label='Null value')

    # Add value labels
    for i, (step, rr) in enumerate(zip(steps, rrs)):
        ax.text(rr + 0.05, i, f'{rr:.3f}', va='center', fontsize=10)

    ax.set_xlabel('Risk Ratio', fontsize=12)
    ax.set_ylabel('Analysis Step', fontsize=12)
    ax.set_title('Bias Decomposition: Sequential Correction', fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    return fig


def plot_tipping_point(
    tipping_results: pd.DataFrame,
    figsize: Tuple[float, float] = (10, 8)
) -> Figure:
    """
    Plot tipping point analysis results.

    Parameters
    ----------
    tipping_results : pd.DataFrame
        Results from tipping point analysis
    figsize : tuple
        Figure size

    Returns
    -------
    Figure
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Separate tipped and not-tipped
    tipped = tipping_results[tipping_results['Tipped']]
    not_tipped = tipping_results[~tipping_results['Tipped']]

    # Scatter plot
    ax.scatter(
        not_tipped['Prevalence_confounder'],
        not_tipped['RR_confounder'],
        c='steelblue',
        alpha=0.6,
        s=50,
        label='Robust (not tipped)',
        edgecolors='black',
        linewidth=0.5
    )

    ax.scatter(
        tipped['Prevalence_confounder'],
        tipped['RR_confounder'],
        c='red',
        alpha=0.6,
        s=50,
        label='Tipped to null',
        edgecolors='black',
        linewidth=0.5
    )

    ax.set_xlabel('Prevalence of Confounder', fontsize=12)
    ax.set_ylabel('RR of Confounder (Exposure & Outcome)', fontsize=12)
    ax.set_title('Tipping Point Analysis', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)

    # Add summary text
    pct_tipped = 100 * len(tipped) / len(tipping_results)
    ax.text(
        0.05, 0.95,
        f'Tipped: {len(tipped)}/{len(tipping_results)} ({pct_tipped:.1f}%)',
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    )

    plt.tight_layout()
    return fig


def plot_rosenbaum_bounds(
    rosenbaum_results: pd.DataFrame,
    figsize: Tuple[float, float] = (10, 6)
) -> Figure:
    """
    Plot Rosenbaum sensitivity bounds.

    Parameters
    ----------
    rosenbaum_results : pd.DataFrame
        Results from Rosenbaum bounds analysis
    figsize : tuple
        Figure size

    Returns
    -------
    Figure
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    gamma = rosenbaum_results['Gamma'].values
    p_lower = rosenbaum_results['P_value_lower'].values
    p_upper = rosenbaum_results['P_value_upper'].values

    # Plot bounds
    ax.plot(gamma, p_upper, linewidth=2, color='red',
            label='Upper bound (worst case)', marker='o', markersize=4)
    ax.plot(gamma, p_lower, linewidth=2, color='blue',
            label='Lower bound (best case)', marker='o', markersize=4)

    # Fill between
    ax.fill_between(gamma, p_lower, p_upper, alpha=0.3, color='gray')

    # Add significance line
    ax.axhline(0.05, color='green', linestyle='--', linewidth=2,
               alpha=0.7, label='α = 0.05')

    # Find tipping point (where upper bound crosses 0.05)
    tipping_idx = np.where(p_upper >= 0.05)[0]
    if len(tipping_idx) > 0:
        tipping_gamma = gamma[tipping_idx[0]]
        ax.axvline(tipping_gamma, color='orange', linestyle=':',
                   linewidth=2, alpha=0.7)
        ax.text(
            tipping_gamma, 0.5,
            f'Tipping point\nΓ = {tipping_gamma:.2f}',
            fontsize=10,
            ha='center',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7)
        )

    ax.set_xlabel('Sensitivity Parameter (Γ)', fontsize=12)
    ax.set_ylabel('P-value', fontsize=12)
    ax.set_title('Rosenbaum Sensitivity Analysis', fontsize=14, fontweight='bold')
    ax.set_yscale('log')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    return fig


def create_comprehensive_plot(
    observed_rr: float,
    mc_result: MonteCarloResult,
    sensitivity_grid: pd.DataFrame,
    figsize: Tuple[float, float] = (16, 12)
) -> Figure:
    """
    Create comprehensive visualization with multiple plots.

    Parameters
    ----------
    observed_rr : float
        Observed risk ratio
    mc_result : MonteCarloResult
        Monte Carlo results
    sensitivity_grid : pd.DataFrame
        Sensitivity analysis grid
    figsize : tuple
        Figure size

    Returns
    -------
    Figure
        Matplotlib figure with subplots
    """
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # E-value curve
    ax1 = fig.add_subplot(gs[0, 0])
    plot_evalue_curve(observed_rr=observed_rr, ax=ax1)

    # Monte Carlo distribution
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.hist(mc_result.adjusted_rr_samples, bins=50, alpha=0.7,
             color='steelblue', edgecolor='black')
    ax2.axvline(mc_result.median, color='green', linestyle='-',
                linewidth=2, label=f'Median = {mc_result.median:.2f}')
    ax2.axvline(observed_rr, color='red', linestyle='--',
                linewidth=2, label=f'Observed = {observed_rr:.2f}')
    ax2.set_xlabel('Bias-Adjusted RR')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Monte Carlo Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Sensitivity contour
    ax3 = fig.add_subplot(gs[1, :])
    plot_sensitivity_contour(sensitivity_grid, ax=ax3)

    fig.suptitle('Comprehensive Bias Analysis', fontsize=16, fontweight='bold', y=0.995)

    return fig
