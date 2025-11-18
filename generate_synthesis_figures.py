"""
Generate figures for the Synthesis document
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

# Create figures directory if it doesn't exist
os.makedirs('figures', exist_ok=True)

# Set publication-quality defaults
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.dpi'] = 300

#############################################################################
# FIGURE 1: QuantBias Workflow Diagram
#############################################################################

fig1 = plt.figure(figsize=(10, 7))
ax1 = fig1.add_subplot(111)
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 10)
ax1.axis('off')

# Define colors
color_input = '#E8F4F8'
color_method = '#FFF4E6'
color_output = '#E8F5E9'
color_arrow = '#666666'

# Helper function to create boxes
def create_box(ax, x, y, width, height, text, color, fontsize=10, fontweight='normal'):
    box = FancyBboxPatch((x, y), width, height,
                          boxstyle="round,pad=0.1",
                          edgecolor='#333333',
                          facecolor=color,
                          linewidth=1.5)
    ax.add_patch(box)
    ax.text(x + width/2, y + height/2, text,
            ha='center', va='center',
            fontsize=fontsize, fontweight=fontweight,
            wrap=True)

# Helper function to create arrows
def create_arrow(ax, x1, y1, x2, y2):
    arrow = FancyArrowPatch((x1, y1), (x2, y2),
                           arrowstyle='->',
                           mutation_scale=20,
                           linewidth=2,
                           color=color_arrow)
    ax.add_patch(arrow)

# Title
ax1.text(5, 9.5, 'QuantBias Analysis Workflow',
         ha='center', va='top', fontsize=14, fontweight='bold')

# Level 1: Input Data
create_box(ax1, 3.5, 8, 3, 0.8, 'Observed Effect Estimate\n(RR, OR, HR, RD)',
           color_input, fontsize=10, fontweight='bold')

# Arrow down
create_arrow(ax1, 5, 8, 5, 7.2)

# Level 2: Three Analysis Paths
# Path 1: E-values
create_box(ax1, 0.2, 5.8, 2.8, 1.2, 'E-value Analysis\n\n• Unmeasured confounding\n• Quick threshold assessment',
           color_method, fontsize=9)

# Path 2: Deterministic Correction
create_box(ax1, 3.6, 5.8, 2.8, 1.2, 'Deterministic Correction\n\n• Specify bias parameters\n• Single corrected estimate',
           color_method, fontsize=9)

# Path 3: Probabilistic Analysis
create_box(ax1, 7.0, 5.8, 2.8, 1.2, 'Probabilistic Analysis\n\n• Parameter distributions\n• Monte Carlo simulation',
           color_method, fontsize=9)

# Arrows from input to three paths
create_arrow(ax1, 4.3, 7.2, 1.6, 7)
create_arrow(ax1, 5, 7.2, 5, 7)
create_arrow(ax1, 5.7, 7.2, 8.4, 7)

# Level 3: Bias Sources (for paths 2 and 3)
create_box(ax1, 3.0, 4.0, 1.4, 1.0, 'Selection\nBias', '#FFE5E5', fontsize=9)
create_box(ax1, 4.6, 4.0, 1.4, 1.0, 'Measurement\nError', '#FFE5E5', fontsize=9)
create_box(ax1, 6.2, 4.0, 1.4, 1.0, 'Unmeasured\nConfounding', '#FFE5E5', fontsize=9)

# Arrows to bias sources
create_arrow(ax1, 4.5, 5.8, 3.7, 5.0)
create_arrow(ax1, 5.0, 5.8, 5.3, 5.0)
create_arrow(ax1, 5.5, 5.8, 6.9, 5.0)
create_arrow(ax1, 8.0, 5.8, 5.3, 5.0)
create_arrow(ax1, 8.5, 5.8, 6.9, 5.0)

# Level 4: Statistical Inference
create_box(ax1, 3.5, 2.5, 3, 0.9, 'Statistical Inference\n\n• Bootstrap / Delta Method / Simulation\n• Proper CI propagation',
           color_method, fontsize=9)

# Arrows to inference
create_arrow(ax1, 4.2, 4.0, 4.5, 3.4)
create_arrow(ax1, 5.3, 4.0, 5.3, 3.4)
create_arrow(ax1, 6.4, 4.0, 6.1, 3.4)

# Level 5: Outputs
create_box(ax1, 0.5, 0.5, 2.5, 1.3, 'E-value & CI\n\nMinimum confounding\nstrength to explain\naway effect',
           color_output, fontsize=9)
create_box(ax1, 3.5, 0.5, 3, 1.3, 'Bias-Adjusted Estimate & CI\n\nCorrected RR/OR/HR\nProper uncertainty quantification\nVariance Inflation Factor',
           color_output, fontsize=9)
create_box(ax1, 7.0, 0.5, 2.5, 1.3, 'Distribution of\nAdjusted Estimates\n\nFull uncertainty\npropagation',
           color_output, fontsize=9)

# Arrows to outputs
create_arrow(ax1, 1.6, 5.8, 1.75, 1.8)
create_arrow(ax1, 5.0, 2.5, 5.0, 1.8)
create_arrow(ax1, 8.4, 5.8, 8.25, 1.8)

# Add annotation box
create_box(ax1, 0.2, 10, 9.6, 0.3,
           'Input: Observed data → Methods: Multiple analysis approaches → Output: Quantified bias impact',
           '#F5F5F5', fontsize=9, fontweight='normal')

plt.tight_layout()
plt.savefig('figures/Figure1_Workflow.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/Figure1_Workflow.pdf', bbox_inches='tight')
print("Figure 1 saved: figures/Figure1_Workflow.png and .pdf")


#############################################################################
# FIGURE 2: E-value Sensitivity Analysis Contour Plot
#############################################################################

fig2, axes = plt.subplots(1, 2, figsize=(12, 5))

# Example scenario: Observed RR = 2.5
observed_rr = 2.5

# Panel A: Contour plot showing when unmeasured confounding explains away effect
ax2a = axes[0]

# Create grid of confounder associations
rr_eu = np.linspace(1.0, 5.0, 100)  # Exposure-Confounder association
rr_ud = np.linspace(1.0, 5.0, 100)  # Confounder-Outcome association

RR_EU, RR_UD = np.meshgrid(rr_eu, rr_ud)

# Calculate bias-adjusted RR for each combination
# Using simplified confounding formula: RR_adjusted ≈ RR_obs / BF
# where BF ≈ (p1*RR_UD + (1-p1)) / (p0*RR_UD + (1-p0))
# For illustration, assume moderate confounder prevalence difference
p1 = 0.5  # Prevalence among exposed
p0 = 0.2  # Prevalence among unexposed

# More precisely: need to account for both RR_EU and RR_UD
# Using relationship: p1 = p0*RR_EU / (1 - p0 + p0*RR_EU)
P1 = p0 * RR_EU / (1 - p0 + p0 * RR_EU)
BF = (P1 * RR_UD + (1 - P1)) / (p0 * RR_UD + (1 - p0))
RR_adjusted = observed_rr / BF

# Create filled contour plot
levels = [0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0]
contourf = ax2a.contourf(RR_EU, RR_UD, RR_adjusted, levels=20, cmap='RdYlGn_r', alpha=0.8)
contour = ax2a.contour(RR_EU, RR_UD, RR_adjusted, levels=levels, colors='black', linewidths=0.5, alpha=0.5)
ax2a.clabel(contour, inline=True, fontsize=8, fmt='%.1f')

# Highlight the null line (RR = 1.0)
null_contour = ax2a.contour(RR_EU, RR_UD, RR_adjusted, levels=[1.0], colors='red', linewidths=2.5)
ax2a.clabel(null_contour, inline=True, fontsize=10, fmt='RR=%.1f')

# Calculate and plot E-value
evalue = observed_rr + np.sqrt(observed_rr * (observed_rr - 1))
ax2a.plot([1, 5], [evalue, evalue], 'b--', linewidth=2, label=f'E-value = {evalue:.2f}')
ax2a.plot([evalue, evalue], [1, 5], 'b--', linewidth=2)

# Add shaded region where both associations exceed E-value
ax2a.fill_between([evalue, 5], evalue, 5, alpha=0.2, color='blue',
                   label='Both associations\nexceed E-value')

ax2a.set_xlabel('Exposure-Confounder Association (RR)', fontweight='bold')
ax2a.set_ylabel('Confounder-Outcome Association (RR)', fontweight='bold')
ax2a.set_title(f'Panel A: Bias-Adjusted RR for Observed RR = {observed_rr}', fontweight='bold')
ax2a.grid(True, alpha=0.3)
ax2a.legend(loc='upper left', framealpha=0.9)
ax2a.set_xlim(1, 5)
ax2a.set_ylim(1, 5)

# Add colorbar
cbar = plt.colorbar(contourf, ax=ax2a)
cbar.set_label('Bias-Adjusted RR', rotation=270, labelpad=20, fontweight='bold')

# Panel B: Bar chart showing E-values for different observed effects
ax2b = axes[1]

observed_effects = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
evalues = [rr + np.sqrt(rr * (rr - 1)) for rr in observed_effects]
ci_evalues = [(rr * 0.7) + np.sqrt((rr * 0.7) * (rr * 0.7 - 1)) if rr * 0.7 > 1 else 1.0
              for rr in observed_effects]

x = np.arange(len(observed_effects))
width = 0.35

bars1 = ax2b.bar(x - width/2, evalues, width, label='Point Estimate E-value',
                 color='#2E86AB', edgecolor='black', linewidth=1)
bars2 = ax2b.bar(x + width/2, ci_evalues, width, label='CI Lower Bound E-value',
                 color='#A23B72', edgecolor='black', linewidth=1)

# Add value labels on bars
for i, (bar1, bar2) in enumerate(zip(bars1, bars2)):
    height1 = bar1.get_height()
    height2 = bar2.get_height()
    ax2b.text(bar1.get_x() + bar1.get_width()/2., height1,
             f'{height1:.2f}', ha='center', va='bottom', fontsize=8, fontweight='bold')
    if height2 > 1.01:  # Only show if meaningfully above 1
        ax2b.text(bar2.get_x() + bar2.get_width()/2., height2,
                 f'{height2:.2f}', ha='center', va='bottom', fontsize=8, fontweight='bold')

ax2b.set_xlabel('Observed Risk Ratio', fontweight='bold')
ax2b.set_ylabel('E-value', fontweight='bold')
ax2b.set_title('Panel B: E-values for Different Effect Sizes', fontweight='bold')
ax2b.set_xticks(x)
ax2b.set_xticklabels([f'{rr:.1f}' for rr in observed_effects])
ax2b.legend(loc='upper left', framealpha=0.9)
ax2b.grid(True, axis='y', alpha=0.3)
ax2b.set_ylim(0, max(evalues) * 1.15)

# Add horizontal reference lines
ax2b.axhline(y=2.0, color='orange', linestyle='--', linewidth=1, alpha=0.7, label='Modest confounding')
ax2b.axhline(y=4.0, color='red', linestyle='--', linewidth=1, alpha=0.7, label='Strong confounding')
ax2b.legend(loc='upper left', framealpha=0.9)

plt.tight_layout()
plt.savefig('figures/Figure2_Evalue_Sensitivity.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/Figure2_Evalue_Sensitivity.pdf', bbox_inches='tight')
print("Figure 2 saved: figures/Figure2_Evalue_Sensitivity.png and .pdf")

print("\nBoth figures generated successfully!")
print("- Figure 1: QuantBias workflow diagram")
print("- Figure 2: E-value sensitivity analysis")
