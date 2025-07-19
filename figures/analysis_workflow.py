#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, ConnectionPatch
import numpy as np

# Set up the figure
fig, ax = plt.subplots(1, 1, figsize=(16, 12))
ax.set_xlim(0, 16)
ax.set_ylim(0, 12)
ax.axis('off')

# Define colors
colors = {
    'basic': '#4CAF50',
    'advanced': '#2196F3',
    'multi_chamber': '#FF9800',
    'data': '#9C27B0',
    'results': '#F44336'
}

# Data Input
data_box = FancyBboxPatch((1, 9), 3, 2, 
                          boxstyle="round,pad=0.1", 
                          facecolor=colors['data'], 
                          edgecolor='black', 
                          linewidth=2)
ax.add_patch(data_box)
ax.text(2.5, 10, 'Raw Data\n287,269 cells\n33,694 genes', 
        ha='center', va='center', fontsize=10, fontweight='bold', color='white')

# Basic Pipeline
basic_box = FancyBboxPatch((6, 9), 3, 2, 
                           boxstyle="round,pad=0.1", 
                           facecolor=colors['basic'], 
                           edgecolor='black', 
                           linewidth=2)
ax.add_patch(basic_box)
ax.text(7.5, 10, 'Basic Pipeline\n(01-05)\nCell Types & Basic\nCommunication', 
        ha='center', va='center', fontsize=9, fontweight='bold', color='white')

# Advanced Communication
advanced_box = FancyBboxPatch((11, 9), 3, 2, 
                              boxstyle="round,pad=0.1", 
                              facecolor=colors['advanced'], 
                              edgecolor='black', 
                              linewidth=2)
ax.add_patch(advanced_box)
ax.text(12.5, 10, 'Advanced\nCommunication\n(06)\nTemporal Patterns &\nCommunication Hubs', 
        ha='center', va='center', fontsize=9, fontweight='bold', color='white')

# Multi-Chamber Atlas
multi_box = FancyBboxPatch((6, 5.5), 3, 2, 
                           boxstyle="round,pad=0.1", 
                           facecolor=colors['multi_chamber'], 
                           edgecolor='black', 
                           linewidth=2)
ax.add_patch(multi_box)
ax.text(7.5, 6.5, 'Multi-Chamber\nAtlas\nChamber-Specific\nBiology & Cross-\nChamber Networks', 
        ha='center', va='center', fontsize=9, fontweight='bold', color='white')

# Scaled Data
scaled_box = FancyBboxPatch((1, 5.5), 3, 2, 
                            boxstyle="round,pad=0.1", 
                            facecolor=colors['data'], 
                            edgecolor='black', 
                            linewidth=2)
ax.add_patch(scaled_box)
ax.text(2.5, 6.5, 'Scaled Data\n50,000 cells\n5,000 genes\n(Optimized for\nM1 MacBook Pro)', 
        ha='center', va='center', fontsize=9, fontweight='bold', color='white')

# Results
results_box = FancyBboxPatch((11, 5.5), 3, 2, 
                             boxstyle="round,pad=0.1", 
                             facecolor=colors['results'], 
                             edgecolor='black', 
                             linewidth=2)
ax.add_patch(results_box)
ax.text(12.5, 6.5, 'Integrated\nResults\nChamber-Specific\nMarkers &\nTherapeutic Targets', 
        ha='center', va='center', fontsize=9, fontweight='bold', color='white')

# Key Findings
findings_box = FancyBboxPatch((6, 2), 3, 2, 
                              boxstyle="round,pad=0.1", 
                              facecolor='lightgray', 
                              edgecolor='black', 
                              linewidth=2)
ax.add_patch(findings_box)
ax.text(7.5, 3, 'Key Findings\n• Chamber-Specific Markers\n• Cross-Chamber Correlations\n• Therapeutic Targets\n• Clinical Implications', 
        ha='center', va='center', fontsize=8, fontweight='bold')

# Arrows
# Data to Basic
arrow1 = ConnectionPatch((4, 10), (6, 10), "data", "data",
                        arrowstyle="->", shrinkA=5, shrinkB=5,
                        mutation_scale=20, fc="black", linewidth=2)
ax.add_patch(arrow1)

# Basic to Advanced
arrow2 = ConnectionPatch((9, 10), (11, 10), "data", "data",
                        arrowstyle="->", shrinkA=5, shrinkB=5,
                        mutation_scale=20, fc="black", linewidth=2)
ax.add_patch(arrow2)

# Data to Scaled
arrow3 = ConnectionPatch((2.5, 9), (2.5, 7.5), "data", "data",
                        arrowstyle="->", shrinkA=5, shrinkB=5,
                        mutation_scale=20, fc="black", linewidth=2)
ax.add_patch(arrow3)

# Scaled to Multi-Chamber
arrow4 = ConnectionPatch((4, 6.5), (6, 6.5), "data", "data",
                        arrowstyle="->", shrinkA=5, shrinkB=5,
                        mutation_scale=20, fc="black", linewidth=2)
ax.add_patch(arrow4)

# Multi-Chamber to Results
arrow5 = ConnectionPatch((9, 6.5), (11, 6.5), "data", "data",
                        arrowstyle="->", shrinkA=5, shrinkB=5,
                        mutation_scale=20, fc="black", linewidth=2)
ax.add_patch(arrow5)

# Results to Findings
arrow6 = ConnectionPatch((12.5, 5.5), (9, 4), "data", "data",
                        arrowstyle="->", shrinkA=5, shrinkB=5,
                        mutation_scale=20, fc="black", linewidth=2)
ax.add_patch(arrow6)

# Basic to Findings
arrow7 = ConnectionPatch((7.5, 9), (7.5, 4), "data", "data",
                        arrowstyle="->", shrinkA=5, shrinkB=5,
                        mutation_scale=20, fc="black", linewidth=2)
ax.add_patch(arrow7)

# Labels
ax.text(5, 10.5, 'Full Dataset\nAnalysis', ha='center', va='center', fontsize=10, fontweight='bold')
ax.text(5, 7, 'Scaled Dataset\nAnalysis', ha='center', va='center', fontsize=10, fontweight='bold')

# Title
ax.text(8, 11.5, 'GRAIL-heart: Comprehensive Analysis Workflow', 
        ha='center', va='center', fontsize=16, fontweight='bold')

# Subtitle
ax.text(8, 11, 'Integration of Basic Pipeline, Advanced Communication, and Multi-Chamber Atlas', 
        ha='center', va='center', fontsize=12, style='italic')

plt.tight_layout()
plt.savefig('figures/analysis_workflow.png', dpi=300, bbox_inches='tight')
plt.close()

print("Workflow figure generated: figures/analysis_workflow.png")
