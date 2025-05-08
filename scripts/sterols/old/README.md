# Yeast Sterol Analysis

This directory contains scripts for analyzing sterol profiles in relationship to genomic patterns in the yeast (*S. cerevisiae*, W303 strain) adaptation project.

## Overview

The sterol analysis focuses on understanding the biochemical consequences of adaptation to different stresses (temperature, low oxygen) and gene modifications (CAS, STC) on the ergosterol biosynthetic pathway. This work builds on the genomic findings that revealed a hierarchical conservation pattern in the ergosterol pathway genes.

## Scripts

The analysis is implemented in two main scripts:

1. **sterol_analysis.py** - Core analysis of sterol profiles:
   - Exploratory analysis of sterol distributions
   - Differential analysis comparing treatments
   - Pathway flux analysis
   - Statistical pattern recognition using PCA and clustering

2. **genomic_correlation.py** - Integration with genomic patterns:
   - Correlation between variant counts and sterol levels
   - Analysis of adaptation-specific sterol patterns
   - Relationship between sterol ratios and genetic variation
   - Biological interpretation connecting genomic conservation with biochemical adaptation

## Key Findings

1. **Adaptation-Specific Sterol Patterns**:
   - Temperature adaptation is associated with maintained or increased ergosterol levels
   - Low oxygen adaptation shows significantly reduced ergosterol levels (30% of control)
   - Different adaptation conditions have distinct sterol profile signatures

2. **Gene Modification Effects**:
   - Gene modifications (CAS, STC) alter sterol profiles in adaptation-specific ways
   - Both gene-modified strains show decreased ergosterol compared to their wild-type counterparts

3. **Sterol Pathway Regulation**:
   - Ergosterol-to-precursor ratios differ between adaptation types
   - Low oxygen adaptation shows higher sterol ratios than temperature adaptation

4. **Genomic-Sterol Integration**:
   - Despite complete conservation of ergosterol pathway genes, adaptation produces significant sterol profile changes
   - This suggests regulation occurs through mechanisms other than direct genetic modification of pathway enzymes
   - The satellite gene architecture may provide a flexible regulatory layer for adaptation

## Usage

To run the complete sterol analysis workflow:

```bash
# Run the core sterol profile analysis
python3 scripts/sterols/sterol_analysis.py

# Run the genomic correlation analysis
python3 scripts/sterols/genomic_correlation.py
```

## Results

Analysis results are saved to the following directories:

- `results/sterol_analysis/exploratory/` - Basic sterol distribution visualizations
- `results/sterol_analysis/differential/` - Fold change analysis and comparisons
- `results/sterol_analysis/pathway_analysis/` - Pathway flux diagrams and ratios
- `results/sterol_analysis/pattern_analysis/` - PCA plots and clustering results
- `results/sterol_analysis/genomic_correlation/` - Integrated genomic-sterol analysis

Key reports:
- `results/sterol_analysis/sterol_analysis_summary.md` - Summary of core sterol analysis
- `results/sterol_analysis/genomic_correlation/genomic_sterol_integrated_report.md` - Integrated findings connecting genomic patterns with biochemical changes

## Biological Significance

The sterol analysis provides crucial biochemical evidence supporting the genomic findings of purifying selection on the ergosterol pathway. It demonstrates that adaptation to environmental stresses occurs through regulatory mechanisms that preserve essential pathway genes while modulating their output.

The distinct sterol profiles in temperature versus low oxygen adaptation reveal different regulatory strategies, suggesting that satellite genes may play important roles in adaptation-specific regulation. This provides a model for understanding how essential pathways can adapt to environmental challenges without compromising their core functions.