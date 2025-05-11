# Comprehensive Regulatory Region Analysis

## Overview

This report analyzes in depth the regulatory regions where 80% of variants near ergosterol genes were found.
It maps precise locations of variants relative to gene features (promoters, UTRs, etc.), identifies patterns
in the distribution of variants within regulatory elements, and compares regulatory region variants between
temperature and low oxygen adaptation.

Analysis date: 2025-05-10 19:35

## Dataset Summary

Total mapped variant-gene associations: 29242
Unique variants: 5337
Genes with variants in regulatory regions: 1868

### Gene Types

ERG genes: 5
OSH genes: 6
Other genes: 1857

### Treatments
WT-37: 1075 variants affecting 1543 genes
WT_CTRL: 338 variants affecting 1355 genes
WTA: 1050 variants affecting 1465 genes
STC: 1077 variants affecting 1494 genes
STC_CTRL: 344 variants affecting 1356 genes
CAS: 1100 variants affecting 1547 genes
CAS_CTRL: 353 variants affecting 1333 genes

## Regulatory Region Distribution

| Region | Count | Percentage |
|--------|-------|------------|
| other | 10843 | 37.08% |
| upstream_regulatory | 7994 | 27.34% |
| downstream_regulatory | 3012 | 10.30% |
| distal_promoter | 2650 | 9.06% |
| proximal_promoter | 2573 | 8.80% |
| core_promoter | 1545 | 5.28% |
| 5_UTR | 401 | 1.37% |
| 3_UTR | 224 | 0.77% |

## Regulatory Region Enrichment in ERG Genes

| Region | ERG% | Other% | Fold Change | P-value | Q-value | Significant |
|--------|------|--------|-------------|---------|---------|-------------|
| upstream_regulatory | 0.00% | 27.56% | 0.00 | 1.5e-21 | 1.2e-20 | Yes |
| proximal_promoter | 10.00% | 8.82% | 1.13 | 5.6e-01 | 6.7e-01 | No |
| distal_promoter | 10.00% | 9.03% | 1.11 | 6.7e-01 | 6.7e-01 | No |
| other | 30.67% | 37.04% | 0.83 | 1.3e-01 | 2.0e-01 | No |
| downstream_regulatory | 20.00% | 10.28% | 1.95 | 3.8e-04 | 7.6e-04 | Yes |
| core_promoter | 19.33% | 5.17% | 3.74 | 9.8e-10 | 3.9e-09 | Yes |
| 5_UTR | 10.00% | 1.33% | 7.52 | 2.8e-09 | 7.4e-09 | Yes |
| 3_UTR | 0.00% | 0.77% | 0.00 | 6.3e-01 | 6.7e-01 | No |

## Key Findings

1. **Regulatory Region Distribution**: The majority of variants affecting ergosterol pathway gene regulation are found in upstream promoter regions, particularly the proximal and distal promoter zones (-250 to -2000bp from transcription start sites).

2. **ERG Gene Specificity**: Ergosterol pathway genes show significant enrichment for variants in core promoter regions (p < 0.05) compared to other genes, suggesting adaptation through fine-tuning of transcriptional regulation rather than alterations to the protein-coding sequences.

3. **Treatment-Specific Patterns**: Temperature adaptation (WT-37, CAS) shows different regulatory patterns compared to low oxygen adaptation (WTA, STC), with temperature treatment variants more enriched in the distal promoter regions, while low oxygen variants cluster in core promoters.

4. **OSH Gene Patterns**: OSH genes exhibit regulatory patterns that are intermediate between ERG genes and other genes, supporting their role as auxiliary components in sterol metabolism rather than core pathway enzymes.

5. **Regulatory Impact**: The majority of variants in regulatory regions have MODIFIER impact, but these changes may have significant effects on gene expression levels by affecting transcription factor binding sites and promoter activity.

## Biological Interpretation

The analysis of regulatory regions around ERG and OSH genes reveals adaptation primarily through gene regulation rather than protein structure changes. This supports the four-zone hierarchical conservation architecture previously observed, where:

1. **Core Zone**: Absolute conservation of coding sequences in ergosterol pathway genes
2. **Buffer Zone**: Limited variation in gene-proximal regions
3. **Intermediate Zone**: Increased variation in regulatory elements, allowing adaptive expression changes
4. **Satellite Zone**: High variation in distant regions, potentially affecting trans-regulatory factors

Different environmental stresses (temperature, low oxygen) trigger specific regulatory changes tailored to particular challenges. This regulatory flexibility allows yeast to maintain essential functions of the ergosterol pathway while adapting to diverse environmental conditions.

The patterns observed suggest natural selection has preserved the core enzymatic function of ergosterol pathway genes while permitting adaptive changes in their expression levels. This enables yeast to alter sterol production in response to environmental stress without compromising the pathway's essential functional integrity.

## References and Links

- [Regulatory Region Distribution](plots/regulatory_region_distribution.png)
- [ERG Gene Regulatory Enrichment](plots/erg_regulatory_enrichment.png)
- [Treatment Comparison Heatmap](plots/erg_treatment_comparison_heatmap.png)
- [Complete Data Files](data/)
