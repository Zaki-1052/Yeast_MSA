# Comparison Report: Scaffold Variant vs. Gene Variant Analysis

This report compares the findings from two different regulatory analyses:

1. **New Regulatory Analysis**: Based on treatment-specific scaffold variants (`results/filtered_scaffold_variants/treatment_specific_scaffold_variants.tsv`)
2. **Specific Regulatory Analysis**: Based on gene variants (`results/gene_variants/all_gene_variants.tsv`)

## Overview Comparison

| Aspect | New Regulatory Analysis | Specific Regulatory Analysis |
|--------|------------------------|------------------------------|
| Total variants | 45 | 44 |
| Variant types | Mix of coding (33.3%) and regulatory (66.7%) | 100% regulatory |
| Distance range | 3,034 - 44,512 bp | 89 - 332 bp |
| Mean distance | 21,841.8 bp | 180.4 bp |
| Conservation zones | Primarily intermediate zone | 100% buffer zone |

## Key Differences

### 1. Variant Data Source

- **New Regulatory Analysis**: Uses scaffold variants that include both coding and regulatory variants found at both close and distant locations from genes
- **Specific Regulatory Analysis**: Uses gene variants that include only upstream variants very close to genes

### 2. Variant Distribution

- **New Regulatory Analysis**:
  - Coding missense: 24.4%
  - Coding synonymous: 2.2%
  - Coding frameshift: 6.7%
  - Distal regulatory: 66.7%

- **Specific Regulatory Analysis**:
  - Distal regulatory: 100%

### 3. Distance from Genes

The most striking difference is the distance of variants from genes:

- **New Regulatory Analysis**:
  - Minimum distance: 3,034 bp
  - Maximum distance: 44,512 bp
  - Mean distance: 21,841.8 bp

- **Specific Regulatory Analysis**:
  - Minimum distance: 89 bp
  - Maximum distance: 332 bp
  - Mean distance: 180.4 bp

This ~100x difference in mean distance indicates that the two analyses are examining fundamentally different types of genetic variants.

### 4. Treatment-Specific Patterns

- **New Regulatory Analysis**:
  - Shows clear differences between treatments
  - Temperature adaptation (WT-37): 38.9% missense variants
  - Low oxygen adaptation (WTA): 0% missense variants

- **Specific Regulatory Analysis**:
  - No treatment-specific differences in variant types
  - All treatments show 100% regulatory variants

## Biological Interpretation of Differences

The stark differences between these two analyses reveal important insights about the genetic architecture of adaptation:

1. **Two Distinct Adaptation Mechanisms**:
   - **Proximal Regulation**: The gene variants analysis (specific) shows consistent regulatory changes very close to genes (within 332 bp)
   - **Distal Adaptation**: The scaffold variants analysis (new) reveals a broader pattern of both coding and regulatory changes at much greater distances

2. **Conservation Zone Architecture**:
   - Gene variants occur exclusively in the buffer zone (very close to genes)
   - Scaffold variants occur primarily in the intermediate zone (further from genes)
   - This supports a model where different types of genetic changes occur at different distances from essential genes

3. **Treatment-Specific Strategies**:
   - Core promoter regulation (revealed in the specific analysis) appears similar across all treatments
   - Distal and coding changes (revealed in the new analysis) show clear treatment-specific patterns

## Conclusions

These two analyses are not contradictory but rather complementary, revealing different layers of adaptation:

1. **Layer 1: Core Regulation** (89-332 bp)
   - Consistent across treatments
   - Affects promoter elements
   - May represent universal adaptation mechanisms

2. **Layer 2: Distal Adaptation** (3,034-44,512 bp)
   - Varies by treatment
   - Includes both coding and regulatory changes
   - May represent environment-specific adaptations

This two-layer model provides a more comprehensive understanding of how yeast adapts to different stresses, showing that adaptation involves both close regulatory changes and more distant genomic modifications that vary by environmental challenge.

These findings strongly support the hierarchical conservation architecture model, showing that genomic changes follow specific patterns at different distances from essential genes.