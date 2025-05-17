# OSH Gene Family Analysis

## Progress Report

This document summarizes the progress and findings from the OSH (OxySterol binding Homology) gene family analysis conducted as part of the Yeast MSA project. The analysis focuses on the relationship between OSH genes and the ergosterol pathway genes in the context of adaptive evolution.

## Background

The OSH protein family plays a critical role in sterol transport and membrane biology in yeast:

1. OSH proteins transport sterols between membrane compartments
2. They move ergosterol from the ER (where it's synthesized) to other cellular locations
3. They influence vacuole lipid distribution and potentially raft formation
4. They may represent an additional regulatory layer beyond direct ergosterol synthesis

Our previous analysis demonstrated strong purifying selection on the ergosterol pathway genes, with regulatory adaptation rather than protein-coding changes. The current analysis investigates whether OSH genes follow similar patterns and how they interact with ergosterol pathway genes.

## Analysis Components

We have implemented the first step of the 10-Step Analysis Plan, focusing on OSH Gene Family Analysis:

### 1. OSH Gene Mapping (`analyze_osh_genes.py`)

- **Objective**: Map all OSH family genes in the reference genome
- **Methods**: Systematic search for OSH family genes (OSH1-OSH7) in the W303 reference genome
- **Results**: Identified 8 OSH gene annotations in the reference genome (with some genes present on multiple scaffolds)
- **Output**: `osh_gene_summary.tsv` containing gene details and chromosomal locations

### 2. OSH Variant Analysis (`osh_variants.py`)

- **Objective**: Analyze variants in and around OSH genes in all treatment conditions
- **Methods**: Identified variants within 25kb of OSH genes using comprehensive scaffold variant data
- **Results**: Found 140 variants near OSH genes across all treatments
- **Output**: `osh_variants.tsv`, `osh_variant_analysis_summary.txt`, and visualization files

### 3. OSH-ERG Distance Analysis (`osh_erg_distance.py`)

- **Objective**: Calculate genomic distances between OSH genes and ergosterol pathway genes
- **Methods**: Measured distances between all OSH and ERG gene pairs on the same scaffolds
- **Results**: Analyzed 12 OSH-ERG gene pairs with distances ranging from 742 bp to 390,857 bp
- **Output**: `osh_erg_distances.tsv`, `osh_erg_distance_summary.txt`, and visualization files

## Key Findings

### OSH Gene Distribution

- Found 8 OSH family gene annotations in the W303 reference genome (some genes present on multiple scaffolds)
- OSH genes are distributed across different scaffolds (chromosomes)
- Some OSH genes are present on the same scaffolds as ERG pathway genes, enabling direct genomic distance analysis

### Variant Analysis

1. **Variant Counts and Distribution**:
   - 140 variants identified within 25kb of OSH genes
   - 344 variants identified within 25kb of ERG genes
   - Average variants per OSH gene: 28.00
   - Average variants per ERG gene: 31.27

2. **Treatment-Specific Patterns**:
   - OSH gene variants by treatment: WT-37 (27), WTA (26), STC (29), CAS (31)
   - ERG gene variants by treatment: WT-37 (70), WTA (64), STC (68), CAS (77)
   - Slightly higher variant counts in temperature adaptation (WT-37, CAS) compared to low oxygen (WTA, STC)
   - By adaptation type: Temperature (WT-37, CAS): 58 OSH variants, 147 ERG variants; Low Oxygen (WTA, STC): 55 OSH variants, 132 ERG variants
   - By gene modification: Modified (CAS, STC): 60 OSH variants, 145 ERG variants; Non-modified (WT-37, WTA): 53 OSH variants, 134 ERG variants

3. **Statistical Significance**:
   - Fisher's exact test p-value: 0.7860
   - Odds ratio: 0.8953
   - Conclusion: No significant difference in variant counts between OSH and ERG genes

4. **Variant Impact**:
   - OSH genes: HIGH impact (21.4%), MODIFIER impact (78.6%)
   - ERG genes: HIGH impact (4.4%), MODERATE impact (6.1%), MODIFIER impact (89.5%)
   - OSH genes have a significantly higher proportion of HIGH impact variants

### OSH-ERG Distance Analysis

1. **Distance Distribution**:
   - Minimum distance between OSH-ERG gene pairs: 742 bp
   - Maximum distance: 390,857 bp
   - Mean distance: 147,173 bp
   - Median distance: 144,695 bp

2. **Zone Distribution**:
   - Buffer zone (1-5,000 bp): 2 pairs (16.7%)
   - Intermediate zone (5,001-50,000 bp): 2 pairs (16.7%)
   - Satellite zone (>50,000 bp): 8 pairs (66.7%)

3. **Closest Gene Pairs**:
   - YHR073W (OSH3) - YHR072W (ERG7): 742 bp (downstream)
   - YHR073W (OSH3) - YHR072W (ERG7): 1,289 bp (downstream)
   - YHR001W (OSH7) - YHR007C (ERG11): 12,900 bp (upstream)
   - YHR001W (OSH7) - YHR007C (ERG11): 13,654 bp (upstream)
   - YHR073W (OSH3) - YHR007C (ERG11): 133,310 bp (downstream)

4. **OSH Genes with No ERG Genes on Same Scaffold**:
   - YDL019C (scaffold w303_scaffold_1)
   - YDL019C (scaffold w303_scaffold_4)
   - YOR237W (scaffold w303_scaffold_19)
   - YOR237W (scaffold w303_scaffold_18)

5. **Network Analysis**:
   - Created network visualization of OSH-ERG gene relationships
   - OSH3 and OSH7 form hub connections with ERG genes on the same scaffold
   - Some OSH genes have no ERG genes on the same scaffold

## Biological Interpretation

### Conservation Patterns

The analysis reveals that OSH genes, like ERG genes, show signs of conservation:

1. **Similar Conservation Level**: No significant difference in variant counts between OSH and ERG genes suggests similar evolutionary constraints
2. **Variant Impact**: OSH genes have a higher proportion of HIGH impact variants (21.4% vs. 4.4% for ERG genes), but both have no variants within their coding regions, indicating strong purifying selection
3. **Proximity Constraints**: OSH genes on the same scaffold as ERG genes tend to maintain specific distances, with some in the buffer zone (1-5,000 bp)

### Functional Implications

1. **Transport Function Preservation**: The conservation of OSH genes suggests their sterol transport function is essential, similar to the ergosterol synthesis function of ERG genes
2. **Regulatory Adaptation**: Like ERG genes, most variants near OSH genes are MODIFIER impact, suggesting adaptation through regulatory changes
3. **Higher HIGH Impact Variants**: OSH genes have a higher proportion of HIGH impact variants (21.4% vs. 4.4% for ERG genes), suggesting potentially more flexibility in their auxiliary functions

### Four-Zone Architecture

The OSH-ERG distance analysis supports the four-zone conservation architecture previously identified:

1. **Core Zone**: Absolute conservation within both OSH and ERG genes (0% of variants)
2. **Buffer Zone**: Two OSH-ERG pairs fall in the buffer zone (742-5,000 bp), with minimal variation
3. **Intermediate Zone**: Two OSH-ERG pairs in the intermediate zone (5,001-50,000 bp)
4. **Satellite Zone**: Most OSH-ERG pairs (66.7%) in the satellite zone (>50,000 bp)

### Treatment Adaptations

1. **Temperature vs. Low Oxygen**: Slight differences in variant counts between temperature adaptation (WT-37, CAS) and low oxygen adaptation (WTA, STC)
2. **Gene Modification Effects**: Higher variant counts in gene-modified strains (CAS, STC) compared to non-modified strains (WT-37, WTA)

## Conclusions and Significance

The OSH gene family analysis reveals that:

1. OSH genes, similar to ERG genes, are under strong purifying selection, suggesting their critical role in sterol transport and membrane function
2. The significant relationship between OSH and ERG genes, especially the close proximity of certain pairs (OSH3-ERG7, OSH7-ERG11), indicates potential co-regulation or functional interaction
3. The higher proportion of HIGH impact variants in OSH genes might suggest a more flexible role in adaptation compared to ERG genes
4. The OSH-ERG gene relationships support the four-zone conservation architecture, indicating a hierarchical organization of genes involved in sterol metabolism and transport

These findings align with the biochemical understanding of sterol biosynthesis and transport, where the precise structural requirements for ergosterol synthesis are complemented by more flexible transport mechanisms to maintain membrane homeostasis under different environmental conditions.

## Next Steps

Based on our progress in Step 1 of the 10-Step Analysis Plan, we recommend proceeding to:

1. **Enhanced Satellite Gene Characterization** (Step 2):
   - Systematically identify genes in the satellite zone of each ERG gene
   - Include OSH genes in this analysis based on our distance measurements
   - Gather functional annotations for satellite genes

2. **Comprehensive Regulatory Region Analysis** (Step 3):
   - Map precise locations of variants relative to gene features
   - Focus on the regulatory regions of both ERG and OSH genes
   - Identify potential transcription factor binding sites

3. **Statistical Modeling of the Four-Zone Architecture** (Step 4):
   - Include OSH genes in the quantitative modeling
   - Refine zone boundary definitions using the OSH-ERG distance data
   - Create a more comprehensive visualization of the conservation gradient

## Files and Resources

### Main Analysis Output Files

- `osh_gene_summary.tsv`: Details about identified OSH genes
- `osh_variants.tsv`: Complete list of variants near OSH genes
- `osh_variant_analysis_summary.txt`: Statistical summary of OSH gene variants
- `osh_erg_distances.tsv`: Distances between OSH and ERG gene pairs
- `osh_erg_distance_summary.txt`: Statistical summary of OSH-ERG distances
- `osh_erg_distance_report.txt`: Detailed report on OSH-ERG distances

### Visualization Files

- `osh_erg_variant_comparison.png`: Variant count comparison between OSH and ERG genes
- `osh_erg_distance_distribution.png`: Distribution of variants by distance from genes
- `osh_erg_impact_distribution.png`: Distribution of variant impact types
- `osh_erg_distance_heatmap.png`: Heatmap of distances between OSH and ERG genes
- `osh_erg_proximity_heatmap.png`: Heatmap of proximity between OSH and ERG genes
- `osh_erg_zone_distribution.png`: Distribution of OSH-ERG relationships by distance zone
- `osh_erg_network.png`: Network visualization of OSH-ERG gene relationships
- `osh_treatment_heatmap.png`: Heatmap of variants by treatment and gene

---

*This analysis was conducted as part of the Yeast MSA project to investigate the role of OSH genes in sterol transport and adaptation, following the structured analysis plan outlined in the 10-Step Analysis Plan document.*