# Key Findings

This document summarizes the most significant findings from the Yeast MSA project.

## SNP Effect Analysis

The comprehensive analysis of variant effects using snpEff revealed important patterns:

### Variant Impact Distribution

- **90.2% MODIFIER impact**: Primarily intergenic and intronic variants (overall project statistics)
- **5.8% LOW impact**: Synonymous variants and other minor effects (overall project statistics)
- **3.4% MODERATE impact**: Primarily missense variants (overall project statistics)
- **0.6% HIGH impact**: Stop gained, frameshift, splice site disruptors (overall project statistics)

#### Treatment-Specific Scaffold Variant Impact Distribution

- **66.67% MODIFIER impact**: Primarily upstream gene variants
- **24.44% MODERATE impact**: Mostly missense variants
- **6.67% HIGH impact**: Frameshift variants with stop loss
- **2.22% LOW impact**: Synonymous variants

### Gene Region Distribution

- **62.7% Intergenic**
- **18.3% Exonic**
- **12.5% Intronic**
- **6.5% Upstream/Downstream**

### Treatment-Specific Patterns

- Temperature-adapted strains (WT-37: 40.00%, CAS: 33.33%) show higher proportion of variants than low oxygen-adapted strains (STC: 17.78%, WTA: 8.89%)
- Gene-modified strains show more HIGH impact variants in non-pathway genes
- Core ergosterol pathway genes show almost exclusively MODIFIER variants
- Distinct gene-specific patterns observed: ERG1 (22.22% of variants), ERG11 (20.00%), ERG24 (17.78%)

### Codon Changes

- Most common amino acid changes: Ala→Val, Ser→Asn, Gly→Asp
- Temperature adaptation shows enrichment for hydrophobic substitutions
- Low oxygen adaptation shows enrichment for charged residue substitutions

## Mutation Spectrum Findings

### Transition/Transversion Patterns

- Ti/Tv ratios range from 0.22 (WT-37) to 0.58 (CAS)
- Gene-modified strains show higher Ti/Tv ratios than non-modified strains
- Temperature and low oxygen adaptations show distinct patterns

### Mutation Type Preferences

- Temperature-adapted strains (WT-37 and CAS) show preference for C>A mutations
- Low oxygen-adapted strains (WTA and STC) show preference for C>G mutations
- Treatment-specific mutation preferences are statistically significant (p < 0.01)

### Adaptation-Specific Signatures

- Temperature adaptation signature enriched for TCT>TAT mutations
- Low oxygen adaptation signature enriched for GCG>GGG mutations
- Gene modification influences mutation signatures but not as strongly as adaptation type

## Genomic Context Findings

### Sequence Context

- GC content around mutations varies slightly between treatments (0.38-0.40)
- High proportion of mutations (~92%) occur near homopolymer regions
- Temperature adaptation shows enrichment for specific trinucleotide contexts

### Context Specificity

- C>A mutations in temperature adaptation occur preferentially in TCT contexts
- C>G mutations in low oxygen adaptation occur preferentially in GCG contexts
- Context preferences are statistically significant (p < 0.001)

## Population Structure Findings

### Sample Clustering

- Samples cluster primarily by adaptation type in PCA and MDS analyses
- Temperature and low oxygen adaptations form distinct genetic profiles
- Gene modification creates subclusters within adaptation types

### Variant Sharing

- High similarity within treatment replicates (Jaccard indices 0.60-0.73)
- Gene-modified strains show increased variant counts
- Variant sharing follows adaptation type more strongly than gene modification status

### Variant Counts

- Gene-modified strains show more variants than non-modified counterparts
  - CAS: 38 variants vs WT-37: 27 variants
  - STC: 27 variants vs WTA: 20 variants
- Controls show significantly fewer variants than treatments (4:1 ratio)

## Regional Enrichment Findings

### Enriched Regions

- Several genomic regions show significant enrichment for mutations
- Analysis identified 45 treatment-specific variants on 5 target scaffolds containing ergosterol pathway genes
- Scaffolds w303_scaffold_7 and w303_scaffold_8 show the highest number of variants (15 each)
- Some enriched regions are shared across treatments, suggesting common hotspots

### Adaptation-Specific Regions

- Temperature adaptation shows specific enrichment on certain scaffolds (WT-37: 40.00% of variants)
- Low oxygen adaptation shows comparatively fewer variants (WTA: 8.89% of variants)
- Gene-modified strains show intermediate variant counts (CAS: 33.33%, STC: 17.78%)
- Treatment-specific regions correlate with specific functional categories

## Hierarchical Conservation Architecture

One of the most significant findings is the discovery of a four-layered conservation pattern:

### 1. Core Zone (Absolute Conservation)

- Core ergosterol pathway genes show complete absence of HIGH/MODERATE impact variants
- No direct mutations in 11 key ergosterol genes across all treatments
- Strong statistical evidence of purifying selection (p < 0.0001)

### 2. Buffer Zone (Minimal Variation)

- Region extending ~5kb from core genes
- Predominantly regulatory (MODIFIER) variants when present
- 7 out of 11 ergosterol genes show zero variants within 5kb

### 3. Intermediate Zone (Moderate Variation)

- Region 5-50kb from core genes
- Moderate impact variants at consistent distances
- 97.78% of variants occur in this zone (5-10kb: 37.78%, 10-50kb: 60.00%)
- Gene modification effects visible in this zone

### 4. Satellite Zone (Adaptive Flexibility)

- Region 50-100kb from core genes
- HIGH impact variants clustered at specific distances
- Shows strongest treatment-specific patterns
- Contains "satellite genes" with regulatory connections to pathway

## Regulatory Adaptation Mechanism

### Upstream Variant Predominance

- 80% of variants near ergosterol genes are in upstream regulatory regions
- Few coding variants affecting protein structure directly
- Strong evidence that adaptation occurs through regulatory changes

### Satellite Gene Function

- Four "satellite genes" identified at consistent distances from pathway genes
- Annotated as "hypothetical proteins" in reference genome
- Correlate with specific sterol production changes
- Show treatment-specific high impact variants

## Sterol Profile Correlation

### Treatment-Specific Profiles

- Temperature adaptation shows significantly higher ergosterol levels (10.25 mean concentration)
- Low oxygen adaptation shows lower ergosterol (2.73 mean concentration, p=0.0109)
- Gene-modified strains have 1.5× more diverse sterol profiles than non-modified strains

### Unique Sterol Markers

- Temperature adaptation has 7 unique sterols
- Low oxygen adaptation has unique Tetrahymanol marker
- CAS-modified strain shows highest sterol diversity

### Genomic-Sterol Correlation

- Frameshift mutation downstream of ERG7 correlates with Tetrahymanol in STC strain
- High impact mutations 8kb upstream of ERG11 correlate with higher ergosterol in temperature-adapted strains
- Moderate missense mutations near ERG25 correlate with Stigmasta-5,22-dien-3-ol acetate

## Treatment Comparison Findings

### Treatment vs Control

- All treatments show significantly fewer variants compared to controls (p < 1e-70)
- Replicates show high consistency (70-72% of variants present in all three replicates)
- Treatment-specific variants cluster in similar functional categories

### Adaptation Comparison

- Temperature adaptation shows more variants than low oxygen adaptation
- Temperature-adapted strains have higher global variant density (0.0032 variants/kb) than low oxygen-adapted strains (0.0025 variants/kb)
- Temperature adaptation variants show stronger clustering

### Gene Modification Effects

- Gene-modified strains show more variants than non-modified counterparts
- CAS modification has stronger effect on variant patterns than STC modification
- Gene modifications influence variant distribution around pathway genes

## OSH Gene Family Analysis

The OSH (OxySterol binding Homology) gene family plays a critical role in sterol transport and membrane biology in yeast. Our analysis of this gene family has revealed important insights:

### OSH Gene Distribution

- 8 OSH family gene annotations identified in the W303 reference genome
- OSH genes are distributed across different scaffolds
- Some OSH genes are present on the same scaffolds as ERG pathway genes, enabling direct genomic distance analysis

### OSH Variant Analysis

- 140 variants identified within 25kb of OSH genes (compared to 344 variants near ERG genes)
- Average variants per OSH gene: 28.00 (similar to ERG genes: 31.27)
- OSH gene variants by treatment: WT-37 (27), WTA (26), STC (29), CAS (31)
- No significant difference in variant counts between OSH and ERG genes (Fisher's exact test p-value: 0.7860)

### OSH Conservation Patterns

- No variants found within OSH gene coding regions (0%), similar to ERG genes
- OSH genes: 13.6% of variants within 5kb (ERG genes: 43.6%)
- Higher proportion of HIGH impact variants near OSH genes (21.4%) compared to ERG genes (4.4%)
- Both gene families show evidence of strong purifying selection

### OSH-ERG Distance Analysis

- Minimum distance between OSH-ERG gene pairs: 742 bp
- Maximum distance: 390,857 bp
- Mean distance: 147,173 bp
- Some OSH-ERG pairs fall in the buffer zone (1-5,000 bp), supporting the four-zone architecture
- Notable close pairs: OSH3-ERG7 (742 bp), OSH7-ERG11 (12,900 bp)

## Integrated Model

The combined findings support a sophisticated model of adaptation:

1. **Conservation of Essential Function**: Core ergosterol pathway genes are protected from mutation through strong purifying selection
2. **Regulatory Flexibility**: Adaptation occurs primarily through changes in upstream regulatory regions
3. **Satellite Gene Modulation**: Specific "satellite genes" at consistent distances from pathway genes enable fine-tuning of the pathway
4. **OSH-ERG Coordination**: OSH genes show similar conservation patterns to ERG genes, suggesting coordinated evolution of sterol synthesis and transport
5. **Hierarchical Architecture**: A four-zone architecture balances conservation and adaptation
6. **Treatment-Specific Responses**: Different stress conditions trigger distinct genomic and biochemical responses

This model demonstrates how yeast balances the need to maintain essential membrane functions while adapting to environmental stress, providing insights into fundamental evolutionary mechanisms.