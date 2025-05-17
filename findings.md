
## Key Findings

The analyses have yielded several important findings:

### SNP Effect Analysis

Comprehensive SNP effect analysis was performed using snpEff to annotate and predict the functional effects of variants:

- **Variant Effect Distribution**:
  - MODIFIER impact: 90.2% (primarily intergenic and intronic variants)
  - LOW impact: 5.8% (synonymous variants and other minor effects)
  - MODERATE impact: 3.4% (missense variants primarily)
  - HIGH impact: 0.6% (stop gained, frameshift, splice site disruptors)

- **Gene Region Distribution**:
  - Intergenic: 62.7% 
  - Exonic: 18.3%
  - Intronic: 12.5%
  - Upstream/Downstream: 6.5%

- **Treatment-Specific Patterns**:
  - Temperature-adapted strains show higher proportion of MODERATE impact variants
  - Gene-modified strains show more HIGH impact variants in non-pathway genes
  - Core ergosterol pathway genes show almost exclusively MODIFIER variants

- **Codon Changes**:
  - Most common amino acid changes: Ala→Val, Ser→Asn, Gly→Asp
  - Temperature adaptation shows enrichment for hydrophobic substitutions
  - Low oxygen adaptation shows enrichment for charged residue substitutions

### 1. Mutation Spectrum

- Different treatments show distinct mutation spectra, with transition/transversion (Ti/Tv) ratios ranging from 0.22 (WT-37) to 0.58 (CAS)
- Temperature-adapted strains (WT-37 and CAS) show a preference for C>A mutations
- Low oxygen-adapted strains (WTA and STC) show a preference for C>G mutations
- Gene modifications influence mutation spectra, with gene-modified strains showing higher Ti/Tv ratios

### 2. Genomic Context

- GC content around mutations varies slightly between treatments (0.38-0.40)
- A high proportion of mutations occur near homopolymer regions (~92%)
- Temperature adaptation shows enrichment for specific trinucleotide contexts (TCT with C>A mutations)
- Gene-modified strains show different sequence context preferences compared to non-modified strains

### 3. Population Structure

- Samples cluster primarily by adaptation type in PCA and MDS analyses
- Temperature and low oxygen adaptations show distinct genetic profiles
- Within treatments, replicates show high similarity (Jaccard indices 0.60-0.73)
- Gene-modified strains show increased variant counts compared to non-modified strains

### 4. Regional Enrichment

- Several genomic regions show significant enrichment for mutations under specific treatments
- Regions on scaffolds CM007977.1 and CM007967.1 show particularly high fold enrichments
- Some enriched regions are shared across multiple treatments, suggesting common mutational hotspots
- Adaptation-specific regions have been identified for both temperature and low oxygen conditions

### 5. Scaffold Distribution

- Mutations are non-randomly distributed across the genome
- Small scaffolds (CM007980.1, LYZE01000020.1) show high variant density
- Treatment correlation analysis shows highest similarity between WT-37 and STC (ρ = 0.53)
- Temperature-adapted strains show higher global variant density (0.0032 variants/kb) than low oxygen-adapted strains (0.0025 variants/kb)

### 6. Treatment vs Control

- All treatments show significantly fewer variants compared to controls (p < 1e-70)
- Gene-modified strains show more variants than their non-modified counterparts (CAS: 38 vs WT-37: 27; STC: 27 vs WTA: 20)
- Replicates show high consistency, with 70-72% of variants present in all three replicates

### 7. Gene Analysis

- Strong evidence of purifying selection in ergosterol pathway genes
- No HIGH or MODERATE impact variants found in core ergosterol pathway genes
- Gene modifications influence the distribution of variants around pathway genes
- Non-genic regions show higher mutation rates than genic regions

### 8. Network Analysis

- Extended ergosterol network reveals adaptation through changes in neighboring genes
- Gene-gene interactions suggest coordinated adaptation response
- Treatment-specific network patterns identified
- HIGH impact variants consistently found at specific distances from pathway genes

### 9. OSH Gene Family Analysis

- 8 OSH family gene annotations identified in the W303 reference genome
- 140 variants identified within 25kb of OSH genes across all treatments
- No variants found within OSH gene coding regions, indicating strong purifying selection similar to ERG genes
- OSH gene variants exhibit balanced distribution across treatments: WT-37 (27), WTA (26), STC (29), CAS (31)
- Higher proportion of HIGH impact variants near OSH genes (21.4%) compared to ERG genes (4.4%)
- Close proximity of certain OSH-ERG pairs identified (e.g., OSH3-ERG7 at 742 bp, OSH7-ERG11 at ~13kb)
- OSH-ERG gene relationships support the four-zone conservation architecture
- OSH gene analysis suggests coordinated evolution between sterol synthesis and transport mechanisms

### 10. Sterol Profile Analysis

- Temperature adaptation shows significantly higher ergosterol levels (10.25 mean concentration) compared to low oxygen adaptation (2.73 mean concentration, p=0.0109)
- Gene-modified strains have 1.5× more diverse sterol profiles than non-modified strains
- Treatment-specific sterol signatures: Temperature adaptation has 7 unique sterols; low oxygen adaptation has unique Tetrahymanol marker
- Sterol compositions validate the hierarchical conservation pattern found in genomic analysis
- "Satellite genes" at consistent distances (7-50kb) from ergosterol pathway genes correlate with specific sterol production
- Evidence of adaptation through altered sterol profiles despite perfect conservation of pathway genes
- Integration reveals a four-layered conservation architecture from absolute conservation (core) to adaptive flexibility (satellite zone)

## Visualizations and Reports

The project generates a rich variety of reports and visualizations that integrate findings across multiple analysis modules. These outputs provide different views of the data, from raw statistics to interactive visualizations.

### Analysis Reports

Text-based summary reports provide detailed analysis of findings:

1. **Core Analysis Summaries**:
   - `combined_analysis_results.txt`: Master summary of all findings across modules
   - `analysis/mutation_spectrum_results/statistical_test_results.txt`: Statistical analysis of mutation patterns
   - `analysis/mutational_signatures_results/mutational_signatures_summary.txt`: Signature analysis findings
   - `analysis/population_structure_results/population_structure_summary.txt`: Population analysis results
   - `analysis/regional_enrichment_results/regional_enrichment_summary.txt`: Enrichment findings
   - `analysis/statistical_pattern_results/statistical_analysis_summary.txt`: Statistical analysis report

2. **Pipeline Processing Reports**:
   - `results/vcf/analysis_report.txt`: Summary of variant statistics across treatments
   - `results/scaffold_variants/scaffold_variant_summary.txt`: Analysis of variants by scaffold
   - `results/gene_variants/gene_summary.tsv`: Summary of gene-specific variant patterns
   - `results/treatment_analysis/analysis_summary.txt`: Treatment comparison summary
   - `vcf/annotated/annotation_summary.txt`: Summary of variant annotations

3. **Sterol Analysis Reports**:
   - `results/sterol_analysis/basic_stats/sterol_statistics.txt`: Basic sterol composition statistics
   - `results/sterol_analysis/comparative/comparative_analysis_summary.txt`: Treatment comparisons
   - `results/sterol_analysis/correlation/statistical_correlation_results.txt`: Genomic-sterol correlations
   - `results/sterol_analysis/pathway/pathway_analysis_summary.txt`: Sterol pathway flux analysis

### Interactive HTML Reports

Interactive HTML dashboards provide rich visualizations and explorable data:

1. **Ergosterol Variant Analysis Dashboard**:
   - **Path**: `results/reports/ergosterol_variant_analysis.html`
   - **Generator**: `scripts/utils/generate_ergosterol_variant_report.py`
   - **Contents**:
     - Interactive variant distribution visualizations
     - Distance-based analysis from pathway genes
     - Treatment-specific variant patterns
     - Purifying selection evidence
     - Biological significance interpretation
   - **Features**:
     - Interactive image gallery with zoom capability
     - Tabbed visualization sections
     - Dark/light mode toggle
     - Mobile-responsive design

2. **Functional Impact Analysis Dashboard**:
   - **Path**: `results/reports/functional_impact.html`
   - **Generator**: `scripts/utils/generate_functional_impact_report.py`
   - **Contents**:
     - Protein domain impact visualizations
     - Conservation patterns analysis
     - Satellite gene architecture diagrams
     - Adaptation mechanisms models
     - Hierarchical conservation zone visualizations

3. **Sterol Profile Analysis Dashboard**:
   - **Path**: `results/reports/sterols.html`
   - **Generator**: `scripts/sterols/generate_html_report.py`
   - **Contents**:
     - Interactive sterol profile visualizations
     - Treatment comparison heatmaps
     - Pathway visualizations with flux indicators
     - Genomic-sterol integration diagrams
     - Adaptation model visualizations

4. **Variant Analysis Dashboard**:
   - **Path**: `results/reports/variant_analysis.html`
   - **Generator**: `scripts/variants/generate_variant_report.py`
   - **Contents**:
     - Sample comparison visualizations
     - Interactive filtering of variants
     - Annotation statistics and distribution
     - Genome browser-like variant visualization
     - Treatment-specific variant patterns

### Markdown Reports

Detailed narrative analysis with embedded visualizations:

1. **Network Analysis Report**:
   - **Path**: `results/network_analysis/network_analysis_report.md`
   - **Contents**:
     - Extended ergosterol network analysis
     - Subnetwork interactions
     - Treatment-specific network patterns
     - Pathway distance models
     - Network statistics and interpretations

2. **OSH Gene Analysis Report**:
   - **Path**: `results/osh_analysis/OSH_Results.md`
   - **Contents**:
     - OSH gene family mapping and analysis
     - OSH variant patterns and conservation
     - OSH-ERG gene distance analysis
     - Biological interpretation of OSH gene findings
     - Integration with four-zone conservation model
     - Visualizations of OSH gene relationships

3. **Integrated Findings Report**:
   - **Path**: `results/sterol_analysis/correlation/integrated_findings_report.md`
   - **Contents**:
     - Integration of sterol and genomic findings
     - Hierarchical conservation model explanation
     - Satellite gene-sterol connections
     - Four-zone conservation architecture
     - Comprehensive adaptation model

4. **High Impact Variants Report**:
   - **Path**: `results/functional_impact/high_impact/high_impact_variants_report.md`
   - **Contents**:
     - Analysis of high impact variants
     - Functional domain mappings
     - Evolutionary conservation patterns
     - Treatment-specific functional changes
     - Pathway impact assessment

### Gene Mapping and Reference Data

Reference data files providing comprehensive annotations:

1. **Gene Mapping Full**:
   - **Path**: `reference/gene_mapping_full.tsv`
   - **Generator**: `scripts/general_gene_analysis/generate_gene_mapping_full.py`
   - **Contents**:
     - W303 gene ID to standard gene name mapping
     - Ergosterol pathway annotations
     - Gene function categories
     - Chromosome and position information
     - Conservation zone annotations

2. **Genes of Interest Mapping**:
   - **Path**: `reference/genes_of_interest_mapping.tsv`
   - **Contents**:
     - Ergosterol pathway genes
     - Functional categorization
     - Pathway position annotations
     - Distance relationships
     - Regulatory connections

3. **Chromosome Mapping**:
   - **Path**: `reference/chromosome_mapping.tsv`
   - **Contents**:
     - W303 scaffold IDs to standard chromosome names
     - Size and coordinate information
     - GC content and feature density
     - Annotation status
     - Reference cross-mapping

### Visualization Galleries

Organized collections of visualizations by analysis type:

1. **Mutation Spectrum Analysis**:
   - **Directory**: `analysis/mutation_spectrum_results/`
   - **Key Visualizations**:
     - `mutation_spectrum_summary.csv`: Numeric data on mutation types
     - `comparative_mutation_spectrum.png`: Cross-treatment comparison
     - `ti_tv_ratios_by_adaptation.png`: Transition/transversion by adaptation
     - Treatment-specific spectra (`CAS_mutation_spectrum.png`, etc.)

2. **Genomic Context Analysis**:
   - **Directory**: `analysis/genomic_context_results/`
   - **Key Visualizations**:
     - `gc_content_by_adaptation.png`: GC content patterns
     - `homopolymer_by_gene_status.png`: Homopolymer analysis
     - `mutation_type_by_treatment_heatmap.png`: Context heatmap
     - Sequence logos for specific contexts (e.g., `logo_Temperature_C_to_A.png`)

3. **Population Structure Analysis**:
   - **Directory**: `analysis/population_structure_results/`
   - **Key Visualizations**:
     - `pca_by_adaptation.png`: Principal component analysis
     - `mds_by_treatment.png`: Multi-dimensional scaling
     - `dendrogram_by_adaptation.png`: Hierarchical clustering
     - `shared_variants_heatmap.png`: Variant sharing patterns

4. **Regional Enrichment Analysis**:
   - **Directory**: `analysis/regional_enrichment_results/`
   - **Key Visualizations**:
     - `enrichment_heatmap.png`: Regional enrichment patterns
     - `clustered_enrichment_heatmap.png`: Clustered enrichment analysis
     - Treatment-specific enrichment results (CSV files)
     - Adaptation-specific enrichment summaries

5. **Scaffold Distribution Analysis**:
   - **Directory**: `analysis/scaffold_distribution_results/`
   - **Key Visualizations**:
     - `comparative_density_heatmap_top30.png`: Density comparison
     - `treatment_correlation_heatmap.png`: Treatment correlations
     - Treatment-specific bubble charts and variant density plots
     - `adaptation_specific_hotspots.txt`: Hotspot locations

6. **Gene Analysis Visualizations**:
   - **Directory**: `analysis/genes_of_interest/treatment_control_analysis/`
   - **Key Visualizations**:
     - `erg_gene_distribution.png`: Ergosterol gene variant distribution
     - `gene_status_distribution.png`: Gene status distribution
     - Treatment-specific purifying selection plots
     - `fold_change_by_gene_status.png`: Gene expression fold changes

7. **Variant Proximity Impact Analysis**:
   - **Directory**: `results/filtered_scaffold_variants/impact/`
   - **Key Visualizations**:
     - `variant_count_heatmap.png`: Heatmap of variant counts by ERG gene and treatment
     - `distance_distribution_heatmap.png`: Distribution of distances to ERG genes
     - `variant_proximity_impact_summary.html`: Interactive summary of variant proximity and impact
     - `variant_proximity_impact_summary.tsv`: Tabular data of variant proximity and impact

8. **OSH Gene Analysis Visualizations**:
   - **Directory**: `results/osh_analysis/`
   - **Key Visualizations**:
     - `osh_erg_variant_comparison.png`: Comparison of variants near OSH and ERG genes
     - `osh_erg_distance_distribution.png`: Distribution of distances between OSH and ERG genes
     - `osh_erg_network.png`: Network visualization of OSH-ERG relationships
     - `osh_treatment_heatmap.png`: Heatmap of OSH gene variants by treatment

9. **Sterol Profile Visualizations**:
   - **Directory**: `results/sterol_analysis/visualizations/`
   - **Key Visualizations**:
     - `sterol_composition_by_treatment.png`: Treatment comparison
     - `pathway_flux_visualization.png`: Pathway analysis
     - `satellite_gene_sterol_correlation.png`: Gene-sterol correlations
     - `four_zone_conservation_model.png`: Conservation architecture

These reports and visualizations collectively provide a comprehensive view of the project's findings, from basic statistics to complex models of adaptation mechanisms. The interactive HTML dashboards offer user-friendly exploration of the data, while the text and markdown reports provide detailed interpretations and biological significance.

All visualization outputs follow consistent color schemes and formatting to facilitate cross-analysis comparisons. The HTML reports are generated using reusable visualization templates that integrate Bootstrap for responsive design and D3.js for interactive elements.
