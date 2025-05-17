# New Regulatory Analysis Planning Document

## Background and Rationale

The current regulatory analysis approach has provided valuable insights but has several limitations that need to be addressed in a new, more comprehensive framework. The observed phenomenon of 80% of variants near ergosterol genes being found in regulatory regions rather than protein-coding sequences suggests that adaptation occurs primarily through changes in gene expression rather than protein structure. This requires a more sophisticated analytical approach to fully understand the regulatory mechanisms involved in adaptation.

## Objectives of the New Regulatory Analysis

1. **Precise Mapping of Regulatory Features**: Map variants to specific regulatory elements with higher resolution and biological context
2. **Transcription Factor Binding Site (TFBS) Analysis**: Identify affected TFBS and their biological significance
3. **Regulatory Network Modeling**: Construct regulatory networks connecting ERG and OSH genes
4. **Adaptation-Specific Regulatory Patterns**: Identify regulatory changes specific to temperature vs. low oxygen adaptation
5. **Statistical Validation**: Apply rigorous statistical analysis to establish significance of regulatory patterns
6. **Integration with Four-Zone Architecture**: Connect regulatory findings with the hierarchical conservation model

## Components of the New Regulatory Analysis Framework

### 1. Enhanced Regulatory Region Definition Framework

**Files**:
- `new_regulatory_analysis/regulatory_region_definition.py`
- `new_regulatory_analysis/data/regulatory_features.json`

**Purpose**:
- Define regulatory regions with higher precision based on yeast genomic features
- Incorporate known regulatory elements from literature and databases
- Create a comprehensive map of potential regulatory regions around genes

**Key Improvements**:
- Include better definitions for promoters, 5'UTRs, 3'UTRs
- Consider nucleosome positioning data to identify accessible regions
- Account for known transcription factor binding motifs
- Define enhancers and silencers based on current yeast literature

### 2. Advanced Variant-to-Regulatory Feature Mapping System

**Files**:
- `new_regulatory_analysis/variant_regulatory_mapping.py`
- `new_regulatory_analysis/data/variant_regulatory_annotations.tsv`

**Purpose**:
- Map variants to specific regulatory features with higher precision
- Categorize variants by their potential regulatory impact
- Quantify distance relationships to transcription start sites

**Key Improvements**:
- Implement position-weight matrices for transcription factor binding sites
- Calculate disruption scores for variants in binding sites
- Consider evolutionary conservation of regulatory elements
- Account for chromatin accessibility data where available

### 3. Treatment-Specific Regulatory Pattern Analysis

**Files**:
- `new_regulatory_analysis/treatment_regulatory_analysis.py`
- `new_regulatory_analysis/data/treatment_regulatory_patterns.tsv`

**Purpose**:
- Identify regulatory patterns specific to different adaptation conditions
- Compare regulatory changes between temperature and low oxygen adaptation
- Analyze the effect of gene modifications (CAS, STC) on regulatory patterns

**Key Improvements**:
- Statistical testing for treatment-specific enrichment
- Multiple testing correction for robust significance assessment
- Effect size calculations for regulatory changes
- Integration with experimental evidence of adaptation

### 4. ERG-OSH Regulatory Relationship Analysis

**Files**:
- `new_regulatory_analysis/erg_osh_regulation.py`
- `new_regulatory_analysis/data/erg_osh_regulatory_relationships.tsv`

**Purpose**:
- Analyze regulatory relationships between ERG and OSH genes
- Identify shared regulatory elements and potential co-regulation
- Investigate how regulatory changes in one group affect the other

**Key Improvements**:
- Comparison of regulatory variant patterns between gene families
- Identification of shared transcription factors
- Distance-based analysis of regulatory relationships
- Integration with findings from OSH gene analysis

### 5. Statistical Validation Framework

**Files**:
- `new_regulatory_analysis/statistical_validation.py`
- `new_regulatory_analysis/data/regulatory_statistics.json`

**Purpose**:
- Implement rigorous statistical tests for regulatory findings
- Account for multiple testing and potential biases
- Calculate effect sizes and confidence intervals

**Key Improvements**:
- Permutation tests for enrichment analysis
- Bootstrap confidence intervals for effect size estimates
- Bayesian analysis of regulatory variant significance
- Simulation-based validation of observed patterns

### 6. Regulatory Impact Score System

**Files**:
- `new_regulatory_analysis/regulatory_impact_scoring.py`
- `new_regulatory_analysis/data/variant_impact_scores.tsv`

**Purpose**:
- Develop a comprehensive scoring system for regulatory impact
- Integrate position, conservation, and motif disruption into a single score
- Prioritize variants by their potential regulatory significance

**Key Improvements**:
- Position-specific scoring relative to transcription start sites
- Motif disruption scoring based on information content
- Conservation-weighted impact calculations
- Integration of experimental evidence where available

### 7. Four-Zone Integration Analysis

**Files**:
- `new_regulatory_analysis/four_zone_integration.py`
- `new_regulatory_analysis/data/zone_regulatory_distribution.tsv`

**Purpose**:
- Connect regulatory findings with the four-zone conservation architecture
- Analyze regulatory variant distribution across conservation zones
- Test hypotheses about zone-specific regulatory mechanisms

**Key Improvements**:
- Zone-specific regulatory element enrichment
- Comparison of regulatory impact across zones
- Integration of regulatory findings with distance-based conservation model
- Visualization of regulatory patterns in the context of conservation zones

### 8. Visualization and Reporting System

**Files**:
- `new_regulatory_analysis/visualization/regulatory_visualizations.py`
- `new_regulatory_analysis/reports/regulatory_analysis_report.md`

**Purpose**:
- Create comprehensive visualizations of regulatory findings
- Generate integrated reports connecting all aspects of the analysis
- Provide interpretable summaries for biological significance

**Key Improvements**:
- Interactive visualizations of regulatory patterns
- Integration of statistical significance into visualizations
- Comprehensive reporting connecting all analysis components
- Biological interpretation of regulatory findings

## Implementation Strategy and Current Status

### Phase 1: Framework Development - ✅ COMPLETED
- ✅ Develop the enhanced regulatory region definition framework
  - Implemented in `scripts/new_regulatory_analysis/regulatory_region_definition.py`
  - Successfully created detailed regulatory region annotations
  - Generated comprehensive reports and visualizations
  
- ✅ Create the advanced variant-to-regulatory feature mapping system
  - Implemented in `scripts/new_regulatory_analysis/variant_regulatory_mapping.py`
  - Successfully mapped 44 variants to specific regulatory features
  - Generated visualization of regulatory patterns and mapping report
  
- ⏳ Implement the statistical validation framework
  - Partially implemented in treatment-specific pattern analysis
  - Remaining components to be developed in Phase 2
  
- ⏳ Build the regulatory impact score system
  - To be implemented in Phase 2

### Phase 2: Analysis Implementation - 🔄 IN PROGRESS
- ✅ Implement the treatment-specific regulatory pattern analysis
  - Implemented in `scripts/new_regulatory_analysis/treatment_regulatory_analysis.py`
  - Fixed variant calling provides complete dataset with all treatment samples
  - Successfully identified treatment-specific regulatory patterns
  - Generated comparative analysis of temperature vs. low oxygen adaptation
  - Key findings: Core promoter regions (65.9%) are preferentially affected across all treatments
  
- ⏳ Develop the ERG-OSH regulatory relationship analysis
  - Next component to be implemented
  - OSH gene analysis completed, providing foundation for regulatory relationship study
  
- ⏳ Create the four-zone integration analysis
  - To be implemented after ERG-OSH relationship analysis
  
- ⏳ Build the visualization and reporting system
  - Basic reporting implemented
  - Advanced visualizations to be implemented

### Phase 3: Integration and Validation - 🔜 PENDING
- ⏳ Integrate all components into a cohesive analytical pipeline
- ⏳ Validate findings against existing knowledge
- ⏳ Perform sensitivity analyses to ensure robustness
- ⏳ Generate comprehensive reports and visualizations

## Expected Outcomes and Progress

1. **Comprehensive Regulatory Annotation**: ✅ COMPLETED
   - Successfully mapped 44 variants to specific regulatory features
   - Classified variants into core promoter (65.9%) and UAS proximal (34.1%) regions
   - Established mean distance of 180.4 bp from genes with maximum of 332 bp

2. **Treatment-Specific Regulatory Signatures**: ✅ COMPLETED
   - Identified consistent pattern across treatments (~66.7% core promoter variants)
   - Generated comparative analysis showing minimal differences (2.0%) between temperature and low oxygen adaptation
   - All treatments show similar regulatory region usage despite different selective pressures

3. **ERG-OSH Regulatory Network**: 🔄 IN PROGRESS
   - OSH gene family analysis completed, providing foundation for regulatory relationship study
   - Identified close physical proximity between certain OSH-ERG gene pairs
   - Further relationship analysis to be implemented in next phase

4. **Regulatory Impact Prioritization**: ⏳ PENDING
   - Framework for impact scoring to be developed

5. **Four-Zone Regulatory Model**: 🔄 IN PROGRESS
   - Confirmed all regulatory variants are located in buffer zone upstream of genes
   - Further integration with conservation architecture to be developed

6. **Statistical Validation**: 🔄 IN PROGRESS
   - Basic statistical analysis completed for treatment comparisons
   - More extensive validation to be implemented

7. **Biological Interpretation**: 🔄 IN PROGRESS
   - Initial interpretation completed, showing adaptation through regulatory changes
   - Core promoter enrichment suggests adaptation through changes in basal transcription machinery recruitment
   - Further biological interpretation to be developed as additional components are implemented

## Biological Significance

This improved regulatory analysis will provide deeper insights into how yeast adapts to environmental stresses through regulatory changes rather than alterations to protein structure. By connecting regulatory changes to the four-zone conservation architecture, we can understand how organisms maintain essential functions while enabling adaptation through carefully controlled regulatory mechanisms. The analysis of ERG-OSH regulatory relationships will further illuminate the coordination between sterol synthesis and transport in membrane adaptation.

This approach moves beyond simple variant counting to a biologically informed analysis of regulatory mechanisms, providing a more complete understanding of adaptation through gene expression changes rather than protein alterations. The integration with our existing findings will create a comprehensive model of yeast adaptation that connects genomic changes to regulatory mechanisms and ultimately to phenotypic outcomes.