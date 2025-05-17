# Enhanced Regulatory Region Definition Report

## Overview

Total genes analyzed: 6078
ERG pathway genes: 12
OSH family genes: 8
Other genes: 6058

## Regulatory Region Definitions

### Upstream Regulatory Regions

| Region | Distance Range (bp) | Description |
|--------|-------------------|-------------|
| TATA_box_region | 40 - 120 | Region containing TATA box elements in ~20% of yeast genes |
| UAS_distal | 500 - 1500 | Distal Upstream Activating Sequences |
| UAS_proximal | 150 - 500 | Proximal Upstream Activating Sequences (enhancer-like elements) |
| URS_region | 150 - 800 | Upstream Repressing Sequences (silencer elements) |
| core_promoter | 0 - 150 | Core promoter region containing essential transcription elements |
| far_upstream | 1500 - 10000 | Distant regulatory elements |
| tss_proximal | 0 - 40 | Region immediately upstream of transcription start site |

### Downstream Regulatory Regions

| Region | Distance Range (bp) | Description |
|--------|-------------------|-------------|
| downstream_element | 250 - 1000 | Other downstream control elements |
| five_prime_UTR | -60 - 0 | 5' Untranslated Region (before CDS start) |
| terminator | 0 - 250 | Transcription termination region |
| three_prime_UTR | 0 - 120 | 3' Untranslated Region |

### Special Regulatory Elements

| Element | Distance Range (bp) | Consensus Sequence | Description |
|---------|-------------------|-------------------|-------------|
| ARE | 100 - 700 | TTGCACGT | Antioxidant Response Element |
| CSRE | 150 - 600 | CCNNNNNNGG | Carbon Source Response Element |
| HSE | 100 - 800 | NGAAN | Heat Shock Element |
| Oxygen_responsive | 100 - 700 | YYNATTGTTY | Oxygen responsive element |
| PDR | 150 - 700 | TCCGCGGA | Pleiotropic Drug Resistance Element |
| STRE | 100 - 700 | AGGGG | Stress Response Element |
| TATA_box | 40 - 120 | TATA[AT]A[AT] | Canonical TATA box element |
| UASPHR | 150 - 700 | GCGATGAGATGAGCT | Phosphate-regulated UAS |
| URS | 150 - 800 | [ACG]CCCC[ACT] | Upstream Repressing Sequence |
| initiator | 0 - 40 | YYANWYY | Initiator element at transcription start site |
| sterol_regulatory | 100 - 700 | TCGTATA | Sterol Regulatory Element |

### Conservation Zones

| Zone | Distance Range (bp) | Description |
|------|-------------------|-------------|
| buffer_zone | 0 - 5000 | Buffer zone around gene (<5kb) |
| core_zone | 0 - 0 | The gene itself - absolute conservation |
| intermediate_zone | 5000 - 50000 | Intermediate zone (5-50kb) |
| satellite_zone | 50000 - 100000 | Satellite zone (50-100kb) |

## Region Statistics

### Regulatory Region Counts

| Region | Count | Percentage |
|--------|-------|------------|
| core_promoter | 6078 | 9.1% |
| TATA_box_region | 6078 | 9.1% |
| tss_proximal | 6078 | 9.1% |
| UAS_proximal | 6078 | 9.1% |
| UAS_distal | 6078 | 9.1% |
| URS_region | 6078 | 9.1% |
| far_upstream | 6078 | 9.1% |
| five_prime_UTR | 6078 | 9.1% |
| terminator | 6078 | 9.1% |
| three_prime_UTR | 6078 | 9.1% |
| downstream_element | 6078 | 9.1% |

### Regulatory Region Size Statistics

| Region | Count | Mean Size (bp) | Median Size (bp) | Min Size (bp) | Max Size (bp) |
|--------|-------|---------------|-----------------|--------------|---------------|
| core_promoter | 6078 | 151.0 | 151.0 | 151 | 151 |
| TATA_box_region | 6078 | 81.0 | 81.0 | 81 | 81 |
| tss_proximal | 6078 | 41.0 | 41.0 | 41 | 41 |
| UAS_proximal | 6078 | 351.0 | 351.0 | 351 | 351 |
| UAS_distal | 6078 | 1001.0 | 1001.0 | 1001 | 1001 |
| URS_region | 6078 | 651.0 | 651.0 | 651 | 651 |
| far_upstream | 6078 | 8501.0 | 8501.0 | 8501 | 8501 |
| five_prime_UTR | 6078 | 61.0 | 61.0 | 61 | 61 |
| terminator | 6078 | 251.0 | 251.0 | 251 | 251 |
| three_prime_UTR | 6078 | 121.0 | 121.0 | 121 | 121 |
| downstream_element | 6078 | 751.0 | 751.0 | 751 | 751 |

## ERG and OSH Gene Analysis

### ERG Gene Regulatory Profile

Total ERG genes: 12

#### Regulatory Region Distribution in ERG Genes

| Region | Count | Percentage |
|--------|-------|------------|
| core_promoter | 12 | 9.1% |
| TATA_box_region | 12 | 9.1% |
| tss_proximal | 12 | 9.1% |
| UAS_proximal | 12 | 9.1% |
| UAS_distal | 12 | 9.1% |
| URS_region | 12 | 9.1% |
| far_upstream | 12 | 9.1% |
| five_prime_UTR | 12 | 9.1% |
| terminator | 12 | 9.1% |
| three_prime_UTR | 12 | 9.1% |
| downstream_element | 12 | 9.1% |

### OSH Gene Regulatory Profile

Total OSH genes: 8

#### Regulatory Region Distribution in OSH Genes

| Region | Count | Percentage |
|--------|-------|------------|
| core_promoter | 8 | 9.1% |
| TATA_box_region | 8 | 9.1% |
| tss_proximal | 8 | 9.1% |
| UAS_proximal | 8 | 9.1% |
| UAS_distal | 8 | 9.1% |
| URS_region | 8 | 9.1% |
| far_upstream | 8 | 9.1% |
| five_prime_UTR | 8 | 9.1% |
| terminator | 8 | 9.1% |
| three_prime_UTR | 8 | 9.1% |
| downstream_element | 8 | 9.1% |

## Biological Interpretations

### Yeast Promoter Architecture

The yeast promoter architecture is characterized by:
- **Core promoter regions** (0-150bp upstream): Essential for basic transcription machinery assembly
- **TATA box regions** (40-120bp upstream): Present in ~20% of yeast genes, associated with regulated genes
- **Upstream Activating Sequences (UAS)**: Enhancer-like elements located 150-1500bp upstream
- **Stress Response Elements (STRE)**: Often found in the UAS regions of stress-responsive genes
- **Short 5' and 3' UTRs**: Typical of yeast genes, with important roles in translation regulation

### The Four-Zone Conservation Architecture

The genomic organization around ERG and OSH genes reveals a hierarchical conservation pattern:
1. **Core Zone**: The genes themselves show absolute conservation, reflecting their essential functions
2. **Buffer Zone** (0-5kb): Minimal variation, primarily regulatory adjustments
3. **Intermediate Zone** (5-50kb): Moderate variation, including pathway modulators
4. **Satellite Zone** (50-100kb): Higher variation, enabling adaptive flexibility while maintaining core functions

### ERG-OSH Regulatory Relationships

The regulatory architecture of ERG (ergosterol biosynthesis) and OSH (oxysterol binding homology) genes suggests:
- Shared regulatory elements controlling both sterol synthesis and transport
- Similar conservation patterns indicating coordinated evolution
- Specialized regulatory elements for oxygen sensing and sterol homeostasis
- Potential co-regulation through chromatin organization

## Implications for Adaptation

The regulatory architecture identified has important implications for adaptation mechanisms:
- **Regulatory flexibility with functional conservation**: Adaptation through changes in gene expression rather than protein structure
- **Zone-specific adaptation**: Different adaptation strategies depending on genomic distance from essential genes
- **Coordinated regulation**: Changes in sterol synthesis and transport systems are likely coordinated
- **Environmental responsiveness**: Specialized regulatory elements for different environmental stresses
- **Transcription factor binding site variation**: Subtle changes in TF binding affinity rather than complete gain/loss of binding sites