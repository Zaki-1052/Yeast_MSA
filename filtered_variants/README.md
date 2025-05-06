# Filtered Variants for Gene-Mapped Chromosomes

## Background and Issue

During the implementation of gene-specific functionality in the mutation analysis scripts, we discovered that variants could not be mapped to genes due to a mismatch between chromosome IDs. After detailed investigation:

1. The reference gene mapping file (`reference/gene_mapping.tsv`) only contains:
   - 12 genes (+ 1 header line)
   - Spread across 5 chromosomes: CM007970.1, CM007971.1, CM007975.1, CM007976.1, CM007977.1
   - These genes are the ergosterol pathway genes of interest

2. The variant data contains:
   - 18 different chromosomes in the VCF files
   - About 330-370 variants per treatment/control sample
   - Many of these variants are on chromosomes not represented in the gene mapping file

## Solution Approach

To enable gene-specific analysis with the limited gene mapping data available, we:

1. Created a filtering script (`filter_variation.py`) that:
   - Extracts variants from all VCF files
   - Keeps only variants on the 5 chromosomes with gene mapping data
   - Saves filtered variant data to CSV files (this directory)

2. Filtered variant findings:
   - 350-374 total variants per treatment sample, of which 134-151 are on gene-mapped chromosomes (~40%)
   - 341 total variants in the control sample, of which 130 are on gene-mapped chromosomes (~38%)
   - Across all samples, chromosome distribution of variants is similar

## Next Steps

Using this filtered data, we can proceed with two parallel approaches:

1. **Limited Gene Analysis**:
   - Modify the gene analysis scripts to use these filtered variant files
   - Focus analysis specifically on the 5 chromosomes with available gene mapping
   - Generate meaningful results for the subset of variants that can be analyzed

2. **Complete Gene Mapping Dataset**:
   - Obtain comprehensive gene mapping data by parsing the annotated genbank files
   - Create a complete gene mapping file that covers all 18 chromosomes
   - Re-run analysis once a complete mapping is available

## Statistics by Treatment

| Treatment  | Total Variants | Filtered Variants | % Retained |
|------------|----------------|-------------------|------------|
| WT-37-55-1 | 350            | 139               | 39.7%      |
| WT-37-55-2 | 372            | 137               | 36.8%      |
| WT-37-55-3 | 367            | 140               | 38.1%      |
| WTA-55-1   | 353            | 139               | 39.4%      |
| WTA-55-2   | 349            | 134               | 38.4%      |
| WTA-55-3   | 360            | 139               | 38.6%      |
| STC-55-1   | 365            | 137               | 37.5%      |
| STC-55-2   | 364            | 141               | 38.7%      |
| STC-55-3   | 365            | 144               | 39.5%      |
| CAS-55-1   | 370            | 142               | 38.4%      |
| CAS-55-2   | 368            | 151               | 41.0%      |
| CAS-55-3   | 374            | 144               | 38.5%      |
| WT-CTRL    | 341            | 130               | 38.1%      |

## Chromosomes with Gene Mapping

| Chromosome | Variants (average) | Notes |
|------------|--------------------|----------------------------------------------------|
| CM007970.1 | 33-37              | Contains ERG1, ERG4, ERG25 genes                   |
| CM007971.1 | 14-17              | Contains ERG11 gene                                |
| CM007975.1 | 34-39              | Contains multiple ergosterol pathway genes         |
| CM007976.1 | 30-33              | Contains ergosterol pathway genes                  |
| CM007977.1 | 18-29              | Contains multiple ergosterol pathway genes         |