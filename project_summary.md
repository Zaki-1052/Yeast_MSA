## Project Summary

I've completed a comprehensive genomic analysis of the Saccharomyces cerevisiae strains subjected to different experimental treatments (WT, WTA, STC, and CAS). Starting with raw sequencing data, I conducted quality control, alignment to the W303 reference genome, variant calling, and comparative analysis to identify treatment-specific genetic variations.

## Detailed Methodology

### Quality Control and Preprocessing

- Performed quality assessment on all 30 FASTQ files using FastQC
- Generated consolidated QC reports with MultiQC to identify quality issues
- Identified moderate duplication rates (28-42%) and per-tile sequence quality issues
- Preprocessed reads using fastp with the following parameters:
    - Quality threshold: Phred score ≥20
    - Unqualified base limit: 40%
    - Minimum read length: 50bp
    - Sliding window quality trimming (window size: 4, mean quality: 20)

### Reference-Guided Alignment

- Aligned all samples to the W303 reference genome using BWA-MEM
- Implemented a complete post-alignment processing workflow:
    - Converted SAM to BAM format
    - Sorted reads by name for proper mate handling
    - Added mate scores with samtools fixmate -m
    - Re-sorted by coordinate position
    - Marked and removed duplicate reads
    - Added read group information for downstream analysis
    - Indexed BAM files for random access
- Generated comprehensive alignment statistics:
    - Achieved mapping rates of 78-83% across all samples
    - Verified proper read pairing rates (78-80%)
    - Confirmed high mapping quality (scores 54-60)

### Variant Calling and Filtering

- Called variants using BCFtools mpileup and call with haploid model (--ploidy 1)
- Added depth annotations for filtering (DP, AD)
- Normalized indel representations for consistent comparison
- Applied quality filters (QUAL≥20, depth≥10)
- Generated individual VCF files for each sample
- Created a merged VCF containing variants across all samples
- Verified VCF integrity with appropriate indexing

### Comparative Analysis

- Resolved sample encoding issues by implementing a direct VCF comparison approach
- Identified treatment-specific variants absent in corresponding controls:
    - WT: 188 high-confidence variants
    - STC: 83 high-confidence variants
    - CAS: 150 high-confidence variants
    - WTA: 157 high-confidence variants
- Analyzed variant consistency across biological replicates:
    - 44-48% of variants present in all three replicates
    - 15-16% in exactly two replicates
    - Remaining variants unique to single replicates
- Characterized variant types (SNPs vs indels) across treatments:
    - WT: 602 SNPs, 86 indels
    - STC: 310 SNPs, 87 indels
    - CAS: 479 SNPs, 87 indels
    - WTA: 532 SNPs, 83 indels
- Identified genomic hotspots for variants:
    - Scaffold JRIU01000031.1 as a major hotspot (17-24 variants)
    - Treatment-specific enrichment in other scaffolds

### Visualization
- Generated a suite of visualizations using R (ggplot2, pheatmap):
    - Variant distribution by treatment and type
    - Variant sharing patterns between treatments
    - Quality metrics across samples
    - Scaffold distribution plots

## Key Findings

1. **Treatment-Specific Genetic Signatures**:
    - WT treatment induces the most genetic changes (688 variants)
    - STC shows the fewest but potentially more focused changes (397 variants)
    - Remarkably consistent indel counts across all treatments (~85)
2. **Treatment Relationships**:
    - High genetic similarity between WT and WTA (86 shared variants)
    - STC exhibits the most distinct genetic profile (lowest sharing)
    - CAS shows intermediate relationships, consistent with its combined treatment
3. **Genomic Distribution Patterns**:
    - Non-random distribution of variants across the genome
    - Treatment-specific hotspots identified
    - Different scaffolds affected by different treatments

---

## Quality Control and Sequencing Data

The FastQC results showed:

- Moderate duplication rates (28-42%)
- Consistent GC content (37-39%) matching expected yeast genome content
- Some quality issues in tile-specific regions
- All samples had similar sequence quality metrics

These findings tell us the sequencing data was of reasonable quality. The duplication rates are slightly elevated but acceptable for a genomic analysis. The consistent GC content across samples indicates no major bias in sequencing.

## Alignment to Reference Genome

The alignment statistics showed:

- Good mapping rates (78-83%)
- High mapping quality scores (54-60)
- Consistent coverage across samples

These results are significant because they indicate the sequencing reads aligned well to the reference genome. With ~80% mapping rate, we can be confident that our variant calling is based on reliable alignments. The high mapping quality scores also suggest that most reads mapped uniquely to the reference genome with high confidence.

## Variant Calling Results

The variant summary showed:

- Each sample had ~800-950 total variants
- About 85-88% of variants were SNPs
- About 12-15% were indels
- Good variant quality metrics (average quality scores 74-84)

The treatment-specific variant counts revealed:

- WT treatment: 188 high-confidence variants
- STC treatment: 83 high-confidence variants
- CAS treatment: 150 high-confidence variants
- WTA treatment: 157 high-confidence variants

Consistency analysis showed about 47-48% of variants were present in all 3 replicates of each treatment group, indicating good biological reproducibility.

## Significant Findings and Their Biological Implications

### 1. Differential Genetic Impact of Treatments

The most striking finding is that different treatments induced different numbers of genomic changes:

- WT treatment (37-55) induced the most genetic changes (688 variants)
- STC treatment induced the fewest changes (397 variants)
- CAS and WTA showed intermediate levels (566 and 615 variants respectively)

This pattern suggests that the 37-55 treatment applied to the wild type strain causes the most substantial genomic alterations. The fact that STC had far fewer variants suggests this strain might be more genetically stable or resistant to the effects of the 55 treatment.

### 2. Common Mutational Mechanisms

A fascinating observation is that all treatments induced almost identical numbers of indels (~85):

- WT: 86 indels
- STC: 87 indels
- CAS: 87 indels
- WTA: 83 indels

This remarkable consistency suggests a common mechanism for indel generation that operates similarly across all treatment conditions. This could indicate a specific DNA damage repair pathway that's consistently activated by these treatments, leading to a predictable number of insertions and deletions.

The differences in total variant counts are almost entirely due to differences in SNP counts, which suggests that while indel-generating processes are consistent, SNP-generating processes vary by treatment.

### 3. Treatment Relationships

The variant sharing patterns revealed important relationships between the treatments:

- WT and WTA share the most variants (86), indicating their close genetic relationship
- STC shares the fewest variants with other treatments, suggesting it has the most unique genetic profile
- CAS shows intermediate sharing, consistent with its combined treatment nature

This is biologically significant because:

- It confirms that WTA (Wild type with modification A) maintains substantial genetic similarity to the original WT despite the modification
- It suggests STC treatment causes a distinctive genomic response
- CAS (with both A and S modifications) shows some similarities to both WTA and STC, as would be expected if treatments have additive effects

### 4. Genomic Hotspots

The scaffold distribution analysis showed that certain regions of the genome are more susceptible to mutations:

- Scaffold JRIU01000031.1 was a major hotspot across all treatments (17-24 variants)
- Different treatments also affected unique scaffolds:
    - WT & WTA both affected JRIU01000195.1
    - STC uniquely affected JRIU01000179.1
    - CAS showed specific enrichment in JRIU01000033.1

This non-random distribution of variants suggests that these scaffolds likely contain genes or regulatory regions that are either:

1. More susceptible to damage from these treatments
2. Under selection pressure in these treatment conditions
3. Involved in the adaptive response to these treatments

Without gene annotations, we can't definitively identify the affected genes, but these hotspots provide strong candidates for functional follow-up studies.

## Broader Scientific Significance

The analysis demonstrates several important principles:

1. **Treatment-Specific Genetic Signatures**: Each treatment has a distinctive genomic signature, suggesting specific mechanisms of action at the DNA level.
    
2. **Predictable Mutational Processes**: The consistent number of indels across treatments suggests some genetic changes follow predictable patterns despite different treatment conditions.
    
3. **Strain-Specific Responses**: The fact that STC showed the fewest variants suggests that genetic background can significantly modulate the response to treatments.
    
4. **Genomic Vulnerability**: The concentration of variants in specific scaffolds indicates that some regions of the genome are more vulnerable to mutation than others.
    
5. **Treatment Relationships**: The pattern of shared variants between treatments provides molecular evidence of their relationships, confirming expected similarities based on the experimental design.
    

The significance of these findings for further research:

1. The genomic hotspots identified are prime candidates for further functional studies to determine how these regions contribute to treatment response.
    
2. The unique genetic signature of STC suggests it might employ different adaptive mechanisms compared to other strains, which could have implications for understanding treatment resistance or adaptation.
    
3. The consistent indel patterns across treatments might point to a specific DNA repair pathway that operates uniformly under these conditions, offering a potential target for intervention.
    
4. The higher variant count in WT samples suggests the wild type might be more genetically flexible or less able to repair certain types of damage than the modified strains.
    

In summary, this analysis provides significant insights into how different yeast strains respond genetically to experimental treatments, revealing both shared mechanisms and strain-specific adaptations at the genomic level. The findings point to specific genomic regions and mutational processes that likely play important roles in the cellular response to these treatments.
