# summary_document.R

# Load required libraries
library(rmarkdown)
library(knitr)
library(pander)

# Create directory for the summary document
dir.create("results/summary", recursive = TRUE, showWarnings = FALSE)

# Create R Markdown content as a string
rmd_content <- '
---
title: "Yeast Genome Analysis: Summary of Findings"
author: "Bioinformatics Team"
date: "`r format(Sys.time(), \"%B %d, %Y\")`"
output: 
  pdf_document:
    toc: true
    toc_depth: 2
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(knitr)
library(dplyr)
library(ggplot2)
library(kableExtra)
```

# Executive Summary

This report summarizes the genetic effects of four experimental treatments (WT, WTA, STC, and CAS) on the genome of *Saccharomyces cerevisiae* (baker\'s yeast). Through whole-genome sequencing and variant analysis, we have identified distinct mutation patterns that provide insight into how each treatment affects DNA integrity.

**Key Findings:**

1. **Each treatment creates a unique "genetic signature"** with characteristic mutation patterns
2. **Mutations are highly clustered** (71-84%) rather than randomly distributed throughout the genome
3. **Treatment effects target specific genomic regions**, particularly scaffold JRIU01000031.1
4. **CAS treatment shows a remarkably specific mutation mechanism** (100% transition mutations)
5. **WT and WTA treatments show complementary mutation patterns** (G→A vs. A→G) suggesting related but opposite effects

## What This Means

These findings suggest the treatments are affecting specific DNA repair or replication mechanisms rather than causing random damage. The high degree of clustering indicates that certain genomic regions are particularly vulnerable to the effects of these treatments.

# Understanding Mutation Patterns

## Mutation Types Explained

For non-specialists, DNA mutations can be classified into two main categories:

**Transitions:** Changes between similar bases:
  - A → G or G → A (purine to purine)
  - C → T or T → C (pyrimidine to pyrimidine)
  
**Transversions:** Changes between dissimilar bases:
  - A → C, A → T, G → C, or G → T (purine to pyrimidine)
  - C → A, C → G, T → A, or T → G (pyrimidine to purine)

Transitions typically occur more frequently in nature and often result from specific biological processes like deamination. Transversions generally require more substantial chemical changes to the DNA.

```{r transition-transversion, fig.width=8, fig.height=5, fig.cap="Transitions vs. transversions across treatments. CAS shows exclusively transition mutations."}
# Import the transition vs. transversion data
transition_data <- read.table("results/visualizations/signatures/transition_transversion.tsv", header=TRUE)

# Create a simple bar chart
ggplot(transition_data, aes(x=Treatment, y=Percentage, fill=Mutation_Class)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(Percentage, 1), "%")), 
            position=position_stack(vjust=0.5), color="white", fontface="bold") +
  scale_fill_manual(values=c("Transition"="#4daf4a", "Transversion"="#e41a1c")) +
  labs(title="Mutation Types by Treatment",
       y="Percentage of Variants (%)",
       fill="Mutation Type") +
  theme_minimal() +
  theme(legend.position="bottom")
```

## Mutation Clustering

One of our most striking findings is that mutations are not randomly distributed across the genome. Instead, they occur in clusters, with 71-84% of mutations located within 10 base pairs of another mutation.

```{r clustering, fig.width=8, fig.height=4, fig.cap="Percentage of variants found in clusters across treatments."}
# Import clustering data
clustering_data <- read.table("results/visualizations/clustering/clustering_summary.tsv", header=TRUE)

# Create a bar chart of clustering percentages
ggplot(clustering_data, aes(x=Treatment, y=Percent_Clustered, fill=Treatment)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(Percent_Clustered, "%")), vjust=-0.5) +
  ylim(0, 100) +
  labs(title="Mutation Clustering by Treatment",
       y="Percentage of Mutations in Clusters (<10bp)") +
  theme_minimal() +
  theme(legend.position="none")
```

# Treatment-Specific Effects

Each treatment produces a characteristic mutation pattern that provides clues about its biological effects:

```{r treatment-summary, echo=FALSE}
# Create a summary table of treatment characteristics
summary_data <- read.table("results/visualizations/integrative/treatment_summary.tsv", 
                          header=TRUE)

# Format the table for display
summary_table <- summary_data %>%
  select(Treatment, Percent_Clustered, TransitionPct, DominantMutation, Percentage, TopHotspot) %>%
  rename("% Clustered" = Percent_Clustered,
         "% Transitions" = TransitionPct,
         "Top Mutation" = DominantMutation,
         "% of Top Mutation" = Percentage,
         "Hotspot Scaffold" = TopHotspot)

knitr::kable(summary_table, 
            caption = "Summary of treatment-specific effects",
            booktabs = TRUE,
            align = "lccccc") %>%
  kable_styling(latex_options = c("striped", "hold_position"),
                full_width = FALSE)
```

## CAS Treatment

The CAS treatment shows a remarkably specific effect on the yeast genome:

- **100% transition mutations** — a highly unusual pattern indicating a very specific DNA repair defect
- **High clustering** (83.6%) — mutations are tightly grouped
- **C→T bias** (40% of mutations) — consistent with cytosine deamination or defective repair of this specific mutation type

These characteristics suggest CAS might specifically affect pathways that repair cytosine deamination damage, such as components of the base excision repair system.

## WT Treatment

The WT treatment shows:

- **Strong G→A bias** (45.5% of mutations)
- **High clustering** (83.1%)
- Targets a different hotspot (JRIU01000157.1) than other treatments

This pattern is typically associated with oxidative damage to guanine bases, which can form 8-oxoguanine that pairs with adenine instead of cytosine during DNA replication.

## WTA Treatment

Interestingly, WTA treatment shows:

- **Strong A→G bias** (55.6%) — the complement of WT\'s dominant mutation pattern
- **High clustering** (81.5%)
- **Similar transition percentage** to WT (77.8% vs. 77.3%)

The complementary mutation pattern suggests that WTA\'s "modification A" might affect similar biological pathways as WT but with an effect on the opposite DNA strand.

## STC Treatment

STC treatment shows the most diverse mutation profile:

- **Lowest clustering** (70.6%)
- **Most balanced mutation spectrum** with highest transversion rate (40%)
- No strongly dominant mutation type

This suggests STC may disrupt multiple DNA repair pathways simultaneously or cause more generalized genomic instability.

# Hotspot Analysis

All treatments show particular vulnerability in specific genomic regions, especially scaffold JRIU01000031.1.

![Position of variants in JRIU01000031.1 across treatments](../visualizations/clustering/position_distribution.png)

Reasons why certain regions might be more susceptible to mutation include:

1. **DNA secondary structures** that challenge replication machinery
2. **Sequence context** such as repeats or homopolymer runs
3. **Late replication timing** which can reduce repair efficiency
4. **Functional importance** of the region to cellular processes

# Treatment Relationships

Our analysis revealed unexpected relationships between treatments:

- **WT and CAS are most similar** (0.82 similarity) despite having different dominant mutations
- **STC and WTA show strong similarity** (0.77)
- **CAS and WTA show lowest similarity** (0.65)

![Treatment relationship network](../visualizations/integrative/treatment_network.png)

# Biological Implications

Based on our findings, we can make several inferences about how these treatments affect yeast biology:

1. **Specific repair pathways are targeted** rather than causing random DNA damage
2. **CAS likely affects mismatch repair** specifically for cytosine deamination events
3. **WT and WTA appear to target the same pathway** but with strand-specific effects
4. **STC creates more generalized genomic instability** affecting multiple repair systems

# Recommendations for Follow-up Studies

To build on these findings, we recommend:

1. **Sequencing of specific DNA repair genes** in treated yeast to identify mutations in repair machinery
2. **Functional assays** to measure the activity of specific DNA repair pathways in treated cells
3. **RNA-seq analysis** to determine if expression of repair genes is altered by treatments
4. **Targeted sequencing of hotspot regions** in a larger sample to validate mutation patterns
5. **Phenotypic correlation studies** to connect these genetic changes with observable traits

# Conclusion

Our analysis has revealed that each treatment produces a distinct genetic signature with specific mutation patterns. The high degree of clustering and the predominance of transitions suggest these treatments are affecting specific DNA repair or replication mechanisms rather than causing random damage.

Most notably, the CAS treatment shows a remarkably specific effect (100% transition mutations), while the complementary patterns in WT and WTA (G→A vs. A→G) suggest related but opposite effects. The STC treatment appears to cause more generalized genomic instability.

These findings provide valuable insights into the molecular mechanisms of these treatments and suggest specific pathways for further investigation.

'

# Write the R Markdown content to a file
writeLines(rmd_content, "results/summary/yeast_analysis_summary.Rmd")

# Render the R Markdown document to PDF
render("results/summary/yeast_analysis_summary.Rmd", 
       output_format = "pdf_document", 
       output_file = "Yeast_Genome_Analysis_Summary.pdf",
       output_dir = "results/summary/")

cat("Non-technical summary document created successfully at: results/summary/Yeast_Genome_Analysis_Summary.pdf\n")