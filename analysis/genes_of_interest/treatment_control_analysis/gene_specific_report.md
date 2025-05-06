# Gene-Specific Analysis of Treatment vs Control Data

## Overview

This report provides an analysis of gene-specific variation patterns across treatments,
with a focus on detecting purifying selection in ergosterol pathway genes.

## Treatment Summary

| Treatment | Total Variants | ERG Genes | Non-ERG Genes | Non-Genic | Overall Fold Change |
| --- | --- | --- | --- | --- | --- |
| WT-37 | 27 | 0 | 0 | 27 | 0.07 |
| WTA | 20 | 0 | 0 | 20 | 0.05 |
| STC | 27 | 0 | 0 | 27 | 0.07 |
| CAS | 38 | 0 | 0 | 38 | 0.10 |
| STC-vs-STCCTRL | 27 | 0 | 0 | 27 | 0.07 |
| CAS-vs-CASCTRL | 38 | 0 | 0 | 38 | 0.10 |

## Purifying Selection Analysis

Purifying selection can be detected by comparing the observed number of variants in a gene
to the expected number based on gene length and genome-wide mutation rate.
A fold change < 1 indicates fewer variants than expected (purifying selection),
while a fold change > 1 indicates more variants than expected (positive selection).
