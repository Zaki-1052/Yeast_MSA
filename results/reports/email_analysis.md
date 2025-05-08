

Wanted to send the few folders of results of the annotated variants in a separate email before the functional impact and sterol regression; you can see each of the statistics of variant analysis from snpEff [here in the stats folder](https://github.com/Zaki-1052/Yeast_MSA/tree/main/vcf/annotated/stats) -- if you aren't familiar, it's short for SNP (like single nucleotide polymorphism) Effect, but obviously it does every type of mutation. And it cross references the list of variants from my [original analysis report here](https://github.com/Zaki-1052/Yeast_MSA/blob/main/vcf/analysis_report.txt) with the newest reference and its associated genome annotations from WGAP. It then predicts the location, effect (precursor to functional impact), and the gene/transcript info with what proteins would be affected. For example, it will annotate each variant like, "This SNP causes a missense mutation in the ERG11 gene, changing an alanine to a valine at position 123," or whatever.

Also, it might be worth mentioning that all of the [preliminary non-annotated VCF analysis](https://github.com/Zaki-1052/Yeast_MSA/tree/main/analysis) still applies, whether or not on the specific genes of interest, as well as the earlier main analysis, scaffold, and annotation reports as HTML files I sent on 4/20; I can always redo them with the mapping to various non-specific genes, now that I have the annotations and don't need to define per-scaffold. Now, I have the main analyses that resulted from those variants on the genes of interest that you listed on 4/13; [you can find that folder with expanded scope here](https://github.com/Zaki-1052/Yeast_MSA/tree/main/results/gene_variants_expanded). The two nested folders inside those results also list the variants specific to each treatment and sorted by each gene in TSV format. Now, there are a few possibilities for why the results for these mutations seemed lacking compared to previous analysis:

  

One is that the distance from the genes of interest is simply too great, and the initial scope of the variant analysis had a higher false positive rate without accounting for how many BPs away the mutations were away from the specific genes. In my opinion, the methodology was perfectly sound, but of course it's worth mentioning that I also could have made a mistake somewhere. The other main possibility is that with mutations mainly found on YGR175C (ERG1) (29 - 50%), and the other half split between YGR060W (15 - ERG25) and YMR202W (15 - ERG2) that purifying selection is taking place on the ergosterol pathways, and measuring enrichment in positive terms is the wrong approach biologically to take. The last possibilities are that the genes of interest listed were represented by an incorrect model (the other 8 genes are insignificant), or either the data from the sequences given is polluted, or there is significant enrichment, just not on the specified genes of interest.

  

Now, in terms of taking the results at face-value, I also performed a scaffold-level analysis while still accounting for the annotated positions/distance, and you can see [the results of the scaffold variants here](https://github.com/Zaki-1052/Yeast_MSA/tree/main/results/scaffold_variants). With this approach, we found that the twelfth scaffold (ERG12) had the most variants, but confirmed that almost 3/4s of the variants are over 50k base pairs away from our target genes. This is what I meant by purifying selection acting on the ergosterol pathway genes, making them relatively resistant to genetic variation, which might make sense biologically once linked to the sterols -- we'd be looking at a negative correlation then. 

  

Some other specific patterns are that YGR060W (ERG25) had the most variants in close proximity (60 within 5kb), YHR007C (ERG11) had 31 variants within 5kb, all upstream, YGR175C (ERG1) -- 29 variants within 5kb, all downstream 1kb -- and YMR202W (ERG2) and YLR056W (ERG3) had 15 variants each within 5kb. Lastly, 7 genes had zero variants within 5kb: ERG4, ERG7, ERG9, ERG6, ERG5, ERG24 (we had mentioned ERG24 in the preliminary as having a statistically significant lesser number of variants compared to the others; or that there was a low Q probability that it would have that few number of variants given normal variation). 

  

Some other interesting results were that there were far more insertions (50%) than single-nucleotide variants (SNVs and deletions were 25%), and that 80% of the upstream gene variants were primarily regulatory (modifier) region effects, and that there were relatively few coding variants; [you can see the exact numbers here](https://github.com/Zaki-1052/Yeast_MSA/blob/main/results/scaffold_variants/scaffold_variant_summary.txt), but basically very few affecting protein structure directly. There was also an extremely similar distribution among treatment groups, though if you look at the group comparisons, control samples obviously showed far fewer variants (by a magnitude of roughly 4). But generally, this confirms that the pipeline was sound, and there were no dramatic differences in impact distribution across treatments (when focused on the ERG pathway genes). 

  

Biologically, I think that this lack of variants within or near most ergosterol genes would mean that there's some selection for conservation of the ergosterol pathway -- which, I'm sorry, from a quick skim, I just realized that was your hypothesis from the PowerPoint, wasn't it? I don't really want to rewrite the email, but you can probably disregard what I said earlier about there being possible error/getting "bad" results. In that case, all of the zeros are probably significant if it's found that all the mutations were concentrated in regions far from the ERG genes of interest. Even the 80% upstream gene variants could be regulatory changes rather than direct protein modifications, which I can confirm with functional impact analysis and the sterols!

  

Some other things to note are that the predominance of upstream variants definitely suggests that when changes do occur, they're more likely to affect gene regulation than protein structure. Also, ERG25 is probably worth examining more closely given that there are so many more scaffold variants (60) -- it might be more tolerant to nearby genetic variation (though if I narrow the scope to my initial approach, we see 0 variants, [found in this file here](https://github.com/Zaki-1052/Yeast_MSA/blob/main/results/gene_variants/gene_summary.tsv)). Importantly, we found 2,087 variants across all scaffolds, but only 44 in the gene-focused analysis, so the coordinates were correct, but variants simply aren't common in these specific genes. Also, the distribution of variants by treatment is proportional to the number of samples in each group, suggesting no strong treatment-specific effect.

  

We can also confirm the upstream variants in ERG1 and ERG2 we found in the initial analysis, and again, the fact that most genes of interest in the ERG pathway have so few variants within 5kb supports the hypothesis that these ergosterol pathway genes are under purifying/negative selection, which supports the idea that deleterious mutations are selectively removed from the population, maintaining the integrity of essential genes and functions. Not sure if this is already supported by existing research, but when I looked this up, I did see that core enzymes are highly conserved across different yeast strains, as mutations disrupting these enzymes would likely be detrimental to cell membrane integrity and viability. You can see the [scaffold visualizations proving this here](https://github.com/Zaki-1052/Yeast_MSA/tree/main/results/scaffold_variants/visualizations).

  

Anyways, those are all of the conclusions we can currently draw from the data; I also took the liberty of [generating an HTML report with visualizations here](https://github.com/Zaki-1052/Yeast_MSA/blob/main/results/reports/ergosterol_variant_analysis.html), so you can view them more easily than the CSV, I've attached it here as `ergosterol_variant_analysis.html`. If I have time tonight, I'll rerun the previous analyses with the gene coordinates/annotations, and maybe functional impact of some of the more prominent variants/genes. The upstream/regulatory classifications (and obviously the unusual distances of the variants from the ERG genes) definitely seems like it supports your hypothesis so far though! 

  

After the Chem midterm on Wednesday, I can definitely finish up the rest of the biological pathways connection from the functional impact analysis, and if you have the exact data frames from the sterol profiles, I can get those in as well on Thursday. Otherwise, it will mainly be guesstimates defined as more/less, without much statistical significance. Sorry this ended up getting so long, but hope it helped with your committee presentation and conclusions! Please let me know if you have any questions or suggestions for further steps beyond the functional impact and sterol correlations for this week!

  

Thanks,

Zaki

  

(Lastly, I'm attaching an LLM summary of this email report -- though you should definitely still take a look at the attached HTML with visualizations and stats -- so you can get just the key points at a glance along with the detailed analysis I wrote here!)

---


# Yeast Genome Variant Analysis: Summary Notes

## Analysis Overview
- Used **snpEff** (SNP Effect) to annotate variants from original analysis
- Cross-referenced variants with newest reference genome and **WGAP annotations**
- Previous non-annotated VCF analysis still applies to both specific and non-specific genes
- Analysis predicts variant **location**, **effect**, and impact on gene/protein function

## Key Findings: Variants in Ergosterol Pathway Genes
- **Low variant density**: Only **44 variants** found in gene-focused analysis vs. **2,087 variants** across all scaffolds
- **Purifying selection evidence**: Almost 75% of variants located >50kb from target genes
- **Variant distribution by gene**:
  - **YGR175C (ERG1)**: 29 variants within 5kb (all downstream 1kb) - ~50% of gene variants
  - **YGR060W (ERG25)**: 60 variants within 5kb - highest number
  - **YHR007C (ERG11)**: 31 variants within 5kb (all upstream)
  - **YMR202W (ERG2)** and **YLR056W (ERG3)**: 15 variants each within 5kb
  - **7 genes showed zero variants** within 5kb: ERG4, ERG7, ERG9, ERG6, ERG5, ERG24

## Variant Characteristics
- **Mutation types**:
  - **Insertions**: ~50% of all variants
  - **SNVs and deletions**: ~25% each
- **Regulatory emphasis**: 80% of variants are upstream gene variants (modifier region effects)
- **Few coding variants** affecting protein structure directly
- **Treatment distribution**: Similar across treatment groups
  - Control samples showed significantly fewer variants (4× difference)
  - No dramatic differences in impact distribution across treatments

## Biological Interpretation
- **Strong conservation** of ergosterol pathway genes (supports purifying/negative selection hypothesis)
- Mutations predominantly affect **gene regulation** rather than protein structure
- **ERG25** shows higher tolerance for nearby genetic variation compared to other pathway genes
- **No strong treatment-specific effects** observed in variant distribution
- **Conservation pattern** consistent with existing research on core enzyme conservation
- Mutations in ergosterol pathway likely detrimental to cell membrane integrity and viability

## Supporting Evidence
- 12th scaffold ("chromosome" XII) showed highest variant count but most were distant from target genes
- Variant distribution proportional to sample numbers across treatment groups
- Similar conclusions from different analytical approaches (gene-focused vs. scaffold-level)

## Next Steps
- **Functional impact analysis** of prominent variants/genes
- **Biological pathway connections** analysis after Chem midterm
- Integration with **sterol profiles** data
- HTML report with visualizations available as `ergosterol_variant_analysis.html`
