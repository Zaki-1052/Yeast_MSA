#!/usr/bin/env python3
"""
variant_gene_annotator.py - Annotate variants with gene information

This script analyzes VCF files and adds annotations related to genes of interest,
providing biological context to help interpret the variants.
"""

import os
import sys
import csv
import re
import argparse
from collections import defaultdict, Counter
import gzip
from Bio import SeqIO
from Bio.Seq import Seq

# Standard genetic code
GENETIC_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

class VariantGeneAnnotator:
    def __init__(self, scaffold_mapping_file, gene_index_file, genes_of_interest_file, 
                 genes_fasta_file, upstream_bp=1000, downstream_bp=500):
        """
        Initialize the VariantGeneAnnotator.
        
        Args:
            scaffold_mapping_file: TSV file mapping scaffolds to JRIU IDs
            gene_index_file: TSV file with gene locations
            genes_of_interest_file: TSV file with genes of interest
            genes_fasta_file: FASTA file with gene sequences
            upstream_bp: Number of base pairs upstream to consider proximal
            downstream_bp: Number of base pairs downstream to consider proximal
        """
        self.upstream_bp = upstream_bp
        self.downstream_bp = downstream_bp
        
        # Load scaffold mapping
        self.jriu_to_scaffold = {}
        self.scaffold_to_jriu = {}
        self.load_scaffold_mapping(scaffold_mapping_file)
        
        # Load gene information
        self.genes = {}
        self.gene_locations = defaultdict(list)
        self.load_gene_index(gene_index_file)
        
        # Load genes of interest
        self.genes_of_interest = set()
        self.sc_gene_id_to_w303 = {}
        self.load_genes_of_interest(genes_of_interest_file)
        
        # Load gene sequences
        self.gene_sequences = {}
        self.load_gene_sequences(genes_fasta_file)
        
        # Treatment groups
        self.treatment_groups = {
            'WT': ['WT-CTRL'],
            'WT-37': ['WT-37-55-1', 'WT-37-55-2', 'WT-37-55-3'],
            'WTA': ['WTA-55-1', 'WTA-55-2', 'WTA-55-3'],
            'STC': ['STC-55-1', 'STC-55-2', 'STC-55-3'],
            'CAS': ['CAS-55-1', 'CAS-55-2', 'CAS-55-3']
        }
        
        # Reverse lookup for sample to treatment
        self.sample_to_treatment = {}
        for treatment, samples in self.treatment_groups.items():
            for sample in samples:
                self.sample_to_treatment[sample] = treatment
        
        # Initialize statistics
        self.stats = defaultdict(Counter)
    
    def load_scaffold_mapping(self, mapping_file):
        """Load scaffold to JRIU ID mapping from TSV file"""
        with open(mapping_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                scaffold_id = row['scaffold_id']
                jriu_id = row['jriu_id']
                
                if jriu_id and jriu_id != 'unknown':
                    self.jriu_to_scaffold[jriu_id] = scaffold_id
                    self.scaffold_to_jriu[scaffold_id] = jriu_id
    
    def load_gene_index(self, gene_index_file):
        """Load gene information from TSV file"""
        with open(gene_index_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_id = row['gene_id']
                scaffold_id = row['scaffold_id']
                start = int(row['start'])
                end = int(row['end'])
                strand = row['strand']
                product = row.get('product', '')
                
                # Store gene information
                self.genes[gene_id] = {
                    'gene_id': gene_id,
                    'scaffold_id': scaffold_id,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'product': product
                }
                
                # Index gene location by scaffold
                self.gene_locations[scaffold_id].append({
                    'gene_id': gene_id,
                    'start': start,
                    'end': end,
                    'strand': strand
                })
    
    def load_genes_of_interest(self, genes_of_interest_file):
        """Load genes of interest from TSV file"""
        with open(genes_of_interest_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sc_gene_id = row['sc_gene_id']
                w303_gene_id = row['w303_gene_id']
                
                self.genes_of_interest.add(w303_gene_id)
                self.sc_gene_id_to_w303[sc_gene_id] = w303_gene_id
    
    def load_gene_sequences(self, genes_fasta_file):
        """Load gene sequences from FASTA file"""
        for record in SeqIO.parse(genes_fasta_file, "fasta"):
            # Parse the FASTA header format: YHR190W|ERG9|W3030BY00190
            header_parts = record.id.split('|')
            if len(header_parts) >= 3:
                sc_gene_id = header_parts[0]
                gene_name = header_parts[1]
                w303_gene_id = header_parts[2]
                
                self.gene_sequences[w303_gene_id] = {
                    'sc_gene_id': sc_gene_id,
                    'gene_name': gene_name,
                    'sequence': str(record.seq)
                }
    
    def _check_variant_gene_overlap(self, jriu_id, position, gene_info):
        """Check if a variant overlaps with a gene or is proximal"""
        scaffold_id = self.jriu_to_scaffold.get(jriu_id)
        if not scaffold_id:
            return None
        
        # Check if the variant is within a gene
        gene_start = gene_info['start']
        gene_end = gene_info['end']
        strand = gene_info['strand']
        
        if gene_start <= position <= gene_end:
            # Variant is within the gene
            rel_position = position - gene_start
            if strand == '-':
                rel_position = gene_end - position
            
            return {
                'location': 'within',
                'distance': 0,
                'relative_position': rel_position
            }
        
        # Check if variant is upstream or downstream
        if strand == '+':
            if position < gene_start:
                distance = gene_start - position
                if distance <= self.upstream_bp:
                    return {
                        'location': 'upstream',
                        'distance': distance,
                        'relative_position': -distance
                    }
            else:  # position > gene_end
                distance = position - gene_end
                if distance <= self.downstream_bp:
                    return {
                        'location': 'downstream',
                        'distance': distance,
                        'relative_position': gene_end - gene_start + distance
                    }
        else:  # strand == '-'
            if position > gene_end:
                distance = position - gene_end
                if distance <= self.upstream_bp:
                    return {
                        'location': 'upstream',
                        'distance': distance,
                        'relative_position': -distance
                    }
            else:  # position < gene_start
                distance = gene_start - position
                if distance <= self.downstream_bp:
                    return {
                        'location': 'downstream',
                        'distance': distance,
                        'relative_position': gene_end - gene_start + distance
                    }
        
        return None
    
    def _predict_coding_impact(self, gene_info, rel_position, ref, alt, gene_seq):
        """Predict the impact of a coding variant"""
        # Skip if we don't have gene sequence
        if not gene_seq:
            return {
                'impact_type': 'unknown',
                'details': 'No gene sequence available'
            }
        
        # Determine codon position
        codon_pos = rel_position // 3
        base_in_codon = rel_position % 3
        
        # Too close to the end?
        if codon_pos * 3 + 3 > len(gene_seq):
            return {
                'impact_type': 'unknown',
                'details': 'Variant too close to gene end'
            }
        
        # Get the codon
        codon = gene_seq[codon_pos * 3:codon_pos * 3 + 3]
        
        # Create the modified codon
        mod_codon = list(codon)
        mod_codon[base_in_codon] = alt
        mod_codon = ''.join(mod_codon)
        
        # Translate codons
        orig_aa = GENETIC_CODE.get(codon.upper(), '?')
        mod_aa = GENETIC_CODE.get(mod_codon.upper(), '?')
        
        # Determine effect
        if orig_aa == mod_aa:
            impact_type = 'synonymous'
            details = f"Codon change: {codon} → {mod_codon}, AA unchanged: {orig_aa}"
        elif mod_aa == '*' and orig_aa != '*':
            impact_type = 'nonsense'
            details = f"Codon change: {codon} → {mod_codon}, AA change: {orig_aa} → STOP"
        elif orig_aa == '*' and mod_aa != '*':
            impact_type = 'stop_lost'
            details = f"Codon change: {codon} → {mod_codon}, AA change: STOP → {mod_aa}"
        else:
            impact_type = 'missense'
            details = f"Codon change: {codon} → {mod_codon}, AA change: {orig_aa} → {mod_aa}"
        
        return {
            'impact_type': impact_type,
            'details': details,
            'ref_codon': codon,
            'alt_codon': mod_codon,
            'ref_aa': orig_aa,
            'alt_aa': mod_aa,
            'codon_position': codon_pos,
            'base_in_codon': base_in_codon
        }
    
    def _predict_variant_impact(self, jriu_id, position, ref, alt, overlap_info, gene_info):
        """Predict the impact of a variant on a gene"""
        w303_gene_id = gene_info['gene_id']
        gene_sequence = None
        
        if w303_gene_id in self.gene_sequences:
            gene_sequence = self.gene_sequences[w303_gene_id]['sequence']
        
        location = overlap_info['location']
        rel_position = overlap_info['relative_position']
        
        if location == 'within':
            return self._predict_coding_impact(gene_info, rel_position, ref, alt, gene_sequence)
        elif location == 'upstream':
            return {
                'impact_type': 'upstream_variant',
                'details': f"Variant {overlap_info['distance']} bp upstream of gene"
            }
        elif location == 'downstream':
            return {
                'impact_type': 'downstream_variant',
                'details': f"Variant {overlap_info['distance']} bp downstream of gene"
            }
        else:
            return {
                'impact_type': 'unknown',
                'details': 'Unknown location relative to gene'
            }
    
    def find_overlapping_genes(self, jriu_id, position):
        """Find genes that overlap with a given position"""
        scaffold_id = self.jriu_to_scaffold.get(jriu_id)
        if not scaffold_id:
            return []
        
        overlapping_genes = []
        for gene_loc in self.gene_locations.get(scaffold_id, []):
            overlap = self._check_variant_gene_overlap(jriu_id, position, gene_loc)
            if overlap:
                gene_id = gene_loc['gene_id']
                gene_info = self.genes.get(gene_id, {})
                
                # Check if this is a gene of interest
                is_goi = gene_id in self.genes_of_interest
                
                # Add to results
                overlapping_genes.append({
                    'gene_id': gene_id,
                    'scaffold_id': scaffold_id,
                    'gene_start': gene_loc['start'],
                    'gene_end': gene_loc['end'],
                    'strand': gene_loc['strand'],
                    'overlap': overlap,
                    'is_gene_of_interest': is_goi
                })
        
        return overlapping_genes
    
    def annotate_vcf(self, vcf_file, output_file, goi_only=True):
        """
        Annotate variants in a VCF file with gene information
        
        Args:
            vcf_file: Input VCF file
            output_file: Output annotated VCF file
            goi_only: Only output variants related to genes of interest
        """
        # Determine if the file is gzipped
        is_gzipped = vcf_file.endswith('.gz')
        opener = gzip.open if is_gzipped else open
        mode = 'rt' if is_gzipped else 'r'
        
        # Extract sample name from filename
        sample_name = os.path.basename(vcf_file).split('.')[0]
        treatment_group = self.sample_to_treatment.get(sample_name, 'unknown')
        
        # Read VCF file
        variant_count = 0
        gene_variant_count = 0
        goi_variant_count = 0
        
        # Dictionary to store annotated variants
        annotated_variants = []
        
        with opener(vcf_file, mode) as f:
            # Process header lines and identify column indexes
            header_lines = []
            column_names = None
            
            for line in f:
                line = line.strip()
                if line.startswith('##'):
                    header_lines.append(line)
                    continue
                elif line.startswith('#CHROM'):
                    header_lines.append(line)
                    column_names = line.split('\t')
                    break
            
            if not column_names:
                raise ValueError(f"Invalid VCF file: {vcf_file} - no header line found")
            
            # Process variant lines
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                variant_count += 1
                fields = line.split('\t')
                
                # Extract variant information
                chrom = fields[0]  # JRIU ID
                pos = int(fields[1])
                variant_id = fields[2]
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                filter_val = fields[6]
                info = fields[7]
                
                # Find overlapping genes
                overlapping_genes = self.find_overlapping_genes(chrom, pos)
                
                # Skip if no overlapping genes and we only want genes of interest
                if not overlapping_genes:
                    continue
                
                gene_variant_count += 1
                
                # Check if any overlapping gene is a gene of interest
                has_goi = any(g['is_gene_of_interest'] for g in overlapping_genes)
                
                # Skip if no gene of interest and we only want genes of interest
                if goi_only and not has_goi:
                    continue
                
                goi_variant_count += 1
                
                # Get annotations for each overlapping gene
                gene_annotations = []
                for gene_overlap in overlapping_genes:
                    gene_id = gene_overlap['gene_id']
                    gene_info = self.genes.get(gene_id, {})
                    overlap_info = gene_overlap['overlap']
                    
                    # Predict impact
                    impact = self._predict_variant_impact(
                        chrom, pos, ref, alt, overlap_info, gene_info
                    )
                    
                    # Add to gene annotations
                    gene_annotations.append({
                        'gene_id': gene_id,
                        'is_gene_of_interest': gene_overlap['is_gene_of_interest'],
                        'location': overlap_info['location'],
                        'distance': overlap_info['distance'],
                        'impact_type': impact['impact_type'],
                        'impact_details': impact['details']
                    })
                    
                    # Update statistics
                    self.stats['sample_counts'][sample_name] += 1
                    self.stats['treatment_counts'][treatment_group] += 1
                    self.stats['gene_counts'][gene_id] += 1
                    self.stats['impact_counts'][impact['impact_type']] += 1
                    if gene_overlap['is_gene_of_interest']:
                        self.stats['goi_counts'][gene_id] += 1
                        self.stats['goi_sample_counts'][f"{gene_id}:{sample_name}"] += 1
                        self.stats['goi_treatment_counts'][f"{gene_id}:{treatment_group}"] += 1
                        self.stats['goi_impact_counts'][f"{gene_id}:{impact['impact_type']}"] += 1
                
                # Add to annotated variants
                annotated_variants.append({
                    'chrom': chrom,
                    'pos': pos,
                    'id': variant_id,
                    'ref': ref,
                    'alt': alt,
                    'qual': qual,
                    'filter': filter_val,
                    'info': info,
                    'sample': sample_name,
                    'treatment': treatment_group,
                    'genes': gene_annotations
                })
        
        # Write annotated VCF
        with open(output_file, 'w') as f:
            # Write custom header lines
            for line in header_lines[:-1]:  # All but the column names
                f.write(line + '\n')
            
            # Add annotation format headers
            f.write('##INFO=<ID=GeneHit,Number=.,Type=String,Description="Gene overlapped by this variant">\n')
            f.write('##INFO=<ID=GeneImpact,Number=.,Type=String,Description="Predicted impact on gene function">\n')
            f.write('##INFO=<ID=GOI,Number=0,Type=Flag,Description="Overlaps a gene of interest">\n')
            
            # Write column names
            f.write(header_lines[-1] + '\n')
            
            # Write variant lines with annotations
            for variant in annotated_variants:
                # Create annotated INFO field
                info = variant['info']
                
                # Add gene annotations to INFO field
                gene_hits = []
                gene_impacts = []
                has_goi = False
                
                for gene_annot in variant['genes']:
                    gene_id = gene_annot['gene_id']
                    impact = gene_annot['impact_type']
                    location = gene_annot['location']
                    
                    gene_hits.append(f"{gene_id}:{location}")
                    gene_impacts.append(f"{gene_id}:{impact}")
                    
                    if gene_annot['is_gene_of_interest']:
                        has_goi = True
                
                if gene_hits:
                    info += f";GeneHit={','.join(gene_hits)}"
                if gene_impacts:
                    info += f";GeneImpact={','.join(gene_impacts)}"
                if has_goi:
                    info += ";GOI"
                
                # Write the annotated variant line
                f.write(f"{variant['chrom']}\t{variant['pos']}\t{variant['id']}\t"
                       f"{variant['ref']}\t{variant['alt']}\t{variant['qual']}\t"
                       f"{variant['filter']}\t{info}\tGT\t1\n")
        
        # Return statistics
        return {
            'sample': sample_name,
            'treatment': treatment_group,
            'total_variants': variant_count,
            'gene_variants': gene_variant_count,
            'goi_variants': goi_variant_count
        }
    
    def annotate_vcf_batch(self, vcf_files, output_dir, goi_only=True):
        """
        Annotate a batch of VCF files
        
        Args:
            vcf_files: List of VCF files to annotate
            output_dir: Output directory for annotated files
            goi_only: Only output variants related to genes of interest
        """
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Process each VCF file
        results = []
        for vcf_file in vcf_files:
            sample_name = os.path.basename(vcf_file).split('.')[0]
            output_file = os.path.join(output_dir, f"{sample_name}.annotated.vcf")
            
            print(f"Annotating {vcf_file}...")
            result = self.annotate_vcf(vcf_file, output_file, goi_only)
            results.append(result)
            
            print(f"  Found {result['total_variants']} variants, "
                  f"{result['gene_variants']} overlapping genes, "
                  f"{result['goi_variants']} overlapping genes of interest")
        
        return results
    
    def generate_summary_report(self, output_dir):
        """Generate summary reports of the annotation"""
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Write sample counts
        sample_counts_file = os.path.join(output_dir, "sample_variant_counts.tsv")
        with open(sample_counts_file, 'w') as f:
            f.write("Sample\tTreatment\tVariantCount\n")
            for sample, count in sorted(self.stats['sample_counts'].items()):
                treatment = self.sample_to_treatment.get(sample, 'unknown')
                f.write(f"{sample}\t{treatment}\t{count}\n")
        
        # Write treatment counts
        treatment_counts_file = os.path.join(output_dir, "treatment_variant_counts.tsv")
        with open(treatment_counts_file, 'w') as f:
            f.write("Treatment\tVariantCount\tSampleCount\tAvgVariantsPerSample\n")
            for treatment, count in sorted(self.stats['treatment_counts'].items()):
                samples = len([s for s, t in self.sample_to_treatment.items() if t == treatment])
                avg = count / samples if samples > 0 else 0
                f.write(f"{treatment}\t{count}\t{samples}\t{avg:.2f}\n")
        
        # Write gene counts
        gene_counts_file = os.path.join(output_dir, "gene_variant_counts.tsv")
        with open(gene_counts_file, 'w') as f:
            f.write("GeneID\tIsGeneOfInterest\tVariantCount\n")
            for gene, count in sorted(self.stats['gene_counts'].items()):
                is_goi = "Yes" if gene in self.genes_of_interest else "No"
                f.write(f"{gene}\t{is_goi}\t{count}\n")
        
        # Write impact counts
        impact_counts_file = os.path.join(output_dir, "impact_counts.tsv")
        with open(impact_counts_file, 'w') as f:
            f.write("ImpactType\tCount\n")
            for impact, count in sorted(self.stats['impact_counts'].items()):
                f.write(f"{impact}\t{count}\n")
        
        # Write genes of interest detailed report
        goi_report_file = os.path.join(output_dir, "genes_of_interest_report.tsv")
        with open(goi_report_file, 'w') as f:
            f.write("GeneID\tSCGeneID\tGeneName\tTreatment\tVariantCount\t"
                    "SynonymousCount\tMissenseCount\tNonsenseCount\tOtherCount\n")
            
            # Get all genes of interest
            for w303_gene_id in sorted(self.genes_of_interest):
                # Find corresponding S. cerevisiae gene ID
                sc_gene_id = "unknown"
                gene_name = "unknown"
                
                if w303_gene_id in self.gene_sequences:
                    sc_gene_id = self.gene_sequences[w303_gene_id]['sc_gene_id']
                    gene_name = self.gene_sequences[w303_gene_id]['gene_name']
                
                # Get counts for each treatment
                for treatment in sorted(self.treatment_groups.keys()):
                    total_count = sum(1 for k, v in self.stats['goi_treatment_counts'].items() 
                                    if k.startswith(f"{w303_gene_id}:{treatment}"))
                    
                    # Impact counts
                    synonymous_count = sum(v for k, v in self.stats['goi_impact_counts'].items() 
                                         if k.startswith(f"{w303_gene_id}:") and k.endswith(":synonymous"))
                    
                    missense_count = sum(v for k, v in self.stats['goi_impact_counts'].items() 
                                       if k.startswith(f"{w303_gene_id}:") and k.endswith(":missense"))
                    
                    nonsense_count = sum(v for k, v in self.stats['goi_impact_counts'].items() 
                                       if k.startswith(f"{w303_gene_id}:") and k.endswith(":nonsense"))
                    
                    other_count = total_count - (synonymous_count + missense_count + nonsense_count)
                    
                    f.write(f"{w303_gene_id}\t{sc_gene_id}\t{gene_name}\t{treatment}\t"
                           f"{total_count}\t{synonymous_count}\t{missense_count}\t"
                           f"{nonsense_count}\t{other_count}\n")
        
        # Generate HTML report
        self._generate_html_report(output_dir)
        
        return {
            'sample_counts_file': sample_counts_file,
            'treatment_counts_file': treatment_counts_file,
            'gene_counts_file': gene_counts_file,
            'impact_counts_file': impact_counts_file,
            'goi_report_file': goi_report_file
        }
    
    def _generate_html_report(self, output_dir):
        """Generate an HTML report with interactive graphs"""
        html_file = os.path.join(output_dir, "annotation_report.html")
        
        # Generate treatment comparison data
        treatment_data = {}
        for treatment, count in sorted(self.stats['treatment_counts'].items()):
            treatment_data[treatment] = count
        
        # Generate gene of interest data
        goi_data = defaultdict(dict)
        for w303_gene_id in sorted(self.genes_of_interest):
            if w303_gene_id in self.gene_sequences:
                sc_gene_id = self.gene_sequences[w303_gene_id]['sc_gene_id']
                gene_name = self.gene_sequences[w303_gene_id]['gene_name']
                
                for treatment in sorted(self.treatment_groups.keys()):
                    key = f"{w303_gene_id}:{treatment}"
                    count = 0
                    for k, v in self.stats['goi_treatment_counts'].items():
                        if k.startswith(key):
                            count += v
                    
                    goi_data[gene_name][treatment] = count
        
        # Generate impact type data
        impact_data = {}
        for impact, count in sorted(self.stats['impact_counts'].items()):
            impact_data[impact] = count
        
        # Write HTML file
        with open(html_file, 'w') as f:
            f.write("""<!DOCTYPE html>
<html>
<head>
    <title>Variant Annotation Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2 { color: #333366; }
        .chart-container { width: 800px; height: 400px; margin: 20px 0; }
        .grid-container { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
    </style>
</head>
<body>
    <h1>Variant Annotation Report</h1>
    <p>Generated on: """ + f"{os.popen('date').read().strip()}" + """</p>
    
    <div class="grid-container">
        <div>
            <h2>Variants by Treatment</h2>
            <div class="chart-container">
                <canvas id="treatmentChart"></canvas>
            </div>
        </div>
        
        <div>
            <h2>Variant Impact Types</h2>
            <div class="chart-container">
                <canvas id="impactChart"></canvas>
            </div>
        </div>
    </div>
    
    <h2>Genes of Interest - Variant Counts by Treatment</h2>
    <div class="chart-container" style="width: 100%;">
        <canvas id="geneChart"></canvas>
    </div>
    
    <script>
        // Treatment chart
        const treatmentCtx = document.getElementById('treatmentChart').getContext('2d');
        const treatmentChart = new Chart(treatmentCtx, {
            type: 'bar',
            data: {
                labels: """ + str(list(treatment_data.keys())) + """,
                datasets: [{
                    label: 'Variant Count',
                    data: """ + str(list(treatment_data.values())) + """,
                    backgroundColor: [
                        'rgba(54, 162, 235, 0.5)',
                        'rgba(255, 99, 132, 0.5)',
                        'rgba(75, 192, 192, 0.5)',
                        'rgba(255, 206, 86, 0.5)',
                        'rgba(153, 102, 255, 0.5)'
                    ],
                    borderColor: [
                        'rgba(54, 162, 235, 1)',
                        'rgba(255, 99, 132, 1)',
                        'rgba(75, 192, 192, 1)',
                        'rgba(255, 206, 86, 1)',
                        'rgba(153, 102, 255, 1)'
                    ],
                    borderWidth: 1
                }]
            },
            options: {
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Variant Count'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Treatment'
                        }
                    }
                }
            }
        });
        
        // Impact chart
        const impactCtx = document.getElementById('impactChart').getContext('2d');
        const impactChart = new Chart(impactCtx, {
            type: 'pie',
            data: {
                labels: """ + str(list(impact_data.keys())) + """,
                datasets: [{
                    label: 'Impact Types',
                    data: """ + str(list(impact_data.values())) + """,
                    backgroundColor: [
                        'rgba(54, 162, 235, 0.5)',
                        'rgba(255, 99, 132, 0.5)',
                        'rgba(75, 192, 192, 0.5)',
                        'rgba(255, 206, 86, 0.5)',
                        'rgba(153, 102, 255, 0.5)',
                        'rgba(255, 159, 64, 0.5)',
                        'rgba(201, 203, 207, 0.5)'
                    ],
                    borderColor: [
                        'rgba(54, 162, 235, 1)',
                        'rgba(255, 99, 132, 1)',
                        'rgba(75, 192, 192, 1)',
                        'rgba(255, 206, 86, 1)',
                        'rgba(153, 102, 255, 1)',
                        'rgba(255, 159, 64, 1)',
                        'rgba(201, 203, 207, 1)'
                    ],
                    borderWidth: 1
                }]
            }
        });
        
        // Gene of interest chart
        const geneCtx = document.getElementById('geneChart').getContext('2d');
        const geneChart = new Chart(geneCtx, {
            type: 'bar',
            data: {
                labels: """ + str(list(goi_data.keys())) + """,
                datasets: [""")
            
            # Generate datasets for each treatment
            treatment_colors = {
                'WT': {'bg': 'rgba(54, 162, 235, 0.5)', 'border': 'rgba(54, 162, 235, 1)'},
                'WT-37': {'bg': 'rgba(255, 99, 132, 0.5)', 'border': 'rgba(255, 99, 132, 1)'},
                'WTA': {'bg': 'rgba(75, 192, 192, 0.5)', 'border': 'rgba(75, 192, 192, 1)'},
                'STC': {'bg': 'rgba(255, 206, 86, 0.5)', 'border': 'rgba(255, 206, 86, 1)'},
                'CAS': {'bg': 'rgba(153, 102, 255, 0.5)', 'border': 'rgba(153, 102, 255, 1)'}
            }
            
            datasets = []
            for i, treatment in enumerate(sorted(self.treatment_groups.keys())):
                dataset = {
                    'label': treatment,
                    'data': [goi_data[gene].get(treatment, 0) for gene in goi_data.keys()],
                    'backgroundColor': treatment_colors.get(treatment, {}).get('bg', f'rgba(100, 100, 100, 0.5)'),
                    'borderColor': treatment_colors.get(treatment, {}).get('border', f'rgba(100, 100, 100, 1)'),
                    'borderWidth': 1
                }
                datasets.append(dataset)
            
            f.write(',\n'.join([str(dataset).replace("'", '"') for dataset in datasets]))
            
            f.write("""
                ]
            },
            options: {
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Variant Count'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Gene of Interest'
                        }
                    }
                }
            }
        });
    </script>
</body>
</html>""")
        
        print(f"Generated HTML report at {html_file}")

def main():
    parser = argparse.ArgumentParser(description='Annotate variants with gene information')
    parser.add_argument('--vcf', required=True, nargs='+', help='VCF file(s) to annotate')
    parser.add_argument('--scaffold-mapping', required=True, help='Scaffold to JRIU mapping file')
    parser.add_argument('--gene-index', required=True, help='Gene index file')
    parser.add_argument('--genes-of-interest', required=True, help='Genes of interest file')
    parser.add_argument('--genes-fasta', required=True, help='Genes FASTA file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--upstream-bp', type=int, default=1000, help='Upstream base pairs')
    parser.add_argument('--downstream-bp', type=int, default=500, help='Downstream base pairs')
    parser.add_argument('--goi-only', action='store_true', help='Only output variants related to genes of interest')
    
    args = parser.parse_args()
    
    # Create annotator
    annotator = VariantGeneAnnotator(
        args.scaffold_mapping,
        args.gene_index,
        args.genes_of_interest,
        args.genes_fasta,
        args.upstream_bp,
        args.downstream_bp
    )
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Annotate VCF files
    annotated_vcf_dir = os.path.join(args.output_dir, 'annotated_vcfs')
    results = annotator.annotate_vcf_batch(args.vcf, annotated_vcf_dir, args.goi_only)
    
    # Generate summary report
    report_dir = os.path.join(args.output_dir, 'reports')
    report_files = annotator.generate_summary_report(report_dir)
    
    print("\nAnnotation Complete")
    print(f"Annotated VCF files: {annotated_vcf_dir}")
    print(f"Report files: {report_dir}")
    
    # Print summary statistics
    total_variants = sum(result['total_variants'] for result in results)
    total_gene_variants = sum(result['gene_variants'] for result in results)
    total_goi_variants = sum(result['goi_variants'] for result in results)
    
    print(f"\nProcessed {len(args.vcf)} VCF files")
    print(f"Total variants: {total_variants}")
    print(f"Variants overlapping genes: {total_gene_variants} ({total_gene_variants/total_variants*100:.2f}%)")
    print(f"Variants overlapping genes of interest: {total_goi_variants} ({total_goi_variants/total_variants*100:.2f}%)")

if __name__ == '__main__':
    main()