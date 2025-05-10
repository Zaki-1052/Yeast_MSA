#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/satellite_genes/satellite_annotation.py

"""
Script to gather functional annotations and GO terms for satellite genes.
This script extends the satellite gene identification by enriching the data
with functional annotations from the GenBank files and SGD database.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import argparse
import sys
import re
import time
from Bio import SeqIO, Entrez
import requests
import json

# Add the parent directory to the path to access shared modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.tools import ensure_dir, save_tsv, load_tsv, setup_plotting_style, save_plot

# Define functional categories and their keywords
FUNCTION_CATEGORIES = {
    'Metabolism': ['metabolism', 'metabolic', 'biosynthesis', 'degradation', 'fatty acid', 'lipid', 
                 'glucose', 'catabolism', 'glycolysis', 'TCA', 'respiration', 'oxidation', 
                 'reduction', 'synthase', 'synthetase', 'dehygrogenase', 'oxidase', 'reductase',
                 'transferase', 'hydrolase', 'kinase', 'phosphatase'],
    'Transport': ['transport', 'transporter', 'transmembrane', 'uptake', 'export', 'import', 
                'channel', 'pump', 'porter', 'exchange', 'trafficking', 'permease', 'carrier',
                'uptake', 'efflux'],
    'Transcription': ['transcription', 'transcriptional', 'transcription factor', 'RNA polymerase', 
                    'promoter', 'gene expression', 'DNA-binding', 'zinc finger', 'helicase',
                    'transcriptase', 'activator', 'repressor', 'silencing'],
    'Translation': ['translation', 'ribosome', 'ribosomal', 'tRNA', 'protein synthesis', 
                  'elongation', 'initiation', 'termination factor', 'aminoacyl', 'synthetase',
                  'peptidyl', 'translational'],
    'Protein Modification': ['chaperone', 'folding', 'protease', 'peptidase', 'isomerase', 
                           'glycosylation', 'ubiquitin', 'proteolysis', 'degradation', 
                           'proteasome', 'isomerase', 'modification', 'processing'],
    'Signaling': ['signaling', 'signal', 'receptor', 'kinase', 'phosphatase', 'cascade', 
                 'pathway', 'regulation', 'regulator', 'GTPase', 'phosphorylation',
                 'transduction', 'response'],
    'Stress Response': ['stress', 'heat shock', 'oxidative', 'response', 'resistance', 
                      'adaptation', 'temperature', 'oxygen', 'hypoxia', 'anaerobic', 
                      'protection', 'defense', 'shock', 'tolerance'],
    'Cell Cycle': ['cell cycle', 'mitosis', 'meiosis', 'division', 'checkpoint', 'spindle', 
                  'chromosome segregation', 'cytokinesis', 'replication', 'G1', 'S phase',
                  'G2', 'M phase', 'mitotic'],
    'Membrane': ['membrane', 'lipid raft', 'vesicle', 'endosome', 'vacuole', 'endoplasmic reticulum', 
               'golgi', 'ergosterol', 'sterol', 'sphingolipid', 'phospholipid', 'transmembrane',
               'integral membrane', 'plasma membrane', 'cell wall'],
    'DNA Processes': ['DNA', 'replication', 'repair', 'recombination', 'nucleotide', 'helicase', 
                    'telomere', 'chromatin', 'histone', 'nuclease', 'topoisomerase',
                    'polymerase', 'nucleosome', 'chromosome'],
}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Annotate satellite genes with functional information")
    parser.add_argument("--satellite-genes", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/satellite_genes/satellite_genes.tsv",
                      help="Path to satellite genes TSV file from identification step")
    parser.add_argument("--genbank-dir", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/w303_annotations",
                      help="Directory containing GenBank annotation files")
    parser.add_argument("--gene-mapping", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/reference/gene_mapping_full.tsv",
                      help="Path to gene mapping file with annotations")
    parser.add_argument("--output-dir", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/satellite_genes",
                      help="Directory to store output files")
    parser.add_argument("--sgd-query", action='store_true',
                      help="Query SGD database for additional annotations (requires internet)")
    return parser.parse_args()

def load_satellite_genes(satellite_file):
    """Load satellite gene data from identification step"""
    if not os.path.exists(satellite_file):
        print(f"ERROR: Satellite gene file not found: {satellite_file}")
        sys.exit(1)
    
    try:
        satellite_df = pd.read_csv(satellite_file, sep="\t")
        print(f"Loaded {len(satellite_df)} satellite genes from file")
        return satellite_df
    except Exception as e:
        print(f"ERROR: Failed to load satellite gene file: {e}")
        sys.exit(1)

def load_gene_mapping(gene_mapping_file):
    """Load gene mapping information with annotations"""
    if not os.path.exists(gene_mapping_file):
        print(f"ERROR: Gene mapping file not found: {gene_mapping_file}")
        sys.exit(1)
    
    try:
        gene_df = pd.read_csv(gene_mapping_file, sep="\t")
        print(f"Loaded gene mapping with {len(gene_df)} genes")
        return gene_df
    except Exception as e:
        print(f"ERROR: Failed to load gene mapping file: {e}")
        sys.exit(1)

def parse_genbank_files(genbank_dir):
    """Parse GenBank files to extract gene annotations"""
    print(f"Parsing GenBank files from {genbank_dir}...")
    
    gene_annotations = {}
    genbank_files = [f for f in os.listdir(genbank_dir) if f.endswith('.genbank')]
    
    for file in genbank_files:
        try:
            file_path = os.path.join(genbank_dir, file)
            for record in SeqIO.parse(file_path, "genbank"):
                for feature in record.features:
                    if feature.type == "gene":
                        gene_id = feature.qualifiers.get("gene", [""])[0]
                        
                        # Skip if we don't have a gene ID
                        if not gene_id:
                            continue
                        
                        # Find the corresponding CDS feature
                        cds = None
                        for cds_feature in record.features:
                            if cds_feature.type == "CDS" and cds_feature.qualifiers.get("gene", [""])[0] == gene_id:
                                cds = cds_feature
                                break
                        
                        if cds:
                            sc_gene_id = None
                            product = cds.qualifiers.get("product", ["hypothetical protein"])[0]
                            note = "; ".join(cds.qualifiers.get("note", []))
                            
                            # Try to extract the systematic gene ID (e.g., YNL331C) from notes
                            if note:
                                # Look for inference or similarity to known genes
                                sc_match = re.search(r'similar to.+?([Y][A-Z]{2}\d+[WC])', note)
                                if sc_match:
                                    sc_gene_id = sc_match.group(1)
                                
                                # Extract function information if available
                                function = note
                            else:
                                function = ""
                            
                            gene_annotations[gene_id] = {
                                "w303_gene_id": gene_id,
                                "sc_gene_id": sc_gene_id,
                                "product": product,
                                "note": note,
                                "function": function
                            }
                            
                            # If we found a systematic ID, also index it
                            if sc_gene_id:
                                gene_annotations[sc_gene_id] = gene_annotations[gene_id]
        
        except Exception as e:
            print(f"Error parsing GenBank file {file}: {e}")
    
    print(f"Extracted annotations for {len(gene_annotations)} genes from GenBank files")
    return gene_annotations

def categorize_by_function(annotations):
    """Categorize genes by function based on their annotations"""
    if not annotations:
        return "Unknown"
    
    # Create a combined text to search
    search_text = ""
    for key in ["product", "note", "function"]:
        if key in annotations and annotations[key]:
            search_text += " " + str(annotations[key]).lower()
    
    # If still empty or only has generic terms
    if not search_text or search_text.strip() in ["", "hypothetical protein", "unknown"]:
        return "Unknown"
    
    # Check each functional category
    for category, keywords in FUNCTION_CATEGORIES.items():
        for keyword in keywords:
            if keyword.lower() in search_text:
                return category
    
    # If we found annotations but couldn't categorize
    if len(search_text.strip()) > 20:  # Arbitrary length to filter out very short annotations
        return "Other"
    
    return "Unknown"

def query_sgd_for_annotation(gene_id):
    """Query Saccharomyces Genome Database (SGD) for additional annotation"""
    try:
        # Use the SGD API to get gene information
        base_url = "https://yeastgenome.org/backend/locus/"
        response = requests.get(f"{base_url}{gene_id}", timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            # Extract relevant information
            go_annotations = {
                'molecular_function': [],
                'biological_process': [],
                'cellular_component': []
            }
            
            if 'go_annotations' in data:
                for annotation in data['go_annotations']:
                    go_type = annotation.get('go_aspect', '').lower().replace(' ', '_')
                    if go_type in go_annotations:
                        term = annotation.get('go_term', {}).get('display_name', '')
                        if term:
                            go_annotations[go_type].append(term)
            
            # Extract phenotypes
            phenotypes = []
            if 'phenotype_annotations' in data:
                for pheno in data['phenotype_annotations']:
                    if 'phenotype' in pheno and 'display_name' in pheno['phenotype']:
                        phenotypes.append(pheno['phenotype']['display_name'])
            
            # Extract description
            description = data.get('description', '')
            
            # Compile results
            result = {
                'sgd_molecular_function': '; '.join(go_annotations['molecular_function']),
                'sgd_biological_process': '; '.join(go_annotations['biological_process']),
                'sgd_cellular_component': '; '.join(go_annotations['cellular_component']),
                'sgd_phenotypes': '; '.join(phenotypes[:5]),  # Limit to top 5 phenotypes
                'sgd_description': description
            }
            
            return result
        else:
            print(f"Warning: Failed to get SGD data for {gene_id}, status code {response.status_code}")
            return {}
    
    except Exception as e:
        print(f"Error querying SGD for {gene_id}: {e}")
        return {}

def enrich_with_sgd_data(annotated_df):
    """Add SGD data to the annotated dataframe"""
    # Create columns for SGD data
    sgd_columns = ['sgd_molecular_function', 'sgd_biological_process', 'sgd_cellular_component', 
                  'sgd_phenotypes', 'sgd_description']
    for col in sgd_columns:
        annotated_df[col] = ''
    
    print(f"Querying SGD for genes (this may take some time)...")
    
    # Get unique gene IDs to avoid redundant queries
    unique_genes = []
    for _, row in annotated_df.iterrows():
        gene_id = row['satellite_gene_id']
        gene_name = row['satellite_gene_name']
        
        if gene_id and gene_id != 'Unknown' and gene_id not in unique_genes:
            unique_genes.append(gene_id)
        
        if gene_name and gene_name != 'Unknown' and gene_name not in unique_genes:
            unique_genes.append(gene_name)
    
    total_genes = len(unique_genes)
    
    gene_data = {}
    for i, gene_id in enumerate(unique_genes):
        if i % 20 == 0:  # Progress update every 20 genes
            print(f"Processed {i}/{total_genes} genes ({i/total_genes*100:.1f}%)")
        
        gene_data[gene_id] = query_sgd_for_annotation(gene_id)
        time.sleep(0.5)  # Be nice to the API
    
    # Add SGD data to dataframe
    for i, row in annotated_df.iterrows():
        # Try satellite_gene_id first
        gene_id = row['satellite_gene_id']
        if gene_id in gene_data and gene_data[gene_id]:
            for col in sgd_columns:
                annotated_df.at[i, col] = gene_data[gene_id].get(col, '')
            continue
            
        # If no data for satellite_gene_id, try satellite_gene_name
        gene_name = row['satellite_gene_name']
        if gene_name in gene_data and gene_data[gene_name]:
            for col in sgd_columns:
                annotated_df.at[i, col] = gene_data[gene_name].get(col, '')
    
    print(f"Completed SGD data enrichment for {len(gene_data)} unique genes")
    return annotated_df

def annotate_satellite_genes(satellite_df, genbank_annotations, gene_df):
    """Annotate satellite genes with functional information"""
    # Create a copy to avoid modifying the original
    annotated_df = satellite_df.copy()
    
    # Initialize annotation columns
    annotated_df['gene_note'] = ''
    annotated_df['gene_function'] = ''
    annotated_df['function_category'] = 'Unknown'
    
    # Go through each satellite gene and add annotations
    for idx, row in annotated_df.iterrows():
        w303_gene_id = row.get('w303_gene_id', '')
        sc_gene_id = row.get('satellite_gene_id', '')
        
        # Try to get annotations from GenBank data
        annotations = None
        
        # Try using w303_gene_id
        if w303_gene_id and w303_gene_id in genbank_annotations:
            annotations = genbank_annotations[w303_gene_id]
        
        # If not found, try sc_gene_id
        elif sc_gene_id and sc_gene_id in genbank_annotations:
            annotations = genbank_annotations[sc_gene_id]
            
        # If not found in GenBank, try mapping file
        elif sc_gene_id:
            mapping_match = gene_df[gene_df['sc_gene_id'] == sc_gene_id]
            if not mapping_match.empty:
                match = mapping_match.iloc[0]
                annotations = {
                    "w303_gene_id": match.get('w303_gene_id', ''),
                    "sc_gene_id": sc_gene_id,
                    "product": match.get('product', 'hypothetical protein'),
                    "note": "",
                    "function": ""
                }
        
        # If we found annotations, add them to the dataframe
        if annotations:
            # Update satellite_gene_name if empty
            if pd.isna(row['satellite_gene_name']) or row['satellite_gene_name'] == '':
                annotated_df.at[idx, 'satellite_gene_name'] = annotations.get('sc_gene_id', '')
            
            # Update satellite_product if it's a generic value
            current_product = row.get('satellite_product', '')
            if pd.isna(current_product) or current_product in ['', 'hypothetical protein', 'Unknown']:
                annotated_df.at[idx, 'satellite_product'] = annotations.get('product', 'hypothetical protein')
            
            # Add note and function
            annotated_df.at[idx, 'gene_note'] = annotations.get('note', '')
            annotated_df.at[idx, 'gene_function'] = annotations.get('function', '')
            
            # Categorize by function
            annotated_df.at[idx, 'function_category'] = categorize_by_function(annotations)
    
    # Fill missing values
    for col in annotated_df.columns:
        if annotated_df[col].dtype == 'object':
            annotated_df[col] = annotated_df[col].fillna('Unknown')
    
    return annotated_df

def visualize_annotations(annotated_df, output_dir):
    """Create visualizations of the functional annotations"""
    setup_plotting_style()
    
    # Create output directory for visualizations
    viz_dir = os.path.join(output_dir, "visualizations")
    ensure_dir(viz_dir)
    
    # 1. Function category distribution
    plt.figure(figsize=(12, 8))
    category_counts = annotated_df['function_category'].value_counts()
    sns.barplot(x=category_counts.index, y=category_counts.values)
    plt.title('Distribution of Satellite Genes by Functional Category')
    plt.xlabel('Functional Category')
    plt.ylabel('Number of Genes')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_plot(plt, os.path.join(viz_dir, "satellite_function_categories.png"))
    plt.close()
    
    # 2. Function categories by ERG gene
    plt.figure(figsize=(14, 10))
    # Get counts of function_category for each erg_gene_name
    crosstab = pd.crosstab(annotated_df['erg_gene_name'], annotated_df['function_category'])
    # Convert to proportions for better visualization
    category_by_erg_prop = crosstab.div(crosstab.sum(axis=1), axis=0)
    sns.heatmap(category_by_erg_prop, cmap='YlGnBu', annot=False, cbar_kws={'label': 'Proportion'})
    plt.title('Functional Categories of Satellite Genes by ERG Gene')
    plt.xlabel('Functional Category')
    plt.ylabel('ERG Gene')
    plt.tight_layout()
    save_plot(plt, os.path.join(viz_dir, "satellite_function_by_erg.png"))
    plt.close()
    
    # 3. Function categories by relative position
    plt.figure(figsize=(12, 8))
    pos_by_func = pd.crosstab(annotated_df['relative_position'], annotated_df['function_category'])
    pos_by_func_prop = pos_by_func.div(pos_by_func.sum(axis=1), axis=0)
    sns.heatmap(pos_by_func_prop, cmap='YlGnBu', annot=True, fmt='.2f')
    plt.title('Functional Categories of Satellite Genes by Position')
    plt.tight_layout()
    save_plot(plt, os.path.join(viz_dir, "satellite_function_by_position.png"))
    plt.close()
    
    # 4. Distance vs function category
    plt.figure(figsize=(12, 8))
    sns.boxplot(x='function_category', y='distance_kb', data=annotated_df)
    plt.title('Distribution of Distances by Functional Category')
    plt.xlabel('Functional Category')
    plt.ylabel('Distance (kb)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_plot(plt, os.path.join(viz_dir, "satellite_distance_by_function.png"))
    plt.close()

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Ensure output directory exists
    ensure_dir(args.output_dir)
    
    print("======================================================")
    print("Satellite Gene Annotation")
    print("======================================================")
    
    # Load satellite genes from identification step
    satellite_df = load_satellite_genes(args.satellite_genes)
    
    # Load gene mapping data with annotations
    gene_df = load_gene_mapping(args.gene_mapping)
    
    # Parse GenBank files for better annotations
    genbank_annotations = parse_genbank_files(args.genbank_dir)
    
    # Annotate satellite genes with functional information
    print("Annotating satellite genes with functional information...")
    annotated_df = annotate_satellite_genes(satellite_df, genbank_annotations, gene_df)
    
    # Query SGD for additional annotations if requested
    if args.sgd_query:
        print("Querying SGD for additional annotations...")
        annotated_df = enrich_with_sgd_data(annotated_df)
    
    # Save annotated data
    output_file = os.path.join(args.output_dir, "satellite_genes_annotated.tsv")
    save_tsv(annotated_df, output_file)
    print(f"Annotated satellite gene data saved to {output_file}")
    
    # Create visualizations
    print("Creating visualizations...")
    visualize_annotations(annotated_df, args.output_dir)
    
    # Generate summary report
    summary_file = os.path.join(args.output_dir, "satellite_annotation_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Satellite Gene Annotation Summary\n")
        f.write("================================\n\n")
        
        f.write("Analysis Overview:\n")
        f.write(f"- Total satellite genes analyzed: {len(satellite_df)}\n")
        f.write(f"- SGD query performed: {'Yes' if args.sgd_query else 'No'}\n\n")
        
        f.write("Functional Categories:\n")
        category_counts = annotated_df['function_category'].value_counts()
        total = len(annotated_df)
        for category, count in category_counts.items():
            f.write(f"- {category}: {count} genes ({count/total*100:.1f}%)\n")
        
        # Add information about functions by ERG gene
        f.write("\nFunctional Category Distribution by ERG Gene:\n")
        for erg_gene in annotated_df['erg_gene_name'].unique():
            if pd.isna(erg_gene) or erg_gene == '':
                continue
                
            erg_genes = annotated_df[annotated_df['erg_gene_name'] == erg_gene]
            categories = erg_genes['function_category'].value_counts()
            
            # Get the top 3 categories for this ERG gene
            top_categories = []
            for category, count in categories.items():
                if len(top_categories) < 3:
                    top_categories.append(f"{category} ({count})")
            
            category_str = ', '.join(top_categories)
            f.write(f"- {erg_gene}: {category_str}\n")
        
        f.write("\nFor full annotations, see the satellite_genes_annotated.tsv file.")
    
    print(f"Summary report saved to {summary_file}")
    print("Satellite gene annotation complete!")

if __name__ == "__main__":
    main()