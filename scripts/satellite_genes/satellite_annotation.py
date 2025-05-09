#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/satellite_genes/satellite_annotation.py

"""
Script to gather functional annotations and GO terms for satellite genes.
This script extends the satellite gene identification by enriching the data
with functional annotations from reference databases.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import argparse
import sys
import requests
import json
import re
import time
from matplotlib.pyplot import figure

# Add the parent directory to the path to access shared modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.tools import ensure_dir, save_tsv, load_tsv, setup_plotting_style, save_plot

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Annotate satellite genes with functional information")
    parser.add_argument("--satellite-genes", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/satellite_genes/satellite_genes.tsv",
                      help="Path to satellite genes TSV file from identification step")
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

def merge_annotations(satellite_df, gene_df):
    """Merge satellite genes with annotations from gene mapping"""
    # Create a copy to avoid modifying the original
    annotated_df = satellite_df.copy()
    
    # Extract relevant annotation columns from gene_df
    annotation_columns = ['sc_gene_id', 'std_gene_name', 'product', 'note', 'function']
    annotation_columns = [col for col in annotation_columns if col in gene_df.columns]
    
    # Add any additional useful columns
    optional_columns = ['go_terms', 'gene_ontology', 'molecular_function', 'biological_process', 'cellular_component']
    for col in optional_columns:
        if col in gene_df.columns:
            annotation_columns.append(col)
    
    # Create a subset of gene_df with only the annotation columns
    gene_annotations = gene_df[annotation_columns]
    
    # Merge satellite_df with annotations based on satellite_gene_id
    annotated_df = pd.merge(
        annotated_df,
        gene_annotations,
        left_on='satellite_gene_id',
        right_on='sc_gene_id',
        how='left',
        suffixes=('', '_annotation')
    )
    
    # Fill missing values
    for col in annotated_df.columns:
        if annotated_df[col].dtype == 'object':
            annotated_df[col] = annotated_df[col].fillna('Unknown')
    
    return annotated_df

def categorize_functions(annotated_df):
    """Categorize satellite genes by function based on annotations"""
    # Initialize function category column
    annotated_df['function_category'] = 'Unknown'
    
    # Define functional categories and their keywords
    function_categories = {
        'Metabolism': ['metabolism', 'metabolic', 'biosynthesis', 'degradation', 'fatty acid', 'lipid', 'glucose', 
                      'catabolism', 'glycolysis', 'TCA', 'respiration', 'oxidation', 'reduction'],
        'Stress Response': ['stress', 'heat shock', 'oxidative', 'response', 'resistance', 'adaptation', 'temperature',
                           'oxygen', 'hypoxia', 'anaerobic', 'protection'],
        'Transcription': ['transcription', 'transcriptional', 'transcription factor', 'RNA polymerase', 'promoter', 
                        'gene expression', 'DNA-binding'],
        'Translation': ['translation', 'ribosome', 'ribosomal', 'tRNA', 'protein synthesis', 'elongation', 'initiation'],
        'Transport': ['transport', 'transporter', 'transmembrane', 'uptake', 'export', 'import', 'channel', 'pump',
                     'porter', 'exchange', 'trafficking'],
        'Signaling': ['signaling', 'signal', 'receptor', 'kinase', 'phosphatase', 'cascade', 'pathway', 'regulation',
                     'regulator', 'GTPase', 'phosphorylation'],
        'Cell Cycle': ['cell cycle', 'mitosis', 'meiosis', 'division', 'checkpoint', 'spindle', 'chromosome segregation'],
        'Protein Modification': ['modification', 'ubiquitin', 'proteolysis', 'degradation', 'proteasome', 'folding',
                               'chaperone', 'isomerase', 'glycosylation'],
        'Membrane': ['membrane', 'lipid raft', 'vesicle', 'endosome', 'vacuole', 'endoplasmic reticulum', 'golgi',
                    'ergosterol', 'sterol', 'sphingolipid', 'phospholipid'],
        'DNA Processes': ['DNA', 'replication', 'repair', 'recombination', 'nucleotide', 'helicase', 'telomere', 'chromatin'],
    }
    
    # Function to categorize based on keywords
    def categorize_by_keywords(row):
        # Combine relevant text fields for searching
        text_fields = []
        for field in ['product', 'note', 'function', 'molecular_function', 'biological_process']:
            if field in row and not pd.isna(row[field]):
                text_fields.append(str(row[field]).lower())
        
        search_text = ' '.join(text_fields)
        
        # Check each category
        for category, keywords in function_categories.items():
            if any(keyword.lower() in search_text for keyword in keywords):
                return category
        
        return 'Unknown'
    
    # Apply categorization
    annotated_df['function_category'] = annotated_df.apply(categorize_by_keywords, axis=1)
    
    return annotated_df

def analyze_go_terms(annotated_df):
    """Analyze GO terms in satellite genes if available"""
    # Check if GO term columns exist
    go_columns = ['go_terms', 'gene_ontology', 'molecular_function', 'biological_process', 'cellular_component']
    available_go_columns = [col for col in go_columns if col in annotated_df.columns]
    
    if not available_go_columns:
        print("No GO term columns found in the data")
        return annotated_df, None
    
    # Create a consolidated GO term field if multiple columns exist
    if len(available_go_columns) > 1:
        annotated_df['consolidated_go'] = annotated_df[available_go_columns].apply(
            lambda row: ' '.join([str(x) for x in row if not pd.isna(x) and str(x) != 'Unknown']), 
            axis=1
        )
    else:
        # Use the single available column
        annotated_df['consolidated_go'] = annotated_df[available_go_columns[0]]
    
    # Extract and count GO terms
    all_go_terms = []
    for terms in annotated_df['consolidated_go']:
        if pd.isna(terms) or terms == 'Unknown':
            continue
        
        # Extract GO terms - can be in multiple formats
        # Format 1: GO:0005634, GO:0016021
        # Format 2: cellular component (GO:0005634), molecular function (GO:0003677)
        go_matches = re.findall(r'GO:\d+', str(terms))
        all_go_terms.extend(go_matches)
    
    # Count occurrences
    go_counter = Counter(all_go_terms)
    go_counts = pd.DataFrame({
        'go_term': list(go_counter.keys()),
        'count': list(go_counter.values())
    }).sort_values('count', ascending=False)
    
    return annotated_df, go_counts

def visualize_annotations(annotated_df, go_counts, output_dir):
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
    category_by_erg = pd.crosstab(annotated_df['erg_gene_name'], annotated_df['function_category'])
    # Convert to proportions
    category_by_erg_prop = category_by_erg.div(category_by_erg.sum(axis=1), axis=0)
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
    
    # 4. Top GO terms if available
    if go_counts is not None and not go_counts.empty:
        plt.figure(figsize=(12, 8))
        top_go = go_counts.head(20)  # Top 20 GO terms
        sns.barplot(x='count', y='go_term', data=top_go)
        plt.title('Top 20 GO Terms Among Satellite Genes')
        plt.xlabel('Number of Genes')
        plt.ylabel('GO Term')
        plt.tight_layout()
        save_plot(plt, os.path.join(viz_dir, "satellite_top_go_terms.png"))
        plt.close()
    
    # 5. Distance vs function category
    plt.figure(figsize=(12, 8))
    sns.boxplot(x='function_category', y='distance_kb', data=annotated_df)
    plt.title('Distribution of Distances by Functional Category')
    plt.xlabel('Functional Category')
    plt.ylabel('Distance (kb)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_plot(plt, os.path.join(viz_dir, "satellite_distance_by_function.png"))
    plt.close()

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
            
            # Compile results
            result = {
                'sgd_molecular_function': '; '.join(go_annotations['molecular_function']),
                'sgd_biological_process': '; '.join(go_annotations['biological_process']),
                'sgd_cellular_component': '; '.join(go_annotations['cellular_component']),
                'sgd_phenotypes': '; '.join(phenotypes[:5]),  # Limit to top 5 phenotypes
                'sgd_description': data.get('description', '')
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
    
    print(f"Querying SGD for {len(annotated_df)} genes (this may take some time)...")
    
    # Get unique gene IDs to avoid redundant queries
    unique_genes = annotated_df['satellite_gene_id'].unique()
    total_genes = len(unique_genes)
    
    gene_data = {}
    for i, gene_id in enumerate(unique_genes):
        if i % 20 == 0:  # Progress update every 20 genes
            print(f"Processed {i}/{total_genes} genes ({i/total_genes*100:.1f}%)")
        
        gene_data[gene_id] = query_sgd_for_annotation(gene_id)
        time.sleep(0.5)  # Be nice to the API
    
    # Add SGD data to dataframe
    for i, row in annotated_df.iterrows():
        gene_id = row['satellite_gene_id']
        if gene_id in gene_data:
            for col in sgd_columns:
                annotated_df.at[i, col] = gene_data[gene_id].get(col, '')
    
    print(f"Completed SGD data enrichment for {len(gene_data)} unique genes")
    return annotated_df

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
    
    # Merge annotations
    print("Merging annotations from gene mapping...")
    annotated_df = merge_annotations(satellite_df, gene_df)
    
    # Categorize by function
    print("Categorizing satellite genes by function...")
    annotated_df = categorize_functions(annotated_df)
    
    # Analyze GO terms
    print("Analyzing GO terms...")
    annotated_df, go_counts = analyze_go_terms(annotated_df)
    
    # Query SGD for additional annotations if requested
    if args.sgd_query:
        print("Querying SGD for additional annotations...")
        annotated_df = enrich_with_sgd_data(annotated_df)
    
    # Save annotated data
    output_file = os.path.join(args.output_dir, "satellite_genes_annotated.tsv")
    save_tsv(annotated_df, output_file)
    print(f"Annotated satellite gene data saved to {output_file}")
    
    # Save GO term counts if available
    if go_counts is not None:
        go_output = os.path.join(args.output_dir, "satellite_go_terms.tsv")
        save_tsv(go_counts, go_output)
        print(f"GO term counts saved to {go_output}")
    
    # Create visualizations
    print("Creating visualizations...")
    visualize_annotations(annotated_df, go_counts, args.output_dir)
    
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
        
        f.write("\nTop GO Terms:\n")
        if go_counts is not None and not go_counts.empty:
            for _, row in go_counts.head(10).iterrows():
                f.write(f"- {row['go_term']}: {row['count']} genes\n")
        else:
            f.write("- No GO terms available in the data\n")
        
        # Category distribution by ERG gene
        f.write("\nFunctional Category Distribution by ERG Gene:\n")
        category_by_erg = pd.crosstab(annotated_df['erg_gene_name'], annotated_df['function_category'])
        for erg_gene, row in category_by_erg.iterrows():
            top_categories = row.sort_values(ascending=False).head(3)
            category_str = ', '.join([f"{cat} ({count})" for cat, count in top_categories.items()])
            f.write(f"- {erg_gene}: {category_str}\n")
        
        f.write("\nFor full annotations, see the satellite_genes_annotated.tsv file.")
    
    print(f"Summary report saved to {summary_file}")
    print("Satellite gene annotation complete!")

if __name__ == "__main__":
    main()