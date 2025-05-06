#!/usr/bin/env python3

'''
Treatment-Control Visualization with Gene Mapping

This script creates visualizations comparing treatment to control samples,
with specific focus on gene-level analysis. It enhances the original TC_Visualization.py
script by adding gene-specific functionality, allowing for visualization of variant
patterns within genes, particularly the ergosterol biosynthesis pathway genes which
may be under purifying selection.
'''

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch
import os
from collections import defaultdict
import logging

# Set up logging for file only, not console
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("analysis.log", mode='a')
    ]
)

# Define output directories
OUTPUT_DIR = "analysis/treatment_control_analysis"
GENE_OUTPUT_DIR = "analysis/genes_of_interest/treatment_control_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(GENE_OUTPUT_DIR, exist_ok=True)

# Initialize gene data structures
GENE_DATA = {}  # Dictionary mapping gene IDs to their details
SCAFFOLD_GENES = defaultdict(list)  # Dictionary mapping scaffolds to lists of genes
GENES_OF_INTEREST = set()  # Set of gene IDs involved in the ergosterol pathway

# Gene status colors
GENE_COLORS = {
    'ERG': '#2ca02c',    # Ergosterol pathway genes
    'Non-ERG': '#1f77b4', # Non-ergosterol genes
    'No Gene': '#7f7f7f'  # No gene
}

# Define treatment information for better biological context
TREATMENT_INFO = {
    'WT-37': {'description': 'Temperature-adapted wild type', 'adaptation': 'Temperature'},
    'WTA': {'description': 'Low oxygen-adapted wild type', 'adaptation': 'Low Oxygen'},
    'STC': {'description': 'STC gene with low oxygen adaptation', 'adaptation': 'Low Oxygen', 'gene': 'STC'},
    'CAS': {'description': 'CAS gene with temperature adaptation', 'adaptation': 'Temperature', 'gene': 'CAS'}
}

# Define a better color scheme for different treatment groups
TREATMENT_COLORS = {
    'WT-37': '#1b9e77',    # Temperature-adapted wild type
    'WTA': '#d95f02',      # Low oxygen-adapted wild type
    'STC': '#7570b3',      # STC gene with low oxygen adaptation
    'CAS': '#e7298a',      # CAS gene with temperature adaptation
    'STC-vs-STCCTRL': '#66a61e',  # STC with original control
    'CAS-vs-CASCTRL': '#e6ab02'   # CAS with original control
}

# Function to load gene mapping data
def load_gene_mapping():
    """
    Load gene mapping data from reference files.
    
    This function loads gene data from the reference directory, including:
    1. Gene coordinates and information from gene_mapping.tsv
    2. Genes of interest (ergosterol pathway) from genes_of_interest_mapping.tsv
    
    Returns:
        bool: True if data was loaded successfully, False otherwise
    """
    global GENE_DATA, SCAFFOLD_GENES, GENES_OF_INTEREST
    
    # Clear existing data
    GENE_DATA.clear()
    SCAFFOLD_GENES.clear()
    GENES_OF_INTEREST.clear()
    
    # Define possible file paths for gene mapping data
    gene_mapping_paths = [
        "reference/gene_mapping.tsv",
        "reference/w303_annotations/gene_mapping.tsv"
    ]
    
    # Load gene mapping data
    gene_mapping_file = None
    for path in gene_mapping_paths:
        if os.path.exists(path):
            gene_mapping_file = path
            break
    
    if gene_mapping_file:
        try:
            # Load gene mapping data
            gene_df = pd.read_csv(gene_mapping_file, sep='\t')
            print(f"Loaded {len(gene_df)} genes from {gene_mapping_file}")
            logging.info(f"Loaded {len(gene_df)} genes from {gene_mapping_file}")
            
            # Process each gene
            for _, row in gene_df.iterrows():
                gene_id = row['w303_gene_id']
                scaffold = row['w303_scaffold']
                
                # Store gene data
                GENE_DATA[gene_id] = {
                    'gene_id': gene_id,
                    'locus_tag': row['locus_tag'] if 'locus_tag' in row else None,
                    'sc_gene_id': row['sc_gene_id'] if 'sc_gene_id' in row else None,
                    'erg_name': row['erg_name'] if 'erg_name' in row else None,
                    'scaffold': scaffold,
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'strand': row['strand'] if 'strand' in row else None,
                    'product': row['product'] if 'product' in row else None
                }
                
                # Map scaffold to genes
                SCAFFOLD_GENES[scaffold].append(gene_id)
            
            # Load genes of interest (ergosterol pathway genes)
            goi_paths = [
                "reference/genes_of_interest_mapping.tsv",
                "reference/w303_annotations/genes_of_interest.tsv"
            ]
            
            goi_file = None
            for path in goi_paths:
                if os.path.exists(path):
                    goi_file = path
                    break
            
            if goi_file:
                try:
                    goi_df = pd.read_csv(goi_file, sep='\t')
                    print(f"Loaded {len(goi_df)} genes of interest from {goi_file}")
                    logging.info(f"Loaded {len(goi_df)} genes of interest from {goi_file}")
                    
                    # Add to our set of genes of interest
                    for _, row in goi_df.iterrows():
                        if 'w303_gene_id' in row:
                            GENES_OF_INTEREST.add(row['w303_gene_id'])
                except Exception as e:
                    print(f"Error loading genes of interest: {e}")
                    logging.error(f"Error loading genes of interest: {e}")
            else:
                print("No genes of interest file found. Using empty set.")
                logging.warning("No genes of interest file found.")
            
            return True
        except Exception as e:
            print(f"Error loading gene mapping data: {e}")
            logging.error(f"Error loading gene mapping data: {e}")
            return False
    else:
        print("No gene mapping file found. Gene mapping will not be available.")
        logging.warning("No gene mapping file found. Gene mapping will not be available.")
        return False

def map_variants_to_genes(variant_df):
    """
    Map variants to genes based on their genomic coordinates.
    
    Args:
        variant_df (pandas.DataFrame): DataFrame containing variant information
    
    Returns:
        pandas.DataFrame: The original DataFrame with additional gene-related columns
    """
    if len(GENE_DATA) == 0 or len(SCAFFOLD_GENES) == 0:
        print("Gene mapping data not loaded. Cannot map variants to genes.")
        logging.warning("Gene mapping data not loaded. Cannot map variants to genes.")
        return variant_df
    
    # Create a copy of the input dataframe
    result_df = variant_df.copy()
    
    # Initialize gene-related columns
    result_df['in_gene'] = False
    result_df['gene_id'] = None
    result_df['gene_name'] = None
    result_df['gene_type'] = None
    result_df['gene_product'] = None
    
    # Process variant IDs to extract scaffold and position
    if 'Variant_ID' in result_df.columns:
        # Extract scaffold and position from variant IDs (format: scaffold_position_ref_alt)
        result_df[['CHROM', 'POS', 'REF', 'ALT']] = result_df['Variant_ID'].str.split('_', expand=True)
        result_df['POS'] = pd.to_numeric(result_df['POS'], errors='coerce')
    
    # Map each variant to genes
    for idx, row in result_df.iterrows():
        if 'CHROM' not in row or pd.isna(row['CHROM']) or 'POS' not in row or pd.isna(row['POS']):
            continue
            
        scaffold = row['CHROM']
        position = row['POS']
        
        # Skip if scaffold has no mapped genes
        if scaffold not in SCAFFOLD_GENES:
            continue
        
        # Check each gene in this scaffold
        for gene_id in SCAFFOLD_GENES[scaffold]:
            gene_data = GENE_DATA[gene_id]
            
            # Check if position falls within gene coordinates
            if gene_data['start'] <= position <= gene_data['end']:
                result_df.at[idx, 'in_gene'] = True
                result_df.at[idx, 'gene_id'] = gene_id
                
                # Add gene name if available
                if gene_data['erg_name']:
                    result_df.at[idx, 'gene_name'] = gene_data['erg_name']
                elif gene_data['sc_gene_id']:
                    result_df.at[idx, 'gene_name'] = gene_data['sc_gene_id']
                
                # Set gene type based on presence in genes of interest
                if gene_id in GENES_OF_INTEREST:
                    result_df.at[idx, 'gene_type'] = 'ergosterol'
                else:
                    result_df.at[idx, 'gene_type'] = 'other'
                
                # Add gene product description if available
                if gene_data['product']:
                    result_df.at[idx, 'gene_product'] = gene_data['product']
                
                # Break since we found a matching gene
                break
    
    # Log the mapping results
    in_gene_count = sum(result_df['in_gene'])
    ergosterol_count = sum(result_df['gene_type'] == 'ergosterol')
    
    print(f"Mapped {in_gene_count} out of {len(result_df)} variants to genes")
    print(f"Found {ergosterol_count} variants in ergosterol pathway genes")
    
    logging.info(f"Mapped {in_gene_count} out of {len(result_df)} variants to genes")
    logging.info(f"Found {ergosterol_count} variants in ergosterol pathway genes")
    
    return result_df

def create_gene_specific_visualizations(results_df, variant_data=None):
    """
    Create gene-specific visualizations for treatment vs control analysis.
    
    Args:
        results_df (pandas.DataFrame): DataFrame with treatment vs control comparison results
        variant_data (pandas.DataFrame, optional): DataFrame with detailed variant information
                                                 mapped to genes
    """
    if variant_data is None or len(variant_data) == 0:
        print("No gene-specific variant data available. Skipping gene-specific visualizations.")
        return
    
    # Create directory for gene-specific visualizations
    os.makedirs(GENE_OUTPUT_DIR, exist_ok=True)
    
    # Count variants by gene status (ergosterol vs. other)
    print("Creating gene-specific visualizations...")
    
    # Aggregate gene status for each treatment
    gene_status_counts = defaultdict(lambda: defaultdict(int))
    
    # For each treatment, count variants by gene status
    for treatment in results_df['Treatment'].unique():
        # Get variants for this treatment
        treatment_variants = variant_data[variant_data['Treatment'] == treatment]
        
        # Count by gene status
        for status in ['ergosterol', 'other']:
            gene_status_counts[treatment][status] = sum(treatment_variants['gene_type'] == status)
        
        # Count non-genic variants
        gene_status_counts[treatment]['non-genic'] = sum(~treatment_variants['in_gene'])
    
    # Convert to DataFrame for plotting
    gene_data = []
    for treatment, counts in gene_status_counts.items():
        for status, count in counts.items():
            gene_status = 'ERG' if status == 'ergosterol' else 'Non-ERG' if status == 'other' else 'No Gene'
            gene_data.append({
                'Treatment': treatment,
                'Gene Status': gene_status,
                'Variants': count
            })
    
    gene_df = pd.DataFrame(gene_data)
    
    # 1. Gene Status Distribution by Treatment
    plt.figure(figsize=(14, 8))
    
    # Create a grouped bar chart
    ax = sns.barplot(
        data=gene_df,
        x='Treatment',
        y='Variants',
        hue='Gene Status',
        palette=GENE_COLORS
    )
    
    # Add value labels
    for container in ax.containers:
        ax.bar_label(container, fmt='%d')
    
    # Customize plot
    plt.title('Distribution of Variants by Gene Status Across Treatments', fontsize=14)
    plt.xlabel('Treatment')
    plt.ylabel('Number of Variants')
    plt.legend(title='Gene Status')
    
    # Save figure
    plt.tight_layout()
    plt.savefig(f'{GENE_OUTPUT_DIR}/gene_status_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Gene Status Proportion by Treatment (Stacked)
    plt.figure(figsize=(14, 8))
    
    # Pivot data for stacked plot
    pivot_df = gene_df.pivot(index='Treatment', columns='Gene Status', values='Variants').fillna(0)
    
    # Calculate proportions
    pivot_df = pivot_df.div(pivot_df.sum(axis=1), axis=0) * 100
    
    # Plot stacked bars
    pivot_df.plot(
        kind='bar', 
        stacked=True,
        color=[GENE_COLORS[status] for status in pivot_df.columns],
        figsize=(14, 8)
    )
    
    # Customize plot
    plt.title('Proportion of Variants by Gene Status Across Treatments', fontsize=14)
    plt.xlabel('Treatment')
    plt.ylabel('Percentage of Variants')
    plt.legend(title='Gene Status')
    
    # Format percentages
    for container in plt.gca().containers:
        plt.gca().bar_label(container, label_type='center', fmt='%.1f%%')
    
    # Save figure
    plt.tight_layout()
    plt.savefig(f'{GENE_OUTPUT_DIR}/gene_status_proportion.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Ergosterol Gene Variant Distribution
    if sum(variant_data['gene_type'] == 'ergosterol') > 0:
        # Get ergosterol gene variants
        erg_variants = variant_data[variant_data['gene_type'] == 'ergosterol']
        
        # Count variants by gene
        gene_counts = erg_variants.groupby(['gene_name', 'Treatment']).size().reset_index(name='Variants')
        
        # Filter for genes with at least one variant
        gene_counts = gene_counts[gene_counts['Variants'] > 0]
        
        if len(gene_counts) > 0:
            plt.figure(figsize=(14, 10))
            
            # Create a grouped bar chart
            ax = sns.barplot(
                data=gene_counts,
                x='gene_name',
                y='Variants',
                hue='Treatment',
                palette=TREATMENT_COLORS
            )
            
            # Customize plot
            plt.title('Ergosterol Pathway Gene Variant Distribution by Treatment', fontsize=14)
            plt.xlabel('Gene')
            plt.ylabel('Number of Variants')
            plt.xticks(rotation=45, ha='right')
            
            # Add value labels
            for container in ax.containers:
                ax.bar_label(container, fmt='%d')
            
            # Save figure
            plt.tight_layout()
            plt.savefig(f'{GENE_OUTPUT_DIR}/erg_gene_distribution.png', dpi=300, bbox_inches='tight')
            plt.close()
    
    # 4. Fold Change by Gene Status
    if 'Fold_Change' in results_df.columns and len(results_df) > 0:
        # Prepare data
        fc_data = []
        for _, row in results_df.iterrows():
            treatment = row['Treatment']
            fold_change = row['Fold_Change']
            
            # Get treatment variants
            treatment_variants = variant_data[variant_data['Treatment'] == treatment]
            
            # Calculate counts
            erg_count = sum(treatment_variants['gene_type'] == 'ergosterol')
            non_erg_count = sum(treatment_variants['gene_type'] == 'other')
            non_gene_count = sum(~treatment_variants['in_gene'])
            
            # Add to data list
            if erg_count > 0:
                fc_data.append({
                    'Treatment': treatment,
                    'Gene Status': 'ERG',
                    'Fold Change': fold_change,
                    'Variants': erg_count
                })
            
            if non_erg_count > 0:
                fc_data.append({
                    'Treatment': treatment,
                    'Gene Status': 'Non-ERG',
                    'Fold Change': fold_change,
                    'Variants': non_erg_count
                })
            
            if non_gene_count > 0:
                fc_data.append({
                    'Treatment': treatment,
                    'Gene Status': 'No Gene',
                    'Fold Change': fold_change,
                    'Variants': non_gene_count
                })
        
        # Create DataFrame
        fc_df = pd.DataFrame(fc_data)
        
        if len(fc_df) > 0:
            plt.figure(figsize=(14, 8))
            
            # Plot fold change with bubble size proportional to variant count
            for status, color in GENE_COLORS.items():
                status_data = fc_df[fc_df['Gene Status'] == status]
                if len(status_data) > 0:
                    plt.scatter(
                        status_data['Treatment'],
                        status_data['Fold Change'],
                        s=status_data['Variants'] * 10,  # Scale bubble size
                        alpha=0.7,
                        color=color,
                        label=status
                    )
            
            # Add horizontal line at fold change = 1
            plt.axhline(y=1, color='black', linestyle='--', alpha=0.5)
            
            # Customize plot
            plt.title('Fold Change by Gene Status Across Treatments', fontsize=14)
            plt.xlabel('Treatment')
            plt.ylabel('Fold Change (Treatment/Control)')
            plt.legend(title='Gene Status')
            plt.xticks(rotation=45, ha='right')
            
            # Save figure
            plt.tight_layout()
            plt.savefig(f'{GENE_OUTPUT_DIR}/fold_change_by_gene_status.png', dpi=300, bbox_inches='tight')
            plt.close()
    
    # 5. Create a gene-specific summary report
    create_gene_specific_report(results_df, variant_data)

def create_gene_specific_report(results_df, variant_data):
    """
    Create a comprehensive gene-specific summary report.
    
    Args:
        results_df (pandas.DataFrame): DataFrame with treatment vs control comparison results
        variant_data (pandas.DataFrame): DataFrame with detailed variant information
    """
    if variant_data is None or len(variant_data) == 0:
        return
    
    # Prepare report
    report_lines = [
        "# Gene-Specific Analysis of Treatment vs Control Data",
        "\n## Overall Statistics\n"
    ]
    
    # Add overall statistics
    total_variants = len(variant_data)
    in_gene_count = sum(variant_data['in_gene'])
    erg_count = sum(variant_data['gene_type'] == 'ergosterol')
    non_erg_count = sum(variant_data['gene_type'] == 'other')
    non_gene_count = total_variants - in_gene_count
    
    report_lines.extend([
        f"Total variants: {total_variants}",
        f"Variants in genes: {in_gene_count} ({in_gene_count/total_variants*100:.1f}%)",
        f"Variants in ergosterol pathway genes: {erg_count} ({erg_count/total_variants*100:.1f}%)",
        f"Variants in other genes: {non_erg_count} ({non_erg_count/total_variants*100:.1f}%)",
        f"Variants in non-genic regions: {non_gene_count} ({non_gene_count/total_variants*100:.1f}%)",
        "\n## Variants by Treatment and Gene Status\n"
    ])
    
    # Add treatment-specific statistics
    for treatment in results_df['Treatment'].unique():
        treatment_variants = variant_data[variant_data['Treatment'] == treatment]
        if len(treatment_variants) == 0:
            continue
            
        treatment_total = len(treatment_variants)
        treatment_in_gene = sum(treatment_variants['in_gene'])
        treatment_erg = sum(treatment_variants['gene_type'] == 'ergosterol')
        treatment_non_erg = sum(treatment_variants['gene_type'] == 'other')
        treatment_non_gene = treatment_total - treatment_in_gene
        
        report_lines.extend([
            f"### {treatment}",
            f"Total variants: {treatment_total}",
            f"Variants in genes: {treatment_in_gene} ({treatment_in_gene/treatment_total*100:.1f}%)",
            f"Variants in ergosterol pathway genes: {treatment_erg} ({treatment_erg/treatment_total*100:.1f}%)",
            f"Variants in other genes: {treatment_non_erg} ({treatment_non_erg/treatment_total*100:.1f}%)",
            f"Variants in non-genic regions: {treatment_non_gene} ({treatment_non_gene/treatment_total*100:.1f}%)",
            ""
        ])
    
    # Add ergosterol pathway gene statistics if present
    if erg_count > 0:
        report_lines.extend([
            "\n## Variants in Ergosterol Pathway Genes\n"
        ])
        
        # Group by gene name and treatment
        erg_variants = variant_data[variant_data['gene_type'] == 'ergosterol']
        gene_counts = erg_variants.groupby(['gene_name', 'Treatment']).size().reset_index(name='Count')
        
        for gene_name in sorted(gene_counts['gene_name'].unique()):
            gene_data = gene_counts[gene_counts['gene_name'] == gene_name]
            
            report_lines.extend([
                f"### {gene_name}",
                "| Treatment | Variant Count |",
                "| --- | --- |"
            ])
            
            for _, row in gene_data.iterrows():
                report_lines.append(f"| {row['Treatment']} | {row['Count']} |")
            
            report_lines.append("")
    
    # Add purifying selection analysis
    report_lines.extend([
        "\n## Purifying Selection Analysis\n",
        "Evidence of purifying selection can be assessed by comparing the proportion of variants",
        "in ergosterol pathway genes versus non-ergosterol genes and non-genic regions."
    ])
    
    # Calculate overall proportions
    if total_variants > 0 and in_gene_count > 0:
        # Calculate expected proportion based on genome coverage
        erg_genes_count = len(GENES_OF_INTEREST)
        total_genes_count = len(GENE_DATA)
        non_erg_genes_count = total_genes_count - erg_genes_count
        
        # Calculate gene length (assuming we have this information)
        erg_length = sum(GENE_DATA[gene_id]['end'] - GENE_DATA[gene_id]['start'] + 1 
                         for gene_id in GENES_OF_INTEREST if gene_id in GENE_DATA)
        total_length = sum(gene['end'] - gene['start'] + 1 for gene in GENE_DATA.values())
        non_erg_length = total_length - erg_length
        
        # Calculate expected proportions
        expected_erg_prop = erg_length / total_length
        expected_non_erg_prop = non_erg_length / total_length
        
        # Calculate actual proportions
        actual_erg_prop = erg_count / in_gene_count
        actual_non_erg_prop = non_erg_count / in_gene_count
        
        # Calculate fold difference (lower values suggest purifying selection)
        erg_fold_diff = actual_erg_prop / expected_erg_prop if expected_erg_prop > 0 else 0
        
        report_lines.extend([
            f"\nErgosterol genes make up {erg_genes_count} out of {total_genes_count} genes ({erg_genes_count/total_genes_count*100:.1f}%)",
            f"Ergosterol genes cover approximately {erg_length} bases out of {total_length} bases in genes ({erg_length/total_length*100:.1f}%)",
            f"Expected proportion of variants in ergosterol genes: {expected_erg_prop:.3f}",
            f"Actual proportion of variants in ergosterol genes: {actual_erg_prop:.3f}",
            f"Fold difference: {erg_fold_diff:.2f}",
            ""
        ])
        
        # Interpret results
        if erg_fold_diff < 0.75:
            report_lines.append("This suggests STRONG purifying selection in ergosterol pathway genes (significantly fewer mutations than expected).")
        elif erg_fold_diff < 0.9:
            report_lines.append("This suggests MODERATE purifying selection in ergosterol pathway genes (fewer mutations than expected).")
        elif erg_fold_diff < 1.1:
            report_lines.append("This suggests NEUTRAL selection in ergosterol pathway genes (mutations as expected).")
        else:
            report_lines.append("This suggests POSITIVE selection in ergosterol pathway genes (more mutations than expected).")
    
    # Write report to file
    with open(f'{GENE_OUTPUT_DIR}/gene_specific_report.md', 'w') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Gene-specific report written to {GENE_OUTPUT_DIR}/gene_specific_report.md")

def create_treatment_control_visualizations(results_df):
    """Create intuitive visualizations of treatment vs control differences."""
    
    # Debug Point 1: Check input data
    print("\nChecking input data:")
    print(results_df.head())
    print("\nChecking for NaN values:")
    print(results_df.isna().sum())
    
    # Set style parameters
    plt.rcParams['figure.figsize'] = [10, 6]
    plt.rcParams['font.size'] = 12
    sns.set_palette("husl")
    
    # Debug Point 2: Print value ranges
    print("\nValue ranges:")
    for col in results_df.columns:
        if pd.api.types.is_numeric_dtype(results_df[col]):
            print(f"{col}: min={results_df[col].min()}, max={results_df[col].max()}")
    
    # Group treatments by control for better visualization
    results_df['Control_Group'] = results_df['Control'].apply(lambda x: x.split('-')[0])
    
    # Define a better color scheme for different treatment groups
    color_map = {
        'WT-37': '#1b9e77',    # Temperature-adapted wild type
        'WTA': '#d95f02',      # Low oxygen-adapted wild type
        'STC': '#7570b3',      # STC gene with low oxygen adaptation
        'CAS': '#e7298a',      # CAS gene with temperature adaptation
        'STC-vs-STCCTRL': '#66a61e',  # STC with original control
        'CAS-vs-CASCTRL': '#e6ab02'   # CAS with original control
    }
    
    # 1. Bar Plot with Control vs Treatment
    plt.figure(figsize=(14, 8))
    
    # Prepare data
    treatments = results_df['Treatment']
    x = np.arange(len(treatments))
    width = 0.35
    
    # Debug Point 3: Print bar plot data
    print("\nBar plot data:")
    print("Control values:", results_df['Control_Variants'].tolist())
    print("Treatment values:", results_df['Treatment_Variants'].tolist())
    
    # Create bars with error checking
    control_vals = np.nan_to_num(results_df['Control_Variants'].values, 0)
    treatment_vals = np.nan_to_num(results_df['Treatment_Variants'].values, 0)
    
    plt.bar(x - width/2, control_vals, width, label='Control', color='lightgray')
    treatment_bars = plt.bar(x + width/2, treatment_vals, width, label='Treatment')
    
    # Apply custom colors to treatment bars
    for i, bar in enumerate(treatment_bars):
        treatment = treatments.iloc[i]
        bar.set_color(color_map.get(treatment, '#2ecc71'))
    
    # Customize plot
    plt.xlabel('Treatment Group')
    plt.ylabel('Number of Variants')
    plt.title('Comparison of Variants: Treatment vs Control', pad=20)
    
    # Custom x-tick labels that include control names
    x_labels = [f"{t}\n(vs {c})" for t, c in zip(results_df['Treatment'], results_df['Control'])]
    plt.xticks(x, x_labels, rotation=45, ha='right')
    
    plt.legend()
    
    # Add treatment effect arrows
    for i, (control, treatment) in enumerate(zip(control_vals, treatment_vals)):
        plt.arrow(x[i], control, 0, treatment-control,
                head_width=0.1, head_length=min(20, max(5, (treatment-control)/10)),
                fc='#e74c3c', ec='#e74c3c',
                alpha=0.5)
    
    # Add significance stars with error checking
    for i, row in enumerate(results_df.itertuples()):
        if pd.notna(row.Q_Value) and pd.notna(row.Treatment_Variants) and pd.notna(row.Control_Variants):
            stars = '***' if row.Q_Value < 0.001 else '**' if row.Q_Value < 0.01 else '*' if row.Q_Value < 0.05 else 'ns'
            max_height = max(row.Treatment_Variants, row.Control_Variants)
            plt.text(i, max_height * 1.1, stars, ha='center', va='bottom')
    
    # Add fold change annotations with error checking
    for i, row in enumerate(results_df.itertuples()):
        if pd.notna(row.Fold_Change) and pd.notna(row.Treatment_Variants):
            plt.text(i, row.Treatment_Variants * 1.02,
                    f'FC: {row.Fold_Change:.1f}×', 
                    ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('analysis/treatment_control_analysis/variant_comparison_barplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Grouped Bar Plot with Primary vs Secondary Comparisons
    primary_comparisons = results_df[results_df['Control'] == 'WT-CTRL']
    secondary_comparisons = results_df[results_df['Control'] != 'WT-CTRL']
    
    if not primary_comparisons.empty and not secondary_comparisons.empty:
        plt.figure(figsize=(14, 8))
        
        # Plot primary comparisons
        ax1 = plt.subplot(1, 2, 1)
        primary_treatments = primary_comparisons['Treatment']
        x_prim = np.arange(len(primary_treatments))
        
        bars1 = plt.bar(x_prim, primary_comparisons['Treatment_Variants'])
        for i, bar in enumerate(bars1):
            treatment = primary_treatments.iloc[i]
            bar.set_color(color_map.get(treatment, '#2ecc71'))
        
        plt.axhline(y=primary_comparisons['Control_Variants'].iloc[0], color='red', linestyle='--', 
                  label=f'WT-CTRL ({primary_comparisons["Control_Variants"].iloc[0]} variants)')
        
        plt.xlabel('Treatment Group')
        plt.ylabel('Number of Variants')
        plt.title('Primary Comparisons (vs WT-CTRL)', pad=20)
        plt.xticks(x_prim, primary_treatments, rotation=45, ha='right')
        plt.legend()
        
        # Plot secondary comparisons if any
        if len(secondary_comparisons) > 0:
            ax2 = plt.subplot(1, 2, 2)
            sec_treatments = secondary_comparisons['Treatment']
            x_sec = np.arange(len(sec_treatments))
            
            # Create paired bars for treatment and its specific control
            bars2 = plt.bar(x_sec, secondary_comparisons['Treatment_Variants'], width=0.4, label='Treatment')
            ctrl_bars = plt.bar(x_sec + 0.4, secondary_comparisons['Control_Variants'], width=0.4, label='Specific Control')
            
            for i, bar in enumerate(bars2):
                treatment = sec_treatments.iloc[i]
                bar.set_color(color_map.get(treatment, '#2ecc71'))
            
            plt.xlabel('Treatment Group')
            plt.ylabel('Number of Variants')
            plt.title('Secondary Comparisons (with Specific Controls)', pad=20)
            plt.xticks(x_sec + 0.2, sec_treatments, rotation=45, ha='right')
            plt.legend()
        
        plt.tight_layout()
        plt.savefig('analysis/treatment_control_analysis/primary_vs_secondary_comparisons.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Fold Change Visualization
    plt.figure(figsize=(12, 7))
    
    # Debug Point 4: Print fold change data
    print("\nFold change data:")
    print(results_df[['Treatment', 'Control', 'Fold_Change', 'Q_Value']].to_string())
    
    # Create fold change bars with error checking and custom colors
    valid_mask = pd.notna(results_df['Fold_Change'])
    
    # Create bars with custom colors
    bars = plt.bar(np.arange(len(results_df[valid_mask])), 
                  results_df['Fold_Change'][valid_mask])
    
    # Apply treatment-specific colors
    for i, bar in enumerate(bars):
        treatment = results_df['Treatment'].iloc[i]
        bar.set_color(color_map.get(treatment, '#2ecc71'))
        
        # Add significance markers using different shapes/colors
        q_val = results_df['Q_Value'].iloc[i]
        if q_val < 0.001:
            plt.scatter(i, bar.get_height() + 1, marker='*', s=200, color='red')
        elif q_val < 0.01:
            plt.scatter(i, bar.get_height() + 1, marker='*', s=150, color='orange')
        elif q_val < 0.05:
            plt.scatter(i, bar.get_height() + 1, marker='*', s=100, color='yellow')
    
    plt.axhline(y=1, color='black', linestyle='--', alpha=0.3)
    plt.ylabel('Fold Change (Treatment/Control)')
    plt.title('Magnitude of Treatment Effect', pad=20)
    
    # Custom x-tick labels that include control information
    plt.xticks(np.arange(len(results_df[valid_mask])), 
              [f"{t}\n(vs {c})" for t, c in zip(results_df['Treatment'][valid_mask], results_df['Control'][valid_mask])], 
              rotation=45, ha='right')
    
    # Add legend
    legend_elements = [
        Patch(facecolor='red', label='p < 0.001 (***)'),
        Patch(facecolor='orange', label='p < 0.01 (**)'),
        Patch(facecolor='yellow', label='p < 0.05 (*)'),
        Patch(facecolor='gray', label='Not significant')
    ]
    plt.legend(handles=legend_elements, title='Significance Level')
    
    # Add value labels with error checking
    for i, bar in enumerate(bars):
        height = bar.get_height()
        if pd.notna(height):
            plt.text(bar.get_x() + bar.get_width()/2., height - 0.1 * max(results_df['Fold_Change']),
                    f'{height:.1f}×', ha='center', va='top', fontsize=9, 
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=1))
    
    plt.ylim(0, max(results_df['Fold_Change']) * 1.2)  # Add some space at the top
    plt.tight_layout()
    plt.savefig('analysis/treatment_control_analysis/fold_change_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Treatment Types - Biological Classification Visual
    plt.figure(figsize=(14, 8))
    
    # Create categories based on biological characteristics
    categories = {
        'Temperature Adaptation': ['WT-37', 'CAS'],
        'Low Oxygen Adaptation': ['WTA', 'STC']
    }
    
    # Prepare data for grouped bar chart
    cat_data = []
    for category, treatments_list in categories.items():
        for treatment in treatments_list:
            if treatment in results_df['Treatment'].values:
                row = results_df[results_df['Treatment'] == treatment].iloc[0]
                cat_data.append({
                    'Category': category,
                    'Treatment': treatment,
                    'Variants': row['Treatment_Variants'],
                    'Control': row['Control'],
                    'Control_Variants': row['Control_Variants'],
                    'Fold_Change': row['Fold_Change'],
                    'Has_Gene': 'Yes' if treatment in ['STC', 'CAS'] else 'No'
                })
    
    if cat_data:  # Only create this plot if we have categorized data
        df_cat = pd.DataFrame(cat_data)
        
        # Plot grouped by category
        g = sns.catplot(
            data=df_cat, kind="bar",
            x="Category", y="Variants", hue="Treatment",
            palette=color_map, height=6, aspect=1.5, legend=False
        )
        
        plt.title('Biological Classification of Treatments', pad=20)
        plt.ylabel('Number of Variants')
        
        # Add legend with treatment descriptions
        handles, labels = plt.gca().get_legend_handles_labels()
        descriptions = {t: d for t, d in zip(results_df['Treatment'], results_df['Description'])}
        legend_labels = [f"{label}: {descriptions.get(label, '')}" for label in labels]
        plt.legend(handles, legend_labels, title="Treatment", loc='best')
        
        plt.tight_layout()
        plt.savefig('analysis/treatment_control_analysis/biological_classification.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 5. Impact Visualization - Horizontal Bar Chart
    plt.figure(figsize=(14, 8))
    
    # Sort by fold change for better visualization
    sorted_df = results_df.sort_values('Fold_Change', ascending=False)
    
    # Create horizontal bars with custom colors
    bars = plt.barh(np.arange(len(sorted_df)), sorted_df['Fold_Change'])
    
    # Apply treatment-specific colors
    for i, bar in enumerate(bars):
        treatment = sorted_df['Treatment'].iloc[i]
        bar.set_color(color_map.get(treatment, '#2ecc71'))
    
    plt.axvline(x=1, color='black', linestyle='--', alpha=0.3)
    plt.xlabel('Fold Change (Treatment/Control)')
    plt.title('Impact of Treatments Ranked by Effect Size', pad=20)
    
    # Custom y-tick labels with treatment and control info
    plt.yticks(np.arange(len(sorted_df)), 
              [f"{t} vs {c}" for t, c in zip(sorted_df['Treatment'], sorted_df['Control'])])
    
    # Add value labels
    for i, bar in enumerate(bars):
        width = bar.get_width()
        if pd.notna(width):
            plt.text(width + 0.5, bar.get_y() + bar.get_height()/2,
                    f'{width:.1f}×', ha='left', va='center')
    
    plt.tight_layout()
    plt.savefig('analysis/treatment_control_analysis/impact_ranked_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 6. Description Table Plot (visualize treatment descriptions)
    plt.figure(figsize=(12, len(results_df) * 0.8))
    
    # Create a table-like visualization
    table_data = []
    for i, row in results_df.iterrows():
        table_data.append([
            row['Treatment'], 
            row['Description'], 
            row['Control'],
            f"{row['Treatment_Variants']} variants",
            f"{row['Fold_Change']:.1f}×"
        ])
    
    # Create table
    table = plt.table(
        cellText=table_data,
        colLabels=['Treatment', 'Description', 'Control', 'Variants', 'Fold Change'],
        cellLoc='center',
        loc='center'
    )
    
    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 1.5)  # Adjust row height
    
    # Hide axes
    plt.axis('off')
    plt.title('Treatment Descriptions and Key Metrics', pad=20)
    
    plt.tight_layout()
    plt.savefig('analysis/treatment_control_analysis/treatment_descriptions.png', dpi=300, bbox_inches='tight')
    plt.close()

def load_gene_level_variant_data(results_df):
    """
    Load and prepare detailed variant data for gene-specific analysis.
    
    Args:
        results_df (pandas.DataFrame): DataFrame with treatment vs control comparison results
        
    Returns:
        pandas.DataFrame: DataFrame with variant details mapped to genes, or None if data not available
    """
    variant_files = [
        'results/gene_variants/all_gene_variants.tsv',
        'results/scaffold_variants/all_scaffold_variants.tsv'
    ]
    
    # Try to find an appropriate variant file
    variant_file = None
    for file_path in variant_files:
        if os.path.exists(file_path):
            variant_file = file_path
            break
            
    if variant_file is None:
        print("No variant files found for gene-level analysis.")
        return None
        
    try:
        # Load variant data
        print(f"Loading variant data from {variant_file}...")
        variant_data = pd.read_csv(variant_file, sep='\t')
        print(f"Loaded {len(variant_data)} variants")
        
        # Check if we need to load gene mapping data
        if 'in_gene' not in variant_data.columns or 'gene_id' not in variant_data.columns:
            # Load gene mapping data
            if load_gene_mapping():
                # Map variants to genes
                variant_data = map_variants_to_genes(variant_data)
            else:
                print("Could not load gene mapping data. Gene-specific analysis will be limited.")
                return None
        
        # Add Treatment column if not already present
        if 'Treatment' not in variant_data.columns and 'Sample' in variant_data.columns:
            # Extract treatment from sample ID (assuming format like "WT-37-55-1")
            variant_data['Treatment'] = variant_data['Sample'].str.extract(r'(WTA|WT-37|STC|CAS)')[0]
        
        return variant_data
    except Exception as e:
        print(f"Error loading variant data: {e}")
        return None

if __name__ == "__main__":
    # Load results and print initial debug info
    try:
        # Load treatment vs control statistics
        results_df = pd.read_csv('analysis/treatment_control_analysis/treatment_vs_control_statistics.csv')
        print("Initial data shape:", results_df.shape)
        print("Columns:", results_df.columns.tolist())
        
        # Create standard visualizations
        create_treatment_control_visualizations(results_df)
        print("Standard visualizations created successfully!")
        
        # Try to load variant data for gene-specific analysis
        variant_data = load_gene_level_variant_data(results_df)
        
        # Create gene-specific visualizations if data is available
        if variant_data is not None and len(variant_data) > 0:
            create_gene_specific_visualizations(results_df, variant_data)
            print("Gene-specific visualizations created successfully!")
        else:
            print("Skipping gene-specific visualizations due to missing data.")
            
        print("\nAll visualizations completed successfully!")
    except Exception as e:
        print(f"Error: {e}")
        print("\nPlease run the analysis script first to generate the statistics file.")