#!/usr/bin/env python3
"""
generate_variant_report.py - Create a visual report of the scaffold variant analysis results

This script generates an HTML report with visualizations of the scaffold variant analysis
that highlight key findings about the ergosterol pathway genes.

Usage:
    python generate_variant_report.py --data_dir <scaffold_variants_directory> --output <output_html_file>
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from datetime import datetime
import base64
from io import BytesIO
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
import json
import textwrap

# Set plotting styles
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['figure.dpi'] = 100
sns.set_context("notebook", font_scale=1.2)

# Define color palettes
main_colors = ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD"]
impact_colors = {"HIGH": "#D62728", "MODERATE": "#FF7F0E", "LOW": "#2CA02C", "MODIFIER": "#1F77B4"}
variant_type_colors = {"SNV": "#1F77B4", "DELETION": "#D62728", "INSERTION": "#2CA02C", "COMPLEX": "#9467BD"}
gene_colors = plt.cm.tab20(np.linspace(0, 1, 12))

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate a visual report of scaffold variant analysis results')
    parser.add_argument('--data_dir', required=True, help='Directory containing scaffold variant analysis results')
    parser.add_argument('--output', required=True, help='Output HTML file path')
    return parser.parse_args()

def load_data(data_dir):
    """Load data from the scaffold variant analysis."""
    data = {}
    
    # Load main variant file
    data['variants'] = pd.read_csv(os.path.join(data_dir, 'treatment_specific_scaffold_variants.tsv'), sep='\t')
    
    # Load gene proximity summary
    data['gene_proximity'] = pd.read_csv(os.path.join(data_dir, 'gene_proximity_treatment_specific_summary.tsv'), sep='\t')
    
    # Load treatment comparison
    data['treatment'] = pd.read_csv(os.path.join(data_dir, 'treatment_specific_comparison.tsv'), sep='\t')
    
    # Load summary text
    with open(os.path.join(data_dir, 'treatment_specific_scaffold_variant_summary.txt'), 'r') as f:
        data['summary'] = f.read()
    
    # Extract key statistics
    data['stats'] = extract_key_stats(data['summary'])
    
    return data

def extract_key_stats(summary_text):
    """Extract key statistics from the summary text."""
    stats = {}
    
    # Extract total variants
    total_match = re.search(r'Total variants on target scaffolds: (\d+)', summary_text)
    if total_match:
        stats['total_variants'] = int(total_match.group(1))
    
    # Extract scaffold counts
    scaffold_section = re.search(r'Variants by scaffold:(.*?)Variants by distance category:', summary_text, re.DOTALL)
    if scaffold_section:
        scaffold_lines = scaffold_section.group(1).strip().split('\n')
        scaffolds = {}
        for line in scaffold_lines:
            scaffold_match = re.search(r'(\w+): (\d+) variants', line.strip())
            if scaffold_match:
                scaffolds[scaffold_match.group(1)] = int(scaffold_match.group(2))
        stats['scaffold_counts'] = scaffolds
    
    # Extract nearest gene counts
    gene_section = re.search(r'Variants by nearest gene:(.*?)Variants by treatment:', summary_text, re.DOTALL)
    if gene_section:
        gene_lines = gene_section.group(1).strip().split('\n')
        genes = {}
        for line in gene_lines:
            gene_match = re.search(r'(\w+) \((\w+)\): (\d+) \(([0-9.]+)%\)', line.strip())
            if gene_match:
                genes[gene_match.group(1)] = {
                    'erg_name': gene_match.group(2),
                    'count': int(gene_match.group(3)),
                    'percentage': float(gene_match.group(4))
                }
        stats['gene_counts'] = genes
    
    # Extract variant type counts
    type_section = re.search(r'Variants by type:(.*?)Top 10 variant effects:', summary_text, re.DOTALL)
    if type_section:
        type_lines = type_section.group(1).strip().split('\n')
        types = {}
        for line in type_lines:
            type_match = re.search(r'(\w+): (\d+) \(([0-9.]+)%\)', line.strip())
            if type_match:
                types[type_match.group(1)] = {
                    'count': int(type_match.group(2)),
                    'percentage': float(type_match.group(3))
                }
        stats['variant_types'] = types
    
    # Extract impact counts
    impact_section = re.search(r'Variants by impact:(.*?)$', summary_text, re.DOTALL)
    if impact_section:
        impact_lines = impact_section.group(1).strip().split('\n')
        impacts = {}
        for line in impact_lines:
            impact_match = re.search(r'(\w+): (\d+) \(([0-9.]+)%\)', line.strip())
            if impact_match:
                impacts[impact_match.group(1)] = {
                    'count': int(impact_match.group(2)),
                    'percentage': float(impact_match.group(3))
                }
        stats['impacts'] = impacts
    
    return stats

def create_gene_variant_proximity_chart(gene_proximity_df):
    """Create a horizontal bar chart showing variant proximity to genes."""
    # Prepare data
    gene_df = gene_proximity_df.copy()
    gene_df['Total_5kb'] = gene_df['Variants_Within'] + gene_df['Variants_Upstream_1kb'] + \
                           gene_df['Variants_Downstream_1kb'] + gene_df['Variants_Upstream_5kb'] + \
                           gene_df['Variants_Downstream_5kb']
    
    # Sort by total variants
    gene_df = gene_df.sort_values('Total_5kb', ascending=False)
    
    # Set up the figure
    fig = plt.figure(figsize=(12, 8))
    
    # Create horizontal stacked bar chart
    categories = ['Variants_Within', 'Variants_Upstream_1kb', 'Variants_Downstream_1kb', 
                  'Variants_Upstream_5kb', 'Variants_Downstream_5kb']
    
    category_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    category_labels = ['Within Gene', 'Upstream (<1kb)', 'Downstream (<1kb)', 
                       'Upstream (1-5kb)', 'Downstream (1-5kb)']
    
    bars = []
    bottom = np.zeros(len(gene_df))
    
    for i, (category, color) in enumerate(zip(categories, category_colors)):
        values = gene_df[category].values
        bars.append(plt.barh(gene_df['SC_Gene_ID'] + ' (' + gene_df['ERG_Name'] + ')', 
                            values, left=bottom, color=color, height=0.7))
        bottom += values
    
    # Add labels and title
    plt.xlabel('Number of Variants')
    plt.title('Variants in Proximity to Ergosterol Pathway Genes', fontsize=14)
    
    # Add legend
    plt.legend(bars, category_labels, loc='upper right', bbox_to_anchor=(1.1, 1))
    
    # Adjust layout
    plt.tight_layout()
    
    # Save to BytesIO
    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight')
    buffer.seek(0)
    plt.close()
    
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def create_pathway_diagram(gene_proximity_df, variants_df):
    """Create an ergosterol pathway diagram with variant information overlay."""
    # Define the ergosterol pathway steps
    pathway_steps = [
        {"gene": "ERG9", "step": 1, "product": "Squalene"},
        {"gene": "ERG1", "step": 2, "product": "Squalene epoxide"},
        {"gene": "ERG7", "step": 3, "product": "Lanosterol"},
        {"gene": "ERG11", "step": 4, "product": "4,4-dimethylcholesta-8,14,24-trienol"},
        {"gene": "ERG24", "step": 5, "product": "4,4-dimethylzymosterol"},
        {"gene": "ERG25", "step": 6, "product": "4-methylzymosterol"},
        {"gene": "ERG6", "step": 7, "product": "Zymosterol"},
        {"gene": "ERG2", "step": 8, "product": "Fecosterol"},
        {"gene": "ERG3", "step": 9, "product": "Episterol"},
        {"gene": "ERG5", "step": 10, "product": "Ergosta-5,7,24-trienol"},
        {"gene": "ERG4", "step": 11, "product": "Ergosterol"}
    ]
    
    # Prepare gene data
    gene_data = {}
    for i, row in gene_proximity_df.iterrows():
        gene_data[row['ERG_Name']] = {
            'sc_id': row['SC_Gene_ID'],
            'total_nearby': row['Total_Within_5kb'],
            'within': row['Variants_Within']
        }
    
    # Calculate variant density by gene
    gene_variant_counts = variants_df.groupby('Nearest_Gene_ERG').size().to_dict()
    
    # Create the figure
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111)
    
    # Define coordinates for pathway steps (vertical layout)
    num_steps = len(pathway_steps)
    step_y = np.linspace(0.1, 0.9, num_steps)
    step_x = 0.5
    box_width = 0.3
    box_height = 0.06
    arrow_length = (step_y[1] - step_y[0]) * 0.6
    
    # Draw pathway
    for i, step in enumerate(pathway_steps):
        gene = step["gene"]
        product = step["product"]
        
        # Determine box color based on variant density
        variant_count = gene_data.get(gene, {}).get('total_nearby', 0)
        
        # Normalize variant count (0-1 scale)
        max_variants = max([g.get('total_nearby', 0) for g in gene_data.values()])
        if max_variants > 0:
            color_intensity = min(variant_count / max_variants, 1.0)
        else:
            color_intensity = 0
            
        # Create a gradient from white to red
        box_color = (1, 1-color_intensity, 1-color_intensity)
        
        # Draw gene box
        gene_box = plt.Rectangle((step_x - box_width/2, step_y[i] - box_height/2), 
                                box_width, box_height, 
                                facecolor=box_color, 
                                edgecolor='black', 
                                alpha=0.8)
        ax.add_patch(gene_box)
        
        # Add gene name and variant count
        plt.text(step_x, step_y[i], f"{gene}\n({gene_data.get(gene, {}).get('sc_id', '')})",
                 ha='center', va='center', fontsize=10, fontweight='bold')
        
        # Add product name
        plt.text(step_x + box_width/2 + 0.02, step_y[i], product,
                 ha='left', va='center', fontsize=9, fontweight='normal')
        
        # Add variant count
        plt.text(step_x - box_width/2 - 0.02, step_y[i], 
                 f"Variants: {variant_count}",
                 ha='right', va='center', fontsize=9, 
                 color='darkred' if variant_count > 0 else 'black')
        
        # Draw arrow to next step
        if i < num_steps - 1:
            ax.arrow(step_x, step_y[i] + box_height/2 + 0.01, 
                    0, arrow_length - 0.02,
                    head_width=0.02, head_length=0.02, 
                    fc='black', ec='black')
    
    # Set plot limits and remove axes
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.axis('off')
    
    # Add title and legend
    plt.title('Ergosterol Biosynthesis Pathway\nVariant Density by Gene', fontsize=16, pad=20)
    
    # Create color gradient for legend
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    
    # Create a separate axis for the color bar
    cax = fig.add_axes([0.2, 0.05, 0.6, 0.03])
    cax.imshow(gradient, aspect='auto', cmap=plt.cm.Reds_r)
    cax.set_title('Variant Density')
    cax.set_xticks([0, 255])
    cax.set_xticklabels(['Low', 'High'])
    cax.set_yticks([])
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])
    
    # Save to BytesIO
    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight')
    buffer.seek(0)
    plt.close()
    
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def create_chromosome_variant_map(variants_df, gene_proximity_df):
    """Create a map of variant distribution along chromosomes/scaffolds."""
    # Get unique scaffolds
    scaffolds = variants_df['Scaffold'].unique()
    
    # Initialize plot
    fig, axs = plt.subplots(len(scaffolds), 1, figsize=(12, 3*len(scaffolds)), sharex=False)
    
    # Process each scaffold
    for i, scaffold in enumerate(sorted(scaffolds)):
        # Get axis for this scaffold
        ax = axs[i] if len(scaffolds) > 1 else axs
        
        # Filter variants for this scaffold
        scaffold_variants = variants_df[variants_df['Scaffold'] == scaffold]
        
        # Get max position for this scaffold
        max_pos = scaffold_variants['Position'].max()
        
        # Create histogram
        ax.hist(scaffold_variants['Position'], bins=50, alpha=0.7, color='#1f77b4')
        
        # Find genes on this scaffold
        scaffold_genes = gene_proximity_df[gene_proximity_df['SC_Gene_ID'].str.contains(
            '|'.join(scaffold_variants['Nearest_Gene_SC_ID'].unique()))]
        
        # Plot gene positions
        for _, gene in scaffold_genes.iterrows():
            # Find gene position from variants
            gene_variants = scaffold_variants[scaffold_variants['Nearest_Gene_SC_ID'] == gene['SC_Gene_ID']]
            if len(gene_variants) > 0:
                # Use the position of variants closest to the gene as an approximation
                gene_pos = gene_variants.loc[gene_variants['Distance'].idxmin(), 'Position']
                
                # Add gene label
                ax.axvline(x=gene_pos, color='red', linestyle='--', alpha=0.6)
                ax.text(gene_pos, ax.get_ylim()[1]*0.9, f"{gene['SC_Gene_ID']}\n({gene['ERG_Name']})", 
                        rotation=90, ha='right', va='top', fontsize=9)
        
        # Add scaffold info
        ax.set_title(f"Variant Distribution on {scaffold}")
        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("Variant Count")
        
        # Add grid
        ax.grid(True, alpha=0.3)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save to BytesIO
    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight')
    buffer.seek(0)
    plt.close()
    
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def create_variant_type_impact_plot(variants_df):
    """Create visualizations of variant types and impacts."""
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot variant types
    variant_counts = variants_df['Variant_Type'].value_counts()
    explode = [0.1] * len(variant_counts)
    ax1.pie(variant_counts, labels=variant_counts.index, autopct='%1.1f%%',
           explode=explode, colors=[variant_type_colors.get(t, '#CCCCCC') for t in variant_counts.index],
           wedgeprops={'edgecolor': 'w', 'linewidth': 1.5})
    ax1.set_title('Variant Types', fontsize=14)
    
    # Plot variant impacts
    impact_counts = variants_df['Impact'].value_counts()
    explode = [0.1] * len(impact_counts)
    ax2.pie(impact_counts, labels=impact_counts.index, autopct='%1.1f%%',
           explode=explode, colors=[impact_colors.get(i, '#CCCCCC') for i in impact_counts.index],
           wedgeprops={'edgecolor': 'w', 'linewidth': 1.5})
    ax2.set_title('Variant Impacts', fontsize=14)
    
    # Add overall title
    plt.suptitle('Variant Types and Impacts', fontsize=16, y=1.05)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save to BytesIO
    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight')
    buffer.seek(0)
    plt.close()
    
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def create_treatment_comparison_plot(treatment_df):
    """Create visualizations comparing treatments."""
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    # Filter out control treatments for cleaner visualization
    non_control_df = treatment_df[~treatment_df['Treatment'].str.contains('CTRL')]
    
    # Plot total variants by treatment
    colors = [main_colors[i % len(main_colors)] for i in range(len(treatment_df))]
    ax1.bar(treatment_df['Treatment'], treatment_df['Total_Variants'], color=colors)
    ax1.set_title('Total Variants by Treatment', fontsize=14)
    ax1.set_ylabel('Number of Variants')
    ax1.set_xticklabels(treatment_df['Treatment'], rotation=45, ha='right')
    ax1.grid(True, axis='y', alpha=0.3)
    
    # Plot stacked bar of variant proximity by treatment
    categories = ['Variants_Within_Genes', 'Variants_Within_1kb', 'Variants_Within_5kb']
    bottoms = np.zeros(len(non_control_df))
    
    category_colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    category_labels = ['Within Genes', 'Within 1kb', '1-5kb from Genes']
    
    for i, (category, color) in enumerate(zip(categories, category_colors)):
        values = non_control_df[category].values
        ax2.bar(non_control_df['Treatment'], values, bottom=bottoms, 
                label=category_labels[i], color=color)
        bottoms += values
    
    ax2.set_title('Variant Proximity to Genes by Treatment', fontsize=14)
    ax2.set_ylabel('Number of Variants')
    ax2.legend()
    ax2.grid(True, axis='y', alpha=0.3)
    
    # Add overall title
    plt.suptitle('Treatment Comparison', fontsize=16, y=1.05)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save to BytesIO
    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight')
    buffer.seek(0)
    plt.close()
    
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def create_gene_variant_heatmap(variants_df):
    """Create a heatmap showing variant counts by gene and treatment."""
    # Filter out variants that are too far from genes
    close_variants = variants_df[variants_df['Distance'] <= 5000]
    
    # Create pivot table
    pivot = pd.pivot_table(
        close_variants,
        values='Position',
        index='Nearest_Gene_SC_ID',
        columns='Treatment',
        aggfunc='count',
        fill_value=0
    )
    
    # Add ERG names to index
    erg_names = {}
    for _, row in close_variants.iterrows():
        if pd.notna(row['Nearest_Gene_SC_ID']) and pd.notna(row['Nearest_Gene_ERG']):
            erg_names[row['Nearest_Gene_SC_ID']] = row['Nearest_Gene_ERG']
    
    new_index = [f"{idx} ({erg_names.get(idx, '')})" for idx in pivot.index]
    pivot.index = new_index
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Create heatmap
    sns.heatmap(pivot, annot=True, cmap='YlOrRd', fmt='g')
    
    # Add title and labels
    plt.title('Variants Within 5kb by Gene and Treatment', fontsize=14)
    plt.ylabel('Gene')
    plt.xlabel('Treatment')
    
    # Rotate x-axis labels
    plt.xticks(rotation=45, ha='right')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save to BytesIO
    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight')
    buffer.seek(0)
    plt.close()
    
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def create_distance_distribution_plot(variants_df):
    """Create a visualization of variant distance distribution."""
    # Create a copy with capped distances for better visualization
    plot_df = variants_df.copy()
    max_dist_to_show = 10000
    plot_df.loc[plot_df['Distance'] > max_dist_to_show, 'Distance'] = max_dist_to_show
    
    # Create figure
    plt.figure(figsize=(12, 6))
    
    # Create histogram
    plt.hist(plot_df['Distance'], bins=50, alpha=0.7, color='#1f77b4')
    
    # Add vertical lines for important distances
    plt.axvline(x=0, color='red', linestyle='--', label='Gene Boundary')
    plt.axvline(x=1000, color='orange', linestyle='--', label='1kb')
    plt.axvline(x=5000, color='green', linestyle='--', label='5kb')
    
    # Add title and labels
    plt.title('Distribution of Variants by Distance to Nearest Gene', fontsize=14)
    plt.xlabel('Distance (bp)')
    plt.ylabel('Number of Variants')
    plt.legend()
    
    # Add grid
    plt.grid(True, alpha=0.3)
    
    # Add text annotation for capped values
    plt.text(max_dist_to_show*0.95, plt.gca().get_ylim()[1]*0.9, 
             f"Distances > {max_dist_to_show}bp\ntruncated for visualization",
             ha='right', va='top', bbox=dict(facecolor='white', alpha=0.7))
    
    # Adjust layout
    plt.tight_layout()
    
    # Save to BytesIO
    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight')
    buffer.seek(0)
    plt.close()
    
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def create_variant_effect_plot(variants_df):
    """Create visualization of top variant effects."""
    # Get top 10 effects
    effect_counts = variants_df['Effect'].value_counts().head(10)
    
    # Create figure
    plt.figure(figsize=(12, 6))
    
    # Create horizontal bar chart
    bars = plt.barh(effect_counts.index, effect_counts.values, color=main_colors[0])
    
    # Add value labels
    for bar in bars:
        width = bar.get_width()
        plt.text(width * 1.01, bar.get_y() + bar.get_height()/2, 
                f"{width} ({width/sum(effect_counts.values)*100:.1f}%)",
                va='center')
    
    # Add title and labels
    plt.title('Top 10 Variant Effects', fontsize=14)
    plt.xlabel('Number of Variants')
    
    # Add grid
    plt.grid(True, axis='x', alpha=0.3)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save to BytesIO
    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight')
    buffer.seek(0)
    plt.close()
    
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def create_key_findings_infographic():
    """Create an infographic highlighting key findings about ergosterol pathway conservation."""
    # Create figure
    fig = plt.figure(figsize=(12, 8))
    
    # Set background color
    fig.patch.set_facecolor('#f8f9fa')
    
    # Add title
    plt.text(0.5, 0.95, 'Key Findings: Ergosterol Pathway Conservation', 
             ha='center', va='top', fontsize=18, fontweight='bold')
    
    # Add divider
    plt.axhline(y=0.9, xmin=0.05, xmax=0.95, color='#343a40', alpha=0.5)
    
    # Add key finding boxes
    findings = [
        {
            'title': 'Pathway Conservation',
            'text': 'Most ergosterol pathway genes show strong conservation with few variants',
            'icon': '🧬',
            'position': (0.25, 0.75)
        },
        {
            'title': 'Regulatory Variation',
            'text': 'When variants occur, they are predominantly in regulatory regions',
            'icon': '🔄',
            'position': (0.75, 0.75)
        },
        {
            'title': 'Gene-Specific Tolerance',
            'text': 'ERG25 shows more tolerance to nearby variants than other genes',
            'icon': '⚖️',
            'position': (0.25, 0.45)
        },
        {
            'title': 'Adaptation Patterns',
            'text': 'Similar variant patterns across treatments suggest background variation',
            'icon': '🧪',
            'position': (0.75, 0.45)
        },
        {
            'title': 'Insertions Dominant',
            'text': 'Insertions (46%) more common than SNVs (26%) or deletions (28%)',
            'icon': '➕',
            'position': (0.5, 0.15)
        }
    ]
    
    for finding in findings:
        # Draw box
        box = mpatches.FancyBboxPatch((finding['position'][0] - 0.2, finding['position'][1] - 0.1), 
                           0.4, 0.2, facecolor='white', edgecolor='#343a40', alpha=0.8,
                           boxstyle='round,pad=0.5')
        plt.gca().add_patch(box)
        
        # Add icon
        plt.text(finding['position'][0] - 0.18, finding['position'][1] + 0.08, 
                finding['icon'], fontsize=20)
        
        # Add title
        plt.text(finding['position'][0] - 0.1, finding['position'][1] + 0.08, 
                finding['title'], fontweight='bold', fontsize=11)
        # Add text
        wrapped_text = textwrap.fill(finding['text'], width=30)  # 30 characters per line
        plt.text(finding['position'][0] - 0.18, finding['position'][1] - 0.02, 
                wrapped_text, fontsize=9)
    
    # Set plot limits and remove axes
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.axis('off')
    
    # Save to BytesIO
    buffer = BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight')
    buffer.seek(0)
    plt.close()
    
    return base64.b64encode(buffer.getvalue()).decode('utf-8')

def create_interactive_variant_map(variants_df, gene_proximity_df):
    """Create an interactive variant map using Plotly."""
    # Convert distance to a bounded value for visualization
    variants_df['Display_Distance'] = variants_df['Distance'].copy()
    max_dist = 50000
    variants_df.loc[variants_df['Display_Distance'] > max_dist, 'Display_Distance'] = max_dist
    
    # Create color mapping based on distance
    variants_df['Color'] = pd.cut(
        variants_df['Display_Distance'],
        bins=[0, 500, 1000, 5000, 10000, max_dist],
        labels=['Within 500bp', '500bp-1kb', '1kb-5kb', '5kb-10kb', '10kb+']
    )
    
    # Create color scale
    color_scale = {
        'Within 500bp': 'rgba(255, 0, 0, 0.8)',
        '500bp-1kb': 'rgba(255, 127, 14, 0.8)',
        '1kb-5kb': 'rgba(44, 160, 44, 0.8)',
        '5kb-10kb': 'rgba(31, 119, 180, 0.8)',
        '10kb+': 'rgba(128, 128, 128, 0.5)'
    }
    
    # Create figures for each scaffold
    figures = []
    
    for scaffold in sorted(variants_df['Scaffold'].unique()):
        # Filter variants for this scaffold
        scaffold_variants = variants_df[variants_df['Scaffold'] == scaffold].copy()
        
        # Create figure
        fig = go.Figure()
        
        # Add scatter plot for variants
        for category in color_scale:
            category_variants = scaffold_variants[scaffold_variants['Color'] == category]
            
            if not category_variants.empty:
                fig.add_trace(go.Scatter(
                    x=category_variants['Position'],
                    y=[1] * len(category_variants),
                    mode='markers',
                    marker=dict(
                        size=8,
                        color=color_scale[category],
                        line=dict(width=1, color='white')
                    ),
                    name=category,
                    text=category_variants.apply(
                        lambda row: f"Position: {row['Position']}<br>" + 
                                   f"Variant: {row['Ref']} > {row['Alt']}<br>" +
                                   f"Effect: {row['Effect']}<br>" +
                                   f"Distance to gene: {row['Distance']}bp<br>" +
                                   f"Nearest gene: {row['Nearest_Gene_SC_ID']} ({row['Nearest_Gene_ERG']})",
                        axis=1
                    ),
                    hoverinfo='text'
                ))
        
        # Add gene positions
        genes_on_scaffold = []
        for gene_id, gene_erg in zip(
            scaffold_variants['Nearest_Gene_SC_ID'].dropna().unique(),
            scaffold_variants['Nearest_Gene_ERG'].dropna().unique()
        ):
            # Filter variants for this gene to approximate position
            gene_variants = scaffold_variants[
                (scaffold_variants['Nearest_Gene_SC_ID'] == gene_id) & 
                (scaffold_variants['Distance'] < 5000)
            ]
            
            if not gene_variants.empty:
                # Use median position of nearby variants as approximation
                gene_pos = gene_variants['Position'].median()
                genes_on_scaffold.append((gene_id, gene_erg, gene_pos))
        
        # Add gene annotations
        for gene_id, gene_erg, pos in genes_on_scaffold:
            fig.add_annotation(
                x=pos,
                y=1,
                text=f"{gene_id}<br>({gene_erg})",
                showarrow=True,
                arrowhead=2,
                arrowsize=1,
                arrowwidth=2,
                arrowcolor="#636363",
                ax=0,
                ay=-40,
                bgcolor="white",
                opacity=0.8
            )
        
        # Update layout
        fig.update_layout(
            title=f"Variant Distribution on {scaffold}",
            xaxis_title="Position (bp)",
            yaxis_visible=False,
            showlegend=True,
            height=400,
            margin=dict(l=50, r=50, t=50, b=50),
            plot_bgcolor='white',
            hovermode='closest'
        )
        
        figures.append(fig)
    
    # Convert to HTML
    html_outputs = []
    for fig in figures:
        html_outputs.append(fig.to_html(full_html=False, include_plotlyjs='cdn'))
    
    return html_outputs

def create_interactive_gene_heatmap(gene_proximity_df):
    """Create an interactive gene heatmap using Plotly."""
    # Create figure
    fig = go.Figure()
    
    # Prepare data
    df = gene_proximity_df.copy()
    
    # Calculate total
    df['Total_5kb'] = df['Variants_Within'] + df['Variants_Upstream_1kb'] + \
                      df['Variants_Downstream_1kb'] + df['Variants_Upstream_5kb'] + \
                      df['Variants_Downstream_5kb']
    
    # Sort by total
    df = df.sort_values('Total_5kb', ascending=False)
    
    # Create heatmap data
    categories = ['Variants_Within', 'Variants_Upstream_1kb', 'Variants_Downstream_1kb', 
                  'Variants_Upstream_5kb', 'Variants_Downstream_5kb']
    
    category_names = ['Within Gene', 'Upstream (<1kb)', 'Downstream (<1kb)', 
                      'Upstream (1-5kb)', 'Downstream (1-5kb)']
    
    z_data = []
    for category in categories:
        z_data.append(df[category].values)
    
    # Create heatmap
    fig.add_trace(go.Heatmap(
        z=z_data,
        x=df['SC_Gene_ID'] + ' (' + df['ERG_Name'] + ')',
        y=category_names,
        colorscale='YlOrRd',
        showscale=True,
        text=[[f"{val} variants" for val in row] for row in z_data],
        hoverinfo='text+y+x'
    ))
    
    # Update layout
    fig.update_layout(
        title='Variant Distribution by Gene and Position',
        xaxis_title='Gene',
        yaxis_title='Position Relative to Gene',
        height=500,
        margin=dict(l=50, r=50, t=50, b=100),
        xaxis=dict(
            tickangle=-45
        )
    )
    
    # Convert to HTML
    html_output = fig.to_html(full_html=False, include_plotlyjs='cdn')
    
    return html_output

def generate_html_report(data, output_file):
    """Generate an HTML report with all visualizations."""
    # Create image data
    images = {}
    
    print("Generating key findings infographic...")
    images['key_findings'] = create_key_findings_infographic()
    
    print("Generating gene variant proximity chart...")
    images['gene_proximity'] = create_gene_variant_proximity_chart(data['gene_proximity'])
    
    print("Generating pathway diagram...")
    images['pathway'] = create_pathway_diagram(data['gene_proximity'], data['variants'])
    
    print("Generating chromosome variant map...")
    images['chromosome_map'] = create_chromosome_variant_map(data['variants'], data['gene_proximity'])
    
    print("Generating variant type and impact plot...")
    images['variant_type_impact'] = create_variant_type_impact_plot(data['variants'])
    
    print("Generating treatment comparison plot...")
    images['treatment_comparison'] = create_treatment_comparison_plot(data['treatment'])
    
    print("Generating gene variant heatmap...")
    images['gene_heatmap'] = create_gene_variant_heatmap(data['variants'])
    
    print("Generating distance distribution plot...")
    images['distance_distribution'] = create_distance_distribution_plot(data['variants'])
    
    print("Generating variant effect plot...")
    images['variant_effect'] = create_variant_effect_plot(data['variants'])
    
    print("Generating interactive visualizations...")
    interactive_maps = create_interactive_variant_map(data['variants'], data['gene_proximity'])
    interactive_heatmap = create_interactive_gene_heatmap(data['gene_proximity'])
    
    # HTML template
    html_template = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Ergosterol Pathway Variant Analysis Report</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                line-height: 1.6;
                color: #333;
                margin: 0;
                padding: 0;
                background-color: #f8f9fa;
            }
            header {
                background-color: #343a40;
                color: white;
                padding: 1rem;
                text-align: center;
            }
            .container {
                max-width: 1200px;
                margin: 0 auto;
                padding: 1rem;
            }
            .section {
                background-color: white;
                margin-bottom: 2rem;
                padding: 1.5rem;
                border-radius: 5px;
                box-shadow: 0 1px 3px rgba(0,0,0,0.12), 0 1px 2px rgba(0,0,0,0.24);
            }
            h1, h2, h3 {
                color: #343a40;
            }
            h2 {
                border-bottom: 2px solid #f8f9fa;
                padding-bottom: 0.5rem;
                margin-top: 0;
            }
            .image-container {
                text-align: center;
                margin: 1.5rem 0;
            }
            img {
                max-width: 100%;
                height: auto;
                border-radius: 5px;
            }
            .stats {
                display: flex;
                flex-wrap: wrap;
                justify-content: space-around;
                margin: 1rem 0;
            }
            .stat-box {
                background-color: #f8f9fa;
                border-radius: 5px;
                padding: 1rem;
                margin: 0.5rem;
                min-width: 200px;
                text-align: center;
            }
            .stat-box h3 {
                margin-top: 0;
                font-size: 1.2rem;
            }
            .stat-box p {
                font-size: 2rem;
                font-weight: bold;
                margin: 0.5rem 0;
                color: #0275d8;
            }
            .stat-box small {
                font-size: 0.9rem;
                color: #6c757d;
            }
            footer {
                background-color: #343a40;
                color: white;
                text-align: center;
                padding: 1rem;
                margin-top: 2rem;
            }
            .interactive-section {
                margin: 2rem 0;
            }
            .interactive-tabs {
                display: flex;
                border-bottom: 1px solid #dee2e6;
            }
            .tab {
                padding: 0.5rem 1rem;
                cursor: pointer;
                border: 1px solid transparent;
                border-top-left-radius: 0.25rem;
                border-top-right-radius: 0.25rem;
            }
            .tab.active {
                color: #495057;
                background-color: #fff;
                border-color: #dee2e6 #dee2e6 #fff;
            }
            .tab-content {
                display: none;
                padding: 1rem;
                border: 1px solid #dee2e6;
                border-top: none;
                background-color: #fff;
            }
            .tab-content.active {
                display: block;
            }
        </style>
    </head>
    <body>
        <header>
            <h1>Ergosterol Pathway Variant Analysis Report</h1>
            <p>Analysis of variants in ergosterol biosynthesis genes across yeast strains</p>
        </header>
        
        <div class="container">
            <!-- Key Findings Section -->
            <div class="section">
                <h2>Key Findings</h2>
                <div class="image-container">
                    <img src="data:image/png;base64,${key_findings}" alt="Key Findings Infographic">
                </div>
            </div>
            
            <!-- High-Level Statistics -->
            <div class="section">
                <h2>Project Overview</h2>
                <div class="stats">
                    <div class="stat-box">
                        <h3>Total Variants</h3>
                        <p>${total_variants}</p>
                        <small>Across all target scaffolds</small>
                    </div>
                    <div class="stat-box">
                        <h3>Target Genes</h3>
                        <p>11</p>
                        <small>Ergosterol pathway genes</small>
                    </div>
                    <div class="stat-box">
                        <h3>Close Variants</h3>
                        <p>${close_variants}</p>
                        <small>Within 5kb of target genes</small>
                    </div>
                </div>
            </div>
            
            <!-- Pathway Analysis -->
            <div class="section">
                <h2>Ergosterol Pathway Analysis</h2>
                <p>The ergosterol biosynthesis pathway shows strong conservation across samples.</p>
                <div class="image-container">
                    <img src="data:image/png;base64,${pathway}" alt="Ergosterol Pathway Diagram">
                </div>
            </div>
            
            <!-- Gene Proximity Analysis -->
            <div class="section">
                <h2>Gene Variant Proximity Analysis</h2>
                <p>Distribution of variants in proximity to ergosterol pathway genes.</p>
                <div class="image-container">
                    <img src="data:image/png;base64,${gene_proximity}" alt="Gene Variant Proximity Chart">
                </div>
                <div class="image-container">
                    <img src="data:image/png;base64,${gene_heatmap}" alt="Gene Variant Heatmap">
                </div>
            </div>
            
            <!-- Variant Characteristics -->
            <div class="section">
                <h2>Variant Characteristics</h2>
                <div class="image-container">
                    <img src="data:image/png;base64,${variant_type_impact}" alt="Variant Type and Impact">
                </div>
                <div class="image-container">
                    <img src="data:image/png;base64,${variant_effect}" alt="Variant Effect Distribution">
                </div>
            </div>
            
            <!-- Distance Distribution -->
            <div class="section">
                <h2>Variant Distance Distribution</h2>
                <div class="image-container">
                    <img src="data:image/png;base64,${distance_distribution}" alt="Distance Distribution">
                </div>
            </div>
            
            <!-- Treatment Comparison -->
            <div class="section">
                <h2>Treatment Comparison</h2>
                <div class="image-container">
                    <img src="data:image/png;base64,${treatment_comparison}" alt="Treatment Comparison">
                </div>
            </div>
            
            <!-- Chromosome Maps -->
            <div class="section">
                <h2>Chromosome Variant Maps</h2>
                <div class="image-container">
                    <img src="data:image/png;base64,${chromosome_map}" alt="Chromosome Variant Map">
                </div>
            </div>
            
            <!-- Interactive Visualizations -->
            <div class="section">
                <h2>Interactive Gene Proximity Analysis</h2>
                ${interactive_heatmap}
            </div>
            
            <div class="section">
                <h2>Interactive Variant Maps</h2>
                
                <div class="interactive-tabs">
                    ${tab_buttons}
                </div>
                
                ${tab_contents}
            </div>
        </div>
        
        <footer>
            <p>Generated on ${date}</p>
        </footer>
        
        <script>
            function openTab(evt, tabName) {
                var i, tabContent, tabLinks;
                
                // Hide all tab content
                tabContent = document.getElementsByClassName("tab-content");
                for (i = 0; i < tabContent.length; i++) {
                    tabContent[i].style.display = "none";
                }
                
                // Remove "active" class from all tab buttons
                tabLinks = document.getElementsByClassName("tab");
                for (i = 0; i < tabLinks.length; i++) {
                    tabLinks[i].className = tabLinks[i].className.replace(" active", "");
                }
                
                // Show the current tab and add "active" class to the button
                document.getElementById(tabName).style.display = "block";
                evt.currentTarget.className += " active";
            }
            
            // Set the first tab as active by default
            document.addEventListener("DOMContentLoaded", function() {
                const firstTab = document.querySelector(".tab");
                const firstTabContent = document.querySelector(".tab-content");
                
                if (firstTab && firstTabContent) {
                    firstTab.classList.add("active");
                    firstTabContent.style.display = "block";
                }
            });
        </script>
    </body>
    </html>
    """
    
    # Create tab buttons and contents
    tab_buttons = ""
    tab_contents = ""
    
    for i, scaffold_map in enumerate(interactive_maps):
        scaffold_name = f"scaffold_{i+1}"
        active_class = " active" if i == 0 else ""
        
        tab_buttons += f'<div class="tab{active_class}" onclick="openTab(event, \'{scaffold_name}\')">{scaffold_name}</div>\n'
        tab_contents += f'<div id="{scaffold_name}" class="tab-content{active_class}">{scaffold_map}</div>\n'
    
    # Calculate close variant count
    close_variants = len(data['variants'][data['variants']['Distance'] <= 5000])
    
    # Fill in template
    html_content = html_template.replace("${key_findings}", images['key_findings']) \
                              .replace("${gene_proximity}", images['gene_proximity']) \
                              .replace("${pathway}", images['pathway']) \
                              .replace("${chromosome_map}", images['chromosome_map']) \
                              .replace("${variant_type_impact}", images['variant_type_impact']) \
                              .replace("${treatment_comparison}", images['treatment_comparison']) \
                              .replace("${gene_heatmap}", images['gene_heatmap']) \
                              .replace("${distance_distribution}", images['distance_distribution']) \
                              .replace("${variant_effect}", images['variant_effect']) \
                              .replace("${interactive_heatmap}", interactive_heatmap) \
                              .replace("${tab_buttons}", tab_buttons) \
                              .replace("${tab_contents}", tab_contents) \
                              .replace("${total_variants}", str(len(data['variants']))) \
                              .replace("${close_variants}", str(close_variants)) \
                              .replace("${date}", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"Report generated: {output_file}")

def main():
    """Main function to generate the report."""
    args = parse_arguments()
    
    # Load data
    print(f"Loading data from {args.data_dir}")
    data = load_data(args.data_dir)
    
    # Generate report
    print(f"Generating report at {args.output}")
    generate_html_report(data, args.output)

if __name__ == "__main__":
    import re  # Import at top level for extract_key_stats function
    main()