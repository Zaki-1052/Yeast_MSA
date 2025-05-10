#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/satellite_genes/satellite_variant_profiling.py

"""
Script to analyze variant patterns in satellite genes.
This script analyzes the distribution, types, and impacts of variants 
in satellite genes across different treatment conditions.
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
from scipy import stats
import matplotlib.patches as mpatches

# Add the parent directory to the path to access shared modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.tools import ensure_dir, save_tsv, load_tsv, setup_plotting_style, save_plot

# Define treatment groups
TREATMENTS = {
    'WT-37': 'Temperature-adapted wild type',
    'WTA': 'Low oxygen-adapted wild type',
    'STC': 'STC gene-modified with low oxygen adaptation',
    'CAS': 'CAS gene-modified with temperature adaptation'
}

# Group treatments by adaptation type
ADAPTATION_GROUPS = {
    'Temperature': ['WT-37', 'CAS'],
    'Low Oxygen': ['WTA', 'STC']
}

# Define variant impact categories
IMPACT_CATEGORIES = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']

# Define variant effect categories
EFFECT_CATEGORIES = [
    'missense_variant', 'synonymous_variant', 'upstream_gene_variant', 
    'downstream_gene_variant', 'intergenic_variant', 'frameshift_variant',
    'stop_gained', 'stop_lost', 'start_lost', 'inframe_insertion',
    'inframe_deletion', 'splice_region_variant'
]

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Analyze variant patterns in satellite genes")
    parser.add_argument("--satellite-genes", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/satellite_genes/satellite_genes_annotated.tsv",
                      help="Path to annotated satellite genes TSV file")
    parser.add_argument("--variant-dir", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/gene_variants_expanded",
                      help="Directory containing variant data")
    parser.add_argument("--output-dir", default="/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/satellite_genes",
                      help="Directory to store output files")
    return parser.parse_args()

def load_satellite_genes(satellite_file):
    """Load satellite gene data from annotation step"""
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

def load_variants(variant_dir):
    """Load variant data for analysis"""
    variants_file = os.path.join(variant_dir, "all_gene_variants.tsv")
    
    if not os.path.exists(variants_file):
        print(f"ERROR: Variant file not found: {variants_file}")
        sys.exit(1)
    
    try:
        variants_df = pd.read_csv(variants_file, sep="\t")
        print(f"Loaded {len(variants_df)} variants from file")

        # Display information about the variants data
        print(f"Columns in variants file: {variants_df.columns.tolist()}")

        # Show sample of the data
        if not variants_df.empty:
            print(f"Sample row: {variants_df.iloc[0].to_dict()}")

            # List unique scaffolds in the variants data
            scaffolds = variants_df['Scaffold'].unique()
            print(f"Unique scaffolds in variants: {scaffolds}")

            # Count variants per scaffold
            scaffold_counts = variants_df['Scaffold'].value_counts()
            print("Variants per scaffold:")
            for scaffold, count in scaffold_counts.items():
                print(f"  {scaffold}: {count} variants")

        return variants_df
    except Exception as e:
        print(f"ERROR: Failed to load variants file: {e}")
        sys.exit(1)

def identify_satellite_variants(satellite_df, variants_df):
    """Identify variants associated with satellite genes"""
    # Create a DataFrame to store satellite gene variants
    satellite_variants = []

    # Get all unique scaffolds in the variants dataframe
    variant_scaffolds = set(variants_df['Scaffold'].unique())

    # For each satellite gene, find associated variants
    for _, satellite in satellite_df.iterrows():
        # Extract satellite gene information
        satellite_id = satellite.get('satellite_gene_id', '')
        satellite_name = satellite.get('satellite_gene_name', '')
        scaffold = satellite.get('scaffold', '')
        start = satellite.get('satellite_start', 0)
        end = satellite.get('satellite_end', 0)
        erg_gene_id = satellite.get('erg_gene_id', '')
        erg_gene_name = satellite.get('erg_gene_name', '')
        erg_start = satellite.get('erg_start', 0)
        erg_end = satellite.get('erg_end', 0)
        function_category = satellite.get('function_category', 'Unknown')

        # Skip if we don't have a valid scaffold or if it's not in variants
        if not scaffold or scaffold not in variant_scaffolds:
            continue

        # Only print for known genes to reduce noise
        if satellite_name and satellite_name != 'Unknown' and not pd.isna(satellite_name):
            print(f"Looking for satellite gene with ID: {satellite_id} or name: {satellite_name}")

        # Filter variants on the same scaffold
        scaffold_variants = variants_df[variants_df['Scaffold'] == scaffold]
        gene_variants = pd.DataFrame()

        # 1. Look for variants within the satellite gene
        if start > 0 and end > 0:
            within_variants = scaffold_variants[
                (scaffold_variants['Position'] >= start) &
                (scaffold_variants['Position'] <= end)
            ]
            if not within_variants.empty:
                gene_variants = pd.concat([gene_variants, within_variants])

        # 2. Look for variants associated with the satellite gene by name/ID
        if satellite_name and satellite_name != 'Unknown' and not pd.isna(satellite_name):
            try:
                # Direct ID match
                id_matches = scaffold_variants[
                    (scaffold_variants['SC_Gene_ID'] == satellite_name) |
                    (scaffold_variants['Gene_ID'] == satellite_name)
                ]
                if not id_matches.empty:
                    gene_variants = pd.concat([gene_variants, id_matches])

                # Partial ID match for systematic names
                if satellite_name.startswith('Y') and len(satellite_name) >= 7:
                    partial_matches = scaffold_variants[
                        scaffold_variants['SC_Gene_ID'].str.contains(satellite_name,
                                                                    regex=False,
                                                                    na=False)
                    ]
                    if not partial_matches.empty:
                        gene_variants = pd.concat([gene_variants, partial_matches])
            except Exception as e:
                # Skip errors in string operations
                pass

        # 3. Look for variants that refer to the ERG gene this satellite is associated with
        if erg_gene_name:
            erg_variants = scaffold_variants[
                (scaffold_variants['ERG_Name'] == erg_gene_name)
            ]
            if not erg_variants.empty:
                # For each ERG variant, check if it's also associated with this satellite gene
                for _, erg_variant in erg_variants.iterrows():
                    erg_variant_pos = erg_variant['Position']
                    # If the variant position is closer to the satellite gene than to the ERG gene
                    if erg_start > 0 and erg_end > 0 and start > 0 and end > 0:
                        # Calculate distances to satellite gene and ERG gene
                        dist_to_satellite = min(abs(erg_variant_pos - start), abs(erg_variant_pos - end))
                        dist_to_erg = min(abs(erg_variant_pos - erg_start), abs(erg_variant_pos - erg_end))

                        # Add variants that are closer to satellite gene than ERG gene
                        if dist_to_satellite < dist_to_erg:
                            gene_variants = pd.concat([gene_variants, pd.DataFrame([erg_variant])])

        # Remove duplicates
        gene_variants = gene_variants.drop_duplicates()
        
        # Add satellite gene information to variants
        for _, variant in gene_variants.iterrows():
            variant_data = variant.to_dict()

            # Map original column names to our standard names for consistency
            # This ensures later analysis functions will work correctly
            standard_variant = {
                'treatment': variant_data.get('Treatment', ''),
                'scaffold': variant_data.get('Scaffold', ''),
                'position': variant_data.get('Position', 0),
                'ref': variant_data.get('Ref', ''),
                'alt': variant_data.get('Alt', ''),
                'quality': variant_data.get('Quality', 0),
                'gene_id': variant_data.get('Gene_ID', ''),
                'erg_gene_name': variant_data.get('ERG_Name', ''),
                'sc_gene_id': variant_data.get('SC_Gene_ID', ''),
                'effect': variant_data.get('Effect', ''),
                'impact': variant_data.get('Impact', ''),
                'satellite_gene_id': satellite_id,
                'satellite_gene_name': satellite_name,
                'satellite_start': start,
                'satellite_end': end,
                'erg_gene_id': erg_gene_id,
                'erg_gene_name': erg_gene_name,
                'function_category': function_category
            }
            satellite_variants.append(standard_variant)
    
    # Convert to DataFrame
    if satellite_variants:
        satellite_variants_df = pd.DataFrame(satellite_variants)
        print(f"Found {len(satellite_variants_df)} variants in satellite genes")
        return satellite_variants_df
    else:
        print("No variants found in satellite genes")
        return pd.DataFrame()

def analyze_satellite_variants(satellite_variants_df):
    """Analyze patterns of variants in satellite genes"""
    if satellite_variants_df.empty:
        print("ERROR: No satellite variants to analyze")
        return {}
    
    # Initialize results dictionary
    results = {
        'total_variants': len(satellite_variants_df),
        'variants_by_treatment': {},
        'variants_by_function': {},
        'variants_by_impact': {},
        'variants_by_effect': {},
        'variants_by_erg_gene': {},
        'treatment_function_crosstab': None,
        'treatment_impact_crosstab': None,
        'treatment_effect_crosstab': None,
        'treatment_erg_crosstab': None,
        'function_impact_crosstab': None,
        'statistical_tests': {}
    }
    
    # Process variants by treatment
    for treatment in TREATMENTS.keys():
        treatment_variants = satellite_variants_df[satellite_variants_df['treatment'] == treatment]
        results['variants_by_treatment'][treatment] = len(treatment_variants)
    
    # Process variants by function category
    function_counts = satellite_variants_df['function_category'].value_counts()
    for function, count in function_counts.items():
        results['variants_by_function'][function] = count
    
    # Process variants by impact
    impact_counts = satellite_variants_df['impact'].value_counts()
    for impact, count in impact_counts.items():
        results['variants_by_impact'][impact] = count
    
    # Process variants by effect
    effect_counts = satellite_variants_df['effect'].value_counts()
    for effect, count in effect_counts.items():
        results['variants_by_effect'][effect] = count
    
    # Process variants by ERG gene
    erg_counts = satellite_variants_df['erg_gene_name'].value_counts()
    for erg_gene, count in erg_counts.items():
        results['variants_by_erg_gene'][erg_gene] = count
    
    # Create cross-tabulations
    # Treatment vs. Function Category
    results['treatment_function_crosstab'] = pd.crosstab(
        satellite_variants_df['treatment'], 
        satellite_variants_df['function_category']
    )
    
    # Treatment vs. Impact
    results['treatment_impact_crosstab'] = pd.crosstab(
        satellite_variants_df['treatment'], 
        satellite_variants_df['impact']
    )
    
    # Treatment vs. Effect
    results['treatment_effect_crosstab'] = pd.crosstab(
        satellite_variants_df['treatment'], 
        satellite_variants_df['effect']
    )
    
    # Treatment vs. ERG Gene
    results['treatment_erg_crosstab'] = pd.crosstab(
        satellite_variants_df['treatment'], 
        satellite_variants_df['erg_gene_name']
    )
    
    # Function vs. Impact
    results['function_impact_crosstab'] = pd.crosstab(
        satellite_variants_df['function_category'], 
        satellite_variants_df['impact']
    )
    
    # Statistical tests
    # Chi-square test for independence between treatment and function category
    try:
        chi2, p, dof, expected = stats.chi2_contingency(results['treatment_function_crosstab'])
        results['statistical_tests']['treatment_function'] = {
            'test': 'chi2_contingency',
            'chi2': chi2,
            'p_value': p,
            'dof': dof
        }
    except:
        results['statistical_tests']['treatment_function'] = {
            'test': 'chi2_contingency',
            'error': 'Failed to calculate'
        }
    
    # Chi-square test for independence between treatment and impact
    try:
        chi2, p, dof, expected = stats.chi2_contingency(results['treatment_impact_crosstab'])
        results['statistical_tests']['treatment_impact'] = {
            'test': 'chi2_contingency',
            'chi2': chi2,
            'p_value': p,
            'dof': dof
        }
    except:
        results['statistical_tests']['treatment_impact'] = {
            'test': 'chi2_contingency',
            'error': 'Failed to calculate'
        }
    
    # Compare variants between adaptation types (Temperature vs. Low Oxygen)
    temp_variants = satellite_variants_df[satellite_variants_df['treatment'].isin(['WT-37', 'CAS'])]
    oxygen_variants = satellite_variants_df[satellite_variants_df['treatment'].isin(['WTA', 'STC'])]
    
    # Test for difference in HIGH impact variants between adaptation types
    temp_high = len(temp_variants[temp_variants['impact'] == 'HIGH'])
    temp_other = len(temp_variants) - temp_high
    oxygen_high = len(oxygen_variants[oxygen_variants['impact'] == 'HIGH'])
    oxygen_other = len(oxygen_variants) - oxygen_high
    
    try:
        table = [[temp_high, temp_other], [oxygen_high, oxygen_other]]
        odds_ratio, p_value = stats.fisher_exact(table)
        results['statistical_tests']['adaptation_high_impact'] = {
            'test': 'fisher_exact',
            'odds_ratio': odds_ratio,
            'p_value': p_value,
            'temp_high': temp_high,
            'temp_total': len(temp_variants),
            'oxygen_high': oxygen_high,
            'oxygen_total': len(oxygen_variants)
        }
    except:
        results['statistical_tests']['adaptation_high_impact'] = {
            'test': 'fisher_exact',
            'error': 'Failed to calculate'
        }
    
    return results

def enrich_variant_analysis(satellite_variants_df):
    """Add additional analysis of satellite variants"""
    if satellite_variants_df.empty:
        return pd.DataFrame(), {}
    
    # Calculate distance from gene start/end
    variant_distances = []
    for i, variant in satellite_variants_df.iterrows():
        position = variant.get('position', 0)
        gene_start = variant.get('satellite_start', 0)
        gene_end = variant.get('satellite_end', 0)
        
        if gene_start <= position <= gene_end:
            # Variant is within the gene
            rel_position = "within"
            distance = 0
        elif position < gene_start:
            # Variant is upstream of the gene
            rel_position = "upstream"
            distance = gene_start - position
        else:
            # Variant is downstream of the gene
            rel_position = "downstream"
            distance = position - gene_end
        
        variant_distances.append({
            'index': i,
            'relative_position': rel_position,
            'distance_to_gene': distance,
            'distance_kb': round(distance / 1000, 2)
        })
    
    # Convert to DataFrame and merge with variants
    if variant_distances:
        distances_df = pd.DataFrame(variant_distances)
        enriched_df = pd.merge(
            satellite_variants_df, 
            distances_df, 
            left_index=True, 
            right_on='index',
            how='left'
        )
        
        # Calculate additional statistics
        stats = {
            'within_gene': len(enriched_df[enriched_df['relative_position'] == 'within']),
            'upstream': len(enriched_df[enriched_df['relative_position'] == 'upstream']),
            'downstream': len(enriched_df[enriched_df['relative_position'] == 'downstream']),
            'avg_distance': enriched_df['distance_to_gene'].mean(),
            'median_distance': enriched_df['distance_to_gene'].median(),
            'distances_by_impact': {},
            'distances_by_function': {}
        }
        
        # Calculate distances by impact
        for impact in IMPACT_CATEGORIES:
            impact_variants = enriched_df[enriched_df['impact'] == impact]
            if not impact_variants.empty:
                stats['distances_by_impact'][impact] = {
                    'count': len(impact_variants),
                    'avg_distance': impact_variants['distance_to_gene'].mean(),
                    'median_distance': impact_variants['distance_to_gene'].median()
                }
        
        # Calculate distances by function
        for function in enriched_df['function_category'].unique():
            function_variants = enriched_df[enriched_df['function_category'] == function]
            if not function_variants.empty:
                stats['distances_by_function'][function] = {
                    'count': len(function_variants),
                    'avg_distance': function_variants['distance_to_gene'].mean(),
                    'median_distance': function_variants['distance_to_gene'].median()
                }
        
        return enriched_df, stats
    else:
        return satellite_variants_df, {}

def visualize_satellite_variants(satellite_variants_df, analysis_results, distance_stats, output_dir):
    """Create visualizations of satellite variant patterns"""
    if satellite_variants_df.empty:
        print("ERROR: No satellite variants to visualize")
        return
    
    setup_plotting_style()
    
    # Create output directory for visualizations
    viz_dir = os.path.join(output_dir, "visualizations")
    ensure_dir(viz_dir)
    
    # 1. Variants by treatment
    plt.figure(figsize=(10, 6))
    treatment_counts = [analysis_results['variants_by_treatment'].get(t, 0) for t in TREATMENTS.keys()]
    sns.barplot(x=list(TREATMENTS.keys()), y=treatment_counts)
    plt.title('Satellite Gene Variants by Treatment')
    plt.xlabel('Treatment')
    plt.ylabel('Number of Variants')
    plt.xticks(rotation=45)
    plt.tight_layout()
    save_plot(plt, os.path.join(viz_dir, "satellite_variants_by_treatment.png"))
    plt.close()
    
    # 2. Variants by function category
    plt.figure(figsize=(12, 6))
    function_data = pd.DataFrame({
        'Function': list(analysis_results['variants_by_function'].keys()),
        'Variants': list(analysis_results['variants_by_function'].values())
    }).sort_values('Variants', ascending=False)
    
    sns.barplot(x='Function', y='Variants', data=function_data)
    plt.title('Satellite Gene Variants by Functional Category')
    plt.xlabel('Functional Category')
    plt.ylabel('Number of Variants')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_plot(plt, os.path.join(viz_dir, "satellite_variants_by_function.png"))
    plt.close()
    
    # 3. Variants by impact
    plt.figure(figsize=(10, 6))
    impact_order = [imp for imp in IMPACT_CATEGORIES if imp in analysis_results['variants_by_impact']]
    impact_counts = [analysis_results['variants_by_impact'].get(imp, 0) for imp in impact_order]
    
    sns.barplot(x=impact_order, y=impact_counts)
    plt.title('Satellite Gene Variants by Impact')
    plt.xlabel('Impact')
    plt.ylabel('Number of Variants')
    plt.tight_layout()
    save_plot(plt, os.path.join(viz_dir, "satellite_variants_by_impact.png"))
    plt.close()
    
    # 4. Variants by effect (top 10)
    plt.figure(figsize=(12, 6))
    effect_data = pd.DataFrame({
        'Effect': list(analysis_results['variants_by_effect'].keys()),
        'Variants': list(analysis_results['variants_by_effect'].values())
    }).sort_values('Variants', ascending=False).head(10)
    
    sns.barplot(x='Effect', y='Variants', data=effect_data)
    plt.title('Top 10 Effects of Satellite Gene Variants')
    plt.xlabel('Effect')
    plt.ylabel('Number of Variants')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_plot(plt, os.path.join(viz_dir, "satellite_variants_by_effect.png"))
    plt.close()
    
    # 5. Treatment vs. Function Category Heatmap
    if analysis_results['treatment_function_crosstab'] is not None:
        plt.figure(figsize=(14, 8))
        sns.heatmap(analysis_results['treatment_function_crosstab'], annot=True, fmt='d', cmap='YlGnBu')
        plt.title('Satellite Gene Variants by Treatment and Functional Category')
        plt.tight_layout()
        save_plot(plt, os.path.join(viz_dir, "satellite_treatment_function_heatmap.png"))
        plt.close()
    
    # 6. Treatment vs. Impact Heatmap
    if analysis_results['treatment_impact_crosstab'] is not None:
        plt.figure(figsize=(12, 8))
        sns.heatmap(analysis_results['treatment_impact_crosstab'], annot=True, fmt='d', cmap='YlGnBu')
        plt.title('Satellite Gene Variants by Treatment and Impact')
        plt.tight_layout()
        save_plot(plt, os.path.join(viz_dir, "satellite_treatment_impact_heatmap.png"))
        plt.close()
    
    # 7. Adaptation-specific analysis
    # Prepare data by grouping treatments into adaptation types
    adaptation_data = []
    for adaptation, treatments in ADAPTATION_GROUPS.items():
        for treatment in treatments:
            treatment_variants = satellite_variants_df[satellite_variants_df['treatment'] == treatment]
            for _, variant in treatment_variants.iterrows():
                adaptation_data.append({
                    'adaptation': adaptation,
                    'impact': variant.get('impact', 'Unknown'),
                    'function_category': variant.get('function_category', 'Unknown'),
                    'effect': variant.get('effect', 'Unknown')
                })
    
    if adaptation_data:
        adaptation_df = pd.DataFrame(adaptation_data)
        
        # Impact distribution by adaptation type
        plt.figure(figsize=(12, 6))
        impact_crosstab = pd.crosstab(adaptation_df['adaptation'], adaptation_df['impact'])
        impact_crosstab_pct = impact_crosstab.div(impact_crosstab.sum(axis=1), axis=0)
        
        impact_crosstab_pct.plot(kind='bar', stacked=True, colormap='viridis')
        plt.title('Distribution of Variant Impacts by Adaptation Type')
        plt.xlabel('Adaptation Type')
        plt.ylabel('Proportion of Variants')
        plt.legend(title='Impact')
        plt.tight_layout()
        save_plot(plt, os.path.join(viz_dir, "satellite_impact_by_adaptation.png"))
        plt.close()
        
        # Function category distribution by adaptation type
        plt.figure(figsize=(14, 6))
        function_crosstab = pd.crosstab(adaptation_df['adaptation'], adaptation_df['function_category'])
        function_crosstab_pct = function_crosstab.div(function_crosstab.sum(axis=1), axis=0)
        
        function_crosstab_pct.plot(kind='bar', stacked=True, colormap='viridis')
        plt.title('Distribution of Functional Categories by Adaptation Type')
        plt.xlabel('Adaptation Type')
        plt.ylabel('Proportion of Variants')
        plt.legend(title='Function', loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        save_plot(plt, os.path.join(viz_dir, "satellite_function_by_adaptation.png"))
        plt.close()
    
    # 8. Distance analysis
    if 'distance_to_gene' in satellite_variants_df.columns:
        # Distance distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(satellite_variants_df['distance_to_gene'], bins=30, kde=True)
        plt.title('Distribution of Variant Distances to Satellite Genes')
        plt.xlabel('Distance (bp)')
        plt.ylabel('Number of Variants')
        plt.tight_layout()
        save_plot(plt, os.path.join(viz_dir, "satellite_variant_distances.png"))
        plt.close()
        
        # Distance vs. Impact
        plt.figure(figsize=(12, 6))
        sns.boxplot(x='impact', y='distance_to_gene', data=satellite_variants_df)
        plt.title('Variant Distances by Impact')
        plt.xlabel('Impact')
        plt.ylabel('Distance to Gene (bp)')
        plt.tight_layout()
        save_plot(plt, os.path.join(viz_dir, "satellite_distance_by_impact.png"))
        plt.close()
        
        # Distance vs. Function Category
        plt.figure(figsize=(14, 6))
        sns.boxplot(x='function_category', y='distance_to_gene', data=satellite_variants_df)
        plt.title('Variant Distances by Functional Category')
        plt.xlabel('Functional Category')
        plt.ylabel('Distance to Gene (bp)')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        save_plot(plt, os.path.join(viz_dir, "satellite_distance_by_function_box.png"))
        plt.close()
    
    # 9. ERG gene analysis - which satellite genes around which ERG genes have the most variants
    plt.figure(figsize=(12, 6))
    erg_data = pd.DataFrame({
        'ERG Gene': list(analysis_results['variants_by_erg_gene'].keys()),
        'Variants': list(analysis_results['variants_by_erg_gene'].values())
    }).sort_values('Variants', ascending=False)
    
    sns.barplot(x='ERG Gene', y='Variants', data=erg_data)
    plt.title('Satellite Gene Variants by Associated ERG Gene')
    plt.xlabel('ERG Gene')
    plt.ylabel('Number of Variants')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_plot(plt, os.path.join(viz_dir, "satellite_variants_by_erg_gene.png"))
    plt.close()
    
    # 10. Treatment vs. ERG Gene Heatmap
    if analysis_results['treatment_erg_crosstab'] is not None:
        plt.figure(figsize=(14, 8))
        sns.heatmap(analysis_results['treatment_erg_crosstab'], annot=True, fmt='d', cmap='YlGnBu')
        plt.title('Satellite Gene Variants by Treatment and ERG Gene')
        plt.tight_layout()
        save_plot(plt, os.path.join(viz_dir, "satellite_treatment_erg_heatmap.png"))
        plt.close()

def generate_variant_profile_report(satellite_variants_df, analysis_results, distance_stats, output_dir):
    """Generate a detailed report of satellite gene variant profiles"""
    report_file = os.path.join(output_dir, "satellite_variant_profiling_report.txt")
    
    with open(report_file, 'w') as f:
        f.write("Satellite Gene Variant Profiling Report\n")
        f.write("=====================================\n\n")
        
        f.write("1. Overview\n")
        f.write("-----------\n")
        f.write(f"Total variants analyzed: {analysis_results.get('total_variants', 0)}\n\n")
        
        f.write("Variants by Treatment:\n")
        for treatment, count in analysis_results.get('variants_by_treatment', {}).items():
            f.write(f"- {treatment}: {count} variants\n")
        
        f.write("\nVariants by Functional Category:\n")
        for function, count in sorted(analysis_results.get('variants_by_function', {}).items(), key=lambda x: x[1], reverse=True):
            f.write(f"- {function}: {count} variants\n")
        
        f.write("\nVariants by Impact:\n")
        for impact in IMPACT_CATEGORIES:
            count = analysis_results.get('variants_by_impact', {}).get(impact, 0)
            f.write(f"- {impact}: {count} variants\n")
        
        f.write("\n2. Statistical Analysis\n")
        f.write("----------------------\n")
        
        # Chi-square test for treatment and function category
        treatment_function_test = analysis_results.get('statistical_tests', {}).get('treatment_function', {})
        if 'error' not in treatment_function_test:
            chi2 = treatment_function_test.get('chi2', 0)
            p_value = treatment_function_test.get('p_value', 1)
            dof = treatment_function_test.get('dof', 0)
            
            f.write("Association between Treatment and Functional Category:\n")
            f.write(f"- Chi-square statistic: {chi2:.2f}\n")
            f.write(f"- Degrees of freedom: {dof}\n")
            f.write(f"- P-value: {p_value:.4f}\n")
            f.write(f"- Interpretation: {'Significant association' if p_value < 0.05 else 'No significant association'}\n\n")
        
        # Chi-square test for treatment and impact
        treatment_impact_test = analysis_results.get('statistical_tests', {}).get('treatment_impact', {})
        if 'error' not in treatment_impact_test:
            chi2 = treatment_impact_test.get('chi2', 0)
            p_value = treatment_impact_test.get('p_value', 1)
            dof = treatment_impact_test.get('dof', 0)
            
            f.write("Association between Treatment and Impact:\n")
            f.write(f"- Chi-square statistic: {chi2:.2f}\n")
            f.write(f"- Degrees of freedom: {dof}\n")
            f.write(f"- P-value: {p_value:.4f}\n")
            f.write(f"- Interpretation: {'Significant association' if p_value < 0.05 else 'No significant association'}\n\n")
        
        # Fisher's exact test for high impact variants between adaptation types
        adaptation_high_test = analysis_results.get('statistical_tests', {}).get('adaptation_high_impact', {})
        if 'error' not in adaptation_high_test:
            odds_ratio = adaptation_high_test.get('odds_ratio', 1)
            p_value = adaptation_high_test.get('p_value', 1)
            temp_high = adaptation_high_test.get('temp_high', 0)
            temp_total = adaptation_high_test.get('temp_total', 0)
            oxygen_high = adaptation_high_test.get('oxygen_high', 0)
            oxygen_total = adaptation_high_test.get('oxygen_total', 0)
            
            f.write("Comparison of HIGH Impact Variants between Adaptation Types:\n")
            f.write(f"- Temperature adaptation: {temp_high}/{temp_total} high impact variants ({temp_high/temp_total*100:.1f}%)\n")
            f.write(f"- Low oxygen adaptation: {oxygen_high}/{oxygen_total} high impact variants ({oxygen_high/oxygen_total*100:.1f}%)\n")
            f.write(f"- Odds ratio: {odds_ratio:.2f}\n")
            f.write(f"- P-value: {p_value:.4f}\n")
            f.write(f"- Interpretation: {'Significant difference' if p_value < 0.05 else 'No significant difference'}\n\n")
        
        f.write("3. Variant Locations\n")
        f.write("------------------\n")
        
        if distance_stats:
            f.write(f"Variants within satellite genes: {distance_stats.get('within_gene', 0)}\n")
            f.write(f"Variants upstream of satellite genes: {distance_stats.get('upstream', 0)}\n")
            f.write(f"Variants downstream of satellite genes: {distance_stats.get('downstream', 0)}\n")
            f.write(f"Average distance to gene: {distance_stats.get('avg_distance', 0):.1f} bp\n")
            f.write(f"Median distance to gene: {distance_stats.get('median_distance', 0):.1f} bp\n\n")
            
            f.write("Distance by Impact:\n")
            for impact, stats in distance_stats.get('distances_by_impact', {}).items():
                f.write(f"- {impact}: {stats.get('count', 0)} variants, ")
                f.write(f"avg distance {stats.get('avg_distance', 0):.1f} bp, ")
                f.write(f"median distance {stats.get('median_distance', 0):.1f} bp\n")
            
            f.write("\nDistance by Functional Category:\n")
            for function, stats in distance_stats.get('distances_by_function', {}).items():
                f.write(f"- {function}: {stats.get('count', 0)} variants, ")
                f.write(f"avg distance {stats.get('avg_distance', 0):.1f} bp, ")
                f.write(f"median distance {stats.get('median_distance', 0):.1f} bp\n")
        
        f.write("\n4. ERG Gene Association\n")
        f.write("---------------------\n")
        
        f.write("Variants by ERG Gene:\n")
        for erg_gene, count in sorted(analysis_results.get('variants_by_erg_gene', {}).items(), key=lambda x: x[1], reverse=True):
            f.write(f"- {erg_gene}: {count} variants\n")
        
        f.write("\n5. Adaptation-Specific Patterns\n")
        f.write("-----------------------------\n")
        
        # Adaptation-specific patterns
        temp_variants = satellite_variants_df[satellite_variants_df['treatment'].isin(['WT-37', 'CAS'])]
        oxygen_variants = satellite_variants_df[satellite_variants_df['treatment'].isin(['WTA', 'STC'])]
        
        f.write(f"Temperature adaptation (WT-37, CAS): {len(temp_variants)} variants\n")
        f.write(f"Low oxygen adaptation (WTA, STC): {len(oxygen_variants)} variants\n\n")
        
        # Impact distribution by adaptation
        f.write("Impact Distribution by Adaptation:\n")
        for impact in IMPACT_CATEGORIES:
            temp_impact = len(temp_variants[temp_variants['impact'] == impact])
            oxygen_impact = len(oxygen_variants[oxygen_variants['impact'] == impact])
            
            temp_pct = temp_impact / len(temp_variants) * 100 if len(temp_variants) > 0 else 0
            oxygen_pct = oxygen_impact / len(oxygen_variants) * 100 if len(oxygen_variants) > 0 else 0
            
            f.write(f"- {impact}: Temperature {temp_impact} ({temp_pct:.1f}%), ")
            f.write(f"Low Oxygen {oxygen_impact} ({oxygen_pct:.1f}%)\n")
        
        f.write("\nFor detailed analysis, see the satellite_variants.tsv file and visualizations directory.\n")
    
    print(f"Variant profile report saved to {report_file}")

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Ensure output directory exists
    ensure_dir(args.output_dir)
    
    print("======================================================")
    print("Satellite Gene Variant Profiling")
    print("======================================================")
    
    # Load annotated satellite genes
    satellite_df = load_satellite_genes(args.satellite_genes)
    
    # Load variant data
    variants_df = load_variants(args.variant_dir)

    # Check for scaffold overlap between satellite genes and variants
    if not satellite_df.empty and not variants_df.empty:
        satellite_scaffolds = satellite_df['scaffold'].unique()
        variant_scaffolds = variants_df['Scaffold'].unique()

        # Find common scaffolds
        common_scaffolds = set(satellite_scaffolds) & set(variant_scaffolds)

        print(f"\nScaffold analysis:")
        print(f"Satellite genes found on {len(satellite_scaffolds)} unique scaffolds")
        print(f"Variants found on {len(variant_scaffolds)} unique scaffolds")
        print(f"Scaffolds in common: {len(common_scaffolds)}")

        if len(common_scaffolds) > 0:
            print(f"Common scaffolds: {common_scaffolds}")

            # Count satellite genes on common scaffolds
            satellite_on_common = satellite_df[satellite_df['scaffold'].isin(common_scaffolds)]
            print(f"Satellite genes on common scaffolds: {len(satellite_on_common)} ({len(satellite_on_common)/len(satellite_df)*100:.1f}% of total)")

            # Analyze variant positions vs satellite gene positions
            print("\nAnalyzing variant positions relative to satellite genes:")

            # Compute min/max positions of satellite genes in common scaffolds
            for scaffold in common_scaffolds:
                scaffold_satellites = satellite_df[satellite_df['scaffold'] == scaffold]
                scaffold_variants = variants_df[variants_df['Scaffold'] == scaffold]

                satellite_min_pos = scaffold_satellites['satellite_start'].min()
                satellite_max_pos = scaffold_satellites['satellite_end'].max()

                variant_min_pos = scaffold_variants['Position'].min()
                variant_max_pos = scaffold_variants['Position'].max()

                print(f"  Scaffold {scaffold}:")
                print(f"    Satellite genes positions: {satellite_min_pos} - {satellite_max_pos}")
                print(f"    Variants positions: {variant_min_pos} - {variant_max_pos}")

                # Check if there's position overlap
                if (variant_min_pos <= satellite_max_pos and variant_max_pos >= satellite_min_pos):
                    print(f"    OVERLAP: Some variants may fall within satellite gene positions")
                else:
                    print(f"    NO OVERLAP: Variants and satellite genes are in different regions")

                # Look at the ERG genes on this scaffold
                scaffold_ergs = []
                for erg_name in scaffold_satellites['erg_gene_name'].unique():
                    if erg_name in scaffold_ergs:
                        continue
                    scaffold_ergs.append(erg_name)

                    # Get position of this ERG gene
                    erg_rows = scaffold_satellites[scaffold_satellites['erg_gene_name'] == erg_name]
                    if not erg_rows.empty:
                        erg_start = erg_rows['erg_start'].iloc[0]
                        erg_end = erg_rows['erg_end'].iloc[0]

                        if erg_start > 0 and erg_end > 0:
                            # Check variants near this ERG gene
                            near_variants = scaffold_variants[
                                (scaffold_variants['Position'] >= erg_start - 5000) &
                                (scaffold_variants['Position'] <= erg_end + 5000)
                            ]
                            if len(near_variants) > 0:
                                print(f"    ERG gene {erg_name} position: {erg_start} - {erg_end}")
                                print(f"    Variants near this ERG gene: {len(near_variants)}")
                                print(f"    Confirming these variants are in the buffer zone, not satellite zone")
        else:
            print("WARNING: No common scaffolds between satellite genes and variants!")
            print("This explains why no satellite variants will be found.")

    # Identify variants in satellite genes
    print("\nIdentifying variants in satellite genes...")
    satellite_variants_df = identify_satellite_variants(satellite_df, variants_df)
    
    if not satellite_variants_df.empty:
        # Analyze satellite variants
        print("Analyzing variant patterns in satellite genes...")
        analysis_results = analyze_satellite_variants(satellite_variants_df)
        
        # Enrich the variant analysis with additional metrics
        print("Enriching variant analysis with distance metrics...")
        satellite_variants_df, distance_stats = enrich_variant_analysis(satellite_variants_df)
        
        # Save the satellite variants data
        satellite_variants_file = os.path.join(args.output_dir, "satellite_variants.tsv")
        save_tsv(satellite_variants_df, satellite_variants_file)
        print(f"Satellite variants data saved to {satellite_variants_file}")
        
        # Create visualizations
        print("Creating visualizations of variant patterns...")
        visualize_satellite_variants(satellite_variants_df, analysis_results, distance_stats, args.output_dir)
        
        # Generate detailed report
        print("Generating variant profile report...")
        generate_variant_profile_report(satellite_variants_df, analysis_results, distance_stats, args.output_dir)
    
    print("Satellite gene variant profiling complete!")

if __name__ == "__main__":
    main()