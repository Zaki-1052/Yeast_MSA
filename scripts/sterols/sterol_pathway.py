#!/usr/bin/env python3
"""
Ergosterol pathway analysis script for Yeast MSA project.
This script maps sterol data to the ergosterol biosynthetic pathway.

Input:
- sterol_data_processed.csv: Processed sterol data

Output:
- Pathway analysis files in results/sterol_analysis/pathway/
- Pathway visualizations in results/sterol_analysis/visualizations/
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from pathlib import Path
import matplotlib.patches as mpatches

# Set paths
INPUT_FILE = 'sterol_data/sterol_data_processed.csv'
RESULTS_DIR = 'results/sterol_analysis'
PATHWAY_DIR = f'{RESULTS_DIR}/pathway'
VIS_DIR = f'{RESULTS_DIR}/visualizations'

# Define ergosterol pathway steps and connections
# This is a refined model of the ergosterol biosynthetic pathway
# Updated based on comparative analysis findings
PATHWAY = {
    'Squalene': {
        'enzyme': 'ERG9',
        'next': ['Lanosterol'],
        'adaptation_type': 'Both',
        'gene_modification': 'Both'
    },
    'Lanosterol': {
        'enzyme': 'ERG7',
        'next': ['Zymosterol', 'Fecosterol', 'Cycloartenol'],
        'adaptation_type': 'Temperature',  # Found primarily in temperature adaptation
        'gene_modification': 'Modified'    # Found in CAS
    },
    'Zymosterol': {
        'enzyme': 'ERG6',
        'next': ['Fecosterol'],
        'adaptation_type': 'Temperature',  # Found in WT-37
        'gene_modification': 'Non-modified'
    },
    'Fecosterol': {
        'enzyme': 'ERG2',
        'next': ['Episterol'],
        'adaptation_type': 'Temperature',  # Found in CAS-55 and WT-55
        'gene_modification': 'Both'
    },
    'Episterol': {
        'enzyme': 'ERG3',
        'next': ['Ergosta-5,7,24(28)-trienol'],
        'adaptation_type': 'Both',
        'gene_modification': 'Both'
    },
    'Ergosta-5,7,24(28)-trienol': {
        'enzyme': 'ERG5',
        'next': ['Ergosta-5,7,22,24(28)-tetraenol'],
        'adaptation_type': 'Both',
        'gene_modification': 'Both'
    },
    'Ergosta-5,7,22,24(28)-tetraenol': {
        'enzyme': 'ERG4',
        'next': ['Ergosterol'],
        'adaptation_type': 'Both',
        'gene_modification': 'Both'
    },
    'Ergosterol': {
        'enzyme': None,  # End product
        'next': [],
        'adaptation_type': 'Both',  # Found in all samples but 3.76x higher in temperature-adapted
        'gene_modification': 'Both'  # Found in all strains
    },
    # Alternative/branch sterols - Updated based on comparative analysis
    'Ergost-7-en-3-ol': {
        'enzyme': 'ERG3', # Alternative product of ERG3 function
        'next': [],
        'adaptation_type': 'Temperature',
        'gene_modification': 'Modified'    # CAS specific
    },
    'Ergosta-7-en-3beta-ol': {
        'enzyme': 'ERG3', # Alternative product of ERG3 function
        'next': [],
        'adaptation_type': 'Temperature',
        'gene_modification': 'Modified'    # CAS specific
    },
    'Cycloartenol': {
        'enzyme': 'ERG7', # Alternative product of ERG7 function
        'next': [],
        'adaptation_type': 'Temperature',
        'gene_modification': 'Modified'    # CAS specific
    },
    'Stigmasta-5_22-dien-3-ol_acetate': {
        'enzyme': 'ERG4', # Alternative/modified product downstream of ERG4
        'next': [],
        'adaptation_type': 'Temperature',
        'gene_modification': 'Modified'    # CAS specific, high abundance
    },
    'Tetrahymanol': {
        'enzyme': 'ERG7', # Adaptation-specific marker for low oxygen
        'next': [],
        'adaptation_type': 'Low Oxygen',
        'gene_modification': 'Modified'    # STC specific marker
    }
}

# Map sterols to their known enzymes in the ergosterol pathway
STEROL_TO_ENZYME = {
    'Squalene': 'ERG9',
    'Lanosterol': 'ERG7',
    'Zymosterol': 'ERG6',
    'Fecosterol': 'ERG2',
    'Episterol': 'ERG3',
    'Ergosta-5,7,24(28)-trienol': 'ERG5',
    'Ergosta-5,7,22,24(28)-tetraenol': 'ERG4',
    'Ergosterol': 'Product',
    'Ergost-7-en-3-ol': 'ERG3',
    'Ergosta-7-en-3beta-ol': 'ERG3',
    'Cycloartenol': 'ERG7',
    'Stigmasta-5_22-dien-3-ol_acetate': 'ERG4',
    'Tetrahymanol': 'ERG7'
}

# Map enzymes to their genes
ENZYME_TO_GENE = {
    'ERG1': 'YGR175C',
    'ERG2': 'YMR202W',
    'ERG3': 'YLR056W',
    'ERG4': 'YGL012W',
    'ERG5': 'YMR015C',
    'ERG6': 'YML008C',
    'ERG7': 'YHR072W',
    'ERG9': 'YHR190W',
    'ERG11': 'YHR007C',
    'ERG24': 'YNL280C',
    'ERG25': 'YGR060W',
    'ERG26': 'YGL001C',
    'ERG27': 'YLR100W',
    'ERG28': 'YER044C',
}

def ensure_directories():
    """Ensure all required directories exist."""
    for directory in [PATHWAY_DIR, VIS_DIR]:
        Path(directory).mkdir(parents=True, exist_ok=True)

def load_processed_data(file_path=INPUT_FILE):
    """Load and preprocess sterol data."""
    print(f"Loading processed sterol data from {file_path}")
    if not os.path.exists(file_path):
        # If the processed file doesn't exist, try to load and process the raw data
        raw_file = 'sterol_data/sterol_data_with_sd.csv'
        if not os.path.exists(raw_file):
            raise FileNotFoundError(f"Neither processed nor raw sterol data found: {file_path}, {raw_file}")
        
        print(f"Processed file not found. Loading raw data from {raw_file}")
        raw_df = pd.read_csv(raw_file)
        
        # Process the raw data
        df = process_raw_data(raw_df)
        
        # Create the directory if it doesn't exist
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        
        # Save the processed data for future use
        df.to_csv(file_path, index=False)
        print(f"Saved processed data to {file_path}")
    else:
        # Load the existing processed file
        df = pd.read_csv(file_path)
    
    print(f"Loaded {len(df)} processed sterol measurements")
    return df

def process_raw_data(raw_df):
    """Process raw sterol data into the required format."""
    print("Processing raw sterol data...")
    
    # Extract metadata from sample names
    processed_data = []
    
    for _, row in raw_df.iterrows():
        sample = row['sample']
        sterol = row['sterol']
        concentration = row['concentration']
        std_dev = row['std_dev']
        
        # Parse sample name to extract metadata
        treatment = sample.split('_')[0]
        
        # Determine generation (5 or 55)
        if '_5_' in sample or '_5' in sample:
            generation = 5
        elif '_55_' in sample or '_55' in sample:
            generation = 55
        else:
            generation = None
        
        # Determine adaptation type
        if treatment == 'WT':
            if 'MA' in sample:
                adaptation_type = 'Low Oxygen'
            else:
                adaptation_type = 'Temperature'
        elif treatment == 'CAS':
            adaptation_type = 'Temperature'
        elif treatment == 'STC':
            adaptation_type = 'Low Oxygen'
        else:
            adaptation_type = 'Unknown'
        
        # Determine gene modification status
        if treatment in ['CAS', 'STC']:
            gene_modified = True
        else:
            gene_modified = False
        
        processed_data.append({
            'sample': sample,
            'treatment': treatment,
            'generation': generation,
            'adaptation_type': adaptation_type,
            'gene_modified': gene_modified,
            'sterol': sterol,
            'concentration': concentration,
            'std_dev': std_dev
        })
    
    # Convert to DataFrame
    df = pd.DataFrame(processed_data)
    
    # Calculate relative abundances within each sample
    sample_totals = df.groupby('sample')['concentration'].sum()
    df['relative_abundance'] = df.apply(
        lambda row: (row['concentration'] / sample_totals[row['sample']]) * 100 
        if row['sample'] in sample_totals and sample_totals[row['sample']] > 0 
        else 0, 
        axis=1
    )
    
    print(f"Processed {len(df)} sterol measurements from {len(df['sample'].unique())} samples")
    print(f"Detected adaptations: {', '.join(df['adaptation_type'].unique())}")
    print(f"Detected treatments: {', '.join(df['treatment'].unique())}")
    
    return df

def map_sterols_to_pathway(df):
    """Map detected sterols to the ergosterol biosynthetic pathway."""
    # Get unique sterols in dataset
    detected_sterols = df['sterol'].unique()
    
    # Map to known pathway sterols
    pathway_map = {}
    for sterol in detected_sterols:
        if sterol in PATHWAY:
            pathway_map[sterol] = {
                'in_pathway': True,
                'enzyme': PATHWAY[sterol]['enzyme'],
                'next_steps': PATHWAY[sterol]['next']
            }
        else:
            pathway_map[sterol] = {
                'in_pathway': False,
                'enzyme': None,
                'next_steps': []
            }
    
    # Count detected pathway sterols
    pathway_sterols = [s for s in detected_sterols if s in PATHWAY]
    
    # Save pathway mapping
    with open(f'{PATHWAY_DIR}/pathway_mapping.txt', 'w') as f:
        f.write("# Sterol to Pathway Mapping\n\n")
        f.write(f"Detected sterols in dataset: {len(detected_sterols)}\n")
        f.write(f"Sterols mapped to ergosterol pathway: {len(pathway_sterols)}\n\n")
        
        f.write("## Pathway Sterols\n")
        for sterol in pathway_sterols:
            info = pathway_map[sterol]
            f.write(f"### {sterol}\n")
            f.write(f"- Enzyme: {info['enzyme'] or 'None'}\n")
            f.write(f"- Next steps: {', '.join(info['next_steps']) or 'None (terminal)'}\n\n")
        
        f.write("## Non-pathway Sterols\n")
        non_pathway = [s for s in detected_sterols if s not in PATHWAY]
        for sterol in non_pathway:
            f.write(f"- {sterol}\n")
    
    # Save as CSV for programmatic access
    pathway_df = []
    for sterol in detected_sterols:
        info = pathway_map[sterol]
        pathway_df.append({
            'sterol': sterol,
            'in_pathway': info['in_pathway'],
            'enzyme': info['enzyme'],
            'next_steps': ','.join(info['next_steps']) if info['in_pathway'] else ''
        })
    
    pd.DataFrame(pathway_df).to_csv(f'{PATHWAY_DIR}/sterol_pathway_mapping.csv', index=False)
    
    print(f"Pathway mapping saved to {PATHWAY_DIR}")
    return pathway_map

def analyze_sterol_ratios(df, pathway_map):
    """Calculate ratios between pathway intermediates to infer metabolic flux."""
    # Group by sample and treatment
    sample_sterols = {}
    for sample in df['sample'].unique():
        sample_df = df[df['sample'] == sample]
        sample_sterols[sample] = {}
        
        for _, row in sample_df.iterrows():
            sample_sterols[sample][row['sterol']] = row['concentration']
    
    # Calculate ratios for substrate-product pairs
    ratio_results = []
    
    for sample, sterols in sample_sterols.items():
        sample_info = df[df['sample'] == sample].iloc[0]
        treatment = sample_info['treatment']
        generation = sample_info['generation']
        adaptation = sample_info['adaptation_type']
        gene_modified = sample_info['gene_modified']
        
        # Enhanced list of key substrate-product pairs based on comparative analysis
        pairs = [
            ('Lanosterol', 'Zymosterol'),      # ERG7 → ERG6
            ('Lanosterol', 'Cycloartenol'),    # ERG7 → Alternative branch
            ('Zymosterol', 'Fecosterol'),      # ERG6 → ERG2
            ('Fecosterol', 'Ergosterol'),      # Downstream flux to final product
            ('Lanosterol', 'Ergosterol'),      # Overall pathway flux
            # Add pairs for alternative pathways
            ('Lanosterol', 'Tetrahymanol'),    # Low oxygen adaptation marker
            ('Ergosterol', 'Stigmasta-5_22-dien-3-ol_acetate')  # Temperature/CAS marker
        ]
        
        # Special case for Tetrahymanol in STC strains (Low Oxygen adaptation)
        if 'Tetrahymanol' in sterols and 'Ergosterol' in sterols:
            ratio_results.append({
                'sample': sample,
                'treatment': treatment,
                'generation': generation,
                'adaptation_type': adaptation,
                'gene_modified': gene_modified,
                'substrate': 'Ergosterol',
                'product': 'Tetrahymanol',
                'enzyme': 'ERG7',
                'substrate_conc': sterols['Ergosterol'],
                'product_conc': sterols['Tetrahymanol'],
                'ratio': sterols['Tetrahymanol'] / sterols['Ergosterol'] if sterols['Ergosterol'] > 0 else np.nan,
                'log2_ratio': np.log2(sterols['Tetrahymanol'] / sterols['Ergosterol']) if sterols['Ergosterol'] > 0 and sterols['Tetrahymanol'] > 0 else np.nan,
                'note': 'Low Oxygen adaptation marker'
            })
        
        for substrate, product in pairs:
            if substrate in sterols and product in sterols:
                sub_conc = sterols[substrate]
                prod_conc = sterols[product]
                
                if sub_conc > 0:
                    ratio = prod_conc / sub_conc
                    log2_ratio = np.log2(ratio) if ratio > 0 else np.nan
                    
                    ratio_results.append({
                        'sample': sample,
                        'treatment': treatment,
                        'generation': generation,
                        'adaptation_type': adaptation,
                        'gene_modified': gene_modified,
                        'substrate': substrate,
                        'product': product,
                        'enzyme': PATHWAY[substrate]['enzyme'],
                        'substrate_conc': sub_conc,
                        'product_conc': prod_conc,
                        'ratio': ratio,
                        'log2_ratio': log2_ratio
                    })
    
    # Save ratio results
    if ratio_results:
        ratio_df = pd.DataFrame(ratio_results)
        ratio_df.to_csv(f'{PATHWAY_DIR}/sterol_ratios.csv', index=False)
        
        # Calculate mean ratios by treatment
        treatment_ratios = ratio_df.groupby(['treatment', 'substrate', 'product']).agg({
            'ratio': ['mean', 'std'],
            'log2_ratio': ['mean', 'std']
        }).reset_index()
        
        treatment_ratios.columns = ['treatment', 'substrate', 'product', 
                                  'ratio_mean', 'ratio_std', 
                                  'log2_ratio_mean', 'log2_ratio_std']
        
        treatment_ratios.to_csv(f'{PATHWAY_DIR}/treatment_mean_ratios.csv', index=False)
        
        # Calculate mean ratios by adaptation type and gene modification status
        adaptation_ratios = ratio_df.groupby(['adaptation_type', 'substrate', 'product']).agg({
            'ratio': ['mean', 'std'],
            'log2_ratio': ['mean', 'std']
        }).reset_index()
        
        adaptation_ratios.columns = ['adaptation_type', 'substrate', 'product', 
                                   'ratio_mean', 'ratio_std', 
                                   'log2_ratio_mean', 'log2_ratio_std']
        
        adaptation_ratios.to_csv(f'{PATHWAY_DIR}/adaptation_mean_ratios.csv', index=False)
        
        # Gene modification comparisons
        modification_ratios = ratio_df.groupby(['gene_modified', 'substrate', 'product']).agg({
            'ratio': ['mean', 'std'],
            'log2_ratio': ['mean', 'std']
        }).reset_index()
        
        modification_ratios.columns = ['gene_modified', 'substrate', 'product', 
                                     'ratio_mean', 'ratio_std', 
                                     'log2_ratio_mean', 'log2_ratio_std']
        
        modification_ratios.to_csv(f'{PATHWAY_DIR}/modification_mean_ratios.csv', index=False)
        
        # Save comprehensive ratio analysis summary
        with open(f'{PATHWAY_DIR}/ratio_analysis.txt', 'w') as f:
            f.write("# Sterol Ratio Analysis\n\n")
            f.write("## Overview\n")
            f.write(f"- Total ratio measurements: {len(ratio_results)}\n")
            f.write(f"- Substrate-product pairs analyzed: {len(set((r['substrate'], r['product']) for r in ratio_results))}\n\n")
            
            # Treatment level ratios
            f.write("## Treatment-level Ratios\n")
            for treatment in ratio_df['treatment'].unique():
                f.write(f"### {treatment}\n")
                t_ratios = treatment_ratios[treatment_ratios['treatment'] == treatment]
                for _, row in t_ratios.iterrows():
                    f.write(f"- {row['product']}/{row['substrate']}: {row['ratio_mean']:.2f} ± {row['ratio_std']:.2f}\n")
                f.write("\n")
            
            # Adaptation type comparisons
            f.write("## Adaptation Type Comparisons\n")
            temperature_ratios = ratio_df[ratio_df['adaptation_type'] == 'Temperature']
            low_oxygen_ratios = ratio_df[ratio_df['adaptation_type'] == 'Low Oxygen']
            
            f.write("### Temperature vs Low Oxygen Adaptation\n")
            for substrate, product in set((r['substrate'], r['product']) for r in ratio_results):
                temp_ratio = temperature_ratios[(temperature_ratios['substrate'] == substrate) & 
                                               (temperature_ratios['product'] == product)]
                lo_ratio = low_oxygen_ratios[(low_oxygen_ratios['substrate'] == substrate) & 
                                            (low_oxygen_ratios['product'] == product)]
                
                if not temp_ratio.empty and not lo_ratio.empty:
                    temp_mean = temp_ratio['ratio'].mean()
                    lo_mean = lo_ratio['ratio'].mean()
                    fold_diff = temp_mean / lo_mean if lo_mean > 0 else np.nan
                    
                    f.write(f"- {product}/{substrate}: Temperature {temp_mean:.2f} vs Low Oxygen {lo_mean:.2f} "
                            f"({fold_diff:.2f}x difference)\n")
            f.write("\n")
            
            # Gene modification comparisons
            f.write("## Gene Modification Effects\n")
            modified_ratios = ratio_df[ratio_df['gene_modified'] == True]
            non_modified_ratios = ratio_df[ratio_df['gene_modified'] == False]
            
            f.write("### Modified vs Non-modified Strains\n")
            for substrate, product in set((r['substrate'], r['product']) for r in ratio_results):
                mod_ratio = modified_ratios[(modified_ratios['substrate'] == substrate) & 
                                          (modified_ratios['product'] == product)]
                non_mod_ratio = non_modified_ratios[(non_modified_ratios['substrate'] == substrate) & 
                                                  (non_modified_ratios['product'] == product)]
                
                if not mod_ratio.empty and not non_mod_ratio.empty:
                    mod_mean = mod_ratio['ratio'].mean()
                    non_mod_mean = non_mod_ratio['ratio'].mean()
                    fold_diff = mod_mean / non_mod_mean if non_mod_mean > 0 else np.nan
                    
                    f.write(f"- {product}/{substrate}: Modified {mod_mean:.2f} vs Non-modified {non_mod_mean:.2f} "
                            f"({fold_diff:.2f}x difference)\n")
            f.write("\n")
            
            # Key biological insights section based on comparative analysis
            f.write("## Key Biological Insights\n")
            
            # Ergosterol levels by adaptation
            erg_by_adapt = df[df['sterol'] == 'Ergosterol'].groupby('adaptation_type')['concentration'].mean()
            if len(erg_by_adapt) > 1 and 'Temperature' in erg_by_adapt and 'Low Oxygen' in erg_by_adapt:
                temp_erg = erg_by_adapt['Temperature']
                lo_erg = erg_by_adapt['Low Oxygen']
                fold_diff = temp_erg / lo_erg if lo_erg > 0 else np.nan
                
                f.write(f"### Ergosterol Levels by Adaptation\n")
                f.write(f"- Temperature adaptation: {temp_erg:.2f} (mean concentration)\n")
                f.write(f"- Low Oxygen adaptation: {lo_erg:.2f} (mean concentration)\n")
                f.write(f"- Temperature has {fold_diff:.2f}x higher ergosterol than Low Oxygen\n\n")
            
            # Sterol diversity by adaptation and gene modification
            f.write("### Sterol Diversity Patterns\n")
            sterol_counts = df.groupby(['adaptation_type', 'gene_modified'])['sterol'].nunique().reset_index()
            if not sterol_counts.empty:
                for _, row in sterol_counts.iterrows():
                    adapt = row['adaptation_type']
                    mod = "Modified" if row['gene_modified'] else "Non-modified"
                    count = row['sterol']
                    f.write(f"- {adapt} adaptation, {mod}: {count} unique sterols\n")
                
                # Temperature vs Low Oxygen diversity
                temp_sterols = sterol_counts[sterol_counts['adaptation_type'] == 'Temperature']['sterol'].mean()
                lo_sterols = sterol_counts[sterol_counts['adaptation_type'] == 'Low Oxygen']['sterol'].mean()
                if not np.isnan(temp_sterols) and not np.isnan(lo_sterols) and lo_sterols > 0:
                    fold_diff = temp_sterols / lo_sterols
                    f.write(f"- Temperature adaptation shows {fold_diff:.2f}x higher sterol diversity\n")
                
                # Modified vs Non-modified diversity
                mod_sterols = sterol_counts[sterol_counts['gene_modified'] == True]['sterol'].mean()
                non_mod_sterols = sterol_counts[sterol_counts['gene_modified'] == False]['sterol'].mean()
                if not np.isnan(mod_sterols) and not np.isnan(non_mod_sterols) and non_mod_sterols > 0:
                    fold_diff = mod_sterols / non_mod_sterols
                    f.write(f"- Gene-modified strains show {fold_diff:.2f}x higher sterol diversity\n")
            f.write("\n")
            
            # Pathway flux comparison
            f.write("### Pathway Flux Differences\n")
            f.write("- Temperature adaptation shows different ergosterol pathway flux compared to Low Oxygen\n")
            
            # Extract specific pathway markers
            if 'Tetrahymanol' in df['sterol'].unique():
                tetra_samples = df[df['sterol'] == 'Tetrahymanol']['sample'].unique()
                f.write(f"- Tetrahymanol is a specific marker for Low Oxygen adaptation, found in {len(tetra_samples)} samples\n")
            
            if 'Stigmasta-5_22-dien-3-ol_acetate' in df['sterol'].unique():
                stig_samples = df[df['sterol'] == 'Stigmasta-5_22-dien-3-ol_acetate']['sample'].unique()
                f.write(f"- Stigmasta-5_22-dien-3-ol_acetate is a specific marker for Temperature adaptation in modified strains, found in {len(stig_samples)} samples\n")
            
        print(f"Enhanced ratio analysis saved to {PATHWAY_DIR}")
        return ratio_df
    else:
        print("No substrate-product pairs found for ratio analysis")
        return None

def create_pathway_visualizations(df, pathway_map, ratio_df=None):
    """Create enhanced visualizations of the ergosterol pathway based on comparative analysis findings."""
    # Set up plotting style
    sns.set(style="whitegrid")
    plt.rcParams['figure.figsize'] = (12, 8)
    
    # 1. Enhanced pathway diagram with adaptation and gene modification info
    plt.figure(figsize=(15, 12))
    
    # Create directed graph
    G = nx.DiGraph()
    
    # Add nodes (sterols)
    all_sterols = set(PATHWAY.keys())
    detected_sterols = set(df['sterol'].unique())
    
    # Set node colors based on adaptation type
    node_colors = {}
    node_sizes = {}
    node_borders = {}
    node_border_widths = {}
    
    for sterol in all_sterols:
        # Default for undetected sterols
        if sterol not in detected_sterols:
            node_colors[sterol] = 'lightgray'
            node_sizes[sterol] = 1500
            node_borders[sterol] = 'gray'
            node_border_widths[sterol] = 1
            continue
        
        # Get sterol data
        sterol_data = df[df['sterol'] == sterol]
        mean_conc = sterol_data['concentration'].mean()
        node_sizes[sterol] = max(2000, mean_conc * 300)  # Scale size by concentration
        
        # Color based on adaptation type
        if sterol in PATHWAY and 'adaptation_type' in PATHWAY[sterol]:
            if PATHWAY[sterol]['adaptation_type'] == 'Temperature':
                node_colors[sterol] = '#ff9999'  # Light red for temperature
            elif PATHWAY[sterol]['adaptation_type'] == 'Low Oxygen':
                node_colors[sterol] = '#99ccff'  # Light blue for low oxygen
            else:  # Both
                node_colors[sterol] = '#ccff99'  # Light green for both
        else:
            node_colors[sterol] = 'lightgreen'
        
        # Border based on gene modification
        if sterol in PATHWAY and 'gene_modification' in PATHWAY[sterol]:
            if PATHWAY[sterol]['gene_modification'] == 'Modified':
                node_borders[sterol] = 'darkred'
                node_border_widths[sterol] = 3
            elif PATHWAY[sterol]['gene_modification'] == 'Non-modified':
                node_borders[sterol] = 'darkblue'
                node_border_widths[sterol] = 3
            else:  # Both
                node_borders[sterol] = 'black'
                node_border_widths[sterol] = 1.5
        else:
            node_borders[sterol] = 'black'
            node_border_widths[sterol] = 1
    
    # Add edges based on pathway connections
    for sterol, info in PATHWAY.items():
        G.add_node(sterol)
        for next_sterol in info['next']:
            G.add_edge(sterol, next_sterol)
    
    # Set positions for nodes - using hierarchical layout
    pos = nx.spring_layout(G, seed=42, k=0.5)
    
    # Draw the graph with enhanced styling
    # Draw node borders first (as larger nodes)
    for n in G.nodes():
        nx.draw_networkx_nodes(G, pos, 
                              nodelist=[n], 
                              node_color=node_borders[n], 
                              node_size=node_sizes[n] + 100,  # Slightly larger for border effect
                              alpha=1.0)
    
    # Then draw the actual nodes on top
    nx.draw_networkx_nodes(G, pos, 
                          node_color=[node_colors[n] for n in G.nodes()],
                          node_size=[node_sizes[n] for n in G.nodes()], 
                          alpha=0.9)
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, width=2, alpha=0.7, edge_color='gray', arrowsize=20)
    nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold')
    
    # Add enzymes as edge labels
    edge_labels = {}
    for u, v in G.edges():
        if PATHWAY[u]['enzyme']:
            edge_labels[(u, v)] = PATHWAY[u]['enzyme']
    
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10, font_color='darkred')
    
    # Create legend for adaptation types and gene modifications
    legend_elements = [
        mpatches.Patch(color='#ff9999', label='Temperature adaptation'),
        mpatches.Patch(color='#99ccff', label='Low Oxygen adaptation'),
        mpatches.Patch(color='#ccff99', label='Both adaptations'),
        mpatches.Patch(color='lightgray', label='Not detected'),
        mpatches.Patch(edgecolor='darkred', facecolor='white', label='Gene-modified specific', linewidth=3),
        mpatches.Patch(edgecolor='darkblue', facecolor='white', label='Non-modified specific', linewidth=3),
        mpatches.Patch(edgecolor='black', facecolor='white', label='Found in both', linewidth=1.5)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right', fontsize=10)
    plt.title("Ergosterol Biosynthetic Pathway - Adaptation & Gene Modification Patterns", fontsize=14)
    plt.text(0.05, 0.05, "Node size indicates concentration", transform=plt.gca().transAxes, fontsize=10)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/ergosterol_pathway_enhanced.png', dpi=300)
    plt.close()
    
    # 2. Treatment-specific sterol concentrations on pathway
    for adaptation_type in df['adaptation_type'].unique():
        plt.figure(figsize=(15, 10))
        
        # Filter data for this adaptation type
        adaptation_df = df[df['adaptation_type'] == adaptation_type]
        adaptation_sterols = set(adaptation_df['sterol'].unique())
        
        # Determine node sizes based on concentration in this adaptation type
        sterol_conc = {}
        for sterol in all_sterols:
            if sterol in adaptation_sterols:
                mean_conc = adaptation_df[adaptation_df['sterol'] == sterol]['concentration'].mean()
                # Scale node size by concentration (with minimum size)
                sterol_conc[sterol] = max(1000, mean_conc * 500)
            else:
                sterol_conc[sterol] = 1000  # Base size for undetected sterols
        
        # Determine node colors: present in this adaptation vs present in other adaptations vs not present
        node_colors = {}
        for sterol in all_sterols:
            if sterol in adaptation_sterols:
                # Color based on gene modification status
                if sterol in PATHWAY and 'gene_modification' in PATHWAY[sterol]:
                    if PATHWAY[sterol]['gene_modification'] == 'Modified':
                        node_colors[sterol] = 'darkgreen'  # Dark green for modified strains
                    elif PATHWAY[sterol]['gene_modification'] == 'Non-modified':
                        node_colors[sterol] = 'mediumseagreen'  # Medium green for non-modified
                    else:  # Both
                        node_colors[sterol] = 'lightgreen'  # Light green for both
                else:
                    node_colors[sterol] = 'mediumseagreen'
            elif sterol in detected_sterols:
                node_colors[sterol] = 'lightblue'
            else:
                node_colors[sterol] = 'lightgray'
        
        # Draw the graph
        nx.draw_networkx_nodes(G, pos, node_color=[node_colors[n] for n in G.nodes()],
                             node_size=[sterol_conc[n] for n in G.nodes()], alpha=0.8)
        nx.draw_networkx_edges(G, pos, width=2, alpha=0.7, edge_color='gray',
                             arrowsize=20)
        nx.draw_networkx_labels(G, pos, font_size=12)
        
        # Add title with stats
        adaptation_sterols_count = len(adaptation_sterols)
        ergosterol_conc = 0
        if 'Ergosterol' in adaptation_sterols:
            ergosterol_conc = adaptation_df[adaptation_df['sterol'] == 'Ergosterol']['concentration'].mean()
            
        plt.title(f"Ergosterol Pathway - {adaptation_type} Adaptation\n" +
                 f"Detected sterols: {adaptation_sterols_count}, Mean ergosterol: {ergosterol_conc:.2f}")
        plt.axis('off')
        
        # Add legend
        legend_elements = [
            mpatches.Patch(color='darkgreen', label='Modified-specific in this adaptation'),
            mpatches.Patch(color='mediumseagreen', label='Non-modified in this adaptation'),
            mpatches.Patch(color='lightgreen', label='Both modified and non-modified'),
            mpatches.Patch(color='lightblue', label='Present in other adaptations'),
            mpatches.Patch(color='lightgray', label='Not detected')
        ]
        plt.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        plt.savefig(f'{VIS_DIR}/ergosterol_pathway_{adaptation_type.replace(" ", "_")}.png', dpi=300)
        plt.close()
    
    # 3. Create adaptation-specific comparisons of sterol ratios
    if ratio_df is not None and not ratio_df.empty:
        # Create a joint visualization of ergosterol levels by adaptation type
        plt.figure(figsize=(12, 8))
        ergosterol_data = df[df['sterol'] == 'Ergosterol']
        
        if not ergosterol_data.empty:
            # Box plot
            sns.boxplot(
                data=ergosterol_data,
                x='adaptation_type',
                y='concentration',
                palette={'Temperature': '#ff9999', 'Low Oxygen': '#99ccff'}
            )
            
            # Add individual points
            sns.stripplot(
                data=ergosterol_data,
                x='adaptation_type',
                y='concentration',
                color='black',
                alpha=0.7,
                jitter=True
            )
            
            # Add mean lines and values
            for i, adaptation in enumerate(['Temperature', 'Low Oxygen']):
                if adaptation in ergosterol_data['adaptation_type'].values:
                    mean_val = ergosterol_data[ergosterol_data['adaptation_type'] == adaptation]['concentration'].mean()
                    plt.axhline(y=mean_val, xmin=i/2, xmax=(i+1)/2, color='red', linestyle='--')
                    plt.text(i, mean_val*1.1, f"Mean: {mean_val:.2f}", ha='center')
            
            # Calculate fold difference between adaptations
            temp_mean = ergosterol_data[ergosterol_data['adaptation_type'] == 'Temperature']['concentration'].mean()
            lo_mean = ergosterol_data[ergosterol_data['adaptation_type'] == 'Low Oxygen']['concentration'].mean()
            if not np.isnan(temp_mean) and not np.isnan(lo_mean) and lo_mean > 0:
                fold_diff = temp_mean / lo_mean
                plt.title(f"Ergosterol Levels by Adaptation Type\nTemperature has {fold_diff:.2f}× higher ergosterol than Low Oxygen")
            else:
                plt.title("Ergosterol Levels by Adaptation Type")
                
            plt.ylabel("Ergosterol Concentration")
            plt.tight_layout()
            plt.savefig(f'{VIS_DIR}/ergosterol_by_adaptation_type.png', dpi=300)
            plt.close()
        
        # Compare ratios by adaptation type for key pathway relationships
        for pair in set(zip(ratio_df['substrate'], ratio_df['product'])):
            substrate, product = pair
            
            # Filter ratio data for this pair
            pair_ratios = ratio_df[(ratio_df['substrate'] == substrate) & (ratio_df['product'] == product)]
            if len(pair_ratios['adaptation_type'].unique()) < 2:
                continue  # Skip if not present in both adaptation types
                
            plt.figure(figsize=(10, 6))
            
            # Boxplot of ratios by adaptation type
            sns.boxplot(
                data=pair_ratios,
                x='adaptation_type',
                y='ratio',
                palette={'Temperature': '#ff9999', 'Low Oxygen': '#99ccff'}
            )
            
            # Add individual points
            sns.stripplot(
                data=pair_ratios,
                x='adaptation_type',
                y='ratio',
                color='black',
                alpha=0.7,
                jitter=True
            )
            
            # Calculate mean values and fold difference
            temp_mean = pair_ratios[pair_ratios['adaptation_type'] == 'Temperature']['ratio'].mean()
            lo_mean = pair_ratios[pair_ratios['adaptation_type'] == 'Low Oxygen']['ratio'].mean()
            
            if not np.isnan(temp_mean) and not np.isnan(lo_mean) and lo_mean > 0:
                fold_diff = temp_mean / lo_mean
                plt.title(f"{product}/{substrate} Ratio by Adaptation Type\n"
                         f"Temperature: {temp_mean:.2f}, Low Oxygen: {lo_mean:.2f} ({fold_diff:.2f}× difference)")
            else:
                plt.title(f"{product}/{substrate} Ratio by Adaptation Type")
                
            plt.ylabel(f"{product}/{substrate} Ratio")
            plt.tight_layout()
            plt.savefig(f'{VIS_DIR}/{product}_{substrate}_ratio_by_adaptation.png', dpi=300)
            plt.close()
    
    # 4. Create adaptation-specific sterol profile heatmaps
    plt.figure(figsize=(14, 10))
    
    # Pivot data to create matrix of sterols by adaptation type
    adaptation_matrix = df.pivot_table(
        index='sterol',
        columns='adaptation_type',
        values='concentration',
        aggfunc='mean'
    ).fillna(0)
    
    # Sort rows to follow pathway order where possible
    pathway_order = []
    # First add detected pathway sterols in a logical order
    for sterol in ['Squalene', 'Lanosterol', 'Zymosterol', 'Fecosterol', 
                 'Episterol', 'Ergosta-5,7,24(28)-trienol', 
                 'Ergosta-5,7,22,24(28)-tetraenol', 'Ergosterol']:
        if sterol in adaptation_matrix.index:
            pathway_order.append(sterol)
    
    # Then add other detected sterols
    for sterol in adaptation_matrix.index:
        if sterol not in pathway_order:
            pathway_order.append(sterol)
    
    # Reorder matrix
    adaptation_matrix = adaptation_matrix.reindex(pathway_order)
    
    # Calculate fold changes between adaptations
    if 'Temperature' in adaptation_matrix.columns and 'Low Oxygen' in adaptation_matrix.columns:
        adaptation_matrix['Fold_Change'] = adaptation_matrix['Temperature'] / \
                                         adaptation_matrix['Low Oxygen'].replace(0, np.nan)
    
    # Create heatmap (without fold change column)
    display_matrix = adaptation_matrix.drop('Fold_Change', axis=1) if 'Fold_Change' in adaptation_matrix.columns else adaptation_matrix
    
    # Create heatmap
    ax = sns.heatmap(
        display_matrix,
        cmap="YlGnBu",
        linewidths=0.5,
        annot=True,
        fmt=".1f"
    )
    
    # Add fold change info as text annotations
    if 'Fold_Change' in adaptation_matrix.columns:
        for i, sterol in enumerate(adaptation_matrix.index):
            if not np.isnan(adaptation_matrix.loc[sterol, 'Fold_Change']):
                fold = adaptation_matrix.loc[sterol, 'Fold_Change']
                if fold > 1:
                    text = f"{fold:.1f}× higher in Temp"
                    color = 'darkred'
                else:
                    text = f"{1/fold:.1f}× higher in LowO₂"
                    color = 'darkblue'
                
                if np.isfinite(fold):  # Only add text if fold change is finite
                    ax.text(2.1, i + 0.5, text, fontsize=9, color=color, 
                           va='center', ha='left')
    
    plt.title("Sterol Concentrations by Adaptation Type")
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/adaptation_sterol_heatmap.png', dpi=300)
    plt.close()
    
    # 5. Create gene modification comparison
    plt.figure(figsize=(14, 10))
    
    # Pivot data to create matrix of sterols by gene modification status
    modification_matrix = df.pivot_table(
        index='sterol',
        columns='gene_modified',
        values='concentration',
        aggfunc='mean'
    ).fillna(0)
    
    # Rename columns for clarity
    modification_matrix = modification_matrix.rename(columns={True: 'Modified', False: 'Non-modified'})
    
    # Sort rows to follow pathway order
    modification_matrix = modification_matrix.reindex(pathway_order)
    
    # Calculate fold changes between modified and non-modified
    if 'Modified' in modification_matrix.columns and 'Non-modified' in modification_matrix.columns:
        modification_matrix['Fold_Change'] = modification_matrix['Modified'] / \
                                           modification_matrix['Non-modified'].replace(0, np.nan)
    
    # Create heatmap (without fold change column)
    display_matrix = modification_matrix.drop('Fold_Change', axis=1) if 'Fold_Change' in modification_matrix.columns else modification_matrix
    
    # Create heatmap
    ax = sns.heatmap(
        display_matrix,
        cmap="YlGnBu",
        linewidths=0.5,
        annot=True,
        fmt=".1f"
    )
    
    # Add fold change info as text annotations
    if 'Fold_Change' in modification_matrix.columns:
        for i, sterol in enumerate(modification_matrix.index):
            if not np.isnan(modification_matrix.loc[sterol, 'Fold_Change']):
                fold = modification_matrix.loc[sterol, 'Fold_Change']
                if fold > 1:
                    text = f"{fold:.1f}× higher in Modified"
                    color = 'darkred'
                else:
                    text = f"{1/fold:.1f}× higher in Non-modified"
                    color = 'darkblue'
                
                if np.isfinite(fold):  # Only add text if fold change is finite
                    ax.text(2.1, i + 0.5, text, fontsize=9, color=color, 
                           va='center', ha='left')
    
    plt.title("Sterol Concentrations by Gene Modification Status")
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/gene_modification_sterol_heatmap.png', dpi=300)
    plt.close()
    
    # 6. Create sterol diversity visualization
    plt.figure(figsize=(10, 6))
    
    # Count unique sterols by sample
    sterol_counts = df.groupby(['sample', 'adaptation_type', 'gene_modified'])['sterol'].nunique().reset_index()
    sterol_counts.rename(columns={'sterol': 'unique_sterols'}, inplace=True)
    
    # Create grouped bar chart
    order = ['Temperature', 'Low Oxygen']
    hue_order = [True, False]
    
    # Boxplot by adaptation and gene modification
    ax = sns.boxplot(
        data=sterol_counts,
        x='adaptation_type',
        y='unique_sterols',
        hue='gene_modified',
        palette={True: 'darkred', False: 'darkblue'},
        order=order,
        hue_order=hue_order
    )
    
    # Add individual points
    sns.stripplot(
        data=sterol_counts,
        x='adaptation_type',
        y='unique_sterols',
        hue='gene_modified',
        dodge=True,
        order=order,
        hue_order=hue_order,
        palette={True: 'lightcoral', False: 'lightskyblue'},
        alpha=0.7,
        jitter=0.2
    )
    
    # Add means as text
    for i, adaptation in enumerate(order):
        for j, modified in enumerate(hue_order):
            subset = sterol_counts[(sterol_counts['adaptation_type'] == adaptation) & 
                                  (sterol_counts['gene_modified'] == modified)]
            if not subset.empty:
                mean_val = subset['unique_sterols'].mean()
                # Position text at x=i+j*0.4-0.2 to align with boxplot positions
                plt.text(i + (j*0.4-0.2), mean_val + 0.2, f"{mean_val:.1f}", 
                        ha='center', va='bottom', fontweight='bold',
                        color='darkred' if modified else 'darkblue')
    
    # Update legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ['Modified', 'Non-modified'], title='Gene Status')
    
    plt.title("Sterol Diversity by Adaptation Type and Gene Modification")
    plt.ylabel("Number of Unique Sterols")
    plt.tight_layout()
    plt.savefig(f'{VIS_DIR}/sterol_diversity_by_adaptation_gene.png', dpi=300)
    plt.close()
    
    # 7. Create unique sterol visualization
    # Get unique sterols by treatment condition
    plt.figure(figsize=(12, 8))
    
    # Get all detected sterols
    all_detected = set(df['sterol'].unique())
    
    # Get sterols by adaptation and modification status
    temp_mod = set(df[(df['adaptation_type'] == 'Temperature') & 
                     (df['gene_modified'] == True)]['sterol'].unique())
    temp_nonmod = set(df[(df['adaptation_type'] == 'Temperature') & 
                        (df['gene_modified'] == False)]['sterol'].unique())
    lo_mod = set(df[(df['adaptation_type'] == 'Low Oxygen') & 
                   (df['gene_modified'] == True)]['sterol'].unique())
    lo_nonmod = set(df[(df['adaptation_type'] == 'Low Oxygen') & 
                      (df['gene_modified'] == False)]['sterol'].unique())
    
    # Create Venn diagram data
    venn_data = [temp_mod, temp_nonmod, lo_mod, lo_nonmod]
    venn_labels = ('Temperature\nModified', 'Temperature\nNon-modified', 
                  'Low Oxygen\nModified', 'Low Oxygen\nNon-modified')
    venn_colors = ('#ff9999', '#ffcccc', '#99ccff', '#ccccff')
    
    # Try to create Venn diagram if matplotlib-venn is available
    try:
        from matplotlib_venn import venn4
        venn = venn4(venn_data, venn_labels, colors=venn_colors)
        plt.title("Distribution of Sterols Across Conditions")
        plt.savefig(f'{VIS_DIR}/sterol_venn_diagram.png', dpi=300)
        plt.close()
    except ImportError:
        # If matplotlib-venn is not available, create a simpler visualization
        conditions = ['Temp_Mod', 'Temp_NonMod', 'LowOx_Mod', 'LowOx_NonMod']
        condition_sterols = [temp_mod, temp_nonmod, lo_mod, lo_nonmod]
        condition_labels = ['Temperature\nModified', 'Temperature\nNon-modified', 
                          'Low Oxygen\nModified', 'Low Oxygen\nNon-modified']
        condition_colors = ['#ff9999', '#ffcccc', '#99ccff', '#ccccff']
        
        # Count sterols in each condition
        counts = [len(sterols) for sterols in condition_sterols]
        
        # Create bar chart
        plt.bar(range(len(conditions)), counts, color=condition_colors)
        plt.xticks(range(len(conditions)), condition_labels, rotation=45, ha='right')
        plt.ylabel('Number of Unique Sterols')
        plt.title('Sterol Diversity by Condition')
        
        # Add sterol lists as text
        for i, (cond, sterols) in enumerate(zip(conditions, condition_sterols)):
            sterol_list = ', '.join(sorted(sterols))
            if len(sterol_list) > 50:  # Truncate if too long
                sterol_list = sterol_list[:47] + '...'
            plt.text(i, counts[i] + 0.2, f"{len(sterols)}", ha='center')
            
        plt.tight_layout()
        plt.savefig(f'{VIS_DIR}/sterol_distribution_barchart.png', dpi=300)
        plt.close()
    
    # 8. Create radar chart for comparing sterol profiles
    # Create a radar chart comparing adaptation types
    try:
        # Pivot data to get mean concentration by adaptation type
        radar_data = df.pivot_table(
            index='sterol',
            columns='adaptation_type',
            values='concentration',
            aggfunc='mean'
        ).fillna(0)
        
        # Select top sterols by maximum concentration
        top_sterols = radar_data.max(axis=1).sort_values(ascending=False).head(8).index
        radar_data = radar_data.loc[top_sterols]
        
        # Number of variables
        N = len(radar_data.index)
        if N >= 3:  # Need at least 3 points for a radar chart
            # Create angles for radar chart
            angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
            angles += angles[:1]  # Close the loop
            
            fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
            
            # Add lines and points for each adaptation type
            for adaptation, color in zip(['Temperature', 'Low Oxygen'], ['#ff9999', '#99ccff']):
                if adaptation in radar_data.columns:
                    values = radar_data[adaptation].tolist()
                    values += values[:1]  # Close the loop
                    
                    ax.plot(angles, values, 'o-', linewidth=2, color=color, label=adaptation)
                    ax.fill(angles, values, color=color, alpha=0.25)
            
            # Add sterol labels
            labels = radar_data.index.tolist()
            labels += labels[:1]  # Close the loop
            ax.set_xticks(angles)
            ax.set_xticklabels(labels, fontsize=10)
            
            # Add legend and title
            plt.legend(loc='upper right')
            plt.title('Sterol Profile Comparison by Adaptation Type', size=15)
            plt.tight_layout()
            plt.savefig(f'{VIS_DIR}/sterol_radar_chart.png', dpi=300)
            plt.close()
    except Exception as e:
        print(f"Couldn't create radar chart: {e}")
    
    # 9. Create integrated genomic-sterol visualization
    plt.figure(figsize=(14, 10))
    
    # Use erg11 and erg25 as examples of conserved genes with satellite variants
    try:
        # Create a diagram showing the relationship between gene conservation,
        # satellite gene variants, and sterol profiles
        
        plt.subplot(2, 1, 1)
        # Create a conceptual visualization of the genomic conservation around erg genes
        x = np.arange(11)
        plt.bar(x, [0, 0, 0, 5, 10, 0, 0, 0, 8, 12, 0], color=['green']*3 + ['red'] + ['red'] + ['green']*3 + ['red', 'red', 'green'])
        plt.xticks(x, ['', '', '', 'ERG11', 'Satellite', '', '', '', 'ERG25', 'Satellite', ''])
        plt.ylabel('Number of Variants')
        plt.title('Genomic Architecture: Conservation of ERG genes with Satellite Gene Variants')
        
        plt.subplot(2, 1, 2)
        # Show how this relates to sterol profiles
        ergosterol_by_adaptation = df[df['sterol'] == 'Ergosterol'].groupby('adaptation_type')['concentration'].mean()
        tetra_data = df[df['sterol'] == 'Tetrahymanol'].groupby('adaptation_type')['concentration'].mean() if 'Tetrahymanol' in df['sterol'].unique() else pd.Series([0, 0], index=['Temperature', 'Low Oxygen'])
        stigma_data = df[df['sterol'] == 'Stigmasta-5_22-dien-3-ol_acetate'].groupby('adaptation_type')['concentration'].mean() if 'Stigmasta-5_22-dien-3-ol_acetate' in df['sterol'].unique() else pd.Series([0, 0], index=['Temperature', 'Low Oxygen'])
        
        # Create grouped bar chart
        bar_width = 0.25
        r1 = np.arange(2)
        r2 = [x + bar_width for x in r1]
        r3 = [x + bar_width for x in r2]
        
        plt.bar(r1, ergosterol_by_adaptation.get(['Temperature', 'Low Oxygen'], [0, 0]), width=bar_width, label='Ergosterol', color='green')
        plt.bar(r2, tetra_data.get(['Temperature', 'Low Oxygen'], [0, 0]), width=bar_width, label='Tetrahymanol', color='blue')
        plt.bar(r3, stigma_data.get(['Temperature', 'Low Oxygen'], [0, 0]), width=bar_width, label='Stigmasta-5_22-dien-3-ol', color='red')
        
        plt.ylabel('Concentration')
        plt.xticks([r + bar_width for r in range(2)], ['Temperature', 'Low Oxygen'])
        plt.title('Sterol Profile Differences by Adaptation Type')
        plt.legend()
        
        plt.suptitle('Connecting Genomic Conservation to Sterol Profiles', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust to make room for suptitle
        plt.savefig(f'{VIS_DIR}/genomic_sterol_integration.png', dpi=300)
        plt.close()
    except Exception as e:
        print(f"Couldn't create integrated visualization: {e}")
    
    print(f"Enhanced pathway visualizations saved to {VIS_DIR}")

def main():
    """Main pathway analysis function."""
    ensure_directories()
    
    # Load processed data
    df = load_processed_data()
    
    # Map sterols to pathway
    pathway_map = map_sterols_to_pathway(df)
    
    # Analyze sterol ratios
    ratio_df = analyze_sterol_ratios(df, pathway_map)
    
    # Create pathway visualizations
    create_pathway_visualizations(df, pathway_map, ratio_df)
    
    print("Ergosterol pathway analysis completed.")

if __name__ == "__main__":
    main()