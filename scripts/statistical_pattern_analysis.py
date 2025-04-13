#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr, spearmanr, chi2_contingency
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
import statsmodels.api as sm
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings('ignore')

# Set matplotlib style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "statistical_pattern_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Function to load mutation data from previous analyses
def load_mutation_data():
    """Load mutation data from previous analyses."""
    # Try to load data from mutation_spectrum_analysis directory
    treatments = ['WT', 'STC', 'CAS', 'WTA']
    all_data = []
    
    for treatment in treatments:
        file_path = f"analysis/MSA/mutation_spectrum_analysis/{treatment}_mutations.txt"
        if os.path.exists(file_path):
            # Load data
            data = pd.read_csv(file_path, sep='\t', header=None,
                              names=['CHROM', 'POS', 'REF', 'ALT'])
            data['Treatment'] = treatment
            all_data.append(data)
    
    if all_data:
        combined_data = pd.concat(all_data, ignore_index=True)
        print(f"Loaded {len(combined_data)} mutations across {len(treatments)} treatments")
        return combined_data
    else:
        print("No mutation data found.")
        return None

# Function to load scaffold information
def load_scaffold_info():
    """Load scaffold information from reference genome."""
    scaffold_info = {}
    
    # Try to load scaffold lengths from fai file
    fai_path = "reference/yeast_w303.fasta.fai"
    if os.path.exists(fai_path):
        scaffold_df = pd.read_csv(fai_path, sep='\t', header=None,
                                 names=['Scaffold', 'Length', 'Offset', 'Linebases', 'Linewidth'])
        scaffold_info = dict(zip(scaffold_df['Scaffold'], scaffold_df['Length']))
        print(f"Loaded information for {len(scaffold_info)} scaffolds")
    else:
        print("No scaffold information found.")
    
    return scaffold_info

# Function to load genomic context data if available
def load_context_data():
    """Load genomic context data by parsing available files and reconstructing from reference genome if needed."""
    context_data = []
    
    # Step 1: Parse the genomic_context_summary.txt file for GC content and homopolymer statistics
    summary_file = "genomic_context_results/genomic_context_summary.txt"
    if os.path.exists(summary_file):
        print("Found genomic context summary file, extracting data...")
        
        # Define dictionaries to store treatment-level statistics
        gc_content = {}
        homopolymer_presence = {}
        dinucleotide_presence = {}
        
        with open(summary_file, 'r') as f:
            lines = f.readlines()
            
            # Parse GC content section
            gc_section = False
            for i, line in enumerate(lines):
                if "GC Content Analysis:" in line:
                    gc_section = True
                    continue
                
                if gc_section and "GC content by treatment:" in line:
                    # Next few lines contain treatment-specific GC content
                    j = i + 1
                    while j < len(lines) and lines[j].strip() and ":" in lines[j]:
                        parts = lines[j].strip().split(":")
                        if len(parts) == 2:
                            treatment = parts[0].strip()
                            gc = float(parts[1].strip())
                            gc_content[treatment] = gc
                        j += 1
                    break
            
            # Parse homopolymer presence section
            homo_section = False
            for i, line in enumerate(lines):
                if "Homopolymer presence by treatment:" in line:
                    homo_section = True
                    j = i + 1
                    while j < len(lines) and lines[j].strip() and ":" in lines[j]:
                        parts = lines[j].strip().split(":")
                        if len(parts) == 2:
                            treatment = parts[0].strip()
                            freq = float(parts[1].strip())
                            homopolymer_presence[treatment] = freq
                        j += 1
                    break
            
            # Parse dinucleotide repeat presence section
            dinuc_section = False
            for i, line in enumerate(lines):
                if "Dinucleotide repeat presence by treatment:" in line:
                    dinuc_section = True
                    j = i + 1
                    while j < len(lines) and lines[j].strip() and ":" in lines[j]:
                        parts = lines[j].strip().split(":")
                        if len(parts) == 2:
                            treatment = parts[0].strip()
                            freq = float(parts[1].strip())
                            dinucleotide_presence[treatment] = freq
                        j += 1
                    break
        
        print(f"Extracted GC content data for {len(gc_content)} treatments")
        print(f"Extracted homopolymer data for {len(homopolymer_presence)} treatments")
        print(f"Extracted dinucleotide data for {len(dinucleotide_presence)} treatments")
        
        if gc_content or homopolymer_presence or dinucleotide_presence:
            # Success! We have some context data. Now load the mutation data to combine with it.
            mutation_files = [f for f in os.listdir("mutation_spectrum_analysis") if f.endswith("_mutations.txt")]
            
            for file in mutation_files:
                treatment = file.split('_')[0]
                
                # Read mutation data
                mutations_df = pd.read_csv(
                    os.path.join("mutation_spectrum_analysis", file), 
                    sep='\t', 
                    header=None,
                    names=['CHROM', 'POS', 'REF', 'ALT']
                )
                
                # Add treatment and context information to each mutation
                for _, row in mutations_df.iterrows():
                    mutation_type = f"{row['REF']}>{row['ALT']}"
                    
                    entry = {
                        'Treatment': treatment,
                        'CHROM': row['CHROM'],
                        'POS': row['POS'],
                        'REF': row['REF'],
                        'ALT': row['ALT'],
                        'Mutation_Type': mutation_type
                    }
                    
                    # Add context features if available for this treatment
                    if treatment in gc_content:
                        entry['GC_Content'] = gc_content[treatment]
                    
                    if treatment in homopolymer_presence:
                        entry['Homopolymer_Presence'] = homopolymer_presence[treatment]
                    
                    if treatment in dinucleotide_presence:
                        entry['Dinucleotide_Presence'] = dinucleotide_presence[treatment]
                    
                    context_data.append(entry)
            
            context_df = pd.DataFrame(context_data)
            print(f"Created context dataset with {len(context_df)} entries")
            return context_df
    
    # Step 2: If we couldn't extract from summary file, try an alternative approach
    # Since we know we have the mutation data and reference genome, we can extract context sequences directly
    try:
        print("Attempting to reconstruct context data from mutation data and reference genome...")
        
        # Check if we have the reference genome
        reference_genome = {}
        genome_file = "reference/yeast_w303.fasta"
        
        if os.path.exists(genome_file):
            from Bio import SeqIO
            
            # Load the reference genome
            for record in SeqIO.parse(genome_file, "fasta"):
                reference_genome[record.id] = str(record.seq).upper()
            
            print(f"Loaded {len(reference_genome)} sequences from reference genome")
            
            # Now load mutation data and extract contexts
            mutation_files = [f for f in os.listdir("analysis/MSA/mutation_spectrum_analysis") if f.endswith("_mutations.txt")]
            
            for file in mutation_files:
                treatment = file.split('_')[0]
                
                # Read mutation data
                mutations_df = pd.read_csv(
                    os.path.join("analysis/MSA/mutation_spectrum_analysis", file), 
                    sep='\t', 
                    header=None,
                    names=['CHROM', 'POS', 'REF', 'ALT']
                )
                
                # For each mutation, extract context and calculate features
                for _, row in mutations_df.iterrows():
                    chrom = row['CHROM']
                    pos = int(row['POS'])
                    ref = row['REF']
                    alt = row['ALT']
                    
                    if chrom in reference_genome:
                        seq = reference_genome[chrom]
                        
                        # Extract ±100bp context if possible
                        context_size = 100
                        start = max(0, pos - context_size - 1)
                        end = min(len(seq), pos + context_size)
                        
                        if start < end:
                            context = seq[start:end]
                            
                            # Calculate GC content
                            gc_count = sum(1 for base in context if base in 'GC')
                            gc_content = gc_count / len(context) if len(context) > 0 else 0
                            
                            # Check for homopolymers (runs of 3+ identical bases)
                            has_homopolymer = any(base * 3 in context for base in 'ACGT')
                            
                            # Check for dinucleotide repeats
                            has_dinucleotide = False
                            for i in range(len(context) - 3):
                                if context[i:i+2] == context[i+2:i+4]:
                                    has_dinucleotide = True
                                    break
                            
                            # Add entry to context data
                            entry = {
                                'Treatment': treatment,
                                'CHROM': chrom,
                                'POS': pos,
                                'REF': ref,
                                'ALT': alt,
                                'Mutation_Type': f"{ref}>{alt}",
                                'GC_Content': gc_content,
                                'Homopolymer_Presence': 1 if has_homopolymer else 0,
                                'Dinucleotide_Presence': 1 if has_dinucleotide else 0
                            }
                            
                            context_data.append(entry)
            
            if context_data:
                context_df = pd.DataFrame(context_data)
                print(f"Successfully reconstructed context data with {len(context_df)} entries")
                return context_df
    except Exception as e:
        print(f"Error reconstructing context data: {e}")
    
    # Step 3: If all else fails, use mutation data alone
    try:
        print("Attempting to create basic context data from mutation data alone...")
        
        # Load mutation data files
        mutation_files = [f for f in os.listdir("analysis/MSA/mutation_spectrum_analysis") if f.endswith("_mutations.txt")]
        
        for file in mutation_files:
            treatment = file.split('_')[0]
            
            # Read mutation data
            mutations_df = pd.read_csv(
                os.path.join("mutation_spectrum_analysis", file), 
                sep='\t', 
                header=None,
                names=['CHROM', 'POS', 'REF', 'ALT']
            )
            
            # For each mutation, create a basic entry
            for _, row in mutations_df.iterrows():
                entry = {
                    'Treatment': treatment,
                    'CHROM': row['CHROM'],
                    'POS': row['POS'],
                    'REF': row['REF'],
                    'ALT': row['ALT'],
                    'Mutation_Type': f"{row['REF']}>{row['ALT']}"
                }
                
                context_data.append(entry)
        
        if context_data:
            context_df = pd.DataFrame(context_data)
            print(f"Created basic context data with {len(context_df)} entries")
            return context_df
    except Exception as e:
        print(f"Error creating basic context data: {e}")
    
    print("Could not load or reconstruct context data. Analysis will proceed with limited features.")
    return None

# Function to integrate all data sources
def integrate_data(mutation_data, scaffold_info, context_data=None):
    """Integrate data from various sources for statistical analysis."""
    if mutation_data is None:
        print("Cannot integrate data: mutation data is missing")
        return None
    
    # Create a copy of mutation data
    integrated_data = mutation_data.copy()
    
    # Add scaffold length information if available
    if scaffold_info:
        integrated_data['Scaffold_Length'] = integrated_data['CHROM'].map(scaffold_info)
    
    # Prepare basic variant features
    integrated_data['Mutation_Type'] = integrated_data.apply(
        lambda row: f"{row['REF']}>{row['ALT']}", axis=1)
    
    # Standardize mutation types to 6 basic categories
    def standardize_mutation(ref, alt):
        """Standardize mutation to one of 6 basic types."""
        if ref in ['C', 'T']:
            return f"{ref}>{alt}"
        else:
            # Convert purine to pyrimidine
            purine_to_pyrimidine = {'A': 'T', 'G': 'C'}
            pyrimidine_to_purine = {'T': 'A', 'C': 'G'}
            std_ref = purine_to_pyrimidine.get(ref, ref)
            std_alt = purine_to_pyrimidine.get(alt, pyrimidine_to_purine.get(alt, alt))
            return f"{std_ref}>{std_alt}"
    
    integrated_data['Std_Mutation'] = integrated_data.apply(
        lambda row: standardize_mutation(row['REF'], row['ALT']), axis=1)
    
    # Create transition/transversion classification
    transitions = ['C>T', 'T>C', 'G>A', 'A>G']
    integrated_data['Mutation_Class'] = integrated_data['Mutation_Type'].apply(
        lambda x: 'Transition' if x in transitions else 'Transversion')
    
    # Add scaffold-level statistics
    scaffold_counts = integrated_data.groupby('CHROM').size().to_dict()
    integrated_data['Scaffold_Variant_Count'] = integrated_data['CHROM'].map(scaffold_counts)
    
    # Calculate variant density if scaffold length is available
    if 'Scaffold_Length' in integrated_data.columns:
        integrated_data['Variant_Density'] = integrated_data.apply(
            lambda row: (row['Scaffold_Variant_Count'] * 1000 / row['Scaffold_Length']) 
            if row['Scaffold_Length'] > 0 else 0, axis=1)
    
    # Merge with context data if available
    if context_data is not None:
        # Implement the merge based on the structure of context_data
        # This would depend on how context_data is structured
        pass
    
    print(f"Created integrated dataset with {len(integrated_data)} rows and {len(integrated_data.columns)} columns")
    return integrated_data

# Function to generate summary statistics
def generate_summary_statistics(data):
    """Generate summary statistics for the integrated data."""
    if data is None:
        print("Cannot generate statistics: data is missing")
        return None
    
    summary_stats = {}
    
    # Variant count by treatment
    treatment_counts = data.groupby('Treatment').size().to_dict()
    summary_stats['treatment_counts'] = treatment_counts
    
    # Mutation type distribution
    mutation_counts = data.groupby('Std_Mutation').size().to_dict()
    summary_stats['mutation_counts'] = mutation_counts
    
    # Transition/transversion ratio by treatment
    ti_tv_ratio = {}
    for treatment in data['Treatment'].unique():
        treatment_data = data[data['Treatment'] == treatment]
        transitions = treatment_data[treatment_data['Mutation_Class'] == 'Transition'].shape[0]
        transversions = treatment_data[treatment_data['Mutation_Class'] == 'Transversion'].shape[0]
        ratio = transitions / transversions if transversions > 0 else float('inf')
        ti_tv_ratio[treatment] = ratio
    summary_stats['ti_tv_ratio'] = ti_tv_ratio
    
    # Scaffold statistics
    scaffold_stats = data.groupby('CHROM').agg({
        'POS': 'count',
        'Scaffold_Length': 'first',
        'Treatment': lambda x: list(x.unique())
    }).reset_index()
    
    if 'Scaffold_Length' in scaffold_stats.columns:
        scaffold_stats['Variant_Density'] = scaffold_stats.apply(
            lambda row: (row['POS'] * 1000 / row['Scaffold_Length']) 
            if row['Scaffold_Length'] > 0 else 0, axis=1)
    
    summary_stats['scaffold_stats'] = scaffold_stats
    
    # Create a DataFrame for easy reporting
    summary_df = pd.DataFrame({
        'Metric': ['Total Variants', 'Unique Scaffolds', 'Average Ti/Tv Ratio'],
        'Value': [len(data), len(scaffold_stats), np.mean(list(ti_tv_ratio.values()))]
    })
    
    return summary_stats, summary_df

# Function to prepare data for PCA
def prepare_data_for_pca(data):
    """Prepare data for Principal Component Analysis."""
    if data is None:
        print("Cannot prepare data for PCA: data is missing")
        return None, None
    
    # Create feature matrix for PCA
    # We'll use a pivot table to create a matrix of treatments by mutation types
    pivot_data = pd.pivot_table(
        data, 
        values='POS',  # Just a placeholder, we only care about counts
        index='Treatment',
        columns='Std_Mutation',
        aggfunc='count',
        fill_value=0
    )
    
    # Normalize by total mutations per treatment
    row_sums = pivot_data.sum(axis=1)
    normalized_pivot = pivot_data.div(row_sums, axis=0)
    
    # Prepare for PCA
    X = normalized_pivot.values
    feature_names = normalized_pivot.columns
    
    # Standard scale the data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    return X_scaled, feature_names, normalized_pivot

# Function to perform PCA
def perform_pca(X, feature_names, sample_names):
    """Perform Principal Component Analysis."""
    if X is None or len(X) == 0:
        print("Cannot perform PCA: data is missing or empty")
        return None
    
    # Determine number of components (2 or number of samples - 1, whichever is smaller)
    n_components = min(2, X.shape[0] - 1)
    if n_components < 1:
        print("Not enough samples for PCA")
        return None
    
    # Perform PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(X)
    
    # Create DataFrame with PCA results
    pca_df = pd.DataFrame(
        data=pca_result,
        columns=[f'PC{i+1}' for i in range(n_components)]
    )
    pca_df['Treatment'] = sample_names
    
    # Calculate feature contributions
    feature_contributions = pd.DataFrame(
        data=pca.components_.T,
        columns=[f'PC{i+1}' for i in range(n_components)],
        index=feature_names
    )
    
    # Calculate explained variance
    explained_variance = pca.explained_variance_ratio_
    
    return pca_df, feature_contributions, explained_variance

# Function to perform correlation analysis
def perform_correlation_analysis(data):
    """Perform correlation analysis between variants and features."""
    if data is None:
        print("Cannot perform correlation analysis: data is missing")
        return None
    
    # Create treatment-by-mutation type matrix
    mutation_pivot = pd.pivot_table(
        data, 
        values='POS',  # Just a placeholder, we only care about counts
        index='Treatment',
        columns='Std_Mutation',
        aggfunc='count',
        fill_value=0
    )
    
    # Calculate correlation between mutation types
    mutation_corr = mutation_pivot.corr(method='spearman')
    
    # Create scaffold-level dataset
    scaffold_data = data.groupby(['CHROM', 'Treatment']).size().unstack(fill_value=0)
    
    # Calculate correlation between treatments at scaffold level
    treatment_corr = scaffold_data.corr(method='spearman')
    
    return mutation_corr, treatment_corr

# Function to build regression models
def build_regression_models(data):
    """Build regression models to identify predictors of variant frequency."""
    if data is None or 'Variant_Density' not in data.columns:
        print("Cannot build regression models: data missing or incomplete")
        return None
    
    # Prepare scaffold-level dataset
    scaffold_data = data.groupby(['CHROM', 'Treatment']).agg({
        'Variant_Density': 'first',
        'Scaffold_Length': 'first'
    }).reset_index()
    
    # Create some derived features
    scaffold_data['Log_Length'] = np.log10(scaffold_data['Scaffold_Length'].replace(0, 1))
    
    # Build models for each treatment
    model_results = {}
    
    for treatment in scaffold_data['Treatment'].unique():
        treatment_data = scaffold_data[scaffold_data['Treatment'] == treatment]
        
        # Skip if too few data points
        if len(treatment_data) < 5:
            continue
        
        # Simple linear regression
        X = sm.add_constant(treatment_data[['Log_Length']])
        y = treatment_data['Variant_Density']
        
        try:
            model = sm.OLS(y, X).fit()
            
            # Extract key statistics
            model_results[treatment] = {
                'r_squared': model.rsquared,
                'coefficients': model.params.to_dict(),
                'p_values': model.pvalues.to_dict(),
                'n': len(treatment_data)
            }
        except:
            print(f"Could not build model for {treatment}")
    
    return model_results

# Function to perform clustering analysis
def perform_clustering_analysis(data):
    """Perform clustering analysis to identify patterns."""
    if data is None:
        print("Cannot perform clustering: data is missing")
        return None
    
    try:
        # Create a matrix of scaffold-by-treatment variant counts
        print("Creating scaffold-by-treatment matrix...")
        pivot_data = pd.pivot_table(
            data, 
            values='POS', 
            index='CHROM',
            columns='Treatment',
            aggfunc='count',
            fill_value=0
        )
        
        print(f"Created pivot table with {len(pivot_data)} scaffolds and {len(pivot_data.columns)} treatments")
        
        # Only keep scaffolds with at least 3 variants
        active_scaffolds = pivot_data.sum(axis=1) >= 3
        pivot_filtered = pivot_data[active_scaffolds]
        
        print(f"After filtering, {len(pivot_filtered)} scaffolds remain")
        
        if len(pivot_filtered) == 0:
            print("No scaffolds with sufficient variants for clustering")
            return None
        
        # Check for at least 2 scaffolds and 2 treatments (minimum for clustering)
        if len(pivot_filtered) < 2 or len(pivot_filtered.columns) < 2:
            print(f"Insufficient data for clustering: {len(pivot_filtered)} scaffolds, {len(pivot_filtered.columns)} treatments")
            return None
        
        # Normalize the data, carefully handling zero sums
        row_sums = pivot_filtered.sum(axis=1)
        
        # Check for zero row sums
        zero_rows = (row_sums == 0).sum()
        if zero_rows > 0:
            print(f"Warning: {zero_rows} rows have zero sum. Removing these rows.")
            pivot_filtered = pivot_filtered[row_sums > 0]
            row_sums = row_sums[row_sums > 0]
            
            if len(pivot_filtered) < 2:
                print("Insufficient non-zero rows for clustering")
                return None
        
        # Now safely normalize
        normalized_pivot = pivot_filtered.div(row_sums, axis=0)
        
        # Check for constant rows/columns which break correlation distance
        constant_cols = []
        for col in normalized_pivot.columns:
            if normalized_pivot[col].nunique() <= 1:
                constant_cols.append(col)
        
        constant_rows = []
        for idx in normalized_pivot.index:
            if normalized_pivot.loc[idx].nunique() <= 1:
                constant_rows.append(idx)
        
        if constant_cols:
            print(f"Warning: {len(constant_cols)} columns have constant values, removing them.")
            if len(normalized_pivot.columns) - len(constant_cols) < 2:
                print("Insufficient varying columns for clustering")
                return None
            normalized_pivot = normalized_pivot.drop(columns=constant_cols)
        
        if constant_rows:
            print(f"Warning: {len(constant_rows)} rows have constant values, removing them.")
            if len(normalized_pivot) - len(constant_rows) < 2:
                print("Insufficient varying rows for clustering")
                return None
            normalized_pivot = normalized_pivot.drop(index=constant_rows)
        
        print(f"Final data shape for clustering: {normalized_pivot.shape}")
        
        # Perform hierarchical clustering with error handling
        try:
            # Cluster the treatments (columns)
            treatment_linkage = hierarchy.linkage(
                normalized_pivot.T, 
                method='average',
                metric='correlation'
            )
            
            # Cluster the scaffolds (rows)
            scaffold_linkage = hierarchy.linkage(
                normalized_pivot, 
                method='average',
                metric='correlation'
            )
            
            return normalized_pivot, treatment_linkage, scaffold_linkage
            
        except Exception as e:
            print(f"Error with correlation distance: {str(e)}")
            print("Trying Euclidean distance instead...")
            
            try:
                # Try with Euclidean distance as fallback
                treatment_linkage = hierarchy.linkage(
                    normalized_pivot.T, 
                    method='average',
                    metric='euclidean'
                )
                
                scaffold_linkage = hierarchy.linkage(
                    normalized_pivot, 
                    method='average',
                    metric='euclidean'
                )
                
                print("Successfully clustered using Euclidean distance")
                return normalized_pivot, treatment_linkage, scaffold_linkage
                
            except Exception as e:
                print(f"Error with Euclidean distance clustering: {str(e)}")
                return None
    
    except Exception as e:
        print(f"Error performing clustering: {str(e)}")
        import traceback
        traceback.print_exc()  # Print the full stack trace for diagnosis
        return None

# Function to plot PCA results
def plot_pca_results(pca_df, feature_contributions, explained_variance, output_dir):
    """Plot PCA results."""
    if pca_df is None or 'PC1' not in pca_df.columns:
        print("Cannot plot PCA results: PCA failed or insufficient components")
        return
    
    # Plot treatments in PCA space
    plt.figure(figsize=(10, 8))
    
    # Define colors for treatments
    treatment_colors = {
        'WT': '#1b9e77', 
        'STC': '#d95f02', 
        'CAS': '#7570b3', 
        'WTA': '#e7298a'
    }
    
    # Plot each treatment
    for treatment in pca_df['Treatment'].unique():
        treatment_data = pca_df[pca_df['Treatment'] == treatment]
        plt.scatter(
            treatment_data['PC1'], 
            treatment_data['PC2'] if 'PC2' in treatment_data.columns else [0] * len(treatment_data),
            label=treatment,
            color=treatment_colors.get(treatment, 'gray'),
            s=100
        )
    
    # Add treatment labels
    for i, row in pca_df.iterrows():
        plt.annotate(
            row['Treatment'],
            (row['PC1'], row['PC2'] if 'PC2' in row else 0),
            textcoords="offset points",
            xytext=(0, 10),
            ha='center'
        )
    
    # Add axis labels with explained variance
    if len(explained_variance) >= 1:
        plt.xlabel(f'PC1 ({explained_variance[0]:.2%} explained variance)')
        
        if len(explained_variance) >= 2:
            plt.ylabel(f'PC2 ({explained_variance[1]:.2%} explained variance)')
        else:
            plt.ylabel('PC2')
    
    plt.title('PCA of Mutation Patterns by Treatment')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "pca_treatments.png"), dpi=300)
    plt.close()
    
    # Plot feature contributions to PCs
    if feature_contributions is not None and not feature_contributions.empty:
        plt.figure(figsize=(12, 8))
        
        # Sort features by contribution to PC1
        sorted_features = feature_contributions.sort_values('PC1', ascending=False)
        
        # Plot contributions to PC1
        plt.subplot(2, 1, 1)
        plt.bar(range(len(sorted_features)), sorted_features['PC1'], color='skyblue')
        plt.xticks(range(len(sorted_features)), sorted_features.index, rotation=90)
        plt.ylabel('Contribution to PC1')
        plt.title('Feature Contributions to PC1')
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # If we have PC2, plot contributions to PC2
        if 'PC2' in sorted_features.columns:
            # Sort features by contribution to PC2
            sorted_features_pc2 = feature_contributions.sort_values('PC2', ascending=False)
            
            plt.subplot(2, 1, 2)
            plt.bar(range(len(sorted_features_pc2)), sorted_features_pc2['PC2'], color='lightgreen')
            plt.xticks(range(len(sorted_features_pc2)), sorted_features_pc2.index, rotation=90)
            plt.ylabel('Contribution to PC2')
            plt.title('Feature Contributions to PC2')
            plt.grid(True, linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "pca_feature_contributions.png"), dpi=300)
        plt.close()

# Function to plot correlation results
def plot_correlation_results(mutation_corr, treatment_corr, output_dir):
    """Plot correlation analysis results."""
    # Plot mutation type correlation matrix
    if mutation_corr is not None and not mutation_corr.empty:
        # Filter to only include standard 6 mutation types
        standard_types = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        
        # Check if these standard types exist in our data
        available_types = [mt for mt in standard_types if mt in mutation_corr.columns]
        
        if available_types:
            # Filter matrix to only include standard SNV types
            filtered_corr = mutation_corr.loc[available_types, available_types]
            
            plt.figure(figsize=(10, 8))
            sns.heatmap(
                filtered_corr, 
                annot=True, 
                cmap='coolwarm', 
                center=0, 
                vmin=-1, 
                vmax=1,
                fmt='.2f',  # Format to 2 decimal places for clarity
                annot_kws={"size": 12}  # Larger annotation text
            )
            plt.title('Correlation Between Mutation Types')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "mutation_correlation.png"), dpi=300)
            plt.close()
        else:
            # If standard types not found, try to create a cleaner version of existing matrix
            # First determine if the matrix is too large
            if len(mutation_corr) > 15:  # Arbitrary threshold for "too large"
                # Extract only SNVs based on pattern (X>Y where X and Y are single nucleotides)
                snv_pattern = r'^[ACGT]>[ACGT]$'
                snv_types = [col for col in mutation_corr.columns if re.match(snv_pattern, col)]
                
                if snv_types:
                    filtered_corr = mutation_corr.loc[snv_types, snv_types]
                else:
                    # Take the top most frequent variants
                    top_types = list(mutation_corr.columns[:10])  # Take top 10
                    filtered_corr = mutation_corr.loc[top_types, top_types]
            else:
                filtered_corr = mutation_corr
            
            plt.figure(figsize=(12, 10))
            sns.heatmap(
                filtered_corr, 
                annot=True, 
                cmap='coolwarm', 
                center=0, 
                vmin=-1, 
                vmax=1,
                fmt='.2f',
                annot_kws={"size": 10}
            )
            plt.title('Correlation Between Mutation Types')
            plt.xticks(rotation=45, ha='right')
            plt.yticks(rotation=0)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "mutation_correlation.png"), dpi=300)
            plt.close()
    
    # Plot treatment correlation matrix - this part is fine as is
    if treatment_corr is not None and not treatment_corr.empty:
        plt.figure(figsize=(8, 6))
        sns.heatmap(treatment_corr, annot=True, cmap='coolwarm', center=0, vmin=-1, vmax=1)
        plt.title('Correlation Between Treatments')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "treatment_correlation.png"), dpi=300)
        plt.close()

# Function to plot regression results
def plot_regression_results(data, model_results, output_dir):
    """Plot regression analysis results."""
    if data is None or model_results is None:
        return
    
    # Create scaffold-level dataset for plotting
    scaffold_data = data.groupby(['CHROM', 'Treatment']).agg({
        'Variant_Density': 'first',
        'Scaffold_Length': 'first'
    }).reset_index()
    
    # Add log-transformed length
    scaffold_data['Log_Length'] = np.log10(scaffold_data['Scaffold_Length'].replace(0, 1))
    
    # Create regression plots for each treatment
    plt.figure(figsize=(12, 10))
    
    # Define colors for treatments
    treatment_colors = {
        'WT': '#1b9e77', 
        'STC': '#d95f02', 
        'CAS': '#7570b3', 
        'WTA': '#e7298a'
    }
    
    # Set up subplots - 2x2 grid
    treatments = list(model_results.keys())
    n_treatments = len(treatments)
    
    for i, treatment in enumerate(treatments):
        if i < 4:  # Only show up to 4 treatments
            plt.subplot(2, 2, i+1)
            
            # Get treatment data
            treatment_data = scaffold_data[scaffold_data['Treatment'] == treatment]
            
            # Scatter plot
            plt.scatter(
                treatment_data['Log_Length'],
                treatment_data['Variant_Density'],
                color=treatment_colors.get(treatment, 'gray'),
                alpha=0.6,
                s=50,
                label=treatment
            )
            
            # Add regression line
            model_info = model_results[treatment]
            if 'coefficients' in model_info:
                # Get coefficients
                intercept = model_info['coefficients'].get('const', 0)
                slope = model_info['coefficients'].get('Log_Length', 0)
                
                # Generate points for line
                x_range = np.linspace(
                    treatment_data['Log_Length'].min(),
                    treatment_data['Log_Length'].max(),
                    100
                )
                y_range = intercept + slope * x_range
                
                # Plot regression line
                plt.plot(x_range, y_range, color='black', linestyle='--')
                
                # Add R-squared value
                r2 = model_info['r_squared']
                plt.text(
                    0.05, 0.95, 
                    f"R² = {r2:.3f}",
                    transform=plt.gca().transAxes,
                    fontsize=12,
                    verticalalignment='top'
                )
            
            plt.xlabel('Log10(Scaffold Length)')
            plt.ylabel('Variant Density (per kb)')
            plt.title(f'{treatment} Treatment')
            plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "regression_analysis.png"), dpi=300)
    plt.close()
    
    # Create a summary bar plot of R-squared values
    plt.figure(figsize=(10, 6))
    
    treatments = []
    r2_values = []
    
    for treatment, model_info in model_results.items():
        treatments.append(treatment)
        r2_values.append(model_info['r_squared'])
    
    # Sort by R-squared value
    sorted_indices = np.argsort(r2_values)[::-1]  # Descending order
    sorted_treatments = [treatments[i] for i in sorted_indices]
    sorted_r2 = [r2_values[i] for i in sorted_indices]
    
    plt.bar(
        range(len(sorted_treatments)), 
        sorted_r2,
        color=[treatment_colors.get(t, 'gray') for t in sorted_treatments]
    )
    
    plt.xticks(range(len(sorted_treatments)), sorted_treatments)
    plt.ylabel('R-squared')
    plt.title('Model Fit (R²) by Treatment')
    plt.grid(True, linestyle='--', alpha=0.7, axis='y')
    
    # Add value labels
    for i, r2 in enumerate(sorted_r2):
        plt.text(
            i, r2 + 0.01, 
            f"{r2:.3f}", 
            ha='center'
        )
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "r2_by_treatment.png"), dpi=300)
    plt.close()

# Function to plot clustering results
def plot_clustering_results(pivot_data, treatment_linkage, scaffold_linkage, output_dir):
    """Plot clustering analysis results."""
    if pivot_data is None or pivot_data.empty:
        return
    
    # Create clustered heatmap
    plt.figure(figsize=(12, 10))
    
    # Create clustered heatmap with both row and column dendrograms
    sns.clustermap(
        pivot_data,
        row_linkage=scaffold_linkage,
        col_linkage=treatment_linkage,
        cmap='YlGnBu',
        figsize=(12, 10),
        yticklabels=False  # Too many scaffolds to show labels
    )
    
    plt.suptitle('Hierarchical Clustering of Scaffold-Treatment Patterns', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "scaffold_treatment_clustering.png"), dpi=300)
    plt.close()
    
    # Create treatment dendrogram
    plt.figure(figsize=(10, 6))
    
    hierarchy.dendrogram(
        treatment_linkage,
        labels=pivot_data.columns,
        leaf_font_size=12
    )
    
    plt.title('Hierarchical Clustering of Treatments')
    plt.xlabel('Treatment')
    plt.ylabel('Distance')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "treatment_dendrogram.png"), dpi=300)
    plt.close()

# Function to create summary report
def create_summary_report(summary_df, model_results, output_dir):
    """Create summary report of statistical analysis results."""
    with open(os.path.join(output_dir, "statistical_analysis_summary.txt"), 'w') as f:
        f.write("Statistical Pattern Analysis Summary\n")
        f.write("==================================\n\n")
        
        # Write basic statistics
        f.write("Basic Statistics:\n")
        f.write("-----------------\n")
        
        for i, row in summary_df.iterrows():
            f.write(f"{row['Metric']}: {row['Value']}\n")
        
        f.write("\n")
        
        # Write regression model results
        f.write("Regression Analysis Results:\n")
        f.write("--------------------------\n")
        
        for treatment, model_info in model_results.items():
            f.write(f"\n{treatment} Treatment:\n")
            f.write(f"  R-squared: {model_info['r_squared']:.4f}\n")
            f.write("  Coefficients:\n")
            
            for param, value in model_info['coefficients'].items():
                p_value = model_info['p_values'].get(param, 1.0)
                significance = ""
                if p_value < 0.001:
                    significance = "***"
                elif p_value < 0.01:
                    significance = "**"
                elif p_value < 0.05:
                    significance = "*"
                
                f.write(f"    {param}: {value:.4f} {significance}\n")
            
            f.write(f"  Sample size: {model_info['n']}\n")
        
        f.write("\n")
        
        # Write main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis examines statistical patterns in the variant data.\n")
        f.write("2. Principal Component Analysis identifies major axes of variation between treatments.\n")
        f.write("3. Correlation analysis reveals relationships between mutation types and treatments.\n")
        f.write("4. Regression models evaluate the relationship between genomic features and variant density.\n")
        f.write("5. Clustering analysis identifies groups of scaffolds with similar variant patterns.\n")

# Main function to run the analysis
def main():
    # Load data from previous analyses
    mutation_data = load_mutation_data()
    scaffold_info = load_scaffold_info()
    context_data = load_context_data()
    
    # Integrate data sources
    integrated_data = integrate_data(mutation_data, scaffold_info, context_data)
    
    if integrated_data is None:
        print("Cannot continue analysis: integrated data is missing")
        return
    
    # Generate summary statistics
    summary_stats, summary_df = generate_summary_statistics(integrated_data)
    
    # Prepare data for PCA
    X_scaled, feature_names, normalized_pivot = prepare_data_for_pca(integrated_data)
    
    # Perform PCA
    if X_scaled is not None and normalized_pivot is not None:
        pca_results = perform_pca(X_scaled, feature_names, normalized_pivot.index)
        if pca_results:
            pca_df, feature_contributions, explained_variance = pca_results
            # Plot PCA results
            plot_pca_results(pca_df, feature_contributions, explained_variance, OUTPUT_DIR)
    
    # Perform correlation analysis
    if integrated_data is not None:
        mutation_corr, treatment_corr = perform_correlation_analysis(integrated_data)
        # Plot correlation results
        plot_correlation_results(mutation_corr, treatment_corr, OUTPUT_DIR)
    
    # Build regression models
    model_results = build_regression_models(integrated_data)
    if model_results:
        # Plot regression results
        plot_regression_results(integrated_data, model_results, OUTPUT_DIR)
    
    # Perform clustering analysis
    clustering_results = perform_clustering_analysis(integrated_data)
    if clustering_results:
        pivot_data, treatment_linkage, scaffold_linkage = clustering_results
        # Plot clustering results
        plot_clustering_results(pivot_data, treatment_linkage, scaffold_linkage, OUTPUT_DIR)
    
    # Create summary report
    create_summary_report(summary_df, model_results, OUTPUT_DIR)
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")

# Run the analysis
if __name__ == "__main__":
    main()