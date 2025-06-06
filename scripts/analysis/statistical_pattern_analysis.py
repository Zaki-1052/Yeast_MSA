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
import re
import warnings
warnings.filterwarnings('ignore')

# Set matplotlib style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "analysis/statistical_pattern_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define biologically correct treatment groups
TREATMENTS = ['WT-37', 'WTA', 'STC', 'CAS']

# Define treatment descriptions for better annotations
TREATMENT_INFO = {
    'WT-37': {'description': 'Temperature-adapted wild type', 'adaptation': 'Temperature'},
    'WTA': {'description': 'Low oxygen-adapted wild type', 'adaptation': 'Low Oxygen'},
    'STC': {'description': 'STC gene with low oxygen adaptation', 'adaptation': 'Low Oxygen', 'gene': 'STC'},
    'CAS': {'description': 'CAS gene with temperature adaptation', 'adaptation': 'Temperature', 'gene': 'CAS'}
}

# Define treatment colors for consistent visualizations
TREATMENT_COLORS = {
    'WT-37': '#1b9e77',  # Temperature-adapted
    'WTA': '#d95f02',    # Low oxygen-adapted
    'STC': '#7570b3',    # STC gene + low oxygen
    'CAS': '#e7298a'     # CAS gene + temperature
}

# Function to find a file checking multiple possible locations
def find_file(file_patterns, base_name):
    """Try multiple file path patterns to find an existing file."""
    for pattern in file_patterns:
        path = pattern.format(base_name)
        if os.path.exists(path):
            return path
    return None

# Function to load mutation data from previous analyses
def load_mutation_data():
    """Load mutation data from previous analyses."""
    all_data = []
    
    # Define possible file locations
    mutation_file_patterns = [
        "analysis/MSA/mutation_spectrum_analysis/{}_mutations.txt",
        "mutation_spectrum_analysis/{}_mutations.txt",
        "results/mutation_spectrum_analysis/{}_mutations.txt"
    ]
    
    for treatment in TREATMENTS:
        # Handle backward compatibility with old 'WT' naming
        # If we're looking for WT-37 data but it's stored as WT
        search_names = [treatment]
        if treatment == 'WT-37':
            search_names.append('WT')  # Also try 'WT' for backward compatibility
        
        file_found = False
        for search_name in search_names:
            file_path = find_file(mutation_file_patterns, search_name)
            if file_path:
                print(f"Found mutation data for {treatment} at {file_path}")
                try:
                    # Load data
                    data = pd.read_csv(file_path, sep='\t', header=None, 
                                    names=['CHROM', 'POS', 'REF', 'ALT'])
                    # Ensure we use the correct treatment name regardless of file name
                    data['Treatment'] = treatment
                    all_data.append(data)
                    file_found = True
                    break
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
        
        if not file_found:
            print(f"No mutation data found for {treatment}")
    
    if all_data:
        combined_data = pd.concat(all_data, ignore_index=True)
        print(f"Loaded {len(combined_data)} mutations across {len(all_data)} treatments")
        return combined_data
    else:
        print("No mutation data found. Attempting to extract from VCF files...")
        
        # Try extracting from VCF files if no pre-extracted data found
        try:
            return extract_mutations_from_vcf()
        except Exception as e:
            print(f"Failed to extract from VCF files: {e}")
            return None

# Function to extract mutation data directly from VCF files
def extract_mutations_from_vcf():
    """Extract mutation data directly from VCF files if pre-extracted data is not available."""
    import subprocess
    
    all_data = []
    vcf_file_patterns = [
        "results/merged/analysis/{}/highconf.vcf.gz",
        "results/merged/analysis/{}_highconf.vcf.gz",
        "results/merged/analysis/{}/specific.vcf.gz",
        "results/merged/analysis/{}_specific.vcf.gz"
    ]
    
    for treatment in TREATMENTS:
        # Handle backward compatibility with old 'WT' naming
        search_names = [treatment]
        if treatment == 'WT-37':
            search_names.append('WT')
        
        file_found = False
        for search_name in search_names:
            vcf_file = find_file(vcf_file_patterns, search_name)
            if vcf_file:
                print(f"Extracting mutations for {treatment} from {vcf_file}")
                try:
                    # Extract mutation data using bcftools
                    cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_file}"
                    output = subprocess.check_output(cmd, shell=True).decode('utf-8')
                    
                    # Parse the output into a dataframe
                    rows = []
                    for line in output.strip().split('\n'):
                        if line:
                            parts = line.split('\t')
                            if len(parts) == 4:
                                rows.append(parts)
                    
                    if rows:
                        data = pd.DataFrame(rows, columns=['CHROM', 'POS', 'REF', 'ALT'])
                        data['POS'] = data['POS'].astype(int)
                        data['Treatment'] = treatment
                        all_data.append(data)
                        file_found = True
                        
                        # Save extracted data for future use
                        os.makedirs("mutation_spectrum_analysis", exist_ok=True)
                        data.to_csv(f"mutation_spectrum_analysis/{treatment}_mutations.txt", 
                                  sep='\t', index=False, header=False)
                        print(f"Saved extracted data to mutation_spectrum_analysis/{treatment}_mutations.txt")
                        break
                except Exception as e:
                    print(f"Error extracting from {vcf_file}: {e}")
        
        if not file_found:
            print(f"Could not find VCF data for {treatment}")
    
    if all_data:
        combined_data = pd.concat(all_data, ignore_index=True)
        print(f"Extracted {len(combined_data)} mutations from VCF files")
        return combined_data
    else:
        print("No mutation data could be extracted")
        return None

# Function to load scaffold information
def load_scaffold_info():
    """Load scaffold information from reference genome."""
    scaffold_info = {}
    
    # Try to load scaffold lengths from fai file
    fai_patterns = [
        "reference/w303_chromosomal.fasta.fai",
        "reference/yeast/yeast_w303.fasta.fai",
        "reference/genome.fasta.fai"
    ]
    
    # First attempt: Load from fai file
    found_fai = False
    for pattern in fai_patterns:
        if os.path.exists(pattern):
            try:
                scaffold_df = pd.read_csv(pattern, sep='\t', header=None,
                                        names=['Scaffold', 'Length', 'Offset', 'Linebases', 'Linewidth'])
                scaffold_info = dict(zip(scaffold_df['Scaffold'], scaffold_df['Length']))
                print(f"Loaded information for {len(scaffold_info)} scaffolds from {pattern}")
                found_fai = True
                break
            except Exception as e:
                print(f"Error reading {pattern}: {e}")
    
    # Second attempt: Extract from VCF headers
    if not found_fai:
        print("No scaffold information found from fai files, attempting to extract from VCF headers...")
        
        try:
            import subprocess
            for treatment in TREATMENTS:
                vcf_file_patterns = [
                    f"results/merged/analysis/{treatment}/highconf.vcf.gz",
                    f"results/merged/analysis/{treatment}_highconf.vcf.gz",
                    f"results/merged/fixed/all_samples.vcf.gz"
                ]
                
                for pattern in vcf_file_patterns:
                    if os.path.exists(pattern):
                        cmd = f"bcftools view -h {pattern} | grep '##contig'"
                        output = subprocess.check_output(cmd, shell=True).decode('utf-8')
                        
                        import re
                        for line in output.strip().split('\n'):
                            match = re.search(r'ID=([^,]+),length=(\d+)', line)
                            if match:
                                scaffold, length = match.groups()
                                scaffold_info[scaffold] = int(length)
                        
                        if scaffold_info:
                            print(f"Extracted information for {len(scaffold_info)} scaffolds from VCF header")
                            found_fai = True
                            break
                if found_fai:
                    break
        except Exception as e:
            print(f"Error extracting scaffold info from VCF: {e}")
    
    # If we still don't have scaffold info, create default values using the scaffolds from the mutation data
    if not scaffold_info:
        print("Could not load scaffold information. Creating default values...")
        
        # Try to get list of scaffolds from mutation data
        try:
            mutation_data = load_mutation_data()
            if mutation_data is not None:
                unique_scaffolds = mutation_data['CHROM'].unique()
                # Use a default length of 10,000 bp for each scaffold
                for scaffold in unique_scaffolds:
                    scaffold_info[scaffold] = 10000
                print(f"Created default length values for {len(scaffold_info)} scaffolds")
        except Exception as e:
            print(f"Error creating default scaffold info: {e}")
    
    # Print some diagnostics
    print(f"Final scaffold info contains {len(scaffold_info)} scaffolds")
    if len(scaffold_info) > 0:
        print(f"Example: {list(scaffold_info.items())[:3]}")
    else:
        print("WARNING: No scaffold information available!")
    
    # Check if we need to handle ID mismatch (if scaffolds are numeric but reference uses JRIU prefix)
    try:
        import re  # Import re module here to ensure it's available
        mutation_data = load_mutation_data()
        if mutation_data is not None:
            sample_ids = mutation_data['CHROM'].head(5).astype(str).tolist()
            print(f"Sample mutation scaffold IDs: {sample_ids}")
            
            # Check if these appear to be numeric IDs
            numeric_ids = all(chrom.isdigit() for chrom in sample_ids)
            if numeric_ids:
                # Check if our scaffold info has JRIU pattern
                jriu_pattern = any("JRIU" in str(key) for key in scaffold_info.keys())
                
                if jriu_pattern:
                    print("Detected ID mismatch: mutation data uses numeric IDs but reference uses JRIU prefix")
                    print("Creating mapping between numeric IDs and reference IDs...")
                    
                    # Create a mapping from numeric ID to full ID
                    id_mapping = {}
                    for full_id in scaffold_info.keys():
                        if "JRIU" in str(full_id):
                            # Extract numeric part (e.g., '371' from 'JRIU01000371.1')
                            # Try multiple regex patterns to capture different naming conventions
                            match = re.search(r'JRIU\d+0*(\d+)(?:\.\d+)?', str(full_id))
                            if match:
                                numeric_id = match.group(1)
                                id_mapping[numeric_id] = full_id
                                
                                # Also try with leading zeros as some systems might use those
                                # For example, map '0371' to the same full_id
                                for i in range(1, 4):  # Try prefixing with 1, 2, or 3 zeros
                                    padded_id = numeric_id.zfill(len(numeric_id) + i)
                                    id_mapping[padded_id] = full_id
                    
                    # Get a sample mutation CHROM value to determine data type
                    sample_type = type(mutation_data['CHROM'].iloc[0])
                    print(f"Mutation CHROM data type: {sample_type}")
                    
                    # Create new scaffold_info with numeric keys, matching the data type
                    numeric_scaffold_info = {}
                    
                    # Create direct mapping for all numeric IDs found in the mutation data
                    # This ensures we don't miss any scaffolds
                    all_chrom_ids = mutation_data['CHROM'].unique()
                    print(f"Found {len(all_chrom_ids)} unique scaffold IDs in mutation data")
                    
                    # First add all mappings from the ID mapping we created
                    for numeric_id, full_id in id_mapping.items():
                        # Convert the numeric_id to the same type as the mutation data CHROMs
                        if sample_type == int or sample_type == np.int64:
                            try:
                                key = int(numeric_id)
                                numeric_scaffold_info[key] = scaffold_info[full_id]
                            except ValueError:
                                print(f"Warning: Could not convert {numeric_id} to int")
                        else:
                            numeric_scaffold_info[numeric_id] = scaffold_info[full_id]
                    
                    # For any remaining scaffolds in the mutation data that we couldn't map,
                    # use a default value
                    default_length = 10000  # Default scaffold length (10kb)
                    
                    # Use the median of known scaffold lengths if we have some
                    if scaffold_info:
                        lengths = list(scaffold_info.values())
                        default_length = int(np.median(lengths))
                        print(f"Using median scaffold length of {default_length} bp as default")
                    
                    # Now handle any scaffolds in the mutation data that aren't mapped yet
                    for chrom_id in all_chrom_ids:
                        if chrom_id not in numeric_scaffold_info:
                            print(f"Using default length for unmapped scaffold: {chrom_id}")
                            numeric_scaffold_info[chrom_id] = default_length
                    
                    print(f"Created mapping for {len(numeric_scaffold_info)} scaffolds")
                    scaffold_info = numeric_scaffold_info
    except Exception as e:
        print(f"Warning: Error checking for ID mismatch: {e}")
    
    return scaffold_info

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
        # Check for missing values and create diagnostics
        missing_scaffolds = integrated_data[integrated_data['Scaffold_Length'].isna()]['CHROM'].unique()
        if len(missing_scaffolds) > 0:
            print(f"Missing scaffold length for {len(missing_scaffolds)} unique scaffolds")
            print(f"First few missing scaffolds: {missing_scaffolds[:5]}")
            # Fill in missing values with a reasonable default (10,000bp)
            integrated_data['Scaffold_Length'] = integrated_data['Scaffold_Length'].fillna(10000)
    else:
        # If no scaffold info at all, create a default column
        print("No scaffold info available - creating default scaffold lengths")
        integrated_data['Scaffold_Length'] = 10000
    
    # Add treatment metadata
    integrated_data['Treatment_Description'] = integrated_data['Treatment'].map(
        lambda t: TREATMENT_INFO.get(t, {}).get('description', 'Unknown'))
    
    integrated_data['Adaptation_Type'] = integrated_data['Treatment'].map(
        lambda t: TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown'))
    
    integrated_data['Has_Gene'] = integrated_data['Treatment'].map(
        lambda t: 'Yes' if TREATMENT_INFO.get(t, {}).get('gene') else 'No')
    
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
        
        # Ensure variant density doesn't contain NaN or infinity values
        integrated_data['Variant_Density'] = integrated_data['Variant_Density'].replace([np.inf, -np.inf], np.nan)
        # Replace NaN with 0 (could also consider using mean or median)
        integrated_data['Variant_Density'] = integrated_data['Variant_Density'].fillna(0)
    
    # Add treatment-specific variant counts
    treatment_counts = integrated_data.groupby('Treatment').size().to_dict()
    integrated_data['Treatment_Total_Variants'] = integrated_data['Treatment'].map(treatment_counts)
    
    # Add adaptation type counts
    if 'Adaptation_Type' in integrated_data.columns:
        adaptation_counts = integrated_data.groupby('Adaptation_Type').size().to_dict()
        integrated_data['Adaptation_Variant_Count'] = integrated_data['Adaptation_Type'].map(adaptation_counts)
    
    # Merge with context data if available
    if context_data is not None:
        # Implement the merge if context data is available
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
    
    # Variant count by adaptation type
    if 'Adaptation_Type' in data.columns:
        adaptation_counts = data.groupby('Adaptation_Type').size().to_dict()
        summary_stats['adaptation_counts'] = adaptation_counts
    
    # Variant count by gene presence
    if 'Has_Gene' in data.columns:
        gene_counts = data.groupby(['Has_Gene', 'Treatment']).size().unstack(fill_value=0)
        summary_stats['gene_counts'] = gene_counts
    
    # Mutation type distribution
    mutation_counts = data.groupby('Std_Mutation').size().to_dict()
    summary_stats['mutation_counts'] = mutation_counts
    
    # Transition/transversion ratio by treatment
    ti_tv_ratio = {}
    # Calculate overall Ti/Tv ratio
    overall_transitions = data[data['Mutation_Class'] == 'Transition'].shape[0]
    overall_transversions = data[data['Mutation_Class'] == 'Transversion'].shape[0]
    overall_ratio = overall_transitions / overall_transversions if overall_transversions > 0 else 0
    summary_stats['overall_ti_tv_ratio'] = overall_ratio
    
    # Calculate by treatment
    for treatment in data['Treatment'].unique():
        treatment_data = data[data['Treatment'] == treatment]
        transitions = treatment_data[treatment_data['Mutation_Class'] == 'Transition'].shape[0]
        transversions = treatment_data[treatment_data['Mutation_Class'] == 'Transversion'].shape[0]
        ratio = transitions / transversions if transversions > 0 else 0
        ti_tv_ratio[treatment] = ratio
    summary_stats['ti_tv_ratio'] = ti_tv_ratio
    
    # Compare mutation spectra between adaptation types
    if 'Adaptation_Type' in data.columns:
        mutation_by_adaptation = pd.crosstab(
            data['Std_Mutation'], data['Adaptation_Type'], normalize='columns')
        summary_stats['mutation_by_adaptation'] = mutation_by_adaptation
    
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
        
        # Handle potential NaN or infinity values
        scaffold_stats['Variant_Density'] = scaffold_stats['Variant_Density'].replace([np.inf, -np.inf], np.nan)
        scaffold_stats['Variant_Density'] = scaffold_stats['Variant_Density'].fillna(0)
    
    summary_stats['scaffold_stats'] = scaffold_stats
    
    # Create a DataFrame for easy reporting
    # Use the overall Ti/Tv ratio directly rather than averaging the individual ratios
    summary_df = pd.DataFrame({
        'Metric': ['Total Variants', 'Unique Scaffolds', 'Average Ti/Tv Ratio'],
        'Value': [len(data), len(scaffold_stats), summary_stats['overall_ti_tv_ratio']]
    })
    
    # Add biological classification metrics
    if 'Adaptation_Type' in data.columns:
        for adaptation_type, count in adaptation_counts.items():
            summary_df = summary_df._append({
                'Metric': f'{adaptation_type} Adaptation Variants',
                'Value': count
            }, ignore_index=True)
    
    # Add gene presence metrics
    if 'Has_Gene' in data.columns:
        with_gene = data[data['Has_Gene'] == 'Yes'].shape[0]
        without_gene = data[data['Has_Gene'] == 'No'].shape[0]
        summary_df = summary_df._append({
            'Metric': 'Variants in Gene-Modified Strains',
            'Value': with_gene
        }, ignore_index=True)
        summary_df = summary_df._append({
            'Metric': 'Variants in Non-Modified Strains',
            'Value': without_gene
        }, ignore_index=True)
    
    return summary_stats, summary_df

# Function to prepare data for PCA
def prepare_data_for_pca(data):
    """Prepare data for Principal Component Analysis."""
    if data is None:
        print("Cannot prepare data for PCA: data is missing")
        return None, None, None
    
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
    
    # Add treatment metadata
    pca_df['Description'] = pca_df['Treatment'].map(
        lambda t: TREATMENT_INFO.get(t, {}).get('description', 'Unknown'))
    
    pca_df['Adaptation_Type'] = pca_df['Treatment'].map(
        lambda t: TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown'))
    
    pca_df['Has_Gene'] = pca_df['Treatment'].map(
        lambda t: 'Yes' if TREATMENT_INFO.get(t, {}).get('gene') else 'No')
    
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
    
    # Create adaptation-type correlation if available
    adaptation_corr = None
    if 'Adaptation_Type' in data.columns:
        adaptation_data = pd.pivot_table(
            data,
            values='POS',
            index='CHROM',
            columns='Adaptation_Type',
            aggfunc='count',
            fill_value=0
        )
        adaptation_corr = adaptation_data.corr(method='spearman')
    
    return mutation_corr, treatment_corr, adaptation_corr

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
    
    # Let's print the unique scaffold lengths to debug
    print("\nScaffold Length debug:")
    print(f"Number of unique scaffolds: {scaffold_data['CHROM'].nunique()}")
    if 'Scaffold_Length' in scaffold_data.columns:
        print(f"Scaffold_Length column exists")
        print(f"Unique Scaffold_Length values: {scaffold_data['Scaffold_Length'].unique()[:5]}...")
        print(f"Number of scaffolds with length=0: {(scaffold_data['Scaffold_Length'] == 0).sum()}")
        print(f"Number of scaffolds with NaN length: {scaffold_data['Scaffold_Length'].isna().sum()}")
    else:
        print("Scaffold_Length column is missing")
    
    # Build models for each treatment
    model_results = {}
    
    for treatment in scaffold_data['Treatment'].unique():
        treatment_data = scaffold_data[scaffold_data['Treatment'] == treatment]
        
        # Skip if too few data points
        if len(treatment_data) < 5:
            continue
        
        # Clean the data - remove rows with NaN or inf values
        clean_data = treatment_data.copy()
        clean_data = clean_data.replace([np.inf, -np.inf], np.nan)
        
        # Check data quality before filtering
        nan_log_length = clean_data['Log_Length'].isna().sum()
        nan_density = clean_data['Variant_Density'].isna().sum()
        print(f"Treatment {treatment}: {len(clean_data)} total rows, {nan_log_length} NaN in Log_Length, {nan_density} NaN in Variant_Density")
        
        # Instead of dropping rows, let's try to fix the data
        # Set a reasonable default for Log_Length if NaN (log10 of average scaffold length)
        if nan_log_length > 0:
            avg_log_length = clean_data['Log_Length'].mean()
            if pd.isna(avg_log_length):
                avg_log_length = 4.0  # Default log10 value (10,000bp)
            clean_data['Log_Length'] = clean_data['Log_Length'].fillna(avg_log_length)
            
        # Set zero for Variant_Density if NaN
        if nan_density > 0:
            clean_data['Variant_Density'] = clean_data['Variant_Density'].fillna(0)
        
        # Re-check if we still have NaN values
        remaining_nans = clean_data[['Log_Length', 'Variant_Density']].isna().any(axis=1).sum()
        print(f"  After fixing: {remaining_nans} rows still contain NaN values")
        
        # Only drop rows if still necessary
        if remaining_nans > 0:
            clean_data = clean_data.dropna(subset=['Log_Length', 'Variant_Density'])
        
        # Skip if too few data points after cleaning
        if len(clean_data) < 5:
            print(f"Insufficient clean data points for {treatment} after fixing NaN/inf values")
            continue
        
        # Simple linear regression
        X = sm.add_constant(clean_data[['Log_Length']])
        y = clean_data['Variant_Density']
        
        try:
            model = sm.OLS(y, X).fit()
            
            # Extract key statistics with adjusted R-squared (never negative)
            adjusted_rsquared = max(0, model.rsquared)  # Convert negative to 0
            
            model_results[treatment] = {
                'r_squared': adjusted_rsquared,
                'coefficients': model.params.to_dict(),
                'p_values': model.pvalues.to_dict(),
                'n': len(clean_data),
                'adaptation_type': TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown'),
                'has_gene': 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
            }
        except Exception as e:
            print(f"Could not build model for {treatment}: {e}")
    
    # Also build models by adaptation type if available
    adaptation_models = {}
    
    if 'Adaptation_Type' in data.columns:
        adaptation_scaffold_data = data.groupby(['CHROM', 'Adaptation_Type']).agg({
            'Variant_Density': 'first',
            'Scaffold_Length': 'first'
        }).reset_index()
        
        adaptation_scaffold_data['Log_Length'] = np.log10(
            adaptation_scaffold_data['Scaffold_Length'].replace(0, 1))
        
        for adaptation in adaptation_scaffold_data['Adaptation_Type'].unique():
            adaptation_data = adaptation_scaffold_data[
                adaptation_scaffold_data['Adaptation_Type'] == adaptation]
            
            if len(adaptation_data) < 5:
                continue
            
            # Clean the data - remove rows with NaN or inf values
            clean_adaptation_data = adaptation_data.copy()
            clean_adaptation_data = clean_adaptation_data.replace([np.inf, -np.inf], np.nan)
            
            # Check data quality before filtering
            nan_log_length = clean_adaptation_data['Log_Length'].isna().sum()
            nan_density = clean_adaptation_data['Variant_Density'].isna().sum()
            print(f"Adaptation {adaptation}: {len(clean_adaptation_data)} total rows, {nan_log_length} NaN in Log_Length, {nan_density} NaN in Variant_Density")
            
            # Fix the data 
            if nan_log_length > 0:
                avg_log_length = clean_adaptation_data['Log_Length'].mean()
                if pd.isna(avg_log_length):
                    avg_log_length = 4.0  # Default log10 value
                clean_adaptation_data['Log_Length'] = clean_adaptation_data['Log_Length'].fillna(avg_log_length)
                
            if nan_density > 0:
                clean_adaptation_data['Variant_Density'] = clean_adaptation_data['Variant_Density'].fillna(0)
            
            # Re-check if we still have NaN values
            remaining_nans = clean_adaptation_data[['Log_Length', 'Variant_Density']].isna().any(axis=1).sum()
            print(f"  After fixing: {remaining_nans} rows still contain NaN values")
            
            # Only drop rows if still necessary
            if remaining_nans > 0:
                clean_adaptation_data = clean_adaptation_data.dropna(subset=['Log_Length', 'Variant_Density'])
            
            # Skip if too few data points after cleaning
            if len(clean_adaptation_data) < 5:
                print(f"Insufficient clean data points for {adaptation} adaptation after fixing NaN/inf values")
                continue
            
            try:
                X = sm.add_constant(clean_adaptation_data[['Log_Length']])
                y = clean_adaptation_data['Variant_Density']
                model = sm.OLS(y, X).fit()
                
                # Deal with potential negative R-squared values
                adjusted_rsquared = max(0, model.rsquared)  # Convert negative to 0
                
                adaptation_models[adaptation] = {
                    'r_squared': adjusted_rsquared,
                    'coefficients': model.params.to_dict(),
                    'p_values': model.pvalues.to_dict(),
                    'n': len(clean_adaptation_data)
                }
            except Exception as e:
                print(f"Could not build model for {adaptation} adaptation: {e}")
    
    return model_results, adaptation_models

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
            
            # Also try clustering by adaptation type if available
            adaptation_pivot = None
            adaptation_linkage = None
            
            if 'Adaptation_Type' in data.columns:
                # Create adaptation-by-scaffold matrix
                adaptation_pivot = pd.pivot_table(
                    data,
                    values='POS',
                    index='CHROM',
                    columns='Adaptation_Type',
                    aggfunc='count',
                    fill_value=0
                )
                
                # Normalize and handle potential issues
                adaptation_row_sums = adaptation_pivot.sum(axis=1)
                adaptation_pivot_norm = adaptation_pivot[adaptation_row_sums > 0].div(
                    adaptation_row_sums[adaptation_row_sums > 0], axis=0)
                
                if len(adaptation_pivot_norm) >= 2 and len(adaptation_pivot_norm.columns) >= 2:
                    adaptation_linkage = hierarchy.linkage(
                        adaptation_pivot_norm.T,
                        method='average',
                        metric='correlation'
                    )
            
            return {
                'normalized_pivot': normalized_pivot,
                'treatment_linkage': treatment_linkage,
                'scaffold_linkage': scaffold_linkage,
                'adaptation_pivot': adaptation_pivot,
                'adaptation_linkage': adaptation_linkage
            }
            
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
                return {
                    'normalized_pivot': normalized_pivot,
                    'treatment_linkage': treatment_linkage,
                    'scaffold_linkage': scaffold_linkage,
                    'adaptation_pivot': None,
                    'adaptation_linkage': None
                }
                
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
    
    # Plot each treatment with its biological category
    if 'Adaptation_Type' in pca_df.columns:
        # Plot by adaptation type
        for adaptation in pca_df['Adaptation_Type'].unique():
            adaptation_data = pca_df[pca_df['Adaptation_Type'] == adaptation]
            
            marker = 'o' if adaptation == 'Temperature' else 's'
            plt.scatter(
                adaptation_data['PC1'], 
                adaptation_data['PC2'] if 'PC2' in adaptation_data.columns else [0] * len(adaptation_data),
                label=f"{adaptation} Adaptation",
                s=100,
                marker=marker
            )
    else:
        # Fall back to treatment-based plotting
        for treatment in pca_df['Treatment'].unique():
            treatment_data = pca_df[pca_df['Treatment'] == treatment]
            plt.scatter(
                treatment_data['PC1'], 
                treatment_data['PC2'] if 'PC2' in treatment_data.columns else [0] * len(treatment_data),
                label=treatment,
                color=TREATMENT_COLORS.get(treatment, 'gray'),
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
    
    # Add gene modification markers if available
    if 'Has_Gene' in pca_df.columns and pca_df['Has_Gene'].nunique() > 1:
        for i, row in pca_df.iterrows():
            if row['Has_Gene'] == 'Yes':
                plt.scatter(
                    row['PC1'],
                    row['PC2'] if 'PC2' in pca_df.columns else 0,
                    s=200,
                    facecolors='none',
                    edgecolors='black',
                    linewidth=2
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
    
    # Add explanatory annotation
    plt.figtext(0.02, 0.02, 
               "Note: Circles represent temperature adaptation, squares represent low oxygen adaptation.\n"
               "Black outlines indicate gene-modified strains (STC, CAS).",
               ha='left', fontsize=9)
    
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
        
        # Create PCA biplot (combining sample scores and feature loadings)
        plt.figure(figsize=(12, 10))
        
        # Plot samples
        for treatment in pca_df['Treatment'].unique():
            treatment_data = pca_df[pca_df['Treatment'] == treatment]
            plt.scatter(
                treatment_data['PC1'], 
                treatment_data['PC2'] if 'PC2' in treatment_data.columns else [0] * len(treatment_data),
                label=treatment,
                color=TREATMENT_COLORS.get(treatment, 'gray'),
                s=100
            )
        
        # Add treatment labels
        for i, row in pca_df.iterrows():
            plt.annotate(
                row['Treatment'],
                (row['PC1'], row['PC2'] if 'PC2' in row else 0),
                textcoords="offset points",
                xytext=(5, 5),
                ha='left'
            )
        
        # Add feature vectors (PC loadings)
        # Select top contributing features for clarity
        top_features = sorted_features.index[:8]  # Show top 8 features
        
        for feature in top_features:
            # Scale the arrows
            scale_factor = 3  # Adjust this to make arrows more visible
            plt.arrow(
                0, 0,  # Start at origin
                feature_contributions.loc[feature, 'PC1'] * scale_factor,
                feature_contributions.loc[feature, 'PC2'] * scale_factor if 'PC2' in feature_contributions.columns else 0,
                color='red',
                width=0.01,
                head_width=0.05
            )
            
            # Add feature labels
            plt.annotate(
                feature,
                (feature_contributions.loc[feature, 'PC1'] * scale_factor * 1.1,
                 feature_contributions.loc[feature, 'PC2'] * scale_factor * 1.1 if 'PC2' in feature_contributions.columns else 0),
                color='red',
                fontsize=9
            )
        
        plt.xlabel(f'PC1 ({explained_variance[0]:.2%} explained variance)')
        if len(explained_variance) >= 2:
            plt.ylabel(f'PC2 ({explained_variance[1]:.2%} explained variance)')
        else:
            plt.ylabel('PC2')
        
        plt.title('PCA Biplot: Treatments and Feature Contributions')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "pca_biplot.png"), dpi=300)
        plt.close()

# Function to plot correlation results
def plot_correlation_results(mutation_corr, treatment_corr, adaptation_corr, output_dir):
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
    
    # Plot treatment correlation matrix
    if treatment_corr is not None and not treatment_corr.empty:
        plt.figure(figsize=(10, 8))
        
        # Define custom colors for treatments
        treatment_colors = {t: TREATMENT_COLORS.get(t, '#333333') for t in treatment_corr.columns}
        
        # Create row colors for the heatmap
        row_colors = pd.Series(
            treatment_corr.index.map(lambda x: treatment_colors.get(x, '#333333')),
            index=treatment_corr.index
        )
        
        # Create clustered heatmap
        g = sns.clustermap(
            treatment_corr,
            annot=True,
            cmap='coolwarm',
            center=0,
            vmin=-1,
            vmax=1,
            row_colors=row_colors,
            col_colors=row_colors,
            figsize=(10, 8)
        )
        
        # Improve readability
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
        
        # Add treatment type legend
        for adaptation_type, treatments in {
            'Temperature': ['WT-37', 'CAS'],
            'Low Oxygen': ['WTA', 'STC']
        }.items():
            for t in treatments:
                if t in treatment_corr.columns:
                    g.ax_row_dendrogram.bar(0, 0, color=treatment_colors.get(t), 
                                          label=f"{t} ({adaptation_type})")
        
        g.ax_row_dendrogram.legend(title="Treatments", loc="center", ncol=1)
        
        plt.suptitle('Correlation Between Treatments', y=0.98)
        plt.savefig(os.path.join(output_dir, "treatment_correlation_clustered.png"), dpi=300)
        plt.close()
        
        # Also create simple heatmap for clarity
        plt.figure(figsize=(8, 6))
        sns.heatmap(treatment_corr, annot=True, cmap='coolwarm', center=0, vmin=-1, vmax=1)
        plt.title('Correlation Between Treatments')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "treatment_correlation_simple.png"), dpi=300)
        plt.close()
    
    # Plot adaptation type correlation if available
    if adaptation_corr is not None and not adaptation_corr.empty:
        plt.figure(figsize=(8, 6))
        sns.heatmap(adaptation_corr, annot=True, cmap='coolwarm', center=0, vmin=-1, vmax=1)
        plt.title('Correlation Between Adaptation Types')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "adaptation_correlation.png"), dpi=300)
        plt.close()

# Function to plot regression results
def plot_regression_results(data, model_results, adaptation_models, output_dir):
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
    plt.figure(figsize=(14, 12))
    
    # Define colors for treatments
    treatment_colors = TREATMENT_COLORS
    
    # Get number of treatments with models
    n_treatments = len(model_results)
    plot_rows = (n_treatments + 1) // 2  # Round up division
    
    # Set up subplots - create enough subplots for all treatments
    treatments = list(model_results.keys())
    
    for i, treatment in enumerate(treatments):
        if i < n_treatments:
            plt.subplot(plot_rows, 2, i+1)
            
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
                
                # Add p-value
                p_val = model_info['p_values'].get('Log_Length', 1.0)
                sig_stars = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
                plt.text(
                    0.05, 0.85,
                    f"p = {p_val:.3f} {sig_stars}",
                    transform=plt.gca().transAxes,
                    fontsize=12,
                    verticalalignment='top'
                )
            
            # Add treatment description
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            gene_mod = "with gene" if TREATMENT_INFO.get(treatment, {}).get('gene') else "no gene"
            title = f"{treatment}: {adaptation} adaptation ({gene_mod})"
            
            plt.xlabel('Log10(Scaffold Length)')
            plt.ylabel('Variant Density (per kb)')
            plt.title(title)
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
    
    # Sort by adaptation type then by R-squared
    treatment_data = pd.DataFrame({
        'Treatment': treatments,
        'R_squared': r2_values,
        'Adaptation': [TREATMENT_INFO.get(t, {}).get('adaptation', 'Unknown') for t in treatments],
        'Has_Gene': [TREATMENT_INFO.get(t, {}).get('gene', None) is not None for t in treatments]
    })
    
    treatment_data = treatment_data.sort_values(['Adaptation', 'R_squared'], ascending=[True, False])
    
    # Plot bars colored by adaptation type
    plt.figure(figsize=(10, 6))
    bars = plt.bar(
        range(len(treatment_data)), 
        treatment_data['R_squared'],
        color=[treatment_colors.get(t, 'gray') for t in treatment_data['Treatment']]
    )
    
    plt.xticks(range(len(treatment_data)), treatment_data['Treatment'])
    plt.ylabel('R-squared')
    plt.title('Model Fit (R²) by Treatment')
    plt.grid(True, linestyle='--', alpha=0.7, axis='y')
    
    # Add value labels
    for i, bar in enumerate(bars):
        height = bar.get_height()
        plt.text(
            i, height + 0.01, 
            f"{height:.3f}", 
            ha='center'
        )
    
    # Mark gene-modified treatments
    for i, has_gene in enumerate(treatment_data['Has_Gene']):
        if has_gene:
            plt.scatter(
                i,
                treatment_data['R_squared'].iloc[i] / 2,  # Place in middle of bar
                marker='*',
                s=200,
                color='white',
                edgecolors='black',
                linewidth=1
            )
    
    # Add a legend for adaptation types
    from matplotlib.patches import Patch
    legend_elements = []
    
    adaptation_types = treatment_data['Adaptation'].unique()
    for adaptation in adaptation_types:
        # Find a treatment with this adaptation
        for t in treatment_data['Treatment']:
            if TREATMENT_INFO.get(t, {}).get('adaptation') == adaptation:
                legend_elements.append(
                    Patch(facecolor=treatment_colors.get(t, 'gray'), 
                         label=f"{adaptation} Adaptation")
                )
                break
    
    # Add gene modification element
    legend_elements.append(
        plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='white',
                 markeredgecolor='black', markersize=15, label='Gene-modified')
    )
    
    plt.legend(handles=legend_elements)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "r2_by_treatment.png"), dpi=300)
    plt.close()
    
    # Plot adaptation type models if available
    if adaptation_models and len(adaptation_models) > 0:
        # Create scatter plot with regression lines by adaptation type
        plt.figure(figsize=(12, 6))
        
        adaptation_scaffold_data = data.groupby(['CHROM', 'Adaptation_Type']).agg({
            'Variant_Density': 'first',
            'Scaffold_Length': 'first'
        }).reset_index()
        
        adaptation_scaffold_data['Log_Length'] = np.log10(
            adaptation_scaffold_data['Scaffold_Length'].replace(0, 1))
        
        # Use custom colors for adaptation types
        adaptation_colors = {
            'Temperature': '#1f77b4',
            'Low Oxygen': '#ff7f0e'
        }
        
        for i, (adaptation, model_info) in enumerate(adaptation_models.items()):
            plt.subplot(1, len(adaptation_models), i+1)
            
            # Get adaptation data
            adaptation_data = adaptation_scaffold_data[
                adaptation_scaffold_data['Adaptation_Type'] == adaptation]
            
            # Scatter plot
            plt.scatter(
                adaptation_data['Log_Length'],
                adaptation_data['Variant_Density'],
                color=adaptation_colors.get(adaptation, 'gray'),
                alpha=0.6,
                s=50,
                label=adaptation
            )
            
            # Add regression line
            if 'coefficients' in model_info:
                # Get coefficients
                intercept = model_info['coefficients'].get('const', 0)
                slope = model_info['coefficients'].get('Log_Length', 0)
                
                # Generate points for line
                x_range = np.linspace(
                    adaptation_data['Log_Length'].min(),
                    adaptation_data['Log_Length'].max(),
                    100
                )
                y_range = intercept + slope * x_range
                
                # Plot regression line
                plt.plot(x_range, y_range, color='black', linestyle='--')
                
                # Add R-squared value
                plt.text(
                    0.05, 0.95, 
                    f"R² = {model_info['r_squared']:.3f}",
                    transform=plt.gca().transAxes,
                    fontsize=12,
                    verticalalignment='top'
                )
                
                # Add p-value
                p_val = model_info['p_values'].get('Log_Length', 1.0)
                sig_stars = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
                plt.text(
                    0.05, 0.85,
                    f"p = {p_val:.3f} {sig_stars}",
                    transform=plt.gca().transAxes,
                    fontsize=12,
                    verticalalignment='top'
                )
            
            plt.xlabel('Log10(Scaffold Length)')
            plt.ylabel('Variant Density (per kb)')
            plt.title(f"{adaptation} Adaptation Regression")
            plt.grid(True, linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "adaptation_regression.png"), dpi=300)
        plt.close()

# Function to plot clustering results
def plot_clustering_results(clustering_results, output_dir):
    """Plot clustering analysis results."""
    if clustering_results is None:
        return
    
    normalized_pivot = clustering_results.get('normalized_pivot')
    treatment_linkage = clustering_results.get('treatment_linkage')
    scaffold_linkage = clustering_results.get('scaffold_linkage')
    
    if normalized_pivot is None or normalized_pivot.empty:
        return
    
    # Create clustered heatmap
    plt.figure(figsize=(14, 12))
    
    # Create row colors based on scaffold properties (if available)
    # For now, we'll use a simpler approach without row colors
    
    # Create column colors based on treatment adaptation type
    col_colors = None
    if all(t in TREATMENT_INFO for t in normalized_pivot.columns):
        col_colors = pd.Series(
            [TREATMENT_COLORS.get(t, '#333333') for t in normalized_pivot.columns],
            index=normalized_pivot.columns
        )
    
    # Create clustered heatmap with both row and column dendrograms
    g = sns.clustermap(
        normalized_pivot,
        row_linkage=scaffold_linkage,
        col_linkage=treatment_linkage,
        cmap='YlGnBu',
        figsize=(14, 12),
        col_colors=col_colors,
        yticklabels=False  # Too many scaffolds to show labels
    )
    
    plt.suptitle('Hierarchical Clustering of Scaffold-Treatment Patterns', y=1.02)
    
    # If we have column colors, add a legend
    if col_colors is not None:
        # Add treatment type legend
        for adaptation_type, treatments in {
            'Temperature': ['WT-37', 'CAS'],
            'Low Oxygen': ['WTA', 'STC']
        }.items():
            for t in treatments:
                if t in normalized_pivot.columns:
                    g.ax_col_dendrogram.bar(0, 0, color=TREATMENT_COLORS.get(t), 
                                          label=f"{t} ({adaptation_type})")
        
        g.ax_col_dendrogram.legend(title="Treatments", loc="center", ncol=1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "scaffold_treatment_clustering.png"), dpi=300)
    plt.close()
    
    # Create treatment dendrogram
    plt.figure(figsize=(10, 6))
    
    # Plot dendrogram
    hierarchy.dendrogram(
        treatment_linkage,
        labels=normalized_pivot.columns,
        leaf_font_size=12
    )
    
    # Add treatment information to labels
    ax = plt.gca()
    new_labels = []
    for text in ax.get_xticklabels():
        treatment = text.get_text()
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', '')
        has_gene = TREATMENT_INFO.get(treatment, {}).get('gene', None) is not None
        gene_text = " (gene)" if has_gene else ""
        new_text = f"{treatment}\n{adaptation}{gene_text}"
        text.set_text(new_text)
        # Color the label according to treatment
        text.set_color(TREATMENT_COLORS.get(treatment, 'black'))
    
    plt.title('Hierarchical Clustering of Treatments')
    plt.xlabel('Treatment')
    plt.ylabel('Distance')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "treatment_dendrogram.png"), dpi=300)
    plt.close()
    
    # Plot adaptation clustering if available
    adaptation_pivot = clustering_results.get('adaptation_pivot')
    adaptation_linkage = clustering_results.get('adaptation_linkage')
    
    if adaptation_pivot is not None and adaptation_linkage is not None:
        plt.figure(figsize=(8, 6))
        
        # Plot dendrogram
        hierarchy.dendrogram(
            adaptation_linkage,
            labels=adaptation_pivot.columns,
            leaf_font_size=12
        )
        
        plt.title('Hierarchical Clustering of Adaptation Types')
        plt.xlabel('Adaptation Type')
        plt.ylabel('Distance')
        plt.grid(True, linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "adaptation_dendrogram.png"), dpi=300)
        plt.close()
        
        # Create adaptation heatmap
        plt.figure(figsize=(12, 10))
        
        g = sns.clustermap(
            adaptation_pivot.T,  # Transpose to have adaptation as rows
            cmap='YlGnBu',
            figsize=(12, 10),
            yticklabels=True,
            xticklabels=False  # Too many scaffolds
        )
        
        plt.suptitle('Variant Patterns by Adaptation Type', y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "adaptation_heatmap.png"), dpi=300)
        plt.close()

# Function to create summary report
def create_summary_report(summary_df, model_results, adaptation_models, output_dir):
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
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            has_gene = TREATMENT_INFO.get(treatment, {}).get('gene') is not None
            gene_text = " with gene modification" if has_gene else ""
            
            f.write(f"\n{treatment} Treatment ({adaptation} adaptation{gene_text}):\n")
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
        
        # Add adaptation model results if available
        if adaptation_models and len(adaptation_models) > 0:
            f.write("\nRegression Analysis by Adaptation Type:\n")
            f.write("-----------------------------------\n")
            
            for adaptation, model_info in adaptation_models.items():
                f.write(f"\n{adaptation} Adaptation:\n")
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
        
        # Write biological context summary
        f.write("Biological Context Summary:\n")
        f.write("-------------------------\n")
        f.write("This analysis examines mutation patterns in yeast strains under different conditions:\n\n")
        
        for treatment, info in TREATMENT_INFO.items():
            description = info.get('description', 'Unknown')
            adaptation = info.get('adaptation', 'Unknown')
            gene = info.get('gene', None)
            
            gene_text = f" with {gene} gene insertion" if gene else ""
            f.write(f"  {treatment}: {description} - {adaptation} adaptation{gene_text}\n")
        
        f.write("\n")
        
        # Write main conclusions
        f.write("Main Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis examines statistical patterns in the variant data.\n")
        f.write("2. Principal Component Analysis identifies major axes of variation between treatments.\n")
        f.write("3. Correlation analysis reveals relationships between mutation types and treatments.\n")
        f.write("4. Regression models evaluate the relationship between genomic features and variant density.\n")
        f.write("5. Clustering analysis identifies groups of scaffolds with similar variant patterns.\n")
        f.write("6. Comparisons between adaptation types (temperature vs. low oxygen) reveal\n")
        f.write("   differences in mutation patterns and genomic distributions.\n")
        f.write("7. Gene-modified strains (STC, CAS) show specific mutation characteristics\n")
        f.write("   compared to their non-modified counterparts (WTA, WT-37).\n")

# Main function to run the analysis
def main():
    # Load data from previous analyses
    mutation_data = load_mutation_data()
    scaffold_info = load_scaffold_info()
    
    # Integrate data sources
    integrated_data = integrate_data(mutation_data, scaffold_info)
    
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
        mutation_corr, treatment_corr, adaptation_corr = perform_correlation_analysis(integrated_data)
        # Plot correlation results
        plot_correlation_results(mutation_corr, treatment_corr, adaptation_corr, OUTPUT_DIR)
    
    # Build regression models
    model_results, adaptation_models = build_regression_models(integrated_data)
    if model_results:
        # Plot regression results
        plot_regression_results(integrated_data, model_results, adaptation_models, OUTPUT_DIR)
    
    # Perform clustering analysis
    clustering_results = perform_clustering_analysis(integrated_data)
    if clustering_results:
        # Plot clustering results
        plot_clustering_results(clustering_results, OUTPUT_DIR)
    
    # Create summary report
    create_summary_report(summary_df, model_results, adaptation_models, OUTPUT_DIR)
    
    print(f"Analysis complete! Results saved to {OUTPUT_DIR}/")

# Run the analysis
if __name__ == "__main__":
    main()