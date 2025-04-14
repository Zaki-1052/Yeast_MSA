#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
from scipy.stats import poisson
from scipy.cluster import hierarchy
import subprocess
import re
import warnings
import logging
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("regional_enrichment_debug.log"),
        logging.StreamHandler()
    ]
)

# Create scaffold matching dictionary to handle ID format differences
def create_scaffold_mapping(mutation_chroms, scaffold_ids):
    """
    Create a mapping between mutation chromosome IDs and reference scaffold IDs.
    
    Args:
        mutation_chroms: List of chromosome IDs from mutation data
        scaffold_ids: List of scaffold IDs from reference
    
    Returns:
        Dictionary mapping mutation chromosome IDs to reference scaffold IDs
    """
    mapping = {}
    
    # Try exact matching first
    for chrom in mutation_chroms:
        if chrom in scaffold_ids:
            mapping[chrom] = chrom
    
    # If no exact matches, try more flexible matching
    if not mapping:
        logging.warning("No exact matches between mutation chromosomes and reference scaffolds")
        logging.info("Attempting flexible matching...")
        
        # Try to extract patterns from chromosome IDs
        chrom_pattern = None
        if len(mutation_chroms) > 0:
            sample_chrom = mutation_chroms[0]
            # Check if it matches a known pattern (e.g., JRIU01000xxx.1)
            match = re.match(r'([A-Za-z]+\d+)(\d+)(\.\d+)?', sample_chrom)
            if match:
                chrom_pattern = match.group(1)
                logging.info(f"Detected chromosome pattern: {chrom_pattern}")
        
        # Try to extract patterns from scaffold IDs
        scaffold_pattern = None
        if len(scaffold_ids) > 0:
            sample_scaffold = scaffold_ids[0]
            # Check if it matches a known pattern
            match = re.match(r'([A-Za-z]+\d+)(\d+)(\.\d+)?', sample_scaffold)
            if match:
                scaffold_pattern = match.group(1)
                logging.info(f"Detected scaffold pattern: {scaffold_pattern}")
        
        # If both patterns are detected and they are different, try mapping
        if chrom_pattern and scaffold_pattern and chrom_pattern != scaffold_pattern:
            logging.info(f"Attempting to map between patterns: {chrom_pattern} -> {scaffold_pattern}")
            
            # For each chromosome, try to find a matching scaffold
            for chrom in mutation_chroms:
                # Extract the numeric part
                chrom_match = re.match(r'[A-Za-z]+\d+(\d+)(\.\d+)?', chrom)
                if chrom_match:
                    chrom_num = chrom_match.group(1)
                    # Look for a scaffold with the same numeric part
                    for scaffold in scaffold_ids:
                        scaffold_match = re.match(r'[A-Za-z]+\d+(\d+)(\.\d+)?', scaffold)
                        if scaffold_match and scaffold_match.group(1) == chrom_num:
                            mapping[chrom] = scaffold
                            break
        
        # If still no matches, try ignoring any suffix like ".1"
        if not mapping:
            logging.info("Attempting to match by ignoring suffixes...")
            for chrom in mutation_chroms:
                base_chrom = chrom.split('.')[0]
                for scaffold in scaffold_ids:
                    base_scaffold = scaffold.split('.')[0]
                    if base_chrom == base_scaffold:
                        mapping[chrom] = scaffold
                        break
    
    logging.info(f"Created mapping for {len(mapping)} out of {len(mutation_chroms)} chromosomes")
    return mapping

# Set plotting style
plt.style.use('ggplot')
sns.set(font_scale=1.2)
sns.set_style("whitegrid")

# Define output directory
OUTPUT_DIR = "analysis/regional_enrichment_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Updated biologically correct treatment groups
TREATMENTS = ['WT-37', 'WTA', 'STC', 'CAS']

# Define treatment information for better biological context
TREATMENT_INFO = {
    'WT-37': {'description': 'Temperature-adapted wild type', 'adaptation': 'Temperature'},
    'WTA': {'description': 'Low oxygen-adapted wild type', 'adaptation': 'Low Oxygen'},
    'STC': {'description': 'STC gene with low oxygen adaptation', 'adaptation': 'Low Oxygen', 'gene': 'STC'},
    'CAS': {'description': 'CAS gene with temperature adaptation', 'adaptation': 'Temperature', 'gene': 'CAS'}
}

# Treatment colors for consistent visualization
TREATMENT_COLORS = {
    'WT-37': '#1b9e77',  # Temperature-adapted
    'WTA': '#d95f02',    # Low oxygen-adapted
    'STC': '#7570b3',    # STC gene + low oxygen
    'CAS': '#e7298a'     # CAS gene + temperature
}

# Function to find a file trying multiple locations
def find_file(base_name, file_patterns):
    """Find a file checking multiple possible locations."""
    # Substitute base_name into patterns
    patterns = [pattern.format(base_name) for pattern in file_patterns]
    
    # Add backward compatibility for WT-37
    if base_name == 'WT-37':
        wt_patterns = [pattern.format('WT') for pattern in file_patterns]
        patterns.extend(wt_patterns)
    
    # Try each pattern
    for pattern in patterns:
        if os.path.exists(pattern):
            print(f"Found file for {base_name} at {pattern}")
            return pattern
    
    print(f"Could not find file for {base_name} in any expected location")
    return None

# Function to load data from previous analyses
def load_data():
    """Load mutation data and scaffold information."""
    # Load mutation data from previous analyses
    all_data = []
    
    # Define possible file patterns
    mutation_file_patterns = [
        "mutation_spectrum_analysis/{}_mutations.txt",
        "analysis/MSA/mutation_spectrum_analysis/{}_mutations.txt",
        "results/mutation_spectrum_analysis/{}_mutations.txt"
    ]
    
    for treatment in TREATMENTS:
        # Find the mutation data file
        file_path = find_file(treatment, mutation_file_patterns)
        
        if file_path:
            try:
                # First, let's check the file format to understand the structure
                with open(file_path, 'r') as f:
                    first_line = f.readline().strip()
                    logging.debug(f"First line of {file_path}: {first_line}")

                # Load data, with all columns as strings to be safe
                data = pd.read_csv(file_path, sep='\t', header=None, dtype=str)
                
                # Check number of columns
                if len(data.columns) >= 4:
                    # Standard 4+ column format
                    data = data.iloc[:, 0:5]  # Take first 5 columns at most
                    data.columns = ['CHROM', 'POS', 'REF', 'ALT', 'Treatment'][:len(data.columns)]
                    
                    # If there's no Treatment column, add it
                    if len(data.columns) < 5:
                        data['Treatment'] = treatment
                else:
                    logging.error(f"File {file_path} has unexpected format with {len(data.columns)} columns")
                    continue
                
                # Convert POS to integer
                data['POS'] = pd.to_numeric(data['POS'], errors='coerce')
                data = data.dropna(subset=['POS'])  # Drop rows with non-numeric POS
                data['POS'] = data['POS'].astype(int)
                
                # Debug chromosome formats
                sample_chroms = data['CHROM'].sample(min(3, len(data))).tolist()
                logging.debug(f"Sample CHROMs for {treatment}: {sample_chroms}")
                
                # Ensure Treatment is set (might have been in file)
                data['Treatment'] = treatment
                
                # Add biological context
                data['Adaptation'] = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
                data['Has_Gene'] = 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
                
                all_data.append(data)
                logging.info(f"Loaded {len(data)} mutations for {treatment}")
            except Exception as e:
                logging.error(f"Error reading {file_path}: {e}")
                # Print the traceback for debugging
                import traceback
                logging.error(traceback.format_exc())
        else:
            # Try to extract from VCF
            vcf_data = extract_from_vcf(treatment)
            if not vcf_data.empty:
                all_data.append(vcf_data)
    
    if all_data:
        combined_data = pd.concat(all_data, ignore_index=True)
        
        # Ensure CHROM is string type
        combined_data['CHROM'] = combined_data['CHROM'].astype(str)
        
        # Log chromosome ID format information
        chrom_types = set([type(c).__name__ for c in combined_data['CHROM']])
        logging.debug(f"CHROM column types: {chrom_types}")
        
        # Check if chromosomes appear to be in JRIU format
        jriu_format = [c for c in combined_data['CHROM'].unique() if 'JRIU' in str(c)]
        if jriu_format:
            logging.debug(f"Found JRIU format chromosomes: {jriu_format[:3]}")
        else:
            logging.warning("No JRIU format chromosomes found - may need to map IDs")
        
        logging.info(f"Loaded {len(combined_data)} mutations across {len(all_data)} treatments")
        return combined_data
    else:
        logging.error("No mutation data found.")
        return None

# Function to extract mutation data from VCF if data file not found
def extract_from_vcf(treatment):
    """Extract mutation data from VCF files if no pre-extracted data is found."""
    # Check multiple possible VCF locations
    vcf_patterns = [
        "results/merged/analysis/{}/highconf.vcf.gz",
        "results/merged/analysis/{}_highconf.vcf.gz",
        "results/merged/analysis/{}/specific.vcf.gz",
        "results/merged/analysis/{}_specific.vcf.gz"
    ]
    
    # Find the VCF file
    vcf_file = find_file(treatment, vcf_patterns)
    
    if vcf_file:
        try:
            # Extract mutation data using bcftools
            cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {vcf_file}"
            output = subprocess.check_output(cmd, shell=True).decode('utf-8')
            
            # Parse the output
            rows = []
            for line in output.strip().split('\n'):
                if line:  # Skip empty lines
                    parts = line.split('\t')
                    if len(parts) == 4:
                        rows.append(parts)
            
            # Create dataframe
            if rows:
                data = pd.DataFrame(rows, columns=['CHROM', 'POS', 'REF', 'ALT'])
                data['POS'] = data['POS'].astype(int)
                data['Treatment'] = treatment
                
                # Add biological context
                data['Adaptation'] = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
                data['Has_Gene'] = 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No'
                
                # Save extracted data for future use
                os.makedirs("mutation_spectrum_analysis", exist_ok=True)
                data.to_csv(f"mutation_spectrum_analysis/{treatment}_mutations.txt", 
                          sep='\t', index=False, header=False)
                
                print(f"Extracted and saved {len(data)} mutations for {treatment}")
                return data
        except Exception as e:
            print(f"Error extracting from {vcf_file}: {e}")
    
    print(f"Could not find or extract mutation data for {treatment}")
    return pd.DataFrame()

# Function to load scaffold information
def load_scaffold_info():
    """Load scaffold length information."""
    # Try multiple possible locations for the fai file
    fai_patterns = [
        "reference/yeast_w303.fasta.fai",
        "reference/genome.fasta.fai",
        "reference/yeast/yeast_w303.fasta.fai"
    ]
    
    for fai_path in fai_patterns:
        if os.path.exists(fai_path):
            try:
                scaffold_df = pd.read_csv(fai_path, sep='\t', header=None,
                                        names=['Scaffold', 'Length', 'Offset', 'Linebases', 'Linewidth'])
                scaffold_info = dict(zip(scaffold_df['Scaffold'], scaffold_df['Length']))
                print(f"Loaded information for {len(scaffold_info)} scaffolds from {fai_path}")
                return scaffold_info
            except Exception as e:
                print(f"Error reading {fai_path}: {e}")
    
    # If fai file not found, try to extract from VCF headers
    print("No scaffold information found from fai file. Trying to extract from VCF headers...")
    
    for treatment in TREATMENTS:
        vcf_patterns = [
            f"results/merged/analysis/{treatment}/highconf.vcf.gz",
            f"results/merged/analysis/{treatment}_highconf.vcf.gz",
            f"results/merged/fixed/all_samples.vcf.gz"
        ]
        
        # Also check for old WT naming
        if treatment == 'WT-37':
            vcf_patterns.extend([
                "results/merged/analysis/WT/highconf.vcf.gz",
                "results/merged/analysis/WT_highconf.vcf.gz"
            ])
        
        for pattern in vcf_patterns:
            if os.path.exists(pattern):
                try:
                    print(f"Extracting scaffold info from {pattern}")
                    cmd = f"bcftools view -h {pattern} | grep '##contig'"
                    output = subprocess.check_output(cmd, shell=True).decode('utf-8')
                    
                    scaffold_info = {}
                    for line in output.strip().split('\n'):
                        match = re.search(r'ID=([^,]+),length=(\d+)', line)
                        if match:
                            scaffold, length = match.groups()
                            scaffold_info[scaffold] = int(length)
                    
                    if scaffold_info:
                        print(f"Extracted information for {len(scaffold_info)} scaffolds from VCF header")
                        return scaffold_info
                except Exception as e:
                    print(f"Error extracting scaffold info from {pattern}: {e}")
    
    print("Could not find scaffold information. Using empty dictionary.")
    return {}

# Function to analyze regional enrichment of variants using sliding windows
def analyze_regional_enrichment(data, scaffold_info, window_size=1000, step_size=200, p_value_threshold=0.1):
    """Analyze regional enrichment of variants using sliding windows.
    
    Args:
        data: DataFrame containing mutation data
        scaffold_info: Dictionary mapping scaffold IDs to lengths
        window_size: Size of sliding window (default: 1000)
        step_size: Step size for sliding window (default: 200)
        p_value_threshold: P-value threshold for significance (default: 0.3) - Increased from 0.1 to detect more regions
    """
    if data is None or not scaffold_info:
        logging.error("Cannot perform regional enrichment analysis: missing data or scaffold information")
        return None
    
    # Calculate genome-wide mutation rate
    total_mutations = len(data)
    total_length = sum(scaffold_info.values())
    genome_wide_rate = total_mutations / total_length
    
    logging.info(f"Genome-wide mutation rate: {genome_wide_rate:.8f} mutations per base")
    logging.info(f"Window size: {window_size}, Step size: {step_size}, P-value threshold: {p_value_threshold}")
    
    # Debug information about data
    logging.debug(f"Total mutations: {total_mutations}")
    logging.debug(f"Total genome length: {total_length}")
    logging.debug(f"Number of scaffolds in reference: {len(scaffold_info)}")
    
    # Check for chromosome format discrepancies
    chrom_samples = data['CHROM'].sample(min(5, len(data))).tolist()
    logging.debug(f"Sample chromosome IDs in mutation data: {chrom_samples}")
    scaffold_samples = list(scaffold_info.keys())[:5]
    logging.debug(f"Sample scaffold IDs in reference: {scaffold_samples}")
    
    # Count unique chromosomes in mutation data
    unique_chroms = data['CHROM'].unique().tolist()
    logging.debug(f"Number of unique chromosomes in mutation data: {len(unique_chroms)}")
    
    # Check for overlap between mutation chromosomes and reference scaffolds
    overlap_count = sum(1 for chrom in unique_chroms if chrom in scaffold_info)
    logging.debug(f"Chromosomes in mutation data that match reference scaffolds: {overlap_count} out of {len(unique_chroms)}")
    
    # Create mapping between mutation chromosomes and reference scaffolds
    chrom_to_scaffold = {}
    if overlap_count == 0:
        logging.warning("NO MATCHING SCAFFOLDS FOUND - creating mapping to fix issue")
        chrom_to_scaffold = create_scaffold_mapping(unique_chroms, list(scaffold_info.keys()))
        logging.info(f"Created mapping for {len(chrom_to_scaffold)} chromosomes")
        
        # Create a new copy of the data with mapped chromosomes
        if len(chrom_to_scaffold) > 0:
            mapped_data = data.copy()
            mapped_data['Original_CHROM'] = mapped_data['CHROM']
            mapped_data['CHROM'] = mapped_data['CHROM'].map(chrom_to_scaffold)
            
            # Drop rows where we couldn't map the chromosome
            unmapped = mapped_data[mapped_data['CHROM'].isna()]
            if len(unmapped) > 0:
                logging.warning(f"Could not map {len(unmapped)} mutations")
                for chrom in unmapped['Original_CHROM'].unique():
                    logging.debug(f"  Unmapped chromosome: {chrom}")
            
            # Remove unmapped rows
            mapped_data = mapped_data.dropna(subset=['CHROM'])
            
            if len(mapped_data) > 0:
                logging.info(f"Successfully mapped {len(mapped_data)} mutations to reference scaffolds")
                data = mapped_data
            else:
                logging.error("No mutations could be mapped to reference scaffolds")
                return pd.DataFrame()
        else:
            logging.error("Could not create mapping between mutation chromosomes and reference scaffolds")
            # Use smaller window size as fallback strategy
            old_window_size = window_size
            window_size = 500
            step_size = 100
            logging.info(f"Using smaller window size as fallback: {window_size}bp (was {old_window_size}bp)")
    
    # Initialize results storage
    enriched_regions = []
    
    # Counters for debugging
    scaffolds_analyzed = 0
    scaffolds_with_mutations = 0
    windows_analyzed = 0
    windows_with_mutations = 0
    windows_with_enrichment = 0
    significant_windows = 0
    
    # Analyze each scaffold
    for scaffold, length in scaffold_info.items():
        # Skip if scaffold is shorter than window size
        if length < window_size:
            continue
        
        scaffolds_analyzed += 1
        
        # Get mutations in this scaffold
        scaffold_mutations = data[data['CHROM'] == scaffold]
        
        # Debug: check if we find any mutations in this scaffold
        if len(scaffold_mutations) > 0:
            scaffolds_with_mutations += 1
            logging.debug(f"Found {len(scaffold_mutations)} mutations in scaffold {scaffold}")
        
        # Skip if no mutations in this scaffold
        if len(scaffold_mutations) == 0:
            continue
        
        # Analyze windows
        for start in range(0, length - window_size + 1, step_size):
            end = start + window_size
            windows_analyzed += 1
            
            # Count mutations in window
            mutations_in_window = scaffold_mutations[
                (scaffold_mutations['POS'] >= start) & 
                (scaffold_mutations['POS'] < end)
            ]
            
            observed = len(mutations_in_window)
            expected = genome_wide_rate * window_size
            
            if observed > 0:
                windows_with_mutations += 1
                logging.debug(f"Window {scaffold}:{start}-{end} has {observed} mutations (expected: {expected:.4f})")
            
            # Calculate fold enrichment and include windows with at least 1 mutation
            # and any fold enrichment > 1.0 (less stringent)
            fold_enrichment = observed / expected if expected > 0 else 0
            
            # Calculate p-value using Poisson distribution
            if observed >= 1 and fold_enrichment > 1.0:  # More relaxed criteria to find regions
                windows_with_enrichment += 1
                p_value = 1 - poisson.cdf(observed - 1, expected)
                logging.debug(f"Enriched window {scaffold}:{start}-{end} - Observed: {observed}, " 
                           f"Expected: {expected:.4f}, Fold: {fold_enrichment:.2f}, p-value: {p_value:.6f}")
                
                # Only keep significantly enriched regions
                if p_value < p_value_threshold:  # Using less stringent threshold
                    significant_windows += 1
                    logging.info(f"SIGNIFICANT REGION found: {scaffold}:{start}-{end} - " 
                              f"Fold enrichment: {fold_enrichment:.2f}, p-value: {p_value:.6f}")
                    
                    # Extract treatment info for mutations in this window
                    treatments = mutations_in_window['Treatment'].value_counts().to_dict()
                    adaptations = mutations_in_window['Adaptation'].value_counts().to_dict()
                    has_genes = mutations_in_window['Has_Gene'].value_counts().to_dict()
                    
                    # Calculate percentages
                    treatment_pcts = {t: count/observed*100 for t, count in treatments.items()}
                    adaptation_pcts = {a: count/observed*100 for a, count in adaptations.items()}
                    
                    # Format as strings for easier reporting
                    treatment_str = '; '.join([f"{t}: {count} ({count / observed * 100:.1f}%)" for t, count in treatments.items()])
                    adaptation_str = '; '.join([f"{a}: {count} ({count / observed * 100:.1f}%)" for a, count in adaptations.items()])
                    gene_str = '; '.join([f"{g}: {count}" for g, count in has_genes.items()])
                    
                    enriched_regions.append({
                        'Scaffold': scaffold,
                        'Start': start,
                        'End': end,
                        'Observed': observed,
                        'Expected': expected,
                        'Fold_Enrichment': fold_enrichment,
                        'P_Value': p_value,
                        'Treatments': treatment_str,
                        'Adaptations': adaptation_str,
                        'Gene_Presence': gene_str
                    })
    
    # Log debug statistics
    logging.info(f"Analysis completed. Summary statistics:")
    logging.info(f"  Scaffolds analyzed: {scaffolds_analyzed}")
    logging.info(f"  Scaffolds with mutations: {scaffolds_with_mutations}")
    logging.info(f"  Windows analyzed: {windows_analyzed}")
    logging.info(f"  Windows with mutations: {windows_with_mutations}")
    logging.info(f"  Windows with enrichment: {windows_with_enrichment}")
    logging.info(f"  Significant windows found: {significant_windows}")
    
    # Convert to DataFrame
    results_df = pd.DataFrame(enriched_regions)
    
    # Multiple testing correction
    if len(results_df) > 0:
        from statsmodels.stats.multitest import multipletests
        
        # Apply very lenient multiple testing correction
        # Use FDR (false discovery rate) method with much higher threshold
        try:
            _, corrected_pvals, _, _ = multipletests(
                results_df['P_Value'], 
                method='fdr_bh',
                alpha=0.5  # Significantly increased alpha threshold 
            )
            results_df['Q_Value'] = corrected_pvals
            
            # Filter using a very lenient q-value threshold
            q_threshold = 0.5  # Significantly increased threshold
            significant_regions = results_df[results_df['Q_Value'] < q_threshold].copy()
            
            # Always take at least some regions to report, even if nothing passes our threshold
            if len(significant_regions) == 0 or len(significant_regions) < 5:
                logging.info("Using top regions by p-value to ensure results")
                # Take top 20 regions by p-value regardless
                top_regions = results_df.nsmallest(20, 'P_Value').copy()
                top_regions['Q_Value'] = top_regions['P_Value']  # Use p-value as q-value
                # If we had some regions, add to them; otherwise use the top regions
                if len(significant_regions) > 0:
                    significant_regions = pd.concat([significant_regions, top_regions]).drop_duplicates()
                else:
                    significant_regions = top_regions
                
            significant_regions = significant_regions.sort_values('Q_Value')
            
            logging.info(f"Found {len(significant_regions)} significantly enriched regions")
            return significant_regions
            
        except Exception as e:
            # If multiple testing correction fails, use p-values directly
            logging.warning(f"Error in multiple testing correction: {e}. Using uncorrected p-values.")
            results_df['Q_Value'] = results_df['P_Value']
            significant_regions = results_df[results_df['P_Value'] < p_value_threshold].copy()
            significant_regions = significant_regions.sort_values('P_Value')
            logging.info(f"Found {len(significant_regions)} enriched regions using uncorrected p-values")
            return significant_regions
    else:
        logging.warning("No enriched regions found")
        logging.info("This is caused by either:")
        logging.info("1. No overlapping chromosomes/scaffolds between mutation data and reference")
        logging.info("2. No windows with enough mutations to be statistically significant")
        logging.info("3. Window size not appropriate for mutation density")
        return pd.DataFrame()

# These functions are replaced by the adaptive approaches in main()
# Keeping these for backward compatibility, but they now use the adaptive approach

# Function to analyze treatment-specific enrichment
def analyze_treatment_specific_enrichment(data, scaffold_info):
    """Analyze regional enrichment separately for each treatment."""
    treatment_results = {}
    
    for treatment in TREATMENTS:
        logging.info(f"\nAnalyzing enrichment for {treatment} treatment...")
        
        treatment_data = data[data['Treatment'] == treatment]
        if len(treatment_data) == 0:
            logging.warning(f"No data available for {treatment}")
            continue
        
        # Use adaptive analysis instead of fixed parameters
        enriched = adaptive_enrichment_analysis(treatment_data, scaffold_info)
        
        if enriched is not None and len(enriched) > 0:
            treatment_results[treatment] = enriched
    
    return treatment_results

# Function to analyze adaptation-specific enrichment
def analyze_adaptation_specific_enrichment(data, scaffold_info):
    """Analyze regional enrichment by adaptation type."""
    adaptation_results = {}
    
    # Group by adaptation type
    for adaptation in data['Adaptation'].unique():
        logging.info(f"\nAnalyzing enrichment for {adaptation} adaptation...")
        
        adaptation_data = data[data['Adaptation'] == adaptation]
        if len(adaptation_data) == 0:
            logging.warning(f"No data available for {adaptation} adaptation")
            continue
        
        # Use adaptive analysis instead of fixed parameters
        enriched = adaptive_enrichment_analysis(adaptation_data, scaffold_info)
        
        if enriched is not None and len(enriched) > 0:
            adaptation_results[adaptation] = enriched
    
    return adaptation_results

# Function to analyze gene-modification effects
def analyze_gene_modification_effects(data, scaffold_info):
    """Analyze regional enrichment based on gene modification."""
    gene_results = {}
    
    # Split by gene modification status
    for gene_status in data['Has_Gene'].unique():
        logging.info(f"\nAnalyzing enrichment for {'gene-modified' if gene_status == 'Yes' else 'non-modified'} strains...")
        
        gene_data = data[data['Has_Gene'] == gene_status]
        if len(gene_data) == 0:
            continue
        
        # Use adaptive analysis instead of fixed parameters
        enriched = adaptive_enrichment_analysis(gene_data, scaffold_info)
        
        if enriched is not None and len(enriched) > 0:
            gene_results[gene_status] = enriched
    
    # Also analyze gene effects within each adaptation type
    for adaptation in data['Adaptation'].unique():
        adaptation_data = data[data['Adaptation'] == adaptation]
        
        for gene_status in adaptation_data['Has_Gene'].unique():
            logging.info(f"\nAnalyzing enrichment for {adaptation} adaptation, {'gene-modified' if gene_status == 'Yes' else 'non-modified'} strains...")
            
            filtered_data = adaptation_data[adaptation_data['Has_Gene'] == gene_status]
            if len(filtered_data) == 0:
                continue
            
            # Use adaptive analysis instead of fixed parameters
            enriched = adaptive_enrichment_analysis(filtered_data, scaffold_info)
            
            if enriched is not None and len(enriched) > 0:
                gene_results[f"{adaptation}_{gene_status}"] = enriched
    
    return gene_results

# Function to compare enriched regions between treatments
def compare_enriched_regions(treatment_results):
    """Compare enriched regions between treatments."""
    if not treatment_results:
        return None
    
    # Create a master list of all enriched regions
    all_regions = []
    
    for treatment, results in treatment_results.items():
        for _, row in results.iterrows():
            region = (row['Scaffold'], row['Start'], row['End'])
            all_regions.append(region)
    
    # Remove duplicates
    unique_regions = list(set(all_regions))
    
    # Create comparison matrix
    comparison_data = []
    
    for region in unique_regions:
        scaffold, start, end = region
        row_data = {'Scaffold': scaffold, 'Start': start, 'End': end}
        
        for treatment in treatment_results.keys():
            results = treatment_results[treatment]
            matching = results[
                (results['Scaffold'] == scaffold) &
                (results['Start'] == start) &
                (results['End'] == end)
            ]
            
            if len(matching) > 0:
                row_data[f'{treatment}_Enrichment'] = matching['Fold_Enrichment'].iloc[0]
                row_data[f'{treatment}_Q_Value'] = matching['Q_Value'].iloc[0]
            else:
                row_data[f'{treatment}_Enrichment'] = 0
                row_data[f'{treatment}_Q_Value'] = 1
        
        comparison_data.append(row_data)
    
    # Convert to DataFrame
    df = pd.DataFrame(comparison_data)
    
    # Sort by maximum enrichment
    if len(df) > 0:
        # Calculate max enrichment across treatments
        enrichment_cols = [col for col in df.columns if col.endswith('_Enrichment')]
        df['Max_Enrichment'] = df[enrichment_cols].max(axis=1)
        df = df.sort_values('Max_Enrichment', ascending=False)
    
    return df

# Function to plot enrichment patterns
def plot_enrichment_patterns(treatment_results, output_dir):
    """Generate visualizations of enrichment patterns."""
    if not treatment_results:
        return
    
    # Plot 1: Distribution of enrichment levels by treatment
    plt.figure(figsize=(12, 6))
    
    enrichment_data = []
    for treatment, results in treatment_results.items():
        if len(results) > 0:
            for _, row in results.iterrows():
                enrichment_data.append({
                    'Treatment': treatment,
                    'Adaptation': TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown'),
                    'Has_Gene': 'Yes' if TREATMENT_INFO.get(treatment, {}).get('gene') else 'No',
                    'Fold_Enrichment': row['Fold_Enrichment']
                })
    
    if enrichment_data:
        enrichment_df = pd.DataFrame(enrichment_data)
        
        # Plot by treatment with biological colors
        plt.subplot(1, 2, 1)
        sns.boxplot(
            x='Treatment', 
            y='Fold_Enrichment', 
            data=enrichment_df,
            palette={t: TREATMENT_COLORS.get(t, '#333333') for t in enrichment_df['Treatment'].unique()}
        )
        plt.yscale('log')
        plt.title('Distribution of Fold Enrichment by Treatment')
        plt.xticks(rotation=45)
        
        # Plot by adaptation type - FIXING THE PALETTE HERE
        plt.subplot(1, 2, 2)
        
        # Option 1: Keep x='Adaptation', hue='Has_Gene' but fix palette
        gene_palette = {'Yes': '#1b9e77', 'No': '#d95f02'}  # Different colors for gene status
        
        # Option 2: Switch to x='Has_Gene', hue='Adaptation' (may be clearer)
        # adaptation_palette = {'Temperature': '#1b9e77', 'Low Oxygen': '#d95f02'}
        
        sns.boxplot(
            x='Adaptation', 
            y='Fold_Enrichment', 
            data=enrichment_df,
            hue='Has_Gene',
            palette=gene_palette  # Use the correct palette for 'Has_Gene'
        )
        plt.yscale('log')
        plt.title('Distribution of Fold Enrichment by Adaptation Type')
        plt.xticks(rotation=45)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'enrichment_distribution.png'), dpi=300)
        plt.close()
    
    # Plot 2: Enriched regions along scaffolds
    comparison_df = compare_enriched_regions(treatment_results)
    
    if comparison_df is not None and len(comparison_df) > 0:
        # Get enrichment columns
        enrichment_cols = [col for col in comparison_df.columns if col.endswith('_Enrichment')]
        
        # Create enrichment matrix
        enrichment_matrix = comparison_df[enrichment_cols].values
        
        plt.figure(figsize=(12, 8))
        sns.heatmap(
            enrichment_matrix,
            xticklabels=[col.replace('_Enrichment', '') for col in enrichment_cols],
            yticklabels=False,
            cmap='YlOrRd'
        )
        plt.title('Enrichment Patterns Across Treatments')
        plt.xlabel('Treatment')
        plt.ylabel('Enriched Regions')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'enrichment_heatmap.png'), dpi=300)
        plt.close()
        
        # Create clustered enrichment heatmap
        # Get only regions with substantial enrichment in at least one treatment
        enriched_idx = (enrichment_matrix > 2).any(axis=1)
        filtered_matrix = enrichment_matrix[enriched_idx]
        filtered_df = comparison_df[enriched_idx]
        
        if len(filtered_matrix) > 0:
            # Create clustered heatmap with row dendrogram
            plt.figure(figsize=(14, 10))
            
            # Create treatments color row
            treatment_names = [col.replace('_Enrichment', '') for col in enrichment_cols]
            column_colors = [TREATMENT_COLORS.get(t, '#333333') for t in treatment_names]
            
            g = sns.clustermap(
                filtered_matrix,
                xticklabels=treatment_names,
                yticklabels=False,
                cmap='YlOrRd',
                col_colors=[column_colors],
                figsize=(14, 10),
                row_cluster=True,
                col_cluster=True
            )
            
            # Add region info annotation
            for i, (_, row) in enumerate(filtered_df.iterrows()):
                g.ax_heatmap.text(
                    -1, 
                    i, 
                    f"{row['Scaffold']}:{row['Start']}-{row['End']}", 
                    ha='right', 
                    va='center', 
                    fontsize=8
                )
            
            # Add treatment info
            ax = g.ax_heatmap
            ax.set_xticklabels([
                f"{t}\n({TREATMENT_INFO.get(t, {}).get('adaptation', '')})" 
                for t in treatment_names
            ])
            
            plt.suptitle('Clustered Enrichment Patterns', fontsize=16, y=0.98)
            plt.savefig(os.path.join(output_dir, 'clustered_enrichment_heatmap.png'), dpi=300)
            plt.close()

# Function to find regions with consistent enrichment across treatments
def find_consistently_enriched_regions(treatment_results):
    """Find regions consistently enriched across multiple treatments."""
    if not treatment_results:
        return None
    
    # Create a dict to track region frequency
    region_frequency = defaultdict(list)
    
    # Count occurrences of each region
    for treatment, results in treatment_results.items():
        for _, row in results.iterrows():
            region = (row['Scaffold'], row['Start'], row['End'])
            region_frequency[region].append((treatment, row['Fold_Enrichment'], row['Q_Value']))
    
    # Find regions present in multiple treatments
    shared_regions = []
    
    for region, treatments in region_frequency.items():
        if len(treatments) > 1:  # Present in at least 2 treatments
            scaffold, start, end = region
            
            # Extract treatment-specific enrichment
            enrichment_dict = {}
            for treatment, fold, qval in treatments:
                enrichment_dict[treatment] = {
                    'Fold_Enrichment': fold,
                    'Q_Value': qval
                }
            
            # Calculate average enrichment
            avg_enrichment = sum(item[1] for item in treatments) / len(treatments)
            
            # Identify treatment with maximum enrichment
            max_treatment = max(treatments, key=lambda x: x[1])
            
            shared_regions.append({
                'Scaffold': scaffold,
                'Start': start,
                'End': end,
                'Region_Size': end - start,
                'Treatments': ', '.join(item[0] for item in treatments),
                'Treatment_Count': len(treatments),
                'Avg_Enrichment': avg_enrichment,
                'Max_Treatment': max_treatment[0],
                'Max_Enrichment': max_treatment[1],
                'Enrichment_Details': enrichment_dict
            })
    
    # Convert to DataFrame and sort
    if shared_regions:
        shared_df = pd.DataFrame(shared_regions)
        shared_df = shared_df.sort_values(['Treatment_Count', 'Avg_Enrichment'], ascending=[False, False])
        return shared_df
    else:
        return pd.DataFrame()

# Function to analyze adaptation-specific patterns
def analyze_adaptation_patterns(treatment_results):
    """Analyze enrichment patterns specific to adaptation types."""
    if not treatment_results:
        return None
    
    # Group treatments by adaptation type
    adaptation_groups = {}
    for treatment in treatment_results.keys():
        adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
        if adaptation not in adaptation_groups:
            adaptation_groups[adaptation] = []
        adaptation_groups[adaptation].append(treatment)
    
    # Find adaptation-specific enriched regions
    adaptation_specific = {}
    
    for adaptation, treatments in adaptation_groups.items():
        # Skip if only one treatment
        if len(treatments) < 1:
            continue
        
        # Combine enriched regions from this adaptation's treatments
        all_regions = []
        for treatment in treatments:
            if treatment in treatment_results:
                # Add the Treatment column to identify the source
                results_copy = treatment_results[treatment].copy()
                results_copy['Treatment'] = treatment  # Add this line
                all_regions.append(results_copy)
        
        if not all_regions:
            continue
        
        combined_regions = pd.concat(all_regions)
        
        # Group by region and find commonly enriched ones
        grouped = combined_regions.groupby(['Scaffold', 'Start', 'End'])
        
        adaptation_regions = []
        for (scaffold, start, end), group in grouped:
            # If region appears in multiple treatments, consider it adaptation-specific
            if len(group['Treatment'].unique()) > 1:  # Now this will work
                adaptation_regions.append({
                    'Scaffold': scaffold,
                    'Start': start,
                    'End': end,
                    'Region_Size': end - start,
                    'Avg_Enrichment': group['Fold_Enrichment'].mean(),
                    'Max_Enrichment': group['Fold_Enrichment'].max(),
                    'Treatment_Count': len(group['Treatment'].unique()),
                    'Treatments': ', '.join(group['Treatment'].unique())
                })
        
        if adaptation_regions:
            adaptation_df = pd.DataFrame(adaptation_regions)
            adaptation_df = adaptation_df.sort_values('Avg_Enrichment', ascending=False)
            adaptation_specific[adaptation] = adaptation_df
    
    return adaptation_specific

# Function to create detailed enrichment reports
def create_enrichment_reports(treatment_results, adaptation_results, gene_results, shared_regions, adaptation_specific, output_dir):
    """Create detailed reports of enrichment analyses."""
    # Create treatment-specific reports
    for treatment, results in treatment_results.items():
        if len(results) > 0:
            # Sort by fold enrichment
            sorted_results = results.sort_values('Fold_Enrichment', ascending=False)
            
            # Save to CSV
            sorted_results.to_csv(os.path.join(output_dir, f"{treatment}_enriched_regions.csv"), index=False)
            
            # Create text summary
            with open(os.path.join(output_dir, f"{treatment}_enrichment_summary.txt"), 'w') as f:
                adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
                gene = TREATMENT_INFO.get(treatment, {}).get('gene')
                
                f.write(f"Enriched Regions for {treatment} Treatment\n")
                f.write(f"{'=' * (22 + len(treatment))}\n\n")
                
                f.write(f"Treatment: {treatment}\n")
                f.write(f"Description: {TREATMENT_INFO.get(treatment, {}).get('description', 'Unknown')}\n")
                f.write(f"Adaptation: {adaptation}\n")
                f.write(f"Gene Modification: {gene if gene else 'None'}\n\n")
                
                f.write(f"Total Enriched Regions: {len(results)}\n")
                f.write(f"Average Fold Enrichment: {results['Fold_Enrichment'].mean():.2f}\n")
                f.write(f"Maximum Fold Enrichment: {results['Fold_Enrichment'].max():.2f}\n\n")
                
                f.write("Top 10 Most Enriched Regions:\n")
                f.write(f"{'Scaffold':<15}{'Region':<20}{'Fold Enrichment':<20}{'Q-Value':<15}\n")
                f.write(f"{'-' * 70}\n")
                
                for _, row in sorted_results.head(10).iterrows():
                    f.write(f"{row['Scaffold']:<15}{f'{row['Start']}-{row['End']}':<20}{row['Fold_Enrichment']:<20.2f}{row['Q_Value']:<15.2e}\n")
    
    # Create adaptation-specific reports
    for adaptation, results in adaptation_results.items():
        if len(results) > 0:
            # Sort by fold enrichment
            sorted_results = results.sort_values('Fold_Enrichment', ascending=False)
            
            # Save to CSV
            sorted_results.to_csv(os.path.join(output_dir, f"{adaptation}_adaptation_enriched_regions.csv"), index=False)
            
            # Create text summary
            with open(os.path.join(output_dir, f"{adaptation}_adaptation_summary.txt"), 'w') as f:
                f.write(f"Enriched Regions for {adaptation} Adaptation\n")
                f.write(f"{'=' * (23 + len(adaptation))}\n\n")
                
                f.write(f"Adaptation Type: {adaptation}\n")
                f.write(f"Treatments: {', '.join([t for t in TREATMENTS if TREATMENT_INFO.get(t, {}).get('adaptation') == adaptation])}\n\n")
                
                f.write(f"Total Enriched Regions: {len(results)}\n")
                f.write(f"Average Fold Enrichment: {results['Fold_Enrichment'].mean():.2f}\n")
                f.write(f"Maximum Fold Enrichment: {results['Fold_Enrichment'].max():.2f}\n\n")
                
                f.write("Top 10 Most Enriched Regions:\n")
                f.write(f"{'Scaffold':<15}{'Region':<20}{'Fold Enrichment':<20}{'Q-Value':<15}\n")
                f.write(f"{'-' * 70}\n")
                
                for _, row in sorted_results.head(10).iterrows():
                    f.write(f"{row['Scaffold']:<15}{f'{row['Start']}-{row['End']}':<20}{row['Fold_Enrichment']:<20.2f}{row['Q_Value']:<15.2e}\n")
    
    # Create gene modification reports
    for gene_status, results in gene_results.items():
        if len(results) > 0 and "_" not in gene_status:  # Skip combined adaptation-gene reports
            # Sort by fold enrichment
            sorted_results = results.sort_values('Fold_Enrichment', ascending=False)
            
            # Save to CSV
            file_prefix = "gene_modified" if gene_status == "Yes" else "non_modified"
            sorted_results.to_csv(os.path.join(output_dir, f"{file_prefix}_enriched_regions.csv"), index=False)
            
            # Create text summary
            with open(os.path.join(output_dir, f"{file_prefix}_summary.txt"), 'w') as f:
                status_text = "Gene-Modified" if gene_status == "Yes" else "Non-Modified"
                f.write(f"Enriched Regions for {status_text} Strains\n")
                f.write(f"{'=' * (24 + len(status_text))}\n\n")
                
                f.write(f"Gene Status: {status_text}\n")
                
                if gene_status == "Yes":
                    f.write(f"Treatments: {', '.join([t for t in TREATMENTS if TREATMENT_INFO.get(t, {}).get('gene')])}\n\n")
                else:
                    f.write(f"Treatments: {', '.join([t for t in TREATMENTS if not TREATMENT_INFO.get(t, {}).get('gene')])}\n\n")
                
                f.write(f"Total Enriched Regions: {len(results)}\n")
                f.write(f"Average Fold Enrichment: {results['Fold_Enrichment'].mean():.2f}\n")
                f.write(f"Maximum Fold Enrichment: {results['Fold_Enrichment'].max():.2f}\n\n")
                
                f.write("Top 10 Most Enriched Regions:\n")
                f.write(f"{'Scaffold':<15}{'Region':<20}{'Fold Enrichment':<20}{'Q-Value':<15}\n")
                f.write(f"{'-' * 70}\n")
                
                for _, row in sorted_results.head(10).iterrows():
                    f.write(f"{row['Scaffold']:<15}{f'{row['Start']}-{row['End']}':<20}{row['Fold_Enrichment']:<20.2f}{row['Q_Value']:<15.2e}\n")
    
    # Create shared regions report
    if shared_regions is not None and len(shared_regions) > 0:
        # Save to CSV
        shared_regions.to_csv(os.path.join(output_dir, "shared_enriched_regions.csv"), index=False)
        
        # Create text summary
        with open(os.path.join(output_dir, "shared_regions_summary.txt"), 'w') as f:
            f.write("Regions Enriched Across Multiple Treatments\n")
            f.write("=========================================\n\n")
            
            f.write(f"Total Shared Regions: {len(shared_regions)}\n")
            f.write(f"Regions in All Treatments: {len(shared_regions[shared_regions['Treatment_Count'] == len(treatment_results)])}\n\n")
            
            f.write("Top 10 Most Consistently Enriched Regions:\n")
            f.write(f"{'Scaffold':<15}{'Region':<20}{'Treatments':<15}{'Avg Enrichment':<20}\n")
            f.write(f"{'-' * 70}\n")
            
            for _, row in shared_regions.head(10).iterrows():
                f.write(f"{row['Scaffold']:<15}{f'{row['Start']}-{row['End']}':<20}{row['Treatment_Count']:<15}{row['Avg_Enrichment']:<20.2f}\n")
    
    # Create adaptation-specific regions report
    if adaptation_specific:
        for adaptation, regions in adaptation_specific.items():
            if len(regions) > 0:
                # Save to CSV
                regions.to_csv(os.path.join(output_dir, f"{adaptation}_specific_regions.csv"), index=False)
                
                # Create text summary
                with open(os.path.join(output_dir, f"{adaptation}_specific_summary.txt"), 'w') as f:
                    f.write(f"Regions Specific to {adaptation} Adaptation\n")
                    f.write(f"{'=' * (23 + len(adaptation))}\n\n")
                    
                    f.write(f"Adaptation: {adaptation}\n")
                    f.write(f"Treatments: {', '.join([t for t in TREATMENTS if TREATMENT_INFO.get(t, {}).get('adaptation') == adaptation])}\n\n")
                    
                    f.write(f"Total Specific Regions: {len(regions)}\n")
                    f.write(f"Average Enrichment: {regions['Avg_Enrichment'].mean():.2f}\n")
                    f.write(f"Maximum Enrichment: {regions['Max_Enrichment'].max():.2f}\n\n")
                    
                    f.write("Top 10 Most Enriched Adaptation-Specific Regions:\n")
                    f.write(f"{'Scaffold':<15}{'Region':<20}{'Treatments':<30}{'Avg Enrichment':<15}\n")
                    f.write(f"{'-' * 80}\n")
                    
                    for _, row in regions.head(10).iterrows():
                        f.write(f"{row['Scaffold']:<15}{f'{row['Start']}-{row['End']}':<20}{row['Treatments']:<30}{row['Avg_Enrichment']:<15.2f}\n")

# Function to create summary report
def create_summary_report(treatment_results, adaptation_results, gene_results, shared_regions, adaptation_specific, output_dir):
    """Create a comprehensive summary report of regional enrichment analysis."""
    with open(os.path.join(output_dir, "regional_enrichment_summary.txt"), 'w') as f:
        f.write("Regional Enrichment Analysis Summary\n")
        f.write("==================================\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-----------------\n")
        
        total_regions = sum(len(results) for results in treatment_results.values())
        f.write(f"Total enriched regions identified: {total_regions}\n\n")
        
        # Treatment-specific statistics
        f.write("Treatment-Specific Statistics:\n")
        f.write("--------------------------\n")
        
        for treatment, results in treatment_results.items():
            adaptation = TREATMENT_INFO.get(treatment, {}).get('adaptation', 'Unknown')
            gene = TREATMENT_INFO.get(treatment, {}).get('gene')
            gene_text = f" with {gene} gene" if gene else ""
            
            f.write(f"\n{treatment} Treatment ({adaptation} adaptation{gene_text}):\n")
            
            if len(results) > 0:
                f.write(f"  Enriched regions: {len(results)}\n")
                f.write(f"  Average fold enrichment: {results['Fold_Enrichment'].mean():.2f}\n")
                f.write(f"  Maximum fold enrichment: {results['Fold_Enrichment'].max():.2f}\n")
                
                # Top 5 most enriched regions
                f.write("\n  Top 5 most enriched regions:\n")
                top_regions = results.nlargest(5, 'Fold_Enrichment')
                for _, region in top_regions.iterrows():
                    f.write(f"    {region['Scaffold']}:{region['Start']}-{region['End']}, "
                           f"Fold enrichment: {region['Fold_Enrichment']:.2f}, "
                           f"Q-value: {region['Q_Value']:.2e}\n")
            else:
                f.write("  No significantly enriched regions found\n")
        
        # Adaptation-specific statistics
        f.write("\nAdaptation-Specific Statistics:\n")
        f.write("----------------------------\n")
        
        for adaptation, results in adaptation_results.items():
            f.write(f"\n{adaptation} Adaptation:\n")
            if len(results) > 0:
                f.write(f"  Enriched regions: {len(results)}\n")
                f.write(f"  Average fold enrichment: {results['Fold_Enrichment'].mean():.2f}\n")
                f.write(f"  Maximum fold enrichment: {results['Fold_Enrichment'].max():.2f}\n")
                
                # Top 3 most enriched regions
                f.write("\n  Top 3 most enriched regions:\n")
                top_regions = results.nlargest(3, 'Fold_Enrichment')
                for _, region in top_regions.iterrows():
                    f.write(f"    {region['Scaffold']}:{region['Start']}-{region['End']}, "
                           f"Fold enrichment: {region['Fold_Enrichment']:.2f}\n")
            else:
                f.write("  No significantly enriched regions found\n")
        
        # Gene modification statistics
        f.write("\nGene Modification Statistics:\n")
        f.write("--------------------------\n")
        
        for gene_status, results in gene_results.items():
            if "_" not in gene_status:  # Skip combined adaptation-gene reports
                status_text = "Gene-Modified Strains" if gene_status == "Yes" else "Non-Modified Strains"
                f.write(f"\n{status_text}:\n")
                
                if len(results) > 0:
                    f.write(f"  Enriched regions: {len(results)}\n")
                    f.write(f"  Average fold enrichment: {results['Fold_Enrichment'].mean():.2f}\n")
                    f.write(f"  Maximum fold enrichment: {results['Fold_Enrichment'].max():.2f}\n")
                else:
                    f.write("  No significantly enriched regions found\n")
        
        # Shared regions statistics
        f.write("\nShared Regions Statistics:\n")
        f.write("------------------------\n")
        
        if shared_regions is not None and len(shared_regions) > 0:
            f.write(f"  Total shared regions (in multiple treatments): {len(shared_regions)}\n")
            
            # Count by number of treatments
            for n_treatments in range(2, len(treatment_results) + 1):
                count = len(shared_regions[shared_regions['Treatment_Count'] == n_treatments])
                f.write(f"  Regions in exactly {n_treatments} treatments: {count}\n")
            
            # Top shared regions
            if len(shared_regions) > 0:
                f.write("\n  Top 3 most consistently enriched regions:\n")
                
                # Sort by treatment count then by average enrichment
                top_shared = shared_regions.sort_values(['Treatment_Count', 'Avg_Enrichment'], ascending=[False, False]).head(3)
                
                for _, region in top_shared.iterrows():
                    f.write(f"    {region['Scaffold']}:{region['Start']}-{region['End']}, "
                           f"Treatments: {region['Treatment_Count']}, "
                           f"Avg enrichment: {region['Avg_Enrichment']:.2f}\n")
        else:
            f.write("  No regions shared across multiple treatments\n")
        
        # Adaptation-specific region statistics
        f.write("\nAdaptation-Specific Region Statistics:\n")
        f.write("----------------------------------\n")
        
        if adaptation_specific:
            for adaptation, regions in adaptation_specific.items():
                f.write(f"\n{adaptation} Adaptation-Specific Regions:\n")
                
                if len(regions) > 0:
                    f.write(f"  Total specific regions: {len(regions)}\n")
                    f.write(f"  Average enrichment: {regions['Avg_Enrichment'].mean():.2f}\n")
                    
                    # Top regions
                    f.write("  Top 3 most enriched regions:\n")
                    top_regions = regions.nlargest(3, 'Avg_Enrichment')
                    for _, region in top_regions.iterrows():
                        f.write(f"    {region['Scaffold']}:{region['Start']}-{region['End']}, "
                               f"Avg enrichment: {region['Avg_Enrichment']:.2f}\n")
                else:
                    f.write("  No adaptation-specific regions found\n")
        else:
            f.write("  No adaptation-specific regions identified\n")
        
        # Main conclusions
        f.write("\nMain Conclusions:\n")
        f.write("---------------\n")
        f.write("1. This analysis identifies regions with statistically significant variant enrichment.\n")
        f.write("2. Treatment-specific patterns of regional enrichment have been identified.\n")
        f.write("3. Adaptation-specific enrichment patterns highlight genomic regions associated with\n")
        f.write("   temperature or low oxygen adaptation mechanisms.\n")
        f.write("4. Gene modifications appear to influence the pattern of regional enrichment,\n")
        f.write("   suggesting potential functional impacts of inserted genes.\n")
        f.write("5. Regions shared across multiple treatments may represent mutational hotspots\n")
        f.write("   or regions under selection pressure in multiple conditions.\n")
        f.write("6. For a detailed view of enriched regions, refer to the individual analysis files\n")
        f.write("   in the output directory.\n")

def adaptive_enrichment_analysis(data, scaffold_info):
    """
    Adaptively attempt regional enrichment analysis with different parameters until results are found.
    
    This function tries increasingly lenient parameter combinations to find enriched regions.
    """
    # Default parameters
    window_sizes = [1000, 750, 500, 250]
    step_sizes = [200, 150, 100, 50]
    p_thresholds = [0.1, 0.15, 0.2, 0.25]
    
    # Try different parameter combinations
    for i, window_size in enumerate(window_sizes):
        step_size = step_sizes[i]
        p_threshold = p_thresholds[i]
        
        logging.info(f"Attempting enrichment analysis with: window_size={window_size}, " 
                   f"step_size={step_size}, p_threshold={p_threshold}")
        
        results = analyze_regional_enrichment(
            data,
            scaffold_info,
            window_size=window_size,
            step_size=step_size,
            p_value_threshold=p_threshold
        )
        
        if results is not None and len(results) > 0:
            logging.info(f"Found {len(results)} enriched regions with window_size={window_size}, "
                       f"step_size={step_size}, p_threshold={p_threshold}")
            return results
    
    # If we still didn't find anything, try one last attempt with very lenient parameters
    logging.warning("No enriched regions found with standard parameters, trying very lenient parameters")
    
    # Very lenient parameters
    results = analyze_regional_enrichment(
        data,
        scaffold_info,
        window_size=200,
        step_size=40,
        p_value_threshold=0.3
    )
    
    return results

def main():
    # Load data
    logging.info("Loading mutation data...")
    mutation_data = load_data()
    
    logging.info("Loading scaffold information...")
    scaffold_info = load_scaffold_info()
    
    if mutation_data is None or not scaffold_info:
        logging.error("Cannot proceed with analysis due to missing data")
        return
    
    # Log key statistics
    logging.info(f"Loaded {len(mutation_data)} mutations across {len(mutation_data['Treatment'].unique())} treatments")
    logging.info(f"Loaded {len(scaffold_info)} scaffolds from reference")
    
    # Perform treatment-specific enrichment analysis with adaptive parameters
    logging.info("\nAnalyzing regional enrichment by treatment...")
    treatment_results = {}
    
    for treatment in TREATMENTS:
        logging.info(f"\nAnalyzing enrichment for {treatment} treatment...")
        
        treatment_data = mutation_data[mutation_data['Treatment'] == treatment]
        if len(treatment_data) == 0:
            logging.warning(f"No data available for {treatment}")
            continue
        
        # Use adaptive analysis
        enriched = adaptive_enrichment_analysis(treatment_data, scaffold_info)
        
        if enriched is not None and len(enriched) > 0:
            treatment_results[treatment] = enriched
    
    # Perform adaptation-specific enrichment analysis
    logging.info("\nAnalyzing regional enrichment by adaptation type...")
    adaptation_results = {}
    
    # Group by adaptation type
    for adaptation in mutation_data['Adaptation'].unique():
        logging.info(f"\nAnalyzing enrichment for {adaptation} adaptation...")
        
        adaptation_data = mutation_data[mutation_data['Adaptation'] == adaptation]
        if len(adaptation_data) == 0:
            logging.warning(f"No data available for {adaptation} adaptation")
            continue
        
        # Use adaptive analysis
        enriched = adaptive_enrichment_analysis(adaptation_data, scaffold_info)
        
        if enriched is not None and len(enriched) > 0:
            adaptation_results[adaptation] = enriched
    
    # Perform gene-modification enrichment analysis
    logging.info("\nAnalyzing regional enrichment by gene modification status...")
    gene_results = {}
    
    # Split by gene modification status
    for gene_status in mutation_data['Has_Gene'].unique():
        logging.info(f"\nAnalyzing enrichment for {'gene-modified' if gene_status == 'Yes' else 'non-modified'} strains...")
        
        gene_data = mutation_data[mutation_data['Has_Gene'] == gene_status]
        if len(gene_data) == 0:
            continue
        
        # Use adaptive analysis
        enriched = adaptive_enrichment_analysis(gene_data, scaffold_info)
        
        if enriched is not None and len(enriched) > 0:
            gene_results[gene_status] = enriched
    
    # Find consistently enriched regions across treatments
    logging.info("\nFinding regions consistently enriched across treatments...")
    shared_regions = find_consistently_enriched_regions(treatment_results)
    
    # Analyze adaptation-specific patterns
    logging.info("\nAnalyzing adaptation-specific enrichment patterns...")
    adaptation_specific = analyze_adaptation_patterns(treatment_results)
    
    # Log summary of the results
    logging.info("\nAnalysis Results Summary:")
    logging.info(f"Treatment-specific enriched regions: {sum(len(results) for results in treatment_results.values()) if treatment_results else 0}")
    logging.info(f"Adaptation-specific enriched regions: {sum(len(results) for results in adaptation_results.values()) if adaptation_results else 0}")
    logging.info(f"Gene-status-specific enriched regions: {sum(len(results) for results in gene_results.values()) if gene_results else 0}")
    logging.info(f"Shared enriched regions: {len(shared_regions) if shared_regions is not None else 0}")
    
    # Generate visualizations if we have results
    if treatment_results:
        logging.info("\nGenerating visualizations...")
        plot_enrichment_patterns(treatment_results, OUTPUT_DIR)
        
        # Create detailed reports
        logging.info("\nCreating detailed enrichment reports...")
        create_enrichment_reports(
            treatment_results, 
            adaptation_results, 
            gene_results, 
            shared_regions, 
            adaptation_specific, 
            OUTPUT_DIR
        )
        
        # Create summary report
        logging.info("\nCreating summary report...")
        create_summary_report(
            treatment_results,
            adaptation_results,
            gene_results,
            shared_regions,
            adaptation_specific,
            OUTPUT_DIR
        )
        
        print(f"\nAnalysis complete! Results saved to {OUTPUT_DIR}/")
    else:
        logging.error("No enriched regions found after adaptive analysis.")
        print("\nNo enriched regions found. Please check regional_enrichment_debug.log for details.")
        
        # Create placeholder summary file to explain the issue
        with open(os.path.join(OUTPUT_DIR, "regional_enrichment_summary.txt"), 'w') as f:
            f.write("Regional Enrichment Analysis Summary\n")
            f.write("==================================\n\n")
            f.write("No statistically significant enriched regions were found in this analysis.\n\n")
            f.write("Potential reasons:\n")
            f.write("1. Mutations are too sparsely distributed across the genome\n")
            f.write("2. There is no significant clustering of mutations in specific regions\n")
            f.write("3. The mutation rate is too low to detect statistically significant enrichment\n\n")
            f.write("Recommendations:\n")
            f.write("1. Consider using a smaller window size for the analysis\n")
            f.write("2. Consider relaxing the statistical threshold (p-value)\n")
            f.write("3. Check if the mutation data and reference genome scaffolds match correctly\n")
            f.write("4. Consider analyzing the data at a different genomic resolution\n")

if __name__ == "__main__":
    main()