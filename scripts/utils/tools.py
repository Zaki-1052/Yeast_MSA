#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/utils/tools.py

"""
Utility functions for the Yeast MSA project
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def ensure_dir(directory):
    """Create directory if it doesn't exist"""
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
            print(f"Created directory: {directory}")
        except Exception as e:
            print(f"Error creating directory {directory}: {e}")
            sys.exit(1)
    return directory

def load_tsv(file_path, sep='\t'):
    """Load a TSV file into a pandas DataFrame with error handling"""
    if not os.path.exists(file_path):
        print(f"ERROR: File not found: {file_path}")
        return None
    
    try:
        df = pd.read_csv(file_path, sep=sep)
        return df
    except Exception as e:
        print(f"ERROR: Failed to load file {file_path}: {e}")
        return None

def save_tsv(df, file_path, sep='\t', index=False):
    """Save a pandas DataFrame to a TSV file with error handling"""
    try:
        ensure_dir(os.path.dirname(file_path))
        df.to_csv(file_path, sep=sep, index=index)
        print(f"Saved file: {file_path}")
        return True
    except Exception as e:
        print(f"ERROR: Failed to save file {file_path}: {e}")
        return False

def create_dir_structure(base_dir, subdirs):
    """Create a directory structure"""
    for subdir in subdirs:
        ensure_dir(os.path.join(base_dir, subdir))
    
    print(f"Created directory structure in {base_dir}")

def setup_plotting_style():
    """Set up a consistent plotting style"""
    # Use a clean, modern style
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Set default figure size
    plt.rcParams['figure.figsize'] = [10, 6]
    
    # Set fonts
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
    
    # Set font sizes
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    
    # Set color palette
    sns.set_palette("colorblind")

def save_plot(fig, output_file, dpi=300):
    """Save a matplotlib figure with error handling"""
    try:
        ensure_dir(os.path.dirname(output_file))
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Saved figure: {output_file}")
        return True
    except Exception as e:
        print(f"ERROR: Failed to save figure {output_file}: {e}")
        return False

def parse_fasta(file_path):
    """Parse a FASTA file and return a dictionary of sequences"""
    if not os.path.exists(file_path):
        print(f"ERROR: FASTA file not found: {file_path}")
        return None
    
    sequences = {}
    
    try:
        with open(file_path, 'r') as f:
            header = None
            seq = []
            
            for line in f:
                line = line.strip()
                
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Save the previous sequence if it exists
                    if header is not None:
                        sequences[header] = ''.join(seq)
                    
                    # Start a new sequence
                    header = line[1:].split()[0]  # Get the first word after '>'
                    seq = []
                else:
                    seq.append(line)
            
            # Save the last sequence
            if header is not None:
                sequences[header] = ''.join(seq)
        
        return sequences
    except Exception as e:
        print(f"ERROR: Failed to parse FASTA file {file_path}: {e}")
        return None

def format_p_value(p_value):
    """Format p-value for display"""
    if p_value < 0.001:
        return f"p < 0.001"
    elif p_value < 0.01:
        return f"p < 0.01"
    elif p_value < 0.05:
        return f"p < 0.05"
    else:
        return f"p = {p_value:.3f}"

def calculate_enrichment(observed, expected, pseudocount=1):
    """Calculate enrichment factor with pseudocount"""
    return (observed + pseudocount) / (expected + pseudocount)

def create_distance_bins(max_distance, num_bins=10, log_scale=False):
    """Create distance bins for analysis"""
    if log_scale:
        # Create logarithmic bins
        return np.logspace(0, np.log10(max_distance), num_bins)
    else:
        # Create linear bins
        return np.linspace(0, max_distance, num_bins)

def calculate_gc_content(sequence):
    """Calculate GC content of a DNA sequence"""
    if not sequence:
        return 0
    
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)

def check_dependencies():
    """Check if required Python packages are installed"""
    required_packages = {
        'pandas': 'DataFrame',
        'numpy': 'array',
        'matplotlib': 'pyplot',
        'seaborn': 'heatmap',
        'scipy': 'stats'
    }
    
    missing_packages = []
    
    for package, attribute in required_packages.items():
        try:
            module = __import__(package)
            if not hasattr(module, attribute):
                missing_packages.append(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        print("ERROR: The following required packages are missing:")
        for package in missing_packages:
            print(f"  - {package}")
        print("\nPlease install them using:")
        print(f"pip install {' '.join(missing_packages)}")
        return False
    
    return True