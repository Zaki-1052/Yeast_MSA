#!/usr/bin/env python3

import os
import pandas as pd
from datetime import datetime

def main():
    # Define output file
    output_file = "combined_analysis_results.txt"
    
    # List of summary files to concatenate
    summary_files = [
        "analysis/genomic_context_results/genomic_context_summary.txt",
        "analysis/mutation_spectrum_results/statistical_test_results.txt",
        "analysis/population_structure_results/population_structure_summary.txt",
        "analysis/mutational_signatures_results/mutational_signatures_summary.txt",
        "analysis/scaffold_distribution_results/scaffold_distribution_summary.txt",
        "analysis/regional_enrichment_results/regional_enrichment_summary.txt",
        "analysis/statistical_pattern_results/statistical_analysis_summary.txt",
        "analysis/treatment_control_analysis/statistical_analysis_report.txt"
    ]
    
    # CSV files to convert and include
    csv_files = [
        "analysis/mutation_spectrum_results/mutation_spectrum_summary.csv",
        "analysis/scaffold_distribution_results/scaffold_distribution_summary.csv"
    ]
    
    # Get current directory
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    
    # Create header with timestamp
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
    
    with open(os.path.join(base_dir, output_file), 'w') as outfile:
        # Write header
        outfile.write(f"# YEAST MSA COMBINED ANALYSIS RESULTS\n")
        outfile.write(f"# Generated: {timestamp}\n")
        outfile.write(f"# ==========================================\n\n")
        
        # Process text files
        for file_path in summary_files:
            full_path = os.path.join(base_dir, file_path)
            if os.path.exists(full_path):
                section_name = os.path.basename(os.path.dirname(file_path))
                outfile.write(f"\n\n{'#' * 80}\n")
                outfile.write(f"# {section_name.upper().replace('_', ' ')}\n")
                outfile.write(f"{'#' * 80}\n\n")
                
                with open(full_path, 'r') as infile:
                    content = infile.read()
                    outfile.write(content)
                    if not content.endswith('\n'):
                        outfile.write('\n')
            else:
                print(f"Warning: File not found - {full_path}")
        
        # Process CSV files
        for csv_path in csv_files:
            full_path = os.path.join(base_dir, csv_path)
            if os.path.exists(full_path):
                section_name = os.path.basename(os.path.dirname(csv_path))
                outfile.write(f"\n\n{'#' * 80}\n")
                outfile.write(f"# {section_name.upper().replace('_', ' ')} - {os.path.basename(csv_path)}\n")
                outfile.write(f"{'#' * 80}\n\n")
                
                try:
                    df = pd.read_csv(full_path)
                    outfile.write(df.to_string(index=False))
                    outfile.write('\n')
                except Exception as e:
                    outfile.write(f"Error processing CSV file: {e}\n")
                    # Fallback to simple read if pandas fails
                    with open(full_path, 'r') as infile:
                        outfile.write(infile.read())
            else:
                print(f"Warning: File not found - {full_path}")
    
    print(f"Combined results saved to {os.path.join(base_dir, output_file)}")

if __name__ == "__main__":
    main()