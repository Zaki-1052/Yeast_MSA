#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency
import statsmodels.stats.multitest as mt
import subprocess

def find_vcf_file(base_name, possible_locations):
    """Find a VCF file by trying multiple possible locations and formats."""
    for location in possible_locations:
        if os.path.exists(location.format(base_name)):
            return location.format(base_name)
    return None

def load_vcf_counts(vcf_file):
    """Extract variant counts from a VCF file."""
    if not vcf_file or not os.path.exists(vcf_file):
        print(f"  Error: VCF file not found: {vcf_file}")
        return 0
    
    try:
        # Count total variants
        cmd = f"bcftools view -H {vcf_file} | wc -l"
        variant_count = int(subprocess.check_output(cmd, shell=True))
        return variant_count
    except Exception as e:
        print(f"  Error processing {vcf_file}: {e}")
        return 0

def analyze_treatment_control_differences():
    """Analyze statistical significance of treatment vs control differences."""
    
    # Define treatments and their controls
    # Updated to reflect the correct biological grouping
    treatments = {
        # Primary comparisons with WT-CTRL as baseline
        'WT-37': {'treatment': 'WT-37-55', 'control': 'WT-CTRL', 'description': 'Temperature-adapted wild type'},
        'WTA': {'treatment': 'WTA-55', 'control': 'WT-CTRL', 'description': 'Low oxygen-adapted wild type'},
        'STC': {'treatment': 'STC-55', 'control': 'WT-CTRL', 'description': 'STC gene with low oxygen adaptation'},
        'CAS': {'treatment': 'CAS-55', 'control': 'WT-CTRL', 'description': 'CAS gene with temperature adaptation'},
        
        # Original control comparisons (preserved for completeness)
        'STC-vs-STCCTRL': {'treatment': 'STC-55', 'control': 'STC-CTRL', 'description': 'STC vs STC control'},
        'CAS-vs-CASCTRL': {'treatment': 'CAS-55', 'control': 'CAS-CTRL', 'description': 'CAS vs CAS control'}
    }
    
    # Define possible treatment VCF locations to check
    treatment_locations = [
        "results/merged/analysis/{}/highconf.vcf.gz",
        "results/merged/analysis/{}_highconf.vcf.gz",
        "results/merged/analysis/{}/specific.vcf.gz",
        "results/merged/analysis/{}_specific.vcf.gz"
    ]
    
    # Define possible control VCF locations to check
    control_locations = [
        # Try individual directory with different extensions
        "results/vcf/individual/{}.vcf.gz",
        "results/vcf/individual/{}.vcf",
        "results/vcf/individual/{}.norm.vcf",
        # Try merged filtered directory
        "results/vcf/merged/filtered/{}.filtered.vcf.gz",
        # Try merged fixed directory
        "results/vcf/merged/fixed/{}.fixed.vcf.gz"
    ]
    
    results = []
    
    # Process each treatment
    for treatment_name, info in treatments.items():
        print(f"Treatment: {treatment_name}")
        
        # Find treatment VCF
        if treatment_name in ['WT-37', 'WTA', 'STC', 'CAS']:
            # For main treatment groups
            treatment_vcf = find_vcf_file(treatment_name, treatment_locations)
        else:
            # For the original control comparisons
            base_treatment = treatment_name.split('-vs-')[0]
            treatment_vcf = find_vcf_file(base_treatment, treatment_locations)
        
        # Find control VCF
        control_vcf = find_vcf_file(info['control'], control_locations)
        
        print(f"  Treatment VCF: {treatment_vcf}")
        print(f"  Control VCF: {control_vcf}")
        
        # Get variant counts
        treatment_count = load_vcf_counts(treatment_vcf)
        control_count = load_vcf_counts(control_vcf)
        
        print(f"  Treatment count: {treatment_count}")
        print(f"  Control count: {control_count}")
        
        # Create contingency table
        # Using genome size as background (total possible positions)
        genome_size = 12000000  # Approximate S. cerevisiae genome size
        
        contingency_table = np.array([
            [treatment_count, genome_size - treatment_count],
            [control_count, genome_size - control_count]
        ])
        
        # Perform Fisher's exact test (stable and doesn't require chi-square assumptions)
        odds_ratio, p_value = fisher_exact(contingency_table)
        
        # Make sure the p-value is not too extreme
        #p_value = max(p_value, 1e-10)  # Limit extremely small p-values
        
        # Store results
        results.append({
            'Treatment': treatment_name,
            'Description': info['description'],
            'Treatment_Variants': treatment_count,
            'Control': info['control'],
            'Control_Variants': control_count,
            'Odds_Ratio': odds_ratio,
            'P_Value': p_value
        })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Perform multiple testing correction
    results_df['Q_Value'] = mt.multipletests(results_df['P_Value'], method='fdr_bh')[1]
    
    # Add fold change
    results_df['Fold_Change'] = results_df['Treatment_Variants'] / results_df['Control_Variants'].replace(0, 1)
    
    # Sort by significance
    results_df = results_df.sort_values('P_Value')
    
    # Create output directory if it doesn't exist
    os.makedirs('analysis/treatment_control_analysis', exist_ok=True)
    
    # Save results
    results_df.to_csv('analysis/treatment_control_analysis/treatment_vs_control_statistics.csv', index=False)
    
    # Create a detailed report
    with open('analysis/treatment_control_analysis/statistical_analysis_report.txt', 'w') as f:
        f.write("Treatment vs Control Statistical Analysis\n")
        f.write("=====================================\n\n")
        
        for _, row in results_df.iterrows():
            f.write(f"{row['Treatment']} Treatment Analysis:\n")
            f.write("-" * (len(row['Treatment']) + 20) + "\n")
            f.write(f"Description: {row['Description']}\n")
            f.write(f"Control: {row['Control']}\n")
            f.write(f"Treatment variants: {row['Treatment_Variants']}\n")
            f.write(f"Control variants: {row['Control_Variants']}\n")
            f.write(f"Fold change: {row['Fold_Change']:.2f}\n")
            f.write(f"Odds ratio: {row['Odds_Ratio']:.2f}\n")
            f.write(f"P-value: {row['P_Value']:.2e}\n")
            f.write(f"Q-value (FDR-corrected): {row['Q_Value']:.2e}\n")
            f.write(f"Statistical significance: {'***' if row['Q_Value'] < 0.001 else '**' if row['Q_Value'] < 0.01 else '*' if row['Q_Value'] < 0.05 else 'ns'}\n\n")
    
    print("Analysis complete! Results saved to treatment_control_analysis/")
    return results_df

if __name__ == "__main__":
    results = analyze_treatment_control_differences()
    print("\nSummary of Results:")
    print(results[['Treatment', 'Control', 'Treatment_Variants', 'Control_Variants', 'Fold_Change', 'Q_Value']])