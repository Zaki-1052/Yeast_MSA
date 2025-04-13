#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency
import statsmodels.stats.multitest as mt
import subprocess

def load_vcf_counts(vcf_file):
    """Extract variant counts from a VCF file."""
    try:
        # Count total variants
        cmd = f"bcftools view -H {vcf_file} | wc -l"
        variant_count = int(subprocess.check_output(cmd, shell=True))
        return variant_count
    except Exception as e:
        print(f"Error processing {vcf_file}: {e}")
        return 0

def analyze_treatment_control_differences():
    """Analyze statistical significance of treatment vs control differences."""
    
    # Define treatments and their controls
    treatments = {
        'WT': {'treatment': 'WT-37-55', 'control': 'WT-CTRL'},
        'STC': {'treatment': 'STC-55', 'control': 'STC-CTRL'},
        'CAS': {'treatment': 'CAS-55', 'control': 'CAS-CTRL'},
        'WTA': {'treatment': 'WTA-55', 'control': 'WT-CTRL'}  # Note: WTA uses WT control
    }
    
    results = []
    
    # Process each treatment
    for treatment_name, info in treatments.items():
        # Get counts from VCF files
        treatment_vcf = f"results/merged/analysis/{treatment_name}_specific.vcf.gz"
        control_vcf = f"results/vcf/individual/{info['control']}.vcf.gz"
        
        treatment_count = load_vcf_counts(treatment_vcf)
        control_count = load_vcf_counts(control_vcf)
        
        # Create contingency table
        # Using genome size as background (total possible positions)
        genome_size = 12000000  # Approximate S. cerevisiae genome size
        
        contingency_table = np.array([
            [treatment_count, genome_size - treatment_count],
            [control_count, genome_size - control_count]
        ])
        
        # Perform Fisher's exact test
        odds_ratio, p_value = fisher_exact(contingency_table)
        
        # Store results
        results.append({
            'Treatment': treatment_name,
            'Treatment_Variants': treatment_count,
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
    os.makedirs('treatment_control_analysis', exist_ok=True)
    
    # Save results
    results_df.to_csv('treatment_control_analysis/treatment_vs_control_statistics.csv', index=False)
    
    # Create a detailed report
    with open('treatment_control_analysis/statistical_analysis_report.txt', 'w') as f:
        f.write("Treatment vs Control Statistical Analysis\n")
        f.write("=====================================\n\n")
        
        for _, row in results_df.iterrows():
            f.write(f"{row['Treatment']} Treatment Analysis:\n")
            f.write("-" * (len(row['Treatment']) + 20) + "\n")
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
    print(results[['Treatment', 'Treatment_Variants', 'Control_Variants', 'Fold_Change', 'Q_Value']])