import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Output directory
OUTDIR = 'results/filtered_scaffold_variants/impact/'
os.makedirs(OUTDIR, exist_ok=True)

# File paths
EXPANDED_TSV = 'results/filtered_scaffold_variants/treatment_specific_variant_annotation_expanded.tsv'
SUMMARY_TSV = os.path.join(OUTDIR, 'variant_proximity_impact_summary.tsv')
HEATMAP_PNG = os.path.join(OUTDIR, 'variant_count_heatmap.png')
DIST_HEATMAP_PNG = os.path.join(OUTDIR, 'distance_distribution_heatmap.png')
SUMMARY_HTML = os.path.join(OUTDIR, 'variant_proximity_impact_summary.html')

# Read expanded annotation table
df = pd.read_csv(EXPANDED_TSV, sep='\t')

# Use ERG name and scaffold for grouping
erg_name_map = df.set_index('Nearest_ERG_Gene')['Nearest_ERG_Gene_ERG'].to_dict()
erg_scaffold_map = df.set_index('Nearest_ERG_Gene')['Scaffold'].to_dict()
df['ERG_Name'] = df['Nearest_ERG_Gene'].map(erg_name_map)
df['ERG_Scaffold'] = df['Nearest_ERG_Gene'].map(erg_scaffold_map)

# Group and aggregate: ERG name, scaffold, treatment, min distance, count, position, impact
summary = df.groupby([
    'ERG_Name', 'ERG_Scaffold', 'Treatment', 'ERG_Position_Relative', 'Impact'
]).agg(
    Variant_Count=('Sample', 'count'),
    Min_Distance_to_ERG=('Distance_to_ERG', lambda x: float(pd.to_numeric(x, errors='coerce').min()))
).reset_index()

# Reorder columns
summary = summary[['ERG_Name', 'ERG_Scaffold', 'Treatment', 'Min_Distance_to_ERG', 'Variant_Count', 'ERG_Position_Relative', 'Impact']]

# Sort by increasing Min_Distance_to_ERG
summary = summary.sort_values('Min_Distance_to_ERG', ascending=True)

# Save summary table
summary.to_csv(SUMMARY_TSV, sep='\t', index=False)

# Heatmap: Variant count for each treatment for each ERG gene (regardless of distance)
count_heatmap_data = df.groupby(['ERG_Name', 'Treatment']).size().unstack(fill_value=0)
plt.figure(figsize=(10, 6))
sns.heatmap(count_heatmap_data, annot=True, fmt='d', cmap='YlOrRd')
plt.title('Variant Count for Each Treatment and ERG Gene')
plt.ylabel('ERG Gene')
plt.xlabel('Treatment')
plt.tight_layout()
plt.savefig(HEATMAP_PNG)
plt.close()

# Distance distribution: For each ERG gene and treatment, show distribution of distances (boxplot)
plt.figure(figsize=(12, 7))
sns.boxplot(
    data=df,
    x='ERG_Name',
    y=pd.to_numeric(df['Distance_to_ERG'], errors='coerce'),
    hue='Treatment',
    showfliers=True
)
plt.title('Distribution of Distances to ERG Genes by Treatment')
plt.ylabel('Distance to ERG Gene (bp)')
plt.xlabel('ERG Gene')
plt.legend(title='Treatment', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(DIST_HEATMAP_PNG)
plt.close()

# Generate HTML report
html = f"""
<!DOCTYPE html>
<html lang='en'>
<head>
    <meta charset='UTF-8'>
    <title>Variant Proximity & Impact Summary</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #333; }}
        table {{ border-collapse: collapse; width: 100%; margin-bottom: 30px; }}
        th, td {{ border: 1px solid #ccc; padding: 6px 10px; text-align: left; }}
        th {{ background: #f2f2f2; }}
        img {{ max-width: 800px; display: block; margin: 20px 0; }}
    </style>
</head>
<body>
    <h1>Variant Proximity & Impact Summary</h1>
    <h2>Summary Table (Sorted by Increasing Distance)</h2>
    <table>
        <tr>{''.join(f'<th>{col}</th>' for col in summary.columns)}</tr>
"""
for _, row in summary.iterrows():
    html += '<tr>' + ''.join(f'<td>{row[col]}</td>' for col in summary.columns) + '</tr>'
html += f"""
    </table>
    <h2>Heatmap: Variant Count for Each Treatment and ERG Gene</h2>
    <img src='variant_count_heatmap.png' alt='Variant Count Heatmap'>
    <h2>Boxplot: Distribution of Distances to ERG Genes by Treatment</h2>
    <img src='distance_distribution_heatmap.png' alt='Distance Distribution Boxplot'>
</body>
</html>
"""

with open(SUMMARY_HTML, 'w') as f:
    f.write(html)

print(f"Summary table written to {SUMMARY_TSV}")
print(f"Variant count heatmap written to {HEATMAP_PNG}")
print(f"Distance distribution boxplot written to {DIST_HEATMAP_PNG}")
print(f"HTML report written to {SUMMARY_HTML}") 