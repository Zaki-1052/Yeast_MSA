#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch

def create_treatment_control_visualizations(results_df):
    """Create intuitive visualizations of treatment vs control differences."""
    
    # Debug Point 1: Check input data
    print("\nChecking input data:")
    print(results_df.head())
    print("\nChecking for NaN values:")
    print(results_df.isna().sum())
    
    # Set style parameters
    plt.rcParams['figure.figsize'] = [10, 6]
    plt.rcParams['font.size'] = 12
    sns.set_palette("husl")
    
    # Debug Point 2: Print value ranges
    print("\nValue ranges:")
    for col in results_df.columns:
        if pd.api.types.is_numeric_dtype(results_df[col]):
            print(f"{col}: min={results_df[col].min()}, max={results_df[col].max()}")
    
    # Group treatments by control for better visualization
    results_df['Control_Group'] = results_df['Control'].apply(lambda x: x.split('-')[0])
    
    # Define a better color scheme for different treatment groups
    color_map = {
        'WT-37': '#1b9e77',    # Temperature-adapted wild type
        'WTA': '#d95f02',      # Low oxygen-adapted wild type
        'STC': '#7570b3',      # STC gene with low oxygen adaptation
        'CAS': '#e7298a',      # CAS gene with temperature adaptation
        'STC-vs-STCCTRL': '#66a61e',  # STC with original control
        'CAS-vs-CASCTRL': '#e6ab02'   # CAS with original control
    }
    
    # 1. Bar Plot with Control vs Treatment
    plt.figure(figsize=(14, 8))
    
    # Prepare data
    treatments = results_df['Treatment']
    x = np.arange(len(treatments))
    width = 0.35
    
    # Debug Point 3: Print bar plot data
    print("\nBar plot data:")
    print("Control values:", results_df['Control_Variants'].tolist())
    print("Treatment values:", results_df['Treatment_Variants'].tolist())
    
    # Create bars with error checking
    control_vals = np.nan_to_num(results_df['Control_Variants'].values, 0)
    treatment_vals = np.nan_to_num(results_df['Treatment_Variants'].values, 0)
    
    plt.bar(x - width/2, control_vals, width, label='Control', color='lightgray')
    treatment_bars = plt.bar(x + width/2, treatment_vals, width, label='Treatment')
    
    # Apply custom colors to treatment bars
    for i, bar in enumerate(treatment_bars):
        treatment = treatments.iloc[i]
        bar.set_color(color_map.get(treatment, '#2ecc71'))
    
    # Customize plot
    plt.xlabel('Treatment Group')
    plt.ylabel('Number of Variants')
    plt.title('Comparison of Variants: Treatment vs Control', pad=20)
    
    # Custom x-tick labels that include control names
    x_labels = [f"{t}\n(vs {c})" for t, c in zip(results_df['Treatment'], results_df['Control'])]
    plt.xticks(x, x_labels, rotation=45, ha='right')
    
    plt.legend()
    
    # Add treatment effect arrows
    for i, (control, treatment) in enumerate(zip(control_vals, treatment_vals)):
        plt.arrow(x[i], control, 0, treatment-control,
                head_width=0.1, head_length=min(20, max(5, (treatment-control)/10)),
                fc='#e74c3c', ec='#e74c3c',
                alpha=0.5)
    
    # Add significance stars with error checking
    for i, row in enumerate(results_df.itertuples()):
        if pd.notna(row.Q_Value) and pd.notna(row.Treatment_Variants) and pd.notna(row.Control_Variants):
            stars = '***' if row.Q_Value < 0.001 else '**' if row.Q_Value < 0.01 else '*' if row.Q_Value < 0.05 else 'ns'
            max_height = max(row.Treatment_Variants, row.Control_Variants)
            plt.text(i, max_height * 1.1, stars, ha='center', va='bottom')
    
    # Add fold change annotations with error checking
    for i, row in enumerate(results_df.itertuples()):
        if pd.notna(row.Fold_Change) and pd.notna(row.Treatment_Variants):
            plt.text(i, row.Treatment_Variants * 1.02,
                    f'FC: {row.Fold_Change:.1f}×', 
                    ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('analysis/treatment_control_analysis/variant_comparison_barplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Grouped Bar Plot with Primary vs Secondary Comparisons
    primary_comparisons = results_df[results_df['Control'] == 'WT-CTRL']
    secondary_comparisons = results_df[results_df['Control'] != 'WT-CTRL']
    
    if not primary_comparisons.empty and not secondary_comparisons.empty:
        plt.figure(figsize=(14, 8))
        
        # Plot primary comparisons
        ax1 = plt.subplot(1, 2, 1)
        primary_treatments = primary_comparisons['Treatment']
        x_prim = np.arange(len(primary_treatments))
        
        bars1 = plt.bar(x_prim, primary_comparisons['Treatment_Variants'])
        for i, bar in enumerate(bars1):
            treatment = primary_treatments.iloc[i]
            bar.set_color(color_map.get(treatment, '#2ecc71'))
        
        plt.axhline(y=primary_comparisons['Control_Variants'].iloc[0], color='red', linestyle='--', 
                  label=f'WT-CTRL ({primary_comparisons["Control_Variants"].iloc[0]} variants)')
        
        plt.xlabel('Treatment Group')
        plt.ylabel('Number of Variants')
        plt.title('Primary Comparisons (vs WT-CTRL)', pad=20)
        plt.xticks(x_prim, primary_treatments, rotation=45, ha='right')
        plt.legend()
        
        # Plot secondary comparisons if any
        if len(secondary_comparisons) > 0:
            ax2 = plt.subplot(1, 2, 2)
            sec_treatments = secondary_comparisons['Treatment']
            x_sec = np.arange(len(sec_treatments))
            
            # Create paired bars for treatment and its specific control
            bars2 = plt.bar(x_sec, secondary_comparisons['Treatment_Variants'], width=0.4, label='Treatment')
            ctrl_bars = plt.bar(x_sec + 0.4, secondary_comparisons['Control_Variants'], width=0.4, label='Specific Control')
            
            for i, bar in enumerate(bars2):
                treatment = sec_treatments.iloc[i]
                bar.set_color(color_map.get(treatment, '#2ecc71'))
            
            plt.xlabel('Treatment Group')
            plt.ylabel('Number of Variants')
            plt.title('Secondary Comparisons (with Specific Controls)', pad=20)
            plt.xticks(x_sec + 0.2, sec_treatments, rotation=45, ha='right')
            plt.legend()
        
        plt.tight_layout()
        plt.savefig('analysis/treatment_control_analysis/primary_vs_secondary_comparisons.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Fold Change Visualization
    plt.figure(figsize=(12, 7))
    
    # Debug Point 4: Print fold change data
    print("\nFold change data:")
    print(results_df[['Treatment', 'Control', 'Fold_Change', 'Q_Value']].to_string())
    
    # Create fold change bars with error checking and custom colors
    valid_mask = pd.notna(results_df['Fold_Change'])
    
    # Create bars with custom colors
    bars = plt.bar(np.arange(len(results_df[valid_mask])), 
                  results_df['Fold_Change'][valid_mask])
    
    # Apply treatment-specific colors
    for i, bar in enumerate(bars):
        treatment = results_df['Treatment'].iloc[i]
        bar.set_color(color_map.get(treatment, '#2ecc71'))
        
        # Add significance markers using different shapes/colors
        q_val = results_df['Q_Value'].iloc[i]
        if q_val < 0.001:
            plt.scatter(i, bar.get_height() + 1, marker='*', s=200, color='red')
        elif q_val < 0.01:
            plt.scatter(i, bar.get_height() + 1, marker='*', s=150, color='orange')
        elif q_val < 0.05:
            plt.scatter(i, bar.get_height() + 1, marker='*', s=100, color='yellow')
    
    plt.axhline(y=1, color='black', linestyle='--', alpha=0.3)
    plt.ylabel('Fold Change (Treatment/Control)')
    plt.title('Magnitude of Treatment Effect', pad=20)
    
    # Custom x-tick labels that include control information
    plt.xticks(np.arange(len(results_df[valid_mask])), 
              [f"{t}\n(vs {c})" for t, c in zip(results_df['Treatment'][valid_mask], results_df['Control'][valid_mask])], 
              rotation=45, ha='right')
    
    # Add legend
    legend_elements = [
        Patch(facecolor='red', label='p < 0.001 (***)'),
        Patch(facecolor='orange', label='p < 0.01 (**)'),
        Patch(facecolor='yellow', label='p < 0.05 (*)'),
        Patch(facecolor='gray', label='Not significant')
    ]
    plt.legend(handles=legend_elements, title='Significance Level')
    
    # Add value labels with error checking
    for i, bar in enumerate(bars):
        height = bar.get_height()
        if pd.notna(height):
            plt.text(bar.get_x() + bar.get_width()/2., height - 0.1 * max(results_df['Fold_Change']),
                    f'{height:.1f}×', ha='center', va='top', fontsize=9, 
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=1))
    
    plt.ylim(0, max(results_df['Fold_Change']) * 1.2)  # Add some space at the top
    plt.tight_layout()
    plt.savefig('analysis/treatment_control_analysis/fold_change_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Treatment Types - Biological Classification Visual
    plt.figure(figsize=(14, 8))
    
    # Create categories based on biological characteristics
    categories = {
        'Temperature Adaptation': ['WT-37', 'CAS'],
        'Low Oxygen Adaptation': ['WTA', 'STC']
    }
    
    # Prepare data for grouped bar chart
    cat_data = []
    for category, treatments_list in categories.items():
        for treatment in treatments_list:
            if treatment in results_df['Treatment'].values:
                row = results_df[results_df['Treatment'] == treatment].iloc[0]
                cat_data.append({
                    'Category': category,
                    'Treatment': treatment,
                    'Variants': row['Treatment_Variants'],
                    'Control': row['Control'],
                    'Control_Variants': row['Control_Variants'],
                    'Fold_Change': row['Fold_Change'],
                    'Has_Gene': 'Yes' if treatment in ['STC', 'CAS'] else 'No'
                })
    
    if cat_data:  # Only create this plot if we have categorized data
        df_cat = pd.DataFrame(cat_data)
        
        # Plot grouped by category
        g = sns.catplot(
            data=df_cat, kind="bar",
            x="Category", y="Variants", hue="Treatment",
            palette=color_map, height=6, aspect=1.5, legend=False
        )
        
        plt.title('Biological Classification of Treatments', pad=20)
        plt.ylabel('Number of Variants')
        
        # Add legend with treatment descriptions
        handles, labels = plt.gca().get_legend_handles_labels()
        descriptions = {t: d for t, d in zip(results_df['Treatment'], results_df['Description'])}
        legend_labels = [f"{label}: {descriptions.get(label, '')}" for label in labels]
        plt.legend(handles, legend_labels, title="Treatment", loc='best')
        
        plt.tight_layout()
        plt.savefig('analysis/treatment_control_analysis/biological_classification.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 5. Impact Visualization - Horizontal Bar Chart
    plt.figure(figsize=(14, 8))
    
    # Sort by fold change for better visualization
    sorted_df = results_df.sort_values('Fold_Change', ascending=False)
    
    # Create horizontal bars with custom colors
    bars = plt.barh(np.arange(len(sorted_df)), sorted_df['Fold_Change'])
    
    # Apply treatment-specific colors
    for i, bar in enumerate(bars):
        treatment = sorted_df['Treatment'].iloc[i]
        bar.set_color(color_map.get(treatment, '#2ecc71'))
    
    plt.axvline(x=1, color='black', linestyle='--', alpha=0.3)
    plt.xlabel('Fold Change (Treatment/Control)')
    plt.title('Impact of Treatments Ranked by Effect Size', pad=20)
    
    # Custom y-tick labels with treatment and control info
    plt.yticks(np.arange(len(sorted_df)), 
              [f"{t} vs {c}" for t, c in zip(sorted_df['Treatment'], sorted_df['Control'])])
    
    # Add value labels
    for i, bar in enumerate(bars):
        width = bar.get_width()
        if pd.notna(width):
            plt.text(width + 0.5, bar.get_y() + bar.get_height()/2,
                    f'{width:.1f}×', ha='left', va='center')
    
    plt.tight_layout()
    plt.savefig('analysis/treatment_control_analysis/impact_ranked_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 6. Description Table Plot (visualize treatment descriptions)
    plt.figure(figsize=(12, len(results_df) * 0.8))
    
    # Create a table-like visualization
    table_data = []
    for i, row in results_df.iterrows():
        table_data.append([
            row['Treatment'], 
            row['Description'], 
            row['Control'],
            f"{row['Treatment_Variants']} variants",
            f"{row['Fold_Change']:.1f}×"
        ])
    
    # Create table
    table = plt.table(
        cellText=table_data,
        colLabels=['Treatment', 'Description', 'Control', 'Variants', 'Fold Change'],
        cellLoc='center',
        loc='center'
    )
    
    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 1.5)  # Adjust row height
    
    # Hide axes
    plt.axis('off')
    plt.title('Treatment Descriptions and Key Metrics', pad=20)
    
    plt.tight_layout()
    plt.savefig('analysis/treatment_control_analysis/treatment_descriptions.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    # Load results and print initial debug info
    try:
        results_df = pd.read_csv('analysis/treatment_control_analysis/treatment_vs_control_statistics.csv')
        print("Initial data shape:", results_df.shape)
        print("Columns:", results_df.columns.tolist())
        
        create_treatment_control_visualizations(results_df)
        print("\nVisualizations created successfully!")
    except Exception as e:
        print(f"Error: {e}")
        print("\nPlease run the analysis script first to generate the statistics file.")