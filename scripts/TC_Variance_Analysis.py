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
    
    # 1. Bar Plot with Control vs Treatment
    plt.figure(figsize=(12, 6))
    
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
    
    # Customize plot
    plt.xlabel('Treatment Group')
    plt.ylabel('Number of Variants')
    plt.title('Comparison of Variants: Treatment vs Control', pad=20)
    plt.xticks(x, treatments)
    plt.legend()
    
    # Add to the first visualization
    # After creating bars but before annotations
    # Add gradient effect to treatment bars
    for bar in treatment_bars:
        bar.set_facecolor('none')
        bar.set_edgecolor('#2ecc71')
        gradient = np.linspace(0, 1, 100)
        z = np.tile(gradient.reshape(-1, 1), (1, 10))
        bar.set_hatch('/')
        bar.set_alpha(0.7)

    # Add treatment effect arrows
    for i, (control, treatment) in enumerate(zip(control_vals, treatment_vals)):
        plt.arrow(x[i], control, 0, treatment-control,
                head_width=0.1, head_length=20,
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
    plt.savefig('treatment_control_analysis/variant_comparison_barplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Fold Change Visualization
    plt.figure(figsize=(10, 6))
    
    # Debug Point 4: Print fold change data
    print("\nFold change data:")
    print(results_df[['Treatment', 'Fold_Change', 'Q_Value']].to_string())
    
    # Create fold change bars with error checking
    valid_mask = pd.notna(results_df['Fold_Change'])
    colors = ['darkred' if q < 0.001 else 'orangered' if q < 0.01 else 'orange' if q < 0.05 else 'gray' 
             for q in results_df['Q_Value']]
    
    bars = plt.bar(treatments[valid_mask], results_df['Fold_Change'][valid_mask], color=colors)
    
    plt.axhline(y=1, color='black', linestyle='--', alpha=0.3)
    plt.ylabel('Fold Change (Treatment/Control)')
    plt.title('Magnitude of Treatment Effect', pad=20)
    plt.xticks(rotation=45)
    
    # Add legend
    legend_elements = [
        Patch(facecolor='darkred', label='p < 0.001 (***)'),
        Patch(facecolor='orangered', label='p < 0.01 (**)'),
        Patch(facecolor='orange', label='p < 0.05 (*)'),
        Patch(facecolor='gray', label='Not significant')
    ]
    plt.legend(handles=legend_elements, title='Significance Level')
    
    # Add value labels with error checking
    for bar in bars:
        height = bar.get_height()
        if pd.notna(height):
            plt.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.1f}×', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig('treatment_control_analysis/fold_change_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create radial plot
    plt.figure(figsize=(10, 10))
    ax = plt.subplot(111, projection='polar')
    
    # Calculate angles for each treatment
    angles = np.linspace(0, 2*np.pi, len(results_df), endpoint=False)
    
    # Plot fold changes
    fold_changes = results_df['Fold_Change'].values
    max_fold = max(fold_changes)
    bars = ax.bar(angles, fold_changes, alpha=0.5)
    
    # Customize bars
    for angle, bar, q_val in zip(angles, bars, results_df['Q_Value']):
        if q_val < 0.001:
            bar.set_facecolor('#e74c3c')
        elif q_val < 0.01:
            bar.set_facecolor('#f39c12')
        else:
            bar.set_facecolor('#3498db')
            
        # Add value labels
        label_radius = bar.get_height() + max_fold * 0.1
        plt.text(angle, label_radius, f'{bar.get_height():.1f}×',
                ha='center', va='center')
    
    # Set treatment labels
    ax.set_xticks(angles)
    ax.set_xticklabels(results_df['Treatment'])
    
    plt.title('Treatment Effect Magnitude\n(Radial View)', pad=20)
    plt.savefig('treatment_control_analysis/radial_effect.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    plt.figure(figsize=(12, 8))
    
    # Create flow-like visualization
    treatments = results_df['Treatment']
    y_base = np.zeros(len(treatments))
    
    # Create flowing effect
    x = np.arange(len(treatments))
    for i in range(20):
        y = results_df['Treatment_Variants'] * np.sin(np.linspace(0, 2*np.pi, len(treatments)) + i/10)
        plt.fill_between(x, y_base, y + results_df['Treatment_Variants'],
                        alpha=0.1, color='#3498db')
    
    # Add main bars
    plt.bar(x, results_df['Treatment_Variants'],
           alpha=0.3, color='#2ecc71')
    
    # Add significance indicators
    for i, row in enumerate(results_df.itertuples()):
        if row.Q_Value < 0.001:
            plt.scatter([i], [row.Treatment_Variants * 1.1],
                       marker='*', s=200, color='#e74c3c')
    
    plt.xlabel('Treatment')
    plt.ylabel('Number of Variants')
    plt.title('Treatment Effect Flow Diagram', pad=20)
    plt.xticks(x, treatments)
    
    plt.savefig('treatment_control_analysis/flow_diagram.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Impact Visualization
    plt.figure(figsize=(12, 6))
    
    # Calculate percentage increase with error checking
    percent_increase = np.where(
        pd.notna(results_df['Control_Variants']) & (results_df['Control_Variants'] != 0),
        (results_df['Treatment_Variants'] - results_df['Control_Variants']) / results_df['Control_Variants'] * 100,
        np.nan
    )
    
    # Debug Point 5: Print percentage increase data
    print("\nPercentage increase data:")
    print(pd.DataFrame({'Treatment': treatments, 'Percent_Increase': percent_increase}).to_string())
    
    # Create horizontal bars for valid data only
    valid_mask = pd.notna(percent_increase)
    bars = plt.barh(treatments[valid_mask], percent_increase[valid_mask])
    
    # Color bars based on significance
    for bar, q_val in zip(bars, results_df['Q_Value'][valid_mask]):
        if pd.notna(q_val):
            if q_val < 0.001:
                bar.set_color('darkred')
            elif q_val < 0.01:
                bar.set_color('orangered')
            elif q_val < 0.05:
                bar.set_color('orange')
            else:
                bar.set_color('gray')
    
    plt.axvline(x=0, color='black', linestyle='--', alpha=0.3)
    plt.xlabel('Percentage Increase in Variants (%)')
    plt.title('Impact of Treatment on Variant Count', pad=20)
    
    # Add value labels with error checking
    for bar in bars:
        width = bar.get_width()
        if pd.notna(width):
            plt.text(width + (5 if width >= 0 else -5), 
                    bar.get_y() + bar.get_height()/2.,
                    f'{width:,.0f}%', 
                    ha='left' if width >= 0 else 'right', 
                    va='center')
    
    plt.legend(handles=legend_elements, title='Significance Level',
              bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig('treatment_control_analysis/percent_increase_plot.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    # Load results and print initial debug info
    results_df = pd.read_csv('treatment_control_analysis/treatment_vs_control_statistics.csv')
    print("Initial data shape:", results_df.shape)
    print("Columns:", results_df.columns.tolist())
    
    create_treatment_control_visualizations(results_df)
    print("\nVisualizations created successfully!")