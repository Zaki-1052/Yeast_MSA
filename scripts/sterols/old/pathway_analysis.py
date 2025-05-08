import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Define ergosterol pathway relationships
def create_pathway_graph():
    """Create a directed graph representing the ergosterol biosynthetic pathway."""
    G = nx.DiGraph()
    
    # Add nodes (sterol compounds)
    nodes = [
        'Lanosterol', 'Cycloartenol', 'Zymosterol', 'Fecosterol', 
        'Ergosta-7-en-3-ol', 'Ergost-7-en-3beta-ol', 'Ergosterol',
        'Tetrahymanol', 'Stigmasta-5_22-dien-3-ol_acetate'
    ]
    
    for node in nodes:
        G.add_node(node)
    
    # Add edges (pathway steps)
    # These are simplified and would need to be verified with actual pathway information
    edges = [
        ('Lanosterol', 'Zymosterol'),
        ('Zymosterol', 'Fecosterol'),
        ('Fecosterol', 'Ergosta-7-en-3-ol'),
        ('Ergosta-7-en-3-ol', 'Ergosterol'),
        ('Ergost-7-en-3beta-ol', 'Ergosterol'),
        ('Cycloartenol', 'Ergosterol')
    ]
    
    for edge in edges:
        G.add_edge(edge[0], edge[1])
    
    return G

def calculate_pathway_ratios(df, pathway_graph):
    """Calculate ratios between connected sterols in the pathway."""
    results = []
    
    for sample in df['sample'].unique():
        sample_data = df[df['sample'] == sample]
        sample_sterols = {row['sterol']: row['concentration'] for _, row in sample_data.iterrows()}
        
        for source, target in pathway_graph.edges():
            if source in sample_sterols and target in sample_sterols:
                if sample_sterols[source] > 0:
                    ratio = sample_sterols[target] / sample_sterols[source]
                    
                    results.append({
                        'sample': sample,
                        'treatment': sample_data['treatment'].iloc[0] if len(sample_data) > 0 else None,
                        'source_sterol': source,
                        'target_sterol': target,
                        'source_concentration': sample_sterols[source],
                        'target_concentration': sample_sterols[target],
                        'ratio': ratio,
                        'log_ratio': np.log2(ratio)
                    })
    
    return pd.DataFrame(results)

def visualize_pathway_flux(df, pathway_graph, ratios_df):
    """Visualize the ergosterol pathway with flux information."""
    plt.figure(figsize=(15, 10))
    
    # Create positions for pathway visualization
    pos = nx.spring_layout(pathway_graph)
    
    # Get node sizes based on concentrations
    max_conc = df['concentration'].max()
    
    # Create node sizes dictionary
    node_sizes = {}
    for node in pathway_graph.nodes():
        node_data = df[df['sterol'] == node]
        if len(node_data) > 0:
            avg_conc = node_data['concentration'].mean()
            node_sizes[node] = 1000 * (avg_conc / max_conc)
        else:
            node_sizes[node] = 300
    
    # Draw the graph
    nx.draw_networkx_nodes(
        pathway_graph, pos,
        node_size=[node_sizes.get(node, 300) for node in pathway_graph.nodes()],
        node_color='skyblue', alpha=0.8
    )
    
    # Edge weights based on ratios
    edge_widths = []
    for source, target in pathway_graph.edges():
        edge_data = ratios_df[
            (ratios_df['source_sterol'] == source) & 
            (ratios_df['target_sterol'] == target)
        ]
        
        if len(edge_data) > 0:
            avg_ratio = edge_data['ratio'].mean()
            edge_widths.append(max(0.5, avg_ratio * 2))
        else:
            edge_widths.append(1.0)
    
    nx.draw_networkx_edges(pathway_graph, pos, width=edge_widths, alpha=0.7, arrows=True)
    nx.draw_networkx_labels(pathway_graph, pos, font_size=12, font_weight='bold')
    
    plt.title('Ergosterol Pathway Flux Diagram')
    plt.axis('off')
    plt.savefig('results/sterol_analysis/pathway_analysis/pathway_flux_diagram.png')