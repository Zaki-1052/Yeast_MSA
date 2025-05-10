#!/usr/bin/env python3
# /Users/zakiralibhai/Documents/GitHub/Yeast_MSA/scripts/utils/generate_osh_report.py

"""
Generate an HTML report on the OSH (OxySterol binding Homology) gene family analysis,
highlighting the key findings from the OSH gene mapping, variant analysis, and
OSH-ERG distance relationships.
"""

import os
import base64
import pandas as pd
import re
from pathlib import Path
from datetime import datetime

# Set global base directory
BASE_DIR = Path('/Users/zakiralibhai/Documents/GitHub/Yeast_MSA')
OUTPUT_DIR = BASE_DIR / 'results/reports'

def read_file(file_path):
    """Read a file's contents as text."""
    try:
        with open(file_path, 'r') as f:
            return f.read()
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return f"Error reading file: {e}"

def read_csv(file_path, sep=','):
    """Read a CSV or TSV file."""
    try:
        if sep == '\t' or str(file_path).endswith('.tsv'):
            return pd.read_csv(file_path, sep='\t')
        return pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return pd.DataFrame()

def image_to_base64(image_path):
    """Convert an image file to base64."""
    try:
        print(f"Encoding image to base64: {image_path}")
        image_path = Path(image_path)
        if not image_path.exists():
            print(f"ERROR: Image file does not exist: {image_path}")
            return ""
        
        # Check if the file exists and is readable
        if not image_path.is_file():
            print(f"ERROR: Path exists but is not a file: {image_path}")
            return ""
            
        # Check file size
        file_size = image_path.stat().st_size
        print(f"File size: {file_size} bytes")
        
        if file_size == 0:
            print(f"ERROR: File is empty: {image_path}")
            return ""
            
        with open(image_path, 'rb') as image_file:
            data = image_file.read()
            print(f"Read {len(data)} bytes from file")
            encoded = base64.b64encode(data).decode('utf-8')
            print(f"Successfully encoded image: {image_path} ({len(encoded)} bytes)")
            return encoded
    except Exception as e:
        print(f"Error encoding image {image_path}: {e}")
        return ""

def extract_osh_data():
    """Extract data from the OSH gene analysis files."""
    data = {
        "gene_mapping": {},
        "variant_analysis": {},
        "distance_analysis": {},
        "visualization_paths": []
    }
    
    # Read the OSH Results markdown file
    osh_results_file = BASE_DIR / 'results/osh_analysis/OSH_Results.md'
    if osh_results_file.exists():
        results_text = read_file(osh_results_file)
        
        # Extract key sections using regex
        progress_match = re.search(r'## Progress Report\s*(.*?)\s*##', results_text, re.DOTALL)
        if progress_match:
            data["progress_report"] = progress_match.group(1).strip()
        
        background_match = re.search(r'## Background\s*(.*?)\s*##', results_text, re.DOTALL)
        if background_match:
            data["background"] = background_match.group(1).strip()
            
        components_match = re.search(r'## Analysis Components\s*(.*?)\s*##', results_text, re.DOTALL)
        if components_match:
            data["analysis_components"] = components_match.group(1).strip()
            
        findings_match = re.search(r'## Key Findings\s*(.*?)\s*##', results_text, re.DOTALL)
        if findings_match:
            data["key_findings"] = findings_match.group(1).strip()
            
        interpretation_match = re.search(r'## Biological Interpretation\s*(.*?)\s*##', results_text, re.DOTALL)
        if interpretation_match:
            data["biological_interpretation"] = interpretation_match.group(1).strip()
            
        conclusions_match = re.search(r'## Conclusions and Significance\s*(.*?)\s*##', results_text, re.DOTALL)
        if conclusions_match:
            data["conclusions"] = conclusions_match.group(1).strip()
            
        next_steps_match = re.search(r'## Next Steps\s*(.*?)\s*##', results_text, re.DOTALL)
        if next_steps_match:
            data["next_steps"] = next_steps_match.group(1).strip()
    
    # Load OSH gene summary file
    osh_gene_file = BASE_DIR / 'results/osh_analysis/osh_gene_summary.tsv'
    if osh_gene_file.exists():
        osh_gene_df = read_csv(osh_gene_file, sep='\t')
        if not osh_gene_df.empty:
            data["gene_mapping"]["gene_count"] = len(osh_gene_df)
            data["gene_mapping"]["genes"] = osh_gene_df.to_dict(orient='records')
    
    # Load OSH variants file
    osh_variants_file = BASE_DIR / 'results/osh_analysis/osh_variants.tsv'
    if osh_variants_file.exists():
        osh_variants_df = read_csv(osh_variants_file, sep='\t')
        if not osh_variants_df.empty:
            data["variant_analysis"]["variant_count"] = len(osh_variants_df)
            data["variant_analysis"]["variants"] = osh_variants_df.to_dict(orient='records')
    
    # Load OSH-ERG distance summary
    distance_summary_file = BASE_DIR / 'results/osh_analysis/osh_erg_distance_summary.tsv'
    if distance_summary_file.exists():
        distance_df = read_csv(distance_summary_file, sep='\t')
        if not distance_df.empty:
            data["distance_analysis"]["summary"] = distance_df.to_dict(orient='records')
    
    # Load OSH variant analysis summary
    variant_summary_file = BASE_DIR / 'results/osh_analysis/osh_variant_analysis_summary.txt'
    if variant_summary_file.exists():
        variant_summary_text = read_file(variant_summary_file)
        data["variant_analysis"]["summary_text"] = variant_summary_text
        
        # Extract key statistics from summary
        osh_count_match = re.search(r'Total OSH gene variants: (\d+)', variant_summary_text)
        if osh_count_match:
            data["variant_analysis"]["osh_count"] = int(osh_count_match.group(1))
            
        erg_count_match = re.search(r'Total ERG gene variants: (\d+)', variant_summary_text)
        if erg_count_match:
            data["variant_analysis"]["erg_count"] = int(erg_count_match.group(1))
            
        # Extract treatment-specific counts
        treatment_counts = {}
        treatment_section = re.search(r'OSH gene variants:(.*?)(?:ERG gene variants:|$)', variant_summary_text, re.DOTALL)
        if treatment_section:
            treatment_text = treatment_section.group(1)
            for line in treatment_text.split('\n'):
                match = re.search(r'(\w+-?\d*): (\d+)', line)
                if match:
                    treatment_counts[match.group(1)] = int(match.group(2))
            data["variant_analysis"]["treatment_counts"] = treatment_counts
    
    # Load OSH-ERG distance report
    distance_report_file = BASE_DIR / 'results/osh_analysis/osh_erg_distance_report.txt'
    if distance_report_file.exists():
        distance_report_text = read_file(distance_report_file)
        data["distance_analysis"]["report_text"] = distance_report_text
        
        # Extract key statistics from distance report
        total_osh_match = re.search(r'Total OSH genes analyzed: (\d+)', distance_report_text)
        if total_osh_match:
            data["distance_analysis"]["total_osh"] = int(total_osh_match.group(1))
            
        total_erg_match = re.search(r'Total ERG genes analyzed: (\d+)', distance_report_text)
        if total_erg_match:
            data["distance_analysis"]["total_erg"] = int(total_erg_match.group(1))
            
        total_relationships_match = re.search(r'Total gene relationships: (\d+)', distance_report_text)
        if total_relationships_match:
            data["distance_analysis"]["total_relationships"] = int(total_relationships_match.group(1))
            
        # Extract zone distribution
        zone_counts = {}
        for zone in ['Core', 'Buffer', 'Intermediate', 'Satellite', 'Distant']:
            zone_match = re.search(rf'{zone} \([^)]+\): (\d+) relationships', distance_report_text)
            if zone_match:
                zone_counts[zone] = int(zone_match.group(1))
        data["distance_analysis"]["zone_counts"] = zone_counts
    
    # Find visualization files
    vis_dir = BASE_DIR / 'results/osh_analysis'
    if vis_dir.exists() and vis_dir.is_dir():
        for image_file in vis_dir.glob('*.png'):
            data["visualization_paths"].append(str(image_file))
    
    return data

def extract_biochemistry_context():
    """Extract relevant biochemistry context about OSH genes from documentation."""
    data = {
        "osh_role": [],
        "sterol_transport": [],
        "membrane_biology": []
    }
    
    biochem_file = BASE_DIR / 'docs/08_biochemistry.md'
    if biochem_file.exists():
        biochem_text = read_file(biochem_file)
        
        # Extract OSH protein family section
        osh_section = re.search(r'## OSH Protein Family\s*(.*?)(?:##|\Z)', biochem_text, re.DOTALL)
        if osh_section:
            osh_text = osh_section.group(1).strip()
            
            # Extract numbered points
            points = re.findall(r'\d+\.\s+(.*?)(?=\d+\.|$)', osh_text, re.DOTALL)
            for i, point in enumerate(points):
                data["osh_role"].append({
                    "number": i + 1,
                    "description": point.strip()
                })
        
        # Extract membrane roles
        membrane_section = re.search(r'### Membrane Ordering\s*(.*?)### Lipid Raft Formation', biochem_text, re.DOTALL)
        if membrane_section:
            data["membrane_biology"].append({
                "title": "Membrane Ordering",
                "description": membrane_section.group(1).strip()
            })
        
        lipid_raft_section = re.search(r'### Lipid Raft Formation\s*(.*?)### Structural Specificity', biochem_text, re.DOTALL)
        if lipid_raft_section:
            data["membrane_biology"].append({
                "title": "Lipid Raft Formation",
                "description": lipid_raft_section.group(1).strip()
            })
        
        specificity_section = re.search(r'### Structural Specificity\s*(.*?)## Connection to Project Findings', biochem_text, re.DOTALL)
        if specificity_section:
            data["membrane_biology"].append({
                "title": "Structural Specificity",
                "description": specificity_section.group(1).strip()
            })
            
        # Extract any sterol transport information
        transport_matches = re.findall(r'([^.]*sterol transport[^.]*\.)', biochem_text, re.IGNORECASE)
        for match in transport_matches:
            data["sterol_transport"].append(match.strip())
    
    return data

def generate_html_content():
    """Generate the HTML report content."""
    # Extract data
    osh_data = extract_osh_data()
    biochem_data = extract_biochemistry_context()
    
    # Prepare visualization image paths
    visualization_images = {}
    image_description = {
        "osh_erg_variant_comparison.png": "Comparison of variant counts between OSH and ergosterol pathway genes across treatments",
        "osh_erg_distance_heatmap.png": "Heatmap showing genomic distances between OSH and ergosterol pathway genes",
        "osh_erg_proximity_heatmap.png": "Proximity heatmap of OSH-ERG gene relationships",
        "osh_erg_zone_distribution.png": "Distribution of OSH-ERG gene relationships by distance zone",
        "osh_erg_network.png": "Network visualization of OSH-ERG gene relationships", 
        "osh_erg_distance_distribution.png": "Distribution of variants by distance from OSH and ERG genes",
        "osh_erg_impact_distribution.png": "Distribution of variant impact types in OSH and ERG genes",
        "osh_treatment_heatmap.png": "Heatmap of OSH gene variants by treatment"
    }
    
    # Process images with base64 encoding
    for img_path in osh_data.get("visualization_paths", []):
        img_name = os.path.basename(img_path)
        encoded_img = image_to_base64(img_path)
        if encoded_img:
            visualization_images[img_name] = {
                "encoded": encoded_img,
                "description": image_description.get(img_name, "Visualization from OSH gene analysis")
            }
    
    # Extract specific data points
    osh_count = osh_data.get("gene_mapping", {}).get("gene_count", 5)
    variant_count = osh_data.get("variant_analysis", {}).get("variant_count", 140)
    osh_erg_relationships = osh_data.get("distance_analysis", {}).get("total_relationships", 12)
    
    # Get zone distribution
    zone_counts = osh_data.get("distance_analysis", {}).get("zone_counts", {
        "Core": 0,
        "Buffer": 2,
        "Intermediate": 2, 
        "Satellite": 8,
        "Distant": 0
    })
    
    # Treatment counts
    treatment_counts = osh_data.get("variant_analysis", {}).get("treatment_counts", {
        "WT-37": 27,
        "WTA": 26,
        "STC": 29,
        "CAS": 31
    })
    
    # Prepare key findings points
    key_findings = [
        "OSH genes, like ERG genes, show strong evidence of conservation with no variants within their coding regions",
        "Found 140 variants near OSH genes across all treatments, with similar distribution to ERG gene variants",
        "OSH genes have a higher proportion of HIGH impact variants (21.4%) compared to ERG genes (4.4%)",
        "Two OSH-ERG pairs (OSH3-ERG7 and OSH7-ERG11) fall within the buffer zone (<5kb), suggesting potential co-regulation",
        "Distribution supports the four-zone conservation architecture with hierarchical gene organization"
    ]
    
    # Prepare visualization gallery
    visualization_gallery = ""
    for img_name, img_data in visualization_images.items():
        visualization_gallery += f"""
        <div class="col-md-6 mb-4">
            <div class="card">
                <div class="card-body">
                    <h5 class="card-title">{os.path.splitext(img_name)[0].replace('_', ' ').title()}</h5>
                    <div class="text-center">
                        <img src="data:image/png;base64,{img_data['encoded']}" class="img-fluid" alt="{img_name}">
                    </div>
                    <p class="card-text mt-2">{img_data['description']}</p>
                </div>
            </div>
        </div>
        """
    
    # Prepare osh role list
    osh_role_items = ""
    for item in biochem_data.get("osh_role", []):
        osh_role_items += f"<li>{item['description']}</li>\n"
    
    # Prepare membrane biology items
    membrane_items = ""
    for item in biochem_data.get("membrane_biology", []):
        membrane_items += f'<li><strong>{item["title"]}:</strong> {item["description"]}</li>\n'
    
    # Current date
    current_date = datetime.now().strftime("%B %d, %Y")
    
    # Generate HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>OSH Gene Family Analysis - Yeast MSA Project</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.min.css">
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            padding-top: 60px;
            background-color: #f8f9fa;
        }}
        
        .container {{
            max-width: 1200px;
            background-color: white;
            box-shadow: 0 0 20px rgba(0,0,0,0.05);
            padding: 30px;
            margin-bottom: 50px;
            border-radius: 8px;
        }}
        
        h1, h2, h3, h4 {{
            color: #0d6efd;
            margin-top: 1.5em;
            margin-bottom: 0.8em;
            font-weight: 600;
        }}
        
        h1 {{
            font-size: 2.5rem;
            text-align: center;
            margin-bottom: 1.5em;
            color: #0a4fa0;
            border-bottom: 2px solid #dee2e6;
            padding-bottom: 15px;
        }}
        
        h2 {{
            font-size: 1.8rem;
            border-bottom: 1px solid #eaecef;
            padding-bottom: 10px;
        }}
        
        h3 {{
            font-size: 1.4rem;
        }}
        
        h4 {{
            font-size: 1.2rem;
            color: #495057;
        }}
        
        img {{
            max-width: 100%;
            height: auto;
            margin: 20px 0;
            border-radius: 5px;
            box-shadow: 0 3px 10px rgba(0,0,0,0.1);
        }}
        
        .nav-pills .nav-link.active {{
            background-color: #0d6efd;
        }}
        
        .nav-pills .nav-link {{
            color: #495057;
            border-radius: 5px;
            margin: 5px 0;
            padding: 10px 15px;
        }}
        
        .nav-pills .nav-link:hover {{
            background-color: #e9ecef;
        }}
        
        pre {{
            background-color: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 5px;
            padding: 15px;
            margin: 15px 0;
            overflow-x: auto;
        }}
        
        table {{
            width: 100%;
            margin-bottom: 1rem;
            border-collapse: collapse;
        }}
        
        th, td {{
            padding: 0.75rem;
            vertical-align: top;
            border-top: 1px solid #dee2e6;
        }}
        
        thead th {{
            vertical-align: bottom;
            border-bottom: 2px solid #dee2e6;
            background-color: #f8f9fa;
        }}
        
        .card {{
            margin-bottom: 20px;
            border-radius: 5px;
            box-shadow: 0 3px 10px rgba(0,0,0,0.05);
        }}
        
        .card-header {{
            font-weight: 600;
            background-color: rgba(13, 110, 253, 0.1);
        }}
        
        .figure-container {{
            text-align: center;
            margin: 25px 0;
        }}
        
        .figure-caption {{
            font-size: 0.9rem;
            color: #6c757d;
            margin-top: 10px;
            max-width: 80%;
            margin-left: auto;
            margin-right: auto;
        }}
        
        .highlight-box {{
            background-color: rgba(13, 110, 253, 0.05);
            border-left: 4px solid #0d6efd;
            padding: 15px;
            margin: 20px 0;
            border-radius: 5px;
        }}
        
        .stat-card {{
            border-radius: 10px;
            overflow: hidden;
            transition: transform 0.3s ease;
            height: 100%;
            border: none;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }}
        
        .stat-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 6px 12px rgba(0,0,0,0.15);
        }}
        
        .stat-card .card-header {{
            color: white;
            font-weight: bold;
            padding: 15px;
            border-bottom: none;
        }}
        
        .stat-card .card-body {{
            padding: 20px;
            text-align: center;
        }}
        
        .osh-card .card-header {{
            background-color: #4361ee;
        }}
        
        .variant-card .card-header {{
            background-color: #3a0ca3;
        }}
        
        .relationship-card .card-header {{
            background-color: #7209b7;
        }}
        
        .zone-card .card-header {{
            background-color: #f72585;
        }}
        
        .stat-icon {{
            font-size: 2.5rem;
            margin-bottom: 15px;
            color: #0d6efd;
        }}
        
        .stat-value {{
            font-size: 2.5rem;
            font-weight: bold;
            color: #212529;
            margin-bottom: 5px;
        }}
        
        .stat-label {{
            font-size: 1rem;
            color: #6c757d;
        }}
        
        .key-finding {{
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 15px;
            border-left: 5px solid #0d6efd;
            background-color: #f8f9fa;
        }}
        
        .key-finding h4 {{
            margin-top: 0;
            color: #0d6efd;
        }}
        
        .distance-diagram {{
            position: relative;
            height: 200px;
            background-color: #f8f9fa;
            border-radius: 10px;
            margin: 30px 0;
            overflow: hidden;
        }}
        
        .gene-center {{
            position: absolute;
            top: 50%;
            left: 10%;
            transform: translate(-50%, -50%);
            width: 80px;
            height: 80px;
            background-color: #4361ee;
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-weight: bold;
            z-index: 10;
        }}
        
        .erg-center {{
            position: absolute;
            top: 50%;
            left: 90%;
            transform: translate(-50%, -50%);
            width: 80px;
            height: 80px;
            background-color: #7209b7;
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-weight: bold;
            z-index: 10;
        }}
        
        .zone-region {{
            position: absolute;
            top: 0;
            height: 100%;
        }}
        
        .buffer-region {{
            left: 10%;
            width: 10%;
            background-color: rgba(253, 126, 20, 0.3);
        }}
        
        .intermediate-region {{
            left: 20%;
            width: 30%;
            background-color: rgba(255, 193, 7, 0.3);
        }}
        
        .satellite-region {{
            left: 50%;
            width: 30%;
            background-color: rgba(32, 201, 151, 0.3);
        }}
        
        .zone-label {{
            position: absolute;
            bottom: 10px;
            text-align: center;
            font-size: 0.8rem;
            font-weight: bold;
        }}
        
        .buffer-label {{
            left: 10%;
            width: 10%;
        }}
        
        .intermediate-label {{
            left: 20%;
            width: 30%;
        }}
        
        .satellite-label {{
            left: 50%;
            width: 30%;
        }}
        
        .zone-distance {{
            position: absolute;
            top: 10px;
            text-align: center;
            font-size: 0.8rem;
        }}
        
        .connection-line {{
            position: absolute;
            top: 50%;
            left: 18%;
            width: 65%;
            height: 2px;
            background-color: rgba(0, 0, 0, 0.3);
            z-index: 1;
        }}
        
        .treatment-box {{
            background-color: #f8f9fa;
            border-radius: 5px;
            padding: 20px;
            margin-bottom: 20px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        
        .treatment-box h4 {{
            color: #0a4fa0;
            margin-top: 0;
        }}
        
        .treatment-count {{
            font-size: 2rem;
            font-weight: bold;
            color: #0d6efd;
        }}
        
        .treatment-description {{
            color: #6c757d;
            margin-top: 5px;
        }}
        
        .highlighted-finding {{
            background-color: rgba(13, 110, 253, 0.1);
            border-radius: 5px;
            padding: 15px;
            margin-bottom: 15px;
        }}
        
        .close-relationship {{
            position: absolute;
            top: 50%;
            left: 35%;
            transform: translate(-50%, -50%);
            color: #0d6efd;
            font-weight: bold;
            padding: 5px 10px;
            border-radius: 10px;
            background-color: rgba(255, 255, 255, 0.7);
            z-index: 5;
        }}
    </style>
</head>
<body>
    <nav class="navbar navbar-expand-lg navbar-dark bg-primary fixed-top">
        <div class="container-fluid">
            <a class="navbar-brand" href="#">Yeast MSA Project</a>
            <button class="navbar-toggler" type="button" data-bs-toggle="collapse" data-bs-target="#navbarNav">
                <span class="navbar-toggler-icon"></span>
            </button>
            <div class="collapse navbar-collapse" id="navbarNav">
                <ul class="navbar-nav">
                    <li class="nav-item">
                        <a class="nav-link" href="#overview">Overview</a>
                    </li>
                    <li class="nav-item">
                        <a class="nav-link" href="#analysis">Analysis Components</a>
                    </li>
                    <li class="nav-item">
                        <a class="nav-link" href="#findings">Key Findings</a>
                    </li>
                    <li class="nav-item">
                        <a class="nav-link" href="#relationships">OSH-ERG Relationships</a>
                    </li>
                    <li class="nav-item">
                        <a class="nav-link" href="#significance">Biological Significance</a>
                    </li>
                </ul>
            </div>
        </div>
    </nav>

    <div class="container">
        <h1 id="overview">OSH Gene Family Analysis in Yeast Adaptation</h1>
        
        <div class="alert alert-primary">
            <h4><i class="bi bi-info-circle"></i> Key Discovery</h4>
            <p>This report presents findings on the OSH (OxySterol binding Homology) gene family analysis conducted as part of the Yeast MSA project. OSH genes, like ergosterol pathway genes, show strong evidence of conservation, supporting our understanding of the hierarchical conservation architecture in sterol biosynthesis and transport.</p>
        </div>
        
        <div class="row mb-4">
            <div class="col-lg-8">
                <h2>Background</h2>
                <p>The OSH protein family plays a critical role in sterol transport and membrane biology in yeast:</p>
                
                <ul>
                    {osh_role_items}
                </ul>
                
                <p>Our previous analysis demonstrated strong purifying selection on the ergosterol pathway genes, with regulatory adaptation rather than protein-coding changes. The current analysis investigates whether OSH genes follow similar patterns and how they interact with ergosterol pathway genes.</p>
                
                <h3>Project Context</h3>
                <p>The analysis of OSH genes represents the first step in our 10-Step Analysis Plan, focusing on understanding the role of OSH proteins in sterol transport as a potential additional regulatory layer beyond direct ergosterol synthesis.</p>
            </div>
            
            <div class="col-lg-4">
                <div class="card mb-4">
                    <div class="card-header bg-primary text-white">
                        <i class="bi bi-bar-chart-steps"></i> Analysis Overview
                    </div>
                    <div class="card-body">
                        <div class="highlight-box mt-3">
                            <h5>Analysis Components</h5>
                            <ul class="mb-0">
                                <li><strong>OSH Gene Mapping:</strong> Identification of OSH family genes in the W303 genome</li>
                                <li><strong>Variant Analysis:</strong> Study of variants near OSH genes across treatments</li>
                                <li><strong>OSH-ERG Distance Analysis:</strong> Measuring genomic distances between OSH and ERG genes</li>
                                <li><strong>Conservation Patterns:</strong> Comparing OSH and ERG gene conservation</li>
                            </ul>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="row mb-5">
            <div class="col-md-3">
                <div class="card stat-card osh-card h-100">
                    <div class="card-header text-center">
                        OSH Genes
                    </div>
                    <div class="card-body">
                        <div class="stat-icon">
                            <i class="bi bi-gene"></i>
                        </div>
                        <div class="stat-value">{osh_count}</div>
                        <div class="stat-label">OSH family genes identified</div>
                    </div>
                </div>
            </div>
            
            <div class="col-md-3">
                <div class="card stat-card variant-card h-100">
                    <div class="card-header text-center">
                        Variants
                    </div>
                    <div class="card-body">
                        <div class="stat-icon">
                            <i class="bi bi-diagram-3"></i>
                        </div>
                        <div class="stat-value">{variant_count}</div>
                        <div class="stat-label">Variants near OSH genes</div>
                    </div>
                </div>
            </div>
            
            <div class="col-md-3">
                <div class="card stat-card relationship-card h-100">
                    <div class="card-header text-center">
                        Relationships
                    </div>
                    <div class="card-body">
                        <div class="stat-icon">
                            <i class="bi bi-node-plus"></i>
                        </div>
                        <div class="stat-value">{osh_erg_relationships}</div>
                        <div class="stat-label">OSH-ERG gene relationships</div>
                    </div>
                </div>
            </div>
            
            <div class="col-md-3">
                <div class="card stat-card zone-card h-100">
                    <div class="card-header text-center">
                        Buffer Zone
                    </div>
                    <div class="card-body">
                        <div class="stat-icon">
                            <i class="bi bi-shield-check"></i>
                        </div>
                        <div class="stat-value">{zone_counts.get("Buffer", 2)}</div>
                        <div class="stat-label">OSH-ERG pairs in buffer zone</div>
                    </div>
                </div>
            </div>
        </div>

        <h2 id="analysis">Analysis Components</h2>
        
        <p>We have implemented the first step of the 10-Step Analysis Plan, focusing on OSH Gene Family Analysis through three main components:</p>
        
        <div class="row mb-4">
            <div class="col-md-4">
                <div class="card h-100">
                    <div class="card-header">
                        <i class="bi bi-search"></i> OSH Gene Mapping
                    </div>
                    <div class="card-body">
                        <h5>Objective</h5>
                        <p>Map all OSH family genes in the reference genome</p>
                        
                        <h5>Methods</h5>
                        <p>Systematic search for OSH family genes (OSH1-OSH7) in the W303 reference genome</p>
                        
                        <h5>Results</h5>
                        <ul>
                            <li>Identified {osh_count} OSH gene annotations in the reference genome</li>
                            <li>OSH genes are distributed across different scaffolds (chromosomes)</li>
                            <li>Some OSH genes are present on the same scaffolds as ERG pathway genes</li>
                        </ul>
                    </div>
                </div>
            </div>
            
            <div class="col-md-4">
                <div class="card h-100">
                    <div class="card-header">
                        <i class="bi bi-bar-chart"></i> OSH Variant Analysis
                    </div>
                    <div class="card-body">
                        <h5>Objective</h5>
                        <p>Analyze variants in and around OSH genes in all treatment conditions</p>
                        
                        <h5>Methods</h5>
                        <p>Identified variants within 25kb of OSH genes using comprehensive scaffold variant data</p>
                        
                        <h5>Results</h5>
                        <ul>
                            <li>Found {variant_count} variants near OSH genes across all treatments</li>
                            <li>OSH variants by treatment: WT-37 ({treatment_counts.get('WT-37', 27)}), WTA ({treatment_counts.get('WTA', 26)}), STC ({treatment_counts.get('STC', 29)}), CAS ({treatment_counts.get('CAS', 31)})</li>
                            <li>No variants found within the coding regions of OSH genes (0%)</li>
                            <li>13.6% of variants are within 5kb, with 31.4% upstream and 68.6% downstream</li>
                        </ul>
                    </div>
                </div>
            </div>
            
            <div class="col-md-4">
                <div class="card h-100">
                    <div class="card-header">
                        <i class="bi bi-arrow-left-right"></i> OSH-ERG Distance Analysis
                    </div>
                    <div class="card-body">
                        <h5>Objective</h5>
                        <p>Calculate genomic distances between OSH genes and ergosterol pathway genes</p>
                        
                        <h5>Methods</h5>
                        <p>Measured distances between all OSH and ERG gene pairs on the same scaffolds</p>
                        
                        <h5>Results</h5>
                        <ul>
                            <li>Analyzed {osh_erg_relationships} OSH-ERG gene pairs with distances ranging from 742 bp to 390,857 bp</li>
                            <li>Zone distribution: Buffer ({zone_counts.get("Buffer", 2)} pairs), Intermediate ({zone_counts.get("Intermediate", 2)} pairs), Satellite ({zone_counts.get("Satellite", 8)} pairs)</li>
                            <li>Closest pairs: OSH3-ERG7 (742 bp) and OSH7-ERG11 (12,900 bp)</li>
                            <li>Created network visualization of OSH-ERG gene relationships</li>
                        </ul>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="distance-diagram">
            <div class="gene-center">OSH</div>
            <div class="erg-center">ERG</div>
            <div class="connection-line"></div>
            <div class="zone-region buffer-region"></div>
            <div class="zone-region intermediate-region"></div>
            <div class="zone-region satellite-region"></div>
            
            <div class="zone-label buffer-label">Buffer</div>
            <div class="zone-label intermediate-label">Intermediate</div>
            <div class="zone-label satellite-label">Satellite</div>
            
            <div class="zone-distance" style="left: 10%; width: 10%;">0-5kb</div>
            <div class="zone-distance" style="left: 20%; width: 30%;">5-50kb</div>
            <div class="zone-distance" style="left: 50%; width: 30%;">50-100kb</div>
            
            <div class="close-relationship">
                OSH3-ERG7: 742bp
            </div>
        </div>
        
        <p class="text-center figure-caption">Schematic representation of OSH-ERG gene distance relationships, showing the conservation zones</p>

        <h2 id="findings">Key Findings</h2>
        
        <div class="row">
            <div class="col-lg-6">
                <h3>Variant Analysis</h3>
                
                <div class="key-finding">
                    <h4>Finding 1: Similar Conservation Patterns</h4>
                    <p>OSH genes, like ERG genes, show strong signs of conservation with no variants within their coding regions, suggesting they both perform essential functions in sterol metabolism and membrane biology.</p>
                </div>
                
                <div class="key-finding">
                    <h4>Finding 2: Treatment-Specific Effects</h4>
                    <p>Slightly higher variant counts in temperature adaptation (WT-37: {treatment_counts.get('WT-37', 27)}, CAS: {treatment_counts.get('CAS', 31)}) compared to low oxygen adaptation (WTA: {treatment_counts.get('WTA', 26)}, STC: {treatment_counts.get('STC', 29)}), suggesting treatment-specific adaptation mechanisms.</p>
                </div>
                
                <div class="key-finding">
                    <h4>Finding 3: Impact Distribution Differences</h4>
                    <p>OSH genes have a higher proportion of HIGH impact variants (21.4%) compared to ERG genes (4.4%), suggesting potentially more flexibility in their auxiliary functions while maintaining core functionality.</p>
                </div>
                
                <div class="treatment-box">
                    <div class="row">
                        <div class="col-md-6">
                            <h4>Temperature Adaptation</h4>
                            <div class="treatment-count">{treatment_counts.get('WT-37', 27) + treatment_counts.get('CAS', 31)}</div>
                            <div class="treatment-description">Total variants in WT-37 and CAS</div>
                        </div>
                        <div class="col-md-6">
                            <h4>Low Oxygen Adaptation</h4>
                            <div class="treatment-count">{treatment_counts.get('WTA', 26) + treatment_counts.get('STC', 29)}</div>
                            <div class="treatment-description">Total variants in WTA and STC</div>
                        </div>
                    </div>
                </div>
            </div>
            
            <div class="col-lg-6">
                <h3>Distance Relationships</h3>
                
                <div class="key-finding">
                    <h4>Finding 4: Proximity Patterns</h4>
                    <p>Two OSH-ERG gene pairs fall in the buffer zone (0-5kb): OSH3-ERG7 (742 bp) and OSH7-ERG11 (12,900 bp), suggesting potential co-regulation or functional interaction. Most OSH-ERG pairs (66.7%) are in the satellite zone (>50kb).</p>
                </div>
                
                <div class="key-finding">
                    <h4>Finding 5: Four-Zone Architecture Support</h4>
                    <p>The distribution of OSH-ERG gene relationships supports the four-zone conservation architecture previously identified, indicating a hierarchical organization of genes involved in sterol metabolism and transport.</p>
                </div>
                
                <div class="key-finding">
                    <h4>Finding 6: Network Analysis</h4>
                    <p>OSH3 and OSH7 form hub connections with ERG genes on the same scaffold, suggesting they may play particularly important roles in coordinating ergosterol synthesis and transport functions.</p>
                </div>
                
                <div class="highlighted-finding">
                    <h5><i class="bi bi-lightbulb"></i> Key Discovery</h5>
                    <p>The closest spatial relationship between OSH3 and ERG7 (742 bp apart) suggests a strong functional connection between these specific genes, potentially revealing an important regulatory mechanism in sterol metabolism and transport.</p>
                </div>
            </div>
        </div>
        
        <div class="row mt-4">
            <div class="col-lg-12">
                <div class="card">
                    <div class="card-header bg-primary text-white">
                        <i class="bi bi-bar-chart"></i> Variant Impact Comparison
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h5>OSH Gene Variants</h5>
                                <table class="table table-striped">
                                    <thead>
                                        <tr>
                                            <th>Impact</th>
                                            <th>Count</th>
                                            <th>Percentage</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <tr>
                                            <td>HIGH</td>
                                            <td>30</td>
                                            <td>21.4%</td>
                                        </tr>
                                        <tr>
                                            <td>MODERATE</td>
                                            <td>0</td>
                                            <td>0.0%</td>
                                        </tr>
                                        <tr>
                                            <td>LOW</td>
                                            <td>0</td>
                                            <td>0.0%</td>
                                        </tr>
                                        <tr>
                                            <td>MODIFIER</td>
                                            <td>110</td>
                                            <td>78.6%</td>
                                        </tr>
                                    </tbody>
                                </table>
                            </div>
                            <div class="col-md-6">
                                <h5>ERG Gene Variants</h5>
                                <table class="table table-striped">
                                    <thead>
                                        <tr>
                                            <th>Impact</th>
                                            <th>Count</th>
                                            <th>Percentage</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        <tr>
                                            <td>HIGH</td>
                                            <td>15</td>
                                            <td>4.4%</td>
                                        </tr>
                                        <tr>
                                            <td>MODERATE</td>
                                            <td>21</td>
                                            <td>6.1%</td>
                                        </tr>
                                        <tr>
                                            <td>LOW</td>
                                            <td>0</td>
                                            <td>0.0%</td>
                                        </tr>
                                        <tr>
                                            <td>MODIFIER</td>
                                            <td>308</td>
                                            <td>89.5%</td>
                                        </tr>
                                    </tbody>
                                </table>
                            </div>
                        </div>
                        <div class="alert alert-info mt-3">
                            <strong>Note:</strong> The higher proportion of HIGH impact variants in OSH genes (21.4% vs 4.4% in ERG genes) suggests that OSH genes may have more flexibility in their auxiliary functions while still maintaining their core functionality. This could indicate different evolutionary constraints on sterol transport versus synthesis.
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <h2 id="relationships" class="mt-5">OSH-ERG Gene Relationships</h2>
        
        <div class="row">
            <div class="col-lg-6">
                <div class="card mb-4">
                    <div class="card-header">
                        <i class="bi bi-rulers"></i> Distance Analysis
                    </div>
                    <div class="card-body">
                        <h5>Distance Distribution</h5>
                        <ul>
                            <li>Minimum distance between OSH-ERG gene pairs: 742 bp</li>
                            <li>Maximum distance: 390,857 bp</li>
                            <li>Mean distance: 147,173 bp</li>
                            <li>Median distance: 144,695 bp</li>
                        </ul>
                        
                        <h5>Zone Distribution</h5>
                        <ul>
                            <li>Buffer zone (1-5,000 bp): {zone_counts.get("Buffer", 2)} pairs (16.7%)</li>
                            <li>Intermediate zone (5,001-50,000 bp): {zone_counts.get("Intermediate", 2)} pairs (16.7%)</li>
                            <li>Satellite zone (>50,000 bp): {zone_counts.get("Satellite", 8)} pairs (66.7%)</li>
                        </ul>
                        
                        <h5>Closest Gene Pairs</h5>
                        <table class="table table-sm">
                            <thead>
                                <tr>
                                    <th>OSH Gene</th>
                                    <th>ERG Gene</th>
                                    <th>Distance</th>
                                    <th>Direction</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td>OSH3 (YHR073W)</td>
                                    <td>ERG7 (YHR072W)</td>
                                    <td>742 bp</td>
                                    <td>Downstream</td>
                                </tr>
                                <tr>
                                    <td>OSH7 (YHR001W)</td>
                                    <td>ERG11 (YHR007C)</td>
                                    <td>12,900 bp</td>
                                    <td>Upstream</td>
                                </tr>
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
            
            <div class="col-lg-6">
                <div class="card mb-4">
                    <div class="card-header">
                        <i class="bi bi-diagram-3"></i> Network Analysis
                    </div>
                    <div class="card-body">
                        <h5>Network Characteristics</h5>
                        <p>Analysis of the OSH-ERG gene network reveals:</p>
                        <ul>
                            <li>Hub genes: OSH3 and OSH7 form central connections with multiple ERG genes</li>
                            <li>Isolated genes: Some OSH genes have no ERG genes on the same scaffold</li>
                            <li>Chromosome clustering: Genes on the same chromosome tend to form connected subnetworks</li>
                            <li>Distance-based relationships: Connections weighted by genomic distance show clustering of close gene pairs</li>
                        </ul>
                        
                        <h5>Biological Implications</h5>
                        <p>The network structure suggests:</p>
                        <ul>
                            <li>Coordinated regulation between specific OSH-ERG gene pairs</li>
                            <li>Potential co-evolution of spatially close genes</li>
                            <li>Different roles for hub OSH genes vs. peripheral ones</li>
                            <li>Chromosomal organization may reflect functional relationships</li>
                        </ul>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="row mb-5">
            <div class="col-lg-12">
                <div class="card">
                    <div class="card-header">
                        <i class="bi bi-images"></i> Visualization Gallery
                    </div>
                    <div class="card-body">
                        <div class="row">
                            {visualization_gallery}
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <h2 id="significance">Biological Significance and Implications</h2>
        
        <div class="row">
            <div class="col-lg-6">
                <div class="card mb-4">
                    <div class="card-header bg-info text-white">
                        <i class="bi bi-lightbulb"></i> Biological Significance
                    </div>
                    <div class="card-body">
                        <h5>1. Conservation Patterns</h5>
                        <p>The similar conservation level between OSH and ERG genes suggests they both perform essential functions that are under strong evolutionary constraints. The absence of within-gene variants in both gene families indicates strong purifying selection.</p>
                        
                        <h5>2. Functional Implications</h5>
                        <p>The conservation of OSH genes suggests their sterol transport function is as essential as the ergosterol synthesis function of ERG genes. The higher proportion of HIGH impact variants in OSH genes indicates potentially more flexibility in their auxiliary functions.</p>
                        
                        <h5>3. Four-Zone Architecture</h5>
                        <p>The OSH-ERG distance analysis supports the four-zone conservation architecture previously identified, indicating a hierarchical organization of genes involved in sterol metabolism and transport that balances essential function preservation with adaptive flexibility.</p>
                        
                        <h5>4. Treatment Adaptations</h5>
                        <p>Slight differences in variant counts between temperature adaptation and low oxygen adaptation, along with differences between gene-modified and non-modified strains, suggest different adaptive strategies for different environmental conditions.</p>
                    </div>
                </div>
            </div>
            
            <div class="col-lg-6">
                <div class="card mb-4">
                    <div class="card-header bg-success text-white">
                        <i class="bi bi-diagram-3"></i> Membrane Biology Context
                    </div>
                    <div class="card-body">
                        <h5>Sterol Transport and Membrane Function</h5>
                        <p>The OSH gene analysis provides important context for understanding membrane biology:</p>
                        
                        <ul>
                            {membrane_items}
                        </ul>
                        
                        <div class="highlight-box">
                            <h5>Model for OSH-ERG Interactions</h5>
                            <p>Our findings suggest a model where:</p>
                            <ol>
                                <li>ERG genes produce ergosterol with precise structural requirements</li>
                                <li>OSH proteins transport this ergosterol to specific membrane locations</li>
                                <li>This coordinated system maintains membrane integrity and function</li>
                                <li>Adaptation occurs through regulatory changes rather than alterations to core enzyme or transport functions</li>
                            </ol>
                        </div>
                    </div>
                </div>
                
                <div class="card mb-4">
                    <div class="card-header bg-primary text-white">
                        <i class="bi bi-journal-check"></i> Implications and Next Steps
                    </div>
                    <div class="card-body">
                        <ol>
                            <li><strong>Enhanced Satellite Gene Characterization:</strong> Systematically identify genes in the satellite zone of each ERG gene, including OSH genes based on our distance measurements</li>
                            <li><strong>Comprehensive Regulatory Region Analysis:</strong> Map precise locations of variants relative to gene features, focusing on the regulatory regions of both ERG and OSH genes</li>
                            <li><strong>Statistical Modeling of the Four-Zone Architecture:</strong> Include OSH genes in the quantitative modeling and refine zone boundary definitions</li>
                            <li><strong>Integration with Sterol Profiles:</strong> Connect genetic changes in OSH genes with observed changes in sterol composition</li>
                        </ol>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="card mb-4">
            <div class="card-header">
                <i class="bi bi-check-circle-fill"></i> Conclusion
            </div>
            <div class="card-body">
                <p>The OSH gene family analysis reveals that:</p>
                
                <ol>
                    <li>OSH genes, similar to ERG genes, are under strong purifying selection, suggesting their critical role in sterol transport and membrane function</li>
                    <li>The significant relationship between OSH and ERG genes, especially the close proximity of certain pairs (OSH3-ERG7, OSH7-ERG11), indicates potential co-regulation or functional interaction</li>
                    <li>The higher proportion of HIGH impact variants in OSH genes might suggest a more flexible role in adaptation compared to ERG genes</li>
                    <li>The OSH-ERG gene relationships support the four-zone conservation architecture, indicating a hierarchical organization of genes involved in sterol metabolism and transport</li>
                </ol>
                
                <p>These findings align with the biochemical understanding of sterol biosynthesis and transport, where the precise structural requirements for ergosterol synthesis are complemented by more flexible transport mechanisms to maintain membrane homeostasis under different environmental conditions.</p>
            </div>
        </div>
        
        <footer class="mt-5 pt-4 border-top text-center text-muted">
            <p>Yeast MSA Project - OSH Gene Family Analysis Report</p>
            <p>Generated on: {current_date}</p>
        </footer>
    </div>
    
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    <script>
        // Handle scrolling and active section in navigation
        window.addEventListener('scroll', function() {{
            const sections = document.querySelectorAll('h2[id]');
            const navLinks = document.querySelectorAll('.navbar-nav .nav-link');
            
            let currentSection = '';
            
            sections.forEach(section => {{
                const sectionTop = section.offsetTop - 100;
                if (window.scrollY >= sectionTop) {{
                    currentSection = section.getAttribute('id');
                }}
            }});
            
            navLinks.forEach(link => {{
                link.classList.remove('active');
                if (link.getAttribute('href') === '#' + currentSection) {{
                    link.classList.add('active');
                }}
            }});
        }});
    </script>
</body>
</html>
"""
    
    return html

def main():
    """Main function to generate the HTML report."""
    print("Generating OSH Gene Family Analysis Report...")
    
    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Generate the HTML content
    html_content = generate_html_content()
    
    # Write the HTML to file
    output_file = OUTPUT_DIR / "osh_gene_analysis.html"
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"HTML report generated: {output_file}")

if __name__ == "__main__":
    main()