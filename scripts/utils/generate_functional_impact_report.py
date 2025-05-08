#!/usr/bin/env python3

import os
import base64
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from pathlib import Path
import re
import json
from datetime import datetime

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
        if sep == '\t' or file_path.endswith('.tsv'):
            return pd.read_csv(file_path, sep='\t')
        return pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return pd.DataFrame()

def image_to_base64(image_path):
    """Convert an image file to base64."""
    try:
        with open(image_path, 'rb') as image_file:
            return base64.b64encode(image_file.read()).decode('utf-8')
    except Exception as e:
        print(f"Error encoding image {image_path}: {e}")
        return ""

def collect_data():
    """Collect all relevant data from the Yeast MSA project."""
    base_dir = Path('/Users/zakiralibhai/Documents/GitHub/Yeast_MSA')
    
    # Function to collect all PNG files from a directory and its subdirectories
    def get_png_files(directory):
        png_dict = {}
        try:
            for path in directory.glob('**/*.png'):
                # Create a relative path from the directory for the key
                rel_path = path.relative_to(directory)
                # Use directory name + file stem as the key
                key = f"{path.parent.name}_{path.stem}" if path.parent.name != directory.name else path.stem
                png_dict[key] = str(path)
        except Exception as e:
            print(f"Error collecting PNG files from {directory}: {e}")
        return png_dict
    
    # Create a dictionary to store all data
    data = {
        "reports": {},
        "tables": {},
        "images": {},
        "summary_stats": {},
    }
    
    # Collect report files
    report_files = [
        "results/scaffold_variants/scaffold_variant_summary.txt",
        "results/network_analysis/network_analysis_report.md",
        "results/functional_impact/high_impact/high_impact_variants_report.md",
        "results/functional_impact/variants_by_distance/variants_by_distance_report.md",
        "results/functional_impact/key_genomic_regions/key_regions_report.md",
        "analysis/gene_mutation_spectrum_results/gene_analysis_report.txt",
        "results/treatment_analysis/analysis_summary.txt",
        "analysis/genomic_context_results/genomic_context_summary.txt",
        "analysis/mutational_signatures_results/mutational_signatures_summary.txt",
        "analysis/regional_enrichment_results/regional_enrichment_summary.txt"
    ]
    
    for file_path in report_files:
        full_path = base_dir / file_path
        if full_path.exists():
            report_name = full_path.stem
            data["reports"][report_name] = read_file(full_path)
    
    # Collect TSV and CSV table files
    key_table_files = [
        "results/functional_impact/high_impact/high_impact_variants.tsv",
        "results/functional_impact/variants_by_distance/high_impact_variants_by_distance.tsv",
        "results/gene_variants_expanded/gene_summary.tsv",
        "results/gene_variants_expanded/impact_summary.tsv",
        "results/gene_variants_expanded/effect_summary.tsv",
        "results/network_analysis/network_edges.tsv",
        "results/network_analysis/network_nodes.tsv",
        "analysis/gene_mutation_spectrum_results/gene_mutation_spectrum_summary.csv"
    ]
    
    for file_path in key_table_files:
        full_path = base_dir / file_path
        if full_path.exists():
            table_name = full_path.stem
            sep = '\t' if file_path.endswith('.tsv') else ','
            data["tables"][table_name] = read_csv(full_path, sep=sep)
    
    # Collect images from results directories
    key_directories = [
        "results/functional_impact/variants_by_distance",
        "results/functional_impact/key_genomic_regions/region_visualizations",
        "results/network_analysis",
        "results/network_analysis/erg_subnetworks",
        "results/treatment_analysis",
        "results/scaffold_variants/visualizations",
        "analysis/gene_mutation_spectrum_results",
        "analysis/genomic_context_results",
        "analysis/mutational_signatures_results",
        "analysis/population_structure_results",
        "analysis/regional_enrichment_results",
        "analysis/scaffold_distribution_results"
    ]
    
    for directory in key_directories:
        dir_path = base_dir / directory
        if dir_path.exists():
            dir_name = dir_path.name
            if dir_name not in data["images"]:
                data["images"][dir_name] = {}
            
            # Add all PNG files from this directory
            png_files = get_png_files(dir_path)
            data["images"][dir_name].update(png_files)
    
    # Process data to extract key statistics and findings
    data["summary_stats"] = extract_key_stats(data)
    
    return data

def extract_key_stats(data):
    """Extract key statistics and findings from the data."""
    stats = {
        "variant_counts": {},
        "pathway_genes": [],
        "satellite_genes": [],
        "key_findings": [],
        "conservation_zones": {},
    }
    
    # Extract variant counts by treatment from treatment_analysis
    if "analysis_summary" in data["reports"]:
        treatment_report = data["reports"]["analysis_summary"]
        # Extract variant counts
        variant_counts_match = re.search(r'1\. Variant Counts by Treatment(.*?)2\.', treatment_report, re.DOTALL)
        if variant_counts_match:
            variant_counts_text = variant_counts_match.group(1)
            for line in variant_counts_text.strip().split('\n'):
                if ':' in line:
                    treatment, count = line.split(':', 1)
                    count = int(count.strip().split()[0])
                    stats["variant_counts"][treatment.strip()] = count
    
    # Extract ergosterol pathway genes and satellite genes
    if "network_analysis_report" in data["reports"]:
        network_report = data["reports"]["network_analysis_report"]
        
        # Extract ergosterol pathway genes
        erg_match = re.search(r'Ergosterol pathway genes: (\d+)', network_report)
        if erg_match:
            # Get the list of ERG genes
            erg_genes = ["ERG1", "ERG2", "ERG3", "ERG4", "ERG5", "ERG6", "ERG7", "ERG9", "ERG11", "ERG24", "ERG25"]
            stats["pathway_genes"] = erg_genes
        
        # Extract satellite genes
        satellite_section = re.search(r'Genes with HIGH Impact Variants(.*?)Treatment-Specific Patterns', network_report, re.DOTALL)
        if satellite_section:
            satellite_text = satellite_section.group(1)
            satellite_pattern = r'- (W\d+) \((W\d+)\): \d+ HIGH impact variants, Near: (ERG\d+), Distance: (\d+) bp'
            satellite_matches = re.findall(satellite_pattern, satellite_text)
            for match in satellite_matches:
                stats["satellite_genes"].append({
                    "gene_id": match[0],
                    "near_gene": match[2],
                    "distance": int(match[3])
                })
    
    # Extract conservation zones from the key regions report
    if "key_regions_report" in data["reports"]:
        regions_report = data["reports"]["key_regions_report"]
        
        # Extract region information
        region_sections = re.findall(r'### \d+\.\d+\. ([^\n]+)[\s\S]*?- ERG Gene: ([^\n]+)[\s\S]*?- Scaffold: ([^\n]+)[\s\S]*?- Region Size: ([^\n]+)[\s\S]*?- Variants: (\d+) \(Expected: ([^\)]+)\)[\s\S]*?- Enrichment: ([^,]+), p-value: ([^\n]+)', regions_report)
        
        for section in region_sections:
            region_name = section[0]
            stats["conservation_zones"][region_name] = {
                "erg_gene": section[1].strip(),
                "scaffold": section[2].strip(),
                "region_size": section[3].strip(),
                "variants": int(section[4]),
                "expected": float(section[5]),
                "enrichment": float(section[6].strip('x')),
                "p_value": float(section[7])
            }
    
    # Extract key findings from the high impact variants report
    if "high_impact_variants_report" in data["reports"]:
        hi_report = data["reports"]["high_impact_variants_report"]
        key_findings_section = re.search(r'Key implications:(.*?)$', hi_report, re.DOTALL)
        if key_findings_section:
            findings_text = key_findings_section.group(1)
            findings = re.findall(r'\d+\.\s+(.*?)$', findings_text, re.MULTILINE)
            stats["key_findings"].extend(findings)
    
    # Extract biological interpretation from the variants by distance report
    if "variants_by_distance_report" in data["reports"]:
        vbd_report = data["reports"]["variants_by_distance_report"]
        bio_interp_section = re.search(r'## Biological Interpretation(.*?)$', vbd_report, re.DOTALL)
        if bio_interp_section:
            interp_text = bio_interp_section.group(1)
            findings = []
            
            # Extract subsections
            subsections = re.findall(r'### ([^\n]+)(.*?)(?=###|\Z)', interp_text, re.DOTALL)
            for title, content in subsections:
                # Extract bullet points
                bullets = re.findall(r'-(.*?)$', content, re.MULTILINE)
                clean_bullets = [bullet.strip() for bullet in bullets if bullet.strip()]
                if clean_bullets:
                    findings.append(f"{title.strip()}: {clean_bullets[0]}")
            
            stats["key_findings"].extend(findings)
    
    return stats

def generate_html_report(data):
    """Generate an HTML report with all functional impact data."""
    
    # Extract some key stats for the overview
    summary_stats = data["summary_stats"]
    variant_counts = summary_stats.get("variant_counts", {})
    pathway_genes = summary_stats.get("pathway_genes", [])
    satellite_genes = summary_stats.get("satellite_genes", [])
    key_findings = summary_stats.get("key_findings", [])
    
    # Begin HTML content
    html = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Yeast MSA Project - Functional Impact Analysis</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.min.css">
        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/dataTables.bootstrap5.min.css">
        <style>
            body {
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                line-height: 1.6;
                color: #333;
                padding-top: 60px;
                background-color: #f8f9fa;
            }
            
            .container {
                max-width: 1200px;
                background-color: white;
                box-shadow: 0 0 20px rgba(0,0,0,0.05);
                padding: 30px;
                margin-bottom: 50px;
                border-radius: 8px;
            }
            
            h1, h2, h3, h4 {
                color: #0d6efd;
                margin-top: 1.5em;
                margin-bottom: 0.8em;
                font-weight: 600;
            }
            
            h1 {
                font-size: 2.5rem;
                text-align: center;
                margin-bottom: 1.5em;
                color: #0a4fa0;
                border-bottom: 2px solid #dee2e6;
                padding-bottom: 15px;
            }
            
            h2 {
                font-size: 1.8rem;
                border-bottom: 1px solid #eaecef;
                padding-bottom: 10px;
            }
            
            h3 {
                font-size: 1.4rem;
            }
            
            h4 {
                font-size: 1.2rem;
                color: #495057;
            }
            
            img {
                max-width: 100%;
                height: auto;
                margin: 20px 0;
                border-radius: 5px;
                box-shadow: 0 3px 10px rgba(0,0,0,0.1);
            }
            
            .nav-pills .nav-link.active {
                background-color: #0d6efd;
            }
            
            .nav-pills .nav-link {
                color: #495057;
                border-radius: 5px;
                margin: 5px 0;
                padding: 10px 15px;
            }
            
            .nav-pills .nav-link:hover {
                background-color: #e9ecef;
            }
            
            pre {
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 5px;
                padding: 15px;
                margin: 15px 0;
                overflow-x: auto;
            }
            
            table {
                width: 100%;
                margin-bottom: 1rem;
                border-collapse: collapse;
            }
            
            th, td {
                padding: 0.75rem;
                vertical-align: top;
                border-top: 1px solid #dee2e6;
            }
            
            thead th {
                vertical-align: bottom;
                border-bottom: 2px solid #dee2e6;
                background-color: #f8f9fa;
            }
            
            .card {
                margin-bottom: 20px;
                border-radius: 5px;
                box-shadow: 0 3px 10px rgba(0,0,0,0.05);
            }
            
            .card-header {
                font-weight: 600;
                background-color: rgba(13, 110, 253, 0.1);
            }
            
            .img-gallery {
                display: flex;
                flex-wrap: wrap;
                gap: 15px;
                margin: 20px 0;
            }
            
            .img-gallery-item {
                flex: 1 0 300px;
                max-width: 500px;
                margin-bottom: 15px;
            }
            
            .img-gallery-item img {
                width: 100%;
                height: auto;
                cursor: pointer;
                transition: transform 0.3s ease;
            }
            
            .img-gallery-item img:hover {
                transform: scale(1.02);
            }
            
            .img-caption {
                font-size: 0.9rem;
                text-align: center;
                margin-top: 5px;
                color: #495057;
            }
            
            #imageModal .modal-dialog {
                max-width: 90%;
                max-height: 90vh;
            }
            
            .text-primary-emphasis {
                color: #084298 !important;
            }
            
            .bg-primary-subtle {
                background-color: #cfe2ff !important;
            }
            
            .text-secondary-emphasis {
                color: #41464b !important;
            }
            
            .bg-secondary-subtle {
                background-color: #e2e3e5 !important;
            }
            
            .border-primary-subtle {
                border-color: #9ec5fe !important;
            }
            
            .markdown-content {
                line-height: 1.8;
            }
            
            .markdown-content p {
                margin-bottom: 1.2em;
            }
            
            .markdown-content ul {
                margin-bottom: 1.2em;
            }
            
            .markdown-content h2 {
                margin-top: 2em;
            }
            
            .markdown-content h3 {
                margin-top: 1.5em;
            }
            
            .metadata {
                color: #6c757d;
                font-size: 0.9rem;
                margin-bottom: 25px;
            }
            
            .summary-card {
                background-color: #f8f9fa;
                border-left: 4px solid #0d6efd;
                padding: 15px;
                margin: 20px 0;
            }
            
            .highlight-box {
                background-color: rgba(13, 110, 253, 0.05);
                border-left: 4px solid #0d6efd;
                padding: 15px;
                margin: 20px 0;
                border-radius: 5px;
            }
            
            .figure-container {
                text-align: center;
                margin: 25px 0;
            }
            
            .figure-caption {
                font-size: 0.9rem;
                color: #6c757d;
                margin-top: 10px;
                max-width: 80%;
                margin-left: auto;
                margin-right: auto;
            }
            
            .alert-compact {
                padding: 0.5rem 1rem;
                margin-bottom: 1rem;
            }
            
            .chart-container {
                height: 400px;
                width: 100%;
                margin: 20px 0;
            }
            
            .tooltip-inner {
                max-width: 300px;
            }
            
            .badge-lg {
                font-size: 0.9rem;
                padding: 0.4em 0.8em;
            }
            
            .loading-spinner {
                text-align: center;
                padding: 40px;
            }
            
            /* Navigation styles */
            .sidebar {
                position: sticky;
                top: 80px;
                height: calc(100vh - 100px);
                overflow-y: auto;
            }
            
            #toc {
                padding-right: 15px;
            }
            
            .section-nav a {
                display: block;
                padding: 0.3rem 0;
                color: #495057;
                text-decoration: none;
            }
            
            .section-nav a:hover {
                color: #0d6efd;
            }
            
            .section-nav a.active {
                color: #0d6efd;
                font-weight: 600;
            }
            
            .section-nav ul {
                padding-left: 1rem;
                list-style: none;
            }
            
            /* DataTables overrides */
            .dataTables_wrapper .dataTables_length, 
            .dataTables_wrapper .dataTables_filter {
                margin-bottom: 15px;
            }
            
            /* Dark mode */
            @media (prefers-color-scheme: dark) {
                .dark-mode-enabled body {
                    background-color: #212529;
                    color: #e9ecef;
                }
                
                .dark-mode-enabled .container {
                    background-color: #2c3034;
                    box-shadow: 0 0 20px rgba(0,0,0,0.3);
                }
                
                .dark-mode-enabled h1, 
                .dark-mode-enabled h2, 
                .dark-mode-enabled h3, 
                .dark-mode-enabled h4 {
                    color: #6ea8fe;
                }
                
                .dark-mode-enabled h1 {
                    border-bottom-color: #495057;
                }
                
                .dark-mode-enabled h2 {
                    border-bottom-color: #495057;
                }
                
                .dark-mode-enabled .nav-pills .nav-link {
                    color: #e9ecef;
                }
                
                .dark-mode-enabled .nav-pills .nav-link:hover {
                    background-color: #343a40;
                }
                
                .dark-mode-enabled pre {
                    background-color: #343a40;
                    border-color: #495057;
                }
                
                .dark-mode-enabled th, 
                .dark-mode-enabled td {
                    border-top-color: #495057;
                }
                
                .dark-mode-enabled thead th {
                    border-bottom-color: #495057;
                    background-color: #343a40;
                }
                
                .dark-mode-enabled .card {
                    background-color: #343a40;
                    border-color: #495057;
                }
                
                .dark-mode-enabled .card-header {
                    background-color: rgba(13, 110, 253, 0.2);
                }
                
                .dark-mode-enabled .text-dark {
                    color: #e9ecef !important;
                }
                
                .dark-mode-enabled .bg-light {
                    background-color: #343a40 !important;
                }
                
                .dark-mode-enabled .border-light {
                    border-color: #495057 !important;
                }
                
                .dark-mode-enabled .img-caption {
                    color: #adb5bd;
                }
                
                .dark-mode-enabled .metadata {
                    color: #adb5bd;
                }
                
                .dark-mode-enabled .summary-card {
                    background-color: #343a40;
                }
                
                .dark-mode-enabled .highlight-box {
                    background-color: rgba(13, 110, 253, 0.1);
                }
                
                .dark-mode-enabled .figure-caption {
                    color: #adb5bd;
                }
                
                .dark-mode-enabled .section-nav a {
                    color: #adb5bd;
                }
                
                .dark-mode-enabled .section-nav a:hover,
                .dark-mode-enabled .section-nav a.active {
                    color: #6ea8fe;
                }
                
                .dark-mode-enabled .dropdown-menu {
                    background-color: #343a40;
                    border-color: #495057;
                }
                
                .dark-mode-enabled .dropdown-item {
                    color: #e9ecef;
                }
                
                .dark-mode-enabled .dropdown-item:hover {
                    background-color: #495057;
                    color: #fff;
                }
                
                .dark-mode-enabled .modal-content {
                    background-color: #343a40;
                    border-color: #495057;
                }
                
                .dark-mode-enabled .modal-header,
                .dark-mode-enabled .modal-footer {
                    border-color: #495057;
                }
                
                .dark-mode-enabled .close {
                    color: #e9ecef;
                }
                
                .dark-mode-enabled .text-primary-emphasis {
                    color: #6ea8fe !important;
                }
                
                .dark-mode-enabled .bg-primary-subtle {
                    background-color: rgba(13, 110, 253, 0.2) !important;
                }
                
                .dark-mode-enabled .text-secondary-emphasis {
                    color: #adb5bd !important;
                }
                
                .dark-mode-enabled .bg-secondary-subtle {
                    background-color: #495057 !important;
                }
                
                .dark-mode-enabled .border-primary-subtle {
                    border-color: rgba(13, 110, 253, 0.3) !important;
                }
            }
            
            /* Print styles */
            @media print {
                body {
                    padding-top: 0;
                    font-size: 12pt;
                }
                
                .container {
                    box-shadow: none;
                    max-width: 100%;
                    padding: 0;
                }
                
                nav, .sidebar, .no-print {
                    display: none !important;
                }
                
                h1, h2, h3, h4 {
                    page-break-after: avoid;
                }
                
                p, table, figure {
                    page-break-inside: avoid;
                }
                
                img {
                    max-width: 600px;
                }
            }
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
                    <ul class="navbar-nav me-auto">
                        <li class="nav-item">
                            <a class="nav-link" href="#overview">Overview</a>
                        </li>
                        <li class="nav-item">
                            <a class="nav-link" href="#conservation">Conservation Analysis</a>
                        </li>
                        <li class="nav-item">
                            <a class="nav-link" href="#network">Network Analysis</a>
                        </li>
                        <li class="nav-item">
                            <a class="nav-link" href="#genetic-analysis">Genetic Analysis</a>
                        </li>
                        <li class="nav-item">
                            <a class="nav-link" href="#visualizations">Visualizations</a>
                        </li>
                        <li class="nav-item">
                            <a class="nav-link" href="#conclusions">Conclusions</a>
                        </li>
                    </ul>
                    <ul class="navbar-nav">
                        <li class="nav-item">
                            <a class="nav-link" href="#" id="toggleDarkMode">
                                <i class="bi bi-moon"></i> Toggle Dark Mode
                            </a>
                        </li>
                    </ul>
                </div>
            </div>
        </nav>

        <div class="container mt-5">
    """
    
    # Overview Section
    html += """
            <h1 id="overview">Functional Impact Analysis Report</h1>
            <div class="metadata text-center mb-4">
                <span class="me-3"><i class="bi bi-calendar"></i> Generated on: """ + datetime.now().strftime("%B %d, %Y") + """</span>
                <span class="me-3"><i class="bi bi-tag"></i> Yeast MSA Project</span>
            </div>
            
            <div class="row">
                <div class="col-lg-8">
                    <div class="alert alert-primary">
                        <h4><i class="bi bi-info-circle"></i> Project Context</h4>
                        <p>This report presents a comprehensive analysis of the functional impact of genetic variants in the Yeast MSA project. The analysis reveals a sophisticated hierarchical conservation pattern around the ergosterol pathway genes, with important implications for understanding yeast adaptation mechanisms.</p>
                    </div>
                    
                    <h2>Project Overview</h2>
                    <p>The Yeast Multiple Sequence Alignment (MSA) project investigates how yeast (S. cerevisiae, W303 strain) adapts to different environmental stresses through genetic mutations, focusing on:</p>
                    <ul>
                        <li><strong>Temperature adaptation</strong> (WT-37 and CAS strains)</li>
                        <li><strong>Low oxygen adaptation</strong> (WTA and STC strains)</li>
                        <li><strong>Gene modifications</strong> (CAS and STC strains)</li>
                    </ul>
                    
                    <p>A key focus of this analysis is the ergosterol biosynthetic pathway, which is critical for cell membrane integrity and function. The analysis examines how this essential pathway can be conserved yet allow for adaptation to different environmental stressors.</p>
                    
                    <h3>Key Analysis Components</h3>
                    <div class="row">
                        <div class="col-md-6">
                            <div class="card h-100">
                                <div class="card-header bg-primary text-white">
                                    <i class="bi bi-shield-check"></i> Conservation Analysis
                                </div>
                                <div class="card-body">
                                    <p>Examination of purifying selection on ergosterol pathway genes and the hierarchical conservation gradient extending from these genes.</p>
                                </div>
                            </div>
                        </div>
                        <div class="col-md-6">
                            <div class="card h-100">
                                <div class="card-header bg-success text-white">
                                    <i class="bi bi-diagram-3"></i> Network Analysis
                                </div>
                                <div class="card-body">
                                    <p>Analysis of the extended ergosterol network including pathway genes and affected satellite genes at consistent distances.</p>
                                </div>
                            </div>
                        </div>
                    </div>
                    
                    <div class="row mt-3">
                        <div class="col-md-6">
                            <div class="card h-100">
                                <div class="card-header bg-info text-white">
                                    <i class="bi bi-bar-chart-steps"></i> Variant Impact Analysis
                                </div>
                                <div class="card-body">
                                    <p>Characterization of HIGH and MODERATE impact variants and their relationship to ergosterol pathway genes.</p>
                                </div>
                            </div>
                        </div>
                        <div class="col-md-6">
                            <div class="card h-100">
                                <div class="card-header bg-warning text-dark">
                                    <i class="bi bi-geo-alt"></i> Variant Distribution
                                </div>
                                <div class="card-body">
                                    <p>Analysis of variant distribution across the genome, with a focus on key genomic regions showing significant enrichment.</p>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                
                <div class="col-lg-4">
                    <div class="card border-primary mb-4">
                        <div class="card-header bg-primary text-white">
                            <i class="bi bi-lightbulb"></i> Key Findings
                        </div>
                        <div class="card-body">
                            <ul class="list-group list-group-flush">
    """
    
    # Add variant counts to key findings
    for treatment, count in variant_counts.items():
        html += f"""
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    {treatment} Variants
                                    <span class="badge bg-primary rounded-pill">{count}</span>
                                </li>
        """
    
    # Add number of pathway genes
    html += f"""
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    Ergosterol Pathway Genes
                                    <span class="badge bg-success rounded-pill">{len(pathway_genes)}</span>
                                </li>
    """
    
    # Add number of satellite genes
    html += f"""
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    Satellite Genes with HIGH Impact
                                    <span class="badge bg-danger rounded-pill">{len(satellite_genes)}</span>
                                </li>
    """
    
    html += """
                            </ul>
                        </div>
                    </div>
                    
                    <div class="card border-success mb-4">
                        <div class="card-header bg-success text-white">
                            <i class="bi bi-clipboard-data"></i> Major Insights
                        </div>
                        <div class="card-body">
    """
    
    # Add key findings as alerts
    for i, finding in enumerate(key_findings[:4]):  # Limit to first 4 findings
        html += f"""
                            <div class="alert alert-light">
                                <p class="mb-1">{finding}</p>
                            </div>
        """
    
    html += """
                        </div>
                    </div>
                </div>
            </div>
    """
    
    # Conservation Analysis Section
    html += """
            <h2 id="conservation">Conservation Analysis</h2>
            <p>One of the most striking findings of our analysis is the strong conservation of ergosterol pathway genes despite adaptation to different environmental stressors.</p>
            
            <div class="card mb-4">
                <div class="card-header">
                    <i class="bi bi-shield-check"></i> Hierarchical Conservation Model
                </div>
                <div class="card-body">
                    <div class="row">
                        <div class="col-lg-6">
                            <h4>Conservation Gradient</h4>
                            <p>Our analysis reveals a hierarchical conservation pattern extending outward from ergosterol pathway genes:</p>
                            <ol>
                                <li><strong>Core Zone (0bp):</strong> Ergosterol genes - Complete conservation with no HIGH/MODERATE impact variants</li>
                                <li><strong>Buffer Zone (0-7kb):</strong> Strong conservation with no variants</li>
                                <li><strong>Satellite Zone (7-50kb):</strong> Specific genes with HIGH/MODERATE impact variants at consistent distances</li>
                                <li><strong>Distant Zone (>50kb):</strong> Less constrained, with ~75% of all variants</li>
                            </ol>
                            <p>This architecture balances essential function preservation with adaptive flexibility.</p>
                        </div>
                        
                        <div class="col-lg-6">
                            <div class="figure-container">
    """
    
    # Add conservation zone image if available
    conservation_image_path = None
    for dir_name, images in data["images"].items():
        for img_name, img_path in images.items():
            if "conservation" in img_name.lower() or "enrichment" in img_name.lower():
                conservation_image_path = img_path
                break
        if conservation_image_path:
            break
    
    if conservation_image_path:
        img_base64 = image_to_base64(conservation_image_path)
        html += f"""
                                <img src="data:image/png;base64,{img_base64}" alt="Conservation Zone Model" class="img-fluid">
                                <div class="figure-caption">The hierarchical conservation model showing how variants are distributed around ergosterol pathway genes.</div>
        """
    else:
        html += """
                                <div class="alert alert-warning">
                                    <i class="bi bi-exclamation-triangle"></i> Conservation zone visualization not available.
                                </div>
        """
    
    html += """
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            
            <h3>Key Genomic Regions</h3>
    """
    
    # Add Key Genomic Regions
    conservation_zones = summary_stats.get("conservation_zones", {})
    if conservation_zones:
        html += """
            <div class="table-responsive">
                <table class="table table-striped">
                    <thead>
                        <tr>
                            <th>Region</th>
                            <th>ERG Gene</th>
                            <th>Variants</th>
                            <th>Expected</th>
                            <th>Enrichment</th>
                            <th>p-value</th>
                        </tr>
                    </thead>
                    <tbody>
        """
        
        for region, stats in conservation_zones.items():
            html += f"""
                        <tr>
                            <td>{region}</td>
                            <td>{stats['erg_gene']}</td>
                            <td>{stats['variants']}</td>
                            <td>{stats['expected']}</td>
                            <td>{stats['enrichment']}x</td>
                            <td>{stats['p_value']:.2e}</td>
                        </tr>
            """
        
        html += """
                    </tbody>
                </table>
            </div>
        """
    
    # Add Variant Distribution by Distance
    html += """
            <h3>Variant Distribution by Distance</h3>
            <p>The analysis of HIGH and MODERATE impact variants relative to ergosterol pathway genes reveals specific distance patterns that are consistent across samples.</p>
            
            <div class="row">
                <div class="col-md-6">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-rulers"></i> Distance Categories
                        </div>
                        <div class="card-body">
                            <ul>
                                <li><strong>0-5000 bp:</strong> No HIGH or MODERATE impact variants</li>
                                <li><strong>5000-10000 bp:</strong> ~20% of variants, primarily near ERG11 and ERG24</li>
                                <li><strong>10000-20000 bp:</strong> ~14% of variants, primarily near ERG25</li>
                                <li><strong>20000-50000 bp:</strong> ~66% of variants, distributed across multiple genes</li>
                            </ul>
                            <p>The complete absence of HIGH/MODERATE impact variants within 5kb of pathway genes provides strong evidence for purifying selection.</p>
                        </div>
                    </div>
                </div>
                
                <div class="col-md-6">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-bar-chart"></i> Key Distance Relationships
                        </div>
                        <div class="card-body">
                            <ul>
                                <li><strong>ERG11:</strong> HIGH impact variants at 8,149 bp upstream</li>
                                <li><strong>ERG7:</strong> HIGH impact variants at 47,676 bp downstream</li>
                                <li><strong>ERG25:</strong> MODERATE impact variants at 15,949 bp upstream and 40,586 bp downstream</li>
                                <li><strong>ERG3:</strong> MODERATE impact variants at 47,606 bp upstream</li>
                                <li><strong>ERG4:</strong> MODERATE impact variants at 26,130 bp upstream</li>
                            </ul>
                            <p>These consistent distances suggest specific functional or structural relationships between the satellite genes and the ergosterol pathway.</p>
                        </div>
                    </div>
                </div>
            </div>
    """
    
    # Add Distance Distribution Visualization
    html += """
            <div class="row">
                <div class="col-lg-12">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-graph-up"></i> Distance Distribution Visualizations
                        </div>
                        <div class="card-body">
                            <div class="img-gallery">
    """
    
    # Add distance-related images
    distance_images = []
    for dir_name, images in data["images"].items():
        for img_name, img_path in images.items():
            if "distance" in img_name.lower():
                distance_images.append((img_name, img_path))
    
    # Display up to 4 distance-related images
    for i, (img_name, img_path) in enumerate(distance_images[:4]):
        img_base64 = image_to_base64(img_path)
        # Create a more readable caption
        caption = img_name.replace("_", " ").title()
        html += f"""
                                <div class="img-gallery-item">
                                    <img src="data:image/png;base64,{img_base64}" alt="{caption}" data-bs-toggle="modal" data-bs-target="#imageModal" onclick="showImageInModal(this)">
                                    <div class="img-caption">{caption}</div>
                                </div>
        """
    
    html += """
                            </div>
                        </div>
                    </div>
                </div>
            </div>
    """
    
    # Network Analysis Section
    html += """
            <h2 id="network">Network Analysis</h2>
            <p>The network analysis provides a systems-level view of the relationship between ergosterol pathway genes and the affected satellite genes that harbor variants.</p>
            
            <div class="row">
                <div class="col-lg-7">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-diagram-3"></i> Extended Ergosterol Network
                        </div>
                        <div class="card-body">
                            <p>The extended ergosterol pathway network includes:</p>
                            <ul>
                                <li><strong>Ergosterol pathway genes:</strong> 11 core genes (ERG1-11, ERG24-25), all showing complete conservation</li>
                                <li><strong>Affected satellite genes:</strong> 6 genes harboring HIGH/MODERATE impact variants at specific distances</li>
                                <li><strong>Network connections:</strong> 22 total connections representing genomic proximity and potential functional relationships</li>
                            </ul>
                            <p>The network reveals how adaptation may occur through changes in the broader genomic neighborhood rather than direct modification of essential pathway genes.</p>
                        </div>
                    </div>
                </div>
                
                <div class="col-lg-5">
                    <div class="figure-container">
    """
    
    # Add extended network visualization if available
    network_image_path = None
    for dir_name, images in data["images"].items():
        for img_name, img_path in images.items():
            if "extended_erg_network" in img_name.lower() and "CAS" not in img_name and "STC" not in img_name and "WT" not in img_name:
                network_image_path = img_path
                break
        if network_image_path:
            break
    
    if network_image_path:
        img_base64 = image_to_base64(network_image_path)
        html += f"""
                        <img src="data:image/png;base64,{img_base64}" alt="Extended Ergosterol Network" class="img-fluid">
                        <div class="figure-caption">The extended ergosterol network showing connections between pathway genes and affected satellite genes.</div>
        """
    else:
        html += """
                        <div class="alert alert-warning">
                            <i class="bi bi-exclamation-triangle"></i> Network visualization not available.
                        </div>
        """
    
    html += """
                    </div>
                </div>
            </div>
            
            <h3>Satellite Genes with HIGH Impact Variants</h3>
            <p>The network analysis identified several satellite genes with HIGH impact variants that are consistently located at specific distances from ergosterol pathway genes.</p>
            
            <div class="table-responsive">
                <table class="table table-striped">
                    <thead>
                        <tr>
                            <th>Gene ID</th>
                            <th>Near ERG Gene</th>
                            <th>Distance (bp)</th>
                            <th>Impact</th>
                            <th>Variant Effect</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td>W3030H00610</td>
                            <td>ERG11</td>
                            <td>8,149 (upstream)</td>
                            <td>HIGH</td>
                            <td>Frameshift variant</td>
                        </tr>
                        <tr>
                            <td>W3030H01660</td>
                            <td>ERG7</td>
                            <td>47,676 (downstream)</td>
                            <td>HIGH</td>
                            <td>Frameshift variant</td>
                        </tr>
                        <tr>
                            <td>W3030G02910</td>
                            <td>ERG25</td>
                            <td>15,949 (upstream)</td>
                            <td>MODERATE</td>
                            <td>Missense variant (Arg340Trp)</td>
                        </tr>
                        <tr>
                            <td>W3030G03230</td>
                            <td>ERG25</td>
                            <td>40,586 (downstream)</td>
                            <td>MODERATE</td>
                            <td>Missense variant (Leu336Val)</td>
                        </tr>
                        <tr>
                            <td>W3030G02200</td>
                            <td>ERG4</td>
                            <td>26,130 (upstream)</td>
                            <td>MODERATE</td>
                            <td>Missense variant (Gly485Val)</td>
                        </tr>
                        <tr>
                            <td>W3030L01080</td>
                            <td>ERG3</td>
                            <td>47,606 (upstream)</td>
                            <td>MODERATE</td>
                            <td>Missense variant (Gly535Arg)</td>
                        </tr>
                    </tbody>
                </table>
            </div>
            
            <h3>Treatment-Specific Network Patterns</h3>
            <div class="row">
                <div class="col-lg-12">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-grid-3x3"></i> Network Visualization by Treatment
                        </div>
                        <div class="card-body">
                            <div class="img-gallery">
    """
    
    # Add treatment-specific network visualizations
    treatment_networks = []
    for dir_name, images in data["images"].items():
        for img_name, img_path in images.items():
            if "extended_erg_network_" in img_name.lower():
                treatment_networks.append((img_name, img_path))
    
    # Display treatment-specific network visualizations
    for img_name, img_path in treatment_networks:
        img_base64 = image_to_base64(img_path)
        # Create a more readable caption
        treatment = img_name.split("_")[-1].split(".")[0]
        caption = f"Extended Ergosterol Network - {treatment}"
        html += f"""
                                <div class="img-gallery-item">
                                    <img src="data:image/png;base64,{img_base64}" alt="{caption}" data-bs-toggle="modal" data-bs-target="#imageModal" onclick="showImageInModal(this)">
                                    <div class="img-caption">{caption}</div>
                                </div>
        """
    
    html += """
                            </div>
                        </div>
                    </div>
                </div>
            </div>
    """
    
    # Genetic Analysis Section
    html += """
            <h2 id="genetic-analysis">Genetic Analysis</h2>
            <p>The genetic analysis of variants provides insights into mutation patterns, effects, and their distribution across different treatment conditions.</p>
            
            <div class="row">
                <div class="col-lg-6">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-file-earmark-code"></i> Variant Effects and Impacts
                        </div>
                        <div class="card-body">
                            <h5>Distribution by Effect</h5>
                            <ul>
                                <li><strong>Upstream gene variant:</strong> 80.35% (1677)</li>
                                <li><strong>Missense variant:</strong> 6.61% (138)</li>
                                <li><strong>Frameshift variant:</strong> 6.47% (135)</li>
                                <li><strong>Synonymous variant:</strong> 2.92% (61)</li>
                                <li><strong>Downstream gene variant:</strong> 2.73% (57)</li>
                            </ul>
                            
                            <h5>Distribution by Impact</h5>
                            <ul>
                                <li><strong>MODIFIER:</strong> 83.09% (1734)</li>
                                <li><strong>MODERATE:</strong> 7.33% (153)</li>
                                <li><strong>HIGH:</strong> 6.66% (139)</li>
                                <li><strong>LOW:</strong> 2.92% (61)</li>
                            </ul>
                            <p>The predominance of upstream gene variants (80.35%) suggests that adaptation primarily occurs through changes in gene regulation rather than protein structure.</p>
                        </div>
                    </div>
                </div>
                
                <div class="col-lg-6">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-clipboard-data"></i> Treatment-Specific Patterns
                        </div>
                        <div class="card-body">
                            <h5>Variant Distribution by Treatment</h5>
                            <ul>
                                <li><strong>CAS:</strong> 437 (20.94%)</li>
                                <li><strong>STC:</strong> 422 (20.22%)</li>
                                <li><strong>WT-37:</strong> 416 (19.93%)</li>
                                <li><strong>WTA:</strong> 412 (19.74%)</li>
                                <li><strong>CAS-CTRL:</strong> 137 (6.56%)</li>
                                <li><strong>STC-CTRL:</strong> 133 (6.37%)</li>
                                <li><strong>WT-CTRL:</strong> 130 (6.23%)</li>
                            </ul>
                            
                            <h5>HIGH Impact Variant Pattern</h5>
                            <ul>
                                <li><strong>Gene-modified (CAS, STC):</strong> 16 variants each</li>
                                <li><strong>Non-modified (WT-37, WTA):</strong> 12 variants each</li>
                                <li><strong>Control (WT):</strong> 4 variants</li>
                            </ul>
                            <p>This creates a perfect 4:3:1 ratio maintained across all measurements, suggesting adaptation amplifies pre-existing genomic variation.</p>
                        </div>
                    </div>
                </div>
            </div>
            
            <h3>Mutation Spectrum Analysis</h3>
            <div class="row">
                <div class="col-lg-12">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-pie-chart"></i> Mutation Type Distribution
                        </div>
                        <div class="card-body">
                            <div class="row">
                                <div class="col-md-4">
                                    <h5>Variant Types</h5>
                                    <ul>
                                        <li><strong>INSERTION:</strong> 45.52% (950)</li>
                                        <li><strong>DELETION:</strong> 28.27% (590)</li>
                                        <li><strong>SNV:</strong> 26.21% (547)</li>
                                    </ul>
                                </div>
                                
                                <div class="col-md-8">
                                    <div class="figure-container">
    """
    
    # Add mutation spectrum visualization if available
    spectrum_image_path = None
    for dir_name, images in data["images"].items():
        for img_name, img_path in images.items():
            if any(term in img_name.lower() for term in ["mutation_spectrum", "mutation_by_adaptation", "comparative_mutation_spectrum"]):
                spectrum_image_path = img_path
                break
        if spectrum_image_path:
            break
    
    if spectrum_image_path:
        img_base64 = image_to_base64(spectrum_image_path)
        html += f"""
                                        <img src="data:image/png;base64,{img_base64}" alt="Mutation Spectrum" class="img-fluid">
                                        <div class="figure-caption">Comparative mutation spectrum across different treatment conditions.</div>
        """
    else:
        # Try to use mutation_by_gene_function.png if available
        gene_function_path = None
        for dir_name, images in data["images"].items():
            for img_name, img_path in images.items():
                if "mutation_by_gene_function" in img_name.lower():
                    gene_function_path = img_path
                    break
            if gene_function_path:
                break
        
        if gene_function_path:
            img_base64 = image_to_base64(gene_function_path)
            html += f"""
                                        <img src="data:image/png;base64,{img_base64}" alt="Mutation by Gene Function" class="img-fluid">
                                        <div class="figure-caption">Distribution of mutations by gene function and treatment.</div>
            """
        else:
            html += """
                                        <div class="alert alert-warning">
                                            <i class="bi bi-exclamation-triangle"></i> Mutation spectrum visualization not available.
                                        </div>
            """
    
    html += """
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
    """
    
    # Visualizations Section
    html += """
            <h2 id="visualizations">Visualizations</h2>
            <p>The following visualizations provide additional insights into the genetic variants, their distribution, and functional impact.</p>
            
            <div class="row mb-4">
                <div class="col-lg-3">
                    <div class="list-group" id="viz-tabs" role="tablist">
                        <a class="list-group-item list-group-item-action active" id="genomic-viz-tab" data-bs-toggle="list" href="#genomic-viz-content" role="tab">
                            Genomic Distribution
                        </a>
                        <a class="list-group-item list-group-item-action" id="impact-viz-tab" data-bs-toggle="list" href="#impact-viz-content" role="tab">
                            Functional Impact
                        </a>
                        <a class="list-group-item list-group-item-action" id="network-viz-tab" data-bs-toggle="list" href="#network-viz-content" role="tab">
                            Network Analysis
                        </a>
                        <a class="list-group-item list-group-item-action" id="treatment-viz-tab" data-bs-toggle="list" href="#treatment-viz-content" role="tab">
                            Treatment Comparisons
                        </a>
                    </div>
                </div>
                
                <div class="col-lg-9">
                    <div class="tab-content" id="viz-tabContent">
                        <div class="tab-pane fade show active" id="genomic-viz-content" role="tabpanel">
                            <h4>Genomic Distribution Visualizations</h4>
                            <p>These visualizations show how variants are distributed across the genome and their relationship to key genes.</p>
                            
                            <div class="img-gallery">
    """
    
    # Add genomic distribution visualizations
    genomic_images = []
    for dir_name, images in data["images"].items():
        for img_name, img_path in images.items():
            if any(term in img_name.lower() for term in ["scaffold", "distribution", "enrichment", "region"]):
                if "network" not in img_name.lower():  # Exclude network images
                    genomic_images.append((img_name, img_path))
    
    # Display up to 6 genomic distribution visualizations
    for img_name, img_path in genomic_images[:6]:
        img_base64 = image_to_base64(img_path)
        # Create a more readable caption
        caption = img_name.replace("_", " ").title()
        html += f"""
                                <div class="img-gallery-item">
                                    <img src="data:image/png;base64,{img_base64}" alt="{caption}" data-bs-toggle="modal" data-bs-target="#imageModal" onclick="showImageInModal(this)">
                                    <div class="img-caption">{caption}</div>
                                </div>
        """
    
    html += """
                            </div>
                        </div>
                        
                        <div class="tab-pane fade" id="impact-viz-content" role="tabpanel">
                            <h4>Functional Impact Visualizations</h4>
                            <p>These visualizations highlight the functional impact of variants and their effects on different gene categories.</p>
                            
                            <div class="img-gallery">
    """
    
    # Add functional impact visualizations
    impact_images = []
    for dir_name, images in data["images"].items():
        for img_name, img_path in images.items():
            if any(term in img_name.lower() for term in ["impact", "effect", "mutation", "signature"]):
                if "network" not in img_name.lower():  # Exclude network images
                    impact_images.append((img_name, img_path))
    
    # Display up to 6 functional impact visualizations
    for img_name, img_path in impact_images[:6]:
        img_base64 = image_to_base64(img_path)
        # Create a more readable caption
        caption = img_name.replace("_", " ").title()
        html += f"""
                                <div class="img-gallery-item">
                                    <img src="data:image/png;base64,{img_base64}" alt="{caption}" data-bs-toggle="modal" data-bs-target="#imageModal" onclick="showImageInModal(this)">
                                    <div class="img-caption">{caption}</div>
                                </div>
        """
    
    html += """
                            </div>
                        </div>
                        
                        <div class="tab-pane fade" id="network-viz-content" role="tabpanel">
                            <h4>Network Analysis Visualizations</h4>
                            <p>These visualizations show the extended ergosterol network and subnetworks centered on specific pathway genes.</p>
                            
                            <div class="img-gallery">
    """
    
    # Add network visualizations
    network_images = []
    for dir_name, images in data["images"].items():
        for img_name, img_path in images.items():
            if "network" in img_name.lower() or "subnetwork" in img_name.lower():
                network_images.append((img_name, img_path))
    
    # Display up to 6 network visualizations
    for img_name, img_path in network_images[:6]:
        img_base64 = image_to_base64(img_path)
        # Create a more readable caption
        caption = img_name.replace("_", " ").title()
        html += f"""
                                <div class="img-gallery-item">
                                    <img src="data:image/png;base64,{img_base64}" alt="{caption}" data-bs-toggle="modal" data-bs-target="#imageModal" onclick="showImageInModal(this)">
                                    <div class="img-caption">{caption}</div>
                                </div>
        """
    
    html += """
                            </div>
                        </div>
                        
                        <div class="tab-pane fade" id="treatment-viz-content" role="tabpanel">
                            <h4>Treatment Comparison Visualizations</h4>
                            <p>These visualizations compare variant patterns across different treatment conditions and adaptation types.</p>
                            
                            <div class="img-gallery">
    """
    
    # Add treatment comparison visualizations
    treatment_images = []
    for dir_name, images in data["images"].items():
        for img_name, img_path in images.items():
            if any(term in img_name.lower() for term in ["treatment", "adaptation", "comparison"]):
                treatment_images.append((img_name, img_path))
    
    # Display up to 6 treatment comparison visualizations
    for img_name, img_path in treatment_images[:6]:
        img_base64 = image_to_base64(img_path)
        # Create a more readable caption
        caption = img_name.replace("_", " ").title()
        html += f"""
                                <div class="img-gallery-item">
                                    <img src="data:image/png;base64,{img_base64}" alt="{caption}" data-bs-toggle="modal" data-bs-target="#imageModal" onclick="showImageInModal(this)">
                                    <div class="img-caption">{caption}</div>
                                </div>
        """
    
    html += """
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            
            <!-- Image Modal -->
            <div class="modal fade" id="imageModal" tabindex="-1" aria-hidden="true">
                <div class="modal-dialog modal-lg">
                    <div class="modal-content">
                        <div class="modal-header">
                            <h5 class="modal-title" id="imageModalLabel">Image Preview</h5>
                            <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                        </div>
                        <div class="modal-body text-center">
                            <img id="modalImage" src="" class="img-fluid" alt="Image Preview">
                        </div>
                    </div>
                </div>
            </div>
    """
    
    # Conclusions Section
    html += """
            <h2 id="conclusions">Conclusions</h2>
            
            <div class="alert alert-primary">
                <h4 class="alert-heading"><i class="bi bi-check-circle"></i> Key Finding</h4>
                <p>Our comprehensive analysis of genetic variants in the Yeast MSA project has revealed a sophisticated hierarchical conservation architecture surrounding the ergosterol biosynthetic pathway. This architecture balances essential function preservation with adaptive flexibility, allowing yeast to respond to environmental stresses while maintaining the integrity of critical cellular processes.</p>
            </div>
            
            <div class="row">
                <div class="col-lg-6">
                    <div class="card mb-4">
                        <div class="card-header bg-info text-white">
                            <i class="bi bi-lightbulb"></i> Key Biological Insights
                        </div>
                        <div class="card-body">
                            <h5>1. Hierarchical Conservation Architecture</h5>
                            <p>The four-layered architecture (Core → Buffer → Satellite → Distant) represents an elegant evolutionary strategy that preserves essential functions while allowing adaptation.</p>
                            
                            <h5>2. Regulatory Adaptation Mechanism</h5>
                            <p>Adaptation occurs primarily through regulatory changes mediated by satellite genes rather than direct modification of essential enzymes. This is supported by the predominance of upstream variants (80.35%).</p>
                            
                            <h5>3. Satellite Gene Architecture</h5>
                            <p>HIGH and MODERATE impact variants occur at specific, consistent distances from pathway genes, suggesting functional or regulatory relationships that allow for adaptation without disrupting essential processes.</p>
                            
                            <h5>4. Pattern of Amplification</h5>
                            <p>The perfect 4:3:1 ratio of variants in gene-modified:adapted:control samples suggests that adaptation and genetic modification amplify pre-existing genomic variation rather than generating novel mutations.</p>
                        </div>
                    </div>
                </div>
                
                <div class="col-lg-6">
                    <div class="card mb-4">
                        <div class="card-header bg-success text-white">
                            <i class="bi bi-journal-check"></i> Implications and Applications
                        </div>
                        <div class="card-body">
                            <h5>1. Evolutionary Conservation Model</h5>
                            <p>The hierarchical conservation pattern provides a model for understanding how essential pathways can evolve despite strong functional constraints. This could inform studies of conservation and adaptation in other organisms.</p>
                            
                            <h5>2. Regulatory Network Insights</h5>
                            <p>The identification of satellite genes with specific relationships to ergosterol pathway genes suggests new regulatory connections that could be targeted in studies of sterol metabolism.</p>
                            
                            <h5>3. Adaptation Mechanisms</h5>
                            <p>The finding that adaptation occurs through regulatory changes rather than enzyme modifications provides insight into how organisms can respond to environmental stresses without compromising essential functions.</p>
                            
                            <h5>4. Methodology for Conservation Analysis</h5>
                            <p>The analytical approach used here, combining genomic, network, and functional analyses, provides a template for studying conservation patterns in other essential pathways.</p>
                        </div>
                    </div>
                </div>
            </div>
            
            <div class="card mb-4">
                <div class="card-header">
                    <i class="bi bi-diagram-3"></i> Integrated Model of Yeast Adaptation
                </div>
                <div class="card-body">
                    <p>Our analysis suggests an integrated model of yeast adaptation where:</p>
                    <ol>
                        <li>Core pathway genes remain under strong purifying selection, maintaining essential cellular functions</li>
                        <li>A buffer zone extends ~7kb from each pathway gene, preserving regulatory regions</li>
                        <li>Satellite genes at specific distances (7-50kb) harbor HIGH/MODERATE impact variants that likely affect regulation of the ergosterol pathway</li>
                        <li>These satellite genes mediate adaptation to environmental stresses through regulatory changes rather than direct enzyme modifications</li>
                        <li>Adaptation and genetic modification amplify pre-existing genomic variation in a systematic pattern</li>
                    </ol>
                    <p>This model provides a framework for understanding how essential pathways can be preserved while allowing for adaptation to changing environmental conditions. It highlights the importance of studying not just the genes of interest, but also their broader genomic neighborhood and regulatory context.</p>
                </div>
            </div>
    """
    
    # Footer
    html += """
            <footer class="mt-5 pt-4 border-top text-center text-muted">
                <p>Yeast MSA Project - Functional Impact Analysis Report</p>
                <p>Generated on: """ + datetime.now().strftime("%B %d, %Y") + """</p>
            </footer>
        </div>
        
        <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/jquery@3.7.1/dist/jquery.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.7/js/dataTables.bootstrap5.min.js"></script>
        
        <script>
            // Initialize DataTables
            $(document).ready(function() {
                $('.table-datatable').DataTable({
                    responsive: true,
                    lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]]
                });
                
                // Enable tooltips
                const tooltipTriggerList = document.querySelectorAll('[data-bs-toggle="tooltip"]');
                const tooltipList = [...tooltipTriggerList].map(tooltipTriggerEl => new bootstrap.Tooltip(tooltipTriggerEl));
            });
            
            // Function to show image in modal
            function showImageInModal(img) {
                const modalImg = document.getElementById("modalImage");
                const modalTitle = document.getElementById("imageModalLabel");
                
                modalImg.src = img.src;
                modalTitle.innerText = img.alt;
            }
            
            // Dark mode toggle
            document.getElementById('toggleDarkMode').addEventListener('click', function(e) {
                e.preventDefault();
                document.documentElement.classList.toggle('dark-mode-enabled');
                
                // Update icon
                const icon = this.querySelector('i');
                if (icon.classList.contains('bi-moon')) {
                    icon.classList.remove('bi-moon');
                    icon.classList.add('bi-sun');
                    this.innerHTML = '<i class="bi bi-sun"></i> Toggle Light Mode';
                } else {
                    icon.classList.remove('bi-sun');
                    icon.classList.add('bi-moon');
                    this.innerHTML = '<i class="bi bi-moon"></i> Toggle Dark Mode';
                }
            });
            
            // Handle scrolling and active section in navigation
            window.addEventListener('scroll', function() {
                const sections = document.querySelectorAll('h2[id]');
                const navLinks = document.querySelectorAll('.navbar-nav .nav-link');
                
                let currentSection = '';
                
                sections.forEach(section => {
                    const sectionTop = section.offsetTop - 100;
                    if (window.scrollY >= sectionTop) {
                        currentSection = section.getAttribute('id');
                    }
                });
                
                navLinks.forEach(link => {
                    link.classList.remove('active');
                    if (link.getAttribute('href') === '#' + currentSection) {
                        link.classList.add('active');
                    }
                });
            });
        </script>
    </body>
    </html>
    """
    
    return html

def main():
    """Main function to collect data and generate the HTML report."""
    print("Collecting data from the Yeast MSA project...")
    data = collect_data()
    
    print("Generating HTML report...")
    html = generate_html_report(data)
    
    # Write HTML to file
    output_path = "/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/functional_impact.html"
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"HTML report generated: {output_path}")

if __name__ == "__main__":
    main()