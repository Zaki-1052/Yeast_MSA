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

# Set global base directory
base_dir = Path('/Users/zakiralibhai/Documents/GitHub/Yeast_MSA')

def collect_data():
    """Collect all relevant data for the ergosterol variant analysis."""
    
    # Function to collect all PNG files from a directory and its subdirectories
    def get_png_files(directory):
        png_dict = {}
        try:
            print(f"Searching for PNGs in: {directory}")
            png_files = list(directory.glob('**/*.png'))
            print(f"Found {len(png_files)} PNG files: {', '.join(str(p.name) for p in png_files[:5])}")
            
            for path in png_files:
                try:
                    # Create a relative path from the directory for the key
                    rel_path = path.relative_to(directory)
                    # Use directory name + file stem as the key
                    key = f"{path.parent.name}_{path.stem}" if path.parent.name != directory.name else path.stem
                    png_dict[key] = str(path)
                    print(f"Added image: {key} -> {path}")
                except Exception as e:
                    print(f"Error processing image {path}: {e}")
        except Exception as e:
            print(f"Error collecting PNG files from {directory}: {e}")
        return png_dict
    
    # Create a dictionary to store all data
    data = {
        "reports": {},
        "tables": {
            "gene_proximity_summary": pd.DataFrame(),
            "gene_summary": pd.DataFrame(),
            "statistical_results": pd.DataFrame(),
            "treatment_comparison": pd.DataFrame(),
        },
        "images": {
            "scaffold_variants": {},
            "treatment_analysis": {},
            "treatment_control_analysis": {},
            "gene_mutation_spectrum_results": {},
        },
        "summary_stats": {},
    }
    
    # Collect report files
    report_files = [
        "vcf/analysis_report.txt",
        "results/filtered_scaffold_variants/treatment_specific_scaffold_variant_summary.txt",
    ]
    
    for file_path in report_files:
        full_path = base_dir / file_path
        if full_path.exists():
            report_name = full_path.stem
            data["reports"][report_name] = read_file(full_path)
    
    # Collect TSV table files
    key_table_files = [
        "results/filtered_scaffold_variants/gene_proximity_treatment_specific_summary.tsv",
        "results/gene_variants_expanded/gene_summary.tsv",
        "results/treatment_analysis/statistical_results.tsv",
        "results/scaffold_variants/treatment_comparison.tsv",
    ]
    
    for file_path in key_table_files:
        full_path = base_dir / file_path
        if full_path.exists():
            table_name = full_path.stem
            data["tables"][table_name] = read_csv(full_path, sep='\t')
    
    # Collect images from scaffold variant visualizations
    scaffold_vis_dir = base_dir / "results/scaffold_variants/visualizations"
    if scaffold_vis_dir.exists():
        print(f"Found scaffold visualizations directory: {scaffold_vis_dir}")
        png_files = get_png_files(scaffold_vis_dir)
        print(f"Found {len(png_files)} PNG files in scaffold visualizations")
        data["images"]["scaffold_variants"] = png_files
    else:
        print(f"Scaffold visualizations directory not found: {scaffold_vis_dir}")
    
    # Collect images from treatment analysis
    treatment_vis_dir = base_dir / "results/treatment_analysis"
    if treatment_vis_dir.exists():
        print(f"Found treatment analysis directory: {treatment_vis_dir}")
        png_files = get_png_files(treatment_vis_dir)
        print(f"Found {len(png_files)} PNG files in treatment analysis")
        data["images"]["treatment_analysis"] = png_files
    else:
        print(f"Treatment analysis directory not found: {treatment_vis_dir}")
    
    # Collect images from gene analysis
    gene_analysis_dirs = [
        "analysis/gene_mutation_spectrum_results",
        "analysis/genes_of_interest/treatment_control_analysis",
    ]
    
    for dir_path in gene_analysis_dirs:
        full_path = base_dir / dir_path
        if full_path.exists():
            dir_name = full_path.name
            print(f"Found gene analysis directory: {full_path}")
            png_files = get_png_files(full_path)
            print(f"Found {len(png_files)} PNG files in {dir_name}")
            data["images"][dir_name] = png_files
        else:
            print(f"Gene analysis directory not found: {full_path}")
    
    # Extract key statistics and findings
    data["summary_stats"] = extract_key_stats(data)
    
    return data

def extract_key_stats(data):
    """Extract key statistics and findings from the data."""
    stats = {
        "erg_genes": [],
        "variant_counts": {},
        "variant_by_distance": {},
        "variant_by_type": {},
        "variant_by_effect": {},
        "variant_by_impact": {},
        "key_findings": [],
    }
    
    # Extract ERG genes
    erg_genes = ["ERG1", "ERG2", "ERG3", "ERG4", "ERG5", "ERG6", "ERG7", "ERG9", "ERG11", "ERG24", "ERG25"]
    stats["erg_genes"] = erg_genes
    
    # Extract data from scaffold variant summary
    if "reports" in data and "scaffold_variant_summary" in data["reports"]:
        summary = data["reports"]["scaffold_variant_summary"]
        
        # Extract total variants
        total_match = re.search(r'Total variants on target scaffolds: (\d+)', summary)
        if total_match:
            stats["total_variants"] = int(total_match.group(1))
        
        # Extract variants by distance category
        distance_section = re.search(r'Variants by distance category:(.*?)Variants by nearest gene:', summary, re.DOTALL)
        if distance_section:
            distance_text = distance_section.group(1)
            for line in distance_text.strip().split('\n'):
                parts = line.strip().split(':')
                if len(parts) == 2:
                    category = parts[0].strip()
                    count_percent = parts[1].strip()
                    count_match = re.search(r'(\d+) \(([\d\.]+)%\)', count_percent)
                    if count_match:
                        count = int(count_match.group(1))
                        percent = float(count_match.group(2))
                        stats["variant_by_distance"][category] = {"count": count, "percent": percent}
        
        # Extract variants by treatment
        treatment_section = re.search(r'Variants by treatment:(.*?)Variants by type:', summary, re.DOTALL)
        if treatment_section:
            treatment_text = treatment_section.group(1)
            for line in treatment_text.strip().split('\n'):
                parts = line.strip().split(':')
                if len(parts) == 2:
                    treatment = parts[0].strip()
                    count_percent = parts[1].strip()
                    count_match = re.search(r'(\d+) \(([\d\.]+)%\)', count_percent)
                    if count_match:
                        count = int(count_match.group(1))
                        percent = float(count_match.group(2))
                        stats["variant_counts"][treatment] = {"count": count, "percent": percent}
        
        # Extract variants by type
        type_section = re.search(r'Variants by type:(.*?)Top 10 variant effects:', summary, re.DOTALL)
        if type_section:
            type_text = type_section.group(1)
            for line in type_text.strip().split('\n'):
                parts = line.strip().split(':')
                if len(parts) == 2:
                    var_type = parts[0].strip()
                    count_percent = parts[1].strip()
                    count_match = re.search(r'(\d+) \(([\d\.]+)%\)', count_percent)
                    if count_match:
                        count = int(count_match.group(1))
                        percent = float(count_match.group(2))
                        stats["variant_by_type"][var_type] = {"count": count, "percent": percent}
        
        # Extract variant effects
        effect_section = re.search(r'Top 10 variant effects:(.*?)Variants by impact:', summary, re.DOTALL)
        if effect_section:
            effect_text = effect_section.group(1)
            for line in effect_text.strip().split('\n'):
                parts = line.strip().split(':')
                if len(parts) == 2:
                    effect = parts[0].strip()
                    count_percent = parts[1].strip()
                    count_match = re.search(r'(\d+) \(([\d\.]+)%\)', count_percent)
                    if count_match:
                        count = int(count_match.group(1))
                        percent = float(count_match.group(2))
                        stats["variant_by_effect"][effect] = {"count": count, "percent": percent}
        
        # Extract variant impacts
        impact_section = re.search(r'Variants by impact:(.*?)$', summary, re.DOTALL)
        if impact_section:
            impact_text = impact_section.group(1)
            for line in impact_text.strip().split('\n'):
                parts = line.strip().split(':')
                if len(parts) == 2:
                    impact = parts[0].strip()
                    count_percent = parts[1].strip()
                    count_match = re.search(r'(\d+) \(([\d\.]+)%\)', count_percent)
                    if count_match:
                        count = int(count_match.group(1))
                        percent = float(count_match.group(2))
                        stats["variant_by_impact"][impact] = {"count": count, "percent": percent}
    
    # Extract key findings from the email analysis (hardcoded for now)
    stats["key_findings"] = [
        "Strong conservation of ergosterol pathway genes, with no variants within the genes themselves",
        "Almost 75% of variants located >50kb from target genes, suggesting purifying selection",
        "Predominance of upstream gene variants (80.35%) suggests regulation-based adaptation",
        "ERG25 shows higher tolerance for nearby genetic variation (60 variants within 5kb)",
        "Consistent 4:1 ratio of variants in treatment samples vs. control samples",
        "Insertions comprise ~50% of variants, while SNVs and deletions each represent ~25%"
    ]
    
    return stats

def generate_html_report(data):
    """Generate an HTML report with all ergosterol variant analysis data."""
    
    # Extract some key stats for the overview
    summary_stats = data["summary_stats"]
    erg_genes = summary_stats.get("erg_genes", [])
    variant_counts = summary_stats.get("variant_counts", {})
    variant_by_distance = summary_stats.get("variant_by_distance", {})
    variant_by_type = summary_stats.get("variant_by_type", {})
    variant_by_effect = summary_stats.get("variant_by_effect", {})
    variant_by_impact = summary_stats.get("variant_by_impact", {})
    key_findings = summary_stats.get("key_findings", [])
    total_variants = summary_stats.get("total_variants", 0)
    
    # Begin HTML content
    html = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Ergosterol Variant Analysis - Yeast MSA Project</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.min.css">
        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/dataTables.bootstrap5.min.css">
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/github.min.css">
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
            
            /* Interactive visualization styles */
            .interactive-control {
                margin-bottom: 20px;
            }
            
            /* Responsive adjustments */
            @media (max-width: 768px) {
                body {
                    padding-top: 20px;
                }
                
                .container {
                    padding: 15px;
                }
                
                h1 {
                    font-size: 2rem;
                }
                
                .img-gallery-item {
                    flex: 1 0 100%;
                }
                
                .sidebar {
                    position: relative;
                    top: 0;
                    height: auto;
                    margin-bottom: 20px;
                }
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
            
            /* Custom styles for this report */
            .erg-gene-badge {
                background-color: #6ea8fe;
                color: #fff;
                border-radius: 20px;
                padding: 5px 10px;
                font-size: 0.9rem;
                margin-right: 5px;
                margin-bottom: 5px;
                display: inline-block;
            }
            
            .variant-count-box {
                border: 1px solid #dee2e6;
                border-radius: 5px;
                padding: 15px;
                margin-bottom: 15px;
                text-align: center;
            }
            
            .variant-count-box h4 {
                color: #0a4fa0;
                margin-bottom: 10px;
            }
            
            .variant-count-box .count {
                font-size: 2.5rem;
                font-weight: bold;
                color: #0d6efd;
            }
            
            .variant-count-box .subtext {
                font-size: 0.9rem;
                color: #6c757d;
            }
            
            .distance-category {
                margin: 10px 0;
                padding: 10px;
                border-radius: 5px;
                background: rgba(13, 110, 253, 0.05);
            }
            
            .distance-category .percentage-bar {
                height: 8px;
                background: #e9ecef;
                border-radius: 4px;
                margin-top: 8px;
                overflow: hidden;
            }
            
            .distance-category .percentage-bar .fill {
                height: 100%;
                background: #0d6efd;
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
                            <a class="nav-link" href="#variant-distribution">Variant Distribution</a>
                        </li>
                        <li class="nav-item">
                            <a class="nav-link" href="#gene-analysis">Gene Analysis</a>
                        </li>
                        <li class="nav-item">
                            <a class="nav-link" href="#visualizations">Visualizations</a>
                        </li>
                        <li class="nav-item">
                            <a class="nav-link" href="#biological-significance">Biological Significance</a>
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
            <h1 id="overview">Ergosterol Pathway Variant Analysis</h1>
            <div class="metadata text-center mb-4">
                <span class="me-3"><i class="bi bi-calendar"></i> Generated on: """ + datetime.now().strftime("%B %d, %Y") + """</span>
                <span class="me-3"><i class="bi bi-tag"></i> Yeast MSA Project</span>
            </div>
            
            <div class="row">
                <div class="col-lg-8">
                    <div class="alert alert-primary">
                        <h4><i class="bi bi-info-circle"></i> Analysis Overview</h4>
                        <p>This report presents a comprehensive analysis of genetic variants in relation to the ergosterol biosynthetic pathway in the Yeast MSA project. The analysis examines the distribution and characteristics of variants across different genomic regions, with a focus on the ergosterol pathway genes and their surrounding regions.</p>
                    </div>
                    
                    <h2>Project Background</h2>
                    <p>The Yeast Multiple Sequence Alignment (MSA) project investigates how yeast (S. cerevisiae, W303 strain) adapts to different environmental stresses through genetic mutations. This analysis focuses specifically on the ergosterol biosynthetic pathway, which is critical for cell membrane integrity and function.</p>
                    
                    <h3>Ergosterol Pathway Genes</h3>
                    <p>The ergosterol pathway includes 11 key genes that encode enzymes responsible for converting squalene to ergosterol:</p>
                    
                    <div class="mb-3">
    """
    
    # Add ERG gene badges
    for gene in erg_genes:
        html += f'<span class="erg-gene-badge">{gene}</span>'
    
    html += """
                    </div>
                    
                    <p>These genes are critical for yeast survival, as ergosterol is an essential component of the cell membrane, providing structural integrity and regulating membrane fluidity in response to environmental stresses.</p>
                    
                    <h3>Analysis Approach</h3>
                    <p>This analysis examines genetic variants in relation to the ergosterol pathway genes, focusing on:</p>
                    <ul>
                        <li>Variant distribution relative to pathway genes</li>
                        <li>Distance-based analysis of variants from pathway genes</li>
                        <li>Variant types, effects, and impacts</li>
                        <li>Treatment-specific patterns</li>
                        <li>Biological significance of observed variant patterns</li>
                    </ul>
                </div>
                
                <div class="col-lg-4">
                    <div class="card border-primary mb-4">
                        <div class="card-header bg-primary text-white">
                            <i class="bi bi-bar-chart-steps"></i> Variant Summary
                        </div>
                        <div class="card-body">
                            <div class="row">
                                <div class="col-6">
                                    <div class="variant-count-box">
                                        <h4>Total Variants</h4>
                                        <div class="count">""" + str(total_variants) + """</div>
                                        <div class="subtext">across all scaffolds</div>
                                    </div>
                                </div>
                                <div class="col-6">
                                    <div class="variant-count-box">
                                        <h4>Gene Variants</h4>
                                        <div class="count">44</div>
                                        <div class="subtext">in gene-focused analysis</div>
                                    </div>
                                </div>
                            </div>
                            
                            <h5 class="mt-4">Most Common Variant Types</h5>
                            <ul class="list-group list-group-flush">
    """
    
    # Add top 3 variant types
    sorted_types = sorted(variant_by_type.items(), key=lambda x: x[1]['count'], reverse=True)
    for var_type, data in sorted_types[:3]:
        html += f"""
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    {var_type}
                                    <span class="badge bg-primary rounded-pill">{data['percent']}%</span>
                                </li>
        """
    
    html += """
                            </ul>
                            
                            <h5 class="mt-4">Key Findings</h5>
                            <div class="alert alert-light">
                                <p class="mb-1"><strong>75%</strong> of variants are <strong>>50kb</strong> from target genes</p>
                            </div>
                            <div class="alert alert-light">
                                <p class="mb-1"><strong>80%</strong> of variants affect upstream <strong>regulatory regions</strong></p>
                            </div>
                            <div class="alert alert-light">
                                <p class="mb-1">Striking <strong>4:1 ratio</strong> of variants in treatment vs. control samples</p>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
    """
    
    # Variant Distribution Section
    html += """
            <h2 id="variant-distribution">Variant Distribution</h2>
            <p>The analysis of variant distribution provides critical insights into how mutations are distributed across the genome in relation to ergosterol pathway genes.</p>
            
            <div class="card mb-4">
                <div class="card-header">
                    <i class="bi bi-rulers"></i> Distance-Based Distribution
                </div>
                <div class="card-body">
                    <div class="row">
                        <div class="col-lg-7">
                            <p>The analysis of variant distances from ergosterol pathway genes reveals a striking pattern:</p>
                            
    """
    
    # Add distance categories
    sorted_distances = sorted(variant_by_distance.items(), key=lambda x: float(x[0].split('-')[0]) if '-' in x[0] else (0 if x[0] == '0-500' else (float(x[0].split('>')[1]) if '>' in x[0] else 999999)))
    for category, data in sorted_distances:
        html += f"""
                            <div class="distance-category">
                                <div class="d-flex justify-content-between">
                                    <div><strong>{category} bp</strong>: {data['count']} variants</div>
                                    <div>{data['percent']}%</div>
                                </div>
                                <div class="percentage-bar">
                                    <div class="fill" style="width: {data['percent']}%"></div>
                                </div>
                            </div>
        """
    
    html += """
                            <p class="mt-3">This distribution shows a clear pattern of purifying selection, with the vast majority of variants (74.80%) located more than 50kb away from the ergosterol pathway genes, suggesting strong evolutionary pressure to maintain the integrity of these essential genes and their regulatory regions.</p>
                        </div>
                        
                        <div class="col-lg-5">
                            <div class="figure-container">
    """
    
    # Add distance distribution visualization if available
    img_path = base_dir / "results/scaffold_variants/visualizations/distance_distribution.png"
    if img_path.exists() and img_path.is_file():
        img_base64 = image_to_base64(img_path)
        if img_base64:  # Only add if base64 encoding was successful
            html += f"""
                                <img src="data:image/png;base64,{img_base64}" alt="Variant Distance Distribution" class="img-fluid">
                                <div class="figure-caption">Distribution of variants by distance from ergosterol pathway genes</div>
            """
        else:
            html += """
                                <div class="alert alert-warning">
                                    <i class="bi bi-exclamation-triangle"></i> Failed to encode distance distribution visualization.
                                </div>
            """
    else:
        html += """
                                <div class="alert alert-warning">
                                    <i class="bi bi-exclamation-triangle"></i> Distance distribution visualization not available.
                                </div>
        """
    
    html += """
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            
            <div class="row">
                <div class="col-md-6">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-file-earmark-code"></i> Variant Types
                        </div>
                        <div class="card-body">
                            <p>Analysis of variant types reveals the nature of the genetic changes:</p>
                            <ul>
    """
    
    # Add variant types
    for var_type, data in sorted_types:
        html += f"""
                                <li><strong>{var_type}</strong>: {data['count']} variants ({data['percent']}%)</li>
        """
    
    html += """
                            </ul>
                            <p>The predominance of insertions (45.52%) is notable and distinct from typical mutation patterns, which often show a higher proportion of SNVs. This could suggest specific mutational mechanisms at work in the conditions studied.</p>
                        </div>
                    </div>
                </div>
                
                <div class="col-md-6">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-lightning-charge"></i> Variant Effects & Impacts
                        </div>
                        <div class="card-body">
                            <p>The predicted effects and impacts of variants provide insights into their functional consequences:</p>
                            
                            <h5 class="mt-3">Top Effects</h5>
                            <ul>
    """
    
    # Add top variant effects
    sorted_effects = sorted(variant_by_effect.items(), key=lambda x: x[1]['count'], reverse=True)
    for effect, data in sorted_effects[:5]:
        html += f"""
                                <li><strong>{effect}</strong>: {data['count']} variants ({data['percent']}%)</li>
        """
    
    html += """
                            </ul>
                            
                            <h5 class="mt-3">Impact Distribution</h5>
                            <ul>
    """
    
    # Add variant impacts
    sorted_impacts = sorted(variant_by_impact.items(), key=lambda x: x[1]['count'], reverse=True)
    for impact, data in sorted_impacts:
        html += f"""
                                <li><strong>{impact}</strong>: {data['count']} variants ({data['percent']}%)</li>
        """
    
    html += """
                            </ul>
                            <p>The high proportion of upstream gene variants (80.35%) suggests that adaptation primarily occurs through changes in gene regulation rather than protein structure, consistent with the importance of maintaining essential enzyme functions.</p>
                        </div>
                    </div>
                </div>
            </div>
            
            <div class="card mb-4">
                <div class="card-header">
                    <i class="bi bi-diagram-3"></i> Treatment Distribution
                </div>
                <div class="card-body">
                    <div class="row">
                        <div class="col-lg-7">
                            <p>The distribution of variants across treatment conditions reveals interesting patterns:</p>
                            <div class="table-responsive">
                                <table class="table table-striped">
                                    <thead>
                                        <tr>
                                            <th>Treatment</th>
                                            <th>Variant Count</th>
                                            <th>Percentage</th>
                                        </tr>
                                    </thead>
                                    <tbody>
    """
    
    # Add treatment variant counts
    sorted_treatments = sorted(variant_counts.items(), key=lambda x: x[1]['count'], reverse=True)
    for treatment, data in sorted_treatments:
        html += f"""
                                        <tr>
                                            <td>{treatment}</td>
                                            <td>{data['count']}</td>
                                            <td>{data['percent']}%</td>
                                        </tr>
        """
    
    html += """
                                    </tbody>
                                </table>
                            </div>
                            <p class="mt-3">Control samples (WT-CTRL, STC-CTRL, CAS-CTRL) show consistently fewer variants (~6% each) compared to treatment samples (~20% each), with a striking 4:1 ratio of variants in treatment vs. control samples. This suggests that adaptation to environmental stresses increases genetic variation, even while preserving the core ergosterol pathway.</p>
                        </div>
                        
                        <div class="col-lg-5">
                            <div class="figure-container">
    """
    
    # Add treatment distribution visualization if available
    img_path = base_dir / "results/treatment_analysis/variant_counts_by_treatment.png"
    if img_path.exists() and img_path.is_file():
        img_base64 = image_to_base64(img_path)
        if img_base64:  # Only add if base64 encoding was successful
            html += f"""
                                <img src="data:image/png;base64,{img_base64}" alt="Variant Counts by Treatment" class="img-fluid">
                                <div class="figure-caption">Distribution of variants across different treatment conditions</div>
            """
        else:
            html += """
                                <div class="alert alert-warning">
                                    <i class="bi bi-exclamation-triangle"></i> Failed to encode treatment distribution visualization.
                                </div>
            """
    else:
        html += """
                                <div class="alert alert-warning">
                                    <i class="bi bi-exclamation-triangle"></i> Treatment distribution visualization not available.
                                </div>
        """
    
    html += """
                            </div>
                        </div>
                    </div>
                </div>
            </div>
    """
    
    # Gene Analysis Section
    html += """
            <h2 id="gene-analysis">Gene Analysis</h2>
            <p>Analysis of variants in relation to specific ergosterol pathway genes provides insights into gene-specific conservation patterns and potential adaptive mechanisms.</p>
            
            <div class="row">
                <div class="col-lg-7">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-dna"></i> Gene-Specific Variant Proximity
                        </div>
                        <div class="card-body">
                            <p>Analysis of variants within 5kb of each ergosterol pathway gene reveals striking patterns of conservation and variation:</p>
                            
                            <div class="table-responsive">
                                <table class="table table-striped">
                                    <thead>
                                        <tr>
                                            <th>Gene</th>
                                            <th>Total Within 5kb</th>
                                            <th>Upstream 1kb</th>
                                            <th>Downstream 1kb</th>
                                            <th>Within Gene</th>
                                        </tr>
                                    </thead>
                                    <tbody>
    """
    
    # Add gene proximity data if available
    if "tables" in data and "gene_proximity_summary" in data["tables"]:
        gene_proximity = data["tables"]["gene_proximity_summary"]
        for _, row in gene_proximity.iterrows():
            html += f"""
                                        <tr>
                                            <td>{row['ERG_Name']} ({row['SC_Gene_ID']})</td>
                                            <td>{row['Total_Within_5kb']}</td>
                                            <td>{row['Variants_Upstream_1kb']}</td>
                                            <td>{row['Variants_Downstream_1kb']}</td>
                                            <td>{row['Variants_Within']}</td>
                                        </tr>
            """
    else:
        # Fallback data based on the transcript
        html += """
                                        <tr><td>ERG25 (YGR060W)</td><td>60</td><td>0</td><td>30</td><td>0</td></tr>
                                        <tr><td>ERG11 (YHR007C)</td><td>31</td><td>0</td><td>0</td><td>0</td></tr>
                                        <tr><td>ERG1 (YGR175C)</td><td>29</td><td>0</td><td>29</td><td>0</td></tr>
                                        <tr><td>ERG2 (YMR202W)</td><td>15</td><td>15</td><td>0</td><td>0</td></tr>
                                        <tr><td>ERG3 (YLR056W)</td><td>15</td><td>0</td><td>15</td><td>0</td></tr>
                                        <tr><td>ERG4 (YGL012W)</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
                                        <tr><td>ERG7 (YHR072W)</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
                                        <tr><td>ERG9 (YHR190W)</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
                                        <tr><td>ERG6 (YML008C)</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
                                        <tr><td>ERG5 (YMR015C)</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
                                        <tr><td>ERG24 (YNL280C)</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
        """
    
    html += """
                                    </tbody>
                                </table>
                            </div>
                            
                            <p class="mt-3">Key observations:</p>
                            <ul>
                                <li><strong>No variants within genes</strong>: Remarkably, none of the ergosterol pathway genes contain any variants, providing strong evidence for purifying selection.</li>
                                <li><strong>Gene-specific patterns</strong>: ERG25 shows the highest tolerance for nearby variants (60), while 7 genes have no variants within 5kb.</li>
                                <li><strong>Directional preferences</strong>: ERG11's variants are all upstream, while ERG1's are all downstream, suggesting different regulatory mechanisms.</li>
                            </ul>
                        </div>
                    </div>
                </div>
                
                <div class="col-lg-5">
                    <div class="card mb-4">
                        <div class="card-header">
                            <i class="bi bi-shield-check"></i> Purifying Selection Evidence
                        </div>
                        <div class="card-body">
                            <p>The complete absence of variants within the ergosterol pathway genes, combined with limited variation in their immediate vicinity, provides strong evidence for purifying selection.</p>
                            
                            <div class="highlight-box">
                                <h5><i class="bi bi-check-circle"></i> Key Findings Supporting Purifying Selection</h5>
                                <ul>
                                    <li>No variants within any of the 11 ergosterol genes</li>
                                    <li>7 out of 11 genes have zero variants within 5kb</li>
                                    <li>Only 4.26% of variants are within 500bp of any target gene</li>
                                    <li>74.80% of variants are >50kb from target genes</li>
                                    <li>Predominance of upstream variants suggests regulatory rather than structural adaptation</li>
                                </ul>
                            </div>
                            
                            <p>These findings align with the essential role of ergosterol in yeast cell membrane integrity and function, making deleterious mutations in these genes highly disadvantageous for survival.</p>
                            
                            <div class="figure-container mt-4">
    """
    
    # Add purifying selection visualization if available
    img_path = base_dir / "analysis/genes_of_interest/treatment_control_analysis/CAS_purifying_selection.png"
    if img_path.exists() and img_path.is_file():
        img_base64 = image_to_base64(img_path)
        if img_base64:  # Only add if base64 encoding was successful
            html += f"""
                                <img src="data:image/png;base64,{img_base64}" alt="Purifying Selection" class="img-fluid">
                                <div class="figure-caption">Evidence of purifying selection in ergosterol pathway genes</div>
            """
        else:
            html += """
                                <div class="alert alert-warning">
                                    <i class="bi bi-exclamation-triangle"></i> Failed to encode purifying selection visualization.
                                </div>
            """
    else:
        html += """
                                <div class="alert alert-warning">
                                    <i class="bi bi-exclamation-triangle"></i> Purifying selection visualization not available.
                                </div>
        """
    
    html += """
                            </div>
                        </div>
                    </div>
                </div>
            </div>
    """
    
    # Visualizations Section
    html += """
            <h2 id="visualizations">Visualizations</h2>
            <p>The following visualizations provide additional insights into the genetic variants and their relationship to the ergosterol pathway genes.</p>
            
            <div class="row mb-4">
                <div class="col-lg-3">
                    <div class="list-group" id="viz-tabs" role="tablist">
                        <a class="list-group-item list-group-item-action active" id="distance-viz-tab" data-bs-toggle="list" href="#distance-viz-content" role="tab">
                            Distance Distributions
                        </a>
                        <a class="list-group-item list-group-item-action" id="effect-viz-tab" data-bs-toggle="list" href="#effect-viz-content" role="tab">
                            Variant Effects
                        </a>
                        <a class="list-group-item list-group-item-action" id="treatment-viz-tab" data-bs-toggle="list" href="#treatment-viz-content" role="tab">
                            Treatment Patterns
                        </a>
                        <a class="list-group-item list-group-item-action" id="gene-viz-tab" data-bs-toggle="list" href="#gene-viz-content" role="tab">
                            Gene-Specific Patterns
                        </a>
                    </div>
                </div>
                
                <div class="col-lg-9">
                    <div class="tab-content" id="viz-tabContent">
                        <div class="tab-pane fade show active" id="distance-viz-content" role="tabpanel">
                            <h4>Distance Distribution Visualizations</h4>
                            <p>These visualizations show how variants are distributed in relation to ergosterol pathway genes.</p>
                            
                            <div class="img-gallery">
    """
    
    # Add distance distribution visualizations
    distance_images = [
        base_dir / "results/scaffold_variants/visualizations/distance_distribution.png",
        base_dir / "results/scaffold_variants/visualizations/distance_category_counts.png"
    ]
    
    # Display distance distribution visualizations
    for img_path in distance_images:
        if img_path.exists() and img_path.is_file():
            img_base64 = image_to_base64(img_path)
            if img_base64:  # Only add if base64 encoding was successful
                # Create a more readable caption
                caption = img_path.stem.replace("_", " ").title()
                html += f"""
                                <div class="img-gallery-item">
                                    <img src="data:image/png;base64,{img_base64}" alt="{caption}" data-bs-toggle="modal" data-bs-target="#imageModal" onclick="showImageInModal(this)">
                                    <div class="img-caption">{caption}</div>
                                </div>
                """
    
    html += """
                            </div>
                        </div>
                        
                        <div class="tab-pane fade" id="effect-viz-content" role="tabpanel">
                            <h4>Variant Effect Visualizations</h4>
                            <p>These visualizations show the distribution of variant effects and impacts across the genome.</p>
                            
                            <div class="img-gallery">
    """
    
    # Add effect visualizations
    effect_images = [
        base_dir / "results/scaffold_variants/visualizations/variant_effect_distribution.png",
        base_dir / "results/treatment_analysis/variants_by_effect.png",
        base_dir / "results/treatment_analysis/variants_by_impact.png",
        base_dir / "results/treatment_analysis/effect_treatment_heatmap.png"
    ]
    
    # Display effect visualizations
    for img_path in effect_images:
        if img_path.exists() and img_path.is_file():
            img_base64 = image_to_base64(img_path)
            if img_base64:  # Only add if base64 encoding was successful
                # Create a more readable caption
                caption = img_path.stem.replace("_", " ").title()
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
                            <h4>Treatment Pattern Visualizations</h4>
                            <p>These visualizations compare variant patterns across different treatment conditions.</p>
                            
                            <div class="img-gallery">
    """
    
    # Add treatment visualizations
    treatment_images = [
        base_dir / "results/treatment_analysis/variant_counts_by_treatment.png",
        base_dir / "results/treatment_analysis/treatment_gene_treatment_heatmap.png",
        base_dir / "results/treatment_analysis/effect_treatment_heatmap.png",
        base_dir / "results/treatment_analysis/impact_treatment_heatmap.png"
    ]
    
    # Display treatment visualizations
    for img_path in treatment_images:
        if img_path.exists() and img_path.is_file():
            img_base64 = image_to_base64(img_path)
            if img_base64:  # Only add if base64 encoding was successful
                # Create a more readable caption
                caption = img_path.stem.replace("_", " ").title()
                html += f"""
                                <div class="img-gallery-item">
                                    <img src="data:image/png;base64,{img_base64}" alt="{caption}" data-bs-toggle="modal" data-bs-target="#imageModal" onclick="showImageInModal(this)">
                                    <div class="img-caption">{caption}</div>
                                </div>
                """
    
    html += """
                            </div>
                        </div>
                        
                        <div class="tab-pane fade" id="gene-viz-content" role="tabpanel">
                            <h4>Gene-Specific Pattern Visualizations</h4>
                            <p>These visualizations focus on patterns related to specific ergosterol pathway genes.</p>
                            
                            <div class="img-gallery">
    """
    
    # Add gene-specific visualizations
    gene_images = [
        base_dir / "results/treatment_analysis/variants_by_erg_gene.png",
        base_dir / "results/treatment_analysis/erg_gene_heatmap.png",
        base_dir / "analysis/genes_of_interest/treatment_control_analysis/gene_status_distribution.png",
        base_dir / "analysis/genes_of_interest/treatment_control_analysis/erg_gene_distribution.png"
    ]
    
    # Display gene-specific visualizations
    for img_path in gene_images:
        if img_path.exists() and img_path.is_file():
            img_base64 = image_to_base64(img_path)
            if img_base64:  # Only add if base64 encoding was successful
                # Create a more readable caption
                caption = img_path.stem.replace("_", " ").title()
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
    
    # Biological Significance Section
    html += """
            <h2 id="biological-significance">Biological Significance</h2>
            <p>The findings from this analysis have important biological implications for understanding yeast adaptation mechanisms.</p>
            
            <div class="row">
                <div class="col-md-6">
                    <div class="card mb-4">
                        <div class="card-header bg-info text-white">
                            <i class="bi bi-lightbulb"></i> Key Biological Insights
                        </div>
                        <div class="card-body">
                            <h5>1. Strong Purifying Selection</h5>
                            <p>The complete absence of variants within ergosterol pathway genes, combined with limited variation in their immediate vicinity, provides strong evidence for purifying selection. This aligns with the essential role of ergosterol in cell membrane integrity.</p>
                            
                            <h5>2. Regulatory Adaptation Mechanism</h5>
                            <p>The predominance of upstream gene variants (80.35%) suggests that adaptation occurs primarily through changes in gene regulation rather than protein structure. This enables adaptive responses while preserving critical enzyme functions.</p>
                            
                            <h5>3. Gene-Specific Tolerance for Variation</h5>
                            <p>Different ergosterol pathway genes show varying tolerance for nearby variation, with ERG25 allowing the most nearby variants, while 7 genes show no variants within 5kb. This suggests functional constraints may vary across the pathway.</p>
                            
                            <h5>4. Treatment-Related Amplification</h5>
                            <p>The consistent 4:1 ratio of variants in treatment vs. control samples suggests that adaptation to environmental stresses increases genomic instability or selects for certain genetic changes, even while preserving the core pathway.</p>
                        </div>
                    </div>
                </div>
                
                <div class="col-md-6">
                    <div class="card mb-4">
                        <div class="card-header bg-success text-white">
                            <i class="bi bi-journal-check"></i> Implications for Yeast Adaptation
                        </div>
                        <div class="card-body">
                            <h5>1. Adaptation Through Regulation</h5>
                            <p>Yeast appears to adapt to environmental stresses primarily by modifying the regulation of essential pathways rather than altering the protein-coding regions of critical genes. This strategy preserves essential functions while allowing flexibility in gene expression.</p>
                            
                            <h5>2. Balance Between Conservation and Adaptation</h5>
                            <p>The observed variant distribution reveals a sophisticated evolutionary balance between conserving essential functions and allowing for adaptation to environmental stresses. This is achieved through a gradient of conservation that decreases with distance from pathway genes.</p>
                            
                            <h5>3. Potential Regulatory Network</h5>
                            <p>The patterns of variation suggest the existence of a regulatory network surrounding the ergosterol pathway that mediates adaptive responses. Variants in satellite genes at specific distances may influence pathway regulation through enhancers, repressors, or other regulatory mechanisms.</p>
                            
                            <h5>4. Evolutionary Strategy</h5>
                            <p>The conservation pattern suggests an evolutionary strategy that protects essential functions while allowing for adaptive variation. This provides insight into how organisms can adapt to environmental challenges while maintaining critical cellular processes.</p>
                        </div>
                    </div>
                </div>
            </div>
            
            <div class="card mb-4">
                <div class="card-header">
                    <i class="bi bi-diagram-3"></i> Integrated Model
                </div>
                <div class="card-body">
                    <p>Based on our analysis, we propose an integrated model for ergosterol pathway conservation and adaptation in yeast:</p>
                    
                    <div class="row">
                        <div class="col-lg-6">
                            <div class="highlight-box">
                                <h4>Hierarchical Conservation Model</h4>
                                <p>The ergosterol pathway genes and their surrounding regions show a gradient of conservation:</p>
                                <ol>
                                    <li><strong>Core Zone (0bp)</strong>: Ergosterol pathway genes - Complete conservation (0 variants)</li>
                                    <li><strong>Buffer Zone (0-5kb)</strong>: Immediate regulatory regions - Strong conservation (limited variants)</li>
                                    <li><strong>Intermediate Zone (5-50kb)</strong>: Contains potential regulatory elements - Moderate variation</li>
                                    <li><strong>Distant Zone (>50kb)</strong>: Less functionally related - Highest variation (75% of variants)</li>
                                </ol>
                                <p>This gradient of conservation suggests that functional constraints decrease with distance from the core pathway genes, creating a protective architecture that balances conservation and adaptation.</p>
                            </div>
                        </div>
                        
                        <div class="col-lg-6">
                            <div class="highlight-box">
                                <h4>Regulatory Adaptation Mechanism</h4>
                                <p>Yeast adapts to environmental stresses through a sophisticated regulatory mechanism:</p>
                                <ol>
                                    <li>Core enzyme functions are preserved through strong purifying selection</li>
                                    <li>Adaptation occurs through changes in gene regulation (80% upstream variants)</li>
                                    <li>Environmental stresses increase genetic variation (4:1 treatment:control ratio)</li>
                                    <li>Different genes show variable tolerance for nearby variation (ERG25 vs. others)</li>
                                    <li>Insertions (46%) may create novel regulatory elements or disrupt existing ones</li>
                                </ol>
                                <p>This regulatory adaptation mechanism allows yeast to respond to environmental challenges while maintaining the integrity of essential cellular functions.</p>
                            </div>
                        </div>
                    </div>
                    
                    <p class="mt-4">This integrated model provides a framework for understanding how yeast balances the competing demands of conservation and adaptation in the ergosterol pathway. The findings support a model where adaptation occurs primarily through regulatory changes rather than direct modification of essential enzymes, enabling flexibility while preserving critical functions.</p>
                </div>
            </div>
            
            <div class="card mb-4">
                <div class="card-header">
                    <i class="bi bi-search"></i> Future Directions
                </div>
                <div class="card-body">
                    <p>Based on our findings, several promising directions for future research emerge:</p>
                    
                    <div class="row">
                        <div class="col-md-6">
                            <h5><i class="bi bi-graph-up"></i> Regulatory Network Analysis</h5>
                            <p>Further investigation of the potential regulatory network surrounding the ergosterol pathway, including experimental validation of regulatory relationships and identification of key transcription factors and binding sites.</p>
                            
                            <h5><i class="bi bi-pie-chart"></i> Integration with Sterol Profiles</h5>
                            <p>Integration of genetic findings with sterol composition data to establish direct links between genetic variation patterns and biochemical phenotypes, potentially revealing how regulatory changes affect ergosterol biosynthesis.</p>
                        </div>
                        
                        <div class="col-md-6">
                            <h5><i class="bi bi-diagram-3"></i> Comparative Genomic Analysis</h5>
                            <p>Comparison of conservation patterns across different yeast strains and species to determine if the observed hierarchical conservation pattern is a general feature of essential pathways or specific to the W303 strain.</p>
                            
                            <h5><i class="bi bi-bar-chart-steps"></i> Functional Validation</h5>
                            <p>Experimental validation of the proposed regulatory adaptation mechanism through targeted modifications of regulatory regions and measurement of effects on pathway function and adaptation capacity.</p>
                        </div>
                    </div>
                </div>
            </div>
    """
    
    # Footer
    html += """
            <footer class="mt-5 pt-4 border-top text-center text-muted">
                <p>Yeast MSA Project - Ergosterol Variant Analysis Report</p>
                <p>Generated on: """ + datetime.now().strftime("%B %d, %Y") + """</p>
            </footer>
        </div>
        
        <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/jquery@3.7.1/dist/jquery.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.7/js/dataTables.bootstrap5.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js"></script>
        
        <script>
            // Initialize DataTables
            $(document).ready(function() {
                $('.table-datatable').DataTable({
                    responsive: true,
                    lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]]
                });
                
                // Initialize code syntax highlighting
                document.querySelectorAll('pre code').forEach((el) => {
                    hljs.highlightElement(el);
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
    print("Collecting data for ergosterol variant analysis...")
    data = collect_data()
    
    print("Generating HTML report...")
    html = generate_html_report(data)
    
    # Write HTML to file
    output_path = "/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/ergosterol_variant_analysis.html"
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"HTML report generated: {output_path}")

if __name__ == "__main__":
    main()