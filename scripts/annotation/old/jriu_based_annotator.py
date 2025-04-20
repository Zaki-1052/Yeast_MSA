#!/usr/bin/env python3
"""
jriu_based_annotator.py - Annotate variants based on JRIU IDs rather than exact positions
"""

import os
import sys
import csv
import argparse
import gzip
from collections import defaultdict, Counter
from Bio import SeqIO

def create_gene_info(genes_of_interest_file):
    """Create a mapping of JRIU IDs to gene information"""
    # Parse genes of interest file
    jriu_to_genes = defaultdict(list)
    sc_gene_map = {}
    
    with open(genes_of_interest_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene_id = row.get('w303_gene_id')
            sc_gene_id = row.get('sc_gene_id')
            jriu_id = row.get('jriu_id')
            
            if gene_id and jriu_id:
                jriu_to_genes[jriu_id].append({
                    'gene_id': gene_id,
                    'sc_gene_id': sc_gene_id
                })
                
                if sc_gene_id:
                    sc_gene_map[sc_gene_id] = gene_id
    
    return jriu_to_genes, sc_gene_map

def process_vcf_file(vcf_file, jriu_to_genes, sc_gene_map, output_file, treatment_map):
    """Process a VCF file to identify variants on JRIUs containing genes of interest"""
    # Extract sample name
    sample_name = os.path.basename(vcf_file).split('.')[0]
    treatment = treatment_map.get(sample_name, 'unknown')
    
    # Determine if file is gzipped
    is_gzipped = vcf_file.endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    # Initialize counters
    total_variants = 0
    jriu_variants = 0
    
    # Statistics by gene and treatment
    gene_variant_counts = defaultdict(int)
    gene_treatment_counts = defaultdict(int)
    
    # Store annotated variants
    annotated_variants = []
    
    # Process VCF file
    with opener(vcf_file, mode) as f:
        # Process header
        header_lines = []
        for line in f:
            if line.startswith('##'):
                header_lines.append(line)
                continue
            elif line.startswith('#CHROM'):
                header_lines.append(line)
                break
        
        # Process variants
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split('\t')
            total_variants += 1
            
            # Extract variant info
            chrom = fields[0]  # JRIU ID
            pos = fields[1]
            variant_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter_val = fields[6]
            info = fields[7]
            
            # Check if this JRIU has genes of interest
            if chrom in jriu_to_genes:
                jriu_variants += 1
                genes = jriu_to_genes[chrom]
                
                # Update counts for each gene
                for gene in genes:
                    gene_id = gene['gene_id']
                    sc_gene_id = gene['sc_gene_id']
                    
                    gene_variant_counts[gene_id] += 1
                    gene_treatment_counts[f"{gene_id}:{treatment}"] += 1
                
                # Add annotation to variant
                gene_info = ";".join([f"{g['gene_id']}({g['sc_gene_id']})" for g in genes])
                
                # Create annotated variant
                annotated_variants.append({
                    'chrom': chrom,
                    'pos': pos,
                    'id': variant_id,
                    'ref': ref,
                    'alt': alt,
                    'qual': qual,
                    'filter': filter_val,
                    'info': f"{info};GOI=true;GenesOnJRIU={gene_info}",
                    'genes': [g['gene_id'] for g in genes]
                })
    
    # Write annotated VCF
    with open(output_file, 'w') as f:
        # Write header
        for line in header_lines[:-1]:
            f.write(line + '\n')
        
        # Add custom INFO headers
        f.write('##INFO=<ID=GOI,Number=0,Type=Flag,Description="Variant is on a JRIU containing genes of interest">\n')
        f.write('##INFO=<ID=GenesOnJRIU,Number=.,Type=String,Description="Genes of interest on this JRIU">\n')
        
        # Write column names
        f.write(header_lines[-1] + '\n')
        
        # Write variants
        for variant in annotated_variants:
            f.write(f"{variant['chrom']}\t{variant['pos']}\t{variant['id']}\t"
                   f"{variant['ref']}\t{variant['alt']}\t{variant['qual']}\t"
                   f"{variant['filter']}\t{variant['info']}\tGT\t1\n")
    
    return {
        'sample': sample_name,
        'treatment': treatment,
        'total_variants': total_variants,
        'jriu_variants': jriu_variants,
        'gene_variant_counts': dict(gene_variant_counts),
        'gene_treatment_counts': dict(gene_treatment_counts)
    }

def process_all_vcfs(vcf_files, jriu_to_genes, sc_gene_map, output_dir, treatment_map):
    """Process all VCF files and generate summary reports"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each VCF file
    results = []
    all_gene_counts = defaultdict(int)
    all_gene_treatment_counts = defaultdict(int)
    
    for vcf_file in vcf_files:
        sample_name = os.path.basename(vcf_file).split('.')[0]
        output_file = os.path.join(output_dir, f"{sample_name}.annotated.vcf")
        
        print(f"Processing {vcf_file}")
        result = process_vcf_file(vcf_file, jriu_to_genes, sc_gene_map, output_file, treatment_map)
        results.append(result)
        
        # Update global counts
        for gene_id, count in result['gene_variant_counts'].items():
            all_gene_counts[gene_id] += count
        
        for gene_treat, count in result['gene_treatment_counts'].items():
            all_gene_treatment_counts[gene_treat] += count
        
        print(f"  Found {result['jriu_variants']}/{result['total_variants']} variants on JRIUs with genes of interest")
    
    # Create summary report
    summary_file = os.path.join(output_dir, "annotation_summary.tsv")
    with open(summary_file, 'w') as f:
        f.write("Sample\tTreatment\tTotalVariants\tGeneJRIUVariants\n")
        for result in results:
            f.write(f"{result['sample']}\t{result['treatment']}\t{result['total_variants']}\t{result['jriu_variants']}\n")
    
    # Create gene-specific report
    gene_report_file = os.path.join(output_dir, "gene_variant_summary.tsv")
    with open(gene_report_file, 'w') as f:
        f.write("GeneID\tSCGeneID\tTotalVariants\n")
        for gene_id, count in sorted(all_gene_counts.items(), key=lambda x: -x[1]):
            sc_gene_id = next((sc for sc, g in sc_gene_map.items() if g == gene_id), "unknown")
            f.write(f"{gene_id}\t{sc_gene_id}\t{count}\n")
    
    # Create treatment-specific report
    treatment_report_file = os.path.join(output_dir, "treatment_variant_summary.tsv")
    with open(treatment_report_file, 'w') as f:
        f.write("GeneID\tSCGeneID\tTreatment\tVariants\n")
        
        # Group by gene ID
        by_gene = defaultdict(dict)
        for gene_treat, count in all_gene_treatment_counts.items():
            gene_id, treatment = gene_treat.split(':', 1)
            by_gene[gene_id][treatment] = count
        
        # Write sorted by gene ID
        for gene_id, treatments in sorted(by_gene.items()):
            sc_gene_id = next((sc for sc, g in sc_gene_map.items() if g == gene_id), "unknown")
            for treatment, count in sorted(treatments.items()):
                f.write(f"{gene_id}\t{sc_gene_id}\t{treatment}\t{count}\n")
    
    # Generate HTML report
    html_report_file = os.path.join(output_dir, "annotation_report.html")
    generate_html_report(html_report_file, results, all_gene_counts, sc_gene_map, all_gene_treatment_counts)
    
    return {
        'summary_file': summary_file,
        'gene_report_file': gene_report_file,
        'treatment_report_file': treatment_report_file,
        'html_report_file': html_report_file
    }

def generate_html_report(html_file, results, gene_counts, sc_gene_map, gene_treatment_counts):
    """Generate an HTML report with charts and tables"""
    # Calculate treatment totals
    treatment_totals = defaultdict(int)
    for result in results:
        treatment_totals[result['treatment']] += result['jriu_variants']
    
    # Group gene treatment counts
    gene_treatment_data = defaultdict(dict)
    for gene_treat, count in gene_treatment_counts.items():
        gene_id, treatment = gene_treat.split(':', 1)
        sc_gene_id = next((sc for sc, g in sc_gene_map.items() if g == gene_id), gene_id)
        gene_name = sc_gene_id.split('Y')[0] if sc_gene_id.startswith('Y') else sc_gene_id
        gene_treatment_data[gene_name][treatment] = count
    
    # Generate HTML
    with open(html_file, 'w') as f:
        f.write("""<!DOCTYPE html>
<html>
<head>
    <title>Gene JRIU Variant Annotation Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2 { color: #333366; }
        .chart-container { width: 800px; height: 400px; margin: 20px 0; }
        .grid-container { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .note { background-color: #ffffdd; padding: 10px; border-left: 5px solid #ffcc00; margin: 20px 0; }
    </style>
</head>
<body>
    <h1>Gene JRIU Variant Annotation Report</h1>
    <p>Generated on: """ + f"{os.popen('date').read().strip()}" + """</p>
    
    <div class="note">
        <strong>Note:</strong> This report uses a JRIU-based approach rather than position-based matching, 
        due to coordinate system differences between the VCF variants and gene annotations.
        Variants are associated with genes if they occur on the same JRIU ID.
    </div>
    
    <div class="grid-container">
        <div>
            <h2>Variants by Treatment</h2>
            <div class="chart-container">
                <canvas id="treatmentChart"></canvas>
            </div>
        </div>
        
        <div>
            <h2>Variants by Gene</h2>
            <div class="chart-container">
                <canvas id="geneChart"></canvas>
            </div>
        </div>
    </div>
    
    <h2>Gene Variants by Treatment</h2>
    <div class="chart-container" style="width: 100%;">
        <canvas id="geneTreatmentChart"></canvas>
    </div>
    
    <h2>Sample Details</h2>
    <table>
        <tr>
            <th>Sample</th>
            <th>Treatment</th>
            <th>Total Variants</th>
            <th>Variants on Gene JRIUs</th>
            <th>Percentage</th>
        </tr>
""")
        
        # Add sample rows
        for result in sorted(results, key=lambda x: x['sample']):
            pct = result['jriu_variants'] / result['total_variants'] * 100 if result['total_variants'] > 0 else 0
            f.write(f"""        <tr>
            <td>{result['sample']}</td>
            <td>{result['treatment']}</td>
            <td>{result['total_variants']}</td>
            <td>{result['jriu_variants']}</td>
            <td>{pct:.2f}%</td>
        </tr>
""")
        
        # Continue with HTML structure
        f.write("""    </table>
    
    <script>
        // Treatment chart
        const treatmentCtx = document.getElementById('treatmentChart').getContext('2d');
        const treatmentChart = new Chart(treatmentCtx, {
            type: 'bar',
            data: {
                labels: """ + str(list(treatment_totals.keys())) + """,
                datasets: [{
                    label: 'Variant Count',
                    data: """ + str(list(treatment_totals.values())) + """,
                    backgroundColor: [
                        'rgba(54, 162, 235, 0.5)',
                        'rgba(255, 99, 132, 0.5)',
                        'rgba(75, 192, 192, 0.5)',
                        'rgba(255, 206, 86, 0.5)',
                        'rgba(153, 102, 255, 0.5)'
                    ],
                    borderColor: [
                        'rgba(54, 162, 235, 1)',
                        'rgba(255, 99, 132, 1)',
                        'rgba(75, 192, 192, 1)',
                        'rgba(255, 206, 86, 1)',
                        'rgba(153, 102, 255, 1)'
                    ],
                    borderWidth: 1
                }]
            },
            options: {
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Variant Count'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Treatment'
                        }
                    }
                }
            }
        });
        
        // Gene chart
        const geneCtx = document.getElementById('geneChart').getContext('2d');
        const geneChart = new Chart(geneCtx, {
            type: 'bar',
            data: {
                labels: """ + str([next((sc for sc, g in sc_gene_map.items() if g == gene_id), gene_id).split('Y')[0] if next((sc for sc, g in sc_gene_map.items() if g == gene_id), gene_id).startswith('Y') else next((sc for sc, g in sc_gene_map.items() if g == gene_id), gene_id) for gene_id in sorted(gene_counts.keys(), key=lambda x: -gene_counts[x])]) + """,
                datasets: [{
                    label: 'Variant Count',
                    data: """ + str([gene_counts[gene_id] for gene_id in sorted(gene_counts.keys(), key=lambda x: -gene_counts[x])]) + """,
                    backgroundColor: 'rgba(75, 192, 192, 0.5)',
                    borderColor: 'rgba(75, 192, 192, 1)',
                    borderWidth: 1
                }]
            },
            options: {
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Variant Count'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Gene'
                        }
                    }
                }
            }
        });
        
        // Gene Treatment chart
        const geneTreatmentCtx = document.getElementById('geneTreatmentChart').getContext('2d');
        const geneTreatmentChart = new Chart(geneTreatmentCtx, {
            type: 'bar',
            data: {
                labels: """ + str(list(gene_treatment_data.keys())) + """,
                datasets: [""")
        
        # Generate datasets for each treatment
        treatment_colors = {
            'WT': {'bg': 'rgba(54, 162, 235, 0.5)', 'border': 'rgba(54, 162, 235, 1)'},
            'WT-37': {'bg': 'rgba(255, 99, 132, 0.5)', 'border': 'rgba(255, 99, 132, 1)'},
            'WTA': {'bg': 'rgba(75, 192, 192, 0.5)', 'border': 'rgba(75, 192, 192, 1)'},
            'STC': {'bg': 'rgba(255, 206, 86, 0.5)', 'border': 'rgba(255, 206, 86, 1)'},
            'CAS': {'bg': 'rgba(153, 102, 255, 0.5)', 'border': 'rgba(153, 102, 255, 1)'}
        }
        
        datasets = []
        for treatment in sorted(treatment_totals.keys()):
            dataset = {
                'label': treatment,
                'data': [gene_treatment_data[gene].get(treatment, 0) for gene in gene_treatment_data.keys()],
                'backgroundColor': treatment_colors.get(treatment, {}).get('bg', f'rgba(100, 100, 100, 0.5)'),
                'borderColor': treatment_colors.get(treatment, {}).get('border', f'rgba(100, 100, 100, 1)'),
                'borderWidth': 1
            }
            datasets.append(dataset)
        
        f.write(',\n'.join([str(dataset).replace("'", '"') for dataset in datasets]))
        
        f.write("""
                ]
            },
            options: {
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Variant Count'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Gene'
                        }
                    }
                }
            }
        });
    </script>
</body>
</html>""")

def main():
    parser = argparse.ArgumentParser(description='Annotate variants based on JRIU IDs')
    parser.add_argument('--vcf', required=True, nargs='+', help='VCF file(s) to annotate')
    parser.add_argument('--genes-of-interest', required=True, help='Genes of interest file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # Define treatment groups
    treatment_map = {
        'WT-CTRL': 'WT',
        'WT-37-55-1': 'WT-37',
        'WT-37-55-2': 'WT-37',
        'WT-37-55-3': 'WT-37',
        'WTA-55-1': 'WTA',
        'WTA-55-2': 'WTA',
        'WTA-55-3': 'WTA',
        'STC-55-1': 'STC',
        'STC-55-2': 'STC',
        'STC-55-3': 'STC',
        'CAS-55-1': 'CAS',
        'CAS-55-2': 'CAS',
        'CAS-55-3': 'CAS',
        'STC-CTRL': 'STC-CTRL',
        'CAS-CTRL': 'CAS-CTRL'
    }
    
    # Create gene info lookup
    jriu_to_genes, sc_gene_map = create_gene_info(args.genes_of_interest)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process all VCF files
    result_files = process_all_vcfs(args.vcf, jriu_to_genes, sc_gene_map, args.output_dir, treatment_map)
    
    print("\nAnnotation Summary:")
    print(f"Summary report: {result_files['summary_file']}")
    print(f"Gene report: {result_files['gene_report_file']}")
    print(f"Treatment report: {result_files['treatment_report_file']}")
    print(f"HTML report: {result_files['html_report_file']}")

if __name__ == '__main__':
    main()