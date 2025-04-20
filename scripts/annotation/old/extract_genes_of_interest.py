#!/usr/bin/env python3
"""
extract_genes_of_interest.py - Extract detailed information and sequences for genes of interest
"""

import os
import sys
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def extract_gene_details(genbank_dir, genes_of_interest_file, output_dir):
    """
    Extract detailed information about genes of interest
    
    Args:
        genbank_dir (str): Directory containing GenBank files
        genes_of_interest_file (str): TSV file with genes of interest info
        output_dir (str): Directory to write output files
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read genes of interest from TSV file
    genes_of_interest = []
    with open(genes_of_interest_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            genes_of_interest.append(row)
    
    print(f"Processing {len(genes_of_interest)} gene matches...")
    
    # Create dictionaries to store GenBank file paths by scaffold ID
    scaffold_to_file = {}
    
    # First scan to map scaffold IDs to file paths
    for filename in os.listdir(genbank_dir):
        if not filename.endswith(('.genbank', '.gb', '.gbk')):
            continue
            
        filepath = os.path.join(genbank_dir, filename)
        
        try:
            record = SeqIO.read(filepath, "genbank")
            scaffold_id = record.id
            
            # Handle "unknown" scaffold IDs
            if scaffold_id == "unknown" and "_scaffold_" in filename:
                scaffold_num = filename.split("_scaffold_")[1].split(".")[0]
                scaffold_id = f"scaffold_{scaffold_num}"
            
            scaffold_to_file[scaffold_id] = filepath
        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
    
    # Extract gene details and sequences
    genes_extracted = 0
    gene_details = []
    
    # Create FASTA file for gene sequences
    fasta_file = os.path.join(output_dir, "genes_of_interest.fasta")
    fasta_records = []
    
    # Process each gene of interest
    for gene in genes_of_interest:
        sc_gene_id = gene['sc_gene_id']
        w303_gene_id = gene['w303_gene_id']
        scaffold_id = gene['scaffold_id']
        start = int(gene['start'])
        end = int(gene['end'])
        strand = gene['strand']
        sc_gene_name = gene['sc_gene_name']
        
        if scaffold_id not in scaffold_to_file:
            print(f"Warning: Could not find GenBank file for {scaffold_id} containing {w303_gene_id}")
            continue
        
        genbank_file = scaffold_to_file[scaffold_id]
        
        try:
            # Parse the GenBank file
            record = SeqIO.read(genbank_file, "genbank")
            
            # Extract the gene sequence
            gene_seq = record.seq[start:end]
            
            # Reverse complement if on negative strand
            if strand == '-':
                gene_seq = gene_seq.reverse_complement()
            
            # Find the corresponding feature in the GenBank file
            gene_feature = None
            cds_feature = None
            
            for feature in record.features:
                if feature.type == "gene" and any(w303_gene_id in q for q in feature.qualifiers.get('locus_tag', [])):
                    gene_feature = feature
                elif feature.type == "CDS" and any(w303_gene_id in q for q in feature.qualifiers.get('locus_tag', [])):
                    cds_feature = feature
            
            # Extract notes and product information
            notes = ""
            product = ""
            translation = ""
            
            if gene_feature:
                notes = "; ".join(gene_feature.qualifiers.get('note', []))
            
            if cds_feature:
                product = cds_feature.qualifiers.get('product', [''])[0]
                if 'translation' in cds_feature.qualifiers:
                    translation = cds_feature.qualifiers['translation'][0]
            
            # Create a FASTA record
            fasta_id = f"{sc_gene_id}|{sc_gene_name}|{w303_gene_id}"
            fasta_desc = f"{scaffold_id}:{start}-{end}:{strand}"
            seq_record = SeqRecord(gene_seq, id=fasta_id, description=fasta_desc)
            fasta_records.append(seq_record)
            
            # Add to gene details
            gene_details.append({
                'sc_gene_id': sc_gene_id,
                'sc_gene_name': sc_gene_name,
                'w303_gene_id': w303_gene_id,
                'scaffold_id': scaffold_id,
                'jriu_id': gene['jriu_id'],
                'start': start,
                'end': end,
                'strand': strand,
                'length': len(gene_seq),
                'product': product,
                'notes': notes,
                'translation_length': len(translation) if translation else 0,
                'gc_content': (gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq) if len(gene_seq) > 0 else 0
            })
            
            genes_extracted += 1
            
        except Exception as e:
            print(f"Error extracting {w303_gene_id} from {genbank_file}: {str(e)}")
    
    # Write FASTA sequences to file
    if fasta_records:
        SeqIO.write(fasta_records, fasta_file, "fasta")
        print(f"Wrote {len(fasta_records)} gene sequences to {fasta_file}")
    
    # Write detailed gene info to file
    if gene_details:
        gene_details_file = os.path.join(output_dir, "genes_of_interest_details.tsv")
        with open(gene_details_file, 'w') as f:
            # Write header
            fieldnames = ['sc_gene_id', 'sc_gene_name', 'w303_gene_id', 'scaffold_id', 'jriu_id', 
                          'start', 'end', 'strand', 'length', 'product', 'notes', 
                          'translation_length', 'gc_content']
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            # Write data
            for gene in gene_details:
                writer.writerow(gene)
        
        print(f"Wrote detailed information for {len(gene_details)} genes to {gene_details_file}")
    
    # Generate HTML report
    html_report_file = os.path.join(output_dir, "genes_of_interest_report.html")
    
    with open(html_report_file, 'w') as f:
        f.write("""<!DOCTYPE html>
<html>
<head>
    <title>Genes of Interest - Yeast Genome Annotation Project</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2 { color: #333366; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .gene-info { margin-bottom: 30px; border: 1px solid #ddd; padding: 15px; border-radius: 5px; }
        .gene-header { background-color: #eeeeff; padding: 10px; margin-bottom: 10px; }
        .note { color: #666; font-style: italic; }
        .sequence { font-family: monospace; word-break: break-all; margin-top: 10px; background-color: #f5f5f5; padding: 10px; }
    </style>
</head>
<body>
    <h1>Genes of Interest - Yeast Genome Annotation Project</h1>
    <p>Generated on: """ + f"{os.popen('date').read().strip()}" + """</p>
    
    <h2>Summary</h2>
    <p>Total genes identified: """ + f"{len(gene_details)}" + """</p>
    
    <table>
        <tr>
            <th>S. cerevisiae ID</th>
            <th>Gene Name</th>
            <th>W303 ID</th>
            <th>Location</th>
            <th>Length (bp)</th>
            <th>GC Content</th>
            <th>Product</th>
        </tr>
""")
        
        # Add table rows for each gene
        for gene in sorted(gene_details, key=lambda g: g['sc_gene_name']):
            location = f"{gene['scaffold_id']}:{gene['start']}-{gene['end']}:{gene['strand']}"
            f.write(f"""        <tr>
            <td>{gene['sc_gene_id']}</td>
            <td>{gene['sc_gene_name']}</td>
            <td>{gene['w303_gene_id']}</td>
            <td>{location}</td>
            <td>{gene['length']}</td>
            <td>{gene['gc_content']:.2f}</td>
            <td>{gene['product']}</td>
        </tr>
""")
        
        f.write("""    </table>
    
    <h2>Detailed Gene Information</h2>
""")
        
        # Add detailed sections for each gene
        for gene in sorted(gene_details, key=lambda g: g['sc_gene_name']):
            location = f"{gene['scaffold_id']}:{gene['start']}-{gene['end']}:{gene['strand']}"
            
            # Get the gene sequence from our FASTA records
            gene_seq = ""
            for record in fasta_records:
                if record.id.startswith(f"{gene['sc_gene_id']}|"):
                    gene_seq = str(record.seq)
                    break
            
            f.write(f"""    <div class="gene-info">
        <div class="gene-header">
            <h3>{gene['sc_gene_name']} ({gene['sc_gene_id']})</h3>
        </div>
        
        <p><strong>W303 ID:</strong> {gene['w303_gene_id']}</p>
        <p><strong>Location:</strong> {location}</p>
        <p><strong>Length:</strong> {gene['length']} bp</p>
        <p><strong>GC Content:</strong> {gene['gc_content']:.2f}</p>
        <p><strong>Product:</strong> {gene['product']}</p>
        
        <p class="note"><strong>Notes:</strong> {gene['notes']}</p>
        
        <p><strong>Sequence:</strong></p>
        <div class="sequence">{gene_seq}</div>
    </div>
""")
        
        f.write("""</body>
</html>""")
    
    print(f"Generated HTML report at {html_report_file}")
    return genes_extracted

def main():
    if len(sys.argv) < 3:
        print("Usage: python extract_genes_of_interest.py <genbank_directory> <genes_of_interest_file> [output_directory]")
        sys.exit(1)
    
    genbank_dir = sys.argv[1]
    genes_of_interest_file = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "annotation/genes_of_interest"
    
    if not os.path.isdir(genbank_dir):
        print(f"Error: {genbank_dir} is not a directory")
        sys.exit(1)
    
    if not os.path.isfile(genes_of_interest_file):
        print(f"Error: {genes_of_interest_file} does not exist")
        sys.exit(1)
    
    num_extracted = extract_gene_details(genbank_dir, genes_of_interest_file, output_dir)
    print(f"\nExtracted detailed information for {num_extracted} genes of interest")
    print(f"Output files are available in {output_dir}")

if __name__ == "__main__":
    main()