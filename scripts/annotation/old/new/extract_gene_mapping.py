# File: scripts/annotation/new/extract_gene_mapping.py

import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_gene_mapping():
    """Extract mapping between standard gene IDs and YGAP IDs from GenBank files."""
    input_dir = os.path.expanduser("~/snpEff/data/w303_cm/genes")
    
    # List of genes of interest
    genes_of_interest = [
        "YHR190W", "YGR175C", "YHR072W", "YHR007C", "YNL280C", 
        "YGR060W", "YML008C", "YMR202W", "YLR056W", "YMR015C", "YGL012W"
    ]
    
    # Dictionary to store gene mappings
    gene_mappings = {}
    
    # Dictionary to track genes of interest
    goi_mappings = {}
    
    # Dictionary to track chromosomes for each gene
    gene_chromosome = {}
    
    # Get list of all GenBank files
    gbk_files = [f for f in os.listdir(input_dir) if f.endswith('.gbk')]
    print(f"Processing {len(gbk_files)} GenBank files")
    
    # Track statistical information
    total_genes = 0
    genes_with_standard_id = 0
    
    # Process each GenBank file
    for gb_file in gbk_files:
        input_path = os.path.join(input_dir, gb_file)
        chromosome = gb_file.replace('.gbk', '')
        
        try:
            # Parse the GenBank file
            record = SeqIO.read(input_path, "genbank")
            
            # Get chromosome/scaffold info from source feature
            scaffold_name = None
            roman_numeral = None
            for feature in record.features:
                if feature.type == "source":
                    for note in feature.qualifiers.get('note', []):
                        if note in ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 
                                   'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'MT', 'p2micron',
                                   'XII_Seq19', 'XII_Seq20', 'XII_Seq21']:
                            roman_numeral = note
                        if note.startswith('scaffold'):
                            scaffold_name = note
            
            # Process gene features
            for feature in record.features:
                if feature.type == "gene":
                    total_genes += 1
                    
                    # Extract YGAP gene ID
                    ygap_id = None
                    if 'gene' in feature.qualifiers:
                        ygap_id = feature.qualifiers['gene'][0]
                    
                    if not ygap_id:
                        continue
                    
                    # Extract standard gene ID, if present
                    standard_id = None
                    product = "hypothetical protein"
                    
                    # Check inference qualifier for standard ID
                    if 'inference' in feature.qualifiers:
                        for inf in feature.qualifiers['inference']:
                            match = re.search(r'INSDC:(Y[A-Z0-9]+)', inf)
                            if match:
                                standard_id = match.group(1)
                                genes_with_standard_id += 1
                    
                    # Check note qualifier for standard ID
                    if not standard_id and 'note' in feature.qualifiers:
                        for note in feature.qualifiers['note']:
                            match = re.search(r'similar to Saccharomyces cerevisiae (\w+) \((Y[A-Z0-9]+)\)', note)
                            if match:
                                gene_name = match.group(1)
                                standard_id = match.group(2)
                                product = gene_name
                                genes_with_standard_id += 1
                    
                    # Store the mapping
                    if standard_id:
                        gene_mappings[ygap_id] = {
                            'standard_id': standard_id,
                            'chromosome': chromosome,
                            'roman_numeral': roman_numeral,
                            'scaffold': scaffold_name,
                            'location': f"{feature.location.start}..{feature.location.end}",
                            'product': product
                        }
                        
                        # Track chromosome for this gene
                        gene_chromosome[standard_id] = {
                            'chromosome': chromosome,
                            'roman_numeral': roman_numeral,
                            'scaffold': scaffold_name
                        }
                        
                        # Check if this is a gene of interest
                        if standard_id in genes_of_interest:
                            goi_mappings[standard_id] = {
                                'ygap_id': ygap_id,
                                'chromosome': chromosome,
                                'roman_numeral': roman_numeral,
                                'scaffold': scaffold_name,
                                'location': f"{feature.location.start}..{feature.location.end}",
                                'product': product
                            }
        
        except Exception as e:
            print(f"Error processing {gb_file}: {str(e)}")
    
    # Write gene mapping table to file
    output_mapping = "reference/gene_mapping.tsv"
    with open(output_mapping, 'w') as f:
        f.write("YGAP_ID\tStandard_ID\tChromosome\tRoman_Numeral\tScaffold\tLocation\tProduct\n")
        for ygap_id, info in gene_mappings.items():
            f.write(f"{ygap_id}\t{info['standard_id']}\t{info['chromosome']}\t{info['roman_numeral']}\t{info['scaffold']}\t{info['location']}\t{info['product']}\n")
    
    # Write genes of interest mapping to file
    goi_output = "reference/genes_of_interest_mapping.tsv"
    with open(goi_output, 'w') as f:
        f.write("Standard_ID\tYGAP_ID\tChromosome\tRoman_Numeral\tScaffold\tLocation\tProduct\n")
        for standard_id in genes_of_interest:
            if standard_id in goi_mappings:
                info = goi_mappings[standard_id]
                f.write(f"{standard_id}\t{info['ygap_id']}\t{info['chromosome']}\t{info['roman_numeral']}\t{info['scaffold']}\t{info['location']}\t{info['product']}\n")
            else:
                f.write(f"{standard_id}\tNOT_FOUND\tN/A\tN/A\tN/A\tN/A\tN/A\n")
    
    # Write chromosome summary to file
    chromosome_output = "reference/chromosome_summary.tsv"
    with open(chromosome_output, 'w') as f:
        f.write("Chromosome\tRoman_Numeral\tScaffold\tIs_Fragment\tCircular\n")
        for gb_file in gbk_files:
            chromosome = gb_file.replace('.gbk', '')
            input_path = os.path.join(input_dir, gb_file)
            
            # Get chromosome/scaffold info
            roman_numeral = "unknown"
            scaffold_name = "unknown"
            is_fragment = "No"
            is_circular = "No"
            
            try:
                record = SeqIO.read(input_path, "genbank")
                for feature in record.features:
                    if feature.type == "source":
                        for note in feature.qualifiers.get('note', []):
                            if note in ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 
                                       'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'MT', 'p2micron',
                                       'XII_Seq19', 'XII_Seq20', 'XII_Seq21']:
                                roman_numeral = note
                            if note.startswith('scaffold'):
                                scaffold_name = note
                
                # Set fragment and circular flags
                if roman_numeral in ['XII_Seq19', 'XII_Seq20', 'XII_Seq21']:
                    is_fragment = "Yes"
                
                if roman_numeral in ['MT', 'p2micron']:
                    is_circular = "Yes"
            
            except Exception as e:
                print(f"Error processing {gb_file} for chromosome summary: {str(e)}")
            
            f.write(f"{chromosome}\t{roman_numeral}\t{scaffold_name}\t{is_fragment}\t{is_circular}\n")
    
    # Print summary statistics
    print(f"Total genes processed: {total_genes}")
    print(f"Genes with standard IDs: {genes_with_standard_id}")
    print(f"Genes of interest found: {len(goi_mappings)} out of {len(genes_of_interest)}")
    print(f"Genes of interest not found: {len(genes_of_interest) - len(goi_mappings)}")
    print(f"Gene mapping table written to: {output_mapping}")
    print(f"Genes of interest mapping written to: {goi_output}")
    print(f"Chromosome summary written to: {chromosome_output}")
    
    return gene_mappings, goi_mappings

if __name__ == "__main__":
    extract_gene_mapping()