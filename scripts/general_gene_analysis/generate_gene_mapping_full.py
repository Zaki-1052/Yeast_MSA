import os
import glob
import re
import pandas as pd
from Bio import SeqIO

# Directory containing GenBank files
GENBANK_DIR = os.path.join(os.path.dirname(__file__), '../../reference/w303_annotations')
OUTPUT_FILE = os.path.join(os.path.dirname(__file__), '../../reference/gene_mapping_full.tsv')
GENES_OF_INTEREST_FILE = os.path.join(os.path.dirname(__file__), '../../reference/genes_of_interest_mapping.tsv')

# Output columns
COLUMNS = [
    'w303_gene_id', 'locus_tag', 'sc_gene_id', 'std_gene_name', 'erg_name',
    'w303_scaffold', 'chromosome_id', 'start', 'end', 'strand', 'product'
]

STD_GENE_REGEX = re.compile(r'Y[A-P][LR][0-9]{3}[CW]')

# Load ergosterol genes from the genes of interest file if it exists
ERG_GENES = {}
if os.path.exists(GENES_OF_INTEREST_FILE):
    try:
        erg_df = pd.read_csv(GENES_OF_INTEREST_FILE, sep='\t')
        # Create a mapping from w303_gene_id to erg_name
        for _, row in erg_df.iterrows():
            if 'w303_gene_id' in row and 'erg_name' in row:
                ERG_GENES[row['w303_gene_id']] = row['erg_name']
        print(f"Loaded {len(ERG_GENES)} ergosterol genes from {GENES_OF_INTEREST_FILE}")
    except Exception as e:
        print(f"Error loading genes of interest: {e}")
        ERG_GENES = {}
SCAFFOLD_NUM_REGEX = re.compile(r'scaffold[ _]?([0-9]+)', re.IGNORECASE)

def extract_std_gene_name(cds_feature):
    # Try db_xref first
    db_xref = cds_feature.qualifiers.get('db_xref', [])
    for xref in db_xref:
        if STD_GENE_REGEX.match(xref.split(':')[-1]):
            return xref.split(':')[-1]
    # Try note
    for field in ['note', 'inference']:
        for val in cds_feature.qualifiers.get(field, []):
            match = STD_GENE_REGEX.search(val)
            if match:
                return match.group(0)
    return ''

def extract_gene_info(gene_feature, cds_feature, w303_scaffold, chromosome_id):
    gene_id = gene_feature.qualifiers.get('gene', [''])[0]
    locus_tag = ''
    product = ''
    sc_gene_id = ''
    std_gene_name = ''
    if cds_feature:
        locus_tag = cds_feature.qualifiers.get('locus_tag', [''])[0]
        product = cds_feature.qualifiers.get('product', [''])[0]
        db_xref = cds_feature.qualifiers.get('db_xref', [])
        for xref in db_xref:
            if xref.startswith('SGD:'):
                sc_gene_id = xref.split(':', 1)[1]
                break
        std_gene_name = extract_std_gene_name(cds_feature)
    
    # Get the locus tag from gene feature if not found in CDS
    if not locus_tag and 'locus_tag' in gene_feature.qualifiers:
        locus_tag = gene_feature.qualifiers['locus_tag'][0]
    
    # Set ergosterol gene name
    erg_name = ''
    # First check if it's in our pre-loaded ERG_GENES mapping
    if locus_tag in ERG_GENES:
        erg_name = ERG_GENES[locus_tag]
    elif gene_id in ERG_GENES:
        erg_name = ERG_GENES[gene_id]
    else:
        # Fallback to heuristic identification:
        # 1. Check gene_id directly (like ERG1, ERG2, etc.)
        # 2. Check standard name (some ERG genes might have different gene_id but standard name starts with ERG)
        # 3. Check product description for ergosterol pathway related terms
        if gene_id and gene_id.upper().startswith('ERG'):
            erg_name = gene_id
        elif std_gene_name and std_gene_name.upper().startswith('ERG'):
            erg_name = std_gene_name
        elif product and ('ergosterol' in product.lower() or 'sterol' in product.lower()):
            # If product mentions ergosterol pathway, use gene_id or std_gene_name
            erg_name = std_gene_name if std_gene_name else gene_id
        
    start = int(gene_feature.location.start) + 1  # GenBank is 0-based, inclusive
    end = int(gene_feature.location.end)
    strand = '+' if gene_feature.location.strand == 1 else '-'
    return gene_id, locus_tag, sc_gene_id, std_gene_name, erg_name, w303_scaffold, chromosome_id, start, end, strand, product

def get_scaffold_and_chromosome(record):
    chromosome_id = record.id
    w303_scaffold = ''
    for feature in record.features:
        if feature.type == 'source':
            note = feature.qualifiers.get('note', [])
            if note:
                chromosome_id = note[0]
                # Try to find scaffold number in note
                for n in note:
                    m = SCAFFOLD_NUM_REGEX.search(n)
                    if m:
                        w303_scaffold = f"w303_scaffold_{m.group(1)}"
                        break
            break
    return w303_scaffold, chromosome_id

def main():
    genbank_files = glob.glob(os.path.join(GENBANK_DIR, '*.genbank'))
    print(f"Found {len(genbank_files)} GenBank files in {GENBANK_DIR}")
    
    # Count statistics
    total_genes = 0
    erg_genes = 0
    scaffolds = set()
    chromosomes = set()
    
    with open(OUTPUT_FILE, 'w') as out:
        out.write('\t'.join(COLUMNS) + '\n')
        for gb_file in genbank_files:
            print(f"Processing {os.path.basename(gb_file)}...")
            gene_count = 0
            
            for record in SeqIO.parse(gb_file, 'genbank'):
                w303_scaffold, chromosome_id = get_scaffold_and_chromosome(record)
                scaffolds.add(w303_scaffold)
                chromosomes.add(chromosome_id)
                
                # Build a mapping from gene name to CDS feature
                gene_to_cds = {}
                for feature in record.features:
                    if feature.type == 'CDS':
                        gene_name = feature.qualifiers.get('gene', [''])[0]
                        if gene_name:
                            gene_to_cds[gene_name] = feature
                
                # Now process gene features
                for feature in record.features:
                    if feature.type == 'gene':
                        gene_id = feature.qualifiers.get('gene', [''])[0]
                        cds_feature = gene_to_cds.get(gene_id)
                        row = extract_gene_info(feature, cds_feature, w303_scaffold, chromosome_id)
                        
                        # Check if this is an ergosterol gene - only count if erg_name is not empty
                        if row[4] and len(row[4].strip()) > 0:  # erg_name field must not be empty
                            erg_genes += 1
                        
                        out.write('\t'.join(str(x) for x in row) + '\n')
                        gene_count += 1
                        total_genes += 1
            
            print(f"  Added {gene_count} genes from {os.path.basename(gb_file)}")
    
    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"Total genes processed: {total_genes}")
    print(f"Ergosterol pathway genes identified: {erg_genes}")
    print(f"Unique scaffolds: {len(scaffolds)}")
    print(f"Unique chromosome IDs: {len(chromosomes)}")
    print(f"\nGene mapping written to {OUTPUT_FILE}")

if __name__ == '__main__':
    main() 