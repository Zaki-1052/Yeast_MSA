import os
import glob
import re
from Bio import SeqIO

# Directory containing GenBank files
GENBANK_DIR = os.path.join(os.path.dirname(__file__), '../../reference/w303_annotations')
OUTPUT_FILE = os.path.join(os.path.dirname(__file__), '../../reference/gene_mapping_full.tsv')

# Output columns
COLUMNS = [
    'w303_gene_id', 'locus_tag', 'sc_gene_id', 'std_gene_name', 'erg_name',
    'w303_scaffold', 'chromosome_id', 'start', 'end', 'strand', 'product'
]

STD_GENE_REGEX = re.compile(r'Y[A-P][LR][0-9]{3}[CW]')
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
    # Erg name: if gene name starts with 'ERG', else blank
    erg_name = gene_id if gene_id.upper().startswith('ERG') else ''
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
    with open(OUTPUT_FILE, 'w') as out:
        out.write('\t'.join(COLUMNS) + '\n')
        for gb_file in genbank_files:
            for record in SeqIO.parse(gb_file, 'genbank'):
                w303_scaffold, chromosome_id = get_scaffold_and_chromosome(record)
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
                        out.write('\t'.join(str(x) for x in row) + '\n')
    print(f"Gene mapping written to {OUTPUT_FILE}")

if __name__ == '__main__':
    main() 