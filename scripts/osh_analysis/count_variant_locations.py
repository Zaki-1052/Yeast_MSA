import csv

# Path to your file
filename = "results/osh_analysis/osh_variants.tsv"

# Initialize counters
total = 0
within_5kb = 0
upstream = 0
downstream = 0
within_gene = 0

# Keywords that might indicate a variant is within the gene
within_gene_keywords = ['genic', 'exonic', 'within_gene', 'intron', 'exon']

with open(filename, newline='') as tsvfile:
    reader = csv.DictReader(tsvfile, delimiter='\t')
    for row in reader:
        total += 1
        # Distance check
        try:
            distance = float(row['distance'])
        except ValueError:
            continue  # skip if distance is not a number
        if distance <= 5000:
            within_5kb += 1
        # Location check
        location = row['location'].lower()
        if location == 'upstream':
            upstream += 1
        elif location == 'downstream':
            downstream += 1
        # Within gene check
        if any(key in location for key in within_gene_keywords):
            within_gene += 1

# Print results
print(f"Total variants: {total}")
print(f"Within 5kb: {within_5kb} ({within_5kb/total*100:.1f}%)")
print(f"Upstream: {upstream} ({upstream/total*100:.1f}%)")
print(f"Downstream: {downstream} ({downstream/total*100:.1f}%)")
print(f"Within gene: {within_gene} ({within_gene/total*100:.1f}%)")
