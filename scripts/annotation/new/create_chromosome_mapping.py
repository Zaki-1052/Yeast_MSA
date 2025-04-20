# File: scripts/annotation/new/create_chromosome_mapping.py

import re

def clean_header(header):
    """Extract just the chromosome ID from a FASTA header."""
    if header.startswith('>CM'):
        # Extract CM007XXX.1
        match = re.search(r'(>)(CM\d+\.\d+)', header)
        if match:
            return match.group(2)
    elif header.startswith('>LYZE'):
        # Extract LYZE ID
        match = re.search(r'(>)(LYZE\d+\.\d+)', header)
        if match:
            return match.group(2)
    elif header.startswith('>'):
        # For Roman numeral headers, just return without '>'
        return header[1:]
    return header

def create_mapping_table():
    """Create a cleaned mapping table from header_mapping.txt."""
    mapping = {}
    reverse_mapping = {}
    
    with open('reference/header_mapping.txt', 'r') as f:
        # Skip header line
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                orig_header = parts[0]
                ygap_header = parts[1]
                
                # Clean headers
                orig_clean = clean_header(orig_header)
                ygap_clean = clean_header(ygap_header)
                
                mapping[orig_clean] = ygap_clean
                reverse_mapping[ygap_clean] = orig_clean
    
    # Write cleaned mapping files
    with open('reference/chromosome_mapping.txt', 'w') as f:
        f.write("CM_ID\tYGAP_ID\n")
        for cm_id, ygap_id in mapping.items():
            f.write(f"{cm_id}\t{ygap_id}\n")
    
    with open('reference/chromosome_mapping_reverse.txt', 'w') as f:
        f.write("YGAP_ID\tCM_ID\n")
        for ygap_id, cm_id in reverse_mapping.items():
            f.write(f"{ygap_id}\t{cm_id}\n")
    
    print(f"Created mapping for {len(mapping)} chromosomes")
    return mapping, reverse_mapping

if __name__ == "__main__":
    mapping, reverse_mapping = create_mapping_table()
    
    # Print sample of the mapping
    print("\nSample of CM → YGAP mapping:")
    items = list(mapping.items())
    for i in range(min(5, len(items))):
        print(f"{items[i][0]} → {items[i][1]}")
    
    print("\nSample of YGAP → CM mapping:")
    items = list(reverse_mapping.items())
    for i in range(min(5, len(items))):
        print(f"{items[i][0]} → {items[i][1]}")