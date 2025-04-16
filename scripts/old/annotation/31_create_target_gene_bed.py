#!/usr/bin/env python3

# File: scripts/annotation/31c_create_bed_only.py
# Purpose: Robustly create a BED file with coordinates for target genes (BED creation only)

import os
import re
import pandas as pd
from datetime import datetime
import logging
import subprocess # Kept for potential future use, but not strictly needed now

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def parse_genbank_for_locations(gbk_file_path, target_genes_sgd):
    """
    Parses a GenBank file to find locations of target SGD genes.

    Args:
        gbk_file_path (str): Path to the GenBank file.
        target_genes_sgd (list): List of SGD gene IDs (e.g., 'YHR190W').

    Returns:
        list: A list of dictionaries, each containing BED format info for found genes.
              Returns None if the file cannot be read.
    """
    if not os.path.exists(gbk_file_path):
        logging.error(f"GenBank file not found: {gbk_file_path}")
        return None

    logging.info(f"Parsing GenBank file: {gbk_file_path}...")
    target_genes_upper = {gene.upper() for gene in target_genes_sgd}
    gene_locations = []
    found_genes_set = set() # Track found upper-case SGD IDs to avoid duplicates

    current_locus = None
    current_feature_type = None # e.g., 'gene', 'CDS'
    feature_qualifiers = {}
    feature_location_str = None
    feature_is_complement = False
    line_num = 0 # Initialize line number

    try:
        with open(gbk_file_path, 'r', errors='ignore') as f:
            for line_num, line in enumerate(f, 1):
                # Detect start of a new entry
                if line.startswith('LOCUS'):
                    # Process the previous feature if it was potentially relevant
                    if current_feature_type and current_locus:
                       process_feature(current_locus, feature_location_str, feature_is_complement,
                                        feature_qualifiers, target_genes_upper, target_genes_sgd, found_genes_set, gene_locations)

                    # Start new entry
                    current_locus = line.split()[1] if len(line.split()) > 1 else None
                    if not current_locus:
                        logging.warning(f"Line {line_num}: Could not parse LOCUS ID from line: {line.strip()}")
                        continue # Skip if LOCUS ID is invalid
                    current_feature_type = None
                    feature_qualifiers = {}
                    feature_location_str = None
                    feature_is_complement = False
                    current_qualifier_key = None # Track the current qualifier being parsed
                    # logging.debug(f"Line {line_num}: Found LOCUS: {current_locus}")

                elif current_locus: # Only process lines if we are inside a valid LOCUS entry
                    # Detect start of a feature (gene or CDS)
                    match = re.match(r'\s{5}(gene|CDS)\s+(complement\()?((?:join|order)\(.*\)|<?\d+\.\.>?\d+)\)?', line) # Require 5 spaces indent
                    if match:
                        # Process the previous feature before starting a new one
                        if current_feature_type:
                           process_feature(current_locus, feature_location_str, feature_is_complement,
                                            feature_qualifiers, target_genes_upper, target_genes_sgd, found_genes_set, gene_locations)

                        # Start new feature
                        current_feature_type = match.group(1) # 'gene' or 'CDS'
                        feature_is_complement = bool(match.group(2)) # Check if 'complement(' exists
                        feature_location_str = match.group(3) # e.g., '123..456' or 'join(...)'
                        feature_qualifiers = {} # Reset qualifiers for the new feature
                        current_qualifier_key = None # Reset current qualifier key
                        # logging.debug(f"Line {line_num}: Found feature {current_feature_type} at {feature_location_str}")

                    # Detect qualifiers for the current feature (require 21 spaces indent)
                    elif current_feature_type and line.startswith(' ' * 21 + '/'):
                        stripped_line = line.strip()
                        # Match key="value" or key=value or /key
                        qualifier_match = re.match(r'/(?P<key>[a-zA-Z_]+)(?:="?(?P<value>.*?)"?)?$', stripped_line)
                        if qualifier_match:
                            current_qualifier_key = qualifier_match.group('key')
                            value = qualifier_match.group('value')
                            # Store value if present, otherwise store True for flags
                            feature_qualifiers[current_qualifier_key] = value.strip() if value is not None else True
                            # logging.debug(f"Line {line_num}: Found qualifier {current_qualifier_key}={str(feature_qualifiers[current_qualifier_key])[:50]}...")
                        else:
                           # Handle cases where qualifier value might be just text without quotes after key=
                           qualifier_match_no_quotes = re.match(r'/(?P<key>\w+)=(?P<value>.*)', stripped_line)
                           if qualifier_match_no_quotes:
                               current_qualifier_key = qualifier_match_no_quotes.group('key')
                               value = qualifier_match_no_quotes.group('value').strip()
                               feature_qualifiers[current_qualifier_key] = value
                           else:
                               # Might be a flag or malformed, store the key part
                               current_qualifier_key = stripped_line.split('=')[0][1:]
                               if current_qualifier_key:
                                   feature_qualifiers[current_qualifier_key] = True # Treat as flag
                               # logging.warning(f"Line {line_num}: Could not parse qualifier: {stripped_line}")


                    # Handle multi-line qualifiers (require 21 spaces indent, no starting '/')
                    elif current_feature_type and line.startswith(' ' * 21) and not line.strip().startswith('/') and current_qualifier_key and isinstance(feature_qualifiers.get(current_qualifier_key), str):
                         feature_qualifiers[current_qualifier_key] += " " + line.strip().replace('"', '')


            # Process the very last feature in the file
            if current_feature_type and current_locus:
                process_feature(current_locus, feature_location_str, feature_is_complement,
                                 feature_qualifiers, target_genes_upper, target_genes_sgd, found_genes_set, gene_locations)

    except Exception as e:
        logging.error(f"Error reading or parsing GenBank file around line {line_num}: {e}", exc_info=True)
        return None

    logging.info(f"Finished parsing. Found locations for {len(gene_locations)} target genes.")
    found_names_log = {entry['name'] for entry in gene_locations}
    logging.info(f"Found genes: {', '.join(sorted(list(found_names_log)))}")
    return gene_locations

def process_feature(locus, location_str, is_complement, qualifiers, target_genes_set_upper, target_genes_original, found_genes_set, results_list):
    """Checks if a parsed feature matches target genes and adds its location to results."""
    # Only process simple locations like 123..456 for now
    coords_match = re.match(r'<?(\d+)\.\.>?(\d+)', location_str if location_str else "")
    if not coords_match:
        # logging.debug(f"Skipping complex or missing location: {location_str}")
        return

    # Check if any qualifier value contains a target SGD ID
    matched_sgd_id_upper = None # Store the upper case version that matched

    # Combine relevant qualifier text for searching
    search_text = ' '.join([
        str(qualifiers.get('gene', '')),        # Ensure value is string
        str(qualifiers.get('product', '')),
        str(qualifiers.get('note', '')),
        str(qualifiers.get('locus_tag', '')),
        str(qualifiers.get('inference', ''))
    ]).upper() # Convert the whole thing to upper case once

    for sgd_gene_upper in target_genes_set_upper:
        # More robust check - look for SGD ID possibly surrounded by non-word characters or as whole word
        # Check inference tag specifically for :SGDID pattern
        if f":{sgd_gene_upper}" in qualifiers.get('inference', '').upper():
             matched_sgd_id_upper = sgd_gene_upper
             break
        # Check other common tags for the SGD ID as a word or part of text
        elif re.search(rf'\b{sgd_gene_upper}\b', search_text):
             matched_sgd_id_upper = sgd_gene_upper
             break

    # If a match was found and we haven't already added this specific SGD ID
    if matched_sgd_id_upper and matched_sgd_id_upper not in found_genes_set:
        try:
            start = int(coords_match.group(1))
            end = int(coords_match.group(2))
            strand = "-" if is_complement else "+"

            # Convert to 0-based start for BED
            bed_start = start - 1

            # Find original case SGD ID from the input list
            original_sgd_id = next((g for g in target_genes_original if g.upper() == matched_sgd_id_upper), matched_sgd_id_upper) # Fallback to upper if not found

            results_list.append({
                'chrom': locus,
                'chromStart': bed_start,
                'chromEnd': end,
                'name': original_sgd_id, # Use the original SGD ID
                'score': '0',
                'strand': strand
            })
            found_genes_set.add(matched_sgd_id_upper) # Mark this upper-case version as found
            # Log the successful mapping including which feature type provided the coordinates
            # feature_type_info = qualifiers.get("feature_type", "Unknown") # Assuming you store this if needed
            # logging.info(f"  -> Mapped {original_sgd_id} to {locus}:{bed_start}-{end} (from {feature_type_info})")

        except ValueError as e:
            logging.warning(f"Coordinate conversion error for '{location_str}': {e}")
        except Exception as e:
             logging.error(f"Unexpected error processing feature for {original_sgd_id}: {e}", exc_info=True)


def main():
    logging.info("=== Creating Target Gene BED File (Robust Method v2) ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")

    # Define directories
    snpeff_dir = "/Users/zakiralibhai/snpEff"
    gene_dir = "annotation/genes_of_interest"
    output_dir = "annotation/gene_coordinates"

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Read target genes
    target_genes_file = f"{gene_dir}/target_genes.txt"
    try:
        with open(target_genes_file, 'r') as f:
            target_genes = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        logging.error(f"Target genes file not found: {target_genes_file}")
        return

    logging.info(f"Extracting locations for {len(target_genes)} target genes...")

    # Path to GenBank file
    genes_gbk = f"{snpeff_dir}/data/w303/genes.gbk"

    # Extract locations using the robust parser
    bed_entries = parse_genbank_for_locations(genes_gbk, target_genes)

    if bed_entries is None:
        logging.error("Failed to parse GenBank file.")
        return

    # Report missing genes
    found_gene_names = {entry['name'] for entry in bed_entries}
    missing_genes = [gene for gene in target_genes if gene not in found_gene_names]

    # Create BED file
    bed_file_path = f"{output_dir}/target_genes_robust.bed"
    written_count = 0
    try:
        # Sort entries by chromosome and start position before writing
        bed_entries.sort(key=lambda x: (x['chrom'], x['chromStart']))
        with open(bed_file_path, 'w') as bed_file:
            # Write header
            bed_file.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n")
            # Write entries
            for entry in bed_entries:
                bed_file.write(f"{entry['chrom']}\t{entry['chromStart']}\t{entry['chromEnd']}\t"
                               f"{entry['name']}\t{entry['score']}\t{entry['strand']}\n")
                written_count += 1
        logging.info(f"Created BED file: {bed_file_path}")
        logging.info(f"Successfully mapped and wrote {written_count} out of {len(target_genes)} target genes.")

    except IOError as e:
         logging.error(f"Failed to write BED file: {e}")
         return # Stop if we can't write the BED file

    if missing_genes:
        logging.warning("Could not find locations for the following genes:")
        for gene in missing_genes:
            logging.warning(f"  - {gene}")
        logging.warning("These genes will be excluded from the BED file.")

    # Removed the shell script generation part

    print("\n=== Robust BED File Creation Complete ===")
    print(f"BED file saved to: {bed_file_path}")
    print("Please review the BED file before proceeding.")
    print("You can now run 'scripts/annotation/32_intersect_variants.sh'")


if __name__ == "__main__":
    main()