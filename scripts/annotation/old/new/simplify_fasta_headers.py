#!/usr/bin/env python3

import re
import os

# Input and output file paths
input_fasta = "reference/w303_chromosomal.fasta"
output_fasta = "reference/w303_for_ygap.fasta"
mapping_file = "reference/header_mapping.txt"

# Dictionary to store the mapping
header_mapping = {}

# Process the FASTA file
with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile, open(mapping_file, 'w') as mapfile:
    # Write header for mapping file
    mapfile.write("Original_Header\tYGAP_Header\n")
    
    current_header = ""
    in_sequence = False
    
    for line in infile:
        line = line.strip()
        
        # If this is a header line
        if line.startswith('>'):
            in_sequence = True
            current_header = line
            
            # Handle the unlocalized chromosome XII sequences specially
            if "chromosome XII Seq19" in line:
                new_header = ">XII_Seq19"
            elif "chromosome XII Seq20" in line:
                new_header = ">XII_Seq20"
            elif "chromosome XII Seq21" in line:
                new_header = ">XII_Seq21"
            # Handle regular chromosomes
            elif "chromosome I," in line:
                new_header = ">I"
            elif "chromosome II," in line:
                new_header = ">II"
            elif "chromosome III," in line:
                new_header = ">III"
            elif "chromosome IV," in line:
                new_header = ">IV"
            elif "chromosome V," in line:
                new_header = ">V"
            elif "chromosome VI," in line:
                new_header = ">VI"
            elif "chromosome VII," in line:
                new_header = ">VII"
            elif "chromosome VIII," in line:
                new_header = ">VIII"
            elif "chromosome IX," in line:
                new_header = ">IX"
            elif "chromosome X," in line:
                new_header = ">X"
            elif "chromosome XI," in line:
                new_header = ">XI"
            elif "chromosome XII," in line and "Seq" not in line:  # Regular chromosome XII
                new_header = ">XII"
            elif "chromosome XIII," in line:
                new_header = ">XIII"
            elif "chromosome XIV," in line:
                new_header = ">XIV"
            elif "chromosome XV," in line:
                new_header = ">XV"
            elif "chromosome XVI," in line:
                new_header = ">XVI"
            elif "mitochondrion" in line:
                new_header = ">MT"
            elif "plasmid p2-micron" in line:
                new_header = ">p2micron"
            else:
                # Handle any truly unknown sequences with a more descriptive name
                match = re.search(r'>(\w+)', line)
                if match:
                    seq_id = match.group(1)
                    new_header = f">{seq_id}"
                else:
                    new_header = ">unknown_" + line[1:10]
            
            # Store the mapping
            header_mapping[current_header] = new_header
            
            # Write the new header
            outfile.write(new_header + "\n")
            
            # Write to mapping file
            mapfile.write(f"{current_header}\t{new_header}\n")
        
        # If this is a sequence line
        elif in_sequence:
            # Write the sequence line unchanged
            outfile.write(line + "\n")

print(f"Converted FASTA file saved to {output_fasta}")
print(f"Header mapping saved to {mapping_file}")