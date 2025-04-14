#!/usr/bin/env python3

# File: scripts/annotation/27_create_gene_mapping_final.py
# Purpose: Create final gene mapping and extract variants in target genes

import os
import subprocess
import pandas as pd
from datetime import datetime
import re
import gzip

def main():
    print("=== Creating Final Gene Mapping and Extracting Variants ===")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # Define directories
    output_dir = "annotation/gene_variants_final"
    gene_dir = "annotation/genes_of_interest"
    vcf_source = "annotation/vcf_ready"
    
    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/by_gene", exist_ok=True)
    
    # Create the gene mapping
    gene_mapping = [
        {"SGD_Gene": "YHR190W", "W303_ID": "W303_0BY00190", "Gene": "W3030BY00190", 
         "Coordinates": "19026..20360", "Description": "similar to ERG9"},
        {"SGD_Gene": "YGR175C", "W303_ID": "W303_0AJ00440", "Gene": "W3030AJ00440", 
         "Coordinates": "complement(71492..72982)", "Description": "similar to ERG1"},
        {"SGD_Gene": "YHR072W", "W303_ID": "W303_0CB00150", "Gene": "W3030CB00150", 
         "Coordinates": "9160..11355", "Description": "similar to ERG7"},
        {"SGD_Gene": "YHR007C", "W303_ID": "W303_0EI00110", "Gene": "W3030EI00110", 
         "Coordinates": "complement(3375..4967)", "Description": "similar to ERG11"},
        {"SGD_Gene": "YNL280C", "W303_ID": "W303_0O00140", "Gene": "W3030O00140", 
         "Coordinates": "complement(6598..7914)", "Description": "similar to ERG24"},
        {"SGD_Gene": "YGR060W", "W303_ID": "W303_0W00270", "Gene": "W3030W00270", 
         "Coordinates": "35900..36829", "Description": "similar to ERG25"},
        {"SGD_Gene": "YML008C", "W303_ID": "W303_0S00270", "Gene": "W3030S00270", 
         "Coordinates": "complement(28556..29707)", "Description": "similar to ERG6"},
        {"SGD_Gene": "YMR202W", "W303_ID": "W303_0AD00260", "Gene": "W3030AD00260", 
         "Coordinates": "33370..34038", "Description": "similar to ERG2"},
        {"SGD_Gene": "YLR056W", "W303_ID": "W303_0E01010", "Gene": "W3030E01010", 
         "Coordinates": "complement(196480..197577)", "Description": "similar to ERG3"},
        {"SGD_Gene": "YMR015C", "W303_ID": "W303_0S00490", "Gene": "W3030S00490", 
         "Coordinates": "complement(77542..79158)", "Description": "similar to ERG5"},
        {"SGD_Gene": "YGL012W", "W303_ID": "W303_0Y00390", "Gene": "W3030Y00390", 
         "Coordinates": "complement(55892..57313)", "Description": "similar to ERG4"}
    ]
    
    # Create mapping dataframe
    mapping_df = pd.DataFrame(gene_mapping)
    mapping_file = f"{output_dir}/gene_mapping.tsv"
    mapping_df.to_csv(mapping_file, sep='\t', index=False)
    
    print(f"Created gene mapping file with {len(mapping_df)} genes")
    print(f"Saved to: {mapping_file}")
    print("")
    
    # Find VCF files
    vcf_files = []
    for root, dirs, files in os.walk(vcf_source):
        for file in files:
            if file.endswith('.sorted.vcf.gz'):
                vcf_files.append(os.path.join(root, file))
    
    if not vcf_files:
        print(f"ERROR: No VCF files found in {vcf_source}")
        return
    
    print(f"Found {len(vcf_files)} VCF files")
    print("")
    
    # Extract variants for each gene
    all_variants = []
    gene_variants = {gene["SGD_Gene"]: [] for gene in gene_mapping}
    
    for vcf_file in vcf_files:
        sample = os.path.basename(vcf_file).replace('.sorted.vcf.gz', '')
        print(f"Processing {sample}...")
        
        # Count variants by gene
        counts = {gene["SGD_Gene"]: 0 for gene in gene_mapping}
        
        # Process each VCF file
        try:
            # Use bcftools to extract all variants
            cmd = f"bcftools view -H {vcf_file} | head -10000 > {output_dir}/temp_variants.txt"
            subprocess.run(cmd, shell=True)
            
            with open(f"{output_dir}/temp_variants.txt", 'r') as f:
                for line in f:
                    for gene_info in gene_mapping:
                        sgd_gene = gene_info["SGD_Gene"]
                        w303_id = gene_info["W303_ID"]
                        gene_id = gene_info["Gene"]
                        
                        # Check if the variant line contains our gene identifier
                        if sgd_gene in line or w303_id in line or gene_id in line:
                            # Parse fields
                            fields = line.strip().split('\t')
                            if len(fields) < 8:
                                continue
                            
                            chrom = fields[0]
                            pos = fields[1]
                            ref = fields[3]
                            alt = fields[4]
                            info = fields[7]
                            
                            # Create a record
                            record = {
                                "Sample": sample,
                                "SGD_Gene": sgd_gene,
                                "W303_ID": w303_id,
                                "Gene": gene_id,
                                "Chromosome": chrom,
                                "Position": pos,
                                "Ref": ref,
                                "Alt": alt,
                                "Info": info
                            }
                            
                            # Add to results
                            all_variants.append(record)
                            gene_variants[sgd_gene].append(record)
                            counts[sgd_gene] += 1
            
            # Delete temp file
            os.remove(f"{output_dir}/temp_variants.txt")
            
            # Report counts
            for gene, count in counts.items():
                if count > 0:
                    print(f"  Found {count} variants affecting {gene}")
        except Exception as e:
            print(f"  Error processing {vcf_file}: {e}")
    
    # Save results
    all_df = pd.DataFrame(all_variants)
    if not all_df.empty:
        all_df.to_csv(f"{output_dir}/all_target_variants.tsv", sep='\t', index=False)
        print(f"\nSaved {len(all_df)} total variants to {output_dir}/all_target_variants.tsv")
        
        # Save by gene
        for sgd_gene, variants in gene_variants.items():
            if variants:
                gene_df = pd.DataFrame(variants)
                gene_df.to_csv(f"{output_dir}/by_gene/{sgd_gene}_variants.tsv", sep='\t', index=False)
                print(f"Saved {len(gene_df)} variants for {sgd_gene}")
    else:
        print("\nNo variants found in target genes")
    
    # Create a fallback search script
    fallback_script = f"scripts/annotation/28_search_gene_variants_directly.sh"
    with open(fallback_script, 'w') as f:
        f.write("""#!/bin/bash

# File: scripts/annotation/28_search_gene_variants_directly.sh
# Purpose: Search directly for gene variants using grep

echo "=== Direct Gene Variant Search ==="
echo "Date: $(date)"
echo ""

# Define directories
OUTPUT_DIR="annotation/gene_variants_grep"
VCF_DIR="annotation/vcf_ready"

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/by_gene"

# Define genes to search for
GENES=(
  "YHR190W" "W303_0BY00190" "W3030BY00190" "ERG9"
  "YGR175C" "W303_0AJ00440" "W3030AJ00440" "ERG1"
  "YHR072W" "W303_0CB00150" "W3030CB00150" "ERG7"
  "YHR007C" "W303_0EI00110" "W3030EI00110" "ERG11"
  "YNL280C" "W303_0O00140" "W3030O00140" "ERG24"
  "YGR060W" "W303_0W00270" "W3030W00270" "ERG25"
  "YML008C" "W303_0S00270" "W3030S00270" "ERG6"
  "YMR202W" "W303_0AD00260" "W3030AD00260" "ERG2"
  "YLR056W" "W303_0E01010" "W3030E01010" "ERG3"
  "YMR015C" "W303_0S00490" "W3030S00490" "ERG5"
  "YGL012W" "W303_0Y00390" "W3030Y00390" "ERG4"
)

# Find VCF files
VCF_FILES=$(find "$VCF_DIR" -name "*.sorted.vcf.gz")

# For each SGD gene
for ((i=0; i<${#GENES[@]}; i+=4)); do
  SGD_GENE="${GENES[$i]}"
  W303_ID="${GENES[$i+1]}"
  GENE_ID="${GENES[$i+2]}"
  ERG_NAME="${GENES[$i+3]}"
  
  echo "Searching for $SGD_GENE ($ERG_NAME)..."
  
  # Output file
  OUTPUT_FILE="$OUTPUT_DIR/by_gene/${SGD_GENE}_variants.txt"
  echo "# Variants for $SGD_GENE ($ERG_NAME)" > "$OUTPUT_FILE"
  echo "# W303 ID: $W303_ID" >> "$OUTPUT_FILE"
  echo "# Gene ID: $GENE_ID" >> "$OUTPUT_FILE"
  echo "# Date: $(date)" >> "$OUTPUT_FILE"
  echo "#" >> "$OUTPUT_FILE"
  echo "# Sample Chromosome Position Ref Alt Info" >> "$OUTPUT_FILE"
  
  # Search for all identifiers in all VCF files
  for VCF_FILE in $VCF_FILES; do
    SAMPLE=$(basename "$VCF_FILE" .sorted.vcf.gz)
    echo "  Processing $SAMPLE..."
    
    # Search for all possible identifiers
    {
      zcat "$VCF_FILE" | grep -v "^#" | grep -e "$SGD_GENE" -e "$W303_ID" -e "$GENE_ID" -e "$ERG_NAME" || true
    } | while read -r LINE; do
      # Extract fields
      CHROM=$(echo "$LINE" | awk '{print $1}')
      POS=$(echo "$LINE" | awk '{print $2}')
      REF=$(echo "$LINE" | awk '{print $4}')
      ALT=$(echo "$LINE" | awk '{print $5}')
      INFO=$(echo "$LINE" | awk '{print $8}')
      
      # Add to output
      echo "$SAMPLE $CHROM $POS $REF $ALT $INFO" >> "$OUTPUT_FILE"
    done
  done
  
  # Count variants
  COUNT=$(wc -l < "$OUTPUT_FILE")
  COUNT=$((COUNT - 6))  # Subtract header lines
  if [ $COUNT -gt 0 ]; then
    echo "  Found $COUNT variants for $SGD_GENE"
  else
    echo "  No variants found for $SGD_GENE"
  fi
  echo ""
done

echo "=== Direct Gene Search Complete ==="
echo "Results saved to $OUTPUT_DIR"
""")
    
    os.chmod(fallback_script, 0o755)
    print(f"\nCreated fallback search script: {fallback_script}")
    
    print("\n=== Gene Mapping and Variant Extraction Complete ===")

if __name__ == "__main__":
    main()
