# File: scripts/annotation/new/update_genbank_chromosome_names.py

import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def load_mapping():
    """Load the reverse mapping (YGAP→CM) from file."""
    mapping = {}
    with open('reference/chromosome_mapping_reverse.txt', 'r') as f:
        # Skip header
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                ygap_id, cm_id = parts
                mapping[ygap_id] = cm_id
    return mapping

def extract_scaffold_info(locus_line):
    """Extract scaffold information from the LOCUS line."""
    match = re.search(r'w303_scaffold_(\d+)', locus_line)
    if match:
        scaffold_num = match.group(1)
        return int(scaffold_num)
    return None

def extract_chromosome_from_note(features):
    """Extract chromosome information from the source feature /note qualifier."""
    if features and len(features) > 0:
        source_feature = features[0]
        if source_feature.type == "source":
            for qualifier in source_feature.qualifiers.get('note', []):
                # Look for Roman numeral or special chromosome name (MT, p2micron)
                if qualifier in ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 
                                'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'MT', 'p2micron',
                                'XII_Seq19', 'XII_Seq20', 'XII_Seq21']:
                    return qualifier
    return None

def update_genbank_files(input_dir, output_dir, mapping):
    """Update GenBank files to use CM007XXX.1 chromosome names."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of all genbank files
    genbank_files = [f for f in os.listdir(input_dir) if f.endswith(('.genbank', '.gb', '.gbk'))]
    print(f"Found {len(genbank_files)} GenBank files to process")
    
    success_count = 0
    error_count = 0
    
    for gb_file in genbank_files:
        input_path = os.path.join(input_dir, gb_file)
        
        try:
            # Parse the GenBank file
            record = SeqIO.read(input_path, "genbank")
            
            # Extract scaffold number from LOCUS line
            scaffold_num = extract_scaffold_info(record.name)
            
            # Extract chromosome info from /note qualifier
            chr_id = extract_chromosome_from_note(record.features)
            
            if chr_id and chr_id in mapping:
                cm_id = mapping[chr_id]
                
                # Update the record name (LOCUS line)
                record.name = cm_id
                
                # Generate output filename
                output_filename = f"{cm_id}.gbk"
                output_path = os.path.join(output_dir, output_filename)
                
                # Write the updated record
                SeqIO.write(record, output_path, "genbank")
                success_count += 1
                
                if success_count % 10 == 0:
                    print(f"Processed {success_count} files...")
            else:
                print(f"Warning: Could not determine chromosome ID for {gb_file}")
                error_count += 1
        
        except Exception as e:
            print(f"Error processing {gb_file}: {str(e)}")
            error_count += 1
    
    print(f"Complete! Successfully processed {success_count} files with {error_count} errors")
    return success_count, error_count

if __name__ == "__main__":
    # Define directories
    input_genbank_dir = "reference/w303_annotations/genbank"
    output_genbank_dir = "reference/w303_annotations/genbank_cm_names"
    
    # Load the mapping
    mapping = load_mapping()
    
    # Process the GenBank files
    success_count, error_count = update_genbank_files(input_genbank_dir, output_genbank_dir, mapping)
    
    if success_count > 0:
        print("\nNext steps:")
        print("1. Copy the modified GenBank files to SnpEff data directory:")
        print("   cp reference/w303_annotations/genbank_cm_names/*.gbk ~/snpEff/data/w303_cm/genes/")
        print("2. Create or update SnpEff config to include the w303_cm genome")
        print("3. Build the SnpEff database with:")
        print("   java -jar snpEff.jar build -genbank -v w303_cm")
    else:
        print("Failed to process any GenBank files successfully.")