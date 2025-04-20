# File: scripts/annotation/new/fix_genbank_coordinates.py

import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def fix_genbank_files():
    """Fix gene coordinates in GenBank files."""
    input_dir = os.path.expanduser("~/snpEff/data/w303_cm/genes")
    output_dir = os.path.expanduser("~/snpEff/data/w303_cm/genes_fixed")
    
    os.makedirs(output_dir, exist_ok=True)
    print(f"Created output directory: {output_dir}")
    
    # Get list of all GenBank files
    gbk_files = [f for f in os.listdir(input_dir) if f.endswith('.gbk')]
    print(f"Found {len(gbk_files)} GenBank files to process")
    
    total_features = 0
    fixed_features = 0
    
    for gb_file in gbk_files:
        input_path = os.path.join(input_dir, gb_file)
        output_path = os.path.join(output_dir, gb_file)
        
        try:
            # Read GenBank file
            record = SeqIO.read(input_path, "genbank")
            seq_length = len(record.seq)
            
            # Track if changes were made
            changes_made = False
            
            # Process each feature
            new_features = []
            for feature in record.features:
                total_features += 1
                
                # Skip source feature (first feature that defines entire sequence)
                if feature.type == "source" and feature == record.features[0]:
                    new_features.append(feature)
                    continue
                
                # Check if feature extends beyond sequence boundaries
                start = feature.location.start
                end = feature.location.end
                
                if start < 0 or end > seq_length:
                    # Fix the coordinates to be within sequence boundaries
                    if start < 0:
                        start = 0
                    if end > seq_length:
                        end = seq_length
                    
                    # Create new feature with fixed coordinates
                    new_location = FeatureLocation(start, end, feature.location.strand)
                    new_feature = SeqFeature(
                        location=new_location,
                        type=feature.type,
                        id=feature.id,
                        qualifiers=feature.qualifiers
                    )
                    new_features.append(new_feature)
                    fixed_features += 1
                    changes_made = True
                else:
                    # Keep the original feature if coordinates are fine
                    new_features.append(feature)
            
            # Update the record features if changes were made
            if changes_made:
                record.features = new_features
                
            # Add circular chromosome annotation if needed
            circular_annotation = False
            for feature in record.features:
                if feature.type == "source":
                    if 'note' in feature.qualifiers:
                        for note in feature.qualifiers['note']:
                            if 'circular' in note.lower():
                                circular_annotation = True
                                break
                    
                    if not circular_annotation:
                        # Add circular note to source feature
                        if 'note' in feature.qualifiers:
                            feature.qualifiers['note'].append('circular chromosome')
                        else:
                            feature.qualifiers['note'] = ['circular chromosome']
                    break
            
            # Write the modified record
            SeqIO.write(record, output_path, "genbank")
            print(f"Processed {gb_file} (length: {seq_length})")
        
        except Exception as e:
            print(f"Error processing {gb_file}: {str(e)}")
    
    print(f"Fixed {fixed_features} features out of {total_features} total features")
    
    # Create a combined file
    combined_path = os.path.expanduser("~/snpEff/data/w303_cm/genes.gbk")
    with open(combined_path, 'w') as outfile:
        for gb_file in sorted(gbk_files):
            file_path = os.path.join(output_dir, gb_file)
            if os.path.exists(file_path):
                with open(file_path, 'r') as infile:
                    outfile.write(infile.read())
                    outfile.write("\n\n")
    
    print(f"Created combined GenBank file: {combined_path}")
    
    # Create a build script
    build_script = "build_snpeff_w303cm_fixed.sh"
    with open(build_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("cd ~/snpEff\n")
        f.write("java -jar snpEff.jar build -genbank -v w303_cm\n")
    
    # Make the script executable
    os.chmod(build_script, 0o755)
    print(f"Created build script: {build_script}")
    
    print("\nNext steps:")
    print(f"1. Run the build script: ./{build_script}")

if __name__ == "__main__":
    fix_genbank_files()