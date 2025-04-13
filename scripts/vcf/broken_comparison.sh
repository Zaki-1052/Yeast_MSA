#!/bin/bash

# Create directory for comparison results
mkdir -p results/vcf/file_comparison
mkdir -p results/vcf/sample_lists

# Create sample list files
echo "WT-CTRL" > results/vcf/sample_lists/wt_ctrl.txt
echo "WT-37-55-1" > results/vcf/sample_lists/wt_treat.txt
echo "WT-37-55-2" >> results/vcf/sample_lists/wt_treat.txt
echo "WT-37-55-3" >> results/vcf/sample_lists/wt_treat.txt

echo "STC-CTRL" > results/vcf/sample_lists/stc_ctrl.txt
echo "STC-55-1" > results/vcf/sample_lists/stc_treat.txt
echo "STC-55-2" >> results/vcf/sample_lists/stc_treat.txt
echo "STC-55-3" >> results/vcf/sample_lists/stc_treat.txt

echo "CAS-CTRL" > results/vcf/sample_lists/cas_ctrl.txt
echo "CAS-55-1" > results/vcf/sample_lists/cas_treat.txt
echo "CAS-55-2" >> results/vcf/sample_lists/cas_treat.txt
echo "CAS-55-3" >> results/vcf/sample_lists/cas_treat.txt

echo "WT-CTRL" > results/vcf/sample_lists/wta_ctrl.txt 
echo "WTA-55-1" > results/vcf/sample_lists/wta_treat.txt
echo "WTA-55-2" >> results/vcf/sample_lists/wta_treat.txt
echo "WTA-55-3" >> results/vcf/sample_lists/wta_treat.txt

# Combine for group extraction
cat results/vcf/sample_lists/wt_ctrl.txt results/vcf/sample_lists/wt_treat.txt > results/vcf/sample_lists/wt_group.txt
cat results/vcf/sample_lists/stc_ctrl.txt results/vcf/sample_lists/stc_treat.txt > results/vcf/sample_lists/stc_group.txt
cat results/vcf/sample_lists/cas_ctrl.txt results/vcf/sample_lists/cas_treat.txt > results/vcf/sample_lists/cas_group.txt
cat results/vcf/sample_lists/wta_ctrl.txt results/vcf/sample_lists/wta_treat.txt > results/vcf/sample_lists/wta_group.txt

# Function to compare a treatment group with its control
compare_group_with_files() {
    GROUP=$1
    TREAT_FILE=$2
    CTRL_FILE=$3
    GROUP_FILE=$4
    
    echo "Comparing $GROUP treatments vs control using sample files"
    
    # Extract just this group using sample file
    bcftools view -S $GROUP_FILE results/vcf/merged/all_samples_fixed.vcf.gz -Oz -o results/vcf/file_comparison/${GROUP}_group.vcf.gz
    bcftools index results/vcf/file_comparison/${GROUP}_group.vcf.gz
    
    # Create treatment-only and control-only VCFs
    bcftools view -S $TREAT_FILE results/vcf/file_comparison/${GROUP}_group.vcf.gz -Oz -o results/vcf/file_comparison/${GROUP}_treat.vcf.gz
    bcftools index results/vcf/file_comparison/${GROUP}_treat.vcf.gz
    
    bcftools view -S $CTRL_FILE results/vcf/file_comparison/${GROUP}_group.vcf.gz -Oz -o results/vcf/file_comparison/${GROUP}_ctrl.vcf.gz
    bcftools index results/vcf/file_comparison/${GROUP}_ctrl.vcf.gz
    
    # Find variants in treatment but not control using isec
    ISEC_DIR="results/vcf/file_comparison/${GROUP}_isec"
    mkdir -p $ISEC_DIR
    
    bcftools isec results/vcf/file_comparison/${GROUP}_treat.vcf.gz \
              results/vcf/file_comparison/${GROUP}_ctrl.vcf.gz \
              -p $ISEC_DIR -n =1 -c any
    
    # Process treatment-specific variants
    if [ -f "${ISEC_DIR}/0000.vcf" ]; then
        bgzip -c "${ISEC_DIR}/0000.vcf" > "results/vcf/file_comparison/${GROUP}_specific.vcf.gz"
        bcftools index "results/vcf/file_comparison/${GROUP}_specific.vcf.gz"
        
        COUNT=$(bcftools view -H "results/vcf/file_comparison/${GROUP}_specific.vcf.gz" | wc -l)
        echo "$GROUP treatment-specific variants: $COUNT"
        
        # Find high-confidence variants (present in most treatment samples)
        THRESHOLD=$(($(wc -l < $TREAT_FILE) - 1))  # At least N-1 samples
        
        # Extract treatment samples from specific variants
        bcftools view -S $TREAT_FILE "results/vcf/file_comparison/${GROUP}_specific.vcf.gz" -Oz -o "results/vcf/file_comparison/${GROUP}_specific_treat.vcf.gz"
        bcftools index "results/vcf/file_comparison/${GROUP}_specific_treat.vcf.gz"
        
        # Filter for high-confidence variants
        bcftools view -c $THRESHOLD "results/vcf/file_comparison/${GROUP}_specific_treat.vcf.gz" -Oz -o "results/vcf/file_comparison/${GROUP}_highconf.vcf.gz"
        bcftools index "results/vcf/file_comparison/${GROUP}_highconf.vcf.gz"
        
        HIGH_COUNT=$(bcftools view -H "results/vcf/file_comparison/${GROUP}_highconf.vcf.gz" | wc -l)
        echo "$GROUP high-confidence treatment-specific variants: $HIGH_COUNT"
    else
        echo "No treatment-specific variants found for $GROUP"
    fi
}

# Run comparisons using file-based approach
compare_group_with_files "WT" results/vcf/sample_lists/wt_treat.txt results/vcf/sample_lists/wt_ctrl.txt results/vcf/sample_lists/wt_group.txt
compare_group_with_files "STC" results/vcf/sample_lists/stc_treat.txt results/vcf/sample_lists/stc_ctrl.txt results/vcf/sample_lists/stc_group.txt
compare_group_with_files "CAS" results/vcf/sample_lists/cas_treat.txt results/vcf/sample_lists/cas_ctrl.txt results/vcf/sample_lists/cas_group.txt
compare_group_with_files "WTA" results/vcf/sample_lists/wta_treat.txt results/vcf/sample_lists/wta_ctrl.txt results/vcf/sample_lists/wta_group.txt

echo "File-based comparison complete!"