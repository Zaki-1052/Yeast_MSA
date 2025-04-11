#!/bin/bash

# Create directories for processed data
mkdir -p data/processed_data
mkdir -p results/fastp

# Process each sample with fastp
for R1 in data/raw_data/*_R1_001.fastq.gz; do
    # Get sample name and R2 file
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    echo "Processing sample: $SAMPLE"
    
    # Run fastp with quality trimming and filtering
    fastp \
        --in1 "$R1" \
        --in2 "$R2" \
        --out1 "data/processed_data/${SAMPLE}_R1_cleaned.fastq.gz" \
        --out2 "data/processed_data/${SAMPLE}_R2_cleaned.fastq.gz" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 50 \
        --cut_right \
        --cut_right_window_size 4 \
        --cut_right_mean_quality 20 \
        --html "results/fastp/${SAMPLE}_fastp.html" \
        --json "results/fastp/${SAMPLE}_fastp.json" \
        --report_title "${SAMPLE}" \
        --thread 4
done

echo "Preprocessing complete. Cleaned data is in data/processed_data directory."