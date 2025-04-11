#!/bin/bash

# Create directories for alignment results (should already exist)
mkdir -p results/bam
mkdir -p results/stats

# Path to reference genome
REF="reference/yeast_w303.fasta"

# Process each sample
for R1 in data/processed_data/*_R1_cleaned.fastq.gz; do
    # Get sample name and R2 file
    SAMPLE=$(basename "$R1" _R1_cleaned.fastq.gz)
    R2="${R1/_R1_cleaned.fastq.gz/_R2_cleaned.fastq.gz}"
    
    # Skip if final BAM already exists and is of proper size (>100MB)
    if [ -f "results/bam/${SAMPLE}.final.bam" ]; then
        FILE_SIZE=$(stat -f %z "results/bam/${SAMPLE}.final.bam")
        if [ $FILE_SIZE -gt 100000000 ]; then
            echo "Skipping ${SAMPLE} - already processed"
            continue
        fi
    fi
    
    echo "Aligning sample: $SAMPLE ($(date))"
    
    # 1. Align with BWA-MEM and pipe to samtools to convert to BAM
    bwa mem -t 4 $REF "$R1" "$R2" | \
    samtools view -b -o "results/bam/${SAMPLE}.raw.bam"
    
    # 2. Sort BAM file by name (required for fixmate)
    samtools sort -n -o "results/bam/${SAMPLE}.namesorted.bam" "results/bam/${SAMPLE}.raw.bam"
    
    # 3. Fix mate information and add MS tags
    samtools fixmate -m "results/bam/${SAMPLE}.namesorted.bam" "results/bam/${SAMPLE}.fixmate.bam"
    
    # 4. Sort by coordinate (required for markdup)
    samtools sort -o "results/bam/${SAMPLE}.sorted.bam" "results/bam/${SAMPLE}.fixmate.bam"
    
    # 5. Mark duplicates
    samtools markdup -r "results/bam/${SAMPLE}.sorted.bam" "results/bam/${SAMPLE}.markdup.bam"
    
    # 6. Add read group information
    samtools addreplacerg -r "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:lib1\tPL:ILLUMINA" \
    -o "results/bam/${SAMPLE}.final.bam" "results/bam/${SAMPLE}.markdup.bam"
    
    # 7. Index final BAM file
    samtools index "results/bam/${SAMPLE}.final.bam"
    
    # 8. Generate alignment statistics
    samtools flagstat "results/bam/${SAMPLE}.final.bam" > "results/stats/${SAMPLE}.flagstat.txt"
    samtools coverage "results/bam/${SAMPLE}.final.bam" > "results/stats/${SAMPLE}.coverage.txt"
    
    # 9. Clean up intermediate files (uncomment if you want to save disk space)
    rm "results/bam/${SAMPLE}.raw.bam"
    rm "results/bam/${SAMPLE}.namesorted.bam"
    rm "results/bam/${SAMPLE}.fixmate.bam"
    rm "results/bam/${SAMPLE}.sorted.bam"
    rm "results/bam/${SAMPLE}.markdup.bam"
    
    echo "Completed alignment for $SAMPLE ($(date))"
done

# Create a summary of alignment statistics
echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapped_Percent\tDuplicates\tDuplicate_Percent" > results/stats/alignment_summary.txt
for STAT in results/stats/*.flagstat.txt; do
    SAMPLE=$(basename "$STAT" .flagstat.txt)
    TOTAL=$(grep "in total" "$STAT" | awk '{print $1}')
    MAPPED=$(grep "mapped (" "$STAT" | awk '{print $1}')
    MAPPED_PCT=$(grep "mapped (" "$STAT" | awk '{print $5}' | tr -d '(%)')
    DUPS=$(grep "duplicates" "$STAT" | head -1 | awk '{print $1}')
    
    if [ -n "$TOTAL" ] && [ "$TOTAL" -gt 0 ]; then
        DUP_PCT=$(echo "scale=2; 100*$DUPS/$TOTAL" | bc)
    else
        DUP_PCT="N/A"
    fi
    
    echo -e "$SAMPLE\t$TOTAL\t$MAPPED\t$MAPPED_PCT\t$DUPS\t$DUP_PCT" >> results/stats/alignment_summary.txt
done

echo "Alignment complete. Final BAM files are in results/bam directory."
echo "Alignment statistics summary saved to results/stats/alignment_summary.txt"