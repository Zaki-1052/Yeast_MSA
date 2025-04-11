#!/bin/bash

# Create directory for FastQC results if it doesn't exist
mkdir -p results/fastqc

# Run FastQC on all FASTQ files
echo "Running FastQC on all samples..."
fastqc data/raw_data/*.fastq.gz -o results/fastqc -t 4

# Optional: Run MultiQC to aggregate results (if installed)
if command -v multiqc &> /dev/null; then
    echo "Generating MultiQC report..."
    multiqc results/fastqc -o results/fastqc
else
    echo "MultiQC not found. Consider installing it for aggregated reports."
fi

echo "FastQC analysis complete. Results are in results/fastqc directory."