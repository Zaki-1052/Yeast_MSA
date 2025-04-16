# Install required tools using Homebrew if not already installed
# Check if Homebrew is installed
if ! command -v brew &> /dev/null; then
    echo "Homebrew not found. Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
fi

# Install FastQC
if ! command -v fastqc &> /dev/null; then
    echo "Installing FastQC..."
    brew install fastqc
fi

# Install MultiQC (optional but recommended)
if ! command -v multiqc &> /dev/null; then
    echo "Installing MultiQC..."
    pip install multiqc
fi

# Check if BWA is installed (we'll need this later)
if ! command -v bwa &> /dev/null; then
    echo "Installing BWA..."
    brew install bwa
fi

# Check if samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "Installing Samtools..."
    brew install samtools
fi

# Check if bcftools is installed
if ! command -v bcftools &> /dev/null; then
    echo "Installing BCFtools..."
    brew install bcftools
fi

# Check if fastp is installed
if ! command -v fastp &> /dev/null; then
    echo "Installing fastp..."
    brew install fastp
fi

echo "All required tools are installed."
