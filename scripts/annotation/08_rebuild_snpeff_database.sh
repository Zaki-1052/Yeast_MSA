#!/bin/bash

# File: scripts/annotation/08_rebuild_snpeff_database.sh
# Purpose: Rebuild the SnpEff database to fix chromosome naming issues

echo "=== Rebuilding SnpEff Database ==="
echo "Date: $(date)"
echo ""

# Define directories
SNPEFF_DIR="/Users/zakiralibhai/snpEff"
GENBANK_DIR="annotation/reference/w303_scaffolds"
BACKUP_DIR="annotation/backup"
W303_DATA_DIR="${SNPEFF_DIR}/data/w303"

# Create backup directory
mkdir -p "$BACKUP_DIR"

# Check for required files
if [ ! -d "$GENBANK_DIR" ]; then
    echo "ERROR: GenBank directory not found: $GENBANK_DIR"
    exit 1
fi

if [ ! -d "$SNPEFF_DIR" ]; then
    echo "ERROR: SnpEff directory not found: $SNPEFF_DIR"
    exit 1
fi

# Check for existing database and backup if it exists
if [ -d "$W303_DATA_DIR" ]; then
    echo "Found existing w303 database, backing up..."
    tar -czf "${BACKUP_DIR}/w303_database_backup_$(date +%Y%m%d%H%M%S).tar.gz" -C "$SNPEFF_DIR/data" "w303"
    rm -rf "$W303_DATA_DIR"
    echo "Backup created and old database removed"
fi

# Create the database directory structure
echo "Creating new database directory structure..."
mkdir -p "$W303_DATA_DIR"
mkdir -p "${W303_DATA_DIR}/genomes"

# Step 1: Copy all GenBank files to a single directory for easier processing
echo "Preparing GenBank files..."
# Convert GenBank files to single file (if needed)
cat "$GENBANK_DIR"/*.genbank > "${W303_DATA_DIR}/genes.gbk"
echo "Created combined genes.gbk file"

# Step 2: Extract chromosome names from VCF files
echo "Extracting chromosome names from VCF files..."
VCF_FILE=$(find "annotation/vcf_ready" -name "*.sorted.vcf.gz" | head -n 1)
if [ -z "$VCF_FILE" ]; then
    echo "ERROR: No VCF files found to extract chromosome names"
    exit 1
fi

echo "Using VCF file: $VCF_FILE"
bcftools view -h "$VCF_FILE" | grep '##contig=' | head -5
echo "..."

# Step 3: Set up snpEff.config
echo "Setting up snpEff configuration..."
CONFIG_FILE="${SNPEFF_DIR}/snpEff.config"

# Backup the config file
cp "$CONFIG_FILE" "${BACKUP_DIR}/snpEff.config.backup"

# Add or update the w303 entry in the config file
if grep -q "^w303." "$CONFIG_FILE"; then
    echo "Updating existing w303 entry in config file..."
    # Use sed to update the entry
    sed -i'.bak' '/^w303\./d' "$CONFIG_FILE"
fi

# Add new entry at the end of the file
echo "" >> "$CONFIG_FILE"
echo "# W303 genome with JRIU chromosomes" >> "$CONFIG_FILE"
echo "w303.genome : Saccharomyces cerevisiae W303" >> "$CONFIG_FILE"
echo "w303.chromosomes : JRIU01000001.1, JRIU01000002.1, JRIU01000003.1" >> "$CONFIG_FILE"
echo "w303.JRIU01000001.1.codonTable : Standard" >> "$CONFIG_FILE"

echo "Updated config file with w303 genome entry"

# Step 4: Build the database
echo "Building SnpEff database..."
cd "$SNPEFF_DIR"
java -jar snpEff.jar build -genbank -v w303
cd - > /dev/null

# Step 5: Verify database was built correctly
echo "Verifying database..."
if [ -f "${W303_DATA_DIR}/snpEffectPredictor.bin" ]; then
    echo "✅ Database built successfully: ${W303_DATA_DIR}/snpEffectPredictor.bin exists"
else
    echo "❌ Database build failed: ${W303_DATA_DIR}/snpEffectPredictor.bin not found"
    exit 1
fi

# Print a list of files in the database directory to verify
echo "Files in w303 database directory:"
ls -la "$W303_DATA_DIR"

echo ""
echo "Next steps:"
echo "1. Run a test annotation using the new database"
echo "2. If successful, proceed with annotating all VCF files"
echo "3. If test fails, check the SnpEff build logs for errors"
echo ""

echo "=== SnpEff Database Rebuild Complete ==="
