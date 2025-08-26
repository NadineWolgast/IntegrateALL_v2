#!/bin/bash
set -euo pipefail

RESOURCES_DIR="$1"
LOG_FILE="$2"

log_info() {
    echo "[$(date)] INFO: $1" | tee -a "$LOG_FILE"
}

cd "$RESOURCES_DIR"

# Download human reference genome (if not exists)
if [ ! -f "genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
    log_info "Downloading human reference genome..."
    mkdir -p genomes
    cd genomes
    
    # Download from Ensembl
    wget -c ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    wget -c ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
    
    # Uncompress
    gunzip -f *.gz
    
    cd ..
    log_info "Reference genome downloaded"
else
    log_info "Reference genome already exists"
fi

# Download Arriba databases
if [ ! -d "databases/arriba" ]; then
    log_info "Downloading Arriba databases..."
    mkdir -p databases/arriba
    cd databases/arriba
    
    wget -c https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz
    tar -xzf arriba_v2.4.0.tar.gz
    mv arriba_v2.4.0/database/* .
    rm -rf arriba_v2.4.0*
    
    cd ../..
    log_info "Arriba databases downloaded"
else
    log_info "Arriba databases already exist"
fi

# Create essential database files
mkdir -p databases
cd databases

# Create driver genes list
cat > driver_genes.txt << 'EOF'
ABL1
ABL2
AFF1
BCL2
BCL6
BCR
CREBBP
CRLF2
DUX4
EBF1
EP300
ETV6
IKZF1
KMT2A
MLL
MYC
NRAS
KRAS
PAX5
PBX1
RUNX1
TCF3
TP53
ZNF384
EOF

# Create oncogenes list
cat > oncogenes.txt << 'EOF'
MYC
BCL2
BCL6
ABL1
EOF

# Create tumor suppressors list
cat > tumor_suppressors.txt << 'EOF'
TP53
RB1
CDKN2A
EOF

log_info "Essential database files created"

