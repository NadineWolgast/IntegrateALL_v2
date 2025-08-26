#!/bin/bash

# IntegrateALL Pipeline - Path Configuration Script
# ================================================
# Allows manual configuration of paths for different environments

set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_header() {
    echo "=================================================="
    echo "  IntegrateALL Pipeline - Path Configuration"  
    echo "=================================================="
    echo ""
    echo "🔧 Konfiguriere Pfade für dein System"
    echo ""
}

configure_paths() {
    local base_path="${1:-}"
    local conda_path="${2:-}"
    
    if [ -z "$base_path" ]; then
        echo "Aktueller Pipeline-Pfad: $PIPELINE_DIR"
        echo ""
        read -p "Basis-Arbeitsverzeichnis (Enter für aktuelles): " base_path
        if [ -z "$base_path" ]; then
            base_path="$PIPELINE_DIR"
        fi
    fi
    
    if [ -z "$conda_path" ]; then
        # Try to auto-detect conda
        local detected_conda=""
        if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
            detected_conda="$HOME/miniconda3"
        elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
            detected_conda="$HOME/anaconda3"
        elif command -v conda &> /dev/null; then
            detected_conda="$(conda info --base)"
        fi
        
        echo "Erkannter Conda-Pfad: $detected_conda"
        read -p "Conda-Installation Pfad (Enter für erkannten): " conda_path
        if [ -z "$conda_path" ]; then
            conda_path="$detected_conda"
        fi
    fi
    
    log_info "Konfiguriere Pfade:"
    echo "   Pipeline: $base_path"
    echo "   Conda: $conda_path"
    echo ""
    
    # Update submit_pipeline.sh
    if [ -f "submit_pipeline.sh" ]; then
        log_info "Aktualisiere submit_pipeline.sh..."
        sed -i "s|^cd \".*IntegrateALL.*\"|cd \"$base_path\"|g" submit_pipeline.sh
        sed -i "s|source \".*conda.sh\"|source \"$conda_path/etc/profile.d/conda.sh\"|g" submit_pipeline.sh
        log_success "submit_pipeline.sh aktualisiert"
    fi
    
    # Update activate_pipeline.sh  
    if [ -f "activate_pipeline.sh" ]; then
        log_info "Aktualisiere activate_pipeline.sh..."
        sed -i "s|source \".*conda.sh\"|source \"$conda_path/etc/profile.d/conda.sh\"|g" activate_pipeline.sh
        log_success "activate_pipeline.sh aktualisiert"
    fi
    
    # Update config.yaml with correct paths
    if [ -f "config/config.yaml" ]; then
        log_info "Aktualisiere config.yaml Pfade..."
        
        # Update reference genome paths
        sed -i "s|^reference_genome:.*|reference_genome: \"$base_path/resources/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa\"|g" config/config.yaml
        sed -i "s|^reference_gtf:.*|reference_gtf: \"$base_path/resources/genomes/Homo_sapiens.GRCh38.109.gtf\"|g" config/config.yaml
        
        # Update database paths
        sed -i "s|arriba_blacklist:.*|arriba_blacklist: \"$base_path/resources/databases/arriba/blacklist_hg38_GRCh38_v2.4.0.tsv.gz\"|g" config/config.yaml
        sed -i "s|arriba_known_fusions:.*|arriba_known_fusions: \"$base_path/resources/databases/arriba/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz\"|g" config/config.yaml
        sed -i "s|arriba_protein_domains:.*|arriba_protein_domains: \"$base_path/resources/databases/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3\"|g" config/config.yaml
        sed -i "s|arriba_cytobands:.*|arriba_cytobands: \"$base_path/resources/databases/arriba/cytobands_hg38_GRCh38_v2.4.0.tsv\"|g" config/config.yaml
        sed -i "s|driver_gene_list:.*|driver_gene_list: \"$base_path/resources/databases/driver_genes.txt\"|g" config/config.yaml
        
        log_success "config.yaml Pfade aktualisiert"
    fi
    
    # Create/update sample template
    log_info "Erstelle Sample-Template..."
    cat > "config/samples_template.tsv" << EOF
# Sample Configuration Template
# ============================
# Ersetze die Pfade mit deinen echten Daten-Pfaden
#
# Beispiel für lokale Daten:
# sample_id	fastq1	fastq2	condition	batch
# sample1	$base_path/data/sample1_R1.fastq.gz	$base_path/data/sample1_R2.fastq.gz	B-ALL	batch1
#
# Beispiel für Cluster-Daten:
# sample_id	fastq1	fastq2	condition	batch  
# sample1	/work_beegfs/user/data/sample1_R1.fastq.gz	/work_beegfs/user/data/sample1_R2.fastq.gz	B-ALL	batch1

sample_id	fastq1	fastq2	condition	batch
test_sample	$base_path/data/test_R1.fastq.gz	$base_path/data/test_R2.fastq.gz	B-ALL	batch1
EOF
    
    log_success "Sample-Template erstellt: config/samples_template.tsv"
    
    echo ""
    log_success "Pfad-Konfiguration abgeschlossen!"
    echo ""
    echo "📋 Nächste Schritte:"
    echo "1. Kopiere und bearbeite das Sample-Template:"
    echo "   cp config/samples_template.tsv config/samples.tsv"
    echo "   nano config/samples.tsv"
    echo ""
    echo "2. Für SLURM-Cluster:"
    echo "   sbatch submit_pipeline.sh"
    echo ""
    echo "3. Für lokale Ausführung:"
    echo "   source activate_pipeline.sh"
    echo "   snakemake --cores 16 --use-conda"
}

main() {
    print_header
    
    # Parse arguments
    local base_path="${1:-}"
    local conda_path="${2:-}"
    
    if [[ "${1:-}" == "--help" ]] || [[ "${1:-}" == "-h" ]]; then
        echo "Usage: $0 [BASE_PATH] [CONDA_PATH]"
        echo ""
        echo "Arguments:"
        echo "  BASE_PATH   Basis-Arbeitsverzeichnis (optional)"
        echo "  CONDA_PATH  Pfad zur Conda-Installation (optional)"
        echo ""
        echo "Examples:"
        echo "  $0                                    # Interaktive Konfiguration"
        echo "  $0 /work_beegfs/user/pipeline         # Mit Basis-Pfad"  
        echo "  $0 /work_beegfs/user/pipeline /opt/conda  # Mit beiden Pfaden"
        echo ""
        exit 0
    fi
    
    cd "$PIPELINE_DIR"
    configure_paths "$base_path" "$conda_path"
}

main "$@"