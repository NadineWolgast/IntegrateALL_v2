#!/bin/bash

# IntegrateALL Pipeline - Simple Snakemake-based Installation
# ==========================================================
# Uses Snakemake rules for installation like the original project

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV_NAME="integrateall"

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
    exit 1
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_header() {
    echo "=================================================="
    echo "  IntegrateALL Pipeline - Installation"  
    echo "=================================================="
    echo ""
    echo "Diese Installation verwendet Snakemake rules wie im"
    echo "ursprünglichen Blast-o-Matic-Fusioninator Projekt."
    echo ""
}

# Check if conda is available
check_conda() {
    if ! command -v conda &> /dev/null; then
        log_error "Conda ist nicht installiert. Bitte installiere Miniconda/Anaconda erst."
    fi
    log_info "Conda gefunden: $(conda --version)"
}

# Setup conda environment if needed
setup_environment() {
    log_info "Prüfe Conda Environment '$CONDA_ENV_NAME'..."
    
    # Source conda
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    fi
    
    # Create environment if it doesn't exist
    if ! conda env list | grep -q "^$CONDA_ENV_NAME\s"; then
        log_info "Erstelle Conda Environment '$CONDA_ENV_NAME'..."
        conda env create -f workflow/envs/python.yaml -n $CONDA_ENV_NAME
    else
        log_info "Environment '$CONDA_ENV_NAME' bereits vorhanden"
    fi
    
    # Activate environment
    conda activate $CONDA_ENV_NAME
    
    # Install/update Snakemake in the environment
    if ! command -v snakemake &> /dev/null; then
        log_info "Installiere Snakemake..."
        conda install -c bioconda -c conda-forge snakemake=7.32.4 -y
    fi
    
    log_success "Environment setup complete"
}

main() {
    print_header
    
    # Parse command line arguments for help
    if [[ "${1:-}" == "--help" ]] || [[ "${1:-}" == "-h" ]]; then
        echo "Usage: $0 [OPTIONS]"
        echo ""
        echo "Options:"
        echo "  --help, -h    Diese Hilfe anzeigen"
        echo ""
        echo "Die Installation läuft vollautomatisch über Snakemake rules:"
        echo "  1. Conda Environment Setup"
        echo "  2. R Packages Installation" 
        echo "  3. Reference Data Download"
        echo "  4. Pipeline Configuration"
        echo "  5. Installation Test"
        echo ""
        exit 0
    fi
    
    cd "$PIPELINE_DIR"
    
    # Check prerequisites
    check_conda
    setup_environment
    
    log_info "Starte Snakemake-basierte Installation..."
    echo ""
    
    # Run installation through Snakemake
    if snakemake --use-conda --cores 4 install_all; then
        log_success "Installation erfolgreich abgeschlossen!"
        echo ""
        echo "🎉 IntegrateALL Pipeline ist bereit!"
        echo ""
        echo "Nächste Schritte:"
        echo "  1. Pipeline aktivieren:    source activate_pipeline.sh"
        echo "  2. Samples konfigurieren:  nano config/samples.tsv"
        echo "  3. Pipeline testen:        snakemake --dry-run --cores 1"
        echo "  4. Pipeline ausführen:     snakemake --cores 16 --use-conda"
        echo ""
    else
        log_error "Installation fehlgeschlagen. Prüfe die Logs für Details."
    fi
}

# Run main function
main "$@"