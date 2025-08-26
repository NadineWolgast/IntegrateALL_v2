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
USER_HOME="$(eval echo ~$USER)"
CONDA_BASE="$(conda info --base)"
WORK_DIR="$PIPELINE_DIR"

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
    echo "🔧 Auto-Konfiguration:"
    echo "   Pipeline Directory: $PIPELINE_DIR"
    echo "   Conda Base: $CONDA_BASE"
    echo "   User Home: $USER_HOME"
    echo ""
    echo "Diese Installation konfiguriert automatisch alle Pfade"
    echo "und verwendet Snakemake rules wie im Original-Projekt."
    echo ""
}

# Auto-configure paths for current environment
auto_configure_paths() {
    log_info "Konfiguriere Pfade automatisch..."
    
    # Create cluster config with auto-detected paths
    local cluster_config="config/cluster.yaml"
    if [ ! -f "$cluster_config" ]; then
        log_info "Erstelle Cluster-Konfiguration..."
        cat > "$cluster_config" << EOF
# Auto-generated SLURM Cluster Configuration
# Pipeline Directory: $PIPELINE_DIR
# Generated: $(date)

__default__:
  cores: 1
  memory: "4G"
  time: "01:00:00"
  partition: "base"

# Quality Control
fastqc:
  cores: 2
  memory: "8G"
  time: "30:00"

multiqc:
  cores: 1
  memory: "4G"
  time: "30:00"

# Alignment (Resource-intensiv)
star_index:
  cores: 16
  memory: "60G"
  time: "02:00:00"
  partition: "base"

star_align:
  cores: 16
  memory: "60G"
  time: "04:00:00"
  partition: "base"

# GATK Variant Calling
gatk_add_read_groups:
  cores: 4
  memory: "16G"
  time: "02:00:00"

gatk_mark_duplicates:
  cores: 4
  memory: "16G"
  time: "02:00:00"

gatk_split_n_cigar_reads:
  cores: 4
  memory: "16G"
  time: "03:00:00"

gatk_haplotype_caller:
  cores: 8
  memory: "32G"
  time: "06:00:00"

gatk_variant_filtration:
  cores: 2
  memory: "8G"
  time: "01:00:00"

# Fusion Detection
arriba:
  cores: 8
  memory: "32G"
  time: "02:00:00"

# Classification & Analysis
allcatchr:
  cores: 4
  memory: "16G"
  time: "01:00:00"

rnaseqcnv:
  cores: 4
  memory: "16G"
  time: "02:00:00"

# Reporting
generate_reports:
  cores: 2
  memory: "8G"
  time: "30:00"
EOF
        log_success "Cluster-Konfiguration erstellt: $cluster_config"
    fi
    
    # Create SLURM submit script with auto-detected paths
    local submit_script="submit_pipeline.sh"
    log_info "Erstelle SLURM Submit-Script..."
    
    cat > "$submit_script" << EOF
#!/bin/bash
#SBATCH --job-name=IntegrateALL_v2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --qos=long
#SBATCH --time=7-00:00:00
#SBATCH --output=integrateall_%j.out
#SBATCH --error=integrateall_%j.err

# Auto-configured paths
echo "Pipeline Directory: $PIPELINE_DIR"
echo "Conda Base: $CONDA_BASE"
echo "User: \$(whoami)"
echo "Host: \$(hostname)"
echo "Date: \$(date)"
echo ""

# Activate Conda (auto-detected)
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
elif [ -f "$USER_HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$USER_HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$USER_HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$USER_HOME/anaconda3/etc/profile.d/conda.sh"
else
    echo "Warning: Could not find conda.sh - trying direct activation"
    export PATH="$CONDA_BASE/bin:\$PATH"
fi

# Activate Pipeline Environment
conda activate $CONDA_ENV_NAME

# Change to Pipeline Directory
cd "$PIPELINE_DIR"

# Create logs directory
mkdir -p logs

# Run Pipeline with SLURM
snakemake \\
    --snakefile workflow/Snakefile \\
    --cluster "sbatch -t {cluster.time} -c {cluster.cores} --mem={cluster.memory} -p {cluster.partition} -o logs/slurm-%j.out -e logs/slurm-%j.err" \\
    --cluster-config config/cluster.yaml \\
    --use-conda \\
    --conda-prefix .snakemake/conda \\
    --cores 32 \\
    --jobs 42 \\
    --keep-going \\
    --rerun-triggers mtime \\
    --latency-wait 60 \\
    --restart-times 3 \\
    --printshellcmds
EOF

    chmod +x "$submit_script"
    log_success "SLURM Submit-Script erstellt: $submit_script"
    
    # Update activate_pipeline.sh with correct paths
    local activate_script="activate_pipeline.sh"
    log_info "Aktualisiere Pipeline-Aktivierungs-Script..."
    
    cat > "$activate_script" << EOF
#!/bin/bash
# Auto-generated Pipeline Activation Script
# Pipeline Directory: $PIPELINE_DIR
# Generated: $(date)

echo "🧬 IntegrateALL Pipeline - Environment Activation"
echo "================================================"
echo "Pipeline Directory: $PIPELINE_DIR"
echo "Conda Environment: $CONDA_ENV_NAME"
echo ""

# Auto-detect and source conda
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
elif [ -f "$USER_HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$USER_HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$USER_HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$USER_HOME/anaconda3/etc/profile.d/conda.sh"
else
    echo "⚠️  Warning: Could not auto-detect conda installation"
    echo "   Please manually activate conda before running this script"
fi

# Activate Pipeline Environment
if conda env list | grep -q "$CONDA_ENV_NAME"; then
    conda activate $CONDA_ENV_NAME
    echo "✅ Pipeline environment '$CONDA_ENV_NAME' activated"
    echo "📍 Working directory: \$(pwd)"
    echo ""
    echo "🚀 Ready to run pipeline!"
    echo "   Local: snakemake --cores 16 --use-conda"
    echo "   SLURM: sbatch submit_pipeline.sh"
else
    echo "❌ Pipeline environment '$CONDA_ENV_NAME' not found"
    echo "   Please run ./install_simple.sh first"
fi
EOF

    chmod +x "$activate_script"
    log_success "Pipeline-Aktivierungs-Script aktualisiert: $activate_script"
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
    
    # Auto-configure paths for current environment
    auto_configure_paths
    
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