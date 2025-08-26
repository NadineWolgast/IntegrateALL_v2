#!/bin/bash
# IntegrateALL Pipeline Aktivierung

# Aktiviere Conda Environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate integrateall

# Wechsle in Pipeline Directory
cd "/media/nadine/InternalMaybe/IntegrateALL_pipeline"

echo "=================================================="
echo "  IntegrateALL Pipeline ist bereit!"
echo "=================================================="
echo ""
echo "Verfügbare Befehle:"
echo "  snakemake --dry-run --cores 1        # Test der Konfiguration"
echo "  snakemake --cores 16 --use-conda     # Pipeline lokal ausführen"
echo "  snakemake --help                     # Hilfe anzeigen"
echo ""
echo "Konfiguration bearbeiten:"
echo "  nano config/config.yaml              # Haupt-Konfiguration"
echo "  nano config/samples.tsv              # Sample-Liste"
echo ""
