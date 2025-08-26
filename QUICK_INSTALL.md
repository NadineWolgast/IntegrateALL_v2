# IntegrateALL Pipeline - Ein-Befehl Installation

## 🚀 Schnellinstallation

```bash
git clone <repository-url> IntegrateALL_pipeline
cd IntegrateALL_pipeline
./install_simple.sh
```

**Das wars!** 🎉

## 🔧 Snakemake-basierte Installation

Die neue Installation verwendet Snakemake rules wie im ursprünglichen Blast-o-Matic-Fusioninator Projekt:

```bash
./install_simple.sh           # Vollautomatische Installation über Snakemake
./install_simple.sh --help    # Hilfe anzeigen
```

**Vorteile:**
- ✅ Keine interaktiven Bestätigungen erforderlich
- ✅ Automatische Parallelisierung 
- ✅ Intelligente Dependency-Behandlung
- ✅ Robuste Fehlerbehandlung
- ✅ Wie im ursprünglichen Projekt bewährt

## 🧠 Alternative: Intelligente Installation

Falls gewünscht, steht auch die smarte bash-basierte Installation zur Verfügung:

```bash
./install.sh --yes            # Installation ohne Bestätigungen
./install.sh --force --yes    # Komplette Neuinstallation ohne Bestätigungen
```

## ⏱️ Erwartete Dauer: 30-60 Minuten

Die Installation läuft vollautomatisch und installiert:
- ✅ Conda Environment mit allen Tools
- ✅ Arriba (Fusion Detection) 
- ✅ STAR, GATK4, FastQC, MultiQC
- ✅ R-Pakete (ALLCatchR, RNAseqCNV)
- ✅ Referenz-Genome (GRCh38)
- ✅ Python Analyse-Skripte

## 🏁 Nach der Installation:

```bash
source activate_pipeline.sh
nano config/samples.tsv              # Deine Samples eintragen
snakemake --dry-run --cores 1        # Test
snakemake --cores 16 --use-conda     # Pipeline starten
```

**Ready for B-ALL RNA-seq Analysis!** 🧬