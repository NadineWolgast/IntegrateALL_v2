# IntegrateALL Pipeline

An improved, high-performance pipeline for integrated analysis of B-cell acute lymphoblastic leukemia (B-ALL) RNA-seq data.

## Overview

IntegrateALL is a comprehensive Snakemake-based pipeline that processes RNA-seq data to perform:

- **Quality Control**: FastQC and MultiQC analysis
- **Alignment**: STAR-based RNA-seq alignment
- **Fusion Detection**: Arriba and FusionCatcher integration
- **Variant Calling**: GATK best practices for RNA-seq
- **Copy Number Analysis**: RNASeqCNV for chromosomal alterations
- **Classification**: ALLCatchR B-ALL subtype prediction
- **Karyotype Prediction**: Machine learning-based karyotype classification
- **Hotspot Analysis**: Targeted mutation analysis
- **Integration**: Comprehensive reporting and visualization

## Key Improvements

### Performance Optimizations
- Modular architecture with optimized resource allocation
- Parallel processing for independent samples
- Efficient temporary file management
- Checkpointing for resumable workflows

### Code Quality
- Unified Python codebase (replacing mixed bash/R/Python scripts)
- Comprehensive error handling and logging
- Standardized configuration management
- Clear separation of concerns

### Usability
- Interactive HTML reports
- Comprehensive logging and benchmarking
- Easy configuration through YAML files
- Docker/Singularity support

## Quick Start

### ⚡ Ein-Befehl Installation

```bash
# Komplette automatische Installation
./install.sh

# Oder Python-basierte Installation
python setup.py
```

**Das war's!** Die Installation konfiguriert automatisch:
- ✅ Conda/Mamba Environment Management
- ✅ Alle Pipeline Dependencies  
- ✅ Referenz-Genome und Datenbanken
- ✅ R Package Installation
- ✅ Pipeline Konfiguration
- ✅ Installations-Test

### 🚀 Pipeline aktivieren und ausführen

```bash
# Pipeline aktivieren
source activate_pipeline.sh

# Sample-Sheet bearbeiten
nano config/samples.tsv

# Pipeline testen
snakemake --dry-run --cores 1

# Pipeline ausführen
snakemake --cores 16 --use-conda
```

### 📝 Sample-Sheet Format

Bearbeite `config/samples.tsv`:
```tsv
sample_id	fastq1	fastq2	condition	batch
sample_001	/path/to/sample_001_R1.fastq.gz	/path/to/sample_001_R2.fastq.gz	B-ALL	batch1
sample_002	/path/to/sample_002_R1.fastq.gz	/path/to/sample_002_R2.fastq.gz	B-ALL	batch1
```

## Directory Structure

```
IntegrateALL_pipeline/
├── config/                 # Configuration files
│   ├── config.yaml         # Main pipeline configuration
│   ├── samples.tsv         # Sample sheet
│   └── resources.yaml      # Resource specifications
├── workflow/
│   ├── rules/              # Modular Snakemake rules
│   ├── scripts/            # Python analysis scripts
│   ├── envs/              # Conda environments
│   └── Snakefile          # Main workflow file
├── resources/              # Reference data
│   ├── genomes/           # Genome references
│   └── databases/         # Tool databases
├── results/               # Output results
│   ├── qc/               # Quality control
│   ├── alignment/        # STAR alignments
│   ├── fusions/          # Fusion detection
│   ├── variants/         # Variant calling
│   ├── classification/   # B-ALL classification
│   └── reports/          # Final reports
├── logs/                  # Execution logs
└── benchmarks/           # Performance benchmarks
```

## Configuration

The pipeline uses YAML-based configuration:

- `config/config.yaml`: Main pipeline parameters
- `config/samples.tsv`: Sample metadata and file paths
- `config/resources.yaml`: Computational resources per rule

## Requirements

- Snakemake ≥ 7.0
- Conda/Mamba
- Singularity (optional)
- Python ≥ 3.9
- R ≥ 4.0 (for specific tools that require it)

## Citation

If you use IntegrateALL in your research, please cite:
[Citation will be added upon publication]

## License

[License information]

## Support

For questions and support, please open an issue on GitHub or contact the developers.