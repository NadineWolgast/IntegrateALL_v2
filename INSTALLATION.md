# IntegrateALL Pipeline Installation Guide

## Overview

IntegrateALL is a comprehensive, high-performance pipeline for integrated analysis of B-cell acute lymphoblastic leukemia (B-ALL) RNA-seq data. This installation guide will help you set up the pipeline on your system.

## System Requirements

### Minimum Hardware Requirements
- **CPU**: 16 cores recommended (minimum 8 cores)
- **RAM**: 64 GB recommended (minimum 32 GB)
- **Storage**: 1 TB free space (for reference data and results)
- **Operating System**: Linux (Ubuntu 18.04+, CentOS 7+, or similar)

### Software Dependencies

#### Core Requirements
- **Conda/Mamba**: For environment management
- **Snakemake**: ≥ 7.0 (workflow management)
- **Python**: ≥ 3.9
- **R**: ≥ 4.0 (for specific tools like ALLCatchR)

#### Optional (but recommended)
- **Singularity**: For containerized tool execution
- **SLURM/PBS**: For cluster computing

## Installation Steps

### 1. Clone the Repository

```bash
git clone <repository-url>
cd IntegrateALL_pipeline
```

### 2. Install Conda/Mamba

If you don't have Conda installed:

```bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Install Mamba (faster than conda)
conda install -c conda-forge mamba
```

### 3. Create Main Environment

```bash
# Create main pipeline environment
mamba env create -f workflow/envs/python.yaml
conda activate integrateall_python

# Install Snakemake
mamba install -c bioconda -c conda-forge snakemake=7.32.4
```

### 4. Download Reference Data

Create a script to download necessary reference files:

```bash
#!/bin/bash
# download_references.sh

# Create reference directory
mkdir -p resources/genomes resources/databases

# Download human reference genome (GRCh38)
cd resources/genomes
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

# Uncompress
gunzip *.gz

cd ../../

# Download Arriba databases
mkdir -p resources/databases/arriba
cd resources/databases/arriba
wget https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz
tar -xzf arriba_v2.4.0.tar.gz
mv arriba_v2.4.0/database/* .

cd ../../../

# Download FusionCatcher database
mkdir -p resources/databases/fusioncatcher
cd resources/databases/fusioncatcher
wget https://github.com/ndaniel/fusioncatcher/archive/refs/heads/master.zip
unzip master.zip
cd fusioncatcher-master/data
# Note: This will take several hours
./download-human-db.sh

echo "Reference download completed!"
```

Make the script executable and run it:

```bash
chmod +x download_references.sh
./download_references.sh
```

### 5. Install R Packages

Some tools require R packages. Install them manually:

```bash
# Activate R environment
conda activate integrateall_classification

# Install ALLCatchR
Rscript -e "
library(devtools)
install_github('ThomasBeder/ALLCatchR_bcrabl1')
"

# Install RNASeqCNV
Rscript -e "
library(devtools)
install_github('honzee/RNAseqCNV')
"
```

### 6. Configure the Pipeline

#### Update Configuration Files

1. **Edit `config/config.yaml`**:
```yaml
# Update paths to your reference files
reference_genome: "/path/to/your/resources/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
reference_gtf: "/path/to/your/resources/genomes/Homo_sapiens.GRCh38.109.gtf"

# Update database paths
databases:
  arriba_blacklist: "resources/databases/arriba/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"
  fusioncatcher_data: "resources/databases/fusioncatcher/fusioncatcher-master/data/human_v102"
  # ... other paths
```

2. **Edit `config/samples.tsv`**:
```tsv
sample_id	fastq1	fastq2	condition	batch
your_sample_001	/path/to/sample_001_R1.fastq.gz	/path/to/sample_001_R2.fastq.gz	B-ALL	batch1
your_sample_002	/path/to/sample_002_R1.fastq.gz	/path/to/sample_002_R2.fastq.gz	B-ALL	batch1
```

### 7. Test Installation

Run a dry-run to test the pipeline:

```bash
# Activate main environment
conda activate integrateall_python

# Test pipeline configuration
snakemake --dry-run --cores 1

# Run on a small test dataset (if available)
snakemake --cores 4 --use-conda
```

## Advanced Configuration

### Cluster Configuration (SLURM)

Create a cluster configuration file `config/slurm.yaml`:

```yaml
# SLURM cluster configuration
restart-times: 3
jobscript: "workflow/scripts/slurm_jobscript.sh"
cluster: "sbatch --partition={cluster.partition} --account={cluster.account} --cpus-per-task={cluster.cpus} --mem={cluster.mem} --time={cluster.time} --job-name={cluster.name}"
cluster-status: "workflow/scripts/slurm_status.py"
max-jobs-per-second: 1
max-status-checks-per-second: 10
local-cores: 1
latency-wait: 60

default-resources:
  - partition: "normal"
  - account: "your_account"
  - cpus: 1
  - mem: "4000MB" 
  - time: "01:00:00"
  - name: "integrateall_{rule}"
```

### Singularity Configuration

If using Singularity containers:

```bash
# Enable Singularity in config.yaml
singularity:
  use_singularity: true
  singularity_prefix: "/path/to/singularity/images"
```

Run with Singularity:

```bash
snakemake --use-singularity --cores 16
```

## Running the Pipeline

### Local Execution

```bash
# Activate environment
conda activate integrateall_python

# Run pipeline locally
snakemake --cores 16 --use-conda

# Run with specific targets
snakemake results/reports/sample_001/sample_001_final_report.html --cores 8 --use-conda
```

### Cluster Execution

```bash
# Submit to SLURM cluster
snakemake --profile config/slurm --jobs 50

# Monitor progress
squeue -u $USER
```

### Resume Failed Jobs

```bash
# Resume from checkpoint
snakemake --cores 16 --use-conda --rerun-incomplete
```

## Directory Structure After Installation

```
IntegrateALL_pipeline/
├── config/
│   ├── config.yaml              # Main configuration
│   ├── samples.tsv              # Your sample sheet
│   └── slurm.yaml              # Cluster configuration
├── workflow/
│   ├── rules/                   # Snakemake rules
│   ├── scripts/                 # Python analysis scripts
│   ├── envs/                   # Conda environments
│   └── Snakefile               # Main workflow
├── resources/
│   ├── genomes/                # Reference genomes
│   └── databases/              # Tool databases
├── results/                    # Analysis results
├── logs/                       # Execution logs
└── benchmarks/                # Performance data
```

## Troubleshooting

### Common Issues

1. **Out of Memory Errors**
   - Increase memory allocation in `config.yaml`
   - Use fewer parallel jobs

2. **Conda Environment Issues**
   - Clear conda cache: `conda clean --all`
   - Recreate environments: `conda env remove -n integrateall_python`

3. **Reference File Issues**
   - Verify file paths in `config.yaml`
   - Check file permissions

4. **SLURM Issues**
   - Verify cluster configuration
   - Check job limits and quotas

### Getting Help

- Check log files in `logs/` directory
- Use Snakemake's built-in debugging: `snakemake --debug`
- Review benchmark files for performance issues

### Performance Optimization

1. **Resource Allocation**
   - Monitor resource usage with `htop`/`squeue`
   - Adjust thread/memory settings per rule

2. **Storage Optimization**
   - Use fast SSD storage for temporary files
   - Clean up intermediate files regularly

3. **Network Optimization**
   - Use local storage for reference data
   - Minimize network I/O during execution

## Validation

After installation, validate the pipeline with test data:

```bash
# Download test dataset (if available)
wget <test-data-url>

# Run validation
snakemake --cores 4 --configfile config/test_config.yaml
```

The pipeline should complete without errors and produce expected outputs.

## Updates and Maintenance

### Updating the Pipeline

```bash
git pull origin main
conda env update -f workflow/envs/python.yaml
```

### Updating Reference Data

Reference databases should be updated periodically:

```bash
# Update Arriba databases
cd resources/databases/arriba
./update_databases.sh

# Update FusionCatcher database
cd ../fusioncatcher/fusioncatcher-master/data
./download-human-db.sh
```

## Support

For installation support:

1. Check the documentation in `README.md`
2. Review existing issues on GitHub
3. Contact the development team

**Note**: Installation may take several hours due to large reference file downloads. Ensure stable internet connection and sufficient storage space.