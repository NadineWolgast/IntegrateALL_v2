"""
Common functions and configuration validation for IntegrateALL pipeline
=====================================================================
"""

import pandas as pd
import os
from pathlib import Path
import yaml

def load_samples(samples_file):
    """
    Load and validate sample sheet
    
    Expected format:
    sample_id,fastq1,fastq2,condition
    """
    try:
        samples_df = pd.read_csv(samples_file, sep='\t')
    except Exception as e:
        print(f"Error loading samples file {samples_file}: {e}")
        sys.exit(1)
    
    # Validate required columns
    required_columns = ['sample_id', 'fastq1', 'fastq2']
    missing_columns = set(required_columns) - set(samples_df.columns)
    if missing_columns:
        print(f"Missing required columns in samples file: {missing_columns}")
        sys.exit(1)
    
    # Check for duplicate sample IDs
    if samples_df['sample_id'].duplicated().any():
        duplicates = samples_df[samples_df['sample_id'].duplicated()]['sample_id'].tolist()
        print(f"Duplicate sample IDs found: {duplicates}")
        sys.exit(1)
    
    # Validate FASTQ files exist
    for _, row in samples_df.iterrows():
        for fastq_col in ['fastq1', 'fastq2']:
            fastq_path = row[fastq_col]
            if not os.path.exists(fastq_path):
                print(f"FASTQ file not found: {fastq_path} for sample {row['sample_id']}")
                sys.exit(1)
    
    print(f"Loaded {len(samples_df)} samples from {samples_file}")
    return samples_df

def validate_config(config):
    """Validate pipeline configuration"""
    
    # Required configuration parameters
    required_params = [
        'samples',
        'reference_genome',
        'reference_gtf',
        'output_dir'
    ]
    
    missing_params = []
    for param in required_params:
        if param not in config:
            missing_params.append(param)
    
    if missing_params:
        print(f"Missing required configuration parameters: {missing_params}")
        sys.exit(1)
    
    # Validate reference files exist
    ref_files = {
        'reference_genome': config.get('reference_genome'),
        'reference_gtf': config.get('reference_gtf')
    }
    
    # Skip validation for installation/download rules
    import sys
    skip_validation_targets = [
        'install_all', 'download_references', 'install_conda_setup',
        'install_r_packages', 'install_references', 'download_genome',
        'download_gtf', 'create_star_index'
    ]
    
    # Check if we're running installation targets
    if len(sys.argv) > 1:
        target_rules = [arg for arg in sys.argv if not arg.startswith('-')]
        skip_validation = any(target in skip_validation_targets for target in target_rules)
    else:
        skip_validation = False
    
    if not skip_validation:
        for name, path in ref_files.items():
            if path and not os.path.exists(path):
                print(f"Reference file not found - {name}: {path}")
                print(f"Run 'snakemake --cores 4 --use-conda download_references' first")
                sys.exit(1)
    
    # Set default values
    config.setdefault('threads', 8)
    config.setdefault('memory_gb', 32)
    config.setdefault('cleanup_temp', True)
    config.setdefault('benchmark', False)
    
    return config

def get_fastq_files(wildcards):
    """Get FASTQ files for a sample"""
    sample_row = samples_df[samples_df['sample_id'] == wildcards.sample].iloc[0]
    return {
        'fastq1': sample_row['fastq1'],
        'fastq2': sample_row['fastq2']
    }

def get_all_outputs():
    """Get all expected output files for the pipeline"""
    outputs = []
    
    for sample in SAMPLES:
        # Quality control outputs
        outputs.extend([
            f"results/qc/fastqc/{sample}_R1_fastqc.html",
            f"results/qc/fastqc/{sample}_R2_fastqc.html",
            f"results/qc/multiqc/{sample}_multiqc_report.html"
        ])
        
        # Alignment outputs
        outputs.extend([
            f"results/alignment/{sample}/{sample}.sorted.bam",
            f"results/alignment/{sample}/{sample}.flagstat",
            f"results/alignment/{sample}/{sample}.gene_counts.tsv"
        ])
        
        # Fusion detection outputs  
        outputs.extend([
            f"results/fusions/{sample}/arriba_fusions.tsv",
            f"results/fusions/{sample}/fusioncatcher_results.txt",
            f"results/fusions/{sample}/fusion_summary.json"
        ])
        
        # Classification outputs
        outputs.extend([
            f"results/classification/{sample}/allcatchr_predictions.tsv", 
            f"results/classification/{sample}/karyotype_prediction.csv"
        ])
        
        # CNV outputs
        outputs.extend([
            f"results/cnv/{sample}/rnaseqcnv_results.tsv",
            f"results/cnv/{sample}/cnv_plot.png"
        ])
        
        # Variant calling outputs
        outputs.extend([
            f"results/variants/{sample}/raw_variants.vcf",
            f"results/variants/{sample}/filtered_variants.vcf", 
            f"results/variants/{sample}/hotspots.csv"
        ])
        
        # Final report
        outputs.extend([
            f"results/reports/{sample}/{sample}_final_report.html",
            f"results/reports/{sample}/{sample}_summary.json"
        ])
    
    # Pipeline-wide outputs
    outputs.append("results/reports/pipeline_summary.html")
    
    return outputs

def create_temp_dir(sample=None):
    """Create sample-specific temporary directory"""
    if sample:
        temp_dir = f"results/temp/{sample}"
    else:
        temp_dir = "results/temp"
    
    os.makedirs(temp_dir, exist_ok=True)
    return temp_dir

def log_rule_start(rule_name, sample=None):
    """Log rule execution start"""
    timestamp = pd.Timestamp.now().isoformat()
    sample_info = f" for sample {sample}" if sample else ""
    print(f"[{timestamp}] Starting rule {rule_name}{sample_info}")

def log_rule_end(rule_name, sample=None):
    """Log rule execution end"""
    timestamp = pd.Timestamp.now().isoformat()
    sample_info = f" for sample {sample}" if sample else ""
    print(f"[{timestamp}] Completed rule {rule_name}{sample_info}")

# Validate configuration on import
config = validate_config(config)

# Load samples
samples_df = load_samples(config["samples"])
SAMPLES = samples_df["sample_id"].tolist()

print(f"IntegrateALL Pipeline initialized with {len(SAMPLES)} samples")