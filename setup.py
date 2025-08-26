#!/usr/bin/env python3
"""
IntegrateALL Pipeline - Python Setup Script
===========================================

Alternative Python-basiertes Setup für erweiterte Konfiguration und Validierung.
"""

import os
import sys
import subprocess
import shutil
import urllib.request
import tarfile
import json
import yaml
from pathlib import Path
from typing import Dict, List, Optional
import argparse

class PipelineSetup:
    """Automatisches Setup der IntegrateALL Pipeline"""
    
    def __init__(self, pipeline_dir: str = None):
        self.pipeline_dir = Path(pipeline_dir) if pipeline_dir else Path(__file__).parent.absolute()
        self.conda_env_name = "integrateall"
        self.log_file = self.pipeline_dir / "setup.log"
        
        # Ensure log directory exists
        self.log_file.parent.mkdir(exist_ok=True)
        
    def log(self, message: str, level: str = "INFO") -> None:
        """Log message to file and console"""
        timestamp = subprocess.check_output(['date'], text=True).strip()
        log_entry = f"[{timestamp}] {level}: {message}"
        
        print(log_entry)
        with open(self.log_file, 'a') as f:
            f.write(log_entry + '\\n')
    
    def run_command(self, command: List[str], check: bool = True, capture_output: bool = False) -> Optional[str]:
        """Run shell command with error handling"""
        self.log(f"Running: {' '.join(command)}")
        
        try:
            if capture_output:
                result = subprocess.run(command, check=check, capture_output=True, text=True)
                return result.stdout.strip()
            else:
                subprocess.run(command, check=check)
                return None
        except subprocess.CalledProcessError as e:
            self.log(f"Command failed: {e}", "ERROR")
            if check:
                sys.exit(1)
            return None
    
    def check_system_requirements(self) -> None:
        """Check system requirements"""
        self.log("Checking system requirements...")
        
        # Check OS
        if sys.platform != "linux":
            self.log(f"This pipeline requires Linux. Current platform: {sys.platform}", "ERROR")
            sys.exit(1)
        
        # Check Python version
        if sys.version_info < (3, 8):
            self.log(f"Python 3.8+ required. Current: {sys.version}", "ERROR")
            sys.exit(1)
        
        # Check disk space (need at least 50GB)
        statvfs = os.statvfs(self.pipeline_dir)
        available_gb = (statvfs.f_frsize * statvfs.f_bavail) / (1024**3)
        
        if available_gb < 50:
            self.log(f"Warning: Less than 50GB available space ({available_gb:.1f}GB)", "WARNING")
        
        # Check required commands
        required_commands = ['curl', 'wget', 'tar', 'gunzip']
        for cmd in required_commands:
            if not shutil.which(cmd):
                self.log(f"Required command not found: {cmd}", "ERROR")
                sys.exit(1)
        
        self.log("System requirements check passed")
    
    def install_conda(self) -> None:
        """Install Miniconda if not present"""
        if shutil.which('conda'):
            self.log("Conda already installed")
            return
        
        self.log("Installing Miniconda...")
        
        # Download Miniconda
        miniconda_url = "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
        installer_path = "/tmp/miniconda_installer.sh"
        
        try:
            urllib.request.urlretrieve(miniconda_url, installer_path)
            
            # Install Miniconda
            subprocess.run([
                'bash', installer_path, '-b', '-p', f'{Path.home()}/miniconda3'
            ], check=True)
            
            # Initialize conda
            conda_sh = f"{Path.home()}/miniconda3/etc/profile.d/conda.sh"
            if Path(conda_sh).exists():
                self.run_command(['bash', '-c', f'source {conda_sh} && conda init bash'])
            
            # Clean up
            os.remove(installer_path)
            
            self.log("Miniconda installed successfully")
            
        except Exception as e:
            self.log(f"Failed to install Miniconda: {e}", "ERROR")
            sys.exit(1)
    
    def setup_conda_environments(self) -> None:
        """Create conda environments for the pipeline"""
        self.log("Setting up conda environments...")
        
        conda_sh = f"{Path.home()}/miniconda3/etc/profile.d/conda.sh"
        
        # Install mamba for faster package management
        self.run_command([
            'bash', '-c', 
            f'source {conda_sh} && conda install -c conda-forge mamba -y'
        ])
        
        # Create main environment
        env_file = self.pipeline_dir / "workflow/envs/python.yaml"
        if env_file.exists():
            self.run_command([
                'bash', '-c',
                f'source {conda_sh} && mamba env create -f {env_file} -n {self.conda_env_name}'
            ])
        
        # Install Snakemake in main environment
        self.run_command([
            'bash', '-c',
            f'source {conda_sh} && conda activate {self.conda_env_name} && '
            f'mamba install -c bioconda -c conda-forge snakemake=7.32.4 -y'
        ])
        
        self.log("Conda environments created successfully")
    
    def install_r_packages(self) -> None:
        """Install required R packages"""
        self.log("Installing R packages...")
        
        conda_sh = f"{Path.home()}/miniconda3/etc/profile.d/conda.sh"
        
        # Create R environment
        r_env_file = self.pipeline_dir / "workflow/envs/classification.yaml"
        if r_env_file.exists():
            self.run_command([
                'bash', '-c',
                f'source {conda_sh} && mamba env create -f {r_env_file} -n integrateall_r'
            ])
        
        # Install R packages
        r_script = '''
        if (!require('devtools', quietly = TRUE)) {
            install.packages('devtools', repos='https://cran.rstudio.com/')
        }
        devtools::install_github('ThomasBeder/ALLCatchR_bcrabl1')
        devtools::install_github('honzee/RNAseqCNV')
        '''
        
        try:
            self.run_command([
                'bash', '-c',
                f'source {conda_sh} && conda activate integrateall_r && '
                f'Rscript -e "{r_script}"'
            ])
        except:
            self.log("R packages installation may have failed - check manually", "WARNING")
        
        self.log("R packages installation completed")
    
    def download_reference_data(self) -> None:
        """Download reference genome and databases"""
        self.log("Downloading reference data...")
        
        resources_dir = self.pipeline_dir / "resources"
        resources_dir.mkdir(exist_ok=True)
        
        # Create subdirectories
        (resources_dir / "genomes").mkdir(exist_ok=True)
        (resources_dir / "databases").mkdir(exist_ok=True)
        
        self._download_genome(resources_dir / "genomes")
        self._download_databases(resources_dir / "databases")
        
        self.log("Reference data download completed")
    
    def _download_genome(self, genome_dir: Path) -> None:
        """Download human reference genome"""
        genome_file = genome_dir / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        gtf_file = genome_dir / "Homo_sapiens.GRCh38.109.gtf"
        
        if genome_file.exists() and gtf_file.exists():
            self.log("Reference genome already exists")
            return
        
        self.log("Downloading human reference genome...")
        
        # Download genome
        if not genome_file.exists():
            genome_url = "ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
            self.run_command(['wget', '-c', '-O', str(genome_file) + '.gz', genome_url])
            self.run_command(['gunzip', '-f', str(genome_file) + '.gz'])
        
        # Download GTF
        if not gtf_file.exists():
            gtf_url = "ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"
            self.run_command(['wget', '-c', '-O', str(gtf_file) + '.gz', gtf_url])
            self.run_command(['gunzip', '-f', str(gtf_file) + '.gz'])
    
    def _download_databases(self, db_dir: Path) -> None:
        """Download tool databases"""
        
        # Download Arriba databases
        arriba_dir = db_dir / "arriba"
        if not arriba_dir.exists():
            self.log("Downloading Arriba databases...")
            arriba_dir.mkdir(exist_ok=True)
            
            arriba_url = "https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz"
            arriba_tar = arriba_dir / "arriba_v2.4.0.tar.gz"
            
            self.run_command(['wget', '-c', '-O', str(arriba_tar), arriba_url])
            
            # Extract and move database files
            with tarfile.open(arriba_tar, 'r:gz') as tar:
                tar.extractall(arriba_dir)
            
            # Move database files to arriba_dir
            source_db_dir = arriba_dir / "arriba_v2.4.0" / "database"
            if source_db_dir.exists():
                for file in source_db_dir.iterdir():
                    shutil.move(str(file), str(arriba_dir))
            
            # Cleanup
            shutil.rmtree(arriba_dir / "arriba_v2.4.0", ignore_errors=True)
            arriba_tar.unlink()
        
        # Create essential database files
        self._create_essential_databases(db_dir)
    
    def _create_essential_databases(self, db_dir: Path) -> None:
        """Create essential database files"""
        
        # Driver genes
        driver_genes = [
            'ABL1', 'ABL2', 'AFF1', 'BCL2', 'BCL6', 'BCR', 'CREBBP', 'CRLF2',
            'DUX4', 'EBF1', 'EP300', 'ETV6', 'IKZF1', 'KMT2A', 'MLL', 'MYC',
            'NRAS', 'KRAS', 'PAX5', 'PBX1', 'RUNX1', 'TCF3', 'TP53', 'ZNF384'
        ]
        
        with open(db_dir / "driver_genes.txt", 'w') as f:
            f.write('\\n'.join(driver_genes))
        
        # Oncogenes
        oncogenes = ['MYC', 'BCL2', 'BCL6', 'ABL1']
        with open(db_dir / "oncogenes.txt", 'w') as f:
            f.write('\\n'.join(oncogenes))
        
        # Tumor suppressors
        tumor_suppressors = ['TP53', 'RB1', 'CDKN2A']
        with open(db_dir / "tumor_suppressors.txt", 'w') as f:
            f.write('\\n'.join(tumor_suppressors))
    
    def update_configuration(self) -> None:
        """Update configuration files with correct paths"""
        self.log("Updating configuration files...")
        
        config_file = self.pipeline_dir / "config/config.yaml"
        
        if config_file.exists():
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
            
            # Update reference paths
            resources_dir = str(self.pipeline_dir / "resources")
            config['reference_genome'] = f"{resources_dir}/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
            config['reference_gtf'] = f"{resources_dir}/genomes/Homo_sapiens.GRCh38.109.gtf"
            
            # Update database paths
            if 'databases' not in config:
                config['databases'] = {}
            
            config['databases'].update({
                'arriba_blacklist': f"{resources_dir}/databases/arriba/blacklist_hg38_GRCh38_v2.4.0.tsv.gz",
                'arriba_known_fusions': f"{resources_dir}/databases/arriba/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz",
                'arriba_protein_domains': f"{resources_dir}/databases/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3",
                'arriba_cytobands': f"{resources_dir}/databases/arriba/cytobands_hg38_GRCh38_v2.4.0.tsv",
                'driver_gene_list': f"{resources_dir}/databases/driver_genes.txt",
                'oncogene_list': f"{resources_dir}/databases/oncogenes.txt",
                'tumor_suppressor_list': f"{resources_dir}/databases/tumor_suppressors.txt"
            })
            
            # Write updated config
            with open(config_file, 'w') as f:
                yaml.dump(config, f, default_flow_style=False)
        
        self.log("Configuration files updated")
    
    def test_installation(self) -> None:
        """Test pipeline installation"""
        self.log("Testing pipeline installation...")
        
        conda_sh = f"{Path.home()}/miniconda3/etc/profile.d/conda.sh"
        
        # Create minimal test sample sheet
        test_samples = self.pipeline_dir / "config/test_samples.tsv"
        with open(test_samples, 'w') as f:
            f.write("sample_id\\tfastq1\\tfastq2\\tcondition\\tbatch\\n")
            f.write("test_sample\\t/dev/null\\t/dev/null\\tB-ALL\\ttest\\n")
        
        try:
            # Test dry run
            self.run_command([
                'bash', '-c',
                f'cd {self.pipeline_dir} && source {conda_sh} && conda activate {self.conda_env_name} && '
                f'snakemake --dry-run --cores 1 --configfile config/config.yaml'
            ])
            
            self.log("Pipeline test successful")
        except:
            self.log("Pipeline test failed - manual verification needed", "WARNING")
        finally:
            # Cleanup test files
            test_samples.unlink(missing_ok=True)
    
    def create_activation_script(self) -> None:
        """Create pipeline activation script"""
        activation_script = self.pipeline_dir / "activate_pipeline.sh"
        
        script_content = f'''#!/bin/bash
# IntegrateALL Pipeline Activation

# Activate Conda Environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate {self.conda_env_name}

# Change to Pipeline Directory
cd "{self.pipeline_dir}"

echo "=================================================="
echo "  IntegrateALL Pipeline is ready!"
echo "=================================================="
echo ""
echo "Available commands:"
echo "  snakemake --dry-run --cores 1        # Test configuration"
echo "  snakemake --cores 16 --use-conda     # Run pipeline locally"
echo "  snakemake --help                     # Show help"
echo ""
echo "Edit configuration:"
echo "  nano config/config.yaml              # Main configuration"
echo "  nano config/samples.tsv              # Sample list"
echo ""
'''
        
        with open(activation_script, 'w') as f:
            f.write(script_content)
        
        activation_script.chmod(0o755)
        self.log("Activation script created")
    
    def run_full_setup(self) -> None:
        """Run complete setup process"""
        self.log("Starting IntegrateALL Pipeline setup...")
        
        try:
            self.check_system_requirements()
            self.install_conda()
            self.setup_conda_environments()
            self.install_r_packages()
            self.download_reference_data()
            self.update_configuration()
            self.test_installation()
            self.create_activation_script()
            
            print("\\n" + "="*50)
            print("  ✅ Setup completed successfully!")
            print("="*50)
            print(f"\\nTo activate the pipeline:")
            print(f"  source {self.pipeline_dir}/activate_pipeline.sh")
            print(f"\\nNext steps:")
            print(f"1. Edit sample sheet: config/samples.tsv")
            print(f"2. Test pipeline: snakemake --dry-run --cores 1")
            print(f"3. Run pipeline: snakemake --cores 16 --use-conda")
            print(f"\\nSetup log: {self.log_file}")
            
        except Exception as e:
            self.log(f"Setup failed: {e}", "ERROR")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='IntegrateALL Pipeline Setup')
    parser.add_argument('--pipeline-dir', help='Pipeline directory path')
    parser.add_argument('--quick', action='store_true', help='Quick setup without reference downloads')
    
    args = parser.parse_args()
    
    setup = PipelineSetup(args.pipeline_dir)
    
    if args.quick:
        print("Quick setup mode - skipping reference data downloads")
        setup.check_system_requirements()
        setup.install_conda()
        setup.setup_conda_environments()
        setup.create_activation_script()
    else:
        setup.run_full_setup()

if __name__ == "__main__":
    main()