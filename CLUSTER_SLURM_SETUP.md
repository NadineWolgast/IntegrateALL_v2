# IntegrateALL_v2 Pipeline - SLURM Cluster Setup
================================================================

## 🚀 **Cluster-Deployment Anleitung**

### 1. **Repository auf Cluster klonen:**
```bash
cd /work_beegfs/suk2m376/
git clone https://github.com/NadineWolgast/IntegrateALL_v2.git
cd IntegrateALL_v2
```

### 2. **Pipeline installieren:**
```bash
# Conda Environment aktivieren (falls vorhanden)
source /work_beegfs/suk2m376/miniconda/bin/activate

# Pipeline installieren
./install_simple.sh

# Pipeline Environment aktivieren
source activate_pipeline.sh
```

### 3. **Sample-Konfiguration anpassen:**
```bash
nano config/samples.tsv
```
**Format:**
```
sample_id	fastq1	fastq2	condition	batch
sample1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz	B-ALL	batch1
sample2	/path/to/sample2_R1.fastq.gz	/path/to/sample2_R2.fastq.gz	B-ALL	batch1
```

### 4. **SLURM Submit Script erstellen:**
```bash
nano submit_integrateall_v2.sh
```

## 📝 **Optimiertes SLURM Submit Script**

```bash
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

# Aktiviere Conda Environment
source /work_beegfs/suk2m376/miniconda/bin/activate
conda activate integrateall_env

# Wechsle in Pipeline Directory
cd /work_beegfs/suk2m376/IntegrateALL_v2

# Pipeline mit SLURM ausführen
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
```

### 5. **Cluster-Konfiguration erstellen:**
```bash
nano config/cluster.yaml
```

```yaml
# SLURM Cluster Configuration
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
```

### 6. **Pipeline starten:**
```bash
# Mache Script ausführbar
chmod +x submit_integrateall_v2.sh

# Job einreichen
sbatch submit_integrateall_v2.sh

# Job Status prüfen
squeue -u suk2m376
```

## 🔍 **Monitoring & Debugging**

### **Pipeline Status prüfen:**
```bash
# Snakemake Status
snakemake --summary

# Fehlgeschlagene Jobs anzeigen
snakemake --list-unfinished

# Dry-run für Debugging
snakemake --dry-run --printshellcmds
```

### **Log Files:**
```bash
# Haupt-Pipeline Log
tail -f integrateall_[JOBID].out
tail -f integrateall_[JOBID].err

# Einzelne Job Logs
ls logs/slurm-*.out
ls logs/slurm-*.err
```

## ⚡ **Performance-Optimierungen**

### **Für große Datasets:**
```bash
# Mehr Ressourcen für STAR
snakemake star_align --cluster "sbatch -c 32 --mem=120G -t 08:00:00"

# Mehr parallele Jobs
snakemake --jobs 64 --cores 128
```

### **Für schnelle Tests:**
```bash
# Nur ein Sample
snakemake --config samples="config/test_sample.tsv" --cores 16
```

## 📊 **Erwartete Laufzeiten**

| Komponente | Kleine Samples | Große Samples |
|------------|----------------|---------------|
| FastQC | 5-10 Min | 20-30 Min |
| STAR Index | 30 Min | 30 Min (einmalig) |
| STAR Align | 10-30 Min | 2-4 Stunden |
| GATK Pipeline | 1-2 Stunden | 4-8 Stunden |
| Arriba | 15-30 Min | 1-2 Stunden |
| Classification | 10-20 Min | 30-60 Min |
| **Gesamt** | **3-5 Stunden** | **8-15 Stunden** |

## 🎯 **Ready für Production!**

Die IntegrateALL_v2 Pipeline ist vollständig für SLURM-Cluster optimiert:

✅ **27-Job Snakemake Workflow**  
✅ **Automatische Conda-Environment Verwaltung**  
✅ **SLURM-native Cluster Integration**  
✅ **Fehlerbehandlung & Restart-Funktionalität**  
✅ **Resource-optimierte Job-Konfiguration**  
✅ **Comprehensive Logging & Monitoring**

**Viel Erfolg bei der B-ALL Analyse auf dem Cluster!** 🧬🚀