# IntegrateALL_v2 Pipeline 🧬

## 🚀 **Smart Installation mit Auto-Konfiguration**

### ⚡ **Ein-Befehl Setup**

```bash
# Repository klonen
git clone https://github.com/NadineWolgast/IntegrateALL_v2.git
cd IntegrateALL_v2

# Pipeline installieren (vollautomatisch mit Pfad-Erkennung!)
./install_simple.sh

# Pipeline aktivieren
source activate_pipeline.sh
```

**Das war's!** Die Pipeline ist vollständig installiert und konfiguriert! 🎉

## ✨ **Automatische Konfiguration**

Die Installation erkennt automatisch:
- 🔍 **Pipeline-Verzeichnis** (aktueller Pfad)
- 🐍 **Conda-Installation** (Miniconda/Anaconda)  
- 🏠 **Benutzer-Home-Verzeichnis**
- 📝 **SLURM Submit-Scripts** werden generiert
- ⚙️ **Cluster-Konfiguration** wird erstellt

### 📁 **Automatisch erstellte Dateien:**
- `submit_pipeline.sh` - SLURM Batch-Script mit korrekten Pfaden
- `config/cluster.yaml` - Optimierte Ressourcen-Konfiguration  
- `activate_pipeline.sh` - Environment-Aktivierung
- `config/samples_template.tsv` - Sample-Konfiguration Template

## 🔧 **Manuelle Pfad-Anpassung (optional)**

Falls du spezielle Pfade brauchst:

```bash
# Interaktive Konfiguration
./setup_paths.sh

# Direkte Pfad-Angabe
./setup_paths.sh /work_beegfs/user/pipeline /opt/conda

# Nur Basis-Pfad setzen
./setup_paths.sh /custom/pipeline/path
```

## 📊 **Sample-Konfiguration**

```bash
# Template kopieren und bearbeiten
cp config/samples_template.tsv config/samples.tsv
nano config/samples.tsv
```

**Format:**
```tsv
sample_id	fastq1	fastq2	condition	batch
sample1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz	B-ALL	batch1
sample2	/path/to/sample2_R1.fastq.gz	/path/to/sample2_R2.fastq.gz	B-ALL	batch1
```

## 🚀 **Pipeline ausführen**

### **SLURM-Cluster:**
```bash
# Job einreichen
sbatch submit_pipeline.sh

# Status prüfen
squeue -u $(whoami)

# Logs verfolgen  
tail -f integrateall_*.out
```

### **Lokal:**
```bash
# Environment aktivieren
source activate_pipeline.sh

# Pipeline starten
snakemake --cores 16 --use-conda

# Dry-run (Test)
snakemake --dry-run --printshellcmds
```

## 🌍 **Multi-Environment Support**

Die Pipeline funktioniert automatisch auf:

### **Workstations/Laptops:**
- ✅ Lokale Conda-Installation
- ✅ Automatische Pfad-Erkennung
- ✅ Resource-Anpassung

### **SLURM-Cluster:**  
- ✅ Auto-generierte Batch-Scripts
- ✅ Cluster-Ressourcen-Konfiguration
- ✅ Job-Dependency-Management

### **Cloud-Systeme:**
- ✅ Pfad-agnostische Installation  
- ✅ Container-ready
- ✅ Skalierbare Konfiguration

## 📋 **Pipeline-Komponenten (27 Jobs)**

### **Quality Control & Alignment**
- FastQC, MultiQC für Input-QC
- STAR alignment mit 2-pass mapping
- BAM indexing und QC-Statistiken

### **Variant Calling (GATK Pipeline)**
- AddOrReplaceReadGroups (für STAR-Output)
- MarkDuplicates  
- SplitNCigarReads (RNA-spezifisch)
- HaplotypeCaller
- VariantFiltration

### **Fusion Detection**
- Arriba via offiziellen Snakemake Wrapper
- Automatische Datenbank-Downloads
- Fusion-Integration und Annotation

### **B-ALL Classification**
- ALLCatchR (Genexpressions-basierte Klassifikation)
- B-ALL Subtyp-Vorhersage  
- Konfidenz-Scores

### **CNV Analysis**
- RNAseqCNV (Copy Number Variation Detection)
- Chromosomale Instabilität
- Ploidie-Schätzung

### **Reporting & Summary**
- Sample-spezifische Reports
- Pipeline-Summary
- Interactive HTML Reports

## 🎯 **Für jeden User geeignet!**

**Anfänger:** Ein-Befehl Installation, alles automatisch  
**Fortgeschrittene:** Flexible Pfad-Konfiguration  
**Cluster-Admins:** Vorkonfigurierte SLURM-Integration  
**Entwickler:** Modulare, erweiterbare Architektur

Die IntegrateALL_v2 Pipeline ist **universell deployment-ready**! 🌟