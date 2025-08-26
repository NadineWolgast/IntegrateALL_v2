#!/bin/bash

# IntegrateALL Pipeline - GitHub Repository Setup Script
echo "🐙 IntegrateALL Pipeline - GitHub Setup"
echo "======================================"

# Check if we're in the right directory
if [ ! -f "install_simple.sh" ]; then
    echo "❌ Error: Run this script from the IntegrateALL_pipeline directory"
    exit 1
fi

echo "📝 Creating .gitignore..."
cat > .gitignore << 'EOF'
# Pipeline outputs and temporary files
results/
logs/
benchmarks/
.snakemake/
**/temp/
**/*.tmp
**/*.log

# Large data files (keep structure, ignore content)
*.fastq.gz
*.fastq
*.fq.gz  
*.fq
*.bam
*.bai
*.sam
*.vcf.gz
*.vcf
*.fa
*.fasta
*.fna
*.gtf.gz
*.gtf
*.gff
*.gff3

# Reference databases (too large for git)
resources/genomes/*
resources/databases/*
# But keep directory structure
!resources/genomes/.gitkeep
!resources/databases/.gitkeep

# Conda environments and caches
.snakemake/conda/
**/__pycache__/
*.pyc
*.pyo

# R files
.RData
.Rhistory
*.Rproj

# System files
.DS_Store
Thumbs.db
*.swp
*~

# IDE files
.vscode/
.idea/
*.sublime-*

# Temporary test files
test_*
*_test.*
EOF

echo "📁 Creating directory structure markers..."
touch resources/genomes/.gitkeep
touch resources/databases/.gitkeep

echo "🔧 Initializing Git repository..."
git init
git branch -M main

echo "📦 Adding files to Git..."
git add .

echo "💬 Creating commit message..."
git commit -m "🚀 IntegrateALL Pipeline v1.0 - Cluster Ready

🧬 Comprehensive B-ALL RNA-seq Analysis Pipeline
===============================================

✅ Features:
• One-command installation (./install_simple.sh)
• 27-job Snakemake pipeline for B-ALL analysis  
• STAR alignment + GATK variant calling
• Arriba fusion detection (via Snakemake wrapper)
• ALLCatchR subtype classification
• RNAseqCNV copy number analysis
• Smart installation with --force/--yes options
• Cluster-ready with SLURM support

📊 Pipeline Components:
• 21 Python analysis scripts
• 7 Conda environment specifications  
• Automated reference data download (8GB)
• Interactive HTML reporting
• Quality control and benchmarking

🎯 Deployment:
• Tested on Intel workstation
• Validated STAR + GATK integration
• Dependency conflicts resolved
• Ready for cluster scaling

🔧 Technology Stack:
• Snakemake workflow management
• Conda/Mamba environment management
• Python 3.9 + R 4.x + Bioconductor
• GATK 4.6.2.0 + STAR 2.7.11b
• Modern bioinformatics best practices

🤖 Generated with Claude Code
https://claude.ai/code"

echo ""
echo "✅ Git repository initialized!"
echo ""
echo "📋 Next Steps:"
echo "1. Create private repository on GitHub:"
echo "   → https://github.com/new"
echo "   → Name: IntegrateALL-Pipeline" 
echo "   → ✅ Private"
echo "   → ✅ Add README"
echo ""
echo "2. Connect to your repository:"
echo "   git remote add origin https://github.com/YOUR_USERNAME/IntegrateALL-Pipeline.git"
echo ""
echo "3. Push to GitHub:"
echo "   git push -u origin main"
echo ""  
echo "4. Clone on cluster:"
echo "   git clone https://github.com/YOUR_USERNAME/IntegrateALL-Pipeline.git"
echo "   cd IntegrateALL-Pipeline"
echo "   ./install_simple.sh"
echo ""
echo "🎉 Repository ready for GitHub!"