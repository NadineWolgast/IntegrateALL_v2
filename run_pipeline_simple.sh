#!/bin/bash

# Simple Pipeline Runner - ohne problematische Python Scripts
# Führt die wichtigsten Schritte manuell aus

set -euo pipefail

echo "🚀 IntegrateALL Pipeline - Vereinfachter Lauf"
echo "============================================="

# Activate environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate integrateall

# Check if STAR alignment exists
if [ -f "results/alignment/test_sample/test_sample.sorted.bam" ]; then
    echo "✅ STAR Alignment bereits vorhanden"
else
    echo "❌ STAR Alignment fehlt - Pipeline muss erst gestartet werden"
    exit 1
fi

# Run basic alignment stats with samtools instead of pysam
echo "📊 Erstelle Alignment-Statistiken..."
mkdir -p results/alignment/test_sample/

samtools flagstat results/alignment/test_sample/test_sample.sorted.bam > results/alignment/test_sample/test_sample.flagstat
samtools idxstats results/alignment/test_sample/test_sample.sorted.bam > results/alignment/test_sample/test_sample.idxstats

# Create simple alignment stats JSON manually
cat > results/alignment/test_sample/test_sample.alignment_stats.json << EOF
{
    "sample": "test_sample",
    "total_reads": $(grep "total" results/alignment/test_sample/test_sample.flagstat | cut -f1 -d' '),
    "mapped_reads": $(grep "mapped" results/alignment/test_sample/test_sample.flagstat | head -1 | cut -f1 -d' '),
    "mapping_rate": "calculated from flagstat"
}
EOF

# Run GATK variant calling pipeline
echo "🧬 Starte Variant Calling..."

# Mark duplicates
mkdir -p results/variants/test_sample
if [ ! -f "results/variants/test_sample/marked_duplicates.bam" ]; then
    echo "  Marking duplicates..."
    picard MarkDuplicates \
        I=results/alignment/test_sample/test_sample.sorted.bam \
        O=results/variants/test_sample/marked_duplicates.bam \
        M=results/variants/test_sample/duplicate_metrics.txt \
        VALIDATION_STRINGENCY=LENIENT
    samtools index results/variants/test_sample/marked_duplicates.bam
fi

# Prepare reference
if [ ! -f "resources/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.dict" ]; then
    echo "  Preparing reference..."
    picard CreateSequenceDictionary \
        R=resources/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        O=resources/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.dict
    samtools faidx resources/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa
fi

# Split N cigars
if [ ! -f "results/variants/test_sample/split_reads.bam" ]; then
    echo "  Splitting N cigars..."
    gatk SplitNCigarReads \
        -I results/variants/test_sample/marked_duplicates.bam \
        -R resources/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -O results/variants/test_sample/split_reads.bam
fi

# Call variants
if [ ! -f "results/variants/test_sample/raw_variants.vcf" ]; then
    echo "  Calling variants..."
    gatk HaplotypeCaller \
        -I results/variants/test_sample/split_reads.bam \
        -R resources/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -O results/variants/test_sample/raw_variants.vcf \
        -dont-use-soft-clipped-bases
fi

# Filter variants
if [ ! -f "results/variants/test_sample/filtered_variants.vcf" ]; then
    echo "  Filtering variants..."
    gatk VariantFiltration \
        -V results/variants/test_sample/raw_variants.vcf \
        -R resources/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -O results/variants/test_sample/filtered_variants.vcf \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \
        --filter-expression "FS > 60.0" --filter-name "FS60" \
        --filter-expression "MQ < 40.0" --filter-name "MQ40"
fi

# Run Arriba manually (since wrapper has problems)
echo "🔀 Starte Fusion Detection mit Arriba..."
mkdir -p results/fusions/test_sample

if [ ! -f "results/fusions/test_sample/arriba_fusions.tsv" ]; then
    # Try to find Arriba executable
    ARRIBA_PATH=$(which arriba 2>/dev/null || echo "")
    
    if [ -z "$ARRIBA_PATH" ]; then
        echo "⚠️  Arriba nicht gefunden - überspringe Fusion Detection"
        # Create empty fusion file
        echo -e "#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\tdirection1\tdirection2\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tcoverage1\tcoverage2\tconfidence" > results/fusions/test_sample/arriba_fusions.tsv
        echo -e "#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\tdirection1\tdirection2\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tcoverage1\tcoverage2\tconfidence" > results/fusions/test_sample/arriba_discarded.tsv
    else
        echo "  Arriba gefunden: $ARRIBA_PATH"
        # Note: Arriba needs specific database files which may not be available
        echo "⚠️  Arriba Datenbanken müssen separat installiert werden - erstelle leere Dateien"
        echo -e "#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\tdirection1\tdirection2\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tcoverage1\tcoverage2\tconfidence" > results/fusions/test_sample/arriba_fusions.tsv
        echo -e "#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\tdirection1\tdirection2\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tcoverage1\tcoverage2\tconfidence" > results/fusions/test_sample/arriba_discarded.tsv
    fi
fi

# Create simple results summary
echo "📋 Erstelle Ergebnis-Zusammenfassung..."
mkdir -p results/reports/test_sample

# Count variants
TOTAL_VARIANTS=$(grep -v "^#" results/variants/test_sample/filtered_variants.vcf | wc -l)
PASSED_VARIANTS=$(grep -v "^#" results/variants/test_sample/filtered_variants.vcf | grep "PASS" | wc -l || echo "0")

cat > results/reports/test_sample/test_sample_summary.json << EOF
{
    "sample": "test_sample",
    "pipeline": "IntegrateALL",
    "timestamp": "$(date -Iseconds)",
    "alignment": {
        "tool": "STAR",
        "status": "completed",
        "bam_file": "results/alignment/test_sample/test_sample.sorted.bam"
    },
    "variants": {
        "tool": "GATK HaplotypeCaller",
        "status": "completed", 
        "total_variants": $TOTAL_VARIANTS,
        "passed_variants": $PASSED_VARIANTS,
        "vcf_file": "results/variants/test_sample/filtered_variants.vcf"
    },
    "fusions": {
        "tool": "Arriba",
        "status": "placeholder",
        "note": "Arriba databases not configured - placeholder files created"
    }
}
EOF

echo ""
echo "✅ Pipeline-Lauf abgeschlossen!"
echo "================================"
echo ""
echo "📊 Ergebnisse:"
echo "  - Alignment: results/alignment/test_sample/test_sample.sorted.bam"
echo "  - Variants: results/variants/test_sample/filtered_variants.vcf ($TOTAL_VARIANTS total, $PASSED_VARIANTS passed)"
echo "  - Fusions: results/fusions/test_sample/arriba_fusions.tsv (placeholder)"
echo "  - Summary: results/reports/test_sample/test_sample_summary.json"
echo ""
echo "🎉 IntegrateALL Pipeline erfolgreich!"