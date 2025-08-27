# Original Classification Rules - Maintains research validity
# ==========================================================
# Based on original Blast-o-Matic-Fusioninator system
# with performance optimizations

# SNV Hotspot Detection using pysamstats
rule detect_snv_hotspots:
    input:
        bam="results/alignment/{sample}/{sample}.sorted.bam",
        bai="results/alignment/{sample}/{sample}.sorted.bam.bai",
        reference=config["reference_genome"]
    output:
        summary="results/variants/{sample}/snv_hotspots_summary.csv"
    params:
        hotspot_dir="results/variants/{sample}/hotspots"
    log:
        "logs/variants/snv_hotspots_{sample}.log"
    conda:
        "../envs/variants.yaml"  
    script:
        "../scripts/detect_snv_hotspots.py"

# Final Classification using original algorithm
rule final_classification:
    input:
        allcatchr_pred="results/classification/{sample}/allcatchr_predictions.tsv",
        karyotype_pred="results/classification/{sample}/karyotype_prediction.csv",
        fusioncatcher_results="results/fusions/{sample}/fusioncatcher_results.txt",
        arriba_results="results/fusions/{sample}/arriba_fusions.tsv",
        hotspot_summary="results/variants/{sample}/snv_hotspots_summary.csv",
        classification_rules="resources/databases/Class_test.csv"
    output:
        final_report="results/classification/{sample}/final_classification_report.csv",
        driver_summary="results/classification/{sample}/driver_fusions.csv", 
        curation_summary="results/classification/{sample}/curation_summary.csv"
    log:
        "logs/classification/final_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/make_final_classification_optimized.py"

# Generate Interactive HTML Report (original format)
rule generate_interactive_report:
    input:
        final_classification="results/classification/{sample}/final_classification_report.csv",
        allcatchr_pred="results/classification/{sample}/allcatchr_predictions.tsv",
        fusion_results="results/fusions/{sample}/arriba_fusions.tsv",
        cnv_results="results/cnv/{sample}/rnaseqcnv_results.tsv",
        variant_results="results/variants/{sample}/filtered_variants.vcf",
        hotspot_summary="results/variants/{sample}/snv_hotspots_summary.csv",
        qc_stats="results/alignment/{sample}/{sample}.flagstat"
    output:
        html_report="results/reports/{sample}/{sample}_interactive_report.html",
        json_summary="results/reports/{sample}/{sample}_classification_summary.json"
    params:
        sample_id="{sample}",
        template="resources/templates/interactive_report_template.html"
    log:
        "logs/reports/interactive_{sample}.log"
    conda:
        "../envs/python.yaml" 
    script:
        "../scripts/generate_interactive_report.py"

# Combined summary for all samples
rule classification_summary:
    input:
        expand("results/classification/{sample}/final_classification_report.csv", sample=SAMPLES)
    output:
        combined_summary="results/reports/all_samples_classification_summary.csv",
        statistics="results/reports/classification_statistics.json"
    log:
        "logs/reports/classification_summary.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_all_classifications.py"