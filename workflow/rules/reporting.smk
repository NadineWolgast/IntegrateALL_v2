"""
Reporting Rules for IntegrateALL Pipeline
==========================================

Generate comprehensive HTML reports and summaries
"""

# Generate sample-specific final report
rule generate_sample_report:
    input:
        qc_data="results/qc/multiqc/{sample}_multiqc_data.json",
        alignment_stats="results/alignment/{sample}/{sample}.alignment_stats.json",
        fusion_summary="results/fusions/{sample}/fusion_summary.json",
        cnv_results="results/cnv/{sample}/rnaseqcnv_results.tsv",
        classification="results/classification/{sample}/integrated_classification.json",
        variant_summary="results/variants/{sample}/hotspot_summary.json"
    output:
        html_report="results/reports/{sample}/{sample}_final_report.html",
        summary_json="results/reports/{sample}/{sample}_summary.json"
    params:
        sample_id="{sample}",
        output_dir="results/reports/{sample}",
        template_dir="workflow/templates"
    log:
        "logs/reports/sample_report_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_final_report.py"

# Create pipeline-wide summary report
rule generate_pipeline_summary:
    input:
        sample_summaries=expand("results/reports/{sample}/{sample}_summary.json", sample=SAMPLES),
        fusion_summary="results/fusions/all_samples_fusion_summary.tsv",
        classification_summary="results/classification/all_samples_classification.tsv",
        cnv_summary="results/cnv/all_samples_cnv_summary.tsv",
        variant_summary="results/variants/all_samples_variants.tsv"
    output:
        html_report="results/reports/pipeline_summary.html",
        json_summary="results/reports/pipeline_summary.json"
    params:
        pipeline_version=config.get("pipeline_version", "1.0.0"),
        analysis_date=lambda wildcards: __import__('datetime').datetime.now().strftime('%Y-%m-%d')
    log:
        "logs/reports/pipeline_summary.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_pipeline_summary.py"

# Quality control summary across all samples
rule qc_summary_report:
    input:
        qc_metrics=expand("results/qc/metrics/{sample}_qc_metrics.json", sample=SAMPLES),
        alignment_stats=expand("results/alignment/{sample}/{sample}.alignment_stats.json", sample=SAMPLES)
    output:
        html_report="results/reports/qc_summary.html",
        csv_summary="results/reports/qc_metrics_summary.csv"
    log:
        "logs/reports/qc_summary.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_qc_summary.py"

# Classification performance report
rule classification_performance:
    input:
        classifications=expand("results/classification/{sample}/integrated_classification.json", sample=SAMPLES),
        confidence_scores=expand("results/classification/{sample}/classification_confidence.json", sample=SAMPLES)
    output:
        performance_report="results/reports/classification_performance.html",
        confidence_distribution="results/reports/confidence_distribution.png"
    log:
        "logs/reports/classification_performance.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/classification_performance_report.py"

# Fusion discovery summary
rule fusion_discovery_report:
    input:
        fusion_summaries=expand("results/fusions/{sample}/fusion_summary.json", sample=SAMPLES),
        driver_fusions=expand("results/fusions/{sample}/driver_fusions.tsv", sample=SAMPLES)
    output:
        fusion_report="results/reports/fusion_discovery.html",
        recurrent_fusions="results/reports/recurrent_fusions.csv"
    log:
        "logs/reports/fusion_discovery.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fusion_discovery_report.py"

# CNV landscape report
rule cnv_landscape_report:
    input:
        cnv_summaries=expand("results/cnv/{sample}/rnaseqcnv_results.tsv", sample=SAMPLES),
        instability_scores=expand("results/cnv/{sample}/chromosomal_instability.json", sample=SAMPLES)
    output:
        cnv_report="results/reports/cnv_landscape.html",
        cnv_heatmap="results/reports/cnv_heatmap.png"
    log:
        "logs/reports/cnv_landscape.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/cnv_landscape_report.py"

# Create interactive dashboard (optional)
if config.get("enable_dashboard", False):
    rule interactive_dashboard:
        input:
            pipeline_summary="results/reports/pipeline_summary.json",
            sample_data=expand("results/reports/{sample}/{sample}_summary.json", sample=SAMPLES)
        output:
            dashboard="results/reports/interactive_dashboard.html"
        params:
            dashboard_config=config.get("dashboard_config", "workflow/templates/dashboard_config.yaml")
        log:
            "logs/reports/interactive_dashboard.log"
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/create_interactive_dashboard.py"

# Generate publication-ready figures
if config.get("generate_publication_figures", False):
    rule publication_figures:
        input:
            classifications="results/classification/all_samples_classification.tsv",
            fusions="results/fusions/all_samples_fusion_summary.tsv",
            cnvs="results/cnv/all_samples_cnv_summary.tsv"
        output:
            figures_dir=directory("results/reports/publication_figures")
        params:
            figure_format=config.get("figure_format", "png"),
            figure_dpi=config.get("figure_dpi", 300)
        log:
            "logs/reports/publication_figures.log"
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/generate_publication_figures.py"

# Export data for external analysis
rule export_data:
    input:
        gene_counts="results/alignment/merged_gene_counts.tsv",
        classifications="results/classification/all_samples_classification.tsv",
        fusions="results/fusions/all_samples_fusion_summary.tsv",
        variants="results/variants/all_samples_variants.tsv",
        cnvs="results/cnv/all_samples_cnv_summary.tsv"
    output:
        export_package="results/reports/integrateall_export.tar.gz",
        metadata="results/reports/export_metadata.json"
    params:
        export_format=config.get("export_format", "standard"),
        include_raw_data=config.get("export_include_raw", False)
    log:
        "logs/reports/export_data.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/export_analysis_data.py"

# Performance benchmarking report
if config.get("benchmark", False):
    rule benchmark_report:
        input:
            benchmarks=expand("benchmarks/{step}_{sample}.benchmark.txt", 
                             step=["star_align", "arriba", "fusioncatcher", "gatk_call", "allcatchr"], 
                             sample=SAMPLES)
        output:
            benchmark_report="results/reports/benchmark_summary.html",
            performance_metrics="results/reports/performance_metrics.json"
        log:
            "logs/reports/benchmark_report.log"
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/summarize_benchmarks.py"

# Email notification (if configured)
if config.get("send_notifications", False) and config.get("notification_recipients", []):
    rule send_notification:
        input:
            pipeline_summary="results/reports/pipeline_summary.html"
        output:
            notification_log="results/reports/notification_sent.log"
        params:
            email_config=config.get("email_config", {}),
            notification_recipients=config.get("notification_recipients", [])
        log:
            "logs/reports/send_notification.log"
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/send_notification.py"

# Archive results for long-term storage
if config.get("auto_archive", False):
    rule archive_results:
        input:
            pipeline_summary="results/reports/pipeline_summary.html",
            export_package="results/reports/integrateall_export.tar.gz"
        output:
            archive="results/reports/archive_{date}.tar.gz"
        params:
            archive_config=config.get("archive_config", {}),
            date=lambda wildcards: wildcards.date if hasattr(wildcards, 'date') else __import__('datetime').datetime.now().strftime('%Y%m%d')
        log:
            "logs/reports/archive_results.log"
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/archive_results.py"