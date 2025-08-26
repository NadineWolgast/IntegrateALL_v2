"""
Classification Rules for IntegrateALL Pipeline
==============================================

B-ALL subtype classification using ALLCatchR and karyotype prediction
"""

# ALLCatchR B-ALL subtype classification
rule allcatchr:
    input:
        counts="results/alignment/{sample}/{sample}.counts_formatted.tsv"
    output:
        predictions="results/classification/{sample}/allcatchr_predictions.tsv",
        scores="results/classification/{sample}/allcatchr_scores.tsv",
        confidence="results/classification/{sample}/allcatchr_confidence.tsv"
    params:
        output_dir="results/classification/{sample}",
        model=config.get("allcatchr_model", "default"),
        extra_args=config.get("allcatchr_args", "")
    resources:
        mem_mb=config.get("allcatchr_memory", 8000),
        runtime=config.get("allcatchr_runtime", 60)
    log:
        "logs/classification/allcatchr_{sample}.log"
    benchmark:
        "benchmarks/classification/allcatchr_{sample}.benchmark.txt"
    conda:
        "../envs/classification.yaml"
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}
        
        # Run ALLCatchR using R script
        Rscript -e "
        library(ALLCatchR)
        
        # Load count data
        counts <- read.table('{input.counts}', header=TRUE, sep='\\t', row.names=1)
        
        # Run ALLCatchR prediction
        result <- ALLCatchR_predict(counts, model='{params.model}')
        
        # Save predictions
        write.table(result\\$predictions, '{output.predictions}', 
                   sep='\\t', quote=FALSE, row.names=FALSE)
        
        # Save scores if available
        if('scores' %in% names(result)) {{
            write.table(result\\$scores, '{output.scores}', 
                       sep='\\t', quote=FALSE, row.names=FALSE)
        }} else {{
            file.create('{output.scores}')
        }}
        
        # Save confidence if available  
        if('confidence' %in% names(result)) {{
            write.table(result\\$confidence, '{output.confidence}', 
                       sep='\\t', quote=FALSE, row.names=FALSE)
        }} else {{
            file.create('{output.confidence}')
        }}
        " &> {log}
        
        # Ensure all output files exist
        touch {output.predictions} {output.scores} {output.confidence}
        """

# Karyotype prediction using machine learning
rule karyotype_prediction:
    input:
        cnv_results="results/cnv/{sample}/rnaseqcnv_results.tsv"
    output:
        prediction="results/classification/{sample}/karyotype_prediction.csv",
        features="results/classification/{sample}/karyotype_features.csv"
    params:
        model_file=config.get("karyotype_model", "resources/models/karyotype_model.pkl"),
        extra_args=config.get("karyotype_args", "")
    log:
        "logs/classification/karyotype_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/predict_karyotype.py"

# Integrative classification combining all results
rule integrated_classification:
    input:
        allcatchr_pred="results/classification/{sample}/allcatchr_predictions.tsv",
        allcatchr_scores="results/classification/{sample}/allcatchr_scores.tsv",
        karyotype_pred="results/classification/{sample}/karyotype_prediction.csv",
        fusion_results="results/fusions/{sample}/driver_fusions.tsv",
        cnv_results="results/cnv/{sample}/rnaseqcnv_results.tsv",
        hotspot_results="results/variants/{sample}/hotspots.csv"
    output:
        classification="results/classification/{sample}/integrated_classification.json",
        confidence="results/classification/{sample}/classification_confidence.json",
        evidence="results/classification/{sample}/classification_evidence.tsv"
    params:
        classification_rules=config.get("classification_rules", "resources/classification_rules.yaml"),
        confidence_threshold=config.get("classification_confidence_threshold", 0.7)
    log:
        "logs/classification/integrated_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/integrated_classification.py"

# Subtype-specific analysis
rule subtype_analysis:
    input:
        classification="results/classification/{sample}/integrated_classification.json",
        gene_expression="results/alignment/{sample}/{sample}.counts_formatted.tsv",
        fusion_results="results/fusions/{sample}/annotated_fusions.tsv"
    output:
        subtype_report="results/classification/{sample}/subtype_analysis.json",
        biomarker_profile="results/classification/{sample}/biomarker_profile.tsv"
    params:
        subtype_markers=config.get("subtype_markers", "resources/databases/subtype_markers.yaml")
    log:
        "logs/classification/subtype_analysis_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/subtype_analysis.py"

# Classification quality control
rule classification_qc:
    input:
        allcatchr_pred="results/classification/{sample}/allcatchr_predictions.tsv",
        allcatchr_scores="results/classification/{sample}/allcatchr_scores.tsv",
        integrated="results/classification/{sample}/integrated_classification.json",
        confidence="results/classification/{sample}/classification_confidence.json"
    output:
        qc_report="results/classification/{sample}/classification_qc.json",
        qc_plot="results/classification/{sample}/classification_qc.png"
    log:
        "logs/classification/qc_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/classification_qc.py"

# Drug response prediction (if enabled)
if config.get("enable_drug_prediction", False):
    rule drug_response_prediction:
        input:
            classification="results/classification/{sample}/integrated_classification.json",
            expression="results/alignment/{sample}/{sample}.counts_formatted.tsv",
            mutations="results/variants/{sample}/filtered_variants.vcf"
        output:
            drug_predictions="results/classification/{sample}/drug_predictions.json",
            drug_report="results/classification/{sample}/drug_response_report.html"
        params:
            drug_db=config.get("drug_database", "resources/databases/drug_targets.db"),
            expression_models=config.get("drug_expression_models", "resources/models/drug_response_models/")
        log:
            "logs/classification/drug_response_{sample}.log"
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/predict_drug_response.py"

# Risk stratification
rule risk_stratification:
    input:
        classification="results/classification/{sample}/integrated_classification.json",
        clinical_data="results/classification/{sample}/clinical_features.json" if config.get("clinical_data") else [],
        molecular_features="results/classification/{sample}/biomarker_profile.tsv"
    output:
        risk_score="results/classification/{sample}/risk_stratification.json",
        risk_report="results/classification/{sample}/risk_report.html"
    params:
        risk_model=config.get("risk_model", "resources/models/risk_model.pkl"),
        risk_thresholds=config.get("risk_thresholds", "resources/risk_thresholds.yaml")
    log:
        "logs/classification/risk_stratification_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/risk_stratification.py"

# Generate classification summary across all samples
rule classification_summary_all:
    input:
        classifications=expand("results/classification/{sample}/integrated_classification.json", sample=SAMPLES),
        confidence_scores=expand("results/classification/{sample}/classification_confidence.json", sample=SAMPLES)
    output:
        summary_table="results/classification/all_samples_classification.tsv",
        summary_plot="results/classification/classification_distribution.png",
        cohort_analysis="results/classification/cohort_analysis.json"
    log:
        "logs/classification/summary_all.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_all_classifications.py"

# Validation against known subtypes (if reference data available)
if config.get("reference_classifications", "") != "":
    rule validate_classification:
        input:
            predictions="results/classification/all_samples_classification.tsv",
            reference=config.get("reference_classifications", "")
        output:
            validation_report="results/classification/validation_report.html",
            confusion_matrix="results/classification/confusion_matrix.png"
        params:
            validation_metrics=config.get("validation_metrics", ["accuracy", "precision", "recall", "f1"])
        log:
            "logs/classification/validation.log"
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/validate_classifications.py"