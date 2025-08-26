"""
Copy Number Variation Analysis Rules for IntegrateALL Pipeline
==============================================================

RNASeqCNV analysis for chromosomal alteration detection in B-ALL
"""

# RNASeqCNV analysis
rule rnaseqcnv:
    input:
        counts="results/alignment/{sample}/{sample}.counts_formatted.tsv",
        vcf="results/variants/{sample}/filtered_variants.vcf"
    output:
        results="results/cnv/{sample}/rnaseqcnv_results.tsv",
        plot="results/cnv/{sample}/cnv_plot.png", 
        log2fc="results/cnv/{sample}/log2_fold_change_per_arm.tsv",
        manual_table="results/cnv/{sample}/manual_annotation_table.tsv",
        estimation_table="results/cnv/{sample}/estimation_table.tsv"
    params:
        output_dir="results/cnv/{sample}",
        dbsnp_file=config.get("dbsnp_rda", "resources/databases/dbSNP_hg38.rda"),
        reference_counts=config.get("rnaseqcnv_reference", "resources/databases/rnaseqcnv_reference_counts.rda"),
        extra_args=config.get("rnaseqcnv_args", "")
    resources:
        mem_mb=config.get("rnaseqcnv_memory", 16000),
        runtime=config.get("rnaseqcnv_runtime", 120)
    log:
        "logs/cnv/rnaseqcnv_{sample}.log"
    benchmark:
        "benchmarks/cnv/rnaseqcnv_{sample}.benchmark.txt"
    conda:
        "../envs/cnv.yaml"
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}
        
        # Run RNASeqCNV analysis using R script
        Rscript -e "
        library(RNAseqCNV)
        library(vcfR)
        
        # Set output directory
        setwd('{params.output_dir}')
        
        # Load gene expression data
        counts <- read.table('{input.counts}', header=TRUE, sep='\\t', row.names=1)
        
        # Process VCF file for SNP data
        if(file.size('{input.vcf}') > 0) {{
            vcf_data <- read.vcfR('{input.vcf}')
            # Convert VCF to format suitable for RNAseqCNV
            snp_data <- vcfR2tidy(vcf_data)
        }} else {{
            snp_data <- NULL
        }}
        
        # Load dbSNP data
        if(file.exists('{params.dbsnp_file}')) {{
            load('{params.dbsnp_file}')
        }}
        
        # Run RNASeqCNV analysis
        sample_id <- '{wildcards.sample}'
        
        # Prepare count data (ensure proper formatting)
        count_matrix <- as.matrix(counts)
        colnames(count_matrix) <- sample_id
        
        # Run CNV estimation
        cnv_result <- RNAseqCNV(
            counts = count_matrix,
            genome = 'hg38',
            sample.name = sample_id,
            snp.data = snp_data,
            plot = TRUE,
            {params.extra_args}
        )
        
        # Save main results
        write.table(cnv_result\\$estimation_table, 
                   '{output.estimation_table}', 
                   sep='\\t', quote=FALSE, row.names=FALSE)
        
        # Save log2 fold change per chromosome arm
        if('log2fc_per_arm' %in% names(cnv_result)) {{
            write.table(cnv_result\\$log2fc_per_arm, 
                       '{output.log2fc}', 
                       sep='\\t', quote=FALSE)
        }}
        
        # Save manual annotation table
        if('manual_an_table' %in% names(cnv_result)) {{
            write.table(cnv_result\\$manual_an_table, 
                       '{output.manual_table}', 
                       sep='\\t', quote=FALSE, row.names=FALSE)
        }}
        
        # Create summary results file
        cnv_summary <- data.frame(
            sample_id = sample_id,
            total_alterations = nrow(cnv_result\\$estimation_table),
            gains = sum(cnv_result\\$estimation_table\\$call == 'gain', na.rm=TRUE),
            losses = sum(cnv_result\\$estimation_table\\$call == 'loss', na.rm=TRUE),
            neutral = sum(cnv_result\\$estimation_table\\$call == 'neutral', na.rm=TRUE)
        )
        
        write.table(cnv_summary, '{output.results}', 
                   sep='\\t', quote=FALSE, row.names=FALSE)
        
        # Copy plot if it exists
        if(file.exists(paste0(sample_id, '_CNV_main_fig.png'))) {{
            file.copy(paste0(sample_id, '_CNV_main_fig.png'), '{output.plot}')
        }}
        " &> {log}
        
        # Ensure all output files exist
        touch {output.results} {output.plot} {output.log2fc} {output.manual_table} {output.estimation_table}
        """

# Annotate CNV results with gene information
rule annotate_cnvs:
    input:
        cnv_results="results/cnv/{sample}/rnaseqcnv_results.tsv",
        estimation_table="results/cnv/{sample}/estimation_table.tsv"
    output:
        annotated_cnvs="results/cnv/{sample}/annotated_cnvs.tsv",
        gene_cnv_map="results/cnv/{sample}/gene_cnv_mapping.tsv"
    params:
        gene_bed=config.get("gene_coordinates", "resources/databases/gene_coordinates.bed"),
        oncogenes=config.get("oncogene_list", "resources/databases/oncogenes.txt"),
        tumor_suppressors=config.get("tumor_suppressor_list", "resources/databases/tumor_suppressors.txt")
    log:
        "logs/cnv/annotate_cnvs_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/annotate_cnvs.py"

# CNV quality control
rule cnv_qc:
    input:
        cnv_results="results/cnv/{sample}/rnaseqcnv_results.tsv",
        estimation_table="results/cnv/{sample}/estimation_table.tsv",
        alignment_stats="results/alignment/{sample}/{sample}.alignment_stats.json"
    output:
        qc_report="results/cnv/{sample}/cnv_qc.json",
        qc_plot="results/cnv/{sample}/cnv_qc.png"
    params:
        min_coverage=config.get("cnv_min_coverage", 10),
        max_noise_level=config.get("cnv_max_noise", 0.3)
    log:
        "logs/cnv/cnv_qc_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/cnv_qc.py"

# Chromosomal instability analysis
rule chromosomal_instability:
    input:
        cnv_results="results/cnv/{sample}/rnaseqcnv_results.tsv",
        log2fc="results/cnv/{sample}/log2_fold_change_per_arm.tsv"
    output:
        instability_score="results/cnv/{sample}/chromosomal_instability.json",
        instability_plot="results/cnv/{sample}/chromosomal_instability.png"
    params:
        instability_threshold=config.get("instability_threshold", 0.2)
    log:
        "logs/cnv/chromosomal_instability_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/calculate_chromosomal_instability.py"

# Focal vs broad CNV classification
rule classify_cnv_type:
    input:
        estimation_table="results/cnv/{sample}/estimation_table.tsv"
    output:
        cnv_classification="results/cnv/{sample}/cnv_type_classification.tsv",
        cnv_summary="results/cnv/{sample}/cnv_type_summary.json"
    params:
        focal_threshold=config.get("focal_cnv_threshold", 3000000),  # 3Mb
        broad_threshold=config.get("broad_cnv_threshold", 50000000)  # 50Mb
    log:
        "logs/cnv/classify_cnv_type_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/classify_cnv_types.py"

# Ploidy estimation
rule estimate_ploidy:
    input:
        log2fc="results/cnv/{sample}/log2_fold_change_per_arm.tsv",
        estimation_table="results/cnv/{sample}/estimation_table.tsv"
    output:
        ploidy_estimate="results/cnv/{sample}/ploidy_estimate.json",
        ploidy_plot="results/cnv/{sample}/ploidy_distribution.png"
    log:
        "logs/cnv/estimate_ploidy_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/estimate_ploidy.py"

# CNV signature analysis
rule cnv_signatures:
    input:
        cnv_results="results/cnv/{sample}/rnaseqcnv_results.tsv",
        classification="results/classification/{sample}/integrated_classification.json"
    output:
        signature_profile="results/cnv/{sample}/cnv_signature.json",
        signature_plot="results/cnv/{sample}/cnv_signature.png"
    params:
        signature_db=config.get("cnv_signatures", "resources/databases/cnv_signatures.yaml")
    log:
        "logs/cnv/cnv_signatures_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/analyze_cnv_signatures.py"

# Combine CNV results across all samples
rule cnv_summary_all:
    input:
        cnv_results=expand("results/cnv/{sample}/rnaseqcnv_results.tsv", sample=SAMPLES),
        instability_scores=expand("results/cnv/{sample}/chromosomal_instability.json", sample=SAMPLES),
        ploidy_estimates=expand("results/cnv/{sample}/ploidy_estimate.json", sample=SAMPLES)
    output:
        summary_table="results/cnv/all_samples_cnv_summary.tsv",
        summary_plot="results/cnv/cohort_cnv_analysis.png",
        cohort_instability="results/cnv/cohort_chromosomal_instability.json"
    log:
        "logs/cnv/cnv_summary_all.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_all_cnvs.py"

# CNV-based clustering
rule cnv_clustering:
    input:
        summary_table="results/cnv/all_samples_cnv_summary.tsv"
    output:
        clusters="results/cnv/cnv_based_clusters.tsv",
        cluster_plot="results/cnv/cnv_clustering.png",
        cluster_heatmap="results/cnv/cnv_cluster_heatmap.png"
    params:
        clustering_method=config.get("cnv_clustering_method", "hierarchical"),
        n_clusters=config.get("cnv_n_clusters", "auto")
    log:
        "logs/cnv/cnv_clustering.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/cluster_by_cnv.py"