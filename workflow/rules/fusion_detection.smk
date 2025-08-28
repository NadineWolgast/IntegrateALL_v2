"""
Fusion Detection Rules for IntegrateALL Pipeline
================================================

Arriba and FusionCatcher integration for comprehensive fusion detection
"""

# Arriba fusion detection
rule arriba:
    input:
        bam="results/alignment/{sample}/{sample}.sorted.bam",
        genome=config["reference_genome"],
        annotation=config["reference_gtf"]
    output:
        fusions="results/fusions/{sample}/arriba_fusions.tsv",
        discarded="results/fusions/{sample}/arriba_discarded.tsv"
    params:
        genome_build="GRCh38",
        default_blacklist=True,
        default_known_fusions=True,
        sv_file="",
        extra=""
    threads:
        config.get("arriba_threads", 4)
    resources:
        mem_mb=config.get("arriba_memory", 16000),
        runtime=config.get("arriba_runtime", 120)
    log:
        "logs/fusions/arriba_{sample}.log"
    benchmark:
        "benchmarks/fusions/arriba_{sample}.benchmark.txt"
    wrapper:
        "v3.10.2/bio/arriba"

# Arriba visualization
rule arriba_visualization:
    input:
        fusions="results/fusions/{sample}/arriba_fusions.tsv",
        bam="results/alignment/{sample}/{sample}.sorted.bam",
        bai="results/alignment/{sample}/{sample}.sorted.bam.bai",
        gtf=config["reference_gtf"]
    output:
        pdf="results/fusions/{sample}/arriba_fusions.pdf"
    params:
        cytobands=config.get("arriba_cytobands", "resources/databases/arriba/cytobands_hg38_GRCh38_v2.4.0.tsv"),
        protein_domains=config.get("arriba_protein_domains", "resources/databases/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3"),
        extra_args=config.get("arriba_draw_args", "")
    log:
        "logs/fusions/arriba_visualization_{sample}.log"
    benchmark:
        "benchmarks/fusions/arriba_visualization_{sample}.benchmark.txt"
    conda:
        "../envs/fusions.yaml"
    shell:
        """
        # Check if fusions file has content (more than just header)
        if [ $(wc -l < {input.fusions}) -gt 1 ]; then
            draw_fusions.R \\
                --fusions={input.fusions} \\
                --alignments={input.bam} \\
                --output={output.pdf} \\
                --annotation={input.gtf} \\
                --cytobands={params.cytobands} \\
                --proteinDomains={params.protein_domains} \\
                {params.extra_args} \\
                &> {log}
        else
            # Create empty PDF if no fusions to visualize
            echo "No fusions found - creating empty PDF" > {log}
            touch {output.pdf}
        fi
        """

# FusionCatcher detection
rule fusioncatcher:
    input:
        fastq1=lambda wildcards: samples_df[samples_df['sample_id'] == wildcards.sample]['fastq1'].iloc[0],
        fastq2=lambda wildcards: samples_df[samples_df['sample_id'] == wildcards.sample]['fastq2'].iloc[0]
    output:
        results="results/fusions/{sample}/fusioncatcher_results.txt",
        dir=directory("results/fusions/{sample}/fusioncatcher_output")
    params:
        data_dir=config.get("fusioncatcher_data", "resources/databases/fusioncatcher/human_v102"),
        extra_args=config.get("fusioncatcher_args", "")
    threads:
        config.get("fusioncatcher_threads", 16)
    resources:
        mem_mb=config.get("fusioncatcher_memory", 64000),
        runtime=config.get("fusioncatcher_runtime", 480)
    log:
        "logs/fusions/fusioncatcher_{sample}.log"
    benchmark:
        "benchmarks/fusions/fusioncatcher_{sample}.benchmark.txt"
    conda:
        "../envs/fusions.yaml"
    shell:
        """
        # Create output directory
        mkdir -p {output.dir}
        
        # Run FusionCatcher (EXACT original parameters)
        fusioncatcher \\
            -d {params.data_dir} \\
            -i {input.fastq1},{input.fastq2} \\
            -o {output.dir} \\
            {params.extra_args} \\
            &> {log}
        
        # Copy main results file to expected output location
        if [ -f "{output.dir}/final-list_candidate-fusion-genes.txt" ]; then
            cp {output.dir}/final-list_candidate-fusion-genes.txt {output.results}
        else
            # Create empty results file if no fusions found
            echo -e "Gene_1_symbol(5end_fusion_partner)\\tGene_2_symbol(3end_fusion_partner)\\tFusion_description\\tCounts_of_common_mapping_reads\\tSpanning_pairs\\tSpanning_unique_reads\\tLongest_anchor_found\\tFusion_finding_method\\tFusion_point_for_gene_1(5end_fusion_partner)\\tFusion_point_for_gene_2(3end_fusion_partner)" > {output.results}
        fi
        """

# Integrate fusion results from both tools
rule integrate_fusions:
    input:
        arriba="results/fusions/{sample}/arriba_fusions.tsv"
        # FusionCatcher disabled due to installation issues
        # fusioncatcher="results/fusions/{sample}/fusioncatcher_results.txt"
    output:
        integrated="results/fusions/{sample}/integrated_fusions.tsv",
        summary="results/fusions/{sample}/fusion_summary.json"
    params:
        min_support=config.get("fusion_min_support", 3),
        driver_genes=config.get("driver_gene_list", "resources/databases/driver_genes.txt")
    log:
        "logs/fusions/integrate_fusions_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/integrate_fusion_results.py"

# Filter and annotate fusions
rule annotate_fusions:
    input:
        integrated="results/fusions/{sample}/integrated_fusions.tsv"
    output:
        annotated="results/fusions/{sample}/annotated_fusions.tsv",
        driver_fusions="results/fusions/{sample}/driver_fusions.tsv"
    params:
        gene_annotations=config.get("gene_annotations", "resources/databases/gene_annotations.gtf"),
        oncogene_list=config.get("oncogene_list", "resources/databases/oncogenes.txt"),
        tumor_suppressor_list=config.get("tumor_suppressor_list", "resources/databases/tumor_suppressors.txt")
    log:
        "logs/fusions/annotate_fusions_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/annotate_fusions.py"

# Quality control for fusion detection
rule fusion_qc:
    input:
        arriba="results/fusions/{sample}/arriba_fusions.tsv",
        # FusionCatcher disabled due to installation issues
        # fusioncatcher="results/fusions/{sample}/fusioncatcher_results.txt", 
        integrated="results/fusions/{sample}/integrated_fusions.tsv",
        bam_stats="results/alignment/{sample}/{sample}.alignment_stats.json"
    output:
        qc_report="results/fusions/{sample}/fusion_qc.json",
        qc_plot="results/fusions/{sample}/fusion_qc.png"
    log:
        "logs/fusions/fusion_qc_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fusion_qc.py"

# Create fusion summary across all samples
rule fusion_summary_all:
    input:
        summaries=expand("results/fusions/{sample}/fusion_summary.json", sample=SAMPLES)
    output:
        summary_table="results/fusions/all_samples_fusion_summary.tsv",
        summary_plot="results/fusions/all_samples_fusion_plot.png"
    log:
        "logs/fusions/fusion_summary_all.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_all_fusions.py"