"""
Quality Control Rules for IntegrateALL Pipeline
===============================================

FastQC and MultiQC analysis with optimized performance
"""

rule fastqc:
    input:
        lambda wildcards: get_fastq_files(wildcards)[f'fastq{wildcards.read}']
    output:
        html="results/qc/fastqc/{sample}_R{read}_fastqc.html",
        zip="results/qc/fastqc/{sample}_R{read}_fastqc.zip"
    params:
        output_dir="results/qc/fastqc",
        extra_args=config.get("fastqc_args", "")
    threads: 
        config.get("threads", 4)
    resources:
        mem_mb=config.get("fastqc_memory", 2000),
        runtime=config.get("fastqc_runtime", 30)
    log:
        "logs/qc/fastqc/{sample}_R{read}.log"
    benchmark:
        "benchmarks/qc/fastqc_{sample}_R{read}.benchmark.txt"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}
        
        # Run FastQC
        fastqc \\
            --threads {threads} \\
            --outdir {params.output_dir} \\
            --extract \\
            {params.extra_args} \\
            "{input}" \\
            &> {log}
        
        # Rename output files to match expected pattern
        BASENAME=$(basename "{input}" | sed 's/\\.gz$//' | sed 's/\\.fq$//' | sed 's/\\.fastq$//')
        mv {params.output_dir}/${{BASENAME}}_fastqc.html {output.html}
        mv {params.output_dir}/${{BASENAME}}_fastqc.zip {output.zip}
        """

rule multiqc_sample:
    input:
        fastqc_r1="results/qc/fastqc/{sample}_R1_fastqc.zip",
        fastqc_r2="results/qc/fastqc/{sample}_R2_fastqc.zip"
    output:
        html="results/qc/multiqc/{sample}_multiqc_report.html",
        data="results/qc/multiqc/{sample}_multiqc_data.json"
    params:
        input_dir="results/qc/fastqc/",
        output_dir="results/qc/multiqc",
        config_file=config.get("multiqc_config", ""),
        extra_args=config.get("multiqc_args", "--force")
    log:
        "logs/qc/multiqc/{sample}.log"
    benchmark:
        "benchmarks/qc/multiqc_{sample}.benchmark.txt"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}
        
        # Create temporary directory for this sample's files
        TEMP_DIR=$(mktemp -d)
        
        # Copy FastQC files for this sample only
        cp {input.fastqc_r1} {input.fastqc_r2} "$TEMP_DIR/"
        
        # Run MultiQC
        multiqc \\
            "$TEMP_DIR" \\
            --outdir {params.output_dir} \\
            --filename {wildcards.sample}_multiqc_report \\
            --data-format json \\
            --export \\
            {params.extra_args} \\
            &> {log}
        
        # Move data file to expected location
        if [ -f "{params.output_dir}/multiqc_data.json" ]; then
            mv {params.output_dir}/multiqc_data.json {output.data}
        else
            echo '{{}}' > {output.data}
        fi
        
        # Cleanup
        rm -rf "$TEMP_DIR"
        """

rule multiqc_all:
    input:
        expand("results/qc/fastqc/{sample}_R{read}_fastqc.zip", 
               sample=SAMPLES, read=[1, 2])
    output:
        html="results/qc/multiqc_all.html",
        data_dir=directory("results/qc/multiqc_data")
    params:
        input_dir="results/qc/fastqc/",
        output_dir="results/qc/",
        config_file=config.get("multiqc_config", ""),
        extra_args=config.get("multiqc_args", "--force")
    log:
        "logs/qc/multiqc_all.log"
    benchmark:
        "benchmarks/qc/multiqc_all.benchmark.txt"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc \\
            {params.input_dir} \\
            --outdir {params.output_dir} \\
            --filename multiqc_all \\
            {params.extra_args} \\
            &> {log}
        """

rule qc_summary:
    input:
        multiqc_data=expand("results/qc/multiqc/{sample}_multiqc_data.json", sample=SAMPLES)
    output:
        summary="results/qc/qc_summary.json",
        plot="results/qc/qc_summary.png"
    log:
        "logs/qc/qc_summary.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_qc.py"

# Rule to extract key QC metrics for downstream analysis
rule extract_qc_metrics:
    input:
        fastqc_r1="results/qc/fastqc/{sample}_R1_fastqc.zip",
        fastqc_r2="results/qc/fastqc/{sample}_R2_fastqc.zip"
    output:
        metrics="results/qc/metrics/{sample}_qc_metrics.json"
    log:
        "logs/qc/extract_metrics_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_qc_metrics.py"