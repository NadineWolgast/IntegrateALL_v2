"""
Alignment Rules for IntegrateALL Pipeline
=========================================

STAR alignment with optimized parameters for fusion detection and RNA-seq analysis
"""

rule star_index:
    input:
        fasta=config["reference_genome"],
        gtf=config["reference_gtf"]
    output:
        directory("resources/star_index")
    params:
        extra_args=config.get("star_index_args", ""),
        sjdb_overhang=config.get("sjdb_overhang", 100),
        sa_index_nbases=config.get("sa_index_nbases", 14)
    threads:
        config.get("star_index_threads", 16)
    resources:
        mem_mb=config.get("star_index_memory", 32000),
        runtime=config.get("star_index_runtime", 120)
    log:
        "logs/alignment/star_index.log"
    benchmark:
        "benchmarks/alignment/star_index.benchmark.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        mkdir -p {output}
        
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN {threads} \\
            --genomeDir {output} \\
            --genomeFastaFiles {input.fasta} \\
            --sjdbGTFfile {input.gtf} \\
            --sjdbOverhang {params.sjdb_overhang} \\
            --genomeSAindexNbases {params.sa_index_nbases} \\
            {params.extra_args} \\
            &> {log}
        """

rule star_align:
    input:
        fastq1=lambda wildcards: samples_df[samples_df['sample_id'] == wildcards.sample]['fastq1'].iloc[0],
        fastq2=lambda wildcards: samples_df[samples_df['sample_id'] == wildcards.sample]['fastq2'].iloc[0],
        index="resources/star_index"
    output:
        bam="results/alignment/{sample}/{sample}.sorted.bam",
        bai="results/alignment/{sample}/{sample}.sorted.bam.bai",
        gene_counts="results/alignment/{sample}/{sample}.gene_counts.tsv",
        log_final="results/alignment/{sample}/{sample}.Log.final.out",
        log_out="results/alignment/{sample}/{sample}.Log.out",
        log_progress="results/alignment/{sample}/{sample}.Log.progress.out",
        sj_tab="results/alignment/{sample}/{sample}.SJ.out.tab"
    params:
        output_prefix="results/alignment/{sample}/{sample}.",
        extra_args=config.get("star_align_args", ""),
        read_files_command=config.get("star_read_files_command", "zcat"),
        two_pass_mode=config.get("star_two_pass_mode", "Basic"),
        chim_segment_min=config.get("star_chim_segment_min", 15),
        chim_junction_overhang_min=config.get("star_chim_junction_overhang_min", 15),
        alignSJDBoverhangMin=config.get("star_alignSJDBoverhangMin", 3),
        outFilterMultimapNmax=config.get("star_outFilterMultimapNmax", 20),
        outFilterMismatchNmax=config.get("star_outFilterMismatchNmax", 999),
        outFilterMismatchNoverReadLmax=config.get("star_outFilterMismatchNoverReadLmax", 0.04),
        alignIntronMin=config.get("star_alignIntronMin", 20),
        alignIntronMax=config.get("star_alignIntronMax", 1000000),
        alignMatesGapMax=config.get("star_alignMatesGapMax", 1000000)
    threads:
        config.get("star_align_threads", 16)
    resources:
        mem_mb=config.get("star_align_memory", 32000),
        runtime=config.get("star_align_runtime", 240)
    log:
        "logs/alignment/star_align_{sample}.log"
    benchmark:
        "benchmarks/alignment/star_align_{sample}.benchmark.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {params.output_prefix})
        
        # Run STAR alignment
        STAR \\
            --runMode alignReads \\
            --runThreadN {threads} \\
            --genomeDir {input.index} \\
            --readFilesIn {input.fastq1} {input.fastq2} \\
            --readFilesCommand {params.read_files_command} \\
            --outFileNamePrefix {params.output_prefix} \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMunmapped Within \\
            --outSAMattributes Standard \\
            --outSAMstrandField intronMotif \\
            --quantMode GeneCounts \\
            --twopassMode {params.two_pass_mode} \\
            --chimSegmentMin {params.chim_segment_min} \\
            --chimJunctionOverhangMin {params.chim_junction_overhang_min} \\
            --chimOutType Junctions WithinBAM SeparateSAMold \\
            --chimMainSegmentMultNmax 1 \\
            --outFilterMultimapNmax {params.outFilterMultimapNmax} \\
            --outFilterMismatchNmax {params.outFilterMismatchNmax} \\
            --outFilterMismatchNoverReadLmax {params.outFilterMismatchNoverReadLmax} \\
            --alignIntronMin {params.alignIntronMin} \\
            --alignIntronMax {params.alignIntronMax} \\
            --alignMatesGapMax {params.alignMatesGapMax} \\
            --alignSJDBoverhangMin {params.alignSJDBoverhangMin} \\
            {params.extra_args} \\
            &> {log}
        
        # Rename output files to match expected names
        mv {params.output_prefix}Aligned.sortedByCoord.out.bam {output.bam}
        mv {params.output_prefix}ReadsPerGene.out.tab {output.gene_counts}
        
        # Index BAM file
        samtools index {output.bam} {output.bai}
        """

rule flagstat:
    input:
        bam="results/alignment/{sample}/{sample}.sorted.bam"
    output:
        flagstat="results/alignment/{sample}/{sample}.flagstat"
    resources:
        mem_mb=config.get("samtools_memory", 4000),
        runtime=30
    log:
        "logs/alignment/flagstat_{sample}.log"
    benchmark:
        "benchmarks/alignment/flagstat_{sample}.benchmark.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        """

rule idxstats:
    input:
        bam="results/alignment/{sample}/{sample}.sorted.bam",
        bai="results/alignment/{sample}/{sample}.sorted.bam.bai"
    output:
        idxstats="results/alignment/{sample}/{sample}.idxstats"
    resources:
        mem_mb=config.get("samtools_memory", 4000),
        runtime=30
    log:
        "logs/alignment/idxstats_{sample}.log"
    benchmark:
        "benchmarks/alignment/idxstats_{sample}.benchmark.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats} 2> {log}
        """

rule alignment_stats:
    input:
        bam="results/alignment/{sample}/{sample}.sorted.bam",
        flagstat="results/alignment/{sample}/{sample}.flagstat",
        idxstats="results/alignment/{sample}/{sample}.idxstats",
        star_log="results/alignment/{sample}/{sample}.Log.final.out"
    output:
        stats="results/alignment/{sample}/{sample}.alignment_stats.json"
    log:
        "logs/alignment/alignment_stats_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/calculate_alignment_stats.py"

rule prepare_gene_counts:
    input:
        gene_counts="results/alignment/{sample}/{sample}.gene_counts.tsv"
    output:
        counts_formatted="results/alignment/{sample}/{sample}.counts_formatted.tsv"
    params:
        sample_id="{sample}"
    log:
        "logs/alignment/prepare_gene_counts_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/format_gene_counts.py"

# Rule to merge gene counts from all samples for downstream analysis
rule merge_gene_counts:
    input:
        counts=expand("results/alignment/{sample}/{sample}.counts_formatted.tsv", sample=SAMPLES)
    output:
        merged_counts="results/alignment/merged_gene_counts.tsv"
    log:
        "logs/alignment/merge_gene_counts.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_gene_counts.py"