"""
Variant Calling Rules for IntegrateALL Pipeline
===============================================

GATK best practices for RNA-seq variant calling
"""

# Prepare reference genome for GATK
rule prepare_reference:
    input:
        fasta=config["reference_genome"]
    output:
        dict=config["reference_genome"].replace(".fa", ".dict").replace(".fasta", ".dict"),
        fai=config["reference_genome"] + ".fai"
    resources:
        mem_mb=4000,
        runtime=60
    log:
        "logs/variants/prepare_reference.log"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        # Create FASTA index
        samtools faidx {input.fasta}
        
        # Create sequence dictionary
        picard CreateSequenceDictionary \\
            R={input.fasta} \\
            O={output.dict} \\
            &> {log}
        """

# Mark duplicates
rule mark_duplicates:
    input:
        bam="results/alignment/{sample}/{sample}.sorted.bam"
    output:
        bam="results/variants/{sample}/marked_duplicates.bam",
        metrics="results/variants/{sample}/duplicate_metrics.txt"
    resources:
        mem_mb=config.get("picard_memory", 16000),
        runtime=config.get("picard_runtime", 120)
    log:
        "logs/variants/mark_duplicates_{sample}.log"
    benchmark:
        "benchmarks/variants/mark_duplicates_{sample}.benchmark.txt"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        
        picard MarkDuplicates \\
            I={input.bam} \\
            O={output.bam} \\
            M={output.metrics} \\
            VALIDATION_STRINGENCY=LENIENT \\
            ASSUME_SORT_ORDER=coordinate \\
            &> {log}
        
        # Index the output BAM
        samtools index {output.bam}
        """

# Split reads that contain Ns in their cigar string
rule split_n_cigar:
    input:
        bam="results/variants/{sample}/marked_duplicates.bam",
        ref=config["reference_genome"],
        dict=config["reference_genome"].replace(".fa", ".dict").replace(".fasta", ".dict")
    output:
        bam="results/variants/{sample}/split_reads.bam"
    resources:
        mem_mb=config.get("gatk_memory", 16000),
        runtime=config.get("gatk_runtime", 180)
    log:
        "logs/variants/split_n_cigar_{sample}.log"
    benchmark:
        "benchmarks/variants/split_n_cigar_{sample}.benchmark.txt"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        gatk SplitNCigarReads \\
            -R {input.ref} \\
            -I {input.bam} \\
            -O {output.bam} \\
            &> {log}
        """

# Call variants with HaplotypeCaller
rule call_variants:
    input:
        bam="results/variants/{sample}/split_reads.bam",
        ref=config["reference_genome"],
        dict=config["reference_genome"].replace(".fa", ".dict").replace(".fasta", ".dict")
    output:
        vcf="results/variants/{sample}/raw_variants.vcf"
    params:
        extra_args=config.get("haplotype_caller_args", ""),
        dont_use_soft_clipped_bases=config.get("gatk_dont_use_soft_clipped", True),
        standard_min_confidence_threshold=config.get("gatk_min_confidence", 20.0)
    resources:
        mem_mb=config.get("gatk_memory", 16000),
        runtime=config.get("gatk_runtime", 240)
    log:
        "logs/variants/call_variants_{sample}.log"
    benchmark:
        "benchmarks/variants/call_variants_{sample}.benchmark.txt"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        gatk HaplotypeCaller \\
            -R {input.ref} \\
            -I {input.bam} \\
            -O {output.vcf} \\
            --dont-use-soft-clipped-bases {params.dont_use_soft_clipped_bases} \\
            --standard-min-confidence-threshold-for-calling {params.standard_min_confidence_threshold} \\
            {params.extra_args} \\
            &> {log}
        """

# Filter variants
rule filter_variants:
    input:
        vcf="results/variants/{sample}/raw_variants.vcf",
        ref=config["reference_genome"]
    output:
        vcf="results/variants/{sample}/filtered_variants.vcf"
    params:
        filter_expression=config.get("variant_filter_expression", "QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"),
        filter_name=config.get("variant_filter_name", "basic_filter")
    resources:
        mem_mb=config.get("gatk_memory", 8000),
        runtime=60
    log:
        "logs/variants/filter_variants_{sample}.log"
    benchmark:
        "benchmarks/variants/filter_variants_{sample}.benchmark.txt"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        gatk VariantFiltration \\
            -R {input.ref} \\
            -V {input.vcf} \\
            -O {output.vcf} \\
            --filter-expression "{params.filter_expression}" \\
            --filter-name {params.filter_name} \\
            &> {log}
        """

# Annotate variants with VEP or SnpEff
rule annotate_variants:
    input:
        vcf="results/variants/{sample}/filtered_variants.vcf"
    output:
        vcf="results/variants/{sample}/annotated_variants.vcf",
        stats="results/variants/{sample}/annotation_stats.html"
    params:
        cache_dir=config.get("vep_cache", "resources/databases/vep_cache"),
        species=config.get("vep_species", "homo_sapiens"),
        assembly=config.get("vep_assembly", "GRCh38"),
        extra_args=config.get("vep_args", "")
    resources:
        mem_mb=config.get("vep_memory", 8000),
        runtime=config.get("vep_runtime", 120)
    log:
        "logs/variants/annotate_variants_{sample}.log"
    benchmark:
        "benchmarks/variants/annotate_variants_{sample}.benchmark.txt"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        vep \\
            --input_file {input.vcf} \\
            --output_file {output.vcf} \\
            --stats_file {output.stats} \\
            --cache \\
            --dir_cache {params.cache_dir} \\
            --species {params.species} \\
            --assembly {params.assembly} \\
            --format vcf \\
            --vcf \\
            --symbol \\
            --gene_phenotype \\
            --regulatory \\
            --biotype \\
            --canonical \\
            --protein \\
            --domains \\
            --numbers \\
            --total_length \\
            --allele_number \\
            {params.extra_args} \\
            &> {log}
        """

# Hotspot analysis for known B-ALL mutations
rule hotspot_analysis:
    input:
        vcf="results/variants/{sample}/annotated_variants.vcf"
    output:
        hotspots="results/variants/{sample}/hotspots.csv",
        summary="results/variants/{sample}/hotspot_summary.json"
    params:
        hotspot_bed=config.get("hotspot_regions", "resources/databases/ball_hotspots.bed"),
        driver_genes=config.get("driver_gene_list", "resources/databases/driver_genes.txt")
    log:
        "logs/variants/hotspot_analysis_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/analyze_hotspots.py"

# Variant quality control
rule variant_qc:
    input:
        raw_vcf="results/variants/{sample}/raw_variants.vcf",
        filtered_vcf="results/variants/{sample}/filtered_variants.vcf",
        annotated_vcf="results/variants/{sample}/annotated_variants.vcf"
    output:
        qc_report="results/variants/{sample}/variant_qc.json",
        qc_plot="results/variants/{sample}/variant_qc.png"
    log:
        "logs/variants/variant_qc_{sample}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/variant_qc.py"

# Create variant summary across all samples
rule variant_summary_all:
    input:
        hotspots=expand("results/variants/{sample}/hotspots.csv", sample=SAMPLES),
        summaries=expand("results/variants/{sample}/hotspot_summary.json", sample=SAMPLES)
    output:
        summary_table="results/variants/all_samples_variants.tsv",
        recurrent_variants="results/variants/recurrent_variants.tsv",
        variant_plot="results/variants/variant_landscape.png"
    log:
        "logs/variants/variant_summary_all.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_all_variants.py"