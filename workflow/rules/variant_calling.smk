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

# Add read groups first (CRITICAL - missing from new pipeline!)
rule replace_rg:
    input:
        "results/alignment/{sample}/{sample}.sorted.bam"
    output:
        temp("results/variants/{sample}/fixed-rg/{sample}.bam")
    log:
        "logs/variants/replace_rg_{sample}.log"
    params:
        extra="--RGLB lib1 --RGPL illumina --RGPU {sample} --RGSM {sample}",
        java_opts=""
    resources:
        mem_mb=2048,
        runtime=60
    wrapper:
        "v3.10.2/bio/picard/addorreplacereadgroups"

# Mark duplicates
rule mark_duplicates:
    input:
        bams="results/variants/{sample}/fixed-rg/{sample}.bam"
    output:
        bam=temp("results/variants/{sample}/deduped_bam/{sample}.bam"),
        metrics=temp("results/variants/{sample}/deduped_bam/{sample}.metrics.txt")
    log:
        "logs/variants/mark_duplicates_{sample}.log"
    resources:
        mem_mb=5000,
        runtime=120
    wrapper:
        "v3.10.2/bio/picard/markduplicates"

# Index BAM file
rule index_bam:
    input:
        "results/variants/{sample}/deduped_bam/{sample}.bam"
    output:
        temp("results/variants/{sample}/deduped_bam/{sample}.bam.bai")
    log:
        "logs/variants/index_bam_{sample}.log"
    params:
        extra=""
    threads: 4
    wrapper:
        "v3.10.2/bio/samtools/index"

# Split reads that contain Ns in their cigar string
rule splitncigarreads:
    input:
        bam="results/variants/{sample}/deduped_bam/{sample}.bam",
        bai="results/variants/{sample}/deduped_bam/{sample}.bam.bai",
        ref=config["reference_genome"]
    output:
        temp("results/variants/{sample}/split/{sample}.bam")
    log:
        "logs/variants/splitncigarreads_{sample}.log"
    params:
        extra="",
        java_opts=""
    resources:
        mem_mb=5000,
        runtime=180
    threads: 4
    wrapper:
        "v3.10.2/bio/gatk/splitncigarreads"

# Base recalibrator
rule gatk_baserecalibrator:
    input:
        bam="results/variants/{sample}/split/{sample}.bam",
        ref=config["reference_genome"],
        dict=config["reference_genome"].replace(".fa", ".dict").replace(".fasta", ".dict")
        # known sites removed - optional in original pipeline
    output:
        recal_table=temp("results/variants/{sample}/recal/{sample}_recal.table")
    log:
        "logs/variants/baserecalibrator_{sample}.log"
    params:
        extra="",
        java_opts=""
    resources:
        mem_mb=5000,
        runtime=120
    threads: 4
    wrapper:
        "v3.10.2/bio/gatk/baserecalibrator"

# Apply base recalibration
rule gatk_applybqsr:
    input:
        bam="results/variants/{sample}/split/{sample}.bam",
        ref=config["reference_genome"],
        dict=config["reference_genome"].replace(".fa", ".dict").replace(".fasta", ".dict"),
        recal_table="results/variants/{sample}/recal/{sample}_recal.table"
    output:
        bam=temp("results/variants/{sample}/recal/{sample}.bam")
    log:
        "logs/variants/applybqsr_{sample}.log"
    params:
        extra="",
        java_opts="",
        embed_ref=True
    resources:
        mem_mb=5000,
        runtime=120
    threads: 4
    wrapper:
        "v3.10.2/bio/gatk/applybqsr"

# Call variants with HaplotypeCaller
rule haplotype_caller:
    input:
        bam="results/variants/{sample}/recal/{sample}.bam",
        ref=config["reference_genome"]
    output:
        vcf=temp("results/variants/{sample}/calls/{sample}.vcf")
    log:
        "logs/variants/haplotypecaller_{sample}.log"
    params:
        extra="-ERC GVCF --output-mode EMIT_ALL_CONFIDENT_SITES --dont-use-soft-clipped-bases -stand-call-conf 20.0",
        java_opts=""
    threads: 4
    resources:
        mem_mb=10000,
        runtime=240
    wrapper:
        "v3.10.2/bio/gatk/haplotypecaller"

# Filter variants
rule gatk_filter:
    input:
        vcf="results/variants/{sample}/calls/{sample}.vcf",
        ref=config["reference_genome"]
    output:
        vcf="results/variants/{sample}/filtered_variants.vcf"
    log:
        "logs/variants/gatk_filter_{sample}.log"
    params:
        filters={"myfilter": "AB < 0.2 || MQ0 > 50"},
        extra="",
        java_opts=""
    resources:
        mem_mb=5000,
        runtime=60
    threads: 4
    wrapper:
        "v3.10.2/bio/gatk/variantfiltration"

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

# Original pysamstats hotspot detection (CRITICAL for B-ALL classification)
rule pysamstat:
    input:
        bam="results/alignment/{sample}/{sample}.sorted.bam",
        bai="results/alignment/{sample}/{sample}.sorted.bam.bai", 
        fa=config["reference_genome"]
    output:
        pysamstats_output_dir=temp(directory("results/variants/{sample}/pysamstats_output_dir/{sample}/")),
        ikzf1="results/variants/{sample}/pysamstats_output_dir/{sample}/{sample}_IKZF1.csv"
    params:
        out_dir="results/variants/{sample}/pysamstats_output_dir/{sample}/"
    conda:
        "../envs/variants.yaml"
    threads: 4
    resources:
        mem_mb=20000,
        runtime=120
    shell:
        """
        mkdir -p {output.pysamstats_output_dir} &&
        pysamstats --type variation --chromosome 7 -u --start 50382593 --end 50382596 -f {input.fa} {input.bam} > {output.ikzf1} &&
        python workflow/scripts/run_pysamstats_original.py {input.bam} {input.fa} {wildcards.sample} {params.out_dir} 
        """

# Get hotspots using original R script
rule get_hotspots:
    input:
        pysamstats_output_dir="results/variants/{sample}/pysamstats_output_dir/{sample}/",
        gatk_file="results/variants/{sample}/filtered_variants.vcf"
    output:
        hotspot_output_dir=directory("results/variants/{sample}/hotspots")
    shell:
        """
        mkdir -p {output.hotspot_output_dir} &&
        cp /media/nadine/InternalMaybe/Blast-o-Matic-Fusioninator_cluster/scripts/Get_Amino_for_Hotspot.R workflow/scripts/ &&
        Rscript workflow/scripts/Get_Amino_for_Hotspot.R {input.pysamstats_output_dir} {input.gatk_file} {output.hotspot_output_dir}
        """

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