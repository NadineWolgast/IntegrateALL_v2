#!/usr/bin/env python3
"""
Calculate Alignment Statistics Script for IntegrateALL Pipeline
===============================================================

Calculates comprehensive alignment statistics from BAM files and STAR logs.
"""

import json
import re
import pandas as pd
import pysam
from pathlib import Path

def parse_star_log(star_log_file):
    """Parse STAR alignment log for key statistics"""
    stats = {}
    
    try:
        with open(star_log_file, 'r') as f:
            content = f.read()
        
        # Extract key statistics using regex
        patterns = {
            'input_reads': r'Number of input reads \|\s+(\d+)',
            'avg_input_read_length': r'Average input read length \|\s+([\d.]+)',
            'uniquely_mapped': r'Uniquely mapped reads number \|\s+(\d+)',
            'uniquely_mapped_percent': r'Uniquely mapped reads % \|\s+([\d.]+)%',
            'avg_mapped_length': r'Average mapped length \|\s+([\d.]+)',
            'splices_total': r'Number of splices: Total \|\s+(\d+)',
            'splices_annotated': r'Number of splices: Annotated \(sjdb\) \|\s+(\d+)',
            'splices_gc_ag': r'Number of splices: GC/AG \|\s+(\d+)',
            'splices_at_ac': r'Number of splices: AT/AC \|\s+(\d+)',
            'splices_non_canonical': r'Number of splices: Non-canonical \|\s+(\d+)',
            'mismatch_rate': r'Mismatch rate per base, % \|\s+([\d.]+)%',
            'deletion_rate': r'Deletion rate per base \|\s+([\d.]+)%',
            'deletion_avg_length': r'Deletion average length \|\s+([\d.]+)',
            'insertion_rate': r'Insertion rate per base \|\s+([\d.]+)%',
            'insertion_avg_length': r'Insertion average length \|\s+([\d.]+)',
            'multimapping': r'Number of reads mapped to multiple loci \|\s+(\d+)',
            'multimapping_percent': r'% of reads mapped to multiple loci \|\s+([\d.]+)%',
            'too_many_loci': r'Number of reads mapped to too many loci \|\s+(\d+)',
            'too_many_loci_percent': r'% of reads mapped to too many loci \|\s+([\d.]+)%',
            'unmapped_mismatches': r'% of reads unmapped: too many mismatches \|\s+([\d.]+)%',
            'unmapped_short': r'% of reads unmapped: too short \|\s+([\d.]+)%',
            'unmapped_other': r'% of reads unmapped: other \|\s+([\d.]+)%',
            'chimeric_reads': r'Number of chimeric reads \|\s+(\d+)',
            'chimeric_percent': r'% of chimeric reads \|\s+([\d.]+)%'
        }
        
        for key, pattern in patterns.items():
            match = re.search(pattern, content)
            if match:
                value = match.group(1)
                try:
                    # Convert to float if contains decimal, otherwise int
                    if '.' in value:
                        stats[key] = float(value)
                    else:
                        stats[key] = int(value)
                except ValueError:
                    stats[key] = value
        
        # Calculate derived statistics
        if 'uniquely_mapped' in stats and 'input_reads' in stats:
            stats['mapping_rate'] = (stats['uniquely_mapped'] / stats['input_reads']) * 100
        
        if 'splices_annotated' in stats and 'splices_total' in stats and stats['splices_total'] > 0:
            stats['annotated_splice_rate'] = (stats['splices_annotated'] / stats['splices_total']) * 100
            
    except FileNotFoundError:
        print(f"STAR log file not found: {star_log_file}")
    except Exception as e:
        print(f"Error parsing STAR log {star_log_file}: {e}")
    
    return stats

def parse_flagstat(flagstat_file):
    """Parse samtools flagstat output"""
    stats = {}
    
    try:
        with open(flagstat_file, 'r') as f:
            content = f.read()
        
        lines = content.strip().split('\n')
        for line in lines:
            # Extract total reads
            if 'in total' in line:
                total_reads = int(line.split()[0])
                stats['total_reads_flagstat'] = total_reads
            
            # Extract mapped reads  
            elif 'mapped (' in line and 'mate mapped' not in line:
                mapped_reads = int(line.split()[0])
                percent_match = re.search(r'\(([\d.]+)%', line)
                stats['mapped_reads_flagstat'] = mapped_reads
                if percent_match:
                    stats['mapped_percent_flagstat'] = float(percent_match.group(1))
            
            # Extract properly paired
            elif 'properly paired' in line:
                properly_paired = int(line.split()[0])
                percent_match = re.search(r'\(([\d.]+)%', line)
                stats['properly_paired'] = properly_paired
                if percent_match:
                    stats['properly_paired_percent'] = float(percent_match.group(1))
            
            # Extract singletons
            elif 'singletons' in line:
                singletons = int(line.split()[0])
                percent_match = re.search(r'\(([\d.]+)%', line)
                stats['singletons'] = singletons
                if percent_match:
                    stats['singletons_percent'] = float(percent_match.group(1))
                    
    except FileNotFoundError:
        print(f"Flagstat file not found: {flagstat_file}")
    except Exception as e:
        print(f"Error parsing flagstat {flagstat_file}: {e}")
    
    return stats

def parse_idxstats(idxstats_file):
    """Parse samtools idxstats output"""
    stats = {}
    
    try:
        df = pd.read_csv(idxstats_file, sep='\t', header=None, 
                        names=['chromosome', 'length', 'mapped', 'unmapped'])
        
        # Calculate chromosome-level statistics
        stats['chromosomes_with_reads'] = len(df[df['mapped'] > 0])
        stats['total_chromosomes'] = len(df[df['chromosome'] != '*'])
        stats['reads_unmapped_idxstats'] = df[df['chromosome'] == '*']['unmapped'].sum() if '*' in df['chromosome'].values else 0
        
        # Top 10 chromosomes by read count
        top_chromosomes = df[df['chromosome'] != '*'].nlargest(10, 'mapped')
        stats['top_chromosomes'] = top_chromosomes[['chromosome', 'mapped']].to_dict('records')
        
        # Mitochondrial reads (if present)
        mito_chromosomes = ['chrM', 'MT', 'chrMT']
        mito_data = df[df['chromosome'].isin(mito_chromosomes)]
        if not mito_data.empty:
            stats['mitochondrial_reads'] = mito_data['mapped'].sum()
            total_mapped = df['mapped'].sum()
            if total_mapped > 0:
                stats['mitochondrial_percent'] = (stats['mitochondrial_reads'] / total_mapped) * 100
        
    except FileNotFoundError:
        print(f"Idxstats file not found: {idxstats_file}")
    except Exception as e:
        print(f"Error parsing idxstats {idxstats_file}: {e}")
    
    return stats

def calculate_bam_stats(bam_file):
    """Calculate additional statistics from BAM file using pysam"""
    stats = {}
    
    try:
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            # Sample a subset of reads for quality statistics (first 100,000 reads)
            read_count = 0
            quality_scores = []
            insert_sizes = []
            
            for read in bam.fetch():
                if read_count >= 100000:  # Sample limit to avoid memory issues
                    break
                
                if not read.is_unmapped and read.is_proper_pair:
                    # Mapping quality
                    quality_scores.append(read.mapping_quality)
                    
                    # Insert size (template length)
                    if read.template_length > 0:
                        insert_sizes.append(abs(read.template_length))
                
                read_count += 1
            
            # Calculate statistics
            if quality_scores:
                stats['avg_mapping_quality'] = sum(quality_scores) / len(quality_scores)
                stats['min_mapping_quality'] = min(quality_scores)
                stats['max_mapping_quality'] = max(quality_scores)
            
            if insert_sizes:
                stats['avg_insert_size'] = sum(insert_sizes) / len(insert_sizes)
                stats['median_insert_size'] = sorted(insert_sizes)[len(insert_sizes)//2]
                stats['min_insert_size'] = min(insert_sizes)
                stats['max_insert_size'] = max(insert_sizes)
                
    except Exception as e:
        print(f"Error analyzing BAM file {bam_file}: {e}")
    
    return stats

def calculate_quality_flags(stats):
    """Calculate quality flags based on alignment statistics"""
    flags = {}
    
    # Mapping rate flag
    mapping_rate = stats.get('mapping_rate', 0)
    if mapping_rate >= 80:
        flags['mapping_rate_flag'] = 'PASS'
    elif mapping_rate >= 60:
        flags['mapping_rate_flag'] = 'WARN'
    else:
        flags['mapping_rate_flag'] = 'FAIL'
    
    # Properly paired flag
    properly_paired_percent = stats.get('properly_paired_percent', 0)
    if properly_paired_percent >= 80:
        flags['properly_paired_flag'] = 'PASS'
    elif properly_paired_percent >= 60:
        flags['properly_paired_flag'] = 'WARN'
    else:
        flags['properly_paired_flag'] = 'FAIL'
    
    # Duplication rate flag (if available)
    if 'singletons_percent' in stats:
        if stats['singletons_percent'] <= 5:
            flags['singletons_flag'] = 'PASS'
        elif stats['singletons_percent'] <= 10:
            flags['singletons_flag'] = 'WARN'
        else:
            flags['singletons_flag'] = 'FAIL'
    
    # Overall alignment quality
    pass_count = sum(1 for flag in flags.values() if flag == 'PASS')
    total_count = len(flags)
    
    if total_count > 0:
        overall_score = pass_count / total_count
        if overall_score >= 0.8:
            flags['overall_alignment_quality'] = 'PASS'
        elif overall_score >= 0.6:
            flags['overall_alignment_quality'] = 'WARN'
        else:
            flags['overall_alignment_quality'] = 'FAIL'
    
    return flags

def main():
    # Get inputs from snakemake
    bam_file = snakemake.input.bam
    flagstat_file = snakemake.input.flagstat
    idxstats_file = snakemake.input.idxstats
    star_log_file = snakemake.input.star_log
    output_stats = snakemake.output.stats
    
    # Extract sample name
    sample_id = Path(bam_file).stem.replace('.sorted', '')
    
    # Collect all statistics
    all_stats = {'sample_id': sample_id}
    
    # Parse STAR log
    star_stats = parse_star_log(star_log_file)
    all_stats.update(star_stats)
    
    # Parse flagstat
    flagstat_stats = parse_flagstat(flagstat_file)
    all_stats.update(flagstat_stats)
    
    # Parse idxstats
    idxstats_stats = parse_idxstats(idxstats_file)
    all_stats.update(idxstats_stats)
    
    # Calculate BAM statistics
    bam_stats = calculate_bam_stats(bam_file)
    all_stats.update(bam_stats)
    
    # Calculate quality flags
    quality_flags = calculate_quality_flags(all_stats)
    all_stats.update(quality_flags)
    
    # Add file paths for reference
    all_stats['input_files'] = {
        'bam_file': str(bam_file),
        'flagstat_file': str(flagstat_file),
        'idxstats_file': str(idxstats_file),
        'star_log_file': str(star_log_file)
    }
    
    # Save statistics
    with open(output_stats, 'w') as f:
        json.dump(all_stats, f, indent=2)
    
    print(f"Calculated alignment statistics for sample {sample_id}")
    print(f"Mapping rate: {all_stats.get('mapping_rate', 'N/A'):.2f}%")
    print(f"Overall quality: {all_stats.get('overall_alignment_quality', 'N/A')}")

if __name__ == "__main__":
    main()