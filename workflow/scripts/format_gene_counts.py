#!/usr/bin/env python3
"""
Format Gene Counts Script for IntegrateALL Pipeline
===================================================

Formats STAR gene count output for downstream analysis.
"""

import pandas as pd
import sys
from pathlib import Path

def format_star_counts(count_file, sample_id):
    """
    Format STAR ReadsPerGene.out.tab file
    
    STAR output format:
    Column 1: gene ID
    Column 2: counts for unstranded RNA-seq
    Column 3: counts for the 1st read strand aligned with RNA
    Column 4: counts for the 2nd read strand aligned with RNA
    """
    
    try:
        # Read STAR count file
        df = pd.read_csv(count_file, sep='\t', header=None, 
                        names=['gene_id', 'unstranded', 'strand1', 'strand2'])
        
        # Remove the first 4 rows which contain summary statistics
        # (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
        summary_rows = df.head(4).copy()
        df = df.iloc[4:].copy()
        
        # Create formatted output with sample ID as column name
        formatted_df = pd.DataFrame({
            'gene_id': df['gene_id'],
            sample_id: df['unstranded']  # Use unstranded counts by default
        })
        
        # Convert counts to integers
        formatted_df[sample_id] = formatted_df[sample_id].astype(int)
        
        # Calculate some basic statistics
        total_reads = formatted_df[sample_id].sum()
        genes_with_counts = (formatted_df[sample_id] > 0).sum()
        
        print(f"Formatted gene counts for {sample_id}:")
        print(f"  Total reads assigned to genes: {total_reads:,}")
        print(f"  Genes with non-zero counts: {genes_with_counts:,}")
        print(f"  Total genes: {len(formatted_df):,}")
        
        # Extract summary statistics for reporting
        summary_stats = {
            'sample_id': sample_id,
            'total_assigned_reads': int(total_reads),
            'genes_with_counts': int(genes_with_counts),
            'total_genes': len(formatted_df)
        }
        
        # Add STAR summary statistics if available
        if len(summary_rows) >= 4:
            summary_stats.update({
                'unmapped_reads': int(summary_rows.iloc[0]['unstranded']),
                'multimapping_reads': int(summary_rows.iloc[1]['unstranded']),
                'no_feature_reads': int(summary_rows.iloc[2]['unstranded']),
                'ambiguous_reads': int(summary_rows.iloc[3]['unstranded'])
            })
        
        return formatted_df, summary_stats
        
    except Exception as e:
        print(f"Error formatting gene counts from {count_file}: {e}")
        sys.exit(1)

def main():
    # Get inputs from snakemake
    count_file = snakemake.input.gene_counts
    output_file = snakemake.output.counts_formatted
    sample_id = snakemake.params.sample_id
    
    # Format the counts
    formatted_df, summary_stats = format_star_counts(count_file, sample_id)
    
    # Save formatted counts
    formatted_df.to_csv(output_file, sep='\t', index=False)
    
    # Save summary statistics to a separate file (optional)
    stats_file = output_file.replace('.tsv', '_stats.json')
    import json
    with open(stats_file, 'w') as f:
        json.dump(summary_stats, f, indent=2)
    
    print(f"Formatted gene counts saved to {output_file}")

if __name__ == "__main__":
    main()