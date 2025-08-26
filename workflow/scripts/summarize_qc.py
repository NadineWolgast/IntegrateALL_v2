#!/usr/bin/env python3
"""
Quality Control Summary Script for IntegrateALL Pipeline
========================================================

Summarizes QC metrics from MultiQC data across all samples.
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import sys

def load_multiqc_data(multiqc_files):
    """Load and combine MultiQC data from multiple samples"""
    combined_data = {}
    
    for file_path in multiqc_files:
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
                
            # Extract sample name from file path
            sample_name = Path(file_path).stem.replace('_multiqc_data', '')
            combined_data[sample_name] = data
            
        except Exception as e:
            print(f"Warning: Could not load {file_path}: {e}")
            continue
    
    return combined_data

def extract_qc_metrics(multiqc_data):
    """Extract key QC metrics from MultiQC data"""
    qc_summary = []
    
    for sample_name, data in multiqc_data.items():
        sample_metrics = {'sample_id': sample_name}
        
        # FastQC metrics
        if 'report_general_stats_data' in data:
            general_stats = data['report_general_stats_data'][0] if data['report_general_stats_data'] else {}
            
            # Extract FastQC metrics
            for key, value in general_stats.items():
                if sample_name in value:
                    sample_data = value[sample_name]
                    if 'fastqc' in key.lower():
                        sample_metrics.update({
                            'total_sequences': sample_data.get('total_sequences', 0),
                            'percent_gc': sample_data.get('percent_gc', 0),
                            'avg_sequence_length': sample_data.get('avg_sequence_length', 0),
                            'percent_duplicates': sample_data.get('percent_duplicates', 0),
                            'percent_fails': sample_data.get('percent_fails', 0)
                        })
        
        # Per base sequence quality
        if 'report_plot_data' in data and 'fastqc_per_base_sequence_quality_plot' in data['report_plot_data']:
            quality_data = data['report_plot_data']['fastqc_per_base_sequence_quality_plot']
            if sample_name in quality_data:
                base_qualities = quality_data[sample_name]['data']
                if base_qualities:
                    avg_quality = np.mean([point[1] for point in base_qualities])
                    sample_metrics['avg_base_quality'] = avg_quality
        
        qc_summary.append(sample_metrics)
    
    return pd.DataFrame(qc_summary)

def create_qc_plots(qc_df, output_file):
    """Create QC summary plots"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Quality Control Summary', fontsize=16)
    
    # Total sequences
    if 'total_sequences' in qc_df.columns:
        axes[0,0].bar(qc_df['sample_id'], qc_df['total_sequences'])
        axes[0,0].set_title('Total Sequences per Sample')
        axes[0,0].set_xlabel('Sample')
        axes[0,0].set_ylabel('Total Sequences')
        axes[0,0].tick_params(axis='x', rotation=45)
    
    # GC content distribution
    if 'percent_gc' in qc_df.columns:
        axes[0,1].hist(qc_df['percent_gc'].dropna(), bins=20, alpha=0.7)
        axes[0,1].set_title('GC Content Distribution')
        axes[0,1].set_xlabel('GC Content (%)')
        axes[0,1].set_ylabel('Number of Samples')
    
    # Average base quality
    if 'avg_base_quality' in qc_df.columns:
        axes[1,0].bar(qc_df['sample_id'], qc_df['avg_base_quality'])
        axes[1,0].set_title('Average Base Quality')
        axes[1,0].set_xlabel('Sample')
        axes[1,0].set_ylabel('Average Quality Score')
        axes[1,0].tick_params(axis='x', rotation=45)
    
    # Duplication rate
    if 'percent_duplicates' in qc_df.columns:
        axes[1,1].bar(qc_df['sample_id'], qc_df['percent_duplicates'])
        axes[1,1].set_title('Duplication Rate')
        axes[1,1].set_xlabel('Sample')
        axes[1,1].set_ylabel('Duplication Rate (%)')
        axes[1,1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Get inputs from snakemake
    multiqc_files = snakemake.input.multiqc_data
    output_summary = snakemake.output.summary
    output_plot = snakemake.output.plot
    
    # Load MultiQC data
    multiqc_data = load_multiqc_data(multiqc_files)
    
    if not multiqc_data:
        print("No valid MultiQC data found")
        # Create empty outputs
        with open(output_summary, 'w') as f:
            json.dump({}, f)
        
        # Create empty plot
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, 'No QC data available', ha='center', va='center')
        plt.savefig(output_plot)
        plt.close()
        return
    
    # Extract QC metrics
    qc_df = extract_qc_metrics(multiqc_data)
    
    # Create summary statistics
    summary_stats = {
        'total_samples': len(qc_df),
        'metrics': {}
    }
    
    # Calculate summary statistics for each metric
    numeric_columns = qc_df.select_dtypes(include=[np.number]).columns
    for col in numeric_columns:
        if col != 'sample_id':
            summary_stats['metrics'][col] = {
                'mean': float(qc_df[col].mean()),
                'median': float(qc_df[col].median()),
                'std': float(qc_df[col].std()),
                'min': float(qc_df[col].min()),
                'max': float(qc_df[col].max())
            }
    
    # Save summary
    with open(output_summary, 'w') as f:
        json.dump(summary_stats, f, indent=2)
    
    # Create plots
    create_qc_plots(qc_df, output_plot)
    
    print(f"QC summary completed for {len(qc_df)} samples")

if __name__ == "__main__":
    main()