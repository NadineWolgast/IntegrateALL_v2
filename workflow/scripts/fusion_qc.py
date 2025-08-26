#!/usr/bin/env python3
"""
Fusion QC Script for IntegrateALL Pipeline
==========================================

Performs quality control analysis on fusion detection results.
"""

import pandas as pd
import matplotlib.pyplot as plt
import json
import numpy as np

def analyze_fusion_quality(arriba_file, fusioncatcher_file, integrated_file):
    """Analyze fusion detection quality metrics"""
    
    metrics = {
        'arriba_fusions': 0,
        'fusioncatcher_fusions': 0,
        'integrated_fusions': 0,
        'consensus_fusions': 0,
        'high_confidence_fusions': 0
    }
    
    try:
        # Load Arriba results
        arriba_df = pd.read_csv(arriba_file, sep='\t')
        metrics['arriba_fusions'] = len(arriba_df)
        
        # Load FusionCatcher results  
        fc_df = pd.read_csv(fusioncatcher_file, sep='\t')
        metrics['fusioncatcher_fusions'] = len(fc_df)
        
        # Load integrated results
        int_df = pd.read_csv(integrated_file, sep='\t')
        metrics['integrated_fusions'] = len(int_df)
        
        if 'callers' in int_df.columns:
            metrics['consensus_fusions'] = len(int_df[int_df['callers'].str.contains(',')])
        
        if 'confidence_score' in int_df.columns:
            metrics['high_confidence_fusions'] = len(int_df[int_df['confidence_score'] >= 0.7])
            
    except Exception as e:
        print(f"Warning: Error analyzing fusion files: {e}")
    
    return metrics

def create_qc_plot(metrics, output_plot):
    """Create fusion QC visualization"""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Fusion caller comparison
    callers = ['Arriba', 'FusionCatcher', 'Integrated']
    counts = [metrics['arriba_fusions'], metrics['fusioncatcher_fusions'], metrics['integrated_fusions']]
    
    ax1.bar(callers, counts)
    ax1.set_title('Fusion Counts by Caller')
    ax1.set_ylabel('Number of Fusions')
    
    # Quality distribution
    quality_labels = ['Total', 'Consensus', 'High Confidence']
    quality_counts = [metrics['integrated_fusions'], metrics['consensus_fusions'], metrics['high_confidence_fusions']]
    
    ax2.bar(quality_labels, quality_counts)
    ax2.set_title('Fusion Quality Distribution')
    ax2.set_ylabel('Number of Fusions')
    
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Get inputs from snakemake
    arriba_file = snakemake.input.arriba
    fusioncatcher_file = snakemake.input.fusioncatcher
    integrated_file = snakemake.input.integrated
    
    output_report = snakemake.output.qc_report
    output_plot = snakemake.output.qc_plot
    
    # Analyze fusion quality
    metrics = analyze_fusion_quality(arriba_file, fusioncatcher_file, integrated_file)
    
    # Add quality flags
    if metrics['integrated_fusions'] == 0:
        metrics['quality_flag'] = 'NO_FUSIONS'
    elif metrics['consensus_fusions'] == 0:
        metrics['quality_flag'] = 'NO_CONSENSUS'
    elif metrics['high_confidence_fusions'] == 0:
        metrics['quality_flag'] = 'LOW_CONFIDENCE'
    else:
        metrics['quality_flag'] = 'PASS'
    
    # Save QC report
    with open(output_report, 'w') as f:
        json.dump(metrics, f, indent=2)
    
    # Create QC plot
    create_qc_plot(metrics, output_plot)
    
    print(f"Fusion QC completed: {metrics['quality_flag']}")

if __name__ == "__main__":
    main()