#!/usr/bin/env python3
"""
Pipeline Summary Script for IntegrateALL Pipeline
=================================================

Generates comprehensive pipeline summary report from all analysis results.
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from datetime import datetime

def load_sample_summaries(sample_files):
    """Load summary data from all samples"""
    all_samples = []
    
    for file_path in sample_files:
        try:
            with open(file_path, 'r') as f:
                sample_data = json.load(f)
                all_samples.append(sample_data)
        except Exception as e:
            print(f"Warning: Could not load {file_path}: {e}")
    
    return all_samples

def load_analysis_summaries(fusion_file, classification_file, cnv_file, variant_file):
    """Load analysis summary files"""
    summaries = {}
    
    # Load fusion summary
    try:
        summaries['fusions'] = pd.read_csv(fusion_file, sep='\t')
    except Exception as e:
        print(f"Warning: Could not load fusion summary: {e}")
        summaries['fusions'] = pd.DataFrame()
    
    # Load classification summary
    try:
        summaries['classifications'] = pd.read_csv(classification_file, sep='\t')
    except Exception as e:
        print(f"Warning: Could not load classification summary: {e}")
        summaries['classifications'] = pd.DataFrame()
    
    # Load CNV summary
    try:
        summaries['cnv'] = pd.read_csv(cnv_file, sep='\t')
    except Exception as e:
        print(f"Warning: Could not load CNV summary: {e}")
        summaries['cnv'] = pd.DataFrame()
    
    # Load variant summary
    try:
        summaries['variants'] = pd.read_csv(variant_file, sep='\t')
    except Exception as e:
        print(f"Warning: Could not load variant summary: {e}")
        summaries['variants'] = pd.DataFrame()
    
    return summaries

def generate_pipeline_statistics(sample_data, analysis_summaries):
    """Generate overall pipeline statistics"""
    
    stats = {
        'run_info': {
            'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'total_samples': len(sample_data),
            'pipeline_version': '1.0.0'
        },
        'sample_stats': {},
        'analysis_stats': {}
    }
    
    if sample_data:
        # Sample-level statistics
        sample_ids = [s.get('sample_id', 'unknown') for s in sample_data]
        stats['sample_stats'] = {
            'sample_count': len(sample_ids),
            'samples_processed': len([s for s in sample_data if s.get('status') == 'completed']),
            'samples_failed': len([s for s in sample_data if s.get('status') == 'failed'])
        }
    
    # Analysis-specific statistics
    if not analysis_summaries['fusions'].empty:
        stats['analysis_stats']['fusions'] = {
            'samples_with_fusions': len(analysis_summaries['fusions']),
            'total_fusions': analysis_summaries['fusions']['total_fusions'].sum() if 'total_fusions' in analysis_summaries['fusions'].columns else 0
        }
    
    if not analysis_summaries['classifications'].empty:
        stats['analysis_stats']['classifications'] = {
            'classified_samples': len(analysis_summaries['classifications']),
            'subtype_distribution': analysis_summaries['classifications']['predicted_subtype'].value_counts().to_dict() if 'predicted_subtype' in analysis_summaries['classifications'].columns else {}
        }
    
    if not analysis_summaries['cnv'].empty:
        stats['analysis_stats']['cnv'] = {
            'samples_with_cnv': len(analysis_summaries['cnv']),
            'avg_alterations': analysis_summaries['cnv']['total_alterations'].mean() if 'total_alterations' in analysis_summaries['cnv'].columns else 0
        }
    
    if not analysis_summaries['variants'].empty:
        stats['analysis_stats']['variants'] = {
            'samples_with_variants': len(analysis_summaries['variants']),
            'total_hotspot_variants': analysis_summaries['variants']['hotspot_variants'].sum() if 'hotspot_variants' in analysis_summaries['variants'].columns else 0
        }
    
    return stats

def create_pipeline_plots(analysis_summaries, output_dir):
    """Create summary plots"""
    
    plots_created = []
    
    # Classification distribution plot
    if not analysis_summaries['classifications'].empty and 'predicted_subtype' in analysis_summaries['classifications'].columns:
        plt.figure(figsize=(10, 6))
        subtype_counts = analysis_summaries['classifications']['predicted_subtype'].value_counts()
        
        plt.pie(subtype_counts.values, labels=subtype_counts.index, autopct='%1.1f%%')
        plt.title('B-ALL Subtype Distribution')
        plot_path = f"{output_dir}/subtype_distribution.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plots_created.append(plot_path)
    
    # Fusion analysis plot
    if not analysis_summaries['fusions'].empty and 'total_fusions' in analysis_summaries['fusions'].columns:
        plt.figure(figsize=(10, 6))
        plt.hist(analysis_summaries['fusions']['total_fusions'], bins=20, alpha=0.7)
        plt.xlabel('Number of Fusions per Sample')
        plt.ylabel('Sample Count')
        plt.title('Fusion Distribution Across Samples')
        plot_path = f"{output_dir}/fusion_distribution.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plots_created.append(plot_path)
    
    # CNV analysis plot
    if not analysis_summaries['cnv'].empty and 'total_alterations' in analysis_summaries['cnv'].columns:
        plt.figure(figsize=(10, 6))
        plt.hist(analysis_summaries['cnv']['total_alterations'], bins=20, alpha=0.7)
        plt.xlabel('Number of CNV Alterations per Sample')
        plt.ylabel('Sample Count')
        plt.title('CNV Alteration Distribution')
        plot_path = f"{output_dir}/cnv_distribution.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        plots_created.append(plot_path)
    
    return plots_created

def generate_html_report(stats, plots, output_html):
    """Generate HTML summary report"""
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>IntegrateALL Pipeline Summary</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
            .section {{ margin: 20px 0; }}
            .stats-table {{ border-collapse: collapse; width: 100%; }}
            .stats-table th, .stats-table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            .stats-table th {{ background-color: #4CAF50; color: white; }}
            .plot {{ text-align: center; margin: 20px 0; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>IntegrateALL Pipeline Summary Report</h1>
            <p><strong>Analysis Date:</strong> {stats['run_info']['analysis_date']}</p>
            <p><strong>Pipeline Version:</strong> {stats['run_info']['pipeline_version']}</p>
            <p><strong>Total Samples:</strong> {stats['run_info']['total_samples']}</p>
        </div>
        
        <div class="section">
            <h2>Sample Processing Summary</h2>
            <table class="stats-table">
                <tr><th>Metric</th><th>Value</th></tr>
                <tr><td>Total Samples</td><td>{stats.get('sample_stats', {}).get('sample_count', 0)}</td></tr>
                <tr><td>Successfully Processed</td><td>{stats.get('sample_stats', {}).get('samples_processed', 0)}</td></tr>
                <tr><td>Failed Samples</td><td>{stats.get('sample_stats', {}).get('samples_failed', 0)}</td></tr>
            </table>
        </div>
        
        <div class="section">
            <h2>Analysis Results</h2>
            <h3>Classification</h3>
            <p>Classified Samples: {stats.get('analysis_stats', {}).get('classifications', {}).get('classified_samples', 0)}</p>
            
            <h3>Fusion Detection</h3>
            <p>Samples with Fusions: {stats.get('analysis_stats', {}).get('fusions', {}).get('samples_with_fusions', 0)}</p>
            <p>Total Fusions: {stats.get('analysis_stats', {}).get('fusions', {}).get('total_fusions', 0)}</p>
            
            <h3>CNV Analysis</h3>
            <p>Samples with CNVs: {stats.get('analysis_stats', {}).get('cnv', {}).get('samples_with_cnv', 0)}</p>
            <p>Average Alterations: {stats.get('analysis_stats', {}).get('cnv', {}).get('avg_alterations', 0):.2f}</p>
            
            <h3>Variant Analysis</h3>
            <p>Samples with Variants: {stats.get('analysis_stats', {}).get('variants', {}).get('samples_with_variants', 0)}</p>
            <p>Total Hotspot Variants: {stats.get('analysis_stats', {}).get('variants', {}).get('total_hotspot_variants', 0)}</p>
        </div>
    """
    
    # Add plots if available
    if plots:
        html_content += '<div class="section"><h2>Analysis Plots</h2>'
        for plot_path in plots:
            plot_name = Path(plot_path).name
            html_content += f'<div class="plot"><img src="{plot_name}" alt="{plot_name}" style="max-width: 800px;"></div>'
        html_content += '</div>'
    
    html_content += """
        </body>
        </html>
    """
    
    with open(output_html, 'w') as f:
        f.write(html_content)

def main():
    # Get inputs from snakemake
    sample_summaries = snakemake.input.sample_summaries
    fusion_summary = snakemake.input.fusion_summary
    classification_summary = snakemake.input.classification_summary
    cnv_summary = snakemake.input.cnv_summary
    variant_summary = snakemake.input.variant_summary
    
    output_html = snakemake.output.html_report
    output_json = snakemake.output.json_summary
    
    pipeline_version = snakemake.params.pipeline_version
    analysis_date = snakemake.params.analysis_date
    
    print("Generating pipeline summary report...")
    
    # Load all data
    sample_data = load_sample_summaries(sample_summaries)
    analysis_summaries = load_analysis_summaries(
        fusion_summary, classification_summary, cnv_summary, variant_summary
    )
    
    # Generate statistics
    stats = generate_pipeline_statistics(sample_data, analysis_summaries)
    stats['run_info']['pipeline_version'] = pipeline_version
    stats['run_info']['analysis_date'] = analysis_date
    
    # Create output directory for plots
    output_dir = Path(output_html).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create plots
    plots = create_pipeline_plots(analysis_summaries, str(output_dir))
    
    # Generate HTML report
    generate_html_report(stats, plots, output_html)
    
    # Save JSON summary
    with open(output_json, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"Pipeline summary completed:")
    print(f"  HTML report: {output_html}")
    print(f"  JSON summary: {output_json}")
    print(f"  Samples processed: {stats.get('sample_stats', {}).get('sample_count', 0)}")
    print(f"  Plots created: {len(plots)}")

if __name__ == "__main__":
    main()