#!/usr/bin/env python3
"""
Final Report Generator for IntegrateALL Pipeline
===============================================

Generates comprehensive HTML reports for each sample with integrated results.
This script replaces the original bash/R reporting scripts.
"""

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import base64
from io import BytesIO
import argparse
import logging
from typing import Dict, List, Any
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ReportGenerator:
    """Generates comprehensive HTML reports for B-ALL analysis"""
    
    def __init__(self, sample_id: str):
        """Initialize report generator"""
        self.sample_id = sample_id
        self.report_data = {}
        
        # Set matplotlib style
        plt.style.use('default')
        sns.set_palette("husl")
    
    def load_all_results(self, file_paths: Dict[str, str]) -> None:
        """Load all analysis results"""
        
        # Load QC results
        self.report_data['qc'] = self._load_json_safe(file_paths.get('qc_metrics', ''))
        
        # Load alignment stats
        self.report_data['alignment'] = self._load_json_safe(file_paths.get('alignment_stats', ''))
        
        # Load fusion results
        self.report_data['fusions'] = self._load_tsv_safe(file_paths.get('fusion_results', ''))
        self.report_data['fusion_summary'] = self._load_json_safe(file_paths.get('fusion_summary', ''))
        
        # Load classification results
        self.report_data['classification'] = self._load_json_safe(file_paths.get('classification', ''))
        self.report_data['classification_confidence'] = self._load_json_safe(file_paths.get('classification_confidence', ''))
        
        # Load CNV results
        self.report_data['cnv'] = self._load_tsv_safe(file_paths.get('cnv_results', ''))
        
        # Load variant results
        self.report_data['hotspots'] = self._load_tsv_safe(file_paths.get('hotspots', ''))
        
    def _load_json_safe(self, file_path: str) -> Dict:
        """Safely load JSON file"""
        try:
            if file_path and Path(file_path).exists():
                with open(file_path, 'r') as f:
                    return json.load(f)
        except Exception as e:
            logger.warning(f"Could not load JSON file {file_path}: {e}")
        return {}
    
    def _load_tsv_safe(self, file_path: str) -> pd.DataFrame:
        """Safely load TSV file"""
        try:
            if file_path and Path(file_path).exists():
                return pd.read_csv(file_path, sep='\\t')
        except Exception as e:
            logger.warning(f"Could not load TSV file {file_path}: {e}")
        return pd.DataFrame()
    
    def create_qc_plot(self) -> str:
        """Create QC summary plot"""
        try:
            qc_data = self.report_data.get('qc', {})
            
            if not qc_data:
                return ""
            
            fig, axes = plt.subplots(2, 2, figsize=(12, 8))
            fig.suptitle(f'Quality Control Summary - {self.sample_id}', fontsize=16)
            
            # Plot 1: Read counts
            if 'total_sequences' in qc_data:
                axes[0, 0].bar(['R1', 'R2'], 
                              [qc_data.get('total_sequences_r1', 0), qc_data.get('total_sequences_r2', 0)])
                axes[0, 0].set_title('Total Sequences')
                axes[0, 0].set_ylabel('Number of reads')
            
            # Plot 2: Quality scores
            if 'mean_quality' in qc_data:
                axes[0, 1].bar(['R1', 'R2'], 
                              [qc_data.get('mean_quality_r1', 0), qc_data.get('mean_quality_r2', 0)])
                axes[0, 1].set_title('Mean Quality Scores')
                axes[0, 1].set_ylabel('Quality score')
                axes[0, 1].axhline(y=30, color='r', linestyle='--', label='Q30')
                axes[0, 1].legend()
            
            # Plot 3: GC content
            if 'gc_content' in qc_data:
                axes[1, 0].bar(['R1', 'R2'], 
                              [qc_data.get('gc_content_r1', 0), qc_data.get('gc_content_r2', 0)])
                axes[1, 0].set_title('GC Content (%)')
                axes[1, 0].set_ylabel('GC %')
                axes[1, 0].axhline(y=50, color='g', linestyle='--', label='Expected')
                axes[1, 0].legend()
            
            # Plot 4: Adapter content
            if 'adapter_content' in qc_data:
                axes[1, 1].bar(['R1', 'R2'], 
                              [qc_data.get('adapter_content_r1', 0), qc_data.get('adapter_content_r2', 0)])
                axes[1, 1].set_title('Adapter Content (%)')
                axes[1, 1].set_ylabel('Adapter %')
            
            plt.tight_layout()
            return self._plot_to_base64(fig)
            
        except Exception as e:
            logger.error(f"Error creating QC plot: {e}")
            return ""
    
    def create_fusion_plot(self) -> str:
        """Create fusion summary plot"""
        try:
            fusion_df = self.report_data.get('fusions', pd.DataFrame())
            
            if fusion_df.empty:
                return ""
            
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            fig.suptitle(f'Fusion Detection Summary - {self.sample_id}', fontsize=16)
            
            # Plot 1: Fusion support
            if 'total_support_merged' in fusion_df.columns:
                top_fusions = fusion_df.nlargest(10, 'total_support_merged')
                axes[0].barh(range(len(top_fusions)), top_fusions['total_support_merged'])
                axes[0].set_yticks(range(len(top_fusions)))
                axes[0].set_yticklabels(top_fusions['canonical_name'])
                axes[0].set_xlabel('Supporting Reads')
                axes[0].set_title('Top Fusions by Support')
            
            # Plot 2: Fusion types
            if 'known_fusion_type' in fusion_df.columns:
                fusion_types = fusion_df['known_fusion_type'].value_counts()
                if len(fusion_types) > 0:
                    axes[1].pie(fusion_types.values, labels=fusion_types.index, autopct='%1.1f%%')
                    axes[1].set_title('Fusion Types')
            
            plt.tight_layout()
            return self._plot_to_base64(fig)
            
        except Exception as e:
            logger.error(f"Error creating fusion plot: {e}")
            return ""
    
    def create_classification_plot(self) -> str:
        """Create classification confidence plot"""
        try:
            classification = self.report_data.get('classification', {})
            confidence_data = self.report_data.get('classification_confidence', {})
            
            if not confidence_data:
                return ""
            
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            fig.suptitle(f'Classification Results - {self.sample_id}', fontsize=16)
            
            # Plot 1: Confidence scores by method
            individual_conf = confidence_data.get('individual_confidences', {})
            if individual_conf:
                methods = list(individual_conf.keys())
                scores = list(individual_conf.values())
                
                bars = axes[0].bar(methods, scores)
                axes[0].set_title('Confidence by Method')
                axes[0].set_ylabel('Confidence Score')
                axes[0].set_ylim(0, 1)
                axes[0].axhline(y=0.7, color='r', linestyle='--', label='Threshold')
                axes[0].legend()
                
                # Color bars based on confidence
                for bar, score in zip(bars, scores):
                    if score >= 0.7:
                        bar.set_color('green')
                    elif score >= 0.5:
                        bar.set_color('orange') 
                    else:
                        bar.set_color('red')
            
            # Plot 2: Final classification
            final_class = classification.get('final_classification', 'Unknown')
            overall_conf = confidence_data.get('overall_confidence', 0.0)
            
            # Create a simple confidence meter
            theta = np.linspace(0, np.pi, 100)
            radius = 1
            x = radius * np.cos(theta)
            y = radius * np.sin(theta)
            
            axes[1].plot(x, y, 'k-', linewidth=2)
            axes[1].fill_between(x[:int(overall_conf*100)], y[:int(overall_conf*100)], 0, 
                               color='green' if overall_conf >= 0.7 else 'orange' if overall_conf >= 0.5 else 'red',
                               alpha=0.6)
            axes[1].set_xlim(-1.2, 1.2)
            axes[1].set_ylim(-0.2, 1.2)
            axes[1].set_aspect('equal')
            axes[1].axis('off')
            axes[1].text(0, -0.1, f'{final_class}\\n{overall_conf:.2f}', 
                        ha='center', va='center', fontsize=14, fontweight='bold')
            axes[1].set_title('Final Classification & Confidence')
            
            plt.tight_layout()
            return self._plot_to_base64(fig)
            
        except Exception as e:
            logger.error(f"Error creating classification plot: {e}")
            return ""
    
    def _plot_to_base64(self, fig) -> str:
        """Convert matplotlib figure to base64 string"""
        buffer = BytesIO()
        fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
        buffer.seek(0)
        plot_data = buffer.getvalue()
        buffer.close()
        plt.close(fig)
        
        return base64.b64encode(plot_data).decode()
    
    def generate_html_report(self, output_file: str) -> None:
        """Generate comprehensive HTML report"""
        
        # Create plots
        qc_plot = self.create_qc_plot()
        fusion_plot = self.create_fusion_plot()
        classification_plot = self.create_classification_plot()
        
        # Get summary data
        classification = self.report_data.get('classification', {})
        fusion_summary = self.report_data.get('fusion_summary', {})
        alignment_stats = self.report_data.get('alignment', {})
        
        html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>IntegrateALL Report - {self.sample_id}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        .section {{
            margin-bottom: 30px;
            padding: 20px;
            border: 1px solid #ddd;
            border-radius: 8px;
        }}
        .section h2 {{
            color: #333;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }}
        .plot {{
            text-align: center;
            margin: 20px 0;
        }}
        .plot img {{
            max-width: 100%;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .summary-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .summary-card h3 {{
            margin-top: 0;
            color: #667eea;
        }}
        .key-value {{
            display: flex;
            justify-content: space-between;
            padding: 5px 0;
            border-bottom: 1px solid #eee;
        }}
        .key-value:last-child {{
            border-bottom: none;
        }}
        .status {{
            padding: 5px 15px;
            border-radius: 20px;
            color: white;
            font-weight: bold;
        }}
        .status.high {{ background-color: #28a745; }}
        .status.medium {{ background-color: #ffc107; color: black; }}
        .status.low {{ background-color: #dc3545; }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #667eea;
            color: white;
        }}
        .footer {{
            text-align: center;
            margin-top: 40px;
            padding: 20px;
            background-color: #f8f9fa;
            border-radius: 8px;
            color: #666;
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <h1>IntegrateALL Analysis Report</h1>
            <h2>Sample: {self.sample_id}</h2>
            <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <!-- Executive Summary -->
        <div class="section">
            <h2>Executive Summary</h2>
            <div class="summary-grid">
                <div class="summary-card">
                    <h3>Classification Result</h3>
                    <div class="key-value">
                        <span>Subtype:</span>
                        <span><strong>{classification.get('final_classification', 'Unknown')}</strong></span>
                    </div>
                    <div class="key-value">
                        <span>Confidence:</span>
                        <span class="status {'high' if classification.get('classification_confidence', 0) >= 0.7 else 'medium' if classification.get('classification_confidence', 0) >= 0.5 else 'low'}">
                            {classification.get('classification_confidence', 0):.2f}
                        </span>
                    </div>
                </div>
                
                <div class="summary-card">
                    <h3>Fusion Detection</h3>
                    <div class="key-value">
                        <span>Total Fusions:</span>
                        <span>{fusion_summary.get('total_fusions', 0)}</span>
                    </div>
                    <div class="key-value">
                        <span>Driver Fusions:</span>
                        <span>{fusion_summary.get('driver_fusions', 0)}</span>
                    </div>
                    <div class="key-value">
                        <span>High Confidence:</span>
                        <span>{fusion_summary.get('high_confidence_fusions', 0)}</span>
                    </div>
                </div>
                
                <div class="summary-card">
                    <h3>Sequencing Quality</h3>
                    <div class="key-value">
                        <span>Total Reads:</span>
                        <span>{alignment_stats.get('total_reads', 'N/A'):,}</span>
                    </div>
                    <div class="key-value">
                        <span>Mapped Reads:</span>
                        <span>{alignment_stats.get('mapped_reads', 'N/A'):,}</span>
                    </div>
                    <div class="key-value">
                        <span>Mapping Rate:</span>
                        <span>{alignment_stats.get('mapping_rate', 'N/A')}</span>
                    </div>
                </div>
            </div>
        </div>
        
        <!-- Quality Control -->
        <div class="section">
            <h2>Quality Control</h2>
            {f'<div class="plot"><img src="data:image/png;base64,{qc_plot}" alt="QC Plot"></div>' if qc_plot else '<p>QC plot not available</p>'}
        </div>
        
        <!-- Classification Results -->
        <div class="section">
            <h2>Classification Results</h2>
            {f'<div class="plot"><img src="data:image/png;base64,{classification_plot}" alt="Classification Plot"></div>' if classification_plot else '<p>Classification plot not available</p>'}
            
            <h3>Classification Details</h3>
            <table>
                <tr><th>Parameter</th><th>Value</th></tr>
                <tr><td>Final Classification</td><td>{classification.get('final_classification', 'Unknown')}</td></tr>
                <tr><td>Confidence Score</td><td>{classification.get('classification_confidence', 0):.3f}</td></tr>
                <tr><td>Method</td><td>{classification.get('classification_method', 'Unknown')}</td></tr>
            </table>
        </div>
        
        <!-- Fusion Detection -->
        <div class="section">
            <h2>Fusion Detection</h2>
            {f'<div class="plot"><img src="data:image/png;base64,{fusion_plot}" alt="Fusion Plot"></div>' if fusion_plot else '<p>Fusion plot not available</p>'}
            
            <h3>Top Detected Fusions</h3>
            {self._create_fusion_table()}
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <p>Report generated by IntegrateALL Pipeline v1.0</p>
            <p>For questions or support, please contact the development team</p>
        </div>
    </div>
</body>
</html>"""

        # Write HTML report
        with open(output_file, 'w') as f:
            f.write(html_template)
        
        logger.info(f"HTML report generated: {output_file}")
    
    def _create_fusion_table(self) -> str:
        """Create HTML table for fusion results"""
        fusion_df = self.report_data.get('fusions', pd.DataFrame())
        
        if fusion_df.empty:
            return "<p>No fusions detected</p>"
        
        # Get top 10 fusions by support
        top_fusions = fusion_df.nlargest(10, 'total_support_merged')
        
        table_html = "<table><tr><th>Fusion</th><th>Support</th><th>Confidence</th><th>Type</th><th>Callers</th></tr>"
        
        for _, fusion in top_fusions.iterrows():
            table_html += f"""
            <tr>
                <td>{fusion.get('canonical_name', 'Unknown')}</td>
                <td>{fusion.get('total_support_merged', 0)}</td>
                <td><span class="status {fusion.get('confidence_merged', 'low')}">{fusion.get('confidence_merged', 'low')}</span></td>
                <td>{fusion.get('known_fusion_type', 'unknown')}</td>
                <td>{fusion.get('callers', 'unknown')}</td>
            </tr>"""
        
        table_html += "</table>"
        return table_html

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Generate comprehensive HTML report')
    parser.add_argument('sample_id', help='Sample identifier')
    parser.add_argument('output_html', help='Output HTML report file')
    parser.add_argument('--qc-metrics', help='QC metrics JSON file')
    parser.add_argument('--alignment-stats', help='Alignment stats JSON file')
    parser.add_argument('--fusion-results', help='Fusion results TSV file')
    parser.add_argument('--fusion-summary', help='Fusion summary JSON file')
    parser.add_argument('--classification', help='Classification results JSON file')
    parser.add_argument('--classification-confidence', help='Classification confidence JSON file')
    parser.add_argument('--cnv-results', help='CNV results TSV file')
    parser.add_argument('--hotspots', help='Hotspot mutations TSV file')
    
    args = parser.parse_args()
    
    # Initialize report generator
    reporter = ReportGenerator(args.sample_id)
    
    # Collect file paths
    file_paths = {
        'qc_metrics': args.qc_metrics,
        'alignment_stats': args.alignment_stats,
        'fusion_results': args.fusion_results,
        'fusion_summary': args.fusion_summary,
        'classification': args.classification,
        'classification_confidence': args.classification_confidence,
        'cnv_results': args.cnv_results,
        'hotspots': args.hotspots
    }
    
    # Load all results
    logger.info("Loading analysis results...")
    reporter.load_all_results(file_paths)
    
    # Generate HTML report
    logger.info("Generating HTML report...")
    reporter.generate_html_report(args.output_html)
    
    # Generate summary JSON
    summary_file = args.output_html.replace('.html', '_summary.json')
    summary_data = {
        'sample_id': args.sample_id,
        'report_generated': datetime.now().isoformat(),
        'classification': reporter.report_data.get('classification', {}),
        'fusion_summary': reporter.report_data.get('fusion_summary', {}),
        'alignment_summary': reporter.report_data.get('alignment', {})
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary_data, f, indent=2)
    
    logger.info("Report generation completed successfully")

if __name__ == "__main__":
    main()