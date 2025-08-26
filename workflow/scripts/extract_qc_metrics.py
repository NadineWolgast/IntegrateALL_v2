#!/usr/bin/env python3
"""
Extract QC Metrics Script for IntegrateALL Pipeline
===================================================

Extracts key quality metrics from FastQC results for downstream analysis.
"""

import json
import zipfile
import re
import sys
from pathlib import Path

def parse_fastqc_data(fastqc_zip_path):
    """Parse FastQC data from ZIP file"""
    metrics = {}
    
    try:
        with zipfile.ZipFile(fastqc_zip_path, 'r') as zip_file:
            # Get the base name (without _fastqc.zip)
            base_name = Path(fastqc_zip_path).stem.replace('_fastqc', '')
            
            # Read fastqc_data.txt
            data_file = f"{base_name}_fastqc/fastqc_data.txt"
            
            if data_file in zip_file.namelist():
                with zip_file.open(data_file) as f:
                    content = f.read().decode('utf-8')
                    metrics.update(parse_fastqc_content(content))
            
            # Read summary.txt for pass/warn/fail status
            summary_file = f"{base_name}_fastqc/summary.txt"
            if summary_file in zip_file.namelist():
                with zip_file.open(summary_file) as f:
                    content = f.read().decode('utf-8')
                    metrics.update(parse_summary_content(content))
                    
    except Exception as e:
        print(f"Error parsing FastQC file {fastqc_zip_path}: {e}")
        
    return metrics

def parse_fastqc_content(content):
    """Parse the main FastQC data content"""
    metrics = {}
    lines = content.strip().split('\n')
    
    current_module = None
    
    for line in lines:
        line = line.strip()
        
        # Skip empty lines and comments
        if not line or line.startswith('#'):
            continue
        
        # Check for module headers
        if line.startswith('>>'):
            if line.endswith('END_MODULE'):
                current_module = None
            else:
                current_module = line[2:].strip()
            continue
        
        # Parse basic statistics
        if current_module == "Basic Statistics":
            parts = line.split('\t')
            if len(parts) >= 2:
                key = parts[0].lower().replace(' ', '_').replace('(', '').replace(')', '').replace('%', 'percent')
                value = parts[1]
                
                # Convert to appropriate type
                try:
                    if '.' in value:
                        value = float(value)
                    elif value.isdigit():
                        value = int(value)
                except ValueError:
                    pass  # Keep as string
                    
                metrics[key] = value
        
        # Parse per base sequence quality (take average)
        elif current_module == "Per base sequence quality":
            parts = line.split('\t')
            if len(parts) >= 2 and not line.startswith('#Base'):
                try:
                    mean_quality = float(parts[1])
                    if 'per_base_qualities' not in metrics:
                        metrics['per_base_qualities'] = []
                    metrics['per_base_qualities'].append(mean_quality)
                except ValueError:
                    continue
        
        # Parse sequence duplication levels
        elif current_module == "Sequence Duplication Levels":
            if 'Total Duplicate Percentage' in line:
                try:
                    dup_percent = float(line.split('\t')[1])
                    metrics['duplicate_percentage'] = dup_percent
                except (ValueError, IndexError):
                    continue
        
        # Parse overrepresented sequences
        elif current_module == "Overrepresented sequences":
            if not line.startswith('#Sequence') and len(line.split('\t')) >= 3:
                if 'overrepresented_sequences' not in metrics:
                    metrics['overrepresented_sequences'] = 0
                metrics['overrepresented_sequences'] += 1
    
    # Calculate average per-base quality if available
    if 'per_base_qualities' in metrics:
        metrics['avg_per_base_quality'] = sum(metrics['per_base_qualities']) / len(metrics['per_base_qualities'])
        del metrics['per_base_qualities']  # Remove the list to save space
    
    return metrics

def parse_summary_content(content):
    """Parse FastQC summary content for pass/warn/fail status"""
    metrics = {}
    lines = content.strip().split('\n')
    
    module_status = {}
    
    for line in lines:
        parts = line.split('\t')
        if len(parts) >= 3:
            status = parts[0]
            module = parts[1]
            module_status[module] = status
    
    # Count pass/warn/fail modules
    metrics['modules_passed'] = sum(1 for status in module_status.values() if status == 'PASS')
    metrics['modules_warned'] = sum(1 for status in module_status.values() if status == 'WARN')
    metrics['modules_failed'] = sum(1 for status in module_status.values() if status == 'FAIL')
    metrics['total_modules'] = len(module_status)
    
    # Calculate overall quality score
    if metrics['total_modules'] > 0:
        quality_score = (metrics['modules_passed'] * 1.0 + 
                        metrics['modules_warned'] * 0.5 + 
                        metrics['modules_failed'] * 0.0) / metrics['total_modules']
        metrics['overall_quality_score'] = quality_score
    
    # Store individual module statuses for key modules
    key_modules = [
        'Per base sequence quality',
        'Per sequence quality scores', 
        'Per base N content',
        'Sequence Duplication Levels',
        'Overrepresented sequences'
    ]
    
    for module in key_modules:
        if module in module_status:
            safe_name = module.lower().replace(' ', '_').replace('-', '_')
            metrics[f"{safe_name}_status"] = module_status[module]
    
    return metrics

def combine_paired_metrics(r1_metrics, r2_metrics, sample_id):
    """Combine metrics from R1 and R2 files"""
    combined = {
        'sample_id': sample_id,
        'read_pair': 'paired'
    }
    
    # Average numeric metrics
    numeric_keys = [
        'total_sequences', 'sequence_length', 'percent_gc',
        'avg_per_base_quality', 'overall_quality_score',
        'duplicate_percentage'
    ]
    
    for key in numeric_keys:
        r1_val = r1_metrics.get(key, 0)
        r2_val = r2_metrics.get(key, 0)
        
        if isinstance(r1_val, (int, float)) and isinstance(r2_val, (int, float)):
            combined[key] = (r1_val + r2_val) / 2
        else:
            combined[key] = r1_val  # Use R1 value if not numeric
    
    # Sum certain metrics
    sum_keys = [
        'modules_passed', 'modules_warned', 'modules_failed', 'total_modules',
        'overrepresented_sequences'
    ]
    
    for key in sum_keys:
        r1_val = r1_metrics.get(key, 0)
        r2_val = r2_metrics.get(key, 0)
        combined[key] = r1_val + r2_val
    
    # Take worst status for key modules
    status_keys = [k for k in r1_metrics.keys() if k.endswith('_status')]
    for key in status_keys:
        r1_status = r1_metrics.get(key, 'PASS')
        r2_status = r2_metrics.get(key, 'PASS')
        
        # FAIL > WARN > PASS priority
        if r1_status == 'FAIL' or r2_status == 'FAIL':
            combined[key] = 'FAIL'
        elif r1_status == 'WARN' or r2_status == 'WARN':
            combined[key] = 'WARN'
        else:
            combined[key] = 'PASS'
    
    return combined

def main():
    # Get inputs from snakemake
    fastqc_r1 = snakemake.input.fastqc_r1
    fastqc_r2 = snakemake.input.fastqc_r2
    output_metrics = snakemake.output.metrics
    sample_id = Path(fastqc_r1).stem.split('_R1_fastqc')[0]
    
    # Parse FastQC files
    r1_metrics = parse_fastqc_data(fastqc_r1)
    r2_metrics = parse_fastqc_data(fastqc_r2)
    
    # Combine metrics from both reads
    combined_metrics = combine_paired_metrics(r1_metrics, r2_metrics, sample_id)
    
    # Add file paths for reference
    combined_metrics['fastqc_r1_file'] = str(fastqc_r1)
    combined_metrics['fastqc_r2_file'] = str(fastqc_r2)
    
    # Save metrics as JSON
    with open(output_metrics, 'w') as f:
        json.dump(combined_metrics, f, indent=2)
    
    print(f"Extracted QC metrics for sample {sample_id}")
    print(f"Overall quality score: {combined_metrics.get('overall_quality_score', 'N/A')}")

if __name__ == "__main__":
    main()