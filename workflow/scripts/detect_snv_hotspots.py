#!/usr/bin/env python3
"""
SNV Hotspot Detection Script for IntegrateALL Pipeline
======================================================
Detects specific B-ALL SNV hotspots using pysamstats.
Maintains compatibility with original classification system.
"""

import os
import pandas as pd
import subprocess
from pathlib import Path
import sys

# Critical SNV hotspots for B-ALL classification
HOTSPOTS = {
    "PAX5_P80R": {
        "chromosome": "9",
        "position": 36838531,
        "ref": "C", 
        "alt": "G",
        "gene": "PAX5",
        "mutation": "P80R"
    },
    "IKZF1_N159Y": {
        "chromosome": "7", 
        "position": 50444501,
        "ref": "A",
        "alt": "T", 
        "gene": "IKZF1",
        "mutation": "N159Y"
    },
    "ZEB2_H1038R": {
        "chromosome": "2",
        "position": 145141514,
        "ref": "A", 
        "alt": "G",
        "gene": "ZEB2", 
        "mutation": "H1038R"
    }
}

def run_pysamstats(bam_file: str, reference: str, chrom: str, pos: int, 
                   window: int = 1) -> pd.DataFrame:
    """Run pysamstats to get variant information at specific position."""
    try:
        cmd = [
            "pysamstats", "--type", "variation",
            "--chromosome", chrom,
            "--start", str(pos - window),
            "--end", str(pos + window),
            "--fasta", reference,
            bam_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse pysamstats output
        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:
            return pd.DataFrame()
            
        header = lines[0].split('\t')
        data_rows = []
        
        for line in lines[1:]:
            data_rows.append(line.split('\t'))
        
        return pd.DataFrame(data_rows, columns=header)
        
    except subprocess.CalledProcessError as e:
        print(f"Warning: pysamstats failed for {chrom}:{pos} - {e}")
        return pd.DataFrame()
    except Exception as e:
        print(f"Error running pysamstats: {e}")
        return pd.DataFrame()

def check_hotspot_mutation(stats_df: pd.DataFrame, hotspot: dict, 
                          min_depth: int = 10, min_alt_freq: float = 0.1) -> bool:
    """Check if hotspot mutation is present based on pysamstats output."""
    if stats_df.empty:
        return False
    
    try:
        # Find the exact position
        pos_data = stats_df[stats_df['pos'] == str(hotspot['position'])]
        if pos_data.empty:
            return False
        
        row = pos_data.iloc[0]
        
        # Get coverage and allele counts
        total_depth = int(row.get('reads_all', 0))
        
        # Check specific nucleotide counts based on mutation
        alt_base = hotspot['alt']
        alt_count = int(row.get(f'{alt_base}', 0))
        
        if total_depth < min_depth:
            return False
        
        alt_frequency = alt_count / total_depth if total_depth > 0 else 0
        
        return alt_frequency >= min_alt_freq
        
    except Exception as e:
        print(f"Error checking hotspot mutation: {e}")
        return False

def detect_hotspots(bam_file: str, reference_fasta: str, output_dir: str, 
                   sample_id: str) -> dict:
    """Detect all B-ALL SNV hotspots and create marker files."""
    
    detected_hotspots = {}
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Detecting SNV hotspots for sample: {sample_id}")
    
    for hotspot_name, hotspot_info in HOTSPOTS.items():
        print(f"Checking {hotspot_name} at {hotspot_info['chromosome']}:{hotspot_info['position']}")
        
        # Run pysamstats for this position
        stats_df = run_pysamstats(
            bam_file, reference_fasta,
            hotspot_info['chromosome'], hotspot_info['position']
        )
        
        # Check if mutation is present
        mutation_found = check_hotspot_mutation(stats_df, hotspot_info)
        
        detected_hotspots[hotspot_name] = mutation_found
        
        # Create marker file if mutation found (original system compatibility)
        if mutation_found:
            marker_file = os.path.join(output_dir, f"{hotspot_name}_{sample_id}.detected")
            with open(marker_file, 'w') as f:
                f.write(f"Hotspot detected: {hotspot_name}\n")
                f.write(f"Sample: {sample_id}\n")
                f.write(f"Position: {hotspot_info['chromosome']}:{hotspot_info['position']}\n")
                f.write(f"Mutation: {hotspot_info['gene']} {hotspot_info['mutation']}\n")
            
            print(f"✓ {hotspot_name} detected and marked")
        else:
            print(f"✗ {hotspot_name} not detected")
    
    return detected_hotspots

def create_summary_report(detected_hotspots: dict, output_file: str, sample_id: str):
    """Create summary report of detected hotspots."""
    
    summary_data = []
    for hotspot_name, detected in detected_hotspots.items():
        hotspot_info = HOTSPOTS[hotspot_name]
        summary_data.append({
            'sample_id': sample_id,
            'hotspot': hotspot_name,
            'gene': hotspot_info['gene'],
            'mutation': hotspot_info['mutation'],
            'chromosome': hotspot_info['chromosome'],
            'position': hotspot_info['position'],
            'detected': detected
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_file, index=False)
    
    print(f"SNV hotspot summary saved to: {output_file}")

def main():
    """Main function using snakemake inputs."""
    # Get inputs from snakemake
    bam_file = snakemake.input.bam
    reference_fasta = snakemake.input.reference
    
    # Outputs
    output_dir = snakemake.params.hotspot_dir
    output_summary = snakemake.output.summary
    
    sample_id = snakemake.wildcards.sample
    
    try:
        # Detect hotspots
        detected_hotspots = detect_hotspots(bam_file, reference_fasta, output_dir, sample_id)
        
        # Create summary report
        create_summary_report(detected_hotspots, output_summary, sample_id)
        
        print(f"SNV hotspot detection completed for sample: {sample_id}")
        
    except Exception as e:
        print(f"Error in SNV hotspot detection: {e}")
        # Create empty summary on error
        empty_df = pd.DataFrame()
        empty_df.to_csv(output_summary, index=False)

if __name__ == "__main__":
    main()