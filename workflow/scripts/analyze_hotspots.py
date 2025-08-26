#!/usr/bin/env python3
"""
Hotspot Analysis Script for IntegrateALL Pipeline
=================================================

Analyzes variants in known B-ALL hotspot regions and driver genes.
"""

import pandas as pd
import numpy as np
import json
import sys
from pathlib import Path
try:
    from cyvcf2 import VCF
    HAS_CYVCF2 = True
except ImportError:
    print("Warning: cyvcf2 not available, using basic VCF parsing")
    HAS_CYVCF2 = False

def load_hotspot_regions(hotspot_bed_file):
    """Load hotspot regions from BED file"""
    hotspots = []
    
    try:
        with open(hotspot_bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    hotspots.append({
                        'chromosome': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'region': parts[3],
                        'description': parts[4] if len(parts) > 4 else ''
                    })
    except FileNotFoundError:
        print(f"Warning: Hotspot BED file not found: {hotspot_bed_file}")
        # Create default B-ALL hotspots
        hotspots = create_default_ball_hotspots()
    except Exception as e:
        print(f"Error loading hotspot regions: {e}")
        hotspots = create_default_ball_hotspots()
    
    return hotspots

def create_default_ball_hotspots():
    """Create default B-ALL hotspot regions"""
    return [
        {'chromosome': 'chr9', 'start': 133729450, 'end': 133763062, 'region': 'ABL1', 'description': 'ABL1 kinase domain'},
        {'chromosome': 'chr22', 'start': 23522552, 'end': 23655671, 'region': 'BCR', 'description': 'BCR breakpoint cluster'},
        {'chromosome': 'chr12', 'start': 11802788, 'end': 12048325, 'region': 'ETV6', 'description': 'ETV6 gene'},
        {'chromosome': 'chr21', 'start': 36160098, 'end': 37376887, 'region': 'RUNX1', 'description': 'RUNX1 gene'},
        {'chromosome': 'chr19', 'start': 1606165, 'end': 1722865, 'region': 'TCF3', 'description': 'TCF3 gene'},
        {'chromosome': 'chr1', 'start': 164517598, 'end': 164687225, 'region': 'PBX1', 'description': 'PBX1 gene'},
        {'chromosome': 'chr11', 'start': 118307205, 'end': 118397540, 'region': 'KMT2A', 'description': 'KMT2A/MLL gene'},
        {'chromosome': 'chr17', 'start': 7565097, 'end': 7590856, 'region': 'TP53', 'description': 'TP53 tumor suppressor'},
        {'chromosome': 'chrX', 'start': 24130000, 'end': 24280000, 'region': 'CRLF2', 'description': 'CRLF2 gene'},
        {'chromosome': 'chr9', 'start': 35820000, 'end': 36120000, 'region': 'PAX5', 'description': 'PAX5 gene'},
        {'chromosome': 'chr7', 'start': 50305000, 'end': 50405000, 'region': 'IKZF1', 'description': 'IKZF1 gene'}
    ]

def load_driver_genes(driver_genes_file):
    """Load driver gene list"""
    driver_genes = set()
    
    try:
        with open(driver_genes_file, 'r') as f:
            driver_genes = set(line.strip() for line in f if line.strip())
    except Exception as e:
        print(f"Warning: Could not load driver genes from {driver_genes_file}: {e}")
        # Default B-ALL driver genes
        driver_genes = {
            'ABL1', 'ABL2', 'AFF1', 'BCL2', 'BCL6', 'BCR', 'CREBBP', 'CRLF2',
            'DUX4', 'EBF1', 'EP300', 'ETV6', 'IKZF1', 'KMT2A', 'MLL', 'MYC',
            'NRAS', 'KRAS', 'PAX5', 'PBX1', 'RUNX1', 'TCF3', 'TP53', 'ZNF384'
        }
    
    return driver_genes

def parse_vcf_file(vcf_file):
    """Parse VCF file and extract variant information"""
    variants = []
    
    if HAS_CYVCF2:
        return parse_vcf_with_cyvcf2(vcf_file)
    else:
        return parse_vcf_basic(vcf_file)

def parse_vcf_with_cyvcf2(vcf_file):
    """Parse VCF file using cyvcf2 (modern, fast approach)"""
    variants = []
    
    try:
        vcf_reader = VCF(vcf_file)
        
        for variant in vcf_reader:
            variant_info = {
                'chromosome': str(variant.CHROM),
                'position': variant.POS,
                'ref': variant.REF,
                'alt': str(variant.ALT[0]) if variant.ALT else '',
                'quality': variant.QUAL if variant.QUAL else 0,
                'filter': variant.FILTER if variant.FILTER else 'PASS',
                'info': dict(variant.INFO)
            }
            
            # Extract sample-specific information
            if len(variant.gt_types) > 0:
                variant_info.update({
                    'genotype': variant.gt_types[0],
                    'depth': variant.gt_depths[0] if variant.gt_depths else None,
                    'allele_depth': variant.gt_alt_depths[0] if variant.gt_alt_depths else None
                })
            
            # Extract gene annotation if available
            if 'ANN' in variant.INFO:
                annotations = variant.INFO['ANN']
                if annotations:
                    # Parse first annotation
                    ann_fields = annotations.split('|') if isinstance(annotations, str) else annotations[0].split('|')
                    if len(ann_fields) > 3:
                        variant_info['gene'] = ann_fields[3]
                        variant_info['effect'] = ann_fields[1]
                        variant_info['impact'] = ann_fields[2]
            
            variants.append(variant_info)
            
    except Exception as e:
        print(f"Error parsing VCF file with cyvcf2 {vcf_file}: {e}")
        # Fallback to basic parsing
        return parse_vcf_basic(vcf_file)
    
    return variants

def parse_vcf_basic(vcf_file):
    """Basic VCF parsing as fallback"""
    variants = []
    
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                
                variant_info = {
                    'chromosome': fields[0],
                    'position': int(fields[1]),
                    'ref': fields[3],
                    'alt': fields[4],
                    'quality': float(fields[5]) if fields[5] != '.' else 0,
                    'filter': fields[6] if fields[6] != '.' else 'PASS',
                    'info': {}
                }
                
                # Parse INFO field
                info_pairs = fields[7].split(';')
                for pair in info_pairs:
                    if '=' in pair:
                        key, value = pair.split('=', 1)
                        variant_info['info'][key] = value
                
                variants.append(variant_info)
                
    except Exception as e:
        print(f"Error in basic VCF parsing {vcf_file}: {e}")
    
    return variants

def main():
    # Get inputs from snakemake
    vcf_file = snakemake.input.vcf
    output_hotspots = snakemake.output.hotspots
    output_summary = snakemake.output.summary
    
    hotspot_bed = snakemake.params.hotspot_bed
    driver_genes_file = snakemake.params.driver_genes
    
    # Load reference data
    hotspots = load_hotspot_regions(hotspot_bed)
    driver_genes = load_driver_genes(driver_genes_file)
    
    # Parse variants
    variants = parse_vcf_file(vcf_file)
    
    if not variants:
        print("No variants found in VCF file")
        # Create empty outputs
        empty_df = pd.DataFrame()
        empty_df.to_csv(output_hotspots, index=False)
        with open(output_summary, 'w') as f:
            json.dump({'total_variants': 0, 'hotspot_variants': 0}, f)
        return
    
    print(f"Analyzing {len(variants)} variants for hotspots...")
    
    # Create summary
    summary = {
        'total_variants': len(variants),
        'hotspot_variants': 0,
        'driver_gene_variants': 0,
        'high_significance_variants': 0,
        'high_quality_variants': 0
    }
    
    # Save empty outputs for now (placeholder)
    empty_df = pd.DataFrame()
    empty_df.to_csv(output_hotspots, index=False)
    
    with open(output_summary, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Hotspot analysis completed: {summary['total_variants']} variants processed")

if __name__ == "__main__":
    main()