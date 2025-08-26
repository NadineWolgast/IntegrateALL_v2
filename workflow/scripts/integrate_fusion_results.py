#!/usr/bin/env python3
"""
Fusion Integration Script for IntegrateALL Pipeline
==================================================

Integrates results from Arriba and FusionCatcher to create a unified fusion callset.
This script replaces the original bash processing with a unified Python approach.
"""

import pandas as pd
import numpy as np
import json
import logging
from pathlib import Path
from typing import List, Dict, Set, Tuple
import argparse
import re

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FusionIntegrator:
    """Integrates fusion calls from multiple tools"""
    
    def __init__(self, min_support: int = 3, driver_genes_file: str = ""):
        """
        Initialize fusion integrator
        
        Args:
            min_support: Minimum number of supporting reads for a fusion
            driver_genes_file: Path to file containing driver genes (optional)
        """
        self.min_support = min_support
        self.driver_genes = self._load_driver_genes(driver_genes_file)
        
        # Known recurrent fusion patterns in B-ALL
        self.known_ball_fusions = {
            'BCR-ABL1': {'priority': 1, 'subtype': 'Ph-positive'},
            'ETV6-RUNX1': {'priority': 1, 'subtype': 'ETV6-RUNX1'},
            'TCF3-PBX1': {'priority': 1, 'subtype': 'TCF3-PBX1'},
            'KMT2A': {'priority': 1, 'subtype': 'KMT2A-rearranged'},  # Any KMT2A fusion
            'DUX4': {'priority': 2, 'subtype': 'DUX4-rearranged'},
            'ZNF384': {'priority': 2, 'subtype': 'ZNF384-rearranged'},
            'PAX5': {'priority': 3, 'subtype': 'PAX5alt'},
            'CRLF2': {'priority': 3, 'subtype': 'Ph-like'},
            'IKZF1': {'priority': 3, 'subtype': 'Ph-like'},
        }
        
        # Fusion filtering criteria
        self.filter_criteria = {
            'min_spanning_reads': 3,
            'min_split_reads': 2,
            'max_distance_same_chr': 100000,  # Filter out very close fusions on same chromosome
            'blacklist_patterns': ['RP11-', 'CTD-', 'AC0', 'AL', 'LINC']  # Common false positives
        }
    
    def _load_driver_genes(self, driver_genes_file: str) -> Set[str]:
        """Load driver genes from file or use default list"""
        
        default_drivers = {
            'ABL1', 'ABL2', 'AFF1', 'BCL2', 'BCL6', 'BCR', 'CREBBP', 'CRLF2', 
            'DUX4', 'EBF1', 'EP300', 'ETV6', 'IKZF1', 'KMT2A', 'MLL', 'MYC', 
            'NRAS', 'KRAS', 'PAX5', 'PBX1', 'RUNX1', 'TCF3', 'TP53', 'ZNF384'
        }
        
        if driver_genes_file and Path(driver_genes_file).exists():
            try:
                with open(driver_genes_file, 'r') as f:
                    loaded_genes = {line.strip() for line in f if line.strip()}
                logger.info(f"Loaded {len(loaded_genes)} driver genes from file")
                return loaded_genes
            except Exception as e:
                logger.warning(f"Error loading driver genes file: {e}. Using defaults.")
        
        logger.info(f"Using default driver gene list ({len(default_drivers)} genes)")
        return default_drivers
    
    def parse_arriba_results(self, arriba_file: str) -> pd.DataFrame:
        """Parse Arriba fusion results"""
        try:
            # Check if file exists and has content
            if not Path(arriba_file).exists() or Path(arriba_file).stat().st_size == 0:
                logger.warning(f"Arriba file is empty or missing: {arriba_file}")
                return pd.DataFrame()
            
            # Read Arriba results
            arriba_df = pd.read_csv(arriba_file, sep='\\t', low_memory=False)
            
            if arriba_df.empty:
                logger.info("No Arriba fusions found")
                return pd.DataFrame()
            
            # Standardize column names and extract key information
            standardized = []
            for _, row in arriba_df.iterrows():
                
                # Parse gene names (handle complex gene names)
                gene1 = str(row.get('#gene1', row.get('gene1', ''))).split(',')[0]
                gene2 = str(row.get('gene2', '')).split(',')[0]
                
                # Skip if essential information is missing
                if not gene1 or not gene2 or gene1 == 'nan' or gene2 == 'nan':
                    continue
                
                # Extract supporting evidence
                split_reads1 = self._safe_int(row.get('split_reads1', 0))
                split_reads2 = self._safe_int(row.get('split_reads2', 0))
                discordant_mates = self._safe_int(row.get('discordant_mates', 0))
                
                total_support = split_reads1 + split_reads2 + discordant_mates
                
                # Extract confidence metrics
                confidence = row.get('confidence', 'medium')
                
                standardized.append({
                    'gene1': gene1,
                    'gene2': gene2,
                    'fusion_name': f"{gene1}-{gene2}",
                    'chr1': row.get('breakpoint1', '').split(':')[0] if ':' in str(row.get('breakpoint1', '')) else '',
                    'pos1': self._extract_position(row.get('breakpoint1', '')),
                    'chr2': row.get('breakpoint2', '').split(':')[0] if ':' in str(row.get('breakpoint2', '')) else '',
                    'pos2': self._extract_position(row.get('breakpoint2', '')),
                    'split_reads1': split_reads1,
                    'split_reads2': split_reads2,
                    'discordant_mates': discordant_mates,
                    'total_support': total_support,
                    'confidence': confidence,
                    'caller': 'arriba',
                    'reading_frame': row.get('reading_frame', ''),
                    'type': row.get('type', ''),
                    'direction1': row.get('direction1', ''),
                    'direction2': row.get('direction2', '')
                })
            
            result_df = pd.DataFrame(standardized)
            logger.info(f"Parsed {len(result_df)} Arriba fusions")
            return result_df
            
        except Exception as e:
            logger.error(f"Error parsing Arriba results: {e}")
            return pd.DataFrame()
    
    def parse_fusioncatcher_results(self, fusioncatcher_file: str) -> pd.DataFrame:
        """Parse FusionCatcher results"""
        try:
            # Check if file exists and has content  
            if not Path(fusioncatcher_file).exists() or Path(fusioncatcher_file).stat().st_size == 0:
                logger.warning(f"FusionCatcher file is empty or missing: {fusioncatcher_file}")
                return pd.DataFrame()
            
            # Read FusionCatcher results
            fc_df = pd.read_csv(fusioncatcher_file, sep='\\t', low_memory=False)
            
            if fc_df.empty:
                logger.info("No FusionCatcher fusions found")
                return pd.DataFrame()
            
            # Standardize column names and extract key information
            standardized = []
            for _, row in fc_df.iterrows():
                
                # Extract gene names
                gene1 = str(row.get('Gene_1_symbol(5end_fusion_partner)', '')).strip()
                gene2 = str(row.get('Gene_2_symbol(3end_fusion_partner)', '')).strip()
                
                # Skip if essential information is missing
                if not gene1 or not gene2 or gene1 == 'nan' or gene2 == 'nan':
                    continue
                
                # Extract supporting evidence
                spanning_pairs = self._safe_int(row.get('Spanning_pairs', 0))
                spanning_unique = self._safe_int(row.get('Spanning_unique_reads', 0))
                common_mapping = self._safe_int(row.get('Counts_of_common_mapping_reads', 0))
                
                total_support = spanning_pairs + spanning_unique
                
                # Extract positions
                fusion_point1 = row.get('Fusion_point_for_gene_1(5end_fusion_partner)', '')
                fusion_point2 = row.get('Fusion_point_for_gene_2(3end_fusion_partner)', '')
                
                standardized.append({
                    'gene1': gene1,
                    'gene2': gene2,
                    'fusion_name': f"{gene1}-{gene2}",
                    'chr1': fusion_point1.split(':')[0] if ':' in fusion_point1 else '',
                    'pos1': self._extract_position(fusion_point1),
                    'chr2': fusion_point2.split(':')[0] if ':' in fusion_point2 else '',
                    'pos2': self._extract_position(fusion_point2),
                    'spanning_pairs': spanning_pairs,
                    'spanning_unique_reads': spanning_unique,
                    'common_mapping_reads': common_mapping,
                    'total_support': total_support,
                    'caller': 'fusioncatcher',
                    'fusion_description': row.get('Fusion_description', ''),
                    'longest_anchor': row.get('Longest_anchor_found', ''),
                    'method': row.get('Fusion_finding_method', '')
                })
            
            result_df = pd.DataFrame(standardized)
            logger.info(f"Parsed {len(result_df)} FusionCatcher fusions")
            return result_df
            
        except Exception as e:
            logger.error(f"Error parsing FusionCatcher results: {e}")
            return pd.DataFrame()
    
    def _safe_int(self, value) -> int:
        """Safely convert value to integer"""
        try:
            if pd.isna(value) or value == '' or value == 'nan':
                return 0
            return int(float(str(value)))
        except (ValueError, TypeError):
            return 0
    
    def _extract_position(self, position_str: str) -> int:
        """Extract genomic position from position string"""
        try:
            if ':' in str(position_str):
                return int(position_str.split(':')[1])
            else:
                return 0
        except (ValueError, TypeError, IndexError):
            return 0
    
    def filter_fusions(self, fusions_df: pd.DataFrame) -> pd.DataFrame:
        """Apply quality filters to fusions"""
        if fusions_df.empty:
            return fusions_df
        
        original_count = len(fusions_df)
        
        # Filter 1: Minimum support reads
        fusions_df = fusions_df[fusions_df['total_support'] >= self.min_support]
        
        # Filter 2: Remove likely false positives based on gene names
        pattern = '|'.join(self.filter_criteria['blacklist_patterns'])
        fusions_df = fusions_df[
            ~fusions_df['gene1'].str.contains(pattern, case=False, na=False) &
            ~fusions_df['gene2'].str.contains(pattern, case=False, na=False)
        ]
        
        # Filter 3: Remove very close fusions on the same chromosome (likely artifacts)
        same_chr_mask = (fusions_df['chr1'] == fusions_df['chr2']) & (fusions_df['chr1'] != '')
        if same_chr_mask.any():
            close_fusions = same_chr_mask & (
                abs(fusions_df['pos1'] - fusions_df['pos2']) < self.filter_criteria['max_distance_same_chr']
            )
            fusions_df = fusions_df[~close_fusions]
        
        # Filter 4: Remove self-fusions (same gene)
        fusions_df = fusions_df[fusions_df['gene1'] != fusions_df['gene2']]
        
        filtered_count = len(fusions_df)
        logger.info(f"Filtered fusions: {original_count} -> {filtered_count}")
        
        return fusions_df
    
    def integrate_fusion_calls(self, arriba_df: pd.DataFrame, fusioncatcher_df: pd.DataFrame) -> pd.DataFrame:
        """Integrate fusion calls from both tools"""
        
        all_fusions = []
        
        # Add Arriba fusions
        for _, fusion in arriba_df.iterrows():
            all_fusions.append(fusion.to_dict())
        
        # Add FusionCatcher fusions
        for _, fusion in fusioncatcher_df.iterrows():
            all_fusions.append(fusion.to_dict())
        
        if not all_fusions:
            logger.info("No fusions to integrate")
            return pd.DataFrame()
        
        # Convert to DataFrame for easier processing
        integrated_df = pd.DataFrame(all_fusions)
        
        # Group similar fusions (same gene pair)
        fusion_groups = {}
        for _, fusion in integrated_df.iterrows():
            # Create canonical fusion name (alphabetically sorted)
            genes = sorted([fusion['gene1'], fusion['gene2']])
            canonical_name = f"{genes[0]}-{genes[1]}"
            
            if canonical_name not in fusion_groups:
                fusion_groups[canonical_name] = []
            fusion_groups[canonical_name].append(fusion.to_dict())
        
        # Merge evidence for each fusion group
        merged_fusions = []
        for canonical_name, group in fusion_groups.items():
            merged_fusion = self._merge_fusion_group(canonical_name, group)
            merged_fusions.append(merged_fusion)
        
        result_df = pd.DataFrame(merged_fusions)
        logger.info(f"Integrated into {len(result_df)} unique fusions")
        
        return result_df
    
    def _merge_fusion_group(self, canonical_name: str, fusion_group: List[Dict]) -> Dict:
        """Merge evidence for fusions targeting the same gene pair"""
        
        # Find the best representative fusion
        best_fusion = max(fusion_group, key=lambda f: f['total_support'])
        
        # Count callers
        callers = list(set(f['caller'] for f in fusion_group))
        caller_count = len(callers)
        
        # Sum support across callers
        total_support = sum(f['total_support'] for f in fusion_group)
        
        # Determine confidence based on caller agreement and support
        if caller_count >= 2:
            confidence = 'high'
        elif total_support >= 10:
            confidence = 'medium'
        else:
            confidence = 'low'
        
        # Check if involves driver genes
        genes = canonical_name.split('-')
        is_driver_fusion = any(gene in self.driver_genes for gene in genes)
        
        # Check for known B-ALL fusions
        known_fusion_type = self._classify_known_fusion(canonical_name, genes)
        
        merged = best_fusion.copy()
        merged.update({
            'canonical_name': canonical_name,
            'callers': ','.join(callers),
            'caller_count': caller_count,
            'total_support_merged': total_support,
            'confidence_merged': confidence,
            'is_driver_fusion': is_driver_fusion,
            'known_fusion_type': known_fusion_type,
            'priority': self._get_fusion_priority(canonical_name, genes)
        })
        
        return merged
    
    def _classify_known_fusion(self, canonical_name: str, genes: List[str]) -> str:
        """Classify fusion based on known B-ALL patterns"""
        
        # Check exact matches first
        if canonical_name in self.known_ball_fusions:
            return self.known_ball_fusions[canonical_name]['subtype']
        
        # Check gene-level matches
        for gene in genes:
            if gene in self.known_ball_fusions:
                return self.known_ball_fusions[gene]['subtype']
        
        return 'unknown'
    
    def _get_fusion_priority(self, canonical_name: str, genes: List[str]) -> int:
        """Get priority score for fusion (lower = higher priority)"""
        
        # Check exact matches
        if canonical_name in self.known_ball_fusions:
            return self.known_ball_fusions[canonical_name]['priority']
        
        # Check gene-level matches
        min_priority = 999
        for gene in genes:
            if gene in self.known_ball_fusions:
                min_priority = min(min_priority, self.known_ball_fusions[gene]['priority'])
        
        if min_priority < 999:
            return min_priority
        
        # Driver gene fusions get medium priority
        if any(gene in self.driver_genes for gene in genes):
            return 5
        
        # All others get low priority
        return 10
    
    def create_summary(self, integrated_df: pd.DataFrame) -> Dict:
        """Create summary of fusion results"""
        
        if integrated_df.empty:
            return {
                'total_fusions': 0,
                'driver_fusions': 0,
                'high_confidence_fusions': 0,
                'known_ball_fusions': 0,
                'fusion_types': {}
            }
        
        summary = {
            'total_fusions': len(integrated_df),
            'driver_fusions': sum(integrated_df['is_driver_fusion']),
            'high_confidence_fusions': sum(integrated_df['confidence_merged'] == 'high'),
            'caller_concordant': sum(integrated_df['caller_count'] >= 2),
            'known_ball_fusions': sum(integrated_df['known_fusion_type'] != 'unknown')
        }
        
        # Count by known fusion types
        fusion_type_counts = integrated_df['known_fusion_type'].value_counts().to_dict()
        summary['fusion_types'] = fusion_type_counts
        
        # Top fusions by support
        top_fusions = integrated_df.nlargest(5, 'total_support_merged')[['canonical_name', 'total_support_merged']].to_dict('records')
        summary['top_fusions'] = top_fusions
        
        return summary

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Integrate fusion detection results')
    parser.add_argument('arriba_file', help='Arriba fusions TSV file')
    parser.add_argument('fusioncatcher_file', help='FusionCatcher results file')
    parser.add_argument('output_integrated', help='Output integrated fusions TSV')
    parser.add_argument('output_summary', help='Output fusion summary JSON')
    parser.add_argument('--min-support', type=int, default=3, help='Minimum supporting reads')
    parser.add_argument('--driver-genes', default='', help='Driver genes file')
    
    args = parser.parse_args()
    
    # Initialize integrator
    integrator = FusionIntegrator(args.min_support, args.driver_genes)
    
    # Parse input files
    logger.info("Parsing Arriba results...")
    arriba_df = integrator.parse_arriba_results(args.arriba_file)
    
    logger.info("Parsing FusionCatcher results...")
    fusioncatcher_df = integrator.parse_fusioncatcher_results(args.fusioncatcher_file)
    
    # Apply quality filters
    logger.info("Applying quality filters...")
    arriba_filtered = integrator.filter_fusions(arriba_df)
    fusioncatcher_filtered = integrator.filter_fusions(fusioncatcher_df)
    
    # Integrate results
    logger.info("Integrating fusion calls...")
    integrated_df = integrator.integrate_fusion_calls(arriba_filtered, fusioncatcher_filtered)
    
    # Create summary
    summary = integrator.create_summary(integrated_df)
    
    # Save results
    logger.info(f"Saving {len(integrated_df)} integrated fusions...")
    
    if not integrated_df.empty:
        # Sort by priority and support
        integrated_df = integrated_df.sort_values(['priority', 'total_support_merged'], ascending=[True, False])
        integrated_df.to_csv(args.output_integrated, sep='\\t', index=False)
    else:
        # Create empty file with header
        empty_df = pd.DataFrame(columns=[
            'gene1', 'gene2', 'fusion_name', 'canonical_name', 'total_support_merged',
            'confidence_merged', 'is_driver_fusion', 'known_fusion_type', 'callers'
        ])
        empty_df.to_csv(args.output_integrated, sep='\\t', index=False)
    
    # Save summary
    with open(args.output_summary, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info("Fusion integration completed successfully")
    logger.info(f"Summary: {summary['total_fusions']} total, {summary['driver_fusions']} driver, {summary['high_confidence_fusions']} high-confidence")

if __name__ == "__main__":
    main()