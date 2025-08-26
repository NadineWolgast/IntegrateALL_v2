#!/usr/bin/env python3
"""
Optimized Final Classification Script for IntegrateALL Pipeline
==============================================================
Based on original make_final_classification_new.py but with performance improvements.
Maintains exact classification logic and output format for research validity.
"""

import csv
import os
import pandas as pd
import numpy as np
from pathlib import Path
import json
from typing import Tuple, List, Dict, Optional

# Original driver fusion gene list - CRITICAL FOR CLASSIFICATION
DRIVER_FUSION_LIST = [
    "ABL1", "ABL2", "ACIN1", "ADAMTSL5", "AFF1", "AMPH", "ANTXR1", "ARID1B", 
    "ATF7IP", "ATXN7L3", "AUTS2", "BCL2", "BCL6", "BCL9", "BCR", "BMP2K", 
    "BRD9", "C7ORF72", "CASC15", "CBFA2T2", "CBFA2T3", "CBL", "CD163", 
    "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CENPC", "CLTC", "CREBBP", "CRLF2", 
    "CSF1R", "CUX1", "DACH1", "DACH2", "DAZAP1", "DBX1", "DCPS", "DDX42", 
    "DGKH", "DMRTA2", "DUX4", "EBF1", "ELMO1", "EP300", "ETV6", "EXOSC2", 
    "FHIT", "FIP1L1", "FOXP1", "GATAD2A", "H3C1", "H3C2", "H3C8", "HHEX", 
    "HOXB4", "ID4", "IGH@", "IKZF1", "IL3", "JAK2", "KAT6A", "KDM6A", 
    "KMT2A", "KMT2D", "KRAS", "LIFR", "LSM14A", "LTK", "LYN", "MED12", 
    "MEF2D", "MEOX1", "MLLT10", "MNX1", "MYC", "NF1", "NRAS", "NSD2", 
    "NUP214", "OTX1", "P2RY8", "PAG1", "PAX5", "PBX1", "PDGFRA", "PDGFRB", 
    "PHF6", "PML", "RANBP2", "RB1", "RCSD1", "RUNX1", "SETD2", "SH2B3", 
    "SNX2", "SSBP2", "TAL1", "TBL1XR1", "TCF3", "TEAD1", "TLX1", "TLX3", 
    "TNIP1", "TP53", "USP7", "VHL", "ZC3HAV1", "ZEB2", "ZNF362", "ZNF384", 
    "ZMIZ1", "ZMYND8"
]

def get_allcatchr_data(allcatchr_file: str) -> Tuple[str, float, str]:
    """Extract ALLCatchR prediction data with error handling."""
    try:
        df = pd.read_csv(allcatchr_file, sep='\t')
        if df.empty:
            return "Unknown", 0.0, "not Ph-pos predicted"
        
        row = df.iloc[0]
        subtype = str(row.get('Prediction', 'Unknown'))
        confidence = float(row.get('Confidence', 0.0))
        
        # Determine Ph-pos status based on subtype
        ph_pos = "lymphoid" if "Ph-pos" in subtype else "not Ph-pos predicted"
        
        return subtype, confidence, ph_pos
    except Exception as e:
        print(f"Warning: Error reading ALLCatchR file {allcatchr_file}: {e}")
        return "Unknown", 0.0, "not Ph-pos predicted"

def get_karyotype_and_probabilities(karyotype_file: str) -> Tuple[str, float]:
    """Extract karyotype prediction with error handling."""
    try:
        if os.path.exists(karyotype_file) and os.path.getsize(karyotype_file) > 0:
            df = pd.read_csv(karyotype_file)
            if not df.empty and 'Output' in df.columns:
                output_text = str(df['Output'].iloc[0])
                # Extract karyotype from output text
                if "Hyperdiploid" in output_text:
                    return "Hyperdiploid", 1.0
                elif "Near haploid" in output_text:
                    return "Near haploid", 1.0
                elif "Low hypodiploid" in output_text:
                    return "Low hypodiploid", 1.0
                elif "Near triploid" in output_text:
                    return "Near triploid", 1.0
        
        return "other", 0.5
    except Exception as e:
        print(f"Warning: Error reading karyotype file {karyotype_file}: {e}")
        return "other", 0.5

def check_hotspot_files(hotspot_dir: str) -> Dict[str, bool]:
    """Check for SNV hotspot files - CRITICAL for classification."""
    relevant_files = {
        "PAX5_P80R": False,
        "IKZF1_N159Y": False,
        "ZEB2_H1038R": False
    }
    
    try:
        if not os.path.exists(hotspot_dir):
            return relevant_files
            
        hotspot_files = os.listdir(hotspot_dir)
        for file in hotspot_files:
            if file.startswith("PAX5_P80R"):
                relevant_files["PAX5_P80R"] = True
            elif file.startswith("IKZF1_N159Y"):
                relevant_files["IKZF1_N159Y"] = True
            elif file.startswith("ZEB2_H1038"):
                relevant_files["ZEB2_H1038R"] = True
    except Exception as e:
        print(f"Warning: Error checking hotspot files: {e}")
    
    return relevant_files

def filter_fusions(fusion_genes: List[Tuple], unique_genes: List[str], 
                  df: pd.DataFrame, subgruppe: str) -> List[Tuple]:
    """Filter fusions based on driver genes and classification rules."""
    filtered_fusions = []
    
    for gene_1, gene_2, caller, unique_spanning_reads in fusion_genes:
        # Handle IGH genes
        if gene_1.startswith("IGH"):
            gene_1 = "IGH@"
        if gene_2.startswith("IGH"):
            gene_2 = "IGH@"
        
        # Check if both genes are in driver list
        if gene_1 in unique_genes and gene_2 in unique_genes:
            # Check if gene pair exists in classification dataframe
            match_df = df[
                (df['Gene_1_symbol(5end_fusion_partner)'] == gene_1) &
                (df['Gene_2_symbol(3end_fusion_partner)'] == gene_2)
            ]
            
            if not match_df.empty and (subgruppe == 'DUX4' or ('DUX4' not in [gene_1, gene_2])):
                filtered_fusions.append((gene_1, gene_2, caller, unique_spanning_reads))
    
    return filtered_fusions

def gather_data(allcatchr_file: str, karyotype_file: str, fusioncatcher_file: str, 
               arriba_file: str, hotspot_dir: str, classification_file: str) -> pd.DataFrame:
    """Gather all classification data - maintains original logic exactly."""
    
    # Extract data from input files
    subgruppe, confidence, bcr_abl1_maincluster_pred = get_allcatchr_data(allcatchr_file)
    karyotype, score = get_karyotype_and_probabilities(karyotype_file)
    relevant_files = check_hotspot_files(hotspot_dir)
    
    fusion_genes = []
    
    # Read FusionCatcher fusions (Performance: use pandas for large files)
    try:
        if os.path.exists(fusioncatcher_file) and os.path.getsize(fusioncatcher_file) > 0:
            fc_df = pd.read_csv(fusioncatcher_file, sep='\t')
            if not fc_df.empty and len(fc_df.columns) >= 6:
                for _, row in fc_df.iterrows():
                    fusion_genes.append((
                        str(row.iloc[0]), str(row.iloc[1]), 
                        'FusionCatcher', str(row.iloc[5])
                    ))
    except Exception as e:
        print(f"Warning: Error reading FusionCatcher file: {e}")
    
    # Read Arriba fusions 
    try:
        if os.path.exists(arriba_file) and os.path.getsize(arriba_file) > 0:
            arriba_df = pd.read_csv(arriba_file, sep='\t')
            if not arriba_df.empty and len(arriba_df.columns) >= 12:
                for _, row in arriba_df.iterrows():
                    fusion_genes.append((
                        str(row.iloc[0]), str(row.iloc[1]),
                        'Arriba', str(row.iloc[11])
                    ))
    except Exception as e:
        print(f"Warning: Error reading Arriba file: {e}")
    
    # Load classification dataframe
    df = pd.read_csv(classification_file, sep=',')
    
    # Filter fusions using original logic
    filtered_fusions = filter_fusions(fusion_genes, DRIVER_FUSION_LIST, df, subgruppe)
    
    # Create result dataframe maintaining exact original structure
    if not filtered_fusions:
        result_df = pd.DataFrame([{
            'ALLCatchR': subgruppe,
            'Ph-pos': bcr_abl1_maincluster_pred,
            'Confidence': confidence,
            'Gene_1_symbol(5end_fusion_partner)': None,
            'Gene_2_symbol(3end_fusion_partner)': None,
            'Fusioncaller': None,
            'Unique_spanning_reads': None,
            'karyotype_classifier': karyotype,
            'PAX5_P80R': relevant_files["PAX5_P80R"],
            'IKZF1_N159Y': relevant_files["IKZF1_N159Y"],
            'ZEB2_H1038R': relevant_files["ZEB2_H1038R"]
        }])
    else:
        # Create dataframe from filtered fusions
        result_df = pd.DataFrame(filtered_fusions, columns=[
            'Gene_1_symbol(5end_fusion_partner)',
            'Gene_2_symbol(3end_fusion_partner)',
            'Fusioncaller',
            'Unique_spanning_reads'
        ])
        
        # Add constant columns
        result_df['ALLCatchR'] = subgruppe
        result_df['Ph-pos'] = bcr_abl1_maincluster_pred
        result_df['Confidence'] = confidence
        result_df['karyotype_classifier'] = karyotype
        result_df['PAX5_P80R'] = relevant_files["PAX5_P80R"]
        result_df['IKZF1_N159Y'] = relevant_files["IKZF1_N159Y"]
        result_df['ZEB2_H1038R'] = relevant_files["ZEB2_H1038R"]
        
        # Reorder columns to match original
        result_df = result_df[[
            'ALLCatchR', 'Ph-pos', 'Confidence',
            'Gene_1_symbol(5end_fusion_partner)',
            'Gene_2_symbol(3end_fusion_partner)',
            'Fusioncaller', 'Unique_spanning_reads',
            'karyotype_classifier', 'PAX5_P80R', 
            'IKZF1_N159Y', 'ZEB2_H1038R'
        ]]
    
    return result_df

def check_conditions(data_df: pd.DataFrame, df_classification: pd.DataFrame) -> List[pd.Series]:
    """Check conditions against classification rules - maintains original logic."""
    comparison_df = df_classification
    matched_rows = []
    
    if not data_df['Gene_1_symbol(5end_fusion_partner)'].isnull().all():
        for _, data_row in data_df.iterrows():
            # Performance: vectorized comparison where possible
            matches = comparison_df[
                (comparison_df['ALLCatchR'] == data_row['ALLCatchR']) &
                (comparison_df['Ph-pos'] == data_row['Ph-pos']) &
                (comparison_df['Gene_1_symbol(5end_fusion_partner)'] == data_row['Gene_1_symbol(5end_fusion_partner)']) &
                (comparison_df['Gene_2_symbol(3end_fusion_partner)'] == data_row['Gene_2_symbol(3end_fusion_partner)']) &
                (comparison_df['karyotype classifier'] == data_row['karyotype_classifier']) &
                (comparison_df['PAX5 P80R'] == data_row['PAX5_P80R']) &
                (comparison_df['IKZF1 N159Y'] == data_row['IKZF1_N159Y']) &
                (comparison_df['ZEB2 H1038R'] == data_row['ZEB2_H1038R'])
            ]
            
            if not matches.empty:
                match_row = matches.iloc[0]
                data_row_copy = data_row.copy()
                data_row_copy['WHO-HAEM5'] = match_row['WHO-HAEM5']
                data_row_copy['ICC'] = match_row['ICC']
                matched_rows.append(data_row_copy)
                break
    else:
        # Handle no fusion case
        for _, data_row in data_df.iterrows():
            matches = comparison_df[
                (comparison_df['ALLCatchR'] == data_row['ALLCatchR']) &
                (comparison_df['Ph-pos'] == data_row['Ph-pos']) &
                (comparison_df['karyotype classifier'] == data_row['karyotype_classifier']) &
                (comparison_df['PAX5 P80R'] == data_row['PAX5_P80R']) &
                (comparison_df['IKZF1 N159Y'] == data_row['IKZF1_N159Y']) &
                (comparison_df['ZEB2 H1038R'] == data_row['ZEB2_H1038R']) &
                (comparison_df['Gene_1_symbol(5end_fusion_partner)'].isnull()) &
                (comparison_df['Gene_2_symbol(3end_fusion_partner)'].isnull())
            ]
            
            if not matches.empty:
                match_row = matches.iloc[0]
                data_row_copy = data_row.copy()
                data_row_copy['WHO-HAEM5'] = match_row['WHO-HAEM5']
                data_row_copy['ICC'] = match_row['ICC']
                matched_rows.append(data_row_copy)
                break
    
    return matched_rows

def main():
    """Main function using snakemake inputs."""
    # Get inputs from snakemake
    allcatchr_file = snakemake.input.allcatchr_pred
    karyotype_file = snakemake.input.karyotype_pred
    fusioncatcher_file = snakemake.input.fusioncatcher_results
    arriba_file = snakemake.input.arriba_results
    hotspot_dir = snakemake.input.hotspot_dir
    classification_file = snakemake.input.classification_rules
    
    # Outputs
    output_report = snakemake.output.final_report
    output_driver = snakemake.output.driver_summary
    output_curation = snakemake.output.curation_summary
    
    sample_id = snakemake.wildcards.sample
    
    print(f"Processing final classification for sample: {sample_id}")
    
    try:
        # Gather all data using original logic
        data_df = gather_data(
            allcatchr_file, karyotype_file, fusioncatcher_file,
            arriba_file, hotspot_dir, classification_file
        )
        
        # Load classification rules
        df_classification = pd.read_csv(classification_file, sep=',')
        
        # Check conditions and get matches
        matched_rows = check_conditions(data_df, df_classification)
        
        # Generate outputs in original format
        if matched_rows:
            final_df = pd.DataFrame(matched_rows)
            
            # Save main report
            final_df.to_csv(output_report, index=False)
            
            # Save driver summary (fusion information)
            driver_df = final_df[[
                'Gene_1_symbol(5end_fusion_partner)',
                'Gene_2_symbol(3end_fusion_partner)',
                'Fusioncaller', 'Unique_spanning_reads'
            ]].copy()
            driver_df.to_csv(output_driver, index=False)
            
            # Save curation summary
            curation_df = final_df[[
                'WHO-HAEM5', 'ICC', 'ALLCatchR', 'Confidence',
                'karyotype_classifier', 'PAX5_P80R', 'IKZF1_N159Y', 'ZEB2_H1038R'
            ]].copy()
            curation_df.to_csv(output_curation, index=False)
            
            print(f"Classification completed: {final_df['WHO-HAEM5'].iloc[0]}")
        else:
            # No match found - create empty outputs
            print("No classification match found")
            empty_df = pd.DataFrame()
            empty_df.to_csv(output_report, index=False)
            empty_df.to_csv(output_driver, index=False) 
            empty_df.to_csv(output_curation, index=False)
            
    except Exception as e:
        print(f"Error in final classification: {e}")
        # Create empty outputs on error
        empty_df = pd.DataFrame()
        empty_df.to_csv(output_report, index=False)
        empty_df.to_csv(output_driver, index=False)
        empty_df.to_csv(output_curation, index=False)

if __name__ == "__main__":
    main()