#!/usr/bin/env python3
"""
Integrated Classification Script for IntegrateALL Pipeline
=========================================================

Combines results from ALLCatchR, karyotype prediction, fusion detection, 
CNV analysis, and hotspot mutations to provide comprehensive B-ALL classification.

This script replaces the original bash/R scripts with a unified Python approach.
"""

import json
import pandas as pd
import numpy as np
import yaml
import logging
from pathlib import Path
from typing import Dict, List, Any, Tuple
import argparse

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class IntegratedClassifier:
    """Integrated B-ALL classifier combining multiple data sources"""
    
    def __init__(self, classification_rules_file: str, confidence_threshold: float = 0.7):
        """
        Initialize the classifier
        
        Args:
            classification_rules_file: Path to YAML file with classification rules
            confidence_threshold: Minimum confidence threshold for classification
        """
        self.confidence_threshold = confidence_threshold
        self.classification_rules = self._load_classification_rules(classification_rules_file)
        
        # Driver gene lists (expanded from original pipeline)
        self.driver_genes = {
            'high_confidence': [
                'BCR-ABL1', 'ETV6-RUNX1', 'TCF3-PBX1', 'KMT2A', 'IKZF1',
                'PAX5', 'CRLF2', 'NRAS', 'KRAS', 'FLT3', 'TP53'
            ],
            'medium_confidence': [
                'BCL2', 'BCL6', 'MYC', 'CREBBP', 'EP300', 'ATM', 'BTK'
            ]
        }
        
        # B-ALL subtypes with their defining characteristics
        self.ball_subtypes = {
            'Ph-positive': {'markers': ['BCR-ABL1'], 'confidence_weight': 1.0},
            'ETV6-RUNX1': {'markers': ['ETV6-RUNX1'], 'confidence_weight': 1.0},
            'TCF3-PBX1': {'markers': ['TCF3-PBX1'], 'confidence_weight': 1.0},
            'KMT2A-rearranged': {'markers': ['KMT2A'], 'confidence_weight': 1.0},
            'Hyperdiploid': {'karyotype': '>50', 'confidence_weight': 0.9},
            'Low-hypodiploid': {'karyotype': '<44', 'confidence_weight': 0.9},
            'Ph-like': {'expression_signature': True, 'confidence_weight': 0.8},
            'PAX5alt': {'markers': ['PAX5'], 'confidence_weight': 0.8},
            'DUX4-rearranged': {'markers': ['DUX4'], 'confidence_weight': 0.8},
            'ZNF384-rearranged': {'markers': ['ZNF384'], 'confidence_weight': 0.8}
        }
    
    def _load_classification_rules(self, rules_file: str) -> Dict:
        """Load classification rules from YAML file"""
        try:
            if Path(rules_file).exists():
                with open(rules_file, 'r') as f:
                    return yaml.safe_load(f)
            else:
                logger.warning(f"Classification rules file not found: {rules_file}")
                return self._get_default_rules()
        except Exception as e:
            logger.error(f"Error loading classification rules: {e}")
            return self._get_default_rules()
    
    def _get_default_rules(self) -> Dict:
        """Provide default classification rules"""
        return {
            'fusion_priority': ['BCR-ABL1', 'ETV6-RUNX1', 'TCF3-PBX1', 'KMT2A'],
            'karyotype_thresholds': {
                'hyperdiploid': {'min_chromosomes': 50, 'max_chromosomes': 67},
                'low_hypodiploid': {'min_chromosomes': 31, 'max_chromosomes': 43},
                'near_haploid': {'min_chromosomes': 24, 'max_chromosomes': 30}
            },
            'confidence_weights': {
                'fusion': 1.0,
                'allcatchr': 0.9,
                'karyotype': 0.8,
                'hotspots': 0.7,
                'cnv': 0.6
            }
        }
    
    def load_allcatchr_results(self, predictions_file: str, scores_file: str) -> Dict:
        """Load and parse ALLCatchR results"""
        try:
            predictions = pd.read_csv(predictions_file, sep='\\t')
            
            # Handle scores file (may be empty)
            scores = {}
            if Path(scores_file).stat().st_size > 0:
                scores_df = pd.read_csv(scores_file, sep='\\t')
                if not scores_df.empty:
                    scores = scores_df.iloc[0].to_dict()
            
            result = {
                'predicted_subtype': predictions.iloc[0]['predicted_subtype'] if not predictions.empty else 'Unknown',
                'confidence': predictions.iloc[0].get('confidence', 0.0) if not predictions.empty else 0.0,
                'scores': scores
            }
            
            logger.info(f"ALLCatchR prediction: {result['predicted_subtype']} (confidence: {result['confidence']:.3f})")
            return result
            
        except Exception as e:
            logger.error(f"Error loading ALLCatchR results: {e}")
            return {'predicted_subtype': 'Unknown', 'confidence': 0.0, 'scores': {}}
    
    def load_karyotype_results(self, karyotype_file: str) -> Dict:
        """Load and parse karyotype prediction results"""
        try:
            karyotype_df = pd.read_csv(karyotype_file)
            
            if karyotype_df.empty:
                return {'karyotype': 'Unknown', 'confidence': 0.0, 'chromosome_count': None}
            
            result = {
                'karyotype': karyotype_df.iloc[0].get('predicted_karyotype', 'Unknown'),
                'confidence': karyotype_df.iloc[0].get('confidence', 0.0),
                'chromosome_count': karyotype_df.iloc[0].get('chromosome_count', None)
            }
            
            logger.info(f"Karyotype prediction: {result['karyotype']} (confidence: {result['confidence']:.3f})")
            return result
            
        except Exception as e:
            logger.error(f"Error loading karyotype results: {e}")
            return {'karyotype': 'Unknown', 'confidence': 0.0, 'chromosome_count': None}
    
    def load_fusion_results(self, fusion_file: str) -> List[Dict]:
        """Load and parse driver fusion results"""
        try:
            fusions_df = pd.read_csv(fusion_file, sep='\\t')
            
            if fusions_df.empty:
                return []
            
            # Filter for high-confidence driver fusions
            driver_fusions = []
            for _, fusion in fusions_df.iterrows():
                fusion_name = f"{fusion.get('gene1', '')}-{fusion.get('gene2', '')}"
                
                # Check if this is a known driver fusion
                is_driver = any(driver in fusion_name for driver in self.driver_genes['high_confidence'])
                
                if is_driver:
                    driver_fusions.append({
                        'fusion': fusion_name,
                        'gene1': fusion.get('gene1', ''),
                        'gene2': fusion.get('gene2', ''),
                        'support_reads': fusion.get('support_reads', 0),
                        'confidence': fusion.get('confidence', 0.0)
                    })
            
            logger.info(f"Found {len(driver_fusions)} driver fusions")
            return driver_fusions
            
        except Exception as e:
            logger.error(f"Error loading fusion results: {e}")
            return []
    
    def load_cnv_results(self, cnv_file: str) -> Dict:
        """Load and parse CNV analysis results"""
        try:
            cnv_df = pd.read_csv(cnv_file, sep='\\t')
            
            if cnv_df.empty:
                return {'alterations': 0, 'gains': 0, 'losses': 0, 'confidence': 0.0}
            
            result = {
                'total_alterations': cnv_df.iloc[0].get('total_alterations', 0),
                'gains': cnv_df.iloc[0].get('gains', 0),
                'losses': cnv_df.iloc[0].get('losses', 0),
                'confidence': min(1.0, cnv_df.iloc[0].get('total_alterations', 0) / 10.0)  # Normalize
            }
            
            logger.info(f"CNV analysis: {result['total_alterations']} alterations")
            return result
            
        except Exception as e:
            logger.error(f"Error loading CNV results: {e}")
            return {'total_alterations': 0, 'gains': 0, 'losses': 0, 'confidence': 0.0}
    
    def load_hotspot_results(self, hotspot_file: str) -> List[Dict]:
        """Load and parse hotspot mutation results"""
        try:
            hotspots_df = pd.read_csv(hotspot_file)
            
            if hotspots_df.empty:
                return []
            
            hotspot_mutations = []
            for _, mutation in hotspots_df.iterrows():
                hotspot_mutations.append({
                    'gene': mutation.get('gene', ''),
                    'mutation': mutation.get('mutation', ''),
                    'vaf': mutation.get('vaf', 0.0),
                    'confidence': mutation.get('confidence', 0.0)
                })
            
            logger.info(f"Found {len(hotspot_mutations)} hotspot mutations")
            return hotspot_mutations
            
        except Exception as e:
            logger.error(f"Error loading hotspot results: {e}")
            return []
    
    def classify_sample(self, allcatchr_results: Dict, karyotype_results: Dict, 
                       fusion_results: List[Dict], cnv_results: Dict, 
                       hotspot_results: List[Dict]) -> Tuple[str, float, Dict]:
        """
        Perform integrated classification
        
        Returns:
            Tuple of (classification, confidence, evidence_dict)
        """
        
        evidence = {
            'allcatchr': allcatchr_results,
            'karyotype': karyotype_results, 
            'fusions': fusion_results,
            'cnv': cnv_results,
            'hotspots': hotspot_results
        }
        
        # Step 1: Check for defining fusions (highest priority)
        if fusion_results:
            for fusion in fusion_results:
                fusion_name = fusion['fusion']
                
                # Check for specific defining fusions
                if 'BCR-ABL1' in fusion_name:
                    return 'Ph-positive', 0.95, evidence
                elif 'ETV6-RUNX1' in fusion_name:
                    return 'ETV6-RUNX1', 0.95, evidence
                elif 'TCF3-PBX1' in fusion_name:
                    return 'TCF3-PBX1', 0.95, evidence
                elif 'KMT2A' in fusion_name:
                    return 'KMT2A-rearranged', 0.90, evidence
                elif 'DUX4' in fusion_name:
                    return 'DUX4-rearranged', 0.85, evidence
                elif 'ZNF384' in fusion_name:
                    return 'ZNF384-rearranged', 0.85, evidence
        
        # Step 2: Check karyotype-based classifications
        if karyotype_results.get('chromosome_count'):
            chrom_count = karyotype_results['chromosome_count']
            karyotype_conf = karyotype_results.get('confidence', 0.0)
            
            if chrom_count >= 50 and chrom_count <= 67:
                return 'Hyperdiploid', 0.85 * karyotype_conf, evidence
            elif chrom_count >= 31 and chrom_count <= 43:
                return 'Low-hypodiploid', 0.80 * karyotype_conf, evidence
            elif chrom_count >= 24 and chrom_count <= 30:
                return 'Near-haploid', 0.80 * karyotype_conf, evidence
        
        # Step 3: Use ALLCatchR prediction if confident enough
        allcatchr_subtype = allcatchr_results.get('predicted_subtype', 'Unknown')
        allcatchr_conf = allcatchr_results.get('confidence', 0.0)
        
        if allcatchr_conf >= self.confidence_threshold and allcatchr_subtype != 'Unknown':
            return allcatchr_subtype, allcatchr_conf, evidence
        
        # Step 4: Check for Ph-like based on gene expression signature and other markers
        if self._is_ph_like_candidate(allcatchr_results, fusion_results, hotspot_results):
            confidence = 0.75  # Lower confidence for Ph-like without definitive markers
            return 'Ph-like', confidence, evidence
        
        # Step 5: Default classification
        if allcatchr_subtype != 'Unknown':
            return allcatchr_subtype, max(0.5, allcatchr_conf), evidence
        else:
            return 'B-other', 0.3, evidence
    
    def _is_ph_like_candidate(self, allcatchr_results: Dict, fusion_results: List[Dict], 
                             hotspot_results: List[Dict]) -> bool:
        """Determine if sample is Ph-like candidate based on available evidence"""
        
        # Check for Ph-like associated alterations
        ph_like_markers = ['CRLF2', 'IKZF1', 'PAX5', 'NRAS', 'KRAS', 'FLT3']
        
        # Check fusions
        for fusion in fusion_results:
            if any(marker in fusion['fusion'] for marker in ph_like_markers):
                return True
        
        # Check hotspot mutations
        for hotspot in hotspot_results:
            if hotspot['gene'] in ph_like_markers:
                return True
        
        # Check ALLCatchR scores for Ph-like signature
        scores = allcatchr_results.get('scores', {})
        ph_like_score = scores.get('Ph-like', 0.0)
        if ph_like_score > 0.6:
            return True
        
        return False
    
    def calculate_overall_confidence(self, classification: str, individual_confidences: Dict) -> float:
        """Calculate overall confidence score based on evidence sources"""
        
        weights = self.classification_rules.get('confidence_weights', {
            'fusion': 1.0, 'allcatchr': 0.9, 'karyotype': 0.8, 
            'hotspots': 0.7, 'cnv': 0.6
        })
        
        total_weight = 0
        weighted_confidence = 0
        
        for source, confidence in individual_confidences.items():
            if source in weights and confidence > 0:
                weight = weights[source]
                weighted_confidence += confidence * weight
                total_weight += weight
        
        if total_weight > 0:
            return min(1.0, weighted_confidence / total_weight)
        else:
            return 0.0

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Integrated B-ALL classification')
    parser.add_argument('allcatchr_pred', help='ALLCatchR predictions file')
    parser.add_argument('allcatchr_scores', help='ALLCatchR scores file') 
    parser.add_argument('karyotype_pred', help='Karyotype predictions file')
    parser.add_argument('fusion_results', help='Driver fusions file')
    parser.add_argument('cnv_results', help='CNV results file')
    parser.add_argument('hotspot_results', help='Hotspot mutations file')
    parser.add_argument('output_classification', help='Output classification JSON')
    parser.add_argument('output_confidence', help='Output confidence JSON')
    parser.add_argument('output_evidence', help='Output evidence TSV')
    parser.add_argument('--rules', default='', help='Classification rules YAML file')
    parser.add_argument('--confidence-threshold', type=float, default=0.7, help='Confidence threshold')
    
    args = parser.parse_args()
    
    # Initialize classifier
    classifier = IntegratedClassifier(args.rules, args.confidence_threshold)
    
    # Load all input data
    logger.info("Loading input data...")
    allcatchr_results = classifier.load_allcatchr_results(args.allcatchr_pred, args.allcatchr_scores)
    karyotype_results = classifier.load_karyotype_results(args.karyotype_pred)
    fusion_results = classifier.load_fusion_results(args.fusion_results)
    cnv_results = classifier.load_cnv_results(args.cnv_results)
    hotspot_results = classifier.load_hotspot_results(args.hotspot_results)
    
    # Perform classification
    logger.info("Performing integrated classification...")
    classification, confidence, evidence = classifier.classify_sample(
        allcatchr_results, karyotype_results, fusion_results, cnv_results, hotspot_results
    )
    
    # Calculate overall confidence
    individual_confidences = {
        'allcatchr': allcatchr_results.get('confidence', 0.0),
        'karyotype': karyotype_results.get('confidence', 0.0),
        'fusion': max([f.get('confidence', 0.0) for f in fusion_results], default=0.0),
        'cnv': cnv_results.get('confidence', 0.0),
        'hotspots': max([h.get('confidence', 0.0) for h in hotspot_results], default=0.0)
    }
    
    overall_confidence = classifier.calculate_overall_confidence(classification, individual_confidences)
    
    # Prepare outputs
    classification_result = {
        'sample_id': Path(args.output_classification).stem.split('_')[0],
        'final_classification': classification,
        'classification_confidence': confidence,
        'overall_confidence': overall_confidence,
        'classification_method': 'integrated',
        'evidence_sources': list(individual_confidences.keys())
    }
    
    confidence_result = {
        'overall_confidence': overall_confidence,
        'individual_confidences': individual_confidences,
        'confidence_threshold': args.confidence_threshold,
        'meets_threshold': overall_confidence >= args.confidence_threshold
    }
    
    # Save results
    logger.info(f"Saving results: {classification} (confidence: {overall_confidence:.3f})")
    
    with open(args.output_classification, 'w') as f:
        json.dump(classification_result, f, indent=2)
    
    with open(args.output_confidence, 'w') as f:
        json.dump(confidence_result, f, indent=2)
    
    # Save evidence as TSV
    evidence_rows = []
    for source, data in evidence.items():
        if isinstance(data, dict):
            for key, value in data.items():
                evidence_rows.append([source, key, str(value)])
        elif isinstance(data, list):
            for i, item in enumerate(data):
                for key, value in item.items():
                    evidence_rows.append([source, f"{key}_{i}", str(value)])
    
    evidence_df = pd.DataFrame(evidence_rows, columns=['evidence_source', 'parameter', 'value'])
    evidence_df.to_csv(args.output_evidence, sep='\\t', index=False)
    
    logger.info("Integrated classification completed successfully")

if __name__ == "__main__":
    main()