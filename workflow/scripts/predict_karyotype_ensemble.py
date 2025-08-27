#!/usr/bin/env python3
"""
Karyotype Prediction using Ensemble Classifier
==============================================
Uses the trained ensemble_classifier_250524.pkl model on RNASeqCNV data.
This is the CORRECT karyotype prediction system.
"""

import csv
import pickle
import sys
import os
from collections import Counter
import numpy as np
import pandas as pd
from pathlib import Path

np.set_printoptions(suppress=True)

def generate_dummy_output(outfile, sample_id):
    """Generate dummy output when insufficient data or model unavailable."""
    with open(outfile, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write headers
        writer.writerow(['Sample', 'Score', 'Prediction', 'Hyperdiploid', 'Low hypodiploid', 'Near haploid', 'iAMP21', 'Other'])
        writer.writerow([sample_id, '0%', 'unclassified', '0%', '0%', '0%', '0%', '0%'])

def generate_prediction_statement(predictions):
    """Generate prediction statement with confidence."""
    # Count frequency of predictions
    prediction_counts = Counter(predictions)
    print("prediction_counts", prediction_counts)
    
    # Find most common predicted subtype
    most_common_subtype = prediction_counts.most_common(1)[0][0]
    print("most_common_subtype", most_common_subtype)
    print(prediction_counts[most_common_subtype])
    
    # Calculate confidence percentage
    confidence_percentage = (prediction_counts[most_common_subtype] / len(predictions)) * 100
    
    # Create statement
    prediction_statement = f"Predicted subtype is {most_common_subtype} with confidence {confidence_percentage:.2f}%"
    return prediction_statement

def load_ensemble_model(model_path):
    """Load the trained ensemble classifier model."""
    try:
        with open(model_path, 'rb') as model_file:
            model_data = pickle.load(model_file)
            rf_classifier = model_data['model']
            scaler = model_data['scaler']
            label_encoder = model_data['label_encoder']
            return rf_classifier, scaler, label_encoder
    except Exception as e:
        print(f"Error loading model from {model_path}: {e}")
        return None, None, None

def predict_karyotype(input_file, model_path, sample_id):
    """Main karyotype prediction function."""
    
    # Load the trained ensemble model
    rf_classifier, scaler, label_encoder = load_ensemble_model(model_path)
    
    if rf_classifier is None:
        print(f"Could not load model. Generating dummy output for {sample_id}")
        return {
            'sample': sample_id,
            'score': '0%',
            'prediction': 'unclassified',
            'hyperdiploid': '0%',
            'low_hypodiploid': '0%',
            'near_haploid': '0%',
            'iamp21': '0%',
            'other': '0%'
        }
    
    try:
        # Read input data
        new = pd.read_csv(input_file)
        
        # Check column count (needs >= 185 features)
        if new.shape[1] < 185:
            print(f"Input file {input_file} has fewer than 185 columns. Generating dummy output.")
            return {
                'sample': sample_id,
                'score': '0%', 
                'prediction': 'unclassified',
                'hyperdiploid': '0%',
                'low_hypodiploid': '0%',
                'near_haploid': '0%',
                'iamp21': '0%',
                'other': '0%'
            }
        
        # Extract features (skip first column, usually sample ID)
        new_data = new.iloc[:, 1:]
        
        # Ensure columns match what the scaler expects
        expected_features = scaler.feature_names_in_
        
        # Add missing columns with 0
        for feature in expected_features:
            if feature not in new_data.columns:
                new_data[feature] = 0
        
        # Remove unknown columns
        new_data = new_data[expected_features]
        
        # Scale data with loaded scaler
        scaled_data = scaler.transform(new_data)
        
        # Make predictions
        pred = rf_classifier.predict(scaled_data)
        pred_pro = rf_classifier.predict_proba(scaled_data)
        
        # Extract predicted label
        predicted_label = label_encoder.inverse_transform([pred[0]])[0]
        index_predicted_label = list(rf_classifier.classes_).index(pred[0])
        probability = pred_pro[0][index_predicted_label]
        
        # Convert probabilities to percentages
        print(f"Probability: {probability}, All probabilities: {pred_pro[0]}")
        pred_pro_percent = [f"{int(prob * 100)}%" for prob in pred_pro[0]]
        
        # Map class probabilities to specific categories
        # This mapping depends on the classes in your trained model
        classes = [label_encoder.inverse_transform([cls])[0] for cls in rf_classifier.classes_]
        
        result = {
            'sample': sample_id,
            'score': f"{int(probability * 100)}%",
            'prediction': predicted_label,
            'hyperdiploid': '0%',
            'low_hypodiploid': '0%', 
            'near_haploid': '0%',
            'iamp21': '0%',
            'other': '0%'
        }
        
        # Map probabilities to categories
        for i, cls in enumerate(classes):
            prob_percent = f"{int(pred_pro[0][i] * 100)}%"
            if cls == 'Hyperdiploid':
                result['hyperdiploid'] = prob_percent
            elif cls == 'Low hypodiploid':
                result['low_hypodiploid'] = prob_percent
            elif cls == 'Near haploid':
                result['near_haploid'] = prob_percent
            elif cls == 'iAMP21':
                result['iamp21'] = prob_percent
            else:
                result['other'] = prob_percent
        
        return result
        
    except Exception as e:
        print(f"Error in prediction: {e}")
        return {
            'sample': sample_id,
            'score': '0%',
            'prediction': 'error',
            'hyperdiploid': '0%',
            'low_hypodiploid': '0%',
            'near_haploid': '0%',
            'iamp21': '0%',
            'other': '0%'
        }

def main():
    """Main function for Snakemake integration."""
    
    # Check if running from command line or snakemake
    import sys
    if len(sys.argv) > 1:
        # Command line mode for testing
        input_file = sys.argv[1]
        model_path = sys.argv[2]
        output_file = sys.argv[3]
        sample_id = sys.argv[4]
    else:
        # Snakemake mode
        input_file = snakemake.input.cnv_features
        model_path = snakemake.input.ensemble_model  
        output_file = snakemake.output.prediction
        sample_id = snakemake.wildcards.sample
    
    print(f"Running karyotype prediction for sample: {sample_id}")
    print(f"Input file: {input_file}")
    print(f"Model path: {model_path}")
    
    try:
        # Run prediction
        result = predict_karyotype(input_file, model_path, sample_id)
        
        # Write output in original format
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            
            # Write headers
            writer.writerow(['Sample', 'Score', 'Prediction', 'Hyperdiploid', 
                           'Low hypodiploid', 'Near haploid', 'iAMP21', 'Other'])
            
            # Write results
            writer.writerow([
                result['sample'],
                result['score'], 
                result['prediction'],
                result['hyperdiploid'],
                result['low_hypodiploid'],
                result['near_haploid'],
                result['iamp21'],
                result['other']
            ])
        
        print(f"Karyotype prediction completed: {result['prediction']} ({result['score']})")
        
    except Exception as e:
        print(f"Error in karyotype prediction: {e}")
        generate_dummy_output(output_file, sample_id)

if __name__ == "__main__":
    main()