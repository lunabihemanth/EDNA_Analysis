
import os
from pathlib import Path

# Application settings
APP_NAME = "DNA Taxonomic Classifier"
VERSION = "1.0.0"

# NCBI settings
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DEFAULT_DATABASE = "nt"  # nucleotide database
DEFAULT_BLAST_PROGRAM = "blastn"

# Search parameters
DEFAULT_MAX_RESULTS = 5
DEFAULT_E_VALUE_THRESHOLD = 0.01
DEFAULT_IDENTITY_THRESHOLD = 70.0
DEFAULT_COVERAGE_THRESHOLD = 50.0

# Cache settings
CACHE_DIR = Path.home() / ".dna_classifier_cache"
CACHE_EXPIRY_DAYS = 7

# Confidence score weights
CONFIDENCE_WEIGHTS = {
    'identity': 0.6,
    'coverage': 0.3,
    'e_value': 0.1
}

# Taxonomic ranks in order of hierarchy
TAXONOMIC_RANKS = [
    'superkingdom',
    'kingdom', 
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species'
]

# Output formats
SUPPORTED_OUTPUT_FORMATS = ['json', 'csv', 'txt', 'xml']

# Rate limiting (requests per second)
NCBI_RATE_LIMIT = 3  # Conservative rate limit for NCBI API
'''

with open('config.py', 'w') as f:
    f.write(config_code)

# 3. Enhanced classifier with machine learning capabilities
ml_classifier_code = '''
"""
Enhanced DNA Classifier with Machine Learning capabilities
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, accuracy_score
from sklearn.preprocessing import LabelEncoder
import pickle
import warnings
warnings.filterwarnings('ignore')

class MLDNAClassifier:
    """
    Machine Learning-based DNA classifier for faster local classification
    This can be trained on known DNA sequences and their taxonomic classifications
    """
    
    def __init__(self, kmer_size=6):
        """
        Initialize ML classifier
        
        Args:
            kmer_size: Size of k-mers for feature extraction
        """
        self.kmer_size = kmer_size
        self.model = RandomForestClassifier(n_estimators=100, random_state=42)
        self.label_encoder = LabelEncoder()
        self.feature_names = []
        self.is_trained = False
    
    def extract_kmer_features(self, sequence):
        """
        Extract k-mer frequency features from DNA sequence
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            dict: k-mer frequency features
        """
        sequence = sequence.upper().replace('N', 'A')  # Replace ambiguous bases
        kmers = {}
        
        # Generate all possible k-mers
        bases = ['A', 'T', 'C', 'G']
        all_kmers = []
        
        def generate_kmers(length, current=""):
            if length == 0:
                all_kmers.append(current)
                return
            for base in bases:
                generate_kmers(length - 1, current + base)
        
        generate_kmers(self.kmer_size)
        
        # Initialize kmer counts
        for kmer in all_kmers:
            kmers[kmer] = 0
        
        # Count k-mers in sequence
        total_kmers = 0
        for i in range(len(sequence) - self.kmer_size + 1):
            kmer = sequence[i:i + self.kmer_size]
            if all(base in bases for base in kmer):
                kmers[kmer] += 1
                total_kmers += 1
        
        # Convert to frequencies
        if total_kmers > 0:
            for kmer in kmers:
                kmers[kmer] = kmers[kmer] / total_kmers
        
        return kmers
    
    def prepare_training_data(self, sequences, labels):
        """
        Prepare training data from sequences and labels
        
        Args:
            sequences: List of DNA sequences
            labels: List of taxonomic labels
            
        Returns:
            tuple: (features_matrix, encoded_labels)
        """
        print(f"Extracting k-mer features from {len(sequences)} sequences...")
        
        features_list = []
        for i, seq in enumerate(sequences):
            if i % 100 == 0:
                print(f"Processed {i}/{len(sequences)} sequences")
            
            kmer_features = self.extract_kmer_features(seq)
            features_list.append(kmer_features)
        
        # Convert to DataFrame for consistent feature ordering
        features_df = pd.DataFrame(features_list)
        features_df = features_df.fillna(0)  # Fill missing k-mers with 0
        
        self.feature_names = list(features_df.columns)
        features_matrix = features_df.values
        
        # Encode labels
        encoded_labels = self.label_encoder.fit_transform(labels)
        
        return features_matrix, encoded_labels
    
    def train(self, sequences, labels, test_size=0.2):
        """
        Train the ML classifier
        
        Args:
            sequences: List of DNA sequences
            labels: List of taxonomic labels
            test_size: Fraction of data to use for testing
        """
        print("Training ML DNA classifier...")
        
        # Prepare data
        X, y = self.prepare_training_data(sequences, labels)
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=42, stratify=y
        )
        
        # Train model
        print("Training Random Forest model...")
        self.model.fit(X_train, y_train)
        
        # Evaluate model
        y_pred = self.model.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)
        
        print(f"Model training completed!")
        print(f"Accuracy on test set: {accuracy:.3f}")
        
        # Print classification report
        target_names = self.label_encoder.classes_
        print("\\nClassification Report:")
        print(classification_report(y_test, y_pred, target_names=target_names))
        
        self.is_trained = True
    
    def predict(self, sequence):
        """
        Predict taxonomic classification for a sequence
        
        Args:
            sequence: DNA sequence to classify
            
        Returns:
            dict: Prediction results with probabilities
        """
        if not self.is_trained:
            raise ValueError("Model must be trained before making predictions")
        
        # Extract features
        kmer_features = self.extract_kmer_features(sequence)
        
        # Convert to feature vector
        feature_vector = np.zeros(len(self.feature_names))
        for i, feature_name in enumerate(self.feature_names):
            if feature_name in kmer_features:
                feature_vector[i] = kmer_features[feature_name]
        
        feature_vector = feature_vector.reshape(1, -1)
        
        # Make prediction
        prediction = self.model.predict(feature_vector)[0]
        prediction_proba = self.model.predict_proba(feature_vector)[0]
        
        # Get predicted class name
        predicted_class = self.label_encoder.inverse_transform([prediction])[0]
        
        # Get top predictions with probabilities
        class_names = self.label_encoder.classes_
        top_indices = np.argsort(prediction_proba)[::-1][:5]  # Top 5 predictions
        
        top_predictions = []
        for idx in top_indices:
            top_predictions.append({
                'class': class_names[idx],
                'probability': prediction_proba[idx],
                'confidence': prediction_proba[idx] * 100
            })
        
        return {
            'predicted_class': predicted_class,
            'confidence': prediction_proba[prediction] * 100,
            'top_predictions': top_predictions
        }
    
    def save_model(self, filepath):
        """Save trained model to file"""
        if not self.is_trained:
            raise ValueError("Model must be trained before saving")
        
        model_data = {
            'model': self.model,
            'label_encoder': self.label_encoder,
            'feature_names': self.feature_names,
            'kmer_size': self.kmer_size,
            'is_trained': self.is_trained
        }
        
        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)
        
        print(f"Model saved to {filepath}")
    
    def load_model(self, filepath):
        """Load trained model from file"""
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)
        
        self.model = model_data['model']
        self.label_encoder = model_data['label_encoder']
        self.feature_names = model_data['feature_names']
        self.kmer_size = model_data['kmer_size']
        self.is_trained = model_data['is_trained']
        
        print(f"Model loaded from {filepath}")

# Example usage and training script
def create_sample_training_data():
    """
    Create sample training data for demonstration
    In practice, you would load real DNA sequences with known classifications
    """
    # This is just example data - replace with real sequences
    sample_data = [
        ("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", "Bacteria"),
        ("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT", "Bacteria"),
        ("TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA", "Archaea"),
        ("CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG", "Archaea"),
        ("AGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT", "Eukaryota"),
        ("TGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA", "Eukaryota"),
    ]
    
    sequences = [item[0] for item in sample_data]
    labels = [item[1] for item in sample_data]
    
    return sequences, labels

if __name__ == "__main__":
    # Example of training and using the ML classifier
    print("DNA ML Classifier Example")
    print("=" * 30)
    
    # Create sample data
    sequences, labels = create_sample_training_data()
    
    # Initialize and train classifier
    ml_classifier = MLDNAClassifier(kmer_size=4)
    
    # Note: This is just example code with minimal data
    # In practice, you need thousands of sequences for effective training
    print("Note: This example uses minimal sample data.")
    print("For real applications, use large datasets of known sequences.")
    
    try:
        ml_classifier.train(sequences, labels)
        
        # Make a prediction
        test_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        result = ml_classifier.predict(test_sequence)
        
        print(f"\\nPrediction for test sequence:")
        print(f"Predicted class: {result['predicted_class']}")
        print(f"Confidence: {result['confidence']:.1f}%")
        
        # Save model
        ml_classifier.save_model("dna_ml_model.pkl")
        
    except Exception as e:
        print(f"Training requires more diverse data: {e}")


with open('ml_classifier.py', 'w') as f:
    f.write(ml_classifier_code)

print("Created additional files:")
print("- requirements.txt")
print("- config.py") 
print("- ml_classifier.py")