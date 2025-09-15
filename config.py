
"""
Configuration settings for DNA Taxonomic Classification Application
"""

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
