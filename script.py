# Create a comprehensive DNA taxonomic classification application
import os

# First, let's create the main application structure
app_code = '''
#!/usr/bin/env python3
"""
DNA Taxonomic Classification Application

This application takes DNA sequences and classifies them to determine:
- Taxa (Kingdom, Phylum, Class, Order, Family, Genus, Species)
- Organism information
- Confidence scores
- Additional biological information

Author: Your Name
Date: September 2024
"""

import sys
import os
import json
import time
import argparse
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from pathlib import Path
import re

# Required imports (install via: pip install biopython requests)
try:
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import requests
    from urllib.parse import quote
    import xml.etree.ElementTree as ET
except ImportError as e:
    print(f"Error: Missing required dependencies. Please install with:")
    print("pip install biopython requests")
    sys.exit(1)

@dataclass
class TaxonomicInfo:
    """Data class to store taxonomic information"""
    taxid: str
    scientific_name: str
    kingdom: str = ""
    phylum: str = ""
    class_name: str = ""
    order: str = ""
    family: str = ""
    genus: str = ""
    species: str = ""
    lineage: List[str] = None
    confidence: float = 0.0
    e_value: float = float('inf')
    identity_percent: float = 0.0
    query_coverage: float = 0.0
    
    def __post_init__(self):
        if self.lineage is None:
            self.lineage = []

class DNATaxonomyClassifier:
    """Main class for DNA taxonomic classification"""
    
    def __init__(self, email: str, api_key: Optional[str] = None):
        """
        Initialize the classifier
        
        Args:
            email: Your email for NCBI Entrez (required)
            api_key: NCBI API key (optional but recommended)
        """
        self.email = email
        self.api_key = api_key
        
        # Configure Entrez
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        
        # Cache for taxonomic information
        self.taxonomy_cache = {}
        
    def validate_dna_sequence(self, sequence: str) -> bool:
        """
        Validate if the input is a valid DNA sequence
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            bool: True if valid DNA sequence
        """
        # Remove whitespace and convert to uppercase
        clean_seq = re.sub(r'\\s+', '', sequence.upper())
        
        # Check if sequence contains only valid DNA bases
        valid_bases = set('ATCGN')
        return all(base in valid_bases for base in clean_seq)
    
    def clean_sequence(self, sequence: str) -> str:
        """
        Clean and prepare DNA sequence for analysis
        
        Args:
            sequence: Raw DNA sequence
            
        Returns:
            str: Cleaned DNA sequence
        """
        # Remove FASTA header if present
        if sequence.startswith('>'):
            lines = sequence.split('\\n')
            sequence = '\\n'.join(lines[1:])
        
        # Remove whitespace and newlines
        clean_seq = re.sub(r'\\s+', '', sequence.upper())
        
        return clean_seq
    
    def blast_sequence(self, sequence: str, database: str = "nt", 
                      hitlist_size: int = 10) -> List[Dict]:
        """
        Perform BLAST search against NCBI database
        
        Args:
            sequence: DNA sequence to search
            database: NCBI database to search against
            hitlist_size: Maximum number of hits to return
            
        Returns:
            List[Dict]: BLAST results
        """
        print(f"Performing BLAST search against {database} database...")
        
        try:
            # Perform BLAST search
            result_handle = NCBIWWW.qblast(
                "blastn", 
                database, 
                sequence,
                hitlist_size=hitlist_size,
                expect=0.01
            )
            
            # Parse results
            blast_records = NCBIXML.parse(result_handle)
            blast_results = []
            
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        # Extract accession number
                        accession = alignment.accession
                        
                        # Calculate metrics
                        identity_percent = (hsp.identities / hsp.align_length) * 100
                        query_coverage = (hsp.align_length / blast_record.query_length) * 100
                        
                        blast_results.append({
                            'accession': accession,
                            'title': alignment.title,
                            'e_value': hsp.expect,
                            'identity_percent': identity_percent,
                            'query_coverage': query_coverage,
                            'align_length': hsp.align_length,
                            'score': hsp.score
                        })
            
            result_handle.close()
            
            # Sort by e-value (best matches first)
            blast_results.sort(key=lambda x: x['e_value'])
            
            return blast_results
            
        except Exception as e:
            print(f"Error during BLAST search: {e}")
            return []
    
    def get_taxonomy_info(self, accession: str) -> Optional[TaxonomicInfo]:
        """
        Get taxonomic information for a given accession number
        
        Args:
            accession: NCBI accession number
            
        Returns:
            TaxonomicInfo: Taxonomic information or None
        """
        if accession in self.taxonomy_cache:
            return self.taxonomy_cache[accession]
        
        try:
            # Get nucleotide record
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            # Extract taxonomy information from features
            organism = None
            taxonomy = []
            
            for feature in record.features:
                if feature.type == "source":
                    organism = feature.qualifiers.get("organism", [None])[0]
                    taxonomy = feature.qualifiers.get("db_xref", [])
                    break
            
            if not organism:
                return None
            
            # Get taxon ID
            taxid = None
            for xref in taxonomy:
                if xref.startswith("taxon:"):
                    taxid = xref.split(":")[1]
                    break
            
            if not taxid:
                return None
            
            # Get detailed taxonomic lineage
            lineage_info = self.get_taxonomic_lineage(taxid)
            
            tax_info = TaxonomicInfo(
                taxid=taxid,
                scientific_name=organism,
                lineage=lineage_info.get('lineage', []),
                kingdom=lineage_info.get('kingdom', ''),
                phylum=lineage_info.get('phylum', ''),
                class_name=lineage_info.get('class', ''),
                order=lineage_info.get('order', ''),
                family=lineage_info.get('family', ''),
                genus=lineage_info.get('genus', ''),
                species=lineage_info.get('species', '')
            )
            
            # Cache the result
            self.taxonomy_cache[accession] = tax_info
            
            return tax_info
            
        except Exception as e:
            print(f"Error getting taxonomy for {accession}: {e}")
            return None
    
    def get_taxonomic_lineage(self, taxid: str) -> Dict[str, str]:
        """
        Get complete taxonomic lineage for a given taxon ID
        
        Args:
            taxid: NCBI taxonomy ID
            
        Returns:
            Dict: Taxonomic lineage information
        """
        try:
            # Fetch taxonomy record
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            lineage_info = {
                'lineage': [],
                'kingdom': '',
                'phylum': '',
                'class': '',
                'order': '',
                'family': '',
                'genus': '',
                'species': ''
            }
            
            if records:
                record = records[0]
                
                # Get lineage
                if 'LineageEx' in record:
                    for item in record['LineageEx']:
                        rank = item.get('Rank', '')
                        name = item.get('ScientificName', '')
                        
                        lineage_info['lineage'].append(f"{name} ({rank})")
                        
                        # Map to standard taxonomic levels
                        if rank == 'kingdom':
                            lineage_info['kingdom'] = name
                        elif rank == 'phylum':
                            lineage_info['phylum'] = name
                        elif rank == 'class':
                            lineage_info['class'] = name
                        elif rank == 'order':
                            lineage_info['order'] = name
                        elif rank == 'family':
                            lineage_info['family'] = name
                        elif rank == 'genus':
                            lineage_info['genus'] = name
                        elif rank == 'species':
                            lineage_info['species'] = name
            
            return lineage_info
            
        except Exception as e:
            print(f"Error getting lineage for taxid {taxid}: {e}")
            return {'lineage': []}
    
    def classify_sequence(self, sequence: str, max_results: int = 5) -> List[TaxonomicInfo]:
        """
        Main method to classify a DNA sequence
        
        Args:
            sequence: DNA sequence to classify
            max_results: Maximum number of classification results to return
            
        Returns:
            List[TaxonomicInfo]: List of taxonomic classifications
        """
        print("Starting DNA sequence classification...")
        
        # Validate and clean sequence
        if not self.validate_dna_sequence(sequence):
            raise ValueError("Invalid DNA sequence provided")
        
        clean_seq = self.clean_sequence(sequence)
        
        if len(clean_seq) < 50:
            print("Warning: Very short sequence may give unreliable results")
        
        # Perform BLAST search
        blast_results = self.blast_sequence(clean_seq, hitlist_size=max_results * 3)
        
        if not blast_results:
            print("No BLAST results found")
            return []
        
        print(f"Found {len(blast_results)} BLAST hits, retrieving taxonomic information...")
        
        # Get taxonomic information for top hits
        classifications = []
        processed = 0
        
        for result in blast_results:
            if processed >= max_results:
                break
            
            accession = result['accession']
            tax_info = self.get_taxonomy_info(accession)
            
            if tax_info:
                # Add BLAST metrics to taxonomic info
                tax_info.e_value = result['e_value']
                tax_info.identity_percent = result['identity_percent']
                tax_info.query_coverage = result['query_coverage']
                
                # Calculate confidence score (simple heuristic)
                tax_info.confidence = min(
                    100.0,
                    (tax_info.identity_percent * 0.6) + 
                    (tax_info.query_coverage * 0.3) + 
                    (max(0, 20 - (-1 * log10(max(tax_info.e_value, 1e-200)))) * 0.1)
                )
                
                classifications.append(tax_info)
                processed += 1
            
            # Add small delay to be respectful to NCBI servers
            time.sleep(0.1)
        
        return classifications
    
    def format_results(self, classifications: List[TaxonomicInfo]) -> str:
        """
        Format classification results for display
        
        Args:
            classifications: List of taxonomic classifications
            
        Returns:
            str: Formatted results
        """
        if not classifications:
            return "No classifications found."
        
        output = []
        output.append("DNA SEQUENCE CLASSIFICATION RESULTS")
        output.append("=" * 50)
        
        for i, tax_info in enumerate(classifications, 1):
            output.append(f"\\nResult #{i}:")
            output.append(f"Scientific Name: {tax_info.scientific_name}")
            output.append(f"Taxonomy ID: {tax_info.taxid}")
            output.append(f"Confidence Score: {tax_info.confidence:.1f}%")
            output.append(f"Identity: {tax_info.identity_percent:.1f}%")
            output.append(f"Query Coverage: {tax_info.query_coverage:.1f}%")
            output.append(f"E-value: {tax_info.e_value:.2e}")
            
            output.append("\\nTaxonomic Classification:")
            if tax_info.kingdom:
                output.append(f"  Kingdom: {tax_info.kingdom}")
            if tax_info.phylum:
                output.append(f"  Phylum: {tax_info.phylum}")
            if tax_info.class_name:
                output.append(f"  Class: {tax_info.class_name}")
            if tax_info.order:
                output.append(f"  Order: {tax_info.order}")
            if tax_info.family:
                output.append(f"  Family: {tax_info.family}")
            if tax_info.genus:
                output.append(f"  Genus: {tax_info.genus}")
            if tax_info.species:
                output.append(f"  Species: {tax_info.species}")
            
            if tax_info.lineage:
                output.append("\\nComplete Lineage:")
                for lineage_item in tax_info.lineage:
                    output.append(f"  {lineage_item}")
            
            output.append("-" * 30)
        
        return "\\n".join(output)
    
    def save_results(self, classifications: List[TaxonomicInfo], filename: str):
        """
        Save classification results to JSON file
        
        Args:
            classifications: List of taxonomic classifications
            filename: Output filename
        """
        results_data = []
        
        for tax_info in classifications:
            result_dict = {
                'taxid': tax_info.taxid,
                'scientific_name': tax_info.scientific_name,
                'kingdom': tax_info.kingdom,
                'phylum': tax_info.phylum,
                'class': tax_info.class_name,
                'order': tax_info.order,
                'family': tax_info.family,
                'genus': tax_info.genus,
                'species': tax_info.species,
                'lineage': tax_info.lineage,
                'confidence': tax_info.confidence,
                'e_value': tax_info.e_value,
                'identity_percent': tax_info.identity_percent,
                'query_coverage': tax_info.query_coverage
            }
            results_data.append(result_dict)
        
        with open(filename, 'w') as f:
            json.dump(results_data, f, indent=2)
        
        print(f"Results saved to {filename}")

def main():
    """Main function for command-line interface"""
    parser = argparse.ArgumentParser(
        description="DNA Taxonomic Classification Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python dna_classifier.py --email your@email.com --sequence "ATCGATCGATCG..."
  python dna_classifier.py --email your@email.com --file sequence.fasta
  python dna_classifier.py --email your@email.com --sequence "ATCG..." --output results.json
        """
    )
    
    parser.add_argument('--email', required=True, 
                       help='Your email address (required by NCBI)')
    parser.add_argument('--api-key', 
                       help='NCBI API key (optional but recommended)')
    parser.add_argument('--sequence', 
                       help='DNA sequence to classify')
    parser.add_argument('--file', 
                       help='FASTA file containing DNA sequence')
    parser.add_argument('--output', 
                       help='Output JSON file for results')
    parser.add_argument('--max-results', type=int, default=5,
                       help='Maximum number of classification results (default: 5)')
    
    args = parser.parse_args()
    
    # Validate input
    if not args.sequence and not args.file:
        print("Error: Either --sequence or --file must be provided")
        sys.exit(1)
    
    # Get DNA sequence
    if args.file:
        try:
            if args.file.endswith(('.fasta', '.fa', '.fas')):
                with open(args.file, 'r') as f:
                    record = SeqIO.read(f, 'fasta')
                    sequence = str(record.seq)
            else:
                with open(args.file, 'r') as f:
                    sequence = f.read().strip()
        except Exception as e:
            print(f"Error reading file {args.file}: {e}")
            sys.exit(1)
    else:
        sequence = args.sequence
    
    try:
        # Initialize classifier
        classifier = DNATaxonomyClassifier(args.email, args.api_key)
        
        # Classify sequence
        results = classifier.classify_sequence(sequence, args.max_results)
        
        # Display results
        print(classifier.format_results(results))
        
        # Save results if requested
        if args.output:
            classifier.save_results(results, args.output)
    
    except Exception as e:
        print(f"Error during classification: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
'''

# Write the main application file
with open('dna_classifier.py', 'w') as f:
    f.write(app_code)

print("Created dna_classifier.py - Main application file")