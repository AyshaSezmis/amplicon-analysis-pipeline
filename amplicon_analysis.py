#!/usr/bin/env python3
"""
Amplicon Sequencing Analysis Pipeline

This script processes Nanopore amplicon sequencing data to calculate species proportions
using Minimap2 for alignment and read classification.

Requirements:
- Python 3.8+
- Biopython
- Minimap2
"""

import os
import gzip
import re
import argparse
import csv
import logging
from typing import List, Dict, Tuple
from multiprocessing import Pool
from Bio import SeqIO
import subprocess

# Configure logging
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

class AmpliconAnalysisPipeline:
    """
    A comprehensive pipeline for processing amplicon sequencing data
    and calculating species proportions.
    """
    def __init__(
        self, 
        input_fastq_dir: str, 
        reference_fasta: str, 
        threads: int = 8,
        min_alignment_score: int = 50
    ):
        """
        Initialize the amplicon analysis pipeline.

        Args:
            input_fastq_dir (str): Directory containing input FASTQ.gz files
            reference_fasta (str): Path to reference sequences FASTA file
            threads (int, optional): Number of processing threads. Defaults to 8.
            min_alignment_score (int, optional): Minimum alignment score. Defaults to 50.
        """
        self.input_fastq_dir = input_fastq_dir
        self.reference_fasta = reference_fasta
        self.threads = threads
        self.min_alignment_score = min_alignment_score
        
        # Validate input paths
        self._validate_input_paths()

        # Load reference sequences
        self.ref_sequences = list(SeqIO.parse(reference_fasta, 'fasta'))
        logger.info(f"Loaded {len(self.ref_sequences)} reference sequences")

    def _validate_input_paths(self):
        """
        Validate existence of input directories and files.
        
        Raises:
            FileNotFoundError: If input paths do not exist
        """
        if not os.path.exists(self.input_fastq_dir):
            raise FileNotFoundError(f"Input directory not found: {self.input_fastq_dir}")
        
        if not os.path.exists(self.reference_fasta):
            raise FileNotFoundError(f"Reference FASTA file not found: {self.reference_fasta}")

    def index_reference_genome(self) -> None:
        """
        Index reference genome using Minimap2.
        
        Raises:
            subprocess.CalledProcessError: If indexing fails
        """
        try:
            subprocess.run(
                ['minimap2', '-d', f'{self.reference_fasta}.mmi', self.reference_fasta], 
                check=True
            )
            logger.info(f"Successfully indexed reference genome: {self.reference_fasta}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to index reference genome: {e}")
            raise

    def run_minimap2(self, reads_file: str, output_sam: str) -> None:
        """
        Run Minimap2 alignment on input reads.

        Args:
            reads_file (str): Path to input reads file
            output_sam (str): Path to output SAM file
        """
        try:
            with open(output_sam, 'w') as out_sam:
                if reads_file.endswith('.gz'):
                    with gzip.open(reads_file, 'rt') as reads:
                        subprocess.run(
                            ['minimap2', '-ax', 'map-ont', '--secondary=no', 
                             f'{self.reference_fasta}.mmi', '-'], 
                            stdin=reads, 
                            stdout=out_sam, 
                            check=True
                        )
                else:
                    subprocess.run(
                        ['minimap2', '-ax', 'map-ont', '--secondary=no', 
                         f'{self.reference_fasta}.mmi', reads_file], 
                        stdout=out_sam, 
                        check=True
                    )
            logger.info(f"Completed Minimap2 alignment for {reads_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Minimap2 alignment failed for {reads_file}: {e}")
            raise

    def parse_sam_file(self, sam_file: str) -> Dict[str, Tuple[str, int]]:
        """
        Parse SAM file to extract best matches for each read.

        Args:
            sam_file (str): Path to input SAM file

        Returns:
            Dict of read IDs to (reference_sequence, alignment_score)
        """
        best_matches = {}
        read_count, aligned_count, tie_count = 0, 0, 0

        with open(sam_file, 'r') as file:
            for line in file:
                if line.startswith('@'):
                    continue  # Skip header lines

                read_count += 1
                fields = line.strip().split('\t')
                read_id = fields[0]
                ref_sequence = fields[2]
                
                # Extract alignment score
                alignment_score = next(
                    (int(value.split(':')[2]) for value in fields if value.startswith('AS:i:')), 
                    0
                )

                if alignment_score >= self.min_alignment_score:
                    aligned_count += 1
                    if read_id in best_matches:
                        # Check for tie
                        if best_matches[read_id][1] == alignment_score:
                            best_matches[read_id] = (None, alignment_score)  # Indicate a tie
                            tie_count += 1
                        elif best_matches[read_id][1] < alignment_score:
                            best_matches[read_id] = (ref_sequence, alignment_score)
                    else:
                        best_matches[read_id] = (ref_sequence, alignment_score)

        # Log summary statistics
        logger.info(f"Total reads processed: {read_count}")
        logger.info(f"Reads with alignment score >= {self.min_alignment_score}: {aligned_count}")
        logger.info(f"Reads with alignment ties: {tie_count}")

        # Filter out reads with ties
        return {k: v for k, v in best_matches.items() if v[0] is not None}

    def process_fastq(self, input_fastq: str) -> Tuple[str, str, List[Tuple[str, int, float]]]:
        """
        Process a single FASTQ file through alignment and species proportion calculation.

        Args:
            input_fastq (str): Path to input FASTQ file

        Returns:
            Tuple containing filename, results string, and species data
        """
        logger.info(f"Processing {input_fastq}")
        output_sam = input_fastq.replace('.fastq.gz', '.sam')
        
        # Run alignment
        self.run_minimap2(input_fastq, output_sam)
        
        # Parse alignment results
        best_matches = self.parse_sam_file(output_sam)

        # Initialize species counts
        species_counts = {ref.id: 0 for ref in self.ref_sequences}

        # Count matches per species
        for read_id, (ref, score) in best_matches.items():
            species_counts[ref] += 1

        # Calculate proportions
        total_counts = sum(species_counts.values())
        if total_counts == 0:
            logger.warning(f"No reads in {input_fastq} could be aligned to any reference sequence.")
            return input_fastq, "No reads could be aligned to any reference sequence.", []

        # Compute species proportions
        species_proportions = {species: count / total_counts for species, count in species_counts.items()}

        # Prepare results
        results_str = "Species counts and proportions:\n"
        species_data = []
        for species, count in species_counts.items():
            proportion = species_proportions[species]
            results_str += f"{species}: {count} reads ({proportion:.2%})\n"
            species_data.append((species, count, proportion))

        return input_fastq, results_str, species_data

    def run_analysis(self) -> None:
        """
        Run the complete amplicon analysis pipeline.
        Process all input FASTQ files and generate results.
        """
        # Index reference genome
        self.index_reference_genome()

        # Find input FASTQ files
        input_fastq_files = [
            os.path.join(self.input_fastq_dir, file) 
            for file in sorted(os.listdir(self.input_fastq_dir)) 
            if file.endswith('.fastq.gz')
        ]

        # Process files in parallel
        results = {}
        with Pool(processes=self.threads) as pool:
            for input_fastq, results_str, species_data in pool.map(self.process_fastq, input_fastq_files):
                results[input_fastq] = (results_str, species_data)

        # Write results to CSV
        output_csv = 'AmpliconMinimap2_results.csv'
        with open(output_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['File', 'Species', 'Count', 'Proportion'])

            # Sort and write results
            for input_fastq in sorted(results.keys(), key=self._sort_filenames):
                logger.info(f"Results for {input_fastq}")
                for sequence, count, proportion in results[input_fastq][1]:
                    logger.info(f"{sequence}: {count} reads ({proportion:.2%})")
                    writer.writerow([input_fastq, sequence, count, f"{proportion:.2%}"])

        logger.info(f"Analysis complete. Results saved to {output_csv}")

    @staticmethod
    def _sort_filenames(filename: str) -> Tuple[str, int]:
        """
        Sort filenames based on their alphanumeric components.

        Args:
            filename (str): Filename to be sorted

        Returns:
            Tuple for sorting (alpha part, numeric part)
        """
        alpha_part = re.sub(r'\d+', '', filename)
        numeric_part = int(re.search(r'\d+', filename).group()) if re.search(r'\d+', filename) else 0
        return (alpha_part, numeric_part)

def main():
    """
    Main function to parse arguments and run the pipeline.
    """
    parser = argparse.ArgumentParser(
        description='Calculate species proportions from amplicon sequencing data.'
    )
    parser.add_argument(
        '--input-dir', 
        required=True, 
        help='Directory containing input FASTQ.gz files'
    )
    parser.add_argument(
        '--reference', 
        required=True, 
        help='Path to reference sequences FASTA file'
    )
    parser.add_argument(
        '--threads', 
        type=int, 
        default=8, 
        help='Number of threads to use for computation (default: 8)'
    )
    parser.add_argument(
        '--min-score', 
        type=int, 
        default=50, 
        help='Minimum alignment score (default: 50)'
    )
    
    args = parser.parse_args()

    try:
        pipeline = AmpliconAnalysisPipeline(
            input_fastq_dir=args.input_dir,
            reference_fasta=args.reference,
            threads=args.threads,
            min_alignment_score=args.min_score
        )
        pipeline.run_analysis()
    except Exception as e:
        logger.error(f"Pipeline execution failed: {e}")
        raise

if __name__ == '__main__':
    main()