# Amplicon Sequencing Analysis Pipeline

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [Arguments](#arguments)
  - [Example](#example)
- [Output](#output)
- [Contributing](#contributing)
- [License](#license)
- [Citation](#citation)
- [Authors](#authors)

## Introduction
The **Amplicon Sequencing Analysis Pipeline** is a comprehensive tool designed to process Nanopore amplicon sequencing data and calculate species proportions using Minimap2 for alignment and read classification. This pipeline is optimized for efficiency and accuracy, supporting parallel processing to handle large datasets seamlessly.

## Features
- **Input Validation:** Ensures valid paths for input files and reference sequences.
- **Minimap2 Integration:** Performs robust alignments for amplicon reads.
- **Parallel Processing:** Leverages multicore systems for faster data processing.
- **Detailed Outputs:** Generates CSV files with species proportions and aligned read counts.

## Prerequisites
Before installing the pipeline, ensure that the following dependencies are met:
- Operating System: Linux or macOS recommended. Windows may require additional configuration.
- Python: Version 3.8 or higher.
- Minimap2: Ensure Minimap2 is installed and accessible in your system's `PATH`. Installation instructions are provided below.

## Installation
### 1. Clone the Repository
Clone the repository to your local machine:
```bash
git clone https://github.com/AyshaSezmis/amplicon-analysis-pipeline.git
cd amplicon-analysis-pipeline
```
### 2. Install Minimap2
Minimap2 is a versatile sequence alignment program used in this pipeline. Install it using one of the following methods:
Using Conda (Recommended)
If you have Conda installed, you can install Minimap2 via the bioconda channel:
```bash
conda install -c bioconda minimap2
```
Using Homebrew (macOS)
For macOS users with Homebrew installed:
```bash
brew install minimap2
```
Verify Installation:
Ensure Minimap2 is installed correctly by checking its version:
```bash
minimap2 --version
```
### 3. Set Up Python Environment
It's recommended to use a virtual environment to manage Python dependencies:
```bash
# Create a virtual environment named 'venv'
python3 -m venv venv
# Activate the virtual environment
source venv/bin/activate  # On Windows: venv\Scripts\activate
# Upgrade pip
pip install --upgrade pip
# Install Python dependencies
pip install -r requirements.txt 
```

## Usage
Run the pipeline using the command-line interface. Below is the basic syntax and detailed explanation of the arguments.
Arguments:
- `--input-dir (required)`: Directory containing input .fastq.gz files.
- `--reference (required)`: Path to the reference sequences FASTA file.
- `--threads (optional)`: Number of threads to use for computation (default: 8).
- `--min-score (optional)`: Minimum alignment score to consider (default: 50).
Example:
Here's an example of how to run the pipeline:
```bash
python3 amplicon_analysis.py \
  --input-dir /home/user/data/fastq_files \
  --reference /home/user/data/reference.fasta \
  --threads 16 \
  --min-score 60
```
Explanation:
- `--input-dir`: Specifies the directory containing your .fastq.gz files.
- `--reference`: Path to your reference FASTA file containing 16S sequences.
- `--threads`: Utilizes 16 CPU threads for parallel processing.
- `--min-score`: Sets the minimum alignment score threshold to 60.

## Output
- CSV File: AmpliconMinimap2_results.csv containing the following columns:
    - File: Input FASTQ file name.
    - Species: Reference species name.
    - Count: Number of reads aligned to the species.
    - Proportion: Proportion of reads aligned to the species.

- SAM Files: Generated for each input FASTQ file, stored in the same directory as the input files.

## Contributing
Contributions are welcome! If you'd like to enhance the pipeline or fix any issues, please follow these steps:
1) Fork the Repository
    Click the "Fork" button at the top-right corner of this repository's page to create a personal copy.
2) Create a Feature Branch
```bash
git checkout -b feature/YourFeature
```
3) Commit Your Changes
```bash
git commit -m "Add your feature"
```
4) Push to the Branch
```bash
git push origin feature/YourFeature
```
5) Open a Pull Request
    Navigate to your forked repository on GitHub and click the "Compare & pull request" button.
Please ensure that your code follows the existing style and includes appropriate tests.

## License
This project is licensed under the terms of the MIT license. For more information, see the [`LICENSE.md`](LICENSE.md) file.

## Citation
If you use this pipeline in your research, please cite:

Gareth Howells, Aysha L. Sezmis, hristopher Blake, and Michael J. McDonald "Co-existence slows diversification in experimental populations of *E. coli* and *P. fluorescens*." In press

## Authors
- Aysha Laila Sezmis
- Christopher Blake
