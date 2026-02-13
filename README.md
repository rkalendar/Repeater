# REPEATER

**A Genome-Wide Tool for Rapid *De Novo* Identification, Profiling, Masking, and Visualization of Interspersed and Tandem Repetitive Sequences**

[![DOI](https://img.shields.io/badge/DOI-10.1177%2F11779322241306391-blue)](https://doi.org/10.1177/11779322241306391)
[![Java](https://img.shields.io/badge/Java-25+-orange.svg)](https://www.oracle.com/java/technologies/downloads/)
[![Platform](https://img.shields.io/badge/Platform-Independent-green.svg)]()
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)

---

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Online Access](#online-access)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Java Installation](#java-installation)
  - [Installing Java with Conda](#installing-java-with-conda)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Basic Syntax](#basic-syntax)
  - [Command Examples](#command-examples)
  - [Large Genome Analysis](#large-genome-analysis)
- [Command-Line Options](#command-line-options)
- [Input Format](#input-format)
- [Output Format](#output-format)
- [Advanced Usage](#advanced-usage)
- [Citation](#citation)
- [Contact](#contact)
- [Troubleshooting](#troubleshooting)

---

## 🔬 Overview

REPEATER is a powerful genome-wide tool designed for sensitive *de novo* identification and visualization of repetitive DNA sequences. It detects both interspersed repeats and tandem repeats across entire genomes, providing comprehensive profiling, masking, and visualization capabilities.

**Key Capabilities:**
- 🧬 Detection of interspersed repeats
- 🔄 Identification of tandem repeats
- 🔍 SSR/microsatellite detection
- 📊 Telomeric sequence identification
- 🎨 Visual representation of repeat landscapes
- 🎯 Sequence masking for downstream analysis

**Developed by:** Ruslan Kalendar  
**Institution:** University of Helsinki  
**Platform:** Cross-platform Java application

---

## ✨ Features

- **De Novo Detection**: No prior knowledge of repeat sequences required
- **Comprehensive Analysis**: Identifies all types of repetitive elements
- **Genome-Wide Scale**: Handles entire genomes efficiently
- **Flexible Parameters**: Customizable k-mer size and repeat length thresholds
- **Visual Output**: Generates publication-quality repeat landscape images
- **Sequence Extraction**: Retrieve repeat sequences with flanking regions
- **Multiple Formats**: Supports FASTA and plain text input
- **GFF3 Output**: Standard format for downstream bioinformatics pipelines
- **SSR-Focused Mode**: Specialized analysis for simple sequence repeats
- **Memory Scalable**: Configurable for genomes of any size

---

## 🌐 Online Access

Try REPEATER online without installation:  
**[https://primerdigital.com/tools/repeats.html](https://primerdigital.com/tools/repeats.html)**

The online version provides a user-friendly interface for analyzing sequences up to moderate genome sizes.

---

## 📦 Installation

### Prerequisites

- **Operating System**: Platform independent (Windows, macOS, Linux)
- **Java Runtime**: Java 25 or higher
- **Memory**: 4GB minimum, 16-64GB recommended for large genomes
- **Disk Space**: Sufficient space for input genomes and output files

### Java Installation

#### Option 1: Download from Oracle

1. Visit [Java Downloads](https://www.oracle.com/java/technologies/downloads/)
2. Download Java 25 or newer for your operating system
3. Follow the installation wizard
4. Verify installation:
   ```bash
   java -version
   ```

#### Option 2: Installing Java with Conda

Conda provides an easy way to manage Java installations in isolated environments.

**Step 1: Configure Conda-Forge Channel**

Add the conda-forge channel and set strict priority:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

**Step 2: Create Java Environment**

Create a dedicated environment with OpenJDK 25:

```bash
conda create -n java25 openjdk=25
```

**Step 3: Activate Environment**

```bash
conda activate java25
```

**Step 4: Verify Installation**

```bash
java -version
```

Expected output:
```
openjdk version "25" ...
```

### Setting Java PATH

If you need to set or change the Java PATH variable:
- **Guide**: [How to set Java PATH](https://www.java.com/en/download/help/path.html)

---

## 🚀 Quick Start

1. **Download or clone this repository**

2. **Locate the executable**  
   The `Repeater2.jar` file is in the `dist` directory

3. **Run on a test file**
   ```bash
   java -jar dist/Repeater2.jar test/4.txt
   ```

4. **Check the output**  
   Results will be saved in the same directory as your input file

---

## 📖 Usage

### Basic Syntax

```bash
java -jar [JAVA_OPTIONS] Repeater2.jar <INPUT_PATH> [OPTIONS]
```

**Components:**
- `[JAVA_OPTIONS]`: Optional JVM settings (memory allocation, etc.)
- `Repeater2.jar`: Path to the executable JAR file
- `<INPUT_PATH>`: Path to input file or directory (required)
- `[OPTIONS]`: Command-line parameters for analysis configuration

### Command Examples

#### Example 1: Basic Analysis
```bash
java -jar dist/Repeater2.jar test/4.txt
```

#### Example 2: Custom Parameters with Image Output
```bash
java -jar dist/Repeater2.jar test/ kmer=18 sln=90 image=5000x300
```

**Parameters used:**
- `kmer=18`: Use 18-mer for detection
- `sln=90`: Minimum string length for clustering
- `image=5000x300`: Generate 5000×300 pixel image

#### Example 3: SSR Analysis with Sequence Extraction
```bash
java -jar dist/Repeater2.jar test/2.txt -ssronly -seqshow flanks=100
```

**Parameters used:**
- `-ssronly`: Analyze only SSR/telomeric loci
- `-seqshow`: Extract repeat sequences
- `flanks=100`: Include 100bp flanking regions

#### Example 4: Full Genome Analysis
```bash
java -jar dist/Repeater2.jar /path/to/genome/ kmer=20 sln=100 image=10000x300
```

#### Example 5: Batch Processing
```bash
java -jar dist/Repeater2.jar /path/to/sequences/
```

Process all sequence files in a directory

### Large Genome Analysis

For large genomes (e.g., human, wheat, pine), increase Java heap memory allocation:

#### Syntax for Large Genomes
```bash
java -Xms<INITIAL_MEMORY> -Xmx<MAX_MEMORY> -jar Repeater2.jar <INPUT_PATH> [OPTIONS]
```

**Memory Flags:**
- `-Xms`: Initial heap size (e.g., `-Xms16g` = start with 16GB)
- `-Xmx`: Maximum heap size (e.g., `-Xmx64g` = up to 64GB)

#### Large Genome Examples

**Example 1: Plant Genome (16-64GB)**
```bash
java -Xms16g -Xmx64g -jar dist/Repeater2.jar genomes/Hydra_vulgaris/ kmer=20 sln=100
```

**Example 2: Mammalian Genome (16-32GB)**
```bash
java -Xms16g -Xmx32g -jar dist/Repeater2.jar genomes/Hordeum_marinum/
```

**Example 3: Human Genome (GRCh38)**
```bash
java -Xms16g -Xmx32g -jar dist/Repeater2.jar genomes/GRCh38.p14/
```

**Memory Recommendations:**
| Genome Size | Min RAM | Recommended RAM | Xms | Xmx |
|-------------|---------|-----------------|-----|-----|
| < 100 Mb    | 2 GB    | 4 GB            | 1g  | 4g  |
| 100-500 Mb  | 4 GB    | 8 GB            | 2g  | 8g  |
| 0.5-1 Gb    | 8 GB    | 16 GB           | 4g  | 16g |
| 1-3 Gb      | 16 GB   | 32 GB           | 8g  | 32g |
| 3-10 Gb     | 32 GB   | 64 GB           | 16g | 64g |
| > 10 Gb     | 64 GB   | 128 GB          | 32g | 128g|

---

## ⚙️ Command-Line Options

### Core Parameters

| Option        | Type    | Default | Description |
|---------------|---------|---------|-------------|
| `kmer=N`      | Integer | 18      | K-mer length for repeat detection (minimum: 12) |
| `min=N`       | Integer | 30      | Minimum repeat length to report (bp) |
| `sln=N`       | Integer | 90      | Minimum string length for clustering analysis (bp) |
| `flanks=N`    | Integer | 0       | Length of flanking sequences to extract (bp) |
| `image=WxH`   | String  | Auto    | Output image dimensions (width×height in pixels) |

### Analysis Modes

| Option       | Description |
|--------------|-------------|
| `-seqshow`   | Extract and output repeat sequences |
| `-ssronly`   | Analyze only SSR/microsatellite and telomeric loci |
| `-maskonly`  | Generate only masked sequence output (skip clustering) |

### Parameter Guidelines

**K-mer Size Selection:**
- **Small genomes (< 100 Mb)**: `kmer=15-18`
- **Medium genomes (100 Mb - 1 Gb)**: `kmer=18-20`
- **Large genomes (> 1 Gb)**: `kmer=20-25`
- Larger k-mers increase specificity but may miss divergent repeats

**Minimum Repeat Length:**
- **SSR detection**: `min=12-30`
- **General repeats**: `min=30-50`
- **Long repeats only**: `min=100+`

**String Length for Clustering:**
- **Fine resolution**: `sln=50-90`
- **Standard**: `sln=90-150`
- **Coarse families**: `sln=200+`

---

## 📄 Input Format

### Supported Formats

1. **Plain Text (.txt)**: Simple text file with DNA sequence
2. **FASTA (.fasta, .fa, .fna)**: Standard FASTA format
3. **No Extension**: Files without extensions are acceptable

### FASTA Format Specification

A valid FASTA file consists of:

**1. Header Line (required):**
- Starts with `>` character
- Followed by sequence identifier
- Optional description after identifier

**2. Sequence Lines:**
- One or more lines containing the sequence
- DNA bases: A, C, G, T (uppercase or lowercase)
- Ambiguous bases supported: N, R, Y, K, M, S, W, B, D, H, V

**Example:**
```
>Chr1 Chromosome 1, complete sequence
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>Chr2 Chromosome 2, complete sequence
AAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT
```

### Multiple Sequences

- Multiple sequences can be included in a single file
- Each sequence must have its own FASTA header
- Sequences are analyzed independently

### Sequence Length

- **No upper limit** on sequence length
- Can analyze entire chromosomes or scaffolds
- Memory allocation should match genome size

---

## 📤 Output Format

### GFF3 Format

REPEATER outputs results in **General Feature Format Version 3 (GFF3)**, a standardized tab-delimited format for genomic annotations.

**Specification:** [GFF3 Documentation](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

### GFF3 Column Description

| Column | Name     | Description |
|--------|----------|-------------|
| 1      | seqid    | Sequence/chromosome identifier (cluster ID) |
| 2      | source   | Program name ("REPEATER") |
| 3      | type     | Feature type ("repeat") |
| 4      | start    | Start position (1-based, inclusive) |
| 5      | end      | End position (1-based, inclusive) |
| 6      | score    | Repeat length (bp) |
| 7      | strand   | Strand orientation ('+', '-', or '.') |
| 8      | phase    | Reading frame (always '.' for repeats) |
| 9      | attributes | Additional information (sequence if `-seqshow` used) |

### Example GFF3 Output

```gff3
##gff-version 3
Chr1	REPEATER	repeat	1000	1500	500	+	.	ID=repeat1;Type=tandem
Chr1	REPEATER	repeat	3500	4200	700	-	.	ID=repeat2;Type=interspersed
Chr2	REPEATER	repeat	800	1100	300	+	.	ID=repeat3;Type=SSR;Motif=AT
```

### Additional Output Files

Depending on options used, REPEATER may generate:

1. **Masked Sequences** (`.masked` extension)
   - Input sequence with repeats masked (lowercase or 'N')

2. **Repeat Sequences** (`.repeats` extension)
   - Extracted repeat sequences with flanking regions

3. **Visualization Images** (`.png` or `.jpg`)
   - Graphical representation of repeat distribution
   - Repeat density plots
   - Chromosome/scaffold overview

4. **Statistics Report** (`.stats` or `.txt`)
   - Summary of repeat content
   - Repeat family statistics
   - Coverage information

---

## 🔧 Advanced Usage

### Workflow Examples

#### Complete Genome Annotation Pipeline

```bash
# Step 1: Initial detection with moderate parameters
java -Xms8g -Xmx32g -jar dist/Repeater2.jar genome.fasta kmer=20 min=50

# Step 2: Fine-grained SSR analysis
java -jar dist/Repeater2.jar genome.fasta -ssronly min=12 -seqshow flanks=500

# Step 3: Generate high-resolution visualization
java -jar dist/Repeater2.jar genome.fasta kmer=20 image=15000x500
```

#### Comparative Analysis Across Species

```bash
# Analyze multiple genomes with consistent parameters
for genome in genomes/*.fasta; do
    java -Xms16g -Xmx32g -jar dist/Repeater2.jar "$genome" kmer=20 sln=100
done
```

#### Masking for Gene Prediction

```bash
# Generate repeat-masked genome for gene finding
java -Xms8g -Xmx16g -jar dist/Repeater2.jar genome.fasta -maskonly kmer=18
```

### Integration with Other Tools

**Using REPEATER output with:**
- **BLAST**: Mask repeats before database searches
- **Gene Predictors**: Provide masked sequences to improve accuracy
- **Repeat Modeler**: Compare *de novo* predictions
- **Genome Browsers**: Load GFF3 files for visualization
- **R/Python**: Parse GFF3 for statistical analysis

### Performance Optimization

**Tips for faster analysis:**

1. **Adjust k-mer size**: Larger k-mers = faster but less sensitive
2. **Increase minimum length**: Filter out short repeats
3. **Use `-maskonly`**: Skip clustering for masking-only tasks
4. **Parallel processing**: Run multiple chromosomes separately
5. **Adequate memory**: Prevent disk swapping with sufficient RAM

---

## 📚 Citation

If you use REPEATER in your research, please cite:

```
Kalendar R, Kairov U. 2024. Genome-wide tool for sensitive de novo 
identification and visualisation of interspersed and tandem repeats. 
Bioinformatics and Biology Insights, 18: 1-11.
DOI: 10.1177/11779322241306391
```

**Full Text:**  
[https://journals.sagepub.com/doi/10.1177/11779322241306391](https://journals.sagepub.com/doi/10.1177/11779322241306391)

**BibTeX:**
```bibtex
@article{kalendar2024genome,
  title={Genome-wide tool for sensitive de novo identification and visualisation of interspersed and tandem repeats},
  author={Kalendar, Ruslan and Kairov, Ulykbek},
  journal={Bioinformatics and Biology Insights},
  volume={18},
  pages={1--11},
  year={2024},
  doi={10.1177/11779322241306391}
}
```

---

## 📧 Contact

**Author:** Ruslan Kalendar  
**Email:** ruslan.kalendar@helsinki.fi  
**Institution:** University of Helsinki  
**Online Tool:** [https://primerdigital.com/tools/repeats.html](https://primerdigital.com/tools/repeats.html)

### Support

For questions, bug reports, or feature requests:
1. Check the [Troubleshooting](#troubleshooting) section
2. Contact via email with detailed description of the issue
3. Include: input file format, command used, error message, system specs

---

## 🔍 Troubleshooting

### Common Issues and Solutions

#### Issue 1: "OutOfMemoryError: Java heap space"

**Cause:** Insufficient memory allocation for genome size

**Solution:**
```bash
# Increase heap memory (example: up to 64GB)
java -Xms16g -Xmx64g -jar dist/Repeater2.jar input.fasta
```

#### Issue 2: "Could not find or load main class"

**Cause:** Incorrect path to JAR file or Java not properly installed

**Solution:**
```bash
# Verify Java installation
java -version

# Use absolute path to JAR file
java -jar /full/path/to/Repeater2.jar input.fasta
```

#### Issue 3: Slow Performance

**Causes & Solutions:**

1. **Insufficient memory:**
   - Increase `-Xmx` value
   - Close other applications

2. **K-mer too small:**
   - Increase `kmer` value (try 20-25)

3. **Processing entire genome:**
   - Split into chromosomes
   - Process separately in parallel

#### Issue 4: No Output Generated

**Check:**
1. Write permissions in output directory
2. Disk space availability
3. Input file format is correct
4. Java version is 25+

**Debug command:**
```bash
# Add verbose output (if available)
java -jar dist/Repeater2.jar input.fasta -verbose
```

#### Issue 5: "Invalid or corrupt JAR file"

**Solution:**
- Re-download `Repeater2.jar`
- Verify file integrity
- Check for complete download

### System Requirements Check

**Verify your system meets requirements:**

```bash
# Check Java version
java -version

# Check available memory (Linux/Mac)
free -h

# Check available memory (Windows)
systeminfo | find "Available Physical Memory"

# Check disk space
df -h
```

### Getting Help

If issues persist:
1. Collect system information (OS, Java version, RAM)
2. Note the exact command used
3. Copy any error messages
4. Describe input file characteristics
5. Contact: ruslan.kalendar@helsinki.fi

---

## 📋 Quick Reference Card

### Most Common Commands

```bash
# Basic analysis
java -jar dist/Repeater2.jar sequence.fasta

# Large genome
java -Xms16g -Xmx32g -jar dist/Repeater2.jar genome.fasta kmer=20

# SSR detection
java -jar dist/Repeater2.jar sequence.fasta -ssronly -seqshow flanks=100

# Custom visualization
java -jar dist/Repeater2.jar sequence.fasta image=8000x400

# Masking only
java -jar dist/Repeater2.jar genome.fasta -maskonly kmer=18
```

### Parameter Quick Guide

| Task | Recommended Parameters |
|------|----------------------|
| Small genome analysis | `kmer=15 min=30` |
| Large genome analysis | `kmer=20 sln=100 -Xmx32g` |
| SSR detection | `-ssronly min=12 flanks=100` |
| Masking for BLAST | `-maskonly kmer=18` |
| Publication figure | `kmer=20 image=10000x500` |

---

## 📄 License

This project is open source. Please contact the author for specific licensing terms.

---

## 🙏 Acknowledgments

We thank the bioinformatics and genomics communities for their valuable feedback and contributions to the development of REPEATER.

---

<div align="center">

**REPEATER** — Advancing Genomic Repeat Analysis Through Computational Innovation

[Online Tool](https://primerdigital.com/tools/repeats.html) • [Publication](https://journals.sagepub.com/doi/10.1177/11779322241306391) • [Contact](mailto:ruslan.kalendar@helsinki.fi)

</div>
