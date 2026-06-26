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
- **Visual Output**: Publication-quality repeat-landscape figures as both raster **PNG** and scalable vector **SVG**
- **Sequence Extraction**: Retrieve repeat sequences with flanking regions
- **Multiple Formats**: Supports FASTA and plain text input
- **GFF-style Report**: Tab-delimited output for downstream bioinformatics pipelines
- **SSR-Focused Mode**: Specialized analysis for simple sequence repeats
- **Configurable Output Folder**: Send all results to any directory (default: current directory)
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
   Results are saved to the **current working directory** by default. Use `out=<path>` to choose another folder (it is created automatically if it does not exist).

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

#### Example 6: Choose Output Folder and Image Format
```bash
java -jar dist/Repeater2.jar genome.fasta out=results/ format=svg
```

**Parameters used:**
- `out=results/`: Write all result files to `results/` (created if missing)
- `format=svg`: Produce only the scalable vector figure (`both` is the default; `none` skips images)

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
| `kmer=N`      | Integer | 18      | K-mer length for repeat detection (range: 12–32; values above 32 are capped) |
| `min=N`       | Integer | 30      | Minimum repeat length to report (bp); raised to `kmer` if smaller |
| `sln=N`       | Integer | 60      | Minimum repeat block length for clustering (bp); raised to `min` if smaller |
| `flanks=N`    | Integer | 0       | Length of flanking sequences to extract (bp; 0–1000) |
| `image=WxH`   | String  | Auto    | Output image dimensions (width×height in pixels; width ≤ 40000) |

### Output Options

| Option         | Type   | Default            | Description |
|----------------|--------|--------------------|-------------|
| `out=<path>`   | String | current directory  | Folder for all result files; created automatically if missing |
| `outdir=<path>`| String | current directory  | Alias for `out=` |
| `format=MODE`  | String | `both`             | Image output: `both`, `png`, `svg`, or `none` (skip image generation) |

### Analysis Modes

| Option       | Description |
|--------------|-------------|
| `-seqshow`   | Extract and output repeat sequences (added to the report's `Sequence` column) |
| `-ssronly`   | Analyze only SSR/microsatellite and telomeric loci |
| `-maskonly`  | Generate only the masked sequence output (skip clustering) |

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

### Where Results Are Written

All result files are written to the **current working directory** by default, or to the folder given by `out=<path>` / `outdir=<path>` (created automatically if it does not exist). File names are derived from the input file name:

- Single-sequence input → `<input>.gff`, `<input>.msk`, `<input>.png`, `<input>.svg`
- Multi-sequence FASTA → one set per record, suffixed `_1`, `_2`, … (e.g. `genome.fasta_1.gff`)

### Generated Files

| File          | Format            | Description |
|---------------|-------------------|-------------|
| `<name>.gff`  | Tab-delimited text | Repeat report (run parameters, coverage, and one row per repeat copy) |
| `<name>.msk`  | FASTA             | Soft-masked sequence — repeat regions in **lowercase**, the rest in **UPPERCASE** |
| `<name>.png`  | PNG image         | Repeat-landscape figure (raster). Written unless `format=svg`/`none` |
| `<name>.svg`  | SVG image         | Repeat-landscape figure (scalable vector). Written unless `format=png`/`none` |

### Report Format (`.gff`)

The report is a tab-delimited, GFF-style text file. It opens with `#` comment lines (run parameters and repeat coverage), followed by a column header line and one row per repeat copy.

**Column header:** `Seqid` · `Repeat` · `ClusterID` · `Start` · `Stop` · `Length` · `Strand` · `Phase` · `Sequence`

| Column | Name       | Description |
|--------|------------|-------------|
| 1      | Seqid      | Sequence identifier (taken from the FASTA header) |
| 2      | Repeat     | Source/type placeholder (always `.`) |
| 3      | ClusterID  | Cluster/family number grouping related repeat copies |
| 4      | Start      | Start position (1-based) |
| 5      | Stop       | End position (1-based, inclusive) |
| 6      | Length     | Repeat length (bp) |
| 7      | Strand     | `+` (forward) or `-` (reverse complement) |
| 8      | Phase      | Reserved placeholder |
| 9      | Sequence   | Repeat sequence — populated only with `-seqshow` (flanking regions in UPPERCASE when `flanks=` is set) |

### Example Report Output

```text
#REPEATER2 (2024) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi) https://github.com/rkalendar/Repeater
#kmer=18
#Minimal repeat=30
#Repeat filter=60
#Quick analysis is false
#Sequence length (bp)=4567
#Sequence coverage by repeats 97.94%
#Repeats search for: special-complex
#Time taken: 0 seconds

Seqid	Repeat	ClusterID	Start	Stop	Length	Strand	Phase	Sequence
special-complex	.	1	62	1028	967	+
special-complex	.	1	3584	4550	967	+
special-complex	.	2	1101	1737	637	+
```

> With `-seqshow`, the repeat nucleotide sequence is appended as the final column of each row.

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
1. Write permissions in the output directory (or the `out=` folder)
2. Disk space availability
3. Input file format is correct — a file containing non-ASCII bytes or no `>` FASTA header reports *"There is no sequence(s)."*
4. Java version is 25+

REPEATER prints progress to the console as it runs (target file, k-mer, coverage, and each file it saves). Re-run and read that log to see where it stopped:
```bash
java -jar dist/Repeater2.jar input.fasta out=results/
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

# Custom visualization into a results folder
java -jar dist/Repeater2.jar sequence.fasta image=8000x400 out=results/

# Vector figure only (SVG)
java -jar dist/Repeater2.jar sequence.fasta format=svg

# Masking only, no images
java -jar dist/Repeater2.jar genome.fasta -maskonly kmer=18 format=none
```

### Parameter Quick Guide

| Task | Recommended Parameters |
|------|----------------------|
| Small genome analysis | `kmer=15 min=30` |
| Large genome analysis | `kmer=20 sln=100 -Xmx32g` |
| SSR detection | `-ssronly min=12 flanks=100` |
| Masking for BLAST | `-maskonly kmer=18 format=none` |
| Publication figure (vector) | `kmer=20 image=10000x500 format=svg` |

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
