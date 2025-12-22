## REPEATER
## A genome-wide tool for the rapid *de novo* identification, profiling, masking and visualisation of interspersed and tandem repetitive sequences.

## Author
Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi

## TotalRepeats online tool: 
https://primerdigital.com/tools/repeats.html

## Citation Reference
Kalendar R, Kairov U 2024. Genome-wide tool for sensitive *de novo* identification and visualisation of interspersed and tandem repeats. Bioinformatics and biology insights, 18: 1-11. DOI: 10.1177/11779322241306391
https://journals.sagepub.com/doi/10.1177/11779322241306391

## Availability and requirements:

Operating system(s): Platform independent

Programming language: Java 24 or higher

Java Downloads: https://www.oracle.com/java/technologies/downloads/

How do I set or change the Java path system variable: https://www.java.com/en/download/help/path.html

## Installing Java using Conda
To install a specific version of OpenJDK using Conda, you need to specify the version number in your installation command and use the conda-forge channel. The latest version is available on the conda-forge channel.
1. Add the conda-forge channel (if not already added). It is recommended to add the conda-forge channel to your configuration and set its priority to strict to ensure packages are preferentially installed from this channel:
   
```conda config --add channels conda-forge```

```conda config --set channel_priority strict```

2. Create a new Conda environment and install the desired OpenJDK version. Creating a dedicated environment helps manage dependencies and avoid conflicts with other projects:

```conda create -n java25 openjdk=25```

3. Activate the new environment:

```conda activate java25```

4. Check if you have Java installed. The output should display information for the installed Java version:

```java -version```

To run the project from the command line. Command-line options, separated by spaces. 
The executive file ```Repeater2.jar``` is in the ```dist``` directory, which can be copied to any location. 
Go to the target folder and type the following; an individual file or a file folder can be specified:

```java -jar Repeater2.jar <target_file_path/Folder_path>```


### Basic Usage:

```java -jar <Repeater2Path>\dist\Repeater2.jar <target_file_path> optional_commands```


### Example Commands:
```
java -jar <Repeater2Path>\dist\Repeater2.jar \test\4.txt  

# Run with k-mer size = 18, minimal repeat length = 90, and image output
java -jar <Repeater2Path>\dist\Repeater2.jar \test\ kmer=18 sln=90 image=5000x300 

# SSR/telomeric loci only, extract sequences with flanks
java -jar <Repeater2Path>\dist\Repeater2.jar \test\2.txt -ssronly -seqshow flanks=100

# Full genome analysis with parameters
java -jar C:\Repeater2\dist\Repeater2.jar D:\Genomes\Hydra_vulgaris\ kmer=20 sln=100 image=10000x300

```

### Large genome usage (you will have to show the program to use more RAM, for example as listed here, up to 64 Gb memory: -Xms16g -Xmx64g)
# For large genomes, increase Java heap memory allocation (up to 64 GB recommended):
```
java -jar -Xms16g -Xmx64g <\dist\Repeater2.jar> <input_Folder_path> kmer=20 sln=100

java -jar -Xms16g -Xmx64g C:\Repeater2\dist\Repeater2.jar E:\Genomes\Hordeum_marinum\
```

### To analyze all sequences in a folder:

```
java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ 

java -jar -Xms16g -Xmx32g C:\Repeater2\dist\Repeater2.jar E:\Genomes\GRCh38.p14\
```
### Common Options
| Option       | Description                                                |
| ------------ | ---------------------------------------------------------- |
| `kmer=19`    | k-mer length (default = 18, minimum = 12)                  |
| `min=30`     | Minimal repeat length (default = 30)                       |
| `sln=90`     | Minimal string length for clustering (default = 90)        |
| `flanks=100` | Extend flanks of repeats by specified length (default = 0) |
| `image=WxH`  | Output image dimensions (default = automatic)              |
| `-seqshow`   | Extract repeat sequences                  |
| `-ssronly`   | Analyze only SSR/telomeric loci            |
| `-maskonly`  | Only generate masked output (skip clustering)    |

## Sequence Entry:

Sequence data files are prepared using a text editor and saved in ASCII as text/plain format (.txt) or in .fasta or without file extensions (a file extension is not obligatory). The program takes a single sequence or accepts multiple DNA sequences in FASTA format. The template length is not limited.

## FASTA format description:
A sequence in FASTA format consists of the following:
One line starts with a ">" sign and a sequence identification code. A textual description of the sequence optionally follows it. Since it is not part of the official format description, software can ignore it when it is present.
One or more lines containing the sequence itself. A file in FASTA format may comprise more than one sequence.


## The output is saved in a GFF3 file: a nine-column, tab-delimited, plain text file. 
 
GFF format General Feature Format describes genes and other features associated with DNA, RNA, and Protein sequences. GFF lines have nine tab-separated fields:
Generic Feature Format Version 3 (GFF3) 
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
1. seqid - a cluster.
2. source - The program that generated this feature.
3. type - "repeat".
4. start - The starting position of the feature in the sequence. The first base is numbered 1.
5. stop - The ending position of the feature (inclusive).
6. score - length 
7. strand - Valid entries include '+', '-', or '.' (for those who don't know/care).
8. phase - If the feature is a coding exon, the frame should be a number between 0 and 2, representing the first base's reading frame. If the feature is not a coding exon, the value should be '.'
9. Sequence (if a parameter used: seqshow=true).
