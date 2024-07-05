## REPEATER2
## Genome-wide tool for rapid de novo identification and visualisation of interspersed and tandem repeats.

## Author
Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi

[Web](http://primerdigital.com/tools/fastpcr.html)

## Availability and requirements:

Operating system(s): Platform independent

Programming language: Java 22 or higher

[Java Downloads](https://www.oracle.com/java/technologies/downloads/)


How do I set or change [the Java path system variable](https://www.java.com/en/download/help/path.html)


To run the project from the command line, go to the target folder and type the following; an individual file or a file folder can be specified:

```java -jar Repeater2.jar <target_file_path/Folder_path>```


### Basic usage:

```java -jar <Repeater2Path>\dist\Repeater2.jar <target_file_path> optional_commands```


### Examples:
```
java -jar C:\Repeater2\dist\Repeater2.jar D:\Genomes\Hydra_vulgaris\ kmer=21 min=25 image=5000x3000

java -jar <Repeater2Path>\dist\Repeater2.jar \test\4.txt

java -jar <Repeater2Path>\dist\Repeater2.jar \test\ kmer=21 min=100 sln=250 image=5000x3000 quick=false mask=false seqshow=true

java -jar <Repeater2Path>\dist\Repeater2.jar \test\2.txt ssr=true seqshow=true flanks=100
```

### Large genome usage:
```
java -jar -Xms16g -Xmx64g \dist\Repeater2.jar input_Folder_path kmer=21 image=10000x5000
```

Analysing all files in the folder:

```
java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ kmer=21 min=30
```


**Common options:**

```
ssr=true analysing only the SSR/telomers loci (default ssr=false)

kmer=5	 minimal kmer=5 (default kmer=12)

min=50	 initial repeat length (default min=25)

sln=100	 string length (default sln=100) used for quick=false

image=5000x3000 (by default, the dimensionality of the image is automatically determined)

flanks=100 extend the flanks of the repeat with an appropriate length (100 nt) (default flanks=0)

mask=true/false generate a new file with masking repeats (default mask=true)

gff=true/false generate a GFF file (default gff=true)

seqshow=true/false extract repeat sequences (default seqshow=false)

quick=true/false quick analysis of repeats, without deep analysis and their clustering (default quick=true)

```

## Sequence Entry:

Sequence data files are prepared using a text editor and saved in ASCII as text/plain format (.txt) or in .fasta or without file extensions (a file extension is not obligatory). The program takes a single sequence or accepts multiple DNA sequences in FASTA format. The template length is not limited.

## FASTA format description:
A sequence in FASTA format consists of the following:
One line starts with a ">" sign and a sequence identification code. A textual description of the sequence optionally follows it. Since it is not part of the official format description, software can ignore it when it is present.
One or more lines containing the sequence itself. A file in FASTA format may comprise more than one sequence.



## The output is saved in a GFF3 file: a nine-column, tab-delimited, plain text file. 
 
GFF format General Feature Format describes genes and other features associated with DNA, RNA and Protein sequences. GFF lines have nine tab-separated fields:
Generic Feature Format Version 3 (GFF3) 
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
1. seqid - a cluster.
2. source - The program that generated this feature.
3. type - "repeat".
4. start - The starting position of the feature in the sequence. The first base is numbered 1.
5. stop - The ending position of the feature (inclusive).
6. score - length 
7. strand - Valid entries include '+', '-', or '.' (for those who don't know/care).
8. phase—If the feature is a coding exon, the frame should be a number between 0 and 2, representing the first base's reading frame. If the feature is not a coding exon, the value should be '.'
9. attributes - All lines with the same group are linked into a single item.
 

