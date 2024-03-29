## REPEATER v2
## Genome-wide tool for sensitive de novo identification of interspersed and tandem repeats.
by Ruslan Kalendar 

email: ruslan.kalendar@helsinki.fi

[Web](https://primerdigital.com/tools/)

## Availability and requirements:

Operating system(s): Platform independent

Programming language: Java 20 or higher

[Java Downloads](https://www.oracle.com/java/technologies/downloads/)


How do I set or change [the Java path system variable](https://www.java.com/en/download/help/path.html)


To run the project from the command line, go to the target folder and type the following; an individual file or a file folder can be specified:

```java -jar Repeater2.jar input_file_path/Folder_path```

### Basic usage:

```java -jar \Repeater2\dist\Repeater2.jar input_file_path optional_commands```


### Examples:
```
java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\4.txt

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ kmer=21 min=100 sln=250 image=5000x3000 quick=false mask=false seqshow=true

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\2.txt ssr=true seqshow=true flanks=100
```

### Large genome usage:
```
java -jar -Xms8g -Xmx32g \Repeater2\dist\Repeater2.jar input_Folder_path kmer=21 image=10000x5000
```

Analyzing all files in the folder:

```
java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ kmer=21 min=100 sln=300
```


**Common options:**

```
ssr=true analyzing only the SSR/telomers loci (default ssr=false)

kmer=5	 minimal kmer=5 (default kmer=12)

min=50	 initial repeat length (default min=50)

sln=150	 string length (default sln=100)

image=5000x3000 (by default, the dimensionality of the image is automatically determined)

flanks=100 extend the flanks of the repeat with an appropriate length (100 nt) (default flanks=0)

mask=true/false generate a new file with masking repeats (default mask=true)

seqshow=true/false extract repeat sequences (default seqshow=false)

quick=true/false quick analysis of repeats, without deep analysis and their clustering (default quick=true)

```

## Sequence Entry:

Sequence data files are prepared using a text editor and saved in ASCII as text/plain format (.txt) or in .fasta or without file extensions (a file extension is not obligatory). The program takes a single sequence or accepts multiple DNA sequences in FASTA format. The template length is not limited.

## FASTA format description:
A sequence in FASTA format consists of:
One line starts with a ">" sign and a sequence identification code. A textual description of the sequence optionally follows it. Since it is not part of the official format description, software can ignore it when it is present.
One or more lines containing the sequence itself. A file in FASTA format may comprise more than one sequence.



## The output is saved in GFF3 file: nine-column, tab-delimited, plain text file. 


