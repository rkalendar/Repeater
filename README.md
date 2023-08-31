REPEATER v2 by Ruslan Kalendar
Genome-wide tool for sensitive de novo identification of interspersed and tandem repeats

To run the project from the command line, go to the dist folder and type the following:

java -jar Repeater2.jar <inputfile>  

Basic usage:

java -jar ..\Repeater2\dist\Repeater2.jar <inputfilepath> optional_commands

Examples:

java -jar Repeater2.jar d:\1.txt 

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ef067844.fasta kmer=21 min=100 

java -jar ..\Repeater2\dist\Repeater2.jar <inputfilepath> kmer=21 min=21 quick=false mask=false seqshow=true



Common options:
kmer= minimal kmer=5, by default 12

min= minimal repeat lenghth is kmer length, by default 50

mask= generate a new file with masking repeats, by default mask=true

seqshow= extract repeat sequences, by default seqshow=false

flangs= extracts sequences around the repeat (50 bp), by default flangs=false

quick= quick analysis of repeats, without deep analysis and their clustering, by default quick=true



