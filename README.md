REPEATER v2 by Ruslan Kalendar
Genome-wide tool for sensitive de novo identification of interspersed and tandem repeats

To run the project from the command line, go to the dist folder and type the following:

java -jar Repeater2.jar <inputfile>  

Basic usage:

java -jar \Repeater2\dist\Repeater2.jar <inputfilepath> optional_commands

Examples:

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ef067844.fasta 

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ef067844.fasta kmer=21 min=100 quick=false mask=false seqshow=true



Common options:

kmer=5-any minimal kmer=5 (default kmer=12)

min=50 minimal repeat length (default min=50)

mask=true/false generate a new file with masking repeats (default mask=true)

seqshow=true/false extract repeat sequences (default seqshow=false)

flangs=true/false extracts sequences around the repeat (50 bp) (default flangs=false)

quick=true/false quick analysis of repeats, without deep analysis and their clustering (default quick=true)



