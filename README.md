REPEATER v2 by Ruslan Kalendar
Genome-wide tool for sensitive de novo identification of interspersed and tandem repeats

To run the project from the command line, go to the dist folder and type the following:

java -jar Repeater2.jar <inputfile>  

Basic usage:

java -jar \Repeater2\dist\Repeater2.jar <inputfilepath> optional_commands

Examples:

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ef067844.fasta 

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ef067844.fasta kmer=21 min=100 quick=false mask=false seqshow=true

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ef067844.fasta ssr=true flangs=100


Common options:

ssr=true analyzing only the SSR/telomers loci (default ssr=false)

kmer=9 minimal kmer=5 (default kmer=12)

min=50 minimal repeat length (default min=50)

flangs=100 extracts sequences around the repeat (100 nt) (default flangs=0)

mask=true/false generate a new file with masking repeats (default mask=true)

seqshow=true/false extract repeat sequences (default seqshow=false)

quick=true/false quick analysis of repeats, without deep analysis and their clustering (default quick=true)
