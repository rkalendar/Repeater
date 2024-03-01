REPEATER v2
"Genome-wide tool for sensitive de novo identification of interspersed and tandem repeats"
by Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi
web: https://primerdigital.com/tools/


To run the project from the command line, go to the dist folder and type the following:

java -jar Repeater2.jar inputfilepath/orFolder 

Basic usage:

java -jar \Repeater2\dist\Repeater2.jar inputfilepath optional_commands


Examples:

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ef067844.fasta 

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ef067844.fasta kmer=21 min=100 sln=250 image=10000x3000 quick=false mask=false seqshow=true

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ef067844.fasta ssr=true seqshow=true flanks=100

Large genome usage:

java -jar -Xms8g -Xmx32g \Repeater2\dist\Repeater2.jar inputfilepath kmer=21 image=20000x5000

Analyzing all files in the folder:

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ kmer=21 min=100 sln=300


Common options:

ssr=true analyzing only the SSR/telomers loci (default ssr=false)

kmer=5	 minimal kmer=5 (default kmer=12)

min=50	 initial repeat length (default min=50)

sln=150	 string length (default sln=100)

image=10000x3000 (by default, the dimensionality of the image is automatically determined)

flanks=100 extend the flanks of the repeat with an appropriate length (100 nt) (default flanks=0)

mask=true/false generate a new file with masking repeats (default mask=true)

seqshow=true/false extract repeat sequences (default seqshow=false)

quick=true/false quick analysis of repeats, without deep analysis and their clustering (default quick=true)
