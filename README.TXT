REPEATER v2 "Genome-wide tool for sensitive de novo identification of interspersed and tandem repeats"
by Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi
web: https://primerdigital.com/tools/


Availability and requirements:
Operating system(s): Platform independent
Programming language: Java
Other requirements: Java 22 or higher
Java Downloads:
https://www.oracle.com/java/technologies/downloads/

How do I set or change the Java path system variable:
https://www.java.com/en/download/help/path.html


To run the project from the command line, go to the target folder and type the following; an individual file or a file folder can be specified:

java -jar Repeater2.jar <inputfile>/<inputfolderpath> <optional_commands> 

Basic usage:

java -jar \Repeater2\dist\Repeater2.jar <inputfilepath> <optional_commands>

Examples:

java -jar \Repeater2\dist\Repeater2.jar \Repeater2\test\ef067844.txt

java -jar \Repeater2\dist\Repeater2.jar <inputfile> kmer=21 quick=false mask=true

java -jar \Repeater2\dist\Repeater2.jar <inputfile> ssr=true seqshow=true flanks=100

Large genome usage:

java -jar -Xms4g -Xmx32g \Repeater2\dist\Repeater2.jar <inputfolderpath> kmer=21 min=100 sln=550


Common options:
kmer=5	 minimal kmer=5 (default kmer=12)

min=50	 initial repeat length (default min=50)

sln=150	 string length (default sln=150)

image=5000x3000 (by default, the dimensionality of the image is automatically determined)

flanks=100	 extend the flanks of the repeat with an appropriate length (100 nt) (default flanks=0)

mask=true/false	 generate a new file with masking repeats (default mask=true)

seqshow=true/false	 extract repeat sequences (default seqshow=false)

quick=true/false	 quick analysis of repeats, without deep analysis and their clustering (default quick=true)

ssr=true	 analyzing only the SSR/telomers loci (default ssr=false) 


 