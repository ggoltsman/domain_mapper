Fabric Genomics - Bioinformatics Scientist Exercise B.1

This program maps protein coordinates to genomic coordinates for
various protein domains found in a gene. The software takes the following inputs:

1. A four-column (tab-separated) file containing information about the protein domains and the
associated genes. The first column is the gene symbol. The second column is the transcript
id. The third column is the genome build. The fourth column contains the protein domain(s)
found in that gene. The format of the fourth column is as follows

“[domain_name]:[aa_start]-[aa_end]"

Multiple domains should be separated by a semicolon, e.g.:

“[domain_name]:[aa_start]-[aa_end];[domain_name]:[aa_start]-[aa_end];"

Example:
"Zincfinger,RING-type:24-658;BRCA1,serine-richdomain:345-507;”



2. A feature file in gff3 format, corresponding to a specific Human build referenced by the domains fil (e.g. GRCh37),
and containing features matching the transcripts in the domains file.  





Output:  tab-delimited file with the following columns

● Column 1: The gene name.
● Column 2: The genome build.
● Column 3: The chromosome.
● Column 4: The protein domain name
● Column 5: The protein coordinates of the domain. (Format: start-end)
● Column 6: The amino acid length of the domain.
● Column 7: The genomic coordinates of the domain. (Format: start-end)
● Column 8: The length of the domain in genomic space.

*If multiple protein domains exist for a single gene, then each domain should be in a separate new entry/line.


