Comparative RNA-Seq: Signal to noise

This repo contains the code, processed data, and IDs used to generate the results in the manuscript "Comparative RNA-Seq: Signal to noise".

1) The work consisted of downloading RNA-Seq experiments, 
   whole genomes, and annotation files.
2) These files were then processed to identify putative sRNAs.
3) The putative sRNAs were analysed with a series of tools
4) The outputs from this was analysed in a number of ways and 
  a classifier was built.


Steps 1-3 were largely carried out in bash (and a little in R) and relied on calling external programs

- HMMER
- Infernal
- RNAalifold
- Alifoldz
- RNACode
- R2R
- Trimmomatic
- sra-tools
- samtools
- curl
- bowtie2 
- R-scape
- PHYLIP

Data was obtained from:
- NCBI Refseq
- SRA (Sequence Read Archive)
- Rfam
- RMFam

Manual curation of the initial datasets was carried out, with checks for the 
quality of the data, identifying suitable clades, and combining of the data
sets.

