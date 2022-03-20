
#!/bin/bash

##Downloads the genome sequence (and gff file if the -g flag is used) when given a 
##GCA accession number.

usage(){
    echo "fetchGenomeGCA.sh is a script for downloading a genome (and GFF file) from a GCA accession.  
Usage:
 fetchGenomeGCA.sh [opts] [input]

Options:
	-h	Display this help

Input	       
	-r	Reference genome accession (required)
	-o	Output name
	-e Fasta file extension
   	-g include the GFF file

"
}

while getopts "r:o:e:gh" arg; do
case $arg in
	r) 
	GENOME=${OPTARG};;
	o) 
	OUTPUT=${OPTARG};;
	e) 
	EXTENSION=${OPTARG};;
	g)
      GFF='y'
      ;;  
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z ${GENOME} ]; then
    echo "Error: No input specified." >&2
    usage
    exit 1
fi

if [ -z ${OUTPUT} ]; then

OUTPUT=${GENOME}

fi

if [ -z ${EXTENSION} ]; then

EXTENSION="fna"

fi


AssemblyName=$(esearch -db assembly -query ${GENOME} | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyName)
refseqID=$(esearch -db assembly -query ${GENOME} | efetch -format docsum | xtract -pattern DocumentSummary -element RefSeq)

refseq1=$(echo $refseqID | head -c 7 | tail -c 3)
refseq2=$(echo $refseqID | head -c 10 | tail -c 3)
refseq3=$(echo $refseqID | head -c 13 | tail -c 3)


##get fasta file
if [ ! -f $OUTPUT.$EXTENSION ];then

fastaLink="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/$refseq1/$refseq2/$refseq3/$refseqID._$AssemblyName/$refseqID._$AssemblyName._genomic.fna.gz"

downloadLink=$(echo $fastaLink | sed 's/\._/_/g')

curl $downloadLink > $OUTPUT.$EXTENSION.gz 
sleep 1
gunzip $OUTPUT.$EXTENSION.gz 

if [ $? -eq 0 ]; then
    echo " "
else
    exit $?
fi


echo "$OUTPUT.$EXTENSION downloaded using $downloadLink"

else

fastaLink="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/$refseq1/$refseq2/$refseq3/$refseqID._$AssemblyName/$refseqID._$AssemblyName._genomic.fna.gz"

downloadLink=$(echo $fastaLink | sed 's/\._/_/g')

echo "$OUTPUT.$EXTENSION already downloaded. To download again use $downloadLink"


fi


##Get gff file
if [[ $GFF = 'y' ]]; then

if [ ! -f $OUTPUT.gff ];then

      
gffLink="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/$refseq1/$refseq2/$refseq3/$refseqID._$AssemblyName/$refseqID._$AssemblyName._genomic.gff.gz"

downloadLink=$(echo $gffLink | sed 's/\._/_/g')

curl $downloadLink > $OUTPUT.gff.gz 
sleep 1
gunzip $OUTPUT.gff.gz 

if [ $? -eq 0 ]; then
    echo " "
else
     exit $?
fi

echo "$OUTPUT.gff downloaded using $downloadLink"

else

gffLink="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/$refseq1/$refseq2/$refseq3/$refseqID._$AssemblyName/$refseqID._$AssemblyName._genomic.gff.gz"

downloadLink=$(echo $gffLink | sed 's/\._/_/g')

echo "$OUTPUT.gff already downloaded. To download again use $downloadLink"


fi

fi

