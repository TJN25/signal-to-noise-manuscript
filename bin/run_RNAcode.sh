#!/bin/bash

##makes alignments and running alifoldz and r-scape
##GCA accession number.

usage(){
    echo "sraAlignAndFold.sh is a script for making a multiple sequence alignment and
    getting secondary structure information.  
Usage:
 fetchGenomeGCA.sh [opts] [input]

Options:
	-h	Display this help

Input	       
	-r	Folder location

"
}

while getopts "i:o:h" arg; do
case $arg in
	i) 
	FOLDER=${OPTARG};;
	o) 
	OUTPUT=${OPTARG};;
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z ${FOLDER} ]; then
	FOLDER="./"
fi




mkdir -p "$FOLDER/rnacode_out"



let "fileNum = 0"
for file in alignments/*.stk; #change to _R for positive control
do

checkname=`basename $file .stk`
if [ -f "./rnacode_out/${checkname}.rnacode" ]; then
echo "Already exists: $file"
continue
else

	
nseqs=`esl-alistat $file | grep "Number of sequences" | cut -d ":" -f2`

length=`esl-alistat $file | grep "Alignment length:" | cut -d ":" -f2`
largest_length=`esl-alistat $file | grep "Largest" | cut -d ":" -f2`


if (( $length < 500 )); then

diffLength=`expr $largest_length - $length`

if (( $diffLength > $length ));then
echo "Alignment is poor: $file"
continue
fi

if (( $nseqs > 5000 )); then
echo "Skipping $file (length: $length, nseqs: $nseqs)"
echo "$file (length: $length, nseqs: $nseqs)" >> skipped_alignments.txt
continue
fi

if (( $nseqs < 3 )); then
echo "Skipping $file (length: $length, nseqs: $nseqs)"
#echo "$file (length: $length, nseqs: $nseqs)" >> skipped_alignments.txt
continue
fi

echo "Running RNAcode on $file (length: $length, nseqs: $nseqs)"
	
esl-reformat clustal $file | RNAcode --outfile ./rnacode_out/$checkname.rnacode

fi
fi

done
