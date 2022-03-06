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




mkdir -p "$FOLDER/alifold/post_script"



let "fileNum = 0"
for file in alignments_rnaalifold/*.stk


do

outname=`basename $file`
if [ -f "$FOLDER/alifold/$outname.alifold" ]; then
	#echo "$FOLDER/alifold/$file.alifold"
	echo "Already exists: $file"
	continue
else
	echo "Checking size: $file"
fi


lines=`wc -l < $file`
if (( $lines < 1 ));then
continue
fi


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



echo "Running alifoldz.pl on $file (length: $length, nseqs: $nseqs)"

time esl-reformat  clustal $file  | alifoldz.pl > ./alifold/$outname.alifold         
 

else

	echo "Skipping: $file"


fi

done







