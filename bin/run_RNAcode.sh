#!/bin/bash

##makes alignments and running alifoldz and r-scape
##GCA accession number.

usage(){
    echo "sraAlignAndFold.sh is a script for making a multiple sequence alignment and
    getting secondary structure information.
Usage:
 fetchGenomeGCA.sh -i <input_folder>

Options:
-i  Input folder where the alignment are in stk format (default is the working
directory)
	-h	Display this help

"
}

while getopts "i:o:h" arg; do
case $arg in
	i)
	FOLDER=${OPTARG};;
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
for file in alignments/*.stk;
do

#Check if the file has already been created. Used for cases where the
#program crashed part way through
checkname=`basename $file .stk`
if [ -f "./rnacode_out/${checkname}.rnacode" ]; then
    echo "Already exists: $file"
    continue
else
    #get number of sequences and length of sequences
    nseqs=`esl-alistat $file | grep "Number of sequences" | cut -d ":" -f2`

    length=`esl-alistat $file | grep "Alignment length:" | cut -d ":" -f2`
    largest_length=`esl-alistat $file | grep "Largest" | cut -d ":" -f2`

    #check the sequences are <500 nt long
    if (( $length < 500 )); then
        diffLength=`expr $largest_length - $length`
        #check that the alignment is not poor
        if (( $diffLength > $length ));then
            echo "Alignment is poor: $file"
            continue
        fi

        echo "Running RNAcode on $file (length: $length, nseqs: $nseqs)"
        esl-reformat clustal $file | RNAcode --outfile ./rnacode_out/$checkname.rnacode

    fi
fi
done
