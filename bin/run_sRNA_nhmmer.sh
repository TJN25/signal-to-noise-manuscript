#!/bin/bash

# Setup Variables #

usage(){
    echo "run_sRNA_nhmmer.sh is a script for running nhmmer and sorting the results.  
Usage:
 run_sRNA_nhmmer.sh [opts] [input]

Required:	       
	-d	<database> The nucleotide database to be searched against as the nhmmer target
	-f	<folder> The folder containing the query file for nhmmer
	-e	<extension> The file extension of the query file in <folder>

Options:
	-h	Display this help
	-c	Check if the sequences are >50 nt (only works on fasta files)
"
}

check_seq_lengths=""
database=""
folder=""
extension=""
outfolder="."
missing=""
exitTrue="F"
evalue="1e-5"

# User Input Options #

while getopts "d:o:f:E:e:ch" arg; do
case $arg in
	d) 
	database=${OPTARG};;
	f)
	folder=${OPTARG};;
	o)
	outfolder=${OPTARG};;
	E)
	evalue=${OPTARG};;
	e)
	extension=${OPTARG};;	
	c)
	check_seq_lengths="T";;
    h)
		usage
		exit
      ;;    
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

# Tests for inputs #

if [[ $database == "" ]]; then
	missing="$missing -d <database> "
	exitTrue="T"
fi

if [[ $folder == "" ]]; then
	missing="$missing -f <folder> "
	exitTrue="T"	
fi

if [[ $extension == "" ]]; then
	missing="$missing -e <extension> "
	exitTrue="T"
fi

if [[ $exitTrue == "T" ]];then
	echo "$missing not found. Use -h for more help."
	exit
fi

if [[ $extension == "stk" ]]; then
	check_seq_lengths=""
fi

echo "extension value: $extension"

if [[ $check_seq_lengths == "T" ]]; then

	echo "checking for short seqs"
	mkdir -p $folder/short_seqs
	
	let "fileNum = 0"
	for file in $folder/*.$extension
	do
	
	seqlength=`grep -v ^">" $file | wc -c`
	if (( $seqlength < 50 )); then
	mv $file  single_seqs/short_seqs/
	fi
	
	done

fi

# Set up folders/files #

echo "making directories in `pwd`"
mkdir -p $outfolder/alignments
mkdir -p $outfolder/hmm
mkdir -p $outfolder/output

# Run nhmmer #

echo "running nhmmer"
let "fileNum = 0"
for file in $folder/*.$extension

do

echo $file
outname=`basename $file`


if [ -f "$outfolder/output/$outname.res" ]; then
	echo "Already exists: $file"
	continue
fi


nseqs=`esl-alistat $file | grep "Number of sequences" | cut -d ":" -f2`
length=`esl-alistat $file | grep "Alignment length:" | cut -d ":" -f2`
echo "Running nhmmer on $file (length: $length, nseqs: $nseqs)"


nhmmer -E $evalue --tblout $outfolder/output/$outname.tbl -A $outfolder/alignments/tmp.stk --tformat FASTA  $file $database > $outfolder/output/$outname.res

lines=`wc -l < $outfolder/alignments/tmp.stk`


if (( $lines > 0 )); then


##esl-weight!!!
esl-alimanip   --lnfract 0.8 --lxfract 1.2 --lmin 50 --lmax 500 --detrunc 30  $outfolder/alignments/tmp.stk | esl-alimask -g --gapthresh 0.8 -p --pfract 0.5 --pthresh 0.5 --keepins - > $outfolder/alignments/$outname



hmmbuild $outfolder/hmm/$outname.hmm $outfolder/alignments/$outname

else

echo "No hits found for $outname"

fi

done



