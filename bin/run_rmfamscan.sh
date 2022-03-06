#!/bin/bash

##Run rmfam over all fasta files

usage(){
    echo "run_rmfam_scan.sh is a script for running rmfam_scan over a set of fasta
    files.  
Usage:
 run_rmfam_scan.sh [opts] [extension]

Options:
	-h	Display this help

Input
	-e extension	       
	-a	align
	-o output
	-k keep tmp files

"
}

align="F"
keep="F"
while getopts "e:ao:hk" arg; do
case $arg in
	e)
		extension=${OPTARG}
		;;
	a) 
		align="T"
		;;
    h)
		usage
		exit
      ;;    
    k)
		keep="T"
      	;;   
	\?) 
	echo "Unknown option: -${OPTARG}" >&2; exit 1;;
    esac
done

COUNTER=0
mkdir -p "rmfam_cmscan"
mkdir -p "rmfam_gff_files"
mkdir -p "rmfam_tblout"


let "fileNum = 0"
if [[ $extension == "fna" ]]; then

	for file in *.$extension
	do
		if [ -f "rmfam_tblout/$file.out.tblout" ]; then
			echo "Already exists: $file"
 			continue
		fi
		length=`grep -v ">" $file | wc -c`
		if (( $length < 500 )); then
			echo "Running rmfam_scan on $file (length: $length)"
			time rmfam_scan.pl -g -f ~/phd/RNASeq/RMfam/scripts/RMfam.cm $file -o $file.out 
			mv *.tblout rmfam_tblout
			mv *.cmscan rmfam_cmscan
			mv *.gff rmfam_gff_files
		else
			echo "Skipping: $file"
		fi
	done
else


for file in *.$extension

do


if [ -f "rmfam_tblout/$file.out.tblout" ]; then
	echo "Already exists: $file"
 	continue
else
echo "Checking size: $file"
fi


lines=`wc -l < $file`
if (( $lines < 1));then
continue
fi

# COUNTER=$((COUNTER+1))
# 
# if (( $COUNTER > 100 )); then
# echo "waiting"
# COUNTER=1
# time wait
# 
# 
# 
# fi


# if [[ $align == "T" ]]; then

# echo "Aligning $file"
# esl-reformat  clustal $file > $file.clustal

nseqs=`grep "#=" $file | cut -d ' ' -f2 | sort | uniq | wc -l`


start=`grep "GCA" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f1`
end=`grep "GCA" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f2`

if [[ $start == "" ]]; then
	start=`grep "NC_" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f1`
	end=`grep "NC_" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f2`

fi

if [[ $start == "" ]]; then
	start=`grep "NZ_" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f1`
	end=`grep "NZ_" $file | head -n 1 | cut -d " " -f2 | cut -d "/" -f2 | cut -d "-" -f2`

fi

if [[ $start == "" ]]; then
	head $file
fi




length=`expr $end - $start`

if (( $length < 0 )); then
	length=$(( -1 * $length ))
fi

if (( $length < 500 )); then

ID=`grep "NC_" $file | head -n 1 | cut -d " " -f2`

if [[ $ID == "" ]]; then
	ID=`grep "NZ_" $file | head -n 1 | cut -d " " -f2`
fi

if [[ $ID == "" ]]; then
	ID=`grep GCA_" $file | head -n 1 | cut -d " " -f2`
fi

if [[ $ID == "" ]]; then
	head $file
else
alignmentLength=`grep $ID $file | grep -v "#" | tr -s ' ' | cut -d ' ' -f2 | wc -c`

diffLength=`expr $alignmentLength - $length`

if (( $diffLength > $length ));then
echo "Alignment is poor: $file"
continue
fi
fi




echo "Running rmfam_scan on $file (length: $length, nseqs: $nseqs)"

 time rmfam_scan.pl -g -f ~/phd/RNASeq/RMfam/scripts/RMfam.cm $file -o $file.out 
 
 
 mv *.tblout rmfam_tblout
mv *.cmscan rmfam_cmscan
mv *.gff rmfam_gff_files
 
else

	echo "Skipping: $file"


fi

done

fi

cd rmfam_tblout

let "fileNum = 0"

for file in *.tblout

do

sed 's/  /__/g' $file | sed 's/ /--/g' | sed 's/__/ /g' | tr -s ' ' | sed 's/ --/ /g' | sed 's/-- / /g' | sed 's/--!/ !/g' | tr ' ' '\t' | sed 's/--/_/g' > $file.formatted.tblout

done

