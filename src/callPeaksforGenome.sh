#!/bin/bash

#Take in list of SRA and a GCA
#
#download the genome and gff
#
#loop through list
#	download SRA
#	map reads
#	remove CDS
#	peak calling
#combine gffs


FILE_PATH=`dirname $0`
number_of_sra="10"
output_path="./"
CPUS='6'
output_log=/dev/stdout
display_available_files="F"
use_sra_file="T"
while getopts "g:n:o:c:sqth" arg; do
  case $arg in
    g)
      gca=$OPTARG
      #echo $file1
      ;;
    n)
      number_of_sra=$OPTARG
      #echo $output_file
      ;;      
    o)
      output_path=$OPTARG
      #echo $output_file
      ;;
	c)
      CPUS=$OPTARG
      #echo $file1
      ;;  
    s)
      use_sra_file="F"
      ;;                
	q)
      output_log=$gca.log
      #echo $file1
      ;;
    t)
    display_available_files="T"
    ;; 
    h)
	echo 'PredVirusHost.sh: compares proteins from a number of contigs to hmm models and then scores'
	echo 'Version 3.0 2018'
	echo '# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
	echo 'usage: [options] <filename> <file type>'
	echo ' '
	echo 'Basic options:'
	echo '  -i : <filename> the protein fasta file with protein identifiers not containing spaces'
	exit      
      ;;
      
    esac
done    


if [[ -z $gca ]]; then
	echo 'Error: GCA needed. Specify with -g <gca>'
	echo ' '
	echo 'Use -h for more help.'
	echo ' '
	exit
fi


#mkdir -p $output_path 

if [[ $display_available_files == "T" ]]; then
	grep $gca ~/phd/RNASeq/SRA_bacteria_RNAseq.txt | grep "PAIRED" | grep "Illumina HiSeq"
	exit
fi



if [[ $use_sra_file == 'T' ]]; then      
	counts=`grep $gca ~/phd/RNASeq/SRA_bacteria_RNAseq.txt | grep "PAIRED" | grep "Illumina HiSeq" | wc -l`
	if (( $counts == 0 )); then
		echo "No valid RNAseq datasets for $gca"
		exit
	fi

	cd $output_path
	mkdir -p "$gca.data"
	cd "$gca.data"
	mkdir gff_files      
	echo "Output to $output_log"
	
	if (( $counts > $number_of_sra )); then
		grep $gca ~/phd/RNASeq/SRA_bacteria_RNAseq.txt | grep "PAIRED" | cut -f1 | head -n $number_of_sra > tmp1
	else
		grep $gca ~/phd/RNASeq/SRA_bacteria_RNAseq.txt | grep "PAIRED" | cut -f1 > tmp1
	fi
	
else

	counts=`cat experiments_list.txt | wc -l`
	
	if (( $counts == 0 )); then
		echo "No valid RNAseq datasets for $gca"
		exit
	fi
	
	
	
	cd $output_path
	mkdir -p "$gca.data"
	cd "$gca.data"
	mkdir gff_files
	cat ../experiments_list.txt > tmp1      
	echo "Output to $output_log"

fi

if [[ -f "${gca}.fna" ]]; then
	echo "$gca.fna already downloaded."	
else
	echo "Downloading $gca Genome and GFF files"
	fetch_genomes_from_GCA.sh -r $gca -g >> $output_log
fi
	


if [ $? -eq 0 ]; then
    echo " "
else
     echo "Error: Downloading $gca Genome and GFF files failed. See fetch_genomes_from_GCA.sh"
     exit $?
fi



file_lines=`cat tmp1`

for line in $file_lines ; 
do
	
	if [[ -f "${line}_sra_calls.gff" ]]; then
	
	echo "$line already downloaded."
	
	else
	
	CURRENT_DIR=`pwd`
	cp $gca.fna ~/phd/RNASeq/fastq_downloads
	cp $gca.gff ~/phd/RNASeq/fastq_downloads

	cd ~/phd/RNASeq/fastq_downloads
	echo "Downloading $line"
    fasterq-dump --split-3 -p $line >> $output_log
	echo "Mapping reads"
    sra2plot.1.0.3.sh -s $line -r $gca -d -n $CPUS  >> $output_log
    
    plot_lenegth=`wc -l $line.plot  | cut -d ' ' -f2`
    rm *.sam    
    	if [ $plot_lenegth -gt 0 ]; then
    	rm ${line}*fwd.plot
    	rm ${line}*.rev.plot
    	rm fastq/${line}*.fastq
    	rm trimmed/${line}*.fastq
    	fi
    rm /Users/thomasnicholson/ncbi/public/sra/*.cache
    mv *.plot $CURRENT_DIR
    cd $CURRENT_DIR
	echo "Removing CDS"
    removeProteinCodingRNA.R -f $line -g $gca >> $output_log
    echo "Calling Peaks"
    run_rnaPeakCalling.R -f $line  -g $gca >> $output_log
    
    fi
    cp ${line}_sra_calls.gff ./gff_files/
done

##Search for rfam models


##Using this function
rfamscan() { counts=$( bc -l <<< "scale=2;$(esl-seqstat $1.fna | grep ^"Total" | tr -s ' ' | cut -d ' ' -f4)*2/1000000"); cmscan -Z $counts  --cut_ga --rfam --nohmmonly --tblout $1.tblout --fmt 2 --clanin ~/Downloads/Rfam.clanin.txt ~/Downloads/Rfam.cm $1.fna; cmscanToGffWrapper.R -f $1.tblout -g $1;}

if [[ -f "${gca}_ncRNA.gff" ]]; then
	echo "${gca}_ncRNA.gff exists"
else
	echo "Running cmscan using rfam models"
	rfamscan $gca  >> $output_log
fi

cp $gca.gff ./gff_files/
cp ${gca}_ncRNA.gff ./gff_files

if [[ ! -f "${gca}_new_calls.txt" ]]; then
combine_gff_files.R -f ./gff_files/ -o $gca
fi


if [[ ! -f "plot_files" ]]; then
mkdir plot_files
fi

mv *.plot plot_files/
mv *.gff gff_files/
if [[ -f "${gca}.fna.sslist" ]]; then
rm "${gca}.fna.sslist"
fi

echo "Finished."
rm tmp1