# Comparative RNA-Seq: Signal to noise

## Description 

This repo contains the code, processed data, and IDs used to generate the 
results in the manuscript "Comparative RNA-Seq: Signal to noise".

## Requirements for replication

### Required software

- R version 3.6.0<sup>1</sup>
  -    Tidyverse version 1.3.1
  -    gplots version 3.1.1
  -    ape version 5.6-1
  -    ggtree version 3.2.1
  -    treeio version 1.18.1
  -    ROSE version 0.0-4
  -    reshape2 version 1.4.4
  -    randomForest version 4.7-1
  -    ggpubr version 0.4.0
  -    UpSetR version 1.4.0
  -    GenomicRanges version 1.36.1
- Python version 3.6
  -   Biopython version 1.77
  -   numpy version 1.19.1
  -   pandas version 0.25.1
- HMMER version 3.1b2
- Infernal version 1.1.2
- RNAalifold version 2.4.14
- Alifoldz
- RNACode version 0.3
- R2R version 1.0.6
- Trimmomatic version 0.36
- sratoolkit version 2.9.2
- samtools version 1.6
- bowtie2 version 2.3.4.2
- R-scape version 1.2.3
- PHYLIP version 3.695

<small><sup>1</sup>A number of R functions are used from the functions.R file. 
This can be included with `source('functions.R')`</small>

### Data sources

- NCBI Refseq version 204
- SRA (IDs available in sra-ids.txt)
- Rfam version 14
- RMFam version 0.3

## Methodoloy

### Data selection

How to get the summary table for the genera for the RNA-seq data:

-   Search the [SRA](https://www.ncbi.nlm.nih.gov/sra/) section of ncbi with 
    the genus name and ['Organism'] tag.

-   Select "Send results to Run selector"

-   Select the "Metadata" button of the "Total" row and "Download" Column.

-   The genome accession for each added genome was looked up on 
    [NCBI](https://www.ncbi.nlm.nih.gov/) using the strain name in the file or 
    by search for the given experiment in the 
    [SRA](https://www.ncbi.nlm.nih.gov/sra/) section.

-   Create a folder for the genus e.g. `mkdir Escherichia`

-   Create a folder using the GenBank assembly 
    accession e.g. `mkdir GCA_002504285.1.data/`

-   Create an `experiments_list.txt` file listing all of the SRA experiments 
    for the given genome and save the file in the `genus` folder.
    
-   Run `callPeaksforGenome.sh -g <GCA Accession> -s` from within 
    the `genus` folder.

-   List of representative genomes selected from RefSeq204
    -   2 genomes from each genera were selected (where possible)
    
Predicting sRNAs:

-   The `callPeaksforGenome.sh` carried out a number of steps
    -   `fasterq-dump <SRA Accession>` downloaded each of the RNA-Seq files
        listed in the `experiments_list.txt` file.
    -   `fetch_genomes_from_GCA.sh` downloaded the nucleotide FASTA file and
        the annotation file in GFF3 format.
    -   Known ncRNAs are identified using `cmscan` with models from Rfam
    -   `sra2plot.1.0.3.sh` took these files and produced a plot file (showing
        read depths at each nucleotide position) for each of the experiments.
    -   `run_rnaPeakCalling.R` took these plot files and identfied regions of
        RNA expression tht occured in intergenic regions
    -   `combine_gff_files.R` took the expressed regions from each of the 
        individual files (including the annotation file) and merged regions
        where there was an overlap.

-   The resulting regions are separated into `Known sRNA` and `Putative RUFs`
    based on whether an expressed region overlapped with a previously 
    identified ncRNA

-   The negative control group (RINCs) was generated 
    using `get_random_srna_sequences.py`
      -   For each putative RUF, a region of the same size was selected at 
          random from an intergenic region in the same genome

```bash
for file in *.txt; do 
  accession=`basename $file _new_calls.txt`; 
  get_random_srna_sequences.py -a $accession; 
  done
```

- The nucleotide sequences for each group were obtained from the nucleotide
  FASTA files using the coordinates of the regions

### Analysis

Homology Search:

-   `run_sRNA_nhmmer.sh -d <group_seqs.fna> -f <path/to/sequences/> -e fna`
    -   The `-e` flag is for the the expected extension of the input files.
    -   This loops through all the individual sequences in each data set and
        carries out a similarity search against all the sequences in the data
        set.
    -   The matched sequences are removed from the remaining search list to 
        prevent multiple identical alignments being produced.
        
-   `run_sRNA_nhmmer.sh -d <path/to/genomes/> -f <path/to/initial_alignment> -o` 
    `<path_to_output>/ -e stk -E 1e-3`
    -   The `-E` flag sets the e-value threshold for inclusion.
    -   This gets repeated two times to generate deep alignments for the 
        sequences from each data set.
        
-   A check was carried out for each of the groups (known sRNAs, predicted 
    RUFs and RINCs) that none of the alignments shared sequences.
      -   RINCs sharing sequences with either of the other data sets were 
          removed.
      -   Predicted RUFs sharing sequences with the known sRNAs were removed.

```bash
group='predicted' #repeat for each group
cd ${group}/alignments/
> ../${group}_contig_pos.txt
for file in *.stk;
do 
  ID=`basename $file .stk`
  # extract the contig name and position for each sequence in each stockholm
  #alignment file
  grep "/" $file | grep -v "//" | cut -d ' ' -f1 | sed 's/\// /g' | 
  sed 's/-/ /g'| cut -d '|' -f2 | sort | uniq | sed -e "s/$/ $ID/" >> 
  ../${group}_contig_pos.txt
done

```

The output of this contains a list of sequences containing:
- Alignment ID
- Contig ID
- Sequence start and stop

This output was then processed to identify alignments with shared sequences 
between groups.

```r
#read in the data
ncdat <- read.table("negative_control_contig_pos.txt", sep = " ", fill = T)
pcdat <- read.table("positive_control_contig_pos.txt", fill = T)
preddat <- read.table("predicted_contig_pos.txt", fill = T)

#reformat data to use in contig overlaps
pcdat <- reformatContigPositionData(dat = pcdat)
ncdat <- reformatContigPositionData(dat = ncdat)
preddat <- reformatContigPositionData(dat = preddat)

nc.pc.Dat <- getOverlapIDs(queryDat =ncdat, targetDat = pcdat)
nc.pred.Dat <- getOverlapIDs(queryDat =ncdat, targetDat = preddat)
pc.pred.Dat <- getOverlapIDs(queryDat =preddat, targetDat = pcdat)

dat <- data.frame(ids = c(nc.pc.Dat$id1, nc.pred.Dat$id1))
write.table(dat, file = "not_negative_control_ids.txt", 
row.names = F, col.names = F, quote = F)
dat2 <- data.frame(ids =  pc.pred.Dat$id1)
write.table(dat2, file = "not_predicted_ids.txt", 
row.names = F, col.names = F, quote = F)
```

-   A similar check was carried out within each group to reduce redundancy
    -   Alignments sharing a sequence were merged into a single alignment.

```r
#read in the data
dat <- read.table("predicted_contig_pos.txt", sep = " ", fill = T, as.is = T)
dat <- reformatContigPositionData(dat = dat)

smallDat <- getOverlapIDs(queryDat=dat, targetDat=dat)

#Swap columns and combine so that every combination is present in both columns
s2 <- smallDat
s2$id1 <- smallDat$id2
s2$id2 <- smallDat$id1
smallDat <- smallDat %>% bind_rows(s2) %>% unique()
ids <- unique(dat$srna)

#Loop over the ids and run a check to see if there are any overlapping ids
#to speed this up check if an id has already been seen and skip it.
#Write a list of overlapping ids to a unqiue file
checked_ids <- c()
for(item in ids){
  if(item %in% checked_ids){
    next
  }
  current_ids <- groupOverlapItems(smallDat, item, current_ids = c())
  checked_ids <- c(checked_ids, current_ids)
    write.table(x = current_ids, 
                file = paste("combined_alignment_ids/", 
                             current_ids[1], 
                             "_combined_list.txt", sep = ""), 
                append = T, quote = F, row.names = F, col.names = F)
}

```

The output is a list of files where each file contains:
- The name of the first alignment found for each group
- A complete list of overlapping alignments by ID

For each of these alignments, all of the sequences are obtained and built into a 
new single alignment. 

```bash
#move into the folder where the lists of overlapping ids are kept
cd combined_alignments_ids/
mkdir -p fasta
mkdir -p hmm
mkdir -p alignments

#loop through the files that contain overlapping ids
for file in *.txt;
do
max_num="0"
outname=`basename $file _combined_list.txt`
#Check if there is only one ID listed in a file. If so there is no merge
#required.
line_count=`wc -l $file | cut -d ' ' -f1`
if (( $line_count == 1 ));then
cp ../large_alignments/alignments/$outname.stk ./alignments
continue
fi

#loop through the file and extract the sequences from each of the stk files
#that correspond to the ids that are overlapping
> fasta/$outname.fa
while read line; 
do 
esl-reformat fasta ../large_alignments/alignments/$line.stk >> fasta/$outname.fa
current_num=`esl-alistat ../large_alignments/alignments/$line.stk | \ 
grep "Number of sequences" | cut -d ' ' -f 4`

#check if the current alignment is the longest
if(( $current_num > $max_num ));
then
max_seq="$line.stk"
max_num="$current_num"
fi
done < $file

#build an hmm using the longest alignment
hmmbuild hmm/$outname.hmm ../large_alignments/alignments/$max_seq

#build an alignment using all the sequences from the overlapping ids
hmmalign --informat fasta hmm/$outname.hmm fasta/$outname.fa | \ 
esl-alimask -g --gapthresh 0.8 -p --pfract 0.5 --pthresh 0.5 - | \ 
esl-alimanip   --lnfract 0.6 --lxfract 1.4 --lmin 50 --lmax 500 --detrunc 50 \ 
- > alignments/$outname.stk
done
```

### File setup

-   There are some files that are consistently used throughout the analysis
  -   `genome_contig_pairs.txt` is a file that maps the genome `GCA_` ids to the 
      contig `NC` ids.
    -   this is necessary due to the alignment reporting contig names, but the 
        predicted RUFs belonging to genomes rather than contigs.
  -   `query_target_pairs.txt` consists of each alignment ID and the contig IDs
      of each sequence within the alignment. 
  
```bash
#loop through all the genomes in the analysis
cd representative_genomes/
> ../genome_contig_pairs.txt
for file in GC*.fna;  
do   
  ID=`echo $file | cut -d '.' -f1,2 | cut -d "_" -f1,2`;   
  grep ^">" $file | cut -d ' ' -f1 | sort | uniq | sed 's/>//g' | sed -e "s/$/ \ 
   $ID/" >> ../genome_contig_pairs.txt;   
done
```

```bash
#for each group loop through all the alignments and extract the contig IDs
#of all the sequences within each alignment
group="positive_control"
cd ${group}/alignments
> ../query_target_pairs.txt
for file in *.stk; 
do 
  ID=`echo $file | cut -d '.' -f1,2 | cut -d "_" -f1,2`; 
  ID_2=`echo $file | cut -d '.' -f1,2 `;  
  grep ^"#=GS" $file | sort | uniq | cut -d "/" -f1 | cut -d ' ' -f2 | rev | \ 
  cut -d '|' -f1 | rev | sed -e "s/$/   $ID \
    $ID_2/" >> ../query_target_pairs.txt;  
done
```


### Evolutionary Distance

-   A `cmscan` of all the genomes using the RF00177 Rfam model (bacterial 
      ribosomal SSU) was carried out.
-   [Phylip](https://evolution.genetics.washington.edu/phylip.html) was used to 
    build a maximum likelihood phyogenetic tree.
      -   The `RF00177.stk` alignment was reformatted to `.phylip` format
      -   `esl-reformat phylip RF00177.stk > RF00177.phylip`
      -   This was used as the input for `dnadist` with an output of 
          `RF00177.dnadist`.
          
-   The file was not formatted in a usable way for reading into R
      -   `cat  RF00177.dnadist | tr '\n' ' ' | sed 's/ N/\nN/g' > RF00177.dists`
      
```r
#open the dnadist file
dat <- read.table("RF00177.dists", sep = "",
                  header = F, fill = T, stringsAsFactors = F, as.is = T)
dat <- dat %>% filter(!is.na(V2))
colnames(dat)[1] <- "phylip_id"
colnames(dat) <- c("names", dat$phylip_id)

#reformat the data to use in subsequent steps
meltDat <- melt(dat)
colnames(meltDat) <- c("phylip_id", "query.name", "distance")

#open list of phylip ids (truncated) and contig ids
load(file = "current_ids.Rda")
current_ids <- current_ids %>% unique()

#open list of genome and contig ids
contig_labels <- read.table("genome_contig_pairs.txt")
colnames(contig_labels) <- c("id", "target.genome")

#combine the phylip ids with the genome ids
contig_labels <- contig_labels %>% 
  left_join(current_ids, by = "id") %>% 
  select(-b)

#combine the genome ids with the distances data for the first column
colnames(contig_labels) <- c("target.id", "target.genome", "phylip_id")
meltDat <- meltDat %>% left_join(contig_labels, by = 'phylip_id')

#combine the genome ids with the distances data for the second column
colnames(meltDat)[1:2] <- c("target.name", "phylip_id")
colnames(contig_labels) <- c("query.id", "query.genome", "phylip_id")
meltDat <- meltDat %>% left_join(contig_labels, by = 'phylip_id')
colnames(meltDat)[2] <- c("query.name")

#take the pairs of genomes and arrange based on distances (to account for
#multiple contigs). Sort these and select a single distance for each pair. 
meltDat <- meltDat %>%  group_by(target.genome, query.genome) %>% 
arrange(as.numeric(distance)) %>% mutate(row_num = row_number()) %>% 
filter(row_num == 1) %>% select(-row_num)
save(meltDat, file = "~/bin/R/r_files/distanceMelt.Rda")

#Reshape the data back to a matrix (not sure this was needed but helped for 
#visualising the problem)
mat <- reshape2::acast(data = meltDat %>% 
                         select(target.genome, query.genome, distance), 
                       formula = target.genome ~ query.genome)
distanceMat <- as.data.frame(mat)
save(distanceMat, file = "~/phd/RNASeq/distanceMat.Rda")

#include a save a load point as the previous steps were time consuming and
#this allows for troubleshooting the next steps
load("~/phd/RNASeq/distanceMat.Rda")

#get the maximum distances between pairs of sequences for each alignment.
max_dists_pred <- calculateMaximumDistances('predicted')
save(max_dists_pred, file = "~/bin/R/r_files/max_dists_pred.Rda")
max_dists_pc <- calculateMaximumDistances('pc')
save(max_dists_pc, file = "~/bin/R/r_files/max_dists_pc.Rda")
max_dists_nc <- calculateMaximumDistances('negative_control')
save(max_dists_nc, file = "~/bin/R/r_files/max_dists_nc.Rda")
```

A set of all distances `allDists.Rda` was obtained by taking the taxonomy 
information of all genomes in RefSeq204 from the assembly summary of the genomes, 
and combining this information with the phylogenetic distances calculated for each
genome in `distanceMat.Rda`.

These distances were used to determine the value where we could expect to see 
genomes from the same family. This was compared with the distribution of distances
for the predicted data, and a threshold of $distance < 0.15$.

### Figure S1 panel B

```r
#TODO make function
load("allDists.Rda")
maxsGenus <- allDists %>% group_by(Genus.x, Genus.y) %>% 
  summarise(max_dist = max(distance)) %>% 
  mutate(same.taxa = ifelse(Genus.x == Genus.y, T, F))
maxsFamily <- allDists %>% group_by(Family.x, Family.y) %>% 
  summarise(max_dist = max(distance)) %>% 
  mutate(same.taxa = ifelse(Family.x == Family.y, T, F))
maxsOrder <- allDists %>% group_by(Order.x, Order.y) %>% 
  summarise(max_dist = max(distance)) %>% 
  mutate(same.taxa = ifelse(Order.x == Order.y, T, F))
maxsClass <- allDists %>% group_by(Class.x, Class.y) %>% 
  summarise(max_dist = max(distance)) %>% 
  mutate(same.taxa = ifelse(Class.x == Class.y, T, F))
maxsPhylum <- allDists %>% group_by(Phylum.x, Phylum.y) %>% 
  summarise(max_dist = max(distance)) %>% 
  mutate(same.taxa = ifelse(Phylum.x == Phylum.y, T, F))
maxsKingdom <- allDists %>% group_by(Kingdom.x, Kingdom.y) %>% 
  summarise(max_dist = max(distance)) %>% 
  mutate(same.taxa = ifelse(Kingdom.x == Kingdom.y, T, F))
genus1 <- maxsGenus %>% ungroup() %>% filter(same.taxa, max_dist > 0) %>% 
select(max_dist) %>% mutate(group = "Genus")
family1 <- maxsFamily %>% ungroup() %>% filter(same.taxa, max_dist > 0) %>% 
select(max_dist) %>% mutate(group = "Family")
order1 <- maxsOrder %>% ungroup() %>% filter(same.taxa, max_dist > 0) %>% 
select(max_dist) %>% mutate(group = "Order")
class1 <- maxsClass %>% ungroup() %>% filter(same.taxa, max_dist > 0) %>% 
select(max_dist) %>% mutate(group = "Class")
phylum1 <- maxsPhylum %>% ungroup() %>% filter(same.taxa, max_dist > 0) %>% 
select(max_dist) %>% mutate(group = "Phylum")
kingdom1 <- maxsKingdom %>% ungroup() %>% filter(same.taxa, max_dist > 0) %>% 
select(max_dist) %>% mutate(group = "Kingdom")

maxs1 <- genus1 %>% bind_rows(family1, order1, class1, phylum1, kingdom1)
maxs1$group <- factor(maxs1$group, levels = c("Kingdom", 
                                              "Phylum", 
                                              "Class", 
                                              "Order", 
                                              "Family", 
                                              "Genus"))

##Figure S1 Panel B
p <- ggplot() +
geom_boxplot(data = maxs1 , aes(x = group, y = max_dist), fill = NA) + 
  coord_flip() +
  geom_hline(yintercept = 0.15)
p
```

RNAcode was used to check for protein coding regions using `run_RNAcode.sh` in
the directory containing the alignment files.

A summary file of output results for each alignment was generated

```bash

```

This was used to generate a list of alignments which were putative coding regions.
These alignments were removed and not used later steps.

```r

```

## Score metrics for RNAs

Within each group a `run_$TOOL.sh` script was used to loop over all the alignments
with the tool. The output from each `run_$TOOL.sh` script was processed into a 
single file (for each dataset) with a summary loop. The data was then processed 
so that it was all in the same format for later steps.

- RNAalifold
  -   `run_RNAalifold.sh`

```bash
group="positive_control"
cd ~/phd/RNASeq/srna_seqs/version_1/${group}/large_alignments/RNAALifold
> ../${group}_alifold_score.txt;
> ../${group}_mfe.txt;
> ../${group}_alifold_covariation.txt; 

for file in *.rnaalifold;
do   
  if [ $file == *"\.stk\.stk\.rnaalifold" ]; then 
    ID=`basename $file .stk.stk.rnaalifold`; 
  else 
    ID=`basename $file .stk.rnaalifold`; 
  fi; 
  grep "=" $file | rev | cut -d "(" -f1 | cut -d '{' -f1 | rev | cut -d "=" -f1 \ 
  | sed -e "s/$/   $ID/"  >> ../${group}_alifold_score.txt; 
  grep "=" $file | rev | cut -d "(" -f1 | grep -v '{'| rev | cut -d "=" -f2 | \
  cut -d "+" -f1 | sed -e "s/$/   $ID/" >> ../${group}_mfe.txt; 
  grep "=" $file | rev | cut -d "(" -f1 | grep -v '{' | rev | cut -d "=" -f2 \
  | cut -d "+" -f2 | tr -d ")" | sed -e "s/$/  \ 
  $ID/" >> ../${group}_alifold_covariation.txt; 
done
```

```r
predMFE <- read.table("chapter_4_files/predicted_mfe.txt")
pcMFE <- read.table("chapter_4_files/positive_control_mfe.txt")
ncMFE <- read.table("chapter_4_files/negative_control_mfe.txt")

colnames(pcMFE) <- c("mfe.score", "ID")
colnames(ncMFE) <- c("mfe.score", "ID")
colnames(predMFE) <- c("mfe.score", "ID")

save(pcMFE, file= "chapter_4_files/pcMFE.Rda")
save(ncMFE, file= "chapter_4_files/ncMFE.Rda")
save(predMFE, file= "chapter_4_files/predMFE.Rda")

```

```r
ncAliCov <- read.table("chapter_4_files/negative_control_alifold_covariation.txt")
pcAliCov <- read.table("chapter_4_files/positive_control_alifold_covariation.txt")
predAliCov <- read.table("chapter_4_files/predicted_alifold_covariation.txt")

colnames(ncAliCov) <- c("alifold_cov_score", "ID")
colnames(pcAliCov) <- c("alifold_cov_score", "ID")
colnames(predAliCov) <- c("alifold_cov_score", "ID")

save(pcAliCov, file = "chapter_4_files/pcCovAli.Rda")
save(ncAliCov, file = "chapter_4_files/ncCovAli.Rda")
save(predAliCov, file = "chapter_4_files/predAliCov.Rda")
```


- R-scape
  -   `run_R-scape.sh`
  
```bash
group="predicted"
cd ~/phd/RNASeq/srna_seqs/version_1/${group}/large_alignments/rscape_out/
> ../${group}.rscape.cov;
for file in *.sorted.cov;
do  
  tmpname=`echo $file | cut -d '_' -f2,3,4`; 
  ID=`basename $tmpname .stk.stk`; 
  cat $file | sed -e "s/$/	$ID/"  >> ../${group}.rscape.cov;  
done
```

```r 
pcCov<- rscapeCovarianceSetup("chapter_4_files/pc.cov")
ncCovRNA <- rscapeCovarianceSetup("chapter_4_files/nc.cov")
predCovRNA <- rscapeCovarianceSetup("chapter_4_files/pred.cov")

save(pcCov, file = "chapter_4_files/pcCovariation.Rda")
save(ncCovRNA, file = "chapter_4_files/ncCovariation.Rda")
save(predCovRNA, file = "chapter_4_files/predCovariation.Rda")
```

- AlifoldZ
  -   `run_AlifoldZ.sh`

```bash
group="negative_control"
cd ~/phd/RNASeq/srna_seqs/version_1/${group}/large_alignments/alifold 
> ../${group}.alifold

for file in *.alifold;  
do 
  ID=`echo $file | cut -d '.' -f1,2`; 
  grep -v "#" $file | grep -v "From" | grep -v "\-\-\-" | tr -s ' ' \
  | sed -e "s/$/   $ID/" >> ../${group}.alifold
done
```

```r
pcAlifold<- alifoldSetup("chapter_4_files/positive_control.alifold")
ncAlifold<- alifoldSetup("chapter_4_files/negative_control.alifold")
predAlifold<- alifoldSetup("chapter_4_files/predicted.alifold")

save(pcAlifold, file = "chapter_4_files/pcAlifold.Rda")
save(ncAlifold, file = "chapter_4_files/ncAlifold.Rda")
save(predAlifold, file = "chapter_4_files/predAlifold.Rda")
```

- RMFam
  -   `run_RMFam.sh`

```bash
group='predicted'
cd ~/phd/RNASeq/srna_seqs/version_1/${group}/large_alignments/alignments_rnaalifold/rmfam_gff_files
> ../../${group}.rmfam;  

for file in *.gff;    
do    
  tmpname=`echo $file | sed 's/\.clustal//g'`;  
  ID=`basename $tmpname .stk.stk.gff`;  
  grep -v "#" $file | sed -e "s/$/   $ID/" >> ../../${group}.rmfam;  
done
```

```r
pcMotif <- motifSetup("chapter_4_files/positive_control.rmfam")
ncMotif <- motifSetup("chapter_4_files/negative_control.rmfam")
predMotif <- motifSetup("chapter_4_files/predicted.rmfam")

save(pcMotif, file = "chapter_4_files/pcMotif.Rda")
save(ncMotif, file = "chapter_4_files/ncMotif.Rda")
save(predMotif, file = "chapter_4_files/predMotif.Rda")
```

- G+C content
  -   `run_GC.sh`

```bash
group="predicted"
cd ~/phd/RNASeq/srna_seqs/version_1/${group}/large_alignments/
> ./${group}_gc_reference.txt

for file in ./alignments_rnaalifold/*;
do
  ID=`basename $file | cut -d '.' -f1,2 | cut -d '_' -f2-`; 
  cat $file | grep "#=GC RF" | cut -d ' ' -f3-  | tr -d ' ' |  grep -o . | sort \
  | uniq -c | sed -e "s/$/ $ID/" >> ./${group}_gc_reference.txt
done
```

```r
predGC <- gcSetup(file_path = "chapter_4_files/predicted_gc_reference.txt")
pcGC <- gcSetup(file_path = "chapter_4_files/positive_control_gc_reference.txt")
ncGC <- gcSetup(file_path = "chapter_4_files/negative_control_gc_reference.txt")

save(pcGC, file= "chapter_4_files/pcGC.Rda")
save(ncGC, file= "chapter_4_files/ncGC.Rda")
save(predGC, file= "chapter_4_files/predGC.Rda")
```

- Evolutionary distance
  -   Used the output from the previous section calculating the distances.

- Transcription
  - Calculates read depths for all sequences in the alignments (where data was 
   available).
  - These differed from the original read depths due to changes in the alignment
  lengths that occured when overlaping sequences were combined, and alignments 
  were built and filtered.

```
#Something here
```

```r
pcRDepth <- readDepthsSetup(file_path = "chapter_4_files/positive_control_read_depths_summary.txt")
ncRDepth <- readDepthsSetup("chapter_4_files/negative_control_read_depths_summary.txt")
predRDepth <- readDepthsSetup("chapter_4_files/predicted_read_depths_summary.txt")

save(ncRDepth, file = "chapter_4_files/ncRDepth.Rda")
save(pcRDepth, file = "chapter_4_files/pcRDepth.Rda")
save(predRDepth, file = "chapter_4_files/predRDepth.Rda")
```

The data for each group was combined into a single file that could be used for 
used in futher steps. Missing data was assumed to be a null value and assigned
the lowest available score. While this is not ideal, the proportion of data
being lost otherwise was impacting the size of the datasets available for the 
random forest classifier.

```r
load("chapter_4_files/max_dists_pc.Rda")
load("chapter_4_files/max_dists_nc.Rda")
load("chapter_4_files/pcRDepth.Rda")
load("chapter_4_files/ncRDepth.Rda")
load("chapter_4_files/pcCovariation.Rda") 
load("chapter_4_files/ncCovariation.Rda")
load("chapter_4_files/pcGC.Rda")
load("chapter_4_files/ncGC.Rda")
load("chapter_4_files/pcAlifold.Rda")
load("chapter_4_files/ncAlifold.Rda")
load("chapter_4_files/pcMFE.Rda")
load("chapter_4_files/ncMFE.Rda")
load("chapter_4_files/pcMotif.Rda")
load("chapter_4_files/ncMotif.Rda")
load("chapter_4_files/pcCovAli.Rda")
load("chapter_4_files/ncCovAli.Rda")

pcDat <- pcMFE %>% 
  full_join(pcGC, by = "ID") %>% 
  full_join(max_dists_pc, by = "ID") %>% 
  full_join(pcRDepth, by = "ID") %>% 
  full_join(pcCov, by = "ID") %>% 
  full_join(pcMotif, by = "ID")%>% 
  full_join(pcAlifold, by = "ID") %>% 
  full_join(pcAliCov, by = "ID") %>% 
  mutate(group = "Positive Control") %>% 
  unique() 

ncDat <- ncMFE %>% 
  full_join(ncGC, by = "ID") %>% 
  full_join(max_dists_nc, by = "ID") %>% 
  full_join(ncRDepth, by = "ID") %>% 
  full_join(ncCovRNA, by = "ID") %>% 
  full_join(ncMotif, by = "ID")%>% 
  full_join(ncAlifold, by = "ID") %>% 
  full_join(ncAliCov, by = "ID") %>% 
  mutate(group = "Negative Control") 
ncDat <- ncDat %>% filter(reads.mean.score < 15)##ensure that the RINCs are not transcribed,  in order to compare untranscribed regions to transcribed sRNAs

dat <- pcDat %>% bind_rows(ncDat) %>% select(-cov.combined.score)
dat <- replaceNAs(dat = dat)
dat <- dat %>% 
  unique() %>% 
  mutate(motif.sum = motif.mean.score*motif_count)
save(dat, file = "chapter_4_files/randomForestDat.Rda")
```

The positive and nagative control data were then used to build a random forest 
classifier. At this stage, the filter for the evolutionary distance was applied,
and all sequences where $distance > 0.15$ were removed.

The classifier was trained on half the data, with the remaining data being kept
as a test set for validation of the classifier.

```r
load("chapter_4_files/randomForestDat.Rda")

#ensure only the feartures of interest are being used
dat <- dat %>% select(-ID) %>% unique() %>% 
  select(read.max.score, read.counts, distance, cov.min.eval, z_max, 
         alifold_cov_score, mfe.score, gc.score, srna.counts, motif.sum, group)

#select a random number for each alignment, as a comparison of the contribution
#of each feature when compared to chance.
set.seed(101)
randomNum <- runif(n = nrow(dat), min = 0, max = 1)
dat$random <- randomNum
randomForestDat <- dat %>% mutate(group = ifelse(group == "Positive Control", 1, 0))
randomForestDat$group <- as.factor(randomForestDat$group)
randomForestDat[is.na(randomForestDat)] <- 0

#Filter for distance, and remove the transcription data (as this is what is 
#being tested.)
randomForestDat <- randomForestDat %>% filter(distance < 0.15) %>% 
                    select(-read.max.score, -srna.counts, -read.counts)
colGroupNum <- match(x = "group", table = colnames(randomForestDat))

#Split the dataset into a training abnd validation dataset.
data_set_size <- floor(nrow(randomForestDat)/2)
indexes <- sample(1:nrow(randomForestDat), size = data_set_size)
training <- randomForestDat[indexes,]
validation <- randomForestDat[-indexes,]
save(indexes, file = "chapter_4_files/indexes_gz.Rda")
save(training, file = "chapter_4_files/training_gz.Rda")
save(validation, file = "chapter_4_files/validation_gz.Rda")

#build the classifier
rf_classifier = randomForest(group ~ ., data=training, ntree=1000, importance=TRUE)
rf_classifier
save(rf_classifier, file = "chapter_4_files/rf_classifier_gz.Rda")
```

The classifier was applied to the validation dataset and the predicted RUFs. 
These results were then combined into a dataset to be used in analysing the 
output of the classifier. 

The probability calculated for each ID was recorded, and data about the 
redundancy of each alignment was added into the dataset. This was done as the
purpose of the analysis was to determine what proportion of expression is 
indistinguishable from random intergenic regions. While redundancy in the data
would affect the quality of the classifier, correlation analysis and the 
usefullness of ROC curves, it is relevant to the final results.

```r
load("chapter_4_files/randomForestDat.Rda")
load("chapter_4_files/indexes_gz.Rda")
load("chapter_4_files/rf_classifier_gz.Rda")

##recreate validation where the reads data is kept
validation <- dat[-indexes,]

##keep reads and repeat previous filter steps
validation <- validation %>% mutate(motif.sum = motif.mean.score*motif_count) %>% filter(distance < 0.15) %>% select(read.max.score, distance, cov.min.eval, z_max, motif.max.score, alifold_cov_score, mfe.score, gc.score, motif.sum, group, ID) ##reads are still kept
set.seed(101)
randomNum <- runif(n = nrow(validation), min = 0, max = 1)
validation$random <- randomNum
colGroupNum <- match(x = "group", table = colnames(validation))
colIDNum <- match(x = "ID", table = colnames(validation))

#regenerate probability scores
prediction_for_roc_curve <- predict(rf_classifier,validation[,-c(colIDNum, colGroupNum)],type="prob")
validation$probability <- prediction_for_roc_curve[,2]

##get the number of times each ID was found expressed before redundancy was removed
redundacy_counts_pc <- read.table("chapter_4_files/pc_counts.txt")
colnames(redundacy_counts_pc) <- c("ID", "srna.counts.2")

validation <- validation %>% left_join(redundacy_counts_pc, by = "ID")
validation <- validation %>% select(read.max.score, distance, cov.min.eval, z_max, motif.max.score, alifold_cov_score, mfe.score, gc.score, srna.counts.2, motif.sum, group, ID, random, probability)

redundacy_counts_nc <- read.table("chapter_4_files/negative_control_counts.txt")
colnames(redundacy_counts_nc) <- c("srna.counts.3", "ID")

##adds the same counts data for the negative controls to validation that was added for the positive controls
validation <- validation %>% left_join(redundacy_counts_nc, by = "ID")
validation2 <- validation %>% mutate(srna.counts.2 = ifelse(is.na(srna.counts.2), srna.counts.3, srna.counts.2)) %>% select(-srna.counts.3)

save(validation2, file = "chapter_4_files/validation2.Rda")
```

```r
load("chapter_4_files/rf_classifier_gz.Rda")
load("chapter_4_files/randomForestDat.Rda")
load("chapter_4_files/max_dists_pred.Rda")
load("chapter_4_files/predRDepth.Rda")
load("chapter_4_files/predCovariation.Rda") 
load("chapter_4_files/predGC.Rda")
load("chapter_4_files/predAlifold.Rda")
load("chapter_4_files/predMFE.Rda")
load("chapter_4_files/predMotif.Rda")
load("chapter_4_files/predAlifoldScore.Rda")
load("chapter_4_files/predAliCov.Rda")

predSRNACounts <- read.table("chapter_4_files/predicted_snra_counts.txt")
colnames(predSRNACounts) <- c("srna.counts", "ID")

predDat <- predMFE %>%
  full_join(predGC, by = "ID") %>%
  full_join(max_dists_pred, by = "ID") %>%
  full_join(predRDepth, by = "ID") %>%
  full_join(predCovRNA, by = "ID") %>%
  full_join(predMotif, by = "ID")%>%
  full_join(predAlifold, by = "ID") %>%
  full_join(predAliCov, by = "ID") %>%
  full_join(predAlifoldScore, by = "ID") %>%
  full_join(predSRNACounts, by = "ID") %>%
  select(-cov.combined.score) %>% 
  mutate(group = "Predicted")

##indicates how many ids have been merged into the current id
redundacy_counts_pred <- read.table("chapter_4_files/predicted_counts.txt")
colnames(redundacy_counts_pred) <- c("srna.counts.2", "ID")

predDat <- predDat %>% left_join(redundacy_counts_pred, by = "ID")
predDat <- predDat %>% mutate(alifold_cov_score = as.numeric(alifold_cov_score))
predDat <- predDat %>% mutate(motif.sum = motif.mean.score*motif_count) %>% 
                      select(distance, cov.min.eval, z_max, motif.max.score, 
                              alifold_cov_score, mfe.score, gc.score, 
                              srna.counts.2, motif.sum, group, ID,  
                              read.max.score, read.counts)

set.seed(101)
randomNum <- runif(n = nrow(predDat), min = 0, max = 1)
predDat$random <- randomNum
predDat$mfe.score[is.na(predDat$mfe.score)] <- 0
predDat$gc.score[is.na(predDat$gc.score)] <- 50
predDat$distance[is.na(predDat$distance)] <- 0
predDat$cov.min.eval[is.na(predDat$cov.min.eval)] <- 10
predDat$motif.max.score[is.na(predDat$motif.max.score)] <- 0
predDat$motif.sum[is.na(predDat$motif.sum)] <- 0
predDat$z_max[is.na(predDat$z_max)] <- 10
predDat <- predDat[predDat$z_max != -Inf,]
predDat$alifold_cov_score[is.na(predDat$alifold_cov_score)] <- 0
predDat$srna.counts[is.na(predDat$srna.counts)] <- 1

colGroupNum <- match(x = "group", table = colnames(predDat))
colIDNum <- match(x = "ID", table = colnames(predDat))
colCountNum <- match(x = "srna.counts.2", table = colnames(predDat))

predRfDat <- predDat %>% select(distance, cov.min.eval, z_max, alifold_cov_score, 
                                mfe.score, gc.score, motif.sum, random)
prediction_for_predcited_data <- predict(rf_classifier,predRfDat,
                                          type = 'response')
prob_for_predcited_data <- predict(rf_classifier,predRfDat, type = 'prob')
predDat$probability <- prob_for_predcited_data[,2]
save(predDat, file = "chapter_4_files/predDat.Rda")

##filters for features that are used in the analysis
featuresSelected <- dat %>% mutate(motif.sum = motif.mean.score*motif_count) %>% 
                            unique() %>% 
                            mutate(motif.sum = motif.mean.score*motif_count) %>% 
                            select(read.max.score, distance, cov.min.eval, z_max, 
                                   motif.max.score, alifold_cov_score, mfe.score, 
                                   gc.score, motif.sum, group)

##adds the predicted data to this dataset so that all data is included.
featuresSelected <- predDat %>% select(read.max.score, distance, cov.min.eval, 
                                       z_max, motif.max.score, alifold_cov_score, 
                                       mfe.score, gc.score, srna.counts, 
                                       motif.sum, group) %>% 
                                bind_rows(featuresSelected) %>%
                                filter(!is.na(group))
save(featuresSelected, file="featuresSelected_gz.Rda")
```


## Results and figures

### UpSetR plot

An UpSet plot shows the co-occurance of elements across a number of different 
groups, similar to what might be shown in a Venn diagram, while allowing clear
visualisation of the relationships between groups, and a comparison of 
co-occurance sizes.

We have used upset plots, alongside a phylogenetic tree to visualise the 
proportion of alignments with sequences coming from genomes of differing 
genera.

### Figure 1
Three panels were produced to generate Figure 1.
- A phylogenetic tree showing the relationship between genera.
  -  This was calculated using maximum likelihood with 16s rRNA
- An upset plot for the known sRNAs.
- An upset plot for the predicted RUFs.  

```r
#ML phylogenetic tree built using 16s rRNA for the Gammaproteobacteria
#The full tree was also built for all bacterial genera, but is not displayed
tree <- read.tree(paste0(data_path, 'upsetr.tree'))
tbl_tree <- as_tibble(tree)
  
#list of extra Gammaproteobacteria that were 
#not included in the RNA-sed analysis
to_drop <- c("Azotobacter", "Marinobacter", "Pseudoalteromonas", "Agarivorans", 
             "Vibrio", "Alishewanella", "Aggregatibacter", "Mannheimia", 
             "Actinobacillus", "Xenorhabdus", "Providencia", "Proteus",
             "Pantoea", "Brenneria", "Lonsdalea", "Buchnera.", "Wigglesworthia", 
             "Sodalis", "Dickeya", "Citrobacter", "Plautiasymbiont", 
             "Shewanella", "Moritella", "Moraxella", "Psychrobacter", 
             "Methylomonas", "Cycloclasticus", "Methylococcus", "Francisella",
             "Pseudoxanthomonas", "Candidatus", "Plautia", "Methylophaga", 
             "Pasteurella", "Salinivibrio")
  
sub_tree <- drop.tip(tree, to_drop)
tbl_sub_tree <- as_tibble(sub_tree)
p <- ggtree(sub_tree) + 
  geom_tiplab(align = T) +
  xlim(0, 0.35)
p
```

```r
#matrix of bool values indicating for each RNA if it was found in a genus
load(paste0(data_path, "upsetSubsetPC.Rda"))
load(paste0(data_path, "upsetSubsetPredicted.Rda"))

#list of extra Gammaproteobacteria that were 
#not included in the RNA-sed analysis
genera_arranged <- c("Lysobacter", "Stenotrophomonas", "Xylella", "Xanthomonas", 
                     "Methylomicrobium", "Pseudomonas", "Acinetobacter", 
                     "Alteromonas", "Photorhabdus", "Yersinia", "Erwinia", 
                     "Edwardsiella", "Serratia", "Klebsiella", "Enterobacter", 
                     "Salmonella", "Shigella", "Escherichia")

upset_pred <- select_columns_by_list(upsetSubsetPredicted, genera_arranged)
upset_pc <- select_columns_by_list(upsetSubsetPC, genera_arranged)

UpSetR::upset(upset_pc, sets = colnames(upset_pc), 
              mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = T)
UpSetR::upset(upset_pred, sets = colnames(upset_pred), 
              mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = T)
```

### Correlation heatmap 

A spearman correlation matrix was produced comparing multiple metrics for each 
feature. The best performing metrics (measured by comparing correlation with
known sRNAs were kept).

### Figure 2

```r
load("~/bin/PhD/Chapter_4/chapter_4_files/randomForestDat.Rda")

dat <- dat %>% select(group, 
                      read.max.score, 
                      distance, 
                      alifold_cov_score,
                      cov.min.eval,
                      motif.sum,
                      z_max,
                      mfe.score,
                      gc.score) %>% 
  unique() %>% 
  mutate(cov.min.eval = -log(cov.min.eval), 
         z_max = -z_max,
         mfe.score = -mfe.score,
         alifold_cov_score = -alifold_cov_score)

set.seed(101)
matNames <- c("Read depth (max)", "Evolutionary distance", 
              "Alifold covariance score", "R-scape covariance score", 
              "Sum of motif scores", "Alifold Z score (max)", "MFE score", 
              "G+C Percentage", "Random", "Known sRNA or RINC")

matrixList <- calculateCorrelationMatrix(dat, matNames)
rhoMatrix <- matrixList$rhoMatrix
sigMatrix <- matrixList$sigMatrix

heatmap.2(rhoMatrix, cellnote=round(rhoMatrix, digits = 2),
          notecex=1.5,notecol="black", col=rev(redblue(40)), 
          density.info="none", trace="none", dendrogram=c("column"), 
          symm=F,symkey=T,symbreaks=T, scale="none", key.title = "", 
          srtRow=45, adjRow=c(0, 1), srtCol=45, adjCol=c(1,1), 
          breaks=(-20:20)/20, margins = c(8, 8), cexRow=1.5, 
          cexCol=1.5,font=1)

```

### Figure 3 panel A

ROC curves were produced for each of the features.

```r
load("~/bin/PhD_Chapters_Code/Chapter 3/chapter_4_files/ncRDepth.Rda")
load("~/bin/PhD_Chapters_Code/Chapter 3/chapter_4_files/pcRDepth.Rda")
pcRDepth <- pcRDepth %>% mutate(response = 1)
ncRDepth <- ncRDepth %>% mutate(response = 0)
rocData <- pcRDepth %>% bind_rows(ncRDepth)
roc.curve(response = rocData$response, predicted = rocData$read.max.score,
          main="ROC curve for Read Depths")
}
```

```r
load(file = "chapter_4_files/max_dists_pc.Rda")
load(file = "chapter_4_files//max_dists_nc.Rda")
max_dists_pc <- max_dists_pc  %>% mutate(response = 1)
max_dists_nc <- max_dists_nc  %>% mutate(response = 0)
rocData <- max_dists_pc%>% bind_rows(max_dists_nc)
roc.curve(response = rocData$response, predicted = rocData$distance,
          main="ROC curve for Maximum Phylogenetic Distance")
```

```r
load("chapter_4_files/pcCovariation.Rda")
load("chapter_4_files/ncCovariation.Rda")
pcCov <- pcCov %>% mutate(response = 1)
ncCov <- ncCovRNA%>% mutate(response = 0)
rocData <- pcCov %>% bind_rows(ncCov) %>% filter(!is.na(cov.min.eval))
roc.curve(response = rocData$response, predicted = rocData$cov.min.eval,
          main="ROC curve for Covariation Scores")
```

```r
load("chapter_4_files/pcGC.Rda")
load("chapter_4_files/ncGC.Rda")
pcGC <- pcGC %>% mutate(response = 1)
ncGC <- ncGC %>% mutate(response = 0)
rocData <- pcGC %>% bind_rows(ncGC)
roc.curve(response = rocData$response, predicted = rocData$gc.score,
          main="ROC curve for GC%")
```

```r
load("chapter_4_files/pcMFE.Rda")
load("chapter_4_files/ncMFE.Rda")
pcMFE <- pcMFE %>% mutate(response = 1)
ncMFE <- ncMFE %>% mutate(response = 0)
rocData <- pcMFE %>% bind_rows(ncMFE)
roc.curve(response = rocData$response, predicted = rocData$mfe.score,
          main="ROC curve for MFE")
```

```r
load("chapter_4_files/pcCovAli.Rda")
load("chapter_4_files/ncCovAli.Rda")
pcAliCov <- pcAliCov %>% mutate(response = 1)
ncAliCov <- ncAliCov %>% mutate(response = 0)
rocData <- pcAliCov %>% bind_rows(ncAliCov)
roc.curve(response = rocData$response, predicted = rocData$alifold_cov_score,
          main="ROC curve for Alifold covariation")
```

```r
load("chapter_4_files/pcMotif.Rda")
load("chapter_4_files/ncMotif.Rda")
pcMotif <- pcMotif %>% mutate(response = 1)
ncMotif <- ncMotif %>% mutate(response = 0)
rocData <- pcMotif %>% bind_rows(ncMotif) %>% mutate(sum.score = motif.mean.score*motif_count)
roc.curve(response = rocData$response, predicted = rocData$sum.score,
          main="ROC curve for sum of Motif scores")
```

```r
load("chapter_4_files/pcAlifold.Rda")
load("chapter_4_files/ncAlifold.Rda")
pcAlifold <- pcAlifold %>% mutate(response = 1)
ncAlifold <- ncAlifold %>% mutate(response = 0)
rocData <- pcAlifold %>% bind_rows(ncAlifold)
roc.curve(response = rocData$response, predicted = rocData$z_max,
          main="ROC curve for Z-score (Alifold)")
```

### Figure 3 panel B

The gini index from the random forest classifier shows the contribution of each
feature to the classifier.

```r
load("~/bin/PhD_Chapters_Code/Chapter 3/chapter_4_files/rf_classifier_gz.Rda")
varImpPlot(rf_classifier)
```

## Random forest interpretation figure 

  -   The output of the RF classifier for the test data and the predicted RUFs

### Figure 4

```r
load("~/bin/manuscript/data/predDat.Rda")
load("~/bin/PhD/Chapter_4/chapter_4_files/validatation2.Rda")
#select the desired columns from the predicted data and validatation data

probDat <- predDat %>% 
  select(probability, ID, group, srna.counts.2) %>% 
  bind_rows(validation2 %>% select(probability, ID, group, srna.counts.2))

#not sure if setting factors will break anything so using another data frame
plotDat <- probDat
plotDat$group <- factor(plotDat$group, 
                        levels = c('Positive Control', 
                                   'Predicted', 
                                   'Negative Control'))

#plot histogram of the probabilities
p <- ggplot() +
  geom_histogram(data = plotDat, 
                 aes(x = probability, 
                     y = ..density.., 
                     group = group, 
                     fill = group), 
                 binwidth = 0.02) +
  geom_vline(xintercept = 0.17) +
  geom_vline(xintercept = 0.5) 
p
```

  -   The predicted RUFs were split into <strong>transcriptional noise</strong>, 
      <strong>possible RNA</strong>, and <strong>putative RNA</strong>.
  -   The categories were then compared to see if there was a difference in the
      level of transcription between the groups.

The statistics for the results were calculated below.

```r
load("predDat.Rda")
load("validation2.Rda")

#select the desired columns from the predicted data and validation data
probDat <- predDat %>% 
  select(probability, ID, group, srna.counts.2) %>% 
  bind_rows(validation2 %>% 
              select(probability, ID, group, srna.counts.2))

probDat$group <- factor(probDat$group, 
                        levels = c('Positive Control', 
                                   'Predicted', 
                                   'Negative Control'))

scores <- scoreProbabities(probDat, 
                           threshold = 0.17, 
                           target_column = 'probability')
printListSubset(scores, 
                vec = c('ppv', 'fnr', 'pred_pos', 'pred_pct'), 
                startText = 'p > 0.17', round_val = 3)

scores <- scoreProbabities(probDat, 
                           threshold = 0.5, 
                           target_column = 'probability')
printListSubset(scores, 
                vec = c('ppv', 'fnr', 'pred_pos', 'pred_pct', 'pc_pct'), 
                startText = 'p > 0.5', round_val = 3)

scores <- scoreProbabities(probDat, 
                           threshold = 0.81, 
                           target_column = 'probability')
printListSubset(scores, 
                vec = c('ppv', 'fnr', 'pred_pos', 'pred_pct', 'pc_pct'), 
                startText = 'p > 0.81', round_val = 3)
```





