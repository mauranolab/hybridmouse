#!/bin/bash
set -e -o pipefail


module load star/2.5.2a
module load rsem/1.2.31


# prepare genomes for STAR and RSEM
# input parameters:
STARgenomeDir=$1 #STAR genome directory
RSEMgenomeDir=$2 #RSEM genome directory
fastaGenome=$3   #fasta file(s),  e.g. "male.hg19.fa"
fastaSpikeins=$4 #fasta file with spike-ins, e.g. "spikes.fixed.fasta"
gtf=$5           #all-inclusive gtf file "gencode.v19.annotation.gtf"

# example
# ./STAR_RSEM_prep.sh  /path/to/STARgenome  /path/to/RSEMgenome  male.hg19.fa  spikes.fixed.fasta gencode.v19.annotation_tRNA_spikeins.gtf


# STAR genome
mkdir -p $STARgenomeDir
STARcommand="STAR --runThreadN $NSLOTS --runMode genomeGenerate --genomeDir $STARgenomeDir --genomeFastaFiles $fastaGenome $fastaSpikeins --sjdbGTFfile $gtf --sjdbOverhang 100 --outFileNamePrefix $STARgenomeDir"
echo $STARcommand
$STARcommand


###Utility files for picard CollectRnaSeqMetrics

#Prep refFlat file
#UCSC genes
#http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz

#gencode
#https://www.snip2code.com/Snippet/77082/Convert-gene-annotations-from-GTF-to-gen
set +e
gtfToGenePred -genePredExt -allErrors -geneNameAsName2 $gtf $TMPDIR/refFlat.tmp.txt
paste <(cut -f 12 $TMPDIR/refFlat.tmp.txt) <(cut -f 1-10 $TMPDIR/refFlat.tmp.txt) | gzip -c > $STARgenomeDir/refFlat.txt.gz
set -e

###Prep rRNA coords from gencode

#https://gist.githubusercontent.com/slowkow/b11c28796508f03cdf4b/raw/38d337698ff1e6915578dfa08826c73631c3e0b5/make_rRNA.sh

gencodeBase=`basename $gtf .annotation.gtf`

rRNA=$STARgenomeDir/rRNA.interval_list

#Faster but depends on chrom.sizes 
#perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:mm10"' /vol/isg/annotation/fasta/mm10/mm10.chrom.sizes > $rRNA

#Easier, but needs bam file
#samtools view -H Aligned.toGenome.out.bam | grep "^@SQ" > $rRNA

#Slow but sure
cat $fastaGenome $fastaSpikeins | awk '$0 ~ ">" {if(NR>1) {print c}; c=0; printf "@SQ\tSN:"substr($0,2,100) "\tLN:"; } $0 !~ ">" {c+=length($0);} END { print c; }' > $rRNA


# Intervals for rRNA transcripts.
grep 'gene_type "rRNA"' $gtf | \
    awk '$3 == "transcript"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '
        /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' | \
    sort -k1V -k2n -k3n \
>> $rRNA


### RSEM

### the command below is for RSEM >=1.2.19
### note, that for RSEM < 1.2.19, --no-polyA should be added


mkdir -p $RSEMgenomeDir
RSEMcommand="rsem-prepare-reference --gtf $gtf $fastaGenome","$fastaSpikeins $RSEMgenomeDir/RSEMref"
echo $RSEMcommand
$RSEMcommand


echo "Done!!!"
date
