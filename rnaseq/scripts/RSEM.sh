#!/bin/bash
set -e -o pipefail

#Run after both mapping and aggregation


TARGET_BAM=$1
RSEMrefDir=$2
dataType=$3 # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE

nThreadsRSEM=$NSLOTS # number of threads for RSEM

TMPDIR="${TMPDIR:-/tmp}"

RSEM=rsem-calculate-expression        


######### RSEM

#### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
trBAMsortRAM=20G


mv $TARGET_BAM Tr.bam 

case "$dataType" in
str_SE|unstr_SE)
      # single-end data
      cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | sort -k1,1 -k2,2n -S $trBAMsortRAM -T $TMPDIR ) | samtools view -@ $nThreadsRSEM -bS - > $TARGET_BAM
      ;;
str_PE|unstr_PE)
      # paired-end data, merge mates into one line before sorting, and un-merge after sorting
      #Used to have sort -l 9 -n '-@' "$nThreadsRSEM" -f - at beginning
      cat \
        <( samtools view -H Tr.bam ) \
        <( samtools view Tr.bam \
            | awk '{printf $0 " "; getline; print}' \
            | sort -k1,1 -k2,2n -S "$trBAMsortRAM" -T "$TMPDIR" \
            | tr ' ' '\n'
        ) \
        | samtools view '-@' "$nThreadsRSEM" -bS - \
        > $TARGET_BAM
      ;;
*)
      echo "Error, dont recognize $dataType"
      exit 1
      ;;
esac

'rm' Tr.bam


echo "Transcriptome BAM statistics"
samtools flagstat $TARGET_BAM


echo "Run RSEM"
# RSEM parameters: common
RSEMparCommon="--bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345 --temporary-folder $TMPDIR"

# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun=" -p $nThreadsRSEM --ci-memory 30000 "

# RSEM parameters: data type dependent

case "$dataType" in
str_SE)
      #OPTION: stranded single end
      RSEMparType="--forward-prob 0"
      ;;
str_PE)
      #OPTION: stranded paired end
      RSEMparType="--paired-end --forward-prob 1.0"
      ;;
unstr_SE)
      #OPTION: unstranded single end
      RSEMparType=""
      ;;
unstr_PE)
      #OPTION: unstranded paired end
      RSEMparType="--paired-end"
      ;;
esac


###### RSEM command
(
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType $TARGET_BAM $RSEMrefDir Quant
     $RSEM $RSEMparCommon $RSEMparRun $RSEMparType $TARGET_BAM $RSEMrefDir Quant
) >& Log.rsem

###### RSEM diagnostic plot creation
# Notes:
# 1. rsem-plot-model requires R (and the Rscript executable)
# 2. This command produces the file Quant.pdf, which contains multiple plots
echo rsem-plot-model Quant Quant.pdf
rsem-plot-model Quant Quant.pdf
