#!/bin/bash
set -e -o pipefail

module load star/2.5.2a

# STAR mapping / RSEM quantification pipeline
# usage: from an empty working directory, run
# ./STAR_RSEM.sh (read1) (read2 or "") (STARgenomeDir) (RSEMrefDir) (dataType) (nThreadsSTAR) (nThreadsRSEM)

# input: gzipped fastq file read1 [read2 for paired-end] 
#        STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script
read1=$1 #gzipped fastq file for read1
read2=$2 #gzipped fastq file for read1, use "" if single-end
STARgenomeDir=$3 
RSEMrefDir=$4
dataType=$5 # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE

nThreadsSTAR=$NSLOTS # number of threads for STAR

# output: all in the working directory, fixed names
# Aligned.out.bam                 # alignments, standard BAM, agreed upon formatting
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Quant.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
# Quant.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
# Quant.pdf                                     # RSEM diagnostic plots
# Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data
# Signal.{Unique,UniqueMultiple}.unstranded.bw  # 2 bigWig files for unstranded data

# executables
bedGraphToBigWig=bedGraphToBigWig              

TMPDIR="${TMPDIR:-/tmp}"

# STAR parameters: common
STARparCommon=" --genomeDir $STARgenomeDir  --readFilesIn $read1 $read2   --outSAMunmapped Within --outFilterType BySJout \
 --outSAMattributes NH HI AS NM MD    --outFilterMultimapNmax 20   --outFilterMismatchNmax 999   \
 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   \
 --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat   --twopassMode Basic"

# STAR parameters: run-time, controlled by DCC
#--genomeLoad LoadAndRemove doesn't work for some reason
#  --limitBAMsortRAM 10000000000
STARparRun=" --runThreadN $nThreadsSTAR --genomeLoad NoSharedMemory"

# STAR parameters: type of BAM output: quantification or sorted BAM or both
#     OPTION: sorted BAM output
## STARparBAM="--outSAMtype BAM SortedByCoordinate"
#     OPTION: transcritomic BAM for quantification
## STARparBAM="--outSAMtype None --quantMode TranscriptomeSAM"
#     OPTION: both
STARparBAM="--outSAMtype BAM Unsorted --quantMode TranscriptomeSAM"


# STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 

case "$dataType" in
str_SE|str_PE)
      #OPTION: stranded data
      STARparStrand=""
      STARparWig="--outWigStrand Stranded"
      ;;
      #OPTION: unstranded data
unstr_SE|unstr_PE)
      STARparStrand="--outSAMstrandField intronMotif"
      STARparWig="--outWigStrand Unstranded"
      ;;
esac

# STAR parameters: metadata
STARparsMeta="--outSAMheaderCommentFile commentsENCODElong.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate"

## not needed ## --outSAMheaderPG @PG ID:Samtools PN:Samtools CL:"$samtoolsCommand" PP:STAR VN:0.1.18"

# ENCODE metadata BAM comments
echo -e '@CO\tLIBID:ENCLB175ZZZ
@CO\tREFID:ENCFF001RGS
@CO\tANNID:gencode.v19.annotation.gtf.gz
@CO\tSPIKEID:ENCFF001RTP VN:Ambion-ERCC Mix, Cat no. 445670' > commentsENCODElong.txt

###### STAR command
echo STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta
STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta


#Process BAM file
#Sort separately -- STAR sorting lead to "NO_COOR reads not in a single block at the end" error upon index
samtools sort -l 1 -@ $NSLOTS -O bam -T $TMPDIR/ Aligned.out.bam -o Aligned.sortedByCoord.out.bam && rm -f Aligned.out.bam

#Add MC tag
java -Xmx20g -jar /cm/shared/apps/picard/1.140/picard.jar FixMateInformation INPUT=Aligned.sortedByCoord.out.bam OUTPUT=Aligned.sortedByCoord.out.bam.new VERBOSITY=ERROR QUIET=TRUE COMPRESSION_LEVEL=9 && mv Aligned.sortedByCoord.out.bam.new Aligned.sortedByCoord.out.bam

samtools index Aligned.sortedByCoord.out.bam


echo
echo "Genome BAM statistics"
samtools flagstat Aligned.sortedByCoord.out.bam


echo "Genome BAM Picard CollectRnaSeqMetrics"
#see https://broadinstitute.github.io/picard/picard-metric-definitions.html
java -jar /cm/shared/apps/picard/1.140/picard.jar CollectRnaSeqMetrics INPUT=Aligned.sortedByCoord.out.bam OUTPUT=rnaSeqMetrics.txt REF_FLAT=$STARgenomeDir/refFlat.txt.gz STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=$STARgenomeDir/rRNA.interval_list CHART_OUTPUT=coverageAlongTranscript.pdf


###### bedGraph generation, now decoupled from STAR alignment step
# working subdirectory for this STAR run
mkdir -p Signal

echo STAR --runMode inputAlignmentsFromBAM   --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
STAR --runMode inputAlignmentsFromBAM   --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr

# move the signal files from the subdirectory
mv Signal/Signal*bg .


###### bigWig conversion commands
# exclude spikeins
grep ^chr $STARgenomeDir/chrNameLength.txt > chrNL.txt

case "$dataType" in
str_SE|str_PE)
      # stranded data
      str[1]=-; str[2]=+;
      for istr in 1 2
      do
      for imult in Unique UniqueMultiple
      do
          grep ^chr Signal.$imult.str$istr.out.bg | sort-bed - > sig.tmp
          $bedGraphToBigWig sig.tmp  chrNL.txt Signal.$imult.strand${str[istr]}.bw
      done
      done
      ;;
unstr_SE|unstr_PE)
      # unstranded data
      for imult in Unique UniqueMultiple
      do
          grep ^chr Signal.$imult.str1.out.bg | sort-bed -  > sig.tmp
          $bedGraphToBigWig sig.tmp chrNL.txt  Signal.$imult.unstranded.bw
      done
      ;;
esac


/home/mauram01/scratch/hybridmice/rnaseq/src/scripts/RSEM.sh Aligned.toTranscriptome.out.bam $RSEMrefDir $dataType
