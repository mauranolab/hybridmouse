#!/bin/bash
set -e -o pipefail

#Aggregate and run RSEM


#source $MODULELOAD
#module load samtools/1.2
#module load gcc/4.7.2     # R dependent on this
#module load R/3.1.0

module load star/2.5.2a   # Just for densities
#module load bedops/2.4.14

module load rsem/1.2.31

STAMPIPES="/home/mauram01/scratch/hybridmice/rnaseq/src"

sample=$1
STARgenomeDir=$2
RSEMrefDir=$3
shift
shift
shift
BAM_FILES=$@
echo "Will aggregate $BAM_FILES for sample $sample"

STARdir=$STARgenomeDir
echo "Using $STARgenomeDir STAR and $RSEMrefDir RSEM ref dirs"

export TARGET_BAM=Aligned.toTranscriptome.out.bam
export GENOME_BAM=Aligned.toGenome.out.bam

numbam=$(wc -w <<< $BAM_FILES)
# Temporary
if [ ! -s "$TARGET_BAM" ] ; then
  if [ $numbam -eq 1 ] ; then
    echo "Copying Transcriptome BAM"
    cp "$BAM_FILES" "$TARGET_BAM"
  else
    echo "Merging Transcriptome BAM"
    samtools merge -@ $NSLOTS -n "$TARGET_BAM" $BAM_FILES
    #will be sorted later
  fi
fi

if [ ! -s "$GENOME_BAM" ] ; then
  echo "Making Genome BAM"
  GENOME_BAM_FILES=$(sed 's/toTranscriptome/sortedByCoord/g' <<< "$BAM_FILES")
  $STAMPIPES/scripts/tophat/merge_or_copy_bam.sh "$GENOME_BAM" $GENOME_BAM_FILES
  samtools index "$GENOME_BAM"
fi

echo
echo "Genome BAM statistics"
samtools flagstat "$GENOME_BAM" > Aligned.toGenome.flagstat.txt
cat Aligned.toGenome.flagstat.txt

samflags="-q 20 -F 524"

#NB relative to PF alignments instead of all reads like for DNase pipeline
PFalignments=`cat Aligned.toGenome.flagstat.txt | grep "in total" | awk '{print $1+$3}'`
#NB used to count both columns of flagstat output ($1 + $3) for remaining metrics but it gives no guarantee of 1 line per tag
uniqMappedTags=`cat Aligned.toGenome.flagstat.txt | grep "mapped (" | awk '{print $1}'`
numMappedTagsMitochondria=`samtools view -F 512 "$GENOME_BAM" chrM | wc -l`
propMappedTagsMitochondria=`echo $numMappedTagsMitochondria/$uniqMappedTags*100 | bc -l -q`
dupTags=`cat Aligned.toGenome.flagstat.txt | grep "duplicates" | awk '{print $1}'`
samflags="-q 20 -F 524"
nonredundantTags=`samtools view $samflags -c "$GENOME_BAM"`
pctUniqMappedTags=`echo $uniqMappedTags/$PFalignments*100 | bc -l -q`
pctDupTags=`echo $dupTags/$uniqMappedTags*100 | bc -l -q`
pctNonredundantTags=`echo $nonredundantTags/$PFalignments*100 | bc -l -q`


echo
echo "*** Overall Stats ***"
echo
#printf "Num_sequenced_tags\t$sequencedTags\t\t$sample\n"
printf "Num_pass_filter_alignments\t$PFalignments\t%.1f%%\t$sample\n" "0"
printf "Num_uniquely_mapped_tags\t$uniqMappedTags\t%.1f%%\t$sample\n" "$pctUniqMappedTags"
printf "Num_mitochondria_tags\t$numMappedTagsMitochondria\t%.1f%%\t$sample\n" "$propMappedTagsMitochondria"
printf "Num_duplicate_tags\t$dupTags\t%.1f%%\t$sample\n" "$pctDupTags"
printf "Num_nonredundant_tags\t$nonredundantTags\t%.1f%%\t$sample\n" "$pctNonredundantTags"


echo "Genome BAM Picard CollectRnaSeqMetrics"
#see https://broadinstitute.github.io/picard/picard-metric-definitions.html
java -jar /cm/shared/apps/picard/1.140/picard.jar CollectRnaSeqMetrics INPUT=$GENOME_BAM OUTPUT=rnaSeqMetrics.txt REF_FLAT=$STARgenomeDir/refFlat.txt.gz STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=$STARgenomeDir/rRNA.interval_list CHART_OUTPUT=coverageAlongTranscript.pdf


if [ ! -s "Signal.UniqueMultiple.str+.starch" ] ; then
  qsub -j y -cwd -V -N "AGG#${AGGREGATION_ID}.den" -S /bin/bash <<'__DEN__'
    # Write starch and bigwig to .tmp files
    function convertBedGraph(){
      in="$1"
      base="$2"
      chrom="$in.onlyChr.bg"
      grep '^chr' "$in" | sort-bed - > $chrom
      bedGraphToBigWig "$chrom" chrNL.txt "$base.bw.tmp"
      starch "$chrom" > "$base.starch.tmp"
    }

    set -x -e -o pipefail

    mkdir -p $TMPDIR/Signal

    echo STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir $TMPDIR/STAR
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir $TMPDIR/STAR

    STARdir="/vol/isg/annotation/STARindex/mm10all/"
    grep '^chr' $STARdir/chrNameLength.txt | sort -k1,1 > chrNL.txt

    convertBedGraph $TMPDIR/Signal/Signal.Unique.str1.out.bg         Signal.Unique.str-
    convertBedGraph $TMPDIR/Signal/Signal.Unique.str2.out.bg         Signal.Unique.str+
    convertBedGraph $TMPDIR/Signal/Signal.UniqueMultiple.str1.out.bg Signal.UniqueMultiple.str-
    convertBedGraph $TMPDIR/Signal/Signal.UniqueMultiple.str2.out.bg Signal.UniqueMultiple.str+

    for i in Signal*.tmp ; do
      mv $i ${i/.tmp/}
    done

__DEN__

fi

if [ ! -s "Quant.genes.results" ] ; then

  qsub -j y -cwd -V -N "AGG#${AGGREGATION_ID}.rsem" -pe threads 4 -S /bin/bash <<__RSEM__
    set -x

    dataType="str_PE"

    /home/mauram01/scratch/hybridmice/rnaseq/src/scripts/RSEM.sh $TARGET_BAM $RSEMrefDir/RSEMref \$dataType
__RSEM__

fi
