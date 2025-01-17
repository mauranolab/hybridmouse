#!/bin/bash
set -e # -o pipefail

alias bedmap='bedmap --ec --sweep-all'


# the function "round()" was taken from 
# http://stempell.com/2009/08/rechnen-in-bash/

# the round function:
round()
{
echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
};


getcolor () {
       case "$1" in
              B6129SF1J*)
                               trackcolor="238,54,36";; #red
              B6C3F1J*)
                               trackcolor="51,128,195";; #blue
              B6SPRETF1J*)
                               trackcolor="13,85,0";; #lt green
              B6CASTF1J*)
                               trackcolor="0,204,0";; #dk green
              B6PWKF1J)
                               trackcolor="120,88,165";; #purple
          *)
                               trackcolor="0,0,0";; #black
       esac
       echo $trackcolor
}


sample=$1
DS=$2
mappedgenome=$3

echo "Making tracks for sample $sample ($DS) against genome $mappedgenome"


#TMPDIR=`pwd`/tmp.makeTracks.$sample
#mkdir -p $TMPDIR
echo "using $TMPDIR as TMPDIR"


samflags="-q 20 -F 524"
#NB chrM being considered in most downstream analyses

date
echo "Making bed file"
samtools view $samflags $sample.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
      readlength = length($10); \
      insertlength = $9; \
#      tagSequence = $10; \
#      color=255; \
#      for(i=12; i<=NF; i++) { \
#            if (match($i, /NM:i:[0-9]+/)) { \
#                  editdistance = substr($i, RSTART+5, RLENGTH-5); \
#            } \
#      } \
      if (and($2, 16)) { \
              strand="-"; \
#            colorString = color ",0,0" ;  \
              chromStart=$4-2+readlength; \
      } else { \
              strand="+"; \
#            colorString = "0,0," color; \
              chromStart=$4-1; \
      } \
#      chromStart=$4-1; \
#      chromEnd=chromStart+readlength; \
      chromEnd=chromStart+1; \
#      print $3, chromStart, chromEnd, tagSequence, editdistance, strand, 0, 0, colorString ; \
      print $3, chromStart, chromEnd, insertlength ; \
}' \
| sort-bed - | 
starch - > $sample.tags.starch


echo
echo "Calculating tag counts"
echo "SAMtools statistics"
samtools flagstat $sample.bam > $TMPDIR/$sample.flagstat.txt
cat $TMPDIR/$sample.flagstat.txt

#BUGBUG breaks for DSall or encode reps
sequencedTags=`cat /vol/mauranolab/mauram01/hybridmice/dnase/inputs.txt /vol/isg/encode/mouseencode/mapped/inputs.txt /vol/mauranolab/mauram01/hybridmice/dnase/mES/inputs.txt | grep $DS | sort | uniq | xargs zcat | awk 'END {print NR/4}'`
PFalignments=`cat $TMPDIR/$sample.flagstat.txt | grep "in total" | awk '{print $1+$3}'`
#NB used to count both columns of flagstat output ($1 + $3) for remaining metrics but it gives no guarantee of 1 line per tag
uniqMappedTags=`cat $TMPDIR/$sample.flagstat.txt | grep "mapped (" | awk '{print $1}'`
numMappedTagsMitochondria=`samtools view -F 512 $sample.bam chrM | wc -l`
propMappedTagsMitochondria=`echo $numMappedTagsMitochondria/$uniqMappedTags*100 | bc -l -q`
#BUGBUG probably exclude chrM reads in duplicate counts
dupTags=`cat $TMPDIR/$sample.flagstat.txt | grep "duplicates" | awk '{print $1}'`
nonredundantTags=`unstarch --elements $sample.tags.starch`
pctPFalignments=`echo $PFalignments/$sequencedTags*100 | bc -l -q`
pctUniqMappedTags=`echo $uniqMappedTags/$sequencedTags*100 | bc -l -q`
pctDupTags=`echo $dupTags/$uniqMappedTags*100 | bc -l -q`
pctNonredundantTags=`echo $nonredundantTags/$sequencedTags*100 | bc -l -q`

#Tally how many reads were recovered from unpaired/SE reads (NB many of these may not even be PF, so are unrepresented)
PFalignmentsSE=`samtools view -F 1 $sample.bam | wc -l`
uniqMappedTagsSE=`samtools view -q 30 -F 1549 $sample.bam | wc -l`
pctUniqMappedTagsSE=`echo $uniqMappedTagsSE/$PFalignmentsSE*100 | bc -l -q`


echo
echo "Making density track"
#Make density track of number of overlapping tags per 150-bp window
#Normalizes the density to 1M tags, ignores enrichment for now
#Note the last awk statement makes the exact intervals conform to Richard's convention that the counts are reported in 20bp windows including tags +/-75 from the center of that window
#Remember wig is 1-indexed (groan)
cat /vol/isg/annotation/fasta/mm10/mm10.chrom.sizes | 
grep -v random | grep -v _hap | grep -v chrM |
awk '{OFS="\t"; $3=$2; $2=0; print}' | sort-bed - | cut -f1,3 | awk 'BEGIN {OFS="\t"} {for(i=0; i<=$2-150; i+=20) {print $1, i, i+150} }' | 
bedmap --bp-ovr 1 --echo --count - $sample.tags.starch | perl -pe 's/\|/\t\t/g;' | awk -F "\t" 'BEGIN {OFS="\t"} {$4="id-" NR; print}' |
awk -F "\t" 'BEGIN {OFS="\t"} {$2+=65; $3-=65; print}' |
awk -v nonredundantTags=$nonredundantTags -F "\t" 'BEGIN {OFS="\t"} {$5=$5/nonredundantTags*1000000; print}' |
tee $TMPDIR/$sample.density.bed |
awk 'lastChr!=$1 {print "fixedStep chrom=" $1 " start=" $2+1 " step=" $3-$2 " span=" $3-$2; lastChr=$1} {print $5}' > $TMPDIR/$sample.wig

starch $TMPDIR/$sample.density.bed > $sample.density.starch

if [[ "$mappedgenome" == "C57B6" ]]; then
       wigToBigWig $TMPDIR/$sample.wig /vol/isg/annotation/fasta/mm10/mm10.chrom.sizes $sample.bw
       
       trackcolor=$(getcolor $sample)

       nonredundantTagsM=`echo $nonredundantTags/1000000 | bc -l -q` 
       nonredundantTagsM=$(round $nonredundantTagsM 1)

       #can't find a way to force autoscale=off. http://genome.ucsc.edu/goldenPath/help/bigWig.html implies it's not a parameter in this context
       echo "track name=$sample description=\"$sample DNase Density (${nonredundantTagsM}M nonredundant tags; normalized to 1M)- BWA alignment\" maxHeightPixels=30 color=$trackcolor viewLimits=0:3 autoScale=off visibility=full type=bigWig bigDataUrl=https://cascade.isg.med.nyu.edu/~mauram01/analysis/scratch/hybridmice/dnase/$sample.bw"
fi

echo "Making cut count track"
samtools view $samflags $sample.bam | sam2bed --do-not-sort | 
awk '{if($6=="+"){s=$2; e=$2+1}else{s=$3; e=$3+1} print $1 "\t"s"\t"e"\tid\t1\t"$6 }' | sort-bed - | tee $TMPDIR/$sample.cuts.bed | 
bedops --chop 1 - | awk -F "\t" 'BEGIN {OFS="\t"} {$4="id-" NR; print}' > $TMPDIR/$sample.cuts.loc.bed
bedmap --delim '\t' --echo --count $TMPDIR/$sample.cuts.loc.bed $TMPDIR/$sample.cuts.bed | 
awk -v nonredundantTags=$nonredundantTags -F "\t" 'BEGIN {OFS="\t"} {$5=$5/nonredundantTags*100000000; print}' |
starch - > $sample.perBase.starch

if [[ "$mappedgenome" == "C57B6" ]]; then
       #Skip chrM since UCSC doesn't like the cut count to the right of the last bp in a chromosome
       gcat $sample.perBase.starch | cut -f1-3,5 | awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM"' > $TMPDIR/$sample.perBase.bedGraph
       
       bedGraphToBigWig $TMPDIR/$sample.perBase.bedGraph /vol/isg/annotation/fasta/mm10/mm10.chrom.sizes $sample.perBase.bw
       
       echo "track name=$sample description=\"$sample cut counts (${nonredundantTagsM}M nonredundant tags- BWA alignment\" maxHeightPixels=30 color=$trackcolor viewLimits=0:3 autoScale=off visibility=full type=bigWig bigDataUrl=https://cascade.isg.med.nyu.edu/~mauram01/analysis/scratch/hybridmice/dnase/$sample.perBase.bw"
fi


if [[ "$mappedgenome" == "C57B6" ]]; then
#BUGBUG hotspots don't work on mapped genomes without mappability file and pointing to right chromInfo.bed
	echo
	echo "Calling hotspots"
	base=`pwd`
	mkdir -p hotspots/$sample


	#Force creation of new density file (ours is normalized)
	hotspotDens=$base/hotspots/$sample/$sample.density.starch
	hotspotBAM=$TMPDIR/${sample}.bam
	#250M is too much, 150M takes 12-24 hrs
	if [ $uniqMappedTags -gt 100000000 ]; then
		echo "$uniqMappedTags uniquely mapped tags. Generating hotspots on subsample of 100M reads"
	
		sampleAsProportionOfUniqMappedTags=`echo 100000000/$uniqMappedTags | bc -l -q`
		samtools view $samflags -b -1 -@ $NSLOTS -s $sampleAsProportionOfUniqMappedTags $sample.bam > $hotspotBAM
	else
		echo "$uniqMappedTags uniquely mapped tags. Generating hotspots on all reads"
	
		sampleAsProportionOfUniqMappedTags=1
		hotspotBAM=$base/$sample.bam
	fi

	cd hotspots/$sample

	/vol/mauranolab/mauram01/hybridmice/dnase/src/callHotspots.sh $hotspotBAM $hotspotDens $base/hotspots/$sample > $base/hotspots/$sample.hotspots.log 2>&1

	cd ../..


	hotspotfile=hotspots/${sample}/${sample}-final/${sample}.fdr0.01.hot.bed
	if [[ "$mappedgenome" == "C57B6" ]]; then
		echo "Hotspots for UCSC browser"
		if [ -f "$hotspotfile" ]; then
			cut -f1-3 $hotspotfile > $TMPDIR/${sample}.fdr0.01.hot.bed
			bedToBigBed -type=bed3 $TMPDIR/${sample}.fdr0.01.hot.bed /vol/isg/annotation/fasta/mm10/mm10.chrom.sizes hotspots/${sample}/${sample}.fdr0.01.hot.bb
		else
			echo "Can't find $hotspotfile to make bigBed"
		fi
	
		peakfile=hotspots/${sample}/${sample}-final/${sample}.fdr0.01.pks.bed
		if [ -f "$peakfile" ]; then
			cut -f1-3 $peakfile > $TMPDIR/${sample}.fdr0.01.pks.bed
			bedToBigBed -type=bed3 $TMPDIR/${sample}.fdr0.01.pks.bed /vol/isg/annotation/fasta/mm10/mm10.chrom.sizes hotspots/${sample}/${sample}.fdr0.01.pks.bb
		else
			echo "Can't find $peakfile to make bigBed"
		fi
	fi


	spotout=hotspots/${sample}/$sample.spot.out

	#subsample for spot 
	if [ $uniqMappedTags -gt 10000000 ]; then
		echo
		echo "$uniqMappedTags uniquely mapped tags. Calculating SPOT score on subsample of 10M reads"
		sampleAsProportionOfUniqMappedTags=`echo 10000000/$uniqMappedTags | bc -l -q`
		samtools view $samflags -b -1 -@ $NSLOTS -s $sampleAsProportionOfUniqMappedTags $sample.bam > $TMPDIR/${sample}.10Mtags.bam
	
		mkdir -p $TMPDIR/$sample.hotspots.10Mtags
		cd $TMPDIR/$sample.hotspots.10Mtags
	
		#NB dens track doesn't exist
		/vol/mauranolab/mauram01/hybridmice/dnase/src/callHotspots.sh $TMPDIR/${sample}.10Mtags.bam $TMPDIR/${sample}.10Mtags.density.starch $TMPDIR/$sample.hotspots.10Mtags > $base/hotspots/$sample.hotspots.10Mtags.log 2>&1
	
		cd - #NB prints pwd
	
		spotout=$TMPDIR/${sample}.hotspots.10Mtags/${sample}.10Mtags.spot.out
	fi
fi


#Stats
echo
echo "*** Overall Stats ***"
echo
printf "Num_sequenced_tags\t$sequencedTags\t\t$sample\n"
printf "Num_pass_filter_alignments\t$PFalignments\t%.1f%%\t$sample\n" "$pctPFalignments"
printf "Num_uniquely_mapped_tags\t$uniqMappedTags\t%.1f%%\t$sample\n" "$pctUniqMappedTags"
printf "Num_mitochondria_tags\t$numMappedTagsMitochondria\t%.1f%%\t$sample\n" "$propMappedTagsMitochondria"
printf "Num_duplicate_tags\t$dupTags\t%.1f%%\t$sample\n" "$pctDupTags"
printf "Num_nonredundant_tags\t$nonredundantTags\t%.1f%%\t$sample\n" "$pctNonredundantTags"

#Don't have denominator of unpaired tags we tried to map, so don't compute % for first
printf "Num_SE_pass_filter_alignments\t$PFalignmentsSE\t\t$sample\n"
#NB denominator is PF tags, not tags sequenced
printf "Num_SE_uniquely_mapped_tags\t$uniqMappedTagsSE\t%.1f%%\t$sample\n" "$pctUniqMappedTagsSE"


if [ -f "$hotspotfile" ]; then
       numhotspots=`cat $hotspotfile | wc -l`
else
       numhotspots="NA"
fi
echo -e "Num_hotspots\t$numhotspots\t\t$sample"


if [ -f "$spotout" ]; then
       spot=`cat $spotout | awk 'BEGIN {OFS="\t"} NR==2 {print $3}'`
else
       spot="NA"
fi
echo -e "SPOT\t$spot\t\t$sample"


echo
echo "Histogram of mismatches"
samtools view $samflags $sample.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
       for(i=12; i<=NF; i++) { \
            if(match($i, /NM:i:/)) { \
                  print substr($i, 6); \
            } \
       } \
}' | sort -g | uniq -c | sort -k2,2 -g


#NB XA tags are computed for the unpaired tags, while MAPQ reflects the final PE location
echo
echo "Histogram of number of best alignments"
samtools view  $samflags $sample.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
       tag="NA"; \
       for(i=12; i<=NF; i++) { \
            if(match($i, /X0:i:/)) { \
                  tag=substr($i, 6); \
#                  if(tag>1) {print} \
            } \
#            if(match($i, /XA:Z:/)) { \
#                  tag=length(split(substr($i, 6), xa, ";")); \
#            } \
       } \
       print tag; \
}' | awk -F "\t" 'BEGIN {OFS="\t"} $0>=10 {$0="10+"} {print}' | sort -g | uniq -c | sort -k2,2 -g


if [[ "$mappedgenome" == "C57B6" ]]; then
	PEtags=`cat $TMPDIR/$sample.flagstat.txt | grep "paired in sequencing" | awk '{print $1+$3}'`
	if [ $PEtags -gt 0 ]; then
		echo
		echo "Template lengths"
		samtools view $samflags $sample.bam | awk -F "\t" 'BEGIN {OFS="\t"} $3!="chrM"' | awk -v sample=$sample -F "\t" 'BEGIN {OFS="\t"} $9>0 {print sample, $9}' | sort -k2,2n | tee $sample.insertlengths.txt | cut -f2 |
		awk -F "\t" 'BEGIN {OFS="\t"} NR==1 {print "Minimum: " $0} {cum+=$0; lengths[NR]=$0; if($0<125) {lastLineUnder125=NR}} END {print "Number of tags: " NR; print "25th percentile: " lengths[int(NR*0.25)]; print "50th percentile: " lengths[int(NR*0.5)]; print "75th percentile: " lengths[int(NR*0.75)]; print "95th percentile: " lengths[int(NR*0.95)]; print "99th percentile: " lengths[int(NR*0.99)]; print "Maximum: " $0; print "Mean: " cum/NR; print "Prop. tags under 125 bp: " lastLineUnder125/NR}'
		gzip -f $sample.insertlengths.txt
	fi

	echo
	echo "Tag lengths:"
	samtools view $samflags $sample.bam | awk -F "\t" 'BEGIN {OFS="\t"} {print length($10)}' | sort -g | uniq -c | sort -k1,1g 
	
	
       #Hack to deal with read names from SRA
       if samtools view $sample.bam | cut -f1 | head -10 | grep -v -q -e "^SRR"; then
              echo
              echo "Tag count by sequencing instrument"
              samtools view $samflags $sample.bam | cut -f1 | cut -d ":" -f1 | sort | uniq -c | sort -k1,1g
	fi
fi


echo
echo -e "\nDone!"
date
