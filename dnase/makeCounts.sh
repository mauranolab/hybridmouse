#!/bin/bash
set -e -o pipefail

#in case this is run with an existing TMPDIR
mkdir -p $TMPDIR


bam=$1
strain=$2

sample=`basename $bam .bam` # | cut -d '.' -f1
mappedgenome=`echo $sample | perl -pe 's/^.+\.//g;'`


PREFIX=$sample

echo "Making counts for sample $sample ($strain mapped to $mappedgenome)"
date


if [[ "$mappedgenome" == "C57B6" ]]; then
       echo "Not supported -- specify non-ref .bam"
       exit 1
else
       refbam=`echo $bam | perl -pe "s/$mappedgenome\.bam\$/C57B6.bam/g;"`
       echo "Will find ref reads in $refbam"
fi
refsample=`basename $refbam .bam`


case "$strain" in
B6129SF1J)
       snpfilestrain=129S1;
       snpfilestrainlong=129S1_SvImJ;;
B6C3F1J)
       snpfilestrain=C3HHeJ;
       snpfilestrainlong=C3H_HeJ;;
B6CASTF1J)
       snpfilestrain=CASTEiJ;
       snpfilestrainlong=CAST_EiJ;;
B6PWKF1J)
       snpfilestrain=PWKPhJ;
       snpfilestrainlong=PWK_PhJ;;
B6SPRETF1J)
       snpfilestrain=SPRETEiJ;
       snpfilestrainlong=SPRET_EiJ;;
*)
       echo "Don't recognize strain $strain";
       exit 1;;
esac


snpfile=/vol/isg/annotation/bwaIndex/mm10all.$snpfilestrain.diploid/snps.$snpfilestrain.bed
echo "SNPs from $snpfile"

indelfile=/home/mauram01/scratch/hybridmice/genotyping/v5/mgp.v5.merged.indels.dbSNP142.normed.$snpfilestrainlong.filtered.mm10.starch
echo "Indels from $indelfile"


#Limit to union of hotspots across all strains for performance, note later filter in doMouse.R for hybridmouse mcv==0
#Unthresh merged hotspots
mergedrefsample=`echo $refsample | perl -pe 's/DS\d+[A-Z]\.C57B6$/DSall.C57B6/g;'`
#Wildcard to get files for all strains
mergedstrainrefsample=`echo $mergedrefsample | perl -pe 's/^B6[^\-]+-/*/g;'`

#original hotspot v1 files
#hotspotfile=/vol/mauranolab/mauram01/hybridmice/dnase/merged/hotspots/${mergedstrainrefsample}/${mergedstrainrefsample}-final/${mergedstrainrefsample}.hot.bed
#hotspot v2 files
hotspotfile=/vol/mauranolab/mauram01/hybridmice/dnase/stam_may2018/hotspots/merged/${mergedstrainrefsample}/${mergedstrainrefsample}.hotspots.fdr0.05.starch

echo "Hotspots from $hotspotfile"


mappabilityfile=/home/mauram01/scratch/hybridmice/dnase/mappability/$snpfilestrain.mappability.bed
echo "Mappabilityfile from $mappabilityfile"


#chrY is in v4 but not v3 for some reason
gcat $snpfile | awk -F "\t" 'BEGIN {OFS="\t"} $1 != "chrY"' |
#Exclude SNVs with more than 1 other SNV too close to others in same strain
bedmap --delim '\t' --range 36 --echo --count - $snpfile | 
awk -F "\t" 'BEGIN {OFS="\t"} $NF <= 2' |
bedops -e -1 - $hotspotfile |
#Since we're not using the indel column for now
bedops --range 72 -n - $indelfile |
#-f 3 keeps only paired-end reads
#-F 1548 to exclude dups, 524 to keep them
/vol/mauranolab/mauram01/hybridmice/dnase/src/countReads.py - $refbam $bam --sample $sample --strain $strain --minReads 1 --require 3 --filter 2828 --minMAPQ 20 --permittedMismatches 2 --maxTemplateLength 200 --minBaseQ 20 --cyclesToSkip 3 --maxReadLength 43 |
bedmap --delim '\t' --range 72 --echo --indicator - $indelfile | perl -pe's/1$/TRUE/;' -e's/0$/FALSE/;' |
closest-features --ec --delim '\t' --dist --closest - $indelfile | awk -F "\t" '{for(i=1; i<=NF-7; i++) {printf $i "\t"} printf $NF "\n"}' |
bedmap --delim '\t' --exact --echo --max - $mappabilityfile | perl -pe 's/\tNAN$/\t0/g;' > $TMPDIR/$PREFIX.counts.txt


awk -F "\t" -v widen=75 'BEGIN {OFS="\t"} {$2=$2 > widen ? $2-widen : 0; $3+=widen; print;}' $TMPDIR/$PREFIX.counts.txt | ~/bin/bed2fasta.pl - /vol/isg/annotation/fasta/mm10 2>/dev/null | grep -v -e "^>" | tr '[a-z]' '[A-Z]' | perl -ne 'chomp; print length($_) != 0 ? ($_ =~ tr/[gcGC]//) / length($_) : "NA"; print "\n"' | paste $TMPDIR/$PREFIX.counts.txt - > $PREFIX.counts.txt


numsites=`cat $PREFIX.counts.txt | wc -l`

echo
echo "Finished processing alignments -- $numsites SNPs"
date


if [[ "$numsites" -gt 0 ]]; then
       R --quiet --no-save << EOF

#Work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

filename <- "$PREFIX.counts.txt"
cat("Reading", filename, "\n")
data <- read(filename)
colnames(data) <- c("chrom", "chromStart", "chromEnd", "rs", "totalReads", "ref.allele", "ref.reads", "nonref.allele", "nonref.reads", "totalReadsInclFailed", "skippedReads", "distToNearestSNP", "indel", "distToNearestIndel", "mappable", "PercentGC")
data <- subset(data, totalReads>=8)
data\$pctRef <- with(data, ref.reads / totalReads)
data\$minReadsPerAllele <- apply(data[,c("ref.reads", "nonref.reads")], FUN=min, MARGIN=1)
data\$distToNearestSNP.bin <- cut(data\$distToNearestSNP, breaks=c(0, 36, 100, 250, 500, 1000, Inf), right=F, include.lowest=T)
data\$distToNearestIndel.bin <- cut(abs(data\$distToNearestIndel), breaks=c(0, 50, 72, 100, 250, 500, 1000, Inf), right=F, include.lowest=T)
data\$totalReads.bin <- cut(data\$totalReads, breaks=c(8,10,12,16,25,50,75,Inf), right=F, include.lowest=T)
data\$pctTotalReads <- data\$totalReads / data\$totalReadsInclFailed
data\$pctTotalReads.bin <- cut(data\$pctTotalReads, breaks=c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1), right=F, include.lowest=T)
data\$mappable.bin <- cut(data\$mappable, breaks=c(0, 0.01, 0.8, 0.95, 1), right=F, include.lowest=T)
data\$PercentGC.bin <- cut.pretty(data\$PercentGC, breaks=c(0, 0.4, 0.5, 0.6, 0.7, 1), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")


print.data.frame(cbind(which="totalReads", summaryBy(pctRef~totalReads.bin, FUN=list(mean, median, length), data=data)), row.names=F, right=F)
print.data.frame(summaryBy(pctRef~1, FUN=list(mean, median, length), data=subset(data)), row.names=F, right=F)


print.data.frame(cbind(which="pctTotalReads", summaryBy(pctRef~pctTotalReads.bin, FUN=list(mean, median, length), data=data)), row.names=F, right=F)
print.data.frame(cbind(which="distToNearestSNP", summaryBy(pctRef~distToNearestSNP.bin, FUN=list(mean, median, length), data=data)), row.names=F, right=F)
print.data.frame(cbind(which="distToNearestIndel", summaryBy(pctRef~distToNearestIndel.bin, FUN=list(mean, median, length), data=data)), row.names=F, right=F)
print.data.frame(cbind(which="mappable", summaryBy(pctRef~mappable.bin, FUN=list(mean, median, length), data=data)), row.names=F, right=F)


if(nrow(data) > 1000) {
       cat("Making plots\n")
       
       png(file="$PREFIX.png", width=800, height=800)
       print(densityplot(~pctRef , data=subset(data), n=1024, plot.points=F, xlab=list("Prop. of tags with reference allele", cex=1), ylab=list("Density", cex=1), scales=list(axs="i", x=list(at=seq.int(0,1,0.1)), cex=1, alternating=F, relation="free"), xlim=c(0,1), panel = function(...) {
              panel.abline(v = 0.4, col="darkgray", lty="dashed", lwd=1)
              panel.abline(v = 0.5, col="darkgray", lty="solid", lwd=1)
              panel.abline(v = 0.6, col="darkgray", lty="dashed", lwd=1)
              panel.densityplot(...)
       }, par.settings=list(superpose.line=list(lwd=2)), 
       main=paste("$PREFIX - $strain\n", nrow(data), "het sites")
       ))
       dev.off()

       png(file="$PREFIX.totalReads.png", width=800, height=800)
       print(densityplot(~pctRef, groups=totalReads.bin[drop=T], data=subset(data), n=1024, plot.points=F, auto.key=list(columns=1, corner=c(0, 1), x = 0.8, y = 0.99, cex=1), xlab=list("Prop. of tags with reference allele", cex=1), ylab=list("Density", cex=1), scales=list(axs="i", x=list(at=seq.int(0,1,0.1)), cex=1, alternating=F, relation="free"), xlim=c(0,1), panel = function(...) {
              panel.abline(v = 0.4, col="darkgray", lty="dashed", lwd=1)
              panel.abline(v = 0.5, col="darkgray", lty="solid", lwd=1)
              panel.abline(v = 0.6, col="darkgray", lty="dashed", lwd=1)
              panel.densityplot(...)
       }, par.settings=list(superpose.line=list(lwd=2)), 
       main=paste("$PREFIX - $strain\n", nrow(data), "het sites - totalReads")
       ))
       dev.off()

       png(file="$PREFIX.pctTotalReads.png", width=800, height=800)
       print(densityplot(~pctRef, groups=pctTotalReads.bin[drop=T], data=subset(data), n=1024, plot.points=F, auto.key=list(columns=1, corner=c(0, 1), x = 0.8, y = 0.99, cex=1), xlab=list("Prop. of tags with reference allele", cex=1), ylab=list("Density", cex=1), scales=list(axs="i", x=list(at=seq.int(0,1,0.1)), cex=1, alternating=F, relation="free"), xlim=c(0,1), panel = function(...) {
              panel.abline(v = 0.4, col="darkgray", lty="dashed", lwd=1)
              panel.abline(v = 0.5, col="darkgray", lty="solid", lwd=1)
              panel.abline(v = 0.6, col="darkgray", lty="dashed", lwd=1)
              panel.densityplot(...)
       }, par.settings=list(superpose.line=list(lwd=2)), 
       main=paste("$PREFIX - $strain\n", nrow(data), "het sites - pctTotalReads")
       ))
       dev.off()

       png(file="$PREFIX.PercentGC.png", width=800, height=800)
       print(densityplot(~pctRef, groups=PercentGC.bin[drop=T], data=subset(data), n=1024, plot.points=F, auto.key=list(columns=1, corner=c(0, 1), x = 0.8, y = 0.99, cex=1), xlab=list("Prop. of tags with reference allele", cex=1), ylab=list("Density", cex=1), scales=list(axs="i", x=list(at=seq.int(0,1,0.1)), cex=1, alternating=F, relation="free"), xlim=c(0,1), panel = function(...) {
              panel.abline(v = 0.4, col="darkgray", lty="dashed", lwd=1)
              panel.abline(v = 0.5, col="darkgray", lty="solid", lwd=1)
              panel.abline(v = 0.6, col="darkgray", lty="dashed", lwd=1)
              panel.densityplot(...)
       }, par.settings=list(superpose.line=list(lwd=2)), 
       main=paste("$PREFIX - $strain\n", nrow(data), "het sites - PercentGC")
       ))
       dev.off()
       
       png(file="$PREFIX.distToNearestSNP.png", width=800, height=800)
       print(densityplot(~pctRef, groups=distToNearestSNP.bin[drop=T], data=subset(data), n=1024, plot.points=F, auto.key=list(columns=1, corner=c(0, 1), x = 0.8, y = 0.99, cex=1), xlab=list("Prop. of tags with reference allele", cex=1), ylab=list("Density", cex=1), scales=list(axs="i", x=list(at=seq.int(0,1,0.1)), cex=1, alternating=F, relation="free"), xlim=c(0,1), panel = function(...) {
              panel.abline(v = 0.4, col="darkgray", lty="dashed", lwd=1)
              panel.abline(v = 0.5, col="darkgray", lty="solid", lwd=1)
              panel.abline(v = 0.6, col="darkgray", lty="dashed", lwd=1)
              panel.densityplot(...)
       }, par.settings=list(superpose.line=list(lwd=2)), 
       main=paste("$PREFIX - $strain\n", nrow(data), "het sites - distToNearestSNP")
       ))
       dev.off()

       png(file="$PREFIX.distToNearestIndel.png", width=800, height=800)
       print(densityplot(~pctRef, groups=distToNearestIndel.bin[drop=T], data=subset(data), n=1024, plot.points=F, auto.key=list(columns=1, corner=c(0, 1), x = 0.8, y = 0.99, cex=1), xlab=list("Prop. of tags with reference allele", cex=1), ylab=list("Density", cex=1), scales=list(axs="i", x=list(at=seq.int(0,1,0.1)), cex=1, alternating=F, relation="free"), xlim=c(0,1), panel = function(...) {
              panel.abline(v = 0.4, col="darkgray", lty="dashed", lwd=1)
              panel.abline(v = 0.5, col="darkgray", lty="solid", lwd=1)
              panel.abline(v = 0.6, col="darkgray", lty="dashed", lwd=1)
              panel.densityplot(...)
       }, par.settings=list(superpose.line=list(lwd=2)), 
       main=paste("$PREFIX - $strain\n", nrow(data), "het sites - distToNearestIndel")
       ))
       dev.off()

       png(file="$PREFIX.mappable.png", width=800, height=800)
       print(densityplot(~pctRef, groups=factor(mappable>0.9), data=subset(data), n=1024, plot.points=F, auto.key=list(columns=1, corner=c(0, 1), x = 0.8, y = 0.99, cex=1), xlab=list("Prop. of tags with reference allele", cex=1), ylab=list("Density", cex=1), scales=list(axs="i", x=list(at=seq.int(0,1,0.1)), cex=1, alternating=F, relation="free"), xlim=c(0,1), panel = function(...) {
              panel.abline(v = 0.4, col="darkgray", lty="dashed", lwd=1)
              panel.abline(v = 0.5, col="darkgray", lty="solid", lwd=1)
              panel.abline(v = 0.6, col="darkgray", lty="dashed", lwd=1)
              panel.densityplot(...)
       }, par.settings=list(superpose.line=list(lwd=2)), 
       main=paste("$PREFIX - $strain\n", nrow(data), "het sites - mappable")
       ))
       dev.off()

       png(file="$PREFIX.mappability.png", width=800, height=800)
       print(densityplot(~pctRef, groups=mappable.bin, data=subset(data), n=1024, plot.points=F, auto.key=list(columns=1, corner=c(0, 1), x = 0.8, y = 0.99, cex=1), xlab=list("Prop. of tags with reference allele", cex=1), ylab=list("Density", cex=1), scales=list(axs="i", x=list(at=seq.int(0,1,0.1)), cex=1, alternating=F, relation="free"), xlim=c(0,1), panel = function(...) {
              panel.abline(v = 0.4, col="darkgray", lty="dashed", lwd=1)
              panel.abline(v = 0.5, col="darkgray", lty="solid", lwd=1)
              panel.abline(v = 0.6, col="darkgray", lty="dashed", lwd=1)
              panel.densityplot(...)
       }, par.settings=list(superpose.line=list(lwd=2)), 
       main=paste("$PREFIX - $strain\n", nrow(data), "het sites - mappable")
       ))
       dev.off()

}
EOF
else
       echo "Skipping QC plots"
fi

echo "Done!!"
date
