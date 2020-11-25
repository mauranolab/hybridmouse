#!/bin/bash
set -e -o pipefail

strain=$1
celltype=$2
DS=$3
mappedgenome=$4

name=${celltype}-${DS}.${mappedgenome}
if [[ "$strain" != "C57B6" ]]; then
       name=${strain}-${name}
fi

echo "output name $name"
date


#Here we generate a list of expected filenames
files=""
#NB originally checked grep R1 but some old solexa lanes didn't include it
for f1 in `grep $DS inputs.txt | grep -v R2`; do
       #echo "adding $f1"
       
       #First get FC
       #NB This must match exactly what is in runAnalysis.sh
       #BUGBUG a bit fragile
       #BUGBUG misses things like UwStam_CH12-DS22536-FCD0TGK-L002_R1_001.fastq.gz and UwStam_CH12-DS22542-FCD0TDB-L001_R1_002.fastq.gz, but don't see any collision
       if [[ `basename $f1` =~ "^s_" ]]; then
              fc=`readlink -f $f1 | perl -pe 's/\/\d\d\d\/s_/\/s_/g;' | xargs dirname | xargs basename | perl -pe 's/_\d+_tag//g;'`
       else
              fc=`readlink -f $f1 | xargs dirname | xargs dirname | xargs dirname | xargs basename | perl -pe 's/_\d+_tag//g;'`
       fi
       if [[ "$fc" == "." ]] ; then
              fc=""
       else
#              echo "Flowcell $fc"
              fc="${fc}."
       fi
       
       
       f2=`echo $f1 | perl -pe 's/_R1(_\d+)?/_R2$1/g;'`
       if echo $f2 | grep -q R2 && grep -q $f2 inputs.txt ; then
              f1=`echo $f1 | perl -pe 's/_R1(_\d+)?/_R1R2\1/g;'`
       fi
       
       curOutputFile=`basename $f1 .fastq.gz`
       curOutputFile="${fc}${curOutputFile}.$mappedgenome.bam"
       
       #echo "cur $files"
       
       if [[ -f "$curOutputFile" ]]; then
              files="$files ${curOutputFile}"
              numfiles=$((numfiles+1))
       else
              echo  "Error: $curOutputFile doesn't exist!"
              #Don't die to work around weirdness in pipeline FQ files
              #exit 1
       fi
done


if [[ "$numfiles" -eq 0 ]]; then
       echo "Error: No files found to merge!"
       exit 2
fi

echo "Will merge $files"


if [[ "$numfiles" -eq 1 ]]; then
       echo "copying file"
       cp $files $name.bam
else
       echo "merging files"
       samtools merge -l 1 -@ $NSLOTS $name.bam $files
fi

echo "Processing $name.bam"
#TODO AddOrReplaceReadGroups or do bwa -r ""


echo
echo "mark dups"
date
samtools sort -@ $NSLOTS -O sam -T $TMPDIR/${name}.sortbyname -l 1 -n $name.bam |
samblaster |
samtools view -Sb - |
samtools sort -@ $NSLOTS -O bam -T $TMPDIR/${name}.sort -l 9 - > $name.markedDups.bam && mv $name.markedDups.bam $name.bam


echo
echo "Indexing"
date
samtools index $name.bam



rm -f $files

echo
echo "Done with merge"
date


echo
/vol/mauranolab/mauram01/hybridmice/dnase/src/makeTracks.sh $name $DS $mappedgenome


echo
echo -e "\nDone!"
date
