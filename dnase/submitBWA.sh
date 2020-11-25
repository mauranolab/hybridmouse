#!/bin/bash
set -e

strain=$1
celltype=$2
DS=$3
name=${strain}-${celltype}-${DS}

base=/vol/mauranolab/mauram01/hybridmice/dnase

mkdir -p tmp

if ! grep -q $DS inputs.txt; then
       echo "Can't find $DS"
       exit 1
fi

numlines=`grep $DS inputs.txt | wc -l`

firstline=`grep -n $DS inputs.txt | head -1 | awk -F ":" '{print $1}'`
lastline=`grep -n $DS inputs.txt | tail -1 | awk -F ":" '{print $1}'`
expectednumlines=`echo "$lastline-$firstline+1" | bc -l -q`
if [ "$expectednumlines" -ne "$numlines" ]; then
       echo "Wrong number of lines (expected $expectednumlines, found $numlines), are all fastq files contiguous in inputs.txt?"
       exit 2
fi


#SGE doesn't accept a complicated -t array, so we'll start R2 jobs that will die instantly rather than prune here
echo "Processing $name (input.txt lines $firstline-$lastline) for strain $strain"
qsub -S /bin/bash -cwd -V -pe threads 6 -terse -j y -b y -t $firstline-$lastline -N map.$name "$base/src/runAnalysis.sh $celltype $DS $strain" | perl -pe 's/[^\d].+$//g;' > sgeid

echo -n "Your job "
cat sgeid | perl -pe 's/\n/ /g;'
echo "has been submitted"


if [[ "$strain" != "C57B6" ]]; then
       echo "C57B6 merge"
       qsub -S /bin/bash -cwd -V -terse -j y -b y -hold_jid `cat sgeid` -N merge.$name.C57B6 "$base/src/BWAmerge.sh $strain $celltype $DS C57B6" | perl -pe 's/[^\d].+$//g;' > sgeid.C57B6
fi


echo "$strain merge"
qsub -S /bin/bash -cwd -V -terse -j y -b y -hold_jid `cat sgeid` -N merge.$name.$strain "$base/src/BWAmerge.sh $strain $celltype $DS $strain" | perl -pe 's/[^\d].+$//g;' > sgeid.$strain

if [[ "$strain" != "C57B6" ]]; then
       mkdir -p counts
       cd counts
       #NB requires merged hotspots
       qsub -cwd -V -hold_jid `cat ../sgeid.C57B6`,`cat ../sgeid.$strain` -N counts.$name -S /bin/bash -b y -j y "$base/src/makeCounts.sh ../$name.$strain.bam $strain"
       cd ..
fi

rm -f sgeid sgeid.C57B6 sgeid.$strain

echo
