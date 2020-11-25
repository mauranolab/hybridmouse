#!/bin/bash
set -e -o pipefail

#mergeSamples.sh: for merging fully mapped .bam files (e.g. for multiple sub-libraries)

jobid=$SGE_TASK_ID
#jobid=1


echo "Running on $HOSTNAME. Using $TMPDIR as tmp"
date


name=`cat inputs.txt | head -n $jobid | tail -1`
strain=`echo $name | cut -d '.' -f1 | cut -d '-' -f1`
celltype=`echo $name | cut -d '.' -f1 | cut -d '-' -f2`
DS=`echo $name | cut -d '.' -f1 | cut -d '-' -f3`
mappedgenome=`echo $name | cut -d '.' -f2`

echo "output name $name"


#Here we generate a list of expected filenames
files=""
for curOutputFile in $(ls /vol/mauranolab/mauram01/hybridmice/dnase/${strain}-${celltype}-DS*.${mappedgenome}.bam); do
       #echo "adding $f1"
       if [[ "$DS" == "DSall" || "$curOutputFile" =~ "${DS}" ]]; then
              if [[ -f "$curOutputFile" ]]; then
                     files="$files ${curOutputFile}"
                     numfiles=$((numfiles+1))
              else
                     echo  "Error: $curOutputFile doesn't exist!"
                     #Don't die to work around weirdness in pipeline FQ files
                     #exit 1
              fi
       fi
done

echo "Will merge $files"


if [[ "$numfiles" -eq 1 ]]; then
       echo "copying file"
       cp $files $TMPDIR/$name.bam
else
       echo "merging files"
       samtools merge -l 1 -@ $NSLOTS $TMPDIR/$name.bam $files
fi

#was -q 30 -F 1548
samflags="-F 524"
samtools view -b -1 -@ $NSLOTS $samflags $TMPDIR/$name.bam > $name.bam
samtools index $name.bam


echo
echo "Done with merge"
date


echo
/vol/mauranolab/mauram01/hybridmice/dnase/src/makeTracks.sh $name $DS $mappedgenome


echo
echo -e "\nDone!"
date
