#!/bin/bash
set -e -o pipefail

#module load samtools/1.2
#module load gcc/4.7.2     # for adapter trimming
#module load R/3.2.4       # for RSEM
#module load coreutils/8.9 # parallel sort

module load star/2.5.2a
module load rsem/1.2.31

PRIORITY=0

SLOTS=4


R1_FASTQ=$1
R2_FASTQ=$2
#eg "/vol/isg/annotation/STARindex/mm10all/"
STARdir=$3
#eg "/vol/isg/annotation/RSEMindex/mm10all/"
RSEMdir=$4
SAMPLE_NAME=$5

ALIGNMENT_ID=`date +%y%m%d`
echo -e "Mapping:\n$R1_FASTQ\n$R2_FASTQ\n"

outdir=$(pwd)

STAMPIPES="/home/mauram01/scratch/hybridmice/rnaseq/src"
scriptdir="$STAMPIPES/scripts"
script="$scriptdir/STAR_RSEM.sh"


TRIMDIR="trimmed"
TRIM_R1=$TRIMDIR/$(basename "$R1_FASTQ")
TRIM_R2=$TRIMDIR/$(basename "$R2_FASTQ")
mkdir -p "$TRIMDIR"

dataType="str_PE"  # 4 types; str_SE str_PE unstr_SE unstr_PE

ADAPTER_FILE=${SAMPLE_NAME}.adapters.txt
VERSION_FILE=${SAMPLE_NAME}.versions.txt

ADAPTER_P7=GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGCGATAGATCTCGTATGCCGTCTTCTGCTTG
ADAPTER_P5=GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG


bash "$scriptdir/versions.sh" &> "$VERSION_FILE"
if [[ ( -n "$ADAPTER_P7" ) && ( -n "$ADAPTER_P5" ) ]] ; then
  echo -e "P7\t$ADAPTER_P7\nP5\t$ADAPTER_P5" > "$ADAPTER_FILE"
fi

# Perform trimming
if [[ ( "$ADAPTER_P7"  == "NOTAVAILABLE" ) || ( "$ADAPTER_P5" == "NOTAVAILABLE" ) ]] ; then
  TRIM_R1=$R1_FASTQ
  TRIM_R2=$R2_FASTQ
else 
  if [[ ( ! -e "$TRIM_R1" ) || ( ! -e "$TRIM_R2" ) ]] ; then
    echo "Trimming..."
    /home/mauram01/src/jvierstra-bio-tools-6fe54fa5a3d9/apps/trim-adapters-illumina/trim-adapters-illumina -f "$ADAPTER_FILE" \
      --threads=2 \
      -1 P5 -2 P7 \
      "$R1_FASTQ" \
      "$R2_FASTQ" \
      "$TRIM_R1.tmp" \
      "$TRIM_R2.tmp" \
      &> "$outdir/adapter_trimming.txt"

    mv "$TRIM_R1.tmp" "$TRIM_R1"
    mv "$TRIM_R2.tmp" "$TRIM_R2"
  fi
fi

jobbase="${SAMPLE_NAME}-ALIGN#${ALIGNMENT_ID}"
starjob="rs.$jobbase"
uploadjob="up.$jobbase"

# Run STAR & RSEM
# TODO: Break these into sub-steps
set +e
$scriptdir/checkcomplete.sh $SAMPLE_NAME
iscomplete=$?
set -e

if [[ $iscomplete != 0 ]]; then
  echo "Run STAR & RSEM"
  qsub -j y -cwd -V -N "$starjob" -pe threads "$SLOTS" -p "$PRIORITY" -S /bin/bash << __RNA-STAR__
    set -x

#Cache for performance
#    STARdir=\$("$STAMPIPES/scripts/cache.sh" "$STARdir")
#    RSEMdir=\$("$STAMPIPES/scripts/cache.sh" "$RSEMdir")
    STARdir=$STARdir
    RSEMdir=$RSEMdir

    cd "$outdir"

    # Run!
    "$script" "$TRIM_R1" "$TRIM_R2" "\$STARdir" "\$RSEMdir/RSEMref" "$dataType"
__RNA-STAR__
else
  echo "Already complete, no further processing!"
fi

# Check for completeness and upload files.
qsub -j y -cwd -V -hold_jid "$starjob" -N "$uploadjob" -S /bin/bash << __UPLOAD__
  set -e
  $scriptdir/checkcomplete.sh $SAMPLE_NAME
__UPLOAD__
