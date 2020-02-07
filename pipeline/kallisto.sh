#!/usr/bin/bash
#SBATCH -p short --mem 24gb -N 1 -n 16  --out logs/kallisto.%a.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

module load kallisto
INDIR=fastq
OUTDIR=results
IDX=Cimm_RS_v2.idx
TX=Coccidioides_immitis_RS.mrna-transcripts.fa
SAMPLEFILE=samples.csv
mkdir -p $OUTDIR
if [ ! -f $IDX ]; then
    gzip -dc data/$TX.gz > /tmp/$TX
    kallisto index -i $IDX /tmp/$TX
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
 N=$1
 if [ -z $N ]; then
     echo "cannot run without a number provided either cmdline or --array in sbatch"
     exit
 fi
fi

IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read SAMPLE CONDITION REP READBASE
do
 OUTNAME=$CONDITION.r${REP}
 if [ ! -f $OUTDIR/$OUTNAME/abundance.h5 ]; then
     kallisto quant -i $IDX -o $OUTDIR/$OUTNAME -t $CPU --bias $INDIR/${READBASE}_R1_001.fastq.gz $INDIR/${READBASE}_R2_001.fastq.gz
 fi
done

