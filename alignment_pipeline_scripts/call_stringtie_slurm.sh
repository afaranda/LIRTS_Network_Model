#!/bin/bash
#SBATCH --job-name=LEC_TS_STRTI
#SBATCH --ntasks=1
#SBATCH --mem=16000


# Run Stringtie transcriptome assembly for DBI Samples
export BAMDIR=/work/abf/LEC_Time_Series/DBI_Alignment_NoTrim
export OUTDIR=/work/abf/LEC_Time_Series/DBI_NoTrim_StringTie
for i in $(ls $BAMDIR | grep sorted_alignment.bam$)
do
    j=$(echo $i | sed 's/_sorted_alignment.bam$//')
    if [ ! -d $OUTDIR"/"$j ]
    then
	sbatch stringtie_slurm.sh "$BAMDIR"/"$i"
    fi
done

# Run Stringtie transcriptome assembly for DNA Link Samples
export BAMDIR=/work/abf/LEC_Time_Series/DNA_Link_Alignment_NoTrim
export OUTDIR=/work/abf/LEC_Time_Series/DNA_Link_NoTrim_StringTie
for i in $(ls $BAMDIR | grep sorted_alignment.bam$)
do
    j=$(echo $i | sed 's/_sorted_alignment.bam//')
    if [ ! -d $OUTDIR"/"$j ]
    then
	sbatch stringtie_slurm.sh "$BAMDIR"/"$i"
    fi
done
