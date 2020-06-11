#!/bin/bash
#SBATCH --job-name=LEC_TS_HTSEQ
#SBATCH --ntasks=1
#SBATCH --mem=16000


# Get Gene Level Counts for DBI Samples
export BAMDIR=/work/abf/LEC_Time_Series/DBI_Alignment_NoTrim
export OUTDIR=/work/abf/LEC_Time_Series/DBI_NoTrim_HTSeq_Count_Gene
for i in $(ls $BAMDIR | grep sorted_alignment.bam$)
do
    j=$(echo $i | sed 's/_sorted_alignment.bam//')
    j=$j"_GeneCount.txt"
    if [ ! -f $OUTDIR"/"$j ]
    then
	echo failed test for $OUTDIR"/"$j
	sbatch htseq_gene_slurm.sh "$BAMDIR"/"$i" "$j"
    fi
done

# Get Transcript Level Counts for DBI Samples
export OUTDIR=/work/abf/LEC_Time_Series/DBI_NoTrim_HTSeq_Count_Xscript
for i in $(ls $BAMDIR | grep sorted_alignment.bam$)
do
    j=$(echo $i | sed 's/_sorted_alignment.bam//')
    j=$j"_XscriptCount.txt"
    if [ ! -f $OUTDIR"/"$j ]
    then
	echo failed test
	sbatch htseq_transcript_slurm.sh "$BAMDIR"/"$i" "$j"
    fi
done


# Get Gene Level Counts for DNA_Link Samples
export BAMDIR=/work/abf/LEC_Time_Series/DNA_Link_Alignment_NoTrim
export OUTDIR=/work/abf/LEC_Time_Series/DNA_Link_NoTrim_HTSeq_Count_Gene
for i in $(ls $BAMDIR | grep sorted_alignment.bam$)
do
    j=$(echo $i | sed 's/_sorted_alignment.bam//')
    j=$j"_GeneCount.txt"
    if [ ! -f $OUTDIR"/"$j ]
    then
	echo failed test
	sbatch htseq_gene_slurm.sh "$BAMDIR"/"$i" "$j"
    fi
done

# Get Transcript Level Counts for DNA_Link Samples
export OUTDIR=/work/abf/LEC_Time_Series/DNA_Link_NoTrim_HTSeq_Count_Xscript
for i in $(ls $BAMDIR | grep sorted_alignment.bam$)
do
    j=$(echo $i | sed 's/_sorted_alignment.bam//')
    j=$j"_XscriptCount.txt"
    if [ ! -f $OUTDIR"/"$j ]
    then
	echo sbatch htseq_transcript_slurm.sh "$BAMDIR"/"$i" "$j"
	sbatch htseq_transcript_slurm.sh "$BAMDIR"/"$i" "$j"
    fi
done
