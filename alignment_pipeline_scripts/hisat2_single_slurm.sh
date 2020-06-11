#!/bin/bash
#SBATCH --job-name=BLASTrun
#SBATCH --ntasks=8
#SBATCH --mem=32000


# Call hisat2 on the target file
export HISAT2_INDEXES=/home/abf/MouseEnsembl96
hisat2 -p8 --verbose --phred33 -x genome_tran -U "$FASTQDIR/$1" -S "$2_aligned_reads.sam"
samtools view -bS "$2_aligned_reads.sam" > "$2_aligned_reads.bam"
samtools sort "$2_aligned_reads.bam" "$2_sorted_alignment"
