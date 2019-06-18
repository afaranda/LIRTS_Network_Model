#!/bin/bash
#SBATCH --job-name=BLASTrun
#SBATCH --ntasks=8
#SBATCH --mem=512000


# Call hisat2 on the target file
export HISAT2_INDEXES=/home/abf/MouseEnsembl96

hisat2 -p8 --verbose --phred33 -x genome_tran -1 "$FASTQDIR/$1" -2 "$FASTQDIR/$2" -S "$OUTDIR/$3_aligned_reads.sam"
samtools view -bS "$OUTDIR/$3_aligned_reads.sam" > "$OUTDIR/$3_aligned_reads.bam"
samtools sort "$OUTDIR/$3_aligned_reads.bam" "$OUTDIR/$3_sorted_alignment"
samtools index "$OUTDIR/$3_sorted_alignment.bam"
rm "$OUTDIR/$3_aligned_reads.bam"
