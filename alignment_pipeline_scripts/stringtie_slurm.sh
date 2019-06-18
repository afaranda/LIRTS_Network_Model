#!/bin/bash
#SBATCH --job-name=StringTie
#SBATCH --ntasks=8
#SBATCH --mem=512000


# Call stringtie on the target file
export GENOME_GTF=/home/abf/MouseEnsembl96/Mus_musculus.GRCm38.96.chr.gtf
S=$(echo $1 | sed 's/_sorted_alignment\.bam$//')
S=$(echo $S | awk '{ n=split($NF, a, "/"); print a[n]}')

echo $OUTDIR
echo $OUTDIR/$S
stringtie $1 -p4 --fr -e -B\
	  -G $GENOME_GTF\
	  -o $OUTDIR"/"$S"/"$S.gtf\
	  -A $OUTDIR"/"$S"/"$S.txt

 
