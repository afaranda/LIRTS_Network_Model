#!/bin/bash
#SBATCH --job-name=BLASTrun
#SBATCH --ntasks=4
#SBATCH --mem=32000


# Call hisat2 on the target file
echo "htseq-count -i transcript_id -f bam -s reverse -m union --type exon $1 /home/abf/MouseEnsembl96/Mus_musculus.GRCm38.96.chr.gtf > $2"
htseq-count -i transcript_id -r pos -f bam -s reverse -m union --type exon $1\
	    /home/abf/MouseEnsembl96/Mus_musculus.GRCm38.96.chr.gtf >\
	    $OUTDIR"/"$2


