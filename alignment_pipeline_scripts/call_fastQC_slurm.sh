#!/bin/bash
#SBATCH --job-name=BLASTrun
#SBATCH --ntasks=8
#SBATCH --mem=16000



# Align 0, 24 and 48 hour data generated by DBI
export FASTQDIR=/home/abf/Fibronectin_Analysis
export OUTDIR=/work/abf/LEC_Time_Series/DBI_FastQC
echo "Processing Files: in $FASTQDIR"
for i in $(ls $FASTQDIR | grep R1_001.fastq.gz)
do
    echo $i

    j=$(echo $i | sed 's/_L001_R1_001\.fastq\.gz//')
    #j=$(echo $j | sed 's/^[0-9]\+_//')
    echo $FASTQDIR"/"$i
    sbatch fastqc_slurm.sh $FASTQDIR"/"$i "$OUTDIR"
done

# Align 0, 6, and 24 hour data generated by DNA_Link
export FASTQDIR=/home/abf/Melinda_Duncan_DNALink/Duncan20180919
export OUTDIR=/work/abf/LEC_Time_Series/DNA_Link_FastQC
echo "Processing Files in: $FASTQDIR"
for i in $(ls $FASTQDIR | grep fastq.gz)
do
    echo $FASTQDIR"/"$i
    sbatch fastqc_slurm.sh $FASTQDIR"/"$i "$OUTDIR"
done