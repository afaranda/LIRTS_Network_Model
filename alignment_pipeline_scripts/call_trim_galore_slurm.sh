#!/bin/bash
#SBATCH --job-name=BLASTrun
#SBATCH --ntasks=8
#SBATCH --mem=16000


for i in $(ls /home/abf/Fibronectin_Analysis | grep fastq.gz)
do
	#p=$(pwd)
	j=$(echo $i | sed 's/_L001_R1_001\.fastq\.gz//')
	j=$(echo $j | sed 's/^[0-9]\+_//')
	#mkdir "HISAT2_$j"
	#cd "HISAT2_$j"
	echo $j
	sbatch trim_galore_slurm.sh $i
	#cd $p
done
