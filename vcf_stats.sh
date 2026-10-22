#!/bin/bash
#SBATCH --job-name=vcf_stats
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=30:00:00
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=NETID@uvm.edu

module load gcc
module load vcftools

##set your path to where the data is, change to your own path!!
INDIR="/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics/mydata"
##set your path to where the data will be stored, change to your own path!!
OUTDIR="/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics/myresults"

##change to directory where data is
cd $INDIR
##change out directory to where the results will go
vcftools --gzvcf Lottia_gigantea_pl_flt3.renamed.vcf.gz --missing-indv --out $OUTDIR/missing_data
vcftools --gzvcf Lottia_gigantea_pl_flt3.renamed.vcf.gz --depth --out $OUTDIR/depth_data

#Sumbit and keep eye on it
sbatch vcf_stats.sh      # submit the job; Slurm replies with a job ID
squeue -u $USER          # check your jobs (PD = pending, R = running)
scancel <jobid>          # cancel a job if you need to
