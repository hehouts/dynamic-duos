#!/bin/bash -l
#SBATCH --job-name=snakemake_vir
#SBATCH --partition=bmm
#SBATCH --output=bashjobs/bash_report_%j.output
#SBATCH --error=bashjobs/bash_err_%j.output
#SBATCH --time=1:00:00

#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 


export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks

srun conda activate dynduo
srun snakemake -n --rerun-incomplete