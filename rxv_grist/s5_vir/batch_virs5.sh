#! /bin/bash -l
#
#SBATCH -p bmh
#SBATCH -J vir_grist
#SBATCH --time=00-04:00:00
#
#SBATCH --mem=28G
#SBATCH -c 1
#
#SBATCH -e outputs.s5_vir/jobs/s5vir.j%j.err
#SBATCH -o outputs.s5_vir/jobs/s5vir.j%j.out
#SBATCH --mail-user=hehouts@ucdavis.edu
#SBATCH --mail-type=ALL


. "/home/hehouts/miniconda3/etc/profile.d/conda.sh"
conda activate grist

#genome-grist run config_s5_vir.yml summarize_mapping -n --rerun-incomplete
genome-grist run config_s5_vir.yml summarize_mapping --rerun-incomplete  -j16 --keep-going --latency-wait 180 


