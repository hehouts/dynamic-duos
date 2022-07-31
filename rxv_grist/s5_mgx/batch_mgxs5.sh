#! /bin/bash -l
#
#SBATCH -p bmh
#SBATCH -J grist
#SBATCH --time=00-12:00:00
#SBATCH --mem=8G
#SBATCH -c 1
#SBATCH -e outputs.s5_mgx/jobs/s5mgx.j%j.err
#SBATCH -o outputs.s5_mgx/jobs/s5mgx.j%j.out
#SBATCH --mail-user=hehouts@ucdavis.edu
#SBATCH --mail-type=ALL


. "/home/hehouts/miniconda3/etc/profile.d/conda.sh"
conda activate grist

#genome-grist run config_s5_mgx.yml summarize_mapping -n --rerun-incomplete
genome-grist run config_s5_mgx.yml summarize_mapping --rerun-incomplete  -j16 --keep-going --latency-wait 180 

#echo "!!!!!!!! Print Me !!!!!!"
#--rerun-incomplete  -n --latency-wait 180 --keep-going 

