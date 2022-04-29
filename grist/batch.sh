#! /bin/bash
#
#SBATCH -p bmm
#SBATCH -J grist
#SBATCH --time=08-00:00:00
#SBATCH --mem=180G
#SBATCH -c 1
#SBATCH -e ~/dynamic-duos-virome/grist/outputs.s5_mgx/jobs/s5mgx.j%j.err
#SBATCH -o ~/dynamic-duos-virome/grist/outputs.s5_mgx/jobs/s5mgx.j%j.out
#SBATCH --mail-user=hehouts@ucdavis.edu
#SBATCH --mail-type=ALL

echo Hello World
sleep 15
date



# initialize conda
#. ~/mambaforge/etc/profile.d/conda.sh

#genome-grist run config_s5_mgx.yml summarize_mapping --rerun-incomplete  -j16 --latency-wait 180 --keep-going 


