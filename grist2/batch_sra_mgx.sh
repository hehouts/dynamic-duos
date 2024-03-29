#! /bin/bash -l
#
#SBATCH -p bmh
#SBATCH -J grist
#SBATCH --time=04-00:00:00
#SBATCH --mem=12G
#SBATCH -c 1
#SBATCH -e jobs/sra_mgx.j%j.err
#SBATCH -o jobs/sra_mgx.j%j.out
#SBATCH --mail-user=hehouts@ucdavis.edu
#SBATCH --mail-type=ALL


. "/home/hehouts/miniconda3/etc/profile.d/conda.sh"
conda activate grist-dev


#genome-grist run config_sra5smpl.yml summarize_mapping -n --rerun-incomplete

genome-grist run config_sra_mgx.yml summarize_mapping --rerun-incomplete  -j1 --latency-wait 180 --keep-going



#genome-grist run config_sra5smpl.yml summarize_mapping -n --rerun-incomplete --unlock
#echo "!!!!!!!! Print Me !!!!!!" 
#--rerun-incomplete  -n --latency-wait 180 --keep-going 

