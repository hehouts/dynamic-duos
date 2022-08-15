#! /bin/bash -l
#
#SBATCH -p bmh
#SBATCH -J grist
#SBATCH --time=04-00:00:00
#SBATCH --mem=12G
#SBATCH -c 1
#SBATCH -e jobs_vir/sra_mgx.j%j.err
#SBATCH -o jobs_vir/sra_mgx.j%j.out
#SBATCH --mail-user=hehouts@ucdavis.edu
#SBATCH --mail-type=ALL


. "/home/hehouts/miniconda3/etc/profile.d/conda.sh"
conda activate grist-dev


#genome-grist run config_sra5smpl.yml summarize_mapping -n --rerun-incomplete

genome-grist run config_sra_vir.yml summarize_mapping --rerun-incomplete  -j1 --latency-wait 180 --keep-going



#genome-grist run config_sra5smpl.yml summarize_mapping -n --rerun-incomplete --unlock
#echo "!!!!!!!! Print Me !!!!!!" 
#--rerun-incomplete  -n --latency-wait 180 --keep-going 

