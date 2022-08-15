#! /bin/bash -l
#
#SBATCH -p bmm
#SBATCH -J grist
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH -c 1
#SBATCH -e jobs/unlock.j%j.err
#SBATCH -o jobs/unlock.j%j.out


. "/home/hehouts/miniconda3/etc/profile.d/conda.sh"
conda activate grist-dev

genome-grist run config_sra_vir.yml -n --rerun-incomplete --unlock 
genome-grist run config_sra_mgx.yml -n --rerun-incomplete --unlock 


