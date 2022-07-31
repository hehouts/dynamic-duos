#! /bin/bash -l
#
#SBATCH --mail-user=hehouts@ucdavis.edu         # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS - ALL, NONE, BEGIN, END, FAIL, REQUEUE
#SBATCH -J movefilestr                           # JOB ID
#SBATCH -e movefilestr.j%j.err                   # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o movefilestr.j%j.out                   # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH -c 1                                    # NUMBER OF PROCESSORS PER TASK
#SBATCH --ntasks=1                              # MINIMUM NUMBER OF NODES TO ALLOCATE TO JOB
#SBATCH --mem=2Gb                               # MEMORY POOL TO ALL CORES
#SBATCH --time=04-00:00:00                      # REQUESTED WALL TIME
#SBATCH -p bmh                                # PARTITION TO SUBMIT TO


mv ~/dynamic-duos-virome/snakecrate/data/raw_data/std/gzips/* /group/ctbrowngrp2/hehouts/ibd_gzips/std/



