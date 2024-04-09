#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicslong
#SBATCH --job-name LCMB_pool1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 150G
#SBATCH --time 200:00:00
#SBATCH --output /projects/p31535/Anne/CellRanger/Logs/%x_oe%j.log
#SBATCH --verbose

date

cd /projects/b1042/Gate_Lab/AN1792-vacc/LCMB_CellRanger

/projects/p31535/thomas/SEA_AD/cellranger-7.2.0/cellranger multi \
--id="LCMB_CellRanger_pool1" \
--csv="/projects/p31535/Anne/CellRanger/Data/Gate20_pool1.csv"