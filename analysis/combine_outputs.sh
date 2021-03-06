#!/bin/bash

#SBATCH -n 1
#SBATCH -t 0-01:00
#SBATCH -p serial_requeue
#SBATCH --mem-per-cpu=4000
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nweir@fas.harvard.edu

csv_dir=$1
source new-modules.sh
source activate PYTO_SEG_ENV

cd $csv_dir
python3 ~/code/batch_segmentation/combine_outputs.py -d $csv_dir
