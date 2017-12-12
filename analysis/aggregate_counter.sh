#!/bin/bash

#SBATCH -n 1
#SBATCH -t 0-04:00
#SBATCH -p serial_requeue
#SBATCH --mem-per-cpu=1500
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrewkane@g.harvard.edu

img_dir=$1
ntasks=$2
expt_type=${3-apc}
source new-modules.sh
source activate PYTO_SEG_ENV

cd $img_dir
python3 ~/code/pyto_seg_slurm/analysis/aggregate_counter.py -d $img_dir \
    -n $SLURM_ARRAY_TASK_ID -a $ntasks -t $expt_type
