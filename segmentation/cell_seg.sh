#!/bin/bash

#SBATCH -n 1
#SBATCH -t 0-12:00
#SBATCH -p serial_requeue
#SBATCH --mem-per-cpu=1500
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrewkane@g.harvard.edu

img_dir=$1
threshold=$2
ntasks=$3
source new-modules.sh
source activate PYTO_SEG_ENV

cd $img_dir
python3 ~/code/pyto_seg_slurm/segmentation/cell_seg.py -d $img_dir -t $threshold \
	 -n $SLURM_ARRAY_TASK_ID -a $ntasks
