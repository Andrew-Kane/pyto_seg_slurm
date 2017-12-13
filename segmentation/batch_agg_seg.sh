#!/bin/bash

#SBATCH -n 1
#SBATCH -t 0-02:00
#SBATCH -p serial_requeue
#SBATCH --mem-per-cpu=1500
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrewkane@g.harvard.edu

img_dir=$1
high_threshold=$2
low_threshold=$3
ntasks=$4
source new-modules.sh
source activate PYTO_SEG_ENV

cd $img_dir
imgs=($(find . -maxdepth 1 -iregex '.*_w[0-9]488.*\.TIF'))
echo ${imgs}
nimgs_per_task=`expr ${#imgs[@]} / $ntasks`
task_chk=$(($ntasks - 1))
first_img=$((${SLURM_ARRAY_TASK_ID} * $nimgs_per_task + 1))
if [ "${SLURM_ARRAY_TASK_ID}" -eq  "$task_chk" ]
    then python3 ~/code/pyto_seg_slurm/segmentation/batch_agg_seg.py -d $img_dir -ht $high_threshold -lt \
        $low_threshold  "${imgs[@]:$first_img}"
    else python3 ~/code/pyto_seg_slurm/segmentation/batch_agg_seg.py -d $img_dir -ht $high_threshold -lt \
        $low_threshold "${imgs[@]:$first_img:$nimgs_per_task}"
fi
