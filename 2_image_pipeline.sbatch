#!/bin/bash
#SBATCH --mail-user=rf2485@nyulangone.org
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --array=1-128
#SBATCH -o ./slurm_output/2_image_pipeline/slurm-%A_%a.out

basedir=/gpfs/data/lazarlab/CamCan995/
projectdir=$basedir/derivatives/mti_whole_wm_gm/

subj_list=$(cut -f1 -d$'\t' $projectdir/mti_over_55_tr50.tsv)
subj_list=($subj_list)
subj_num=$(($SLURM_ARRAY_TASK_ID))
j=${subj_list[$subj_num]}

rawanat=$basedir/raw/$j/anat/${j}
rawmti=$basedir/raw/$j/mti/${j}

freesurferdir=$projectdir/freesurfer
mkdir -p $freesurferdir

module load freesurfer/7.4.1
export SUBJECTS_DIR=$freesurferdir
module load fsl/6.0.7

#calculate magnetization transfer ratio
##register bl and mti
mkdir -p $freesurferdir/$j/mti/tmp
mkdir -p $freesurferdir/$j/stats
cd $freesurferdir/$j/mti/tmp
cp ${rawmti}_mti_bl.nii.gz mti_bl_raw.nii.gz
cp  ${rawmti}_mti.nii.gz mti_raw.nii.gz
flirt -in mti_raw -ref mti_bl_raw -cost mutualinfo -interp sinc -searchrx -30 30 -searchry -30 30 -searchrz -30 30 -dof 6 -out mti_reg
fslmaths mti_bl_raw.nii.gz -sub mti_reg.nii.gz -div mti_bl_raw.nii.gz -thr 0 -uthr 1 $freesurferdir/$j/mti/mtr.nii.gz
cd $freesurferdir/$j/mti
mri_synthstrip -i mtr.nii.gz -o mtr_brain.nii.gz --no-csf --threads 4
mri_synthseg --i mtr_brain.nii.gz --o tmp/mtr_brain_seg_upsampled.nii.gz --threads 4
mri_convert tmp/mtr_brain_seg_upsampled.nii.gz mtr_brain_seg.mgz --rl mtr_brain.nii.gz --rt nearest
cd $freesurferdir/$j
mri_segstats --seg mti/mtr_brain_seg.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i mti/mtr_brain.nii.gz --sum stats/mtr_brain.stats
rm -rf $freesurferdir/$j/mti/tmp