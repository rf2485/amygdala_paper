#!/bin/bash
#SBATCH --mail-user=rf2485@nyulangone.org
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --array=1-306
#SBATCH -o ./slurm_output/2_image_pipeline/slurm-%A_%a.out

basedir=/gpfs/data/lazarlab/CamCan995/
projectdir=$basedir/derivatives/mti_whole_wm_gm/

subj_list=$(cut -f1 -d$'\t' $projectdir/mti_over_55.tsv)
subj_list=($subj_list)
subj_num=$(($SLURM_ARRAY_TASK_ID))
j=${subj_list[$subj_num]}

rawanat=$basedir/raw/$j/anat/${j}
rawmti=$basedir/raw/$j/mti/${j}

freesurferdir=$projectdir/freesurfer
mkdir -p $freesurferdir

module load singularity/3.9.8
module load freesurfer/7.4.1
export SUBJECTS_DIR=$freesurferdir
module load fsl/6.0.7

if [ ! -f $freesurferdir/$j/mri/wmparc.mgz ]; then
	#recon-all-clinical to generate segmentations
	rm $freesurferdir/$j/scripts/IsRunning.lh+rh
	recon-all-clinical.sh ${rawanat}_T1w.nii.gz $j 4 $freesurferdir
fi

if [ ! -f $freesurferdir/$j/label/BA_exvivo.thresh.ctab ]; then
	#create files necessary for bbregister
	mkdir -p $freesurferdir/$j/mri/orig
	recon-all -motioncor -talairach -nuintensitycor -normalization -gcareg -pctsurfcon -segstats -wmparc -balabels -subjid ${j} -sd ${freesurferdir} -threads 4 -no-isrunning
fi

if [ ! -f $freesurferdir/$j/mti/mtr.nii.gz ]; then
	#copy raw mti images
	mkdir -p $freesurferdir/$j/mti/intermediate/
	cd $freesurferdir/$j/mti
	cp ${rawmti}_mti_bl.nii.gz mti_bl_raw.nii.gz
	cp  ${rawmti}_mti.nii.gz mti_raw.nii.gz
	#register mti_bl and mti
	flirt -in mti_raw -ref mti_bl_raw -dof 6 -out intermediate/1_mti_reg
	#degibbs mti_bl and registered mti
	singularity exec --bind $basedir/ $basedir/neurodock_latest.sif mrdegibbs mti_bl_raw.nii.gz intermediate/2_mti_bl_degibbs.nii.gz
	singularity exec --bind $basedir/ $basedir/neurodock_latest.sif mrdegibbs intermediate/1_mti_reg.nii.gz intermediate/2_mti_degibbs.nii.gz
	#calculate mtr
	fslmaths intermediate/2_mti_bl_degibbs.nii.gz -sub intermediate/2_mti_degibbs.nii.gz -div intermediate/2_mti_bl_degibbs.nii.gz -thr 0 -uthr 1 mtr.nii.gz
fi

if [ ! -f $freesurferdir/$j/mti/fs2mtr_field.mgz ]; then
	#register mtr to t1, using settings for T1 images (gray matter darker than white matter)
	cd $freesurferdir/$j
	mri_easyreg --ref mri/T1.mgz --flo mti/mtr.nii.gz --ref_seg mri/synthseg.mgz --flo_seg mti/mtr_synthseg.mgz --bak_field mti/fs2mtr_field.mgz --threads 4 --affine_only
	#register aparc+aseg to mtr, using the inverse of the mtr to t1 registration
	mri_easywarp --i mri/aparc+aseg.mgz --o mti/aparc+aseg2mtr.mgz --field mti/fs2mtr_field.mgz --threads 4 --nearest
fi

if [ ! -f $freesurferdir/$j/stats/mtr_subcort_gm.stats ]; then
	cd $freesurferdir/$j
	#calculate mean MTR in each region of aparc+aseg atlas
	mri_segstats --seg mti/aparc+aseg2mtr.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i mti/mtr.nii.gz --sum stats/mtr_aparc+aseg.stats
	#segment whole brain wm and gm in MTR
	mri_binarize --i mti/aparc+aseg2mtr.mgz --all-wm --o mti/mtr_wm.mgz #all wm
	mri_binarize --i mti/aparc+aseg2mtr.mgz --gm --o mti/mtr_gm.mgz #all gm
	mri_binarize --i mti/aparc+aseg2mtr.mgz --ctx-wm --o mti/mtr_ctx_wm.mgz #cortex wm
	mri_binarize --i mti/aparc+aseg2mtr.mgz --min 999.5 --o mti/mtr_ctx_gm.mgz #cortex gm
	mri_binarize --i mti/aparc+aseg2mtr.mgz --min 999.5 --max 1999.5 --o mti/mtr_lh_ctx_gm.mgz #left cortical gm
	mri_binarize --i mti/aparc+aseg2mtr.mgz --min 1999.5 --o mti/mtr_rh_ctx_gm.mgz #right cortical gm
	mri_binarize --i mti/aparc+aseg2mtr.mgz --subcort-gm --o mti/mtr_subcort_gm.mgz #subcortical gm	
	#right and left cortical WM already included in aparc+aseg

	#calculate mean mtr for each segmentation
	mri_segstats --seg mti/mtr_wm.mgz --mask mti/mtr_wm.mgz --excludeid 0 --i mti/mtr.nii.gz --sum stats/mtr_wm.stats
	mri_segstats --seg mti/mtr_gm.mgz --mask mti/mtr_gm.mgz --excludeid 0 --i mti/mtr.nii.gz --sum stats/mtr_gm.stats
	mri_segstats --seg mti/mtr_ctx_wm.mgz --mask mti/mtr_ctx_wm.mgz --excludeid 0 --i mti/mtr.nii.gz --sum stats/mtr_ctx_wm.stats
	mri_segstats --seg mti/mtr_ctx_gm.mgz --mask mti/mtr_ctx_gm.mgz --excludeid 0 --i mti/mtr.nii.gz --sum stats/mtr_ctx_gm.stats
	mri_segstats --seg mti/mtr_lh_ctx_gm.mgz --mask mti/mtr_lh_ctx_gm.mgz --excludeid 0 --i mti/mtr.nii.gz --sum stats/mtr_lh_ctx_gm.stats
	mri_segstats --seg mti/mtr_rh_ctx_gm.mgz --mask mti/mtr_rh_ctx_gm.mgz --excludeid 0 --i mti/mtr.nii.gz --sum stats/mtr_rh_ctx_gm.stats
	mri_segstats --seg mti/mtr_subcort_gm.mgz --mask mti/mtr_subcort_gm.mgz --excludeid 0 --i mti/mtr.nii.gz --sum stats/mtr_subcort_gm.stats
fi

#AD signature
if [ ! -f $freesurferdir/$j/stats/mtr_AD_sig.stats ]; then
	cd $freesurferdir/$j
	mri_easywarp --i mri/aparc.a2009s+aseg.mgz --o mti/aparc.a2009s+aseg2mtr.mgz --field mti/fs2mtr_field.mgz --threads 4 --nearest
	mri_binarize --i mti/aparc.a2009s+aseg2mtr.mgz --match-ctab $projectdir/AD_signature_LUT.txt --o mti/mtr_AD_sig.mgz
	mri_segstats --seg mti/mtr_AD_sig.mgz --mask mti/mtr_AD_sig.mgz --excludeid 0 --i mti/mtr.nii.gz --sum stats/mtr_AD_sig.stats
fi


