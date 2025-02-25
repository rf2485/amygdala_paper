basedir=/gpfs/data/lazarlab/CamCan995/
# basedir=/Volumes/Research/lazarm03lab/labspace/AD/camcan995/
projectdir=$basedir/derivatives/scd/main/
freesurferdir=$projectdir/freesurfer

export FREESURFER_HOME=/gpfs/share/apps/freesurfer/7.4.1/raw/freesurfer/
module load freesurfer/7.4.1
export SUBJECTS_DIR=$freesurferdir
module load miniconda3/gpu/4.9.2
source activate ~/.conda/envs/fsl_eddy/
export FSLDIR=$CONDA_PREFIX
source $FSLDIR/etc/fslconf/fsl.sh

subj_list=$(cut -f1 -d$'\t' $projectdir/dwi_over_55.tsv)
subj_list=($subj_list)

for j in "${subj_list[@]}"; do
	echo $j
	designerdir=$projectdir/dwi_processed/$j/
	mri_convert $designerdir/residual.nii $designerdir/residual.nii.gz
	rm $designerdir/residual.nii
	mri_convert $designerdir/metrics/dti_V1.nii $designerdir/metrics/dti_V1.nii.gz
	rm $designerdir/metrics/dti_V1.nii
	smi_list=( smi_matlab_f smi_matlab_Da smi_matlab_DePar smi_matlab_DePerp smi_matlab_p2 )
	for meas in "${smi_list[@]}"; do
		mri_convert $designerdir/metrics/$meas.nii $designerdir/metrics/$meas.nii.gz
		rm $designerdir/metrics/$meas.nii
	done
	mri_convert $designerdir/ODI_mask.nii $designerdir/ODI_mask.nii.gz
	rm $designerdir/ODI_mask.nii
	mri_convert $designerdir/FWF_mask.nii $designerdir/FWF_mask.nii.gz
	rm $designerdir/FWF_mask.nii
	mri_convert $freesurferdir/$j/diffusion/B0.nii $freesurferdir/$j/diffusion/B0.nii.gz
	rm $freesurferdir/$j/diffusion/B0.nii
	nii_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd smi_matlab_f smi_matlab_Da smi_matlab_DePar smi_matlab_DePerp smi_matlab_p2 mtr2diff g_ratio )
	for meas in "${nii_list[@]}"; do
		mri_convert $freesurferdir/$j/diffusion/${meas}.nii $freesurferdir/$j/diffusion/${meas}.nii.gz
		rm $freesurferdir/$j/diffusion/${meas}.nii
	done
done
