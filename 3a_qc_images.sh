basedir=/gpfs/data/lazarlab/CamCan995/
projectdir=$basedir/derivatives/scd/gm_roi/
t1dir=$projectdir/freesurfer/
dwidir=$projectdir/dwi_processed/

export FREESURFER_HOME=/gpfs/share/apps/freesurfer/7.4.1/raw/freesurfer/
module load freesurfer/7.4.1
export SUBJECTS_DIR=$t1dir
module load miniconda3/gpu/4.9.2
source activate ~/.conda/envs/fsl_eddy/
export FSLDIR=$CONDA_PREFIX
source $FSLDIR/etc/fslconf/fsl.sh

subj_list=$(cut -f1 $projectdir/dwi_over_55.tsv)
subj_list=($subj_list)

mkdir -p $t1dir/group_qc/1_recon
echo '<HTML><TITLE>recon</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/1_recon/index.html

cd $t1dir/group_qc/1_recon
for subj in "${subj_list[@]}"; do
	freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 102 128 128 -nocursor -screenshot $t1dir/group_qc/1_recon/grota.png
	freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 128 128 128 -nocursor -screenshot $t1dir/group_qc/1_recon/grotb.png
	freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 154 128 128 -nocursor -screenshot $t1dir/group_qc/1_recon/grotc.png
	freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 128 128 102 -nocursor -screenshot $t1dir/group_qc/1_recon/grotd.png
	freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 128 128 128 -nocursor -screenshot $t1dir/group_qc/1_recon/grote.png
	freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 128 128 154 -nocursor -screenshot $t1dir/group_qc/1_recon/grotf.png
	freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 128 102 128 -nocursor -screenshot $t1dir/group_qc/1_recon/grotg.png
	freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 128 128 128 -nocursor -screenshot $t1dir/group_qc/1_recon/groth.png
	freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/aparc+aseg.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 128 154 128 -nocursor -screenshot $t1dir/group_qc/1_recon/groti.png
	pngappend grota.png + grotb.png + grotc.png + grotd.png + grote.png + grotf.png + grotg.png + groth.png + groti.png $subj.png
	echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> $t1dir/group_qc/1_recon/index.html
done

echo '</BODY></HTML>' >> $t1dir/group_qc/1_recon/index.html

mkdir -p $dwidir/group_qc/intermediate_nifti
mkdir -p $dwidir/group_qc/metrics

nii_list=( dwi_raw noisemap residual intermediate_nifti/1_dwi_denoised intermediate_nifti/2_dwi_degibbs intermediate_nifti/2_dwi_undistorted intermediate_nifti/3_dwi_smoothed intermediate_nifti/4_dwi_rician B0 B1000 B2000 metrics/dti_md metrics/dti_rd metrics/dti_ad metrics/dti_fa metrics/dki_mk metrics/dki_rk metrics/dki_ak metrics/dki_kfa metrics/dki_mkt )
for i in "${nii_list[@]}"; do
	slicesdir $dwidir/*/${i}.nii
	rm -rf $dwidir/group_qc/$i
	mv slicesdir $dwidir/group_qc/$i
done

mask_list=( brain_mask csf_mask FWF_mask ODI_mask )
for mask in "${mask_list[@]}"; do
	mkdir -p $dwidir/group_qc/$mask
	for subj in "${subj_list[@]}"; do
		echo $subj
		cp $dwidir/$subj/B0.nii $dwidir/group_qc/$mask/${subj}_1_B0.nii
		cp $dwidir/$subj/$mask.nii $dwidir/group_qc/$mask/${subj}_2_${mask}.nii
	done
	slicesdir -o $dwidir/group_qc/$mask/*.nii
	rm -rf $dwidir/group_qc/$mask/
	mv slicesdir $dwidir/group_qc/$mask
done

noddi_list=( AMICO/NODDI/fit_FWF AMICO/NODDI/fit_NDI AMICO/NODDI/fit_ODI )
mkdir -p $dwidir/group_qc/AMICO/NODDI
for  i in "${noddi_list[@]}"; do
	slicesdir $dwidir/*/${i}.nii.gz
	rm -rf $dwidir/group_qc/$i
	mv slicesdir $dwidir/group_qc/$i
done

mkdir -p $t1dir/group_qc/2_aparc2diff
echo '<HTML><TITLE>aparc2diff</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/2_aparc2diff/index.html

cd $t1dir/group_qc/2_aparc2diff
for subj in "${subj_list[@]}"; do
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 38 48 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grota.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grotb.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 58 48 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grotc.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 38 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grotd.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grote.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 58 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grotf.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 26 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grotg.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/groth.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 40 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/groti.png
	pngappend grota.png + grotb.png + grotc.png + grotd.png + grote.png + grotf.png + grotg.png + groth.png + groti.png $subj.png
	echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> $t1dir/group_qc/2_aparc2diff/index.html
done

echo '</BODY></HTML>' >> $t1dir/group_qc/2_aparc2diff/index.html