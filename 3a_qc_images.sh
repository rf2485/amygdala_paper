# basedir=/gpfs/data/lazarlab/CamCan995/
basedir=/Volumes/Research/lazarm03lab/labspace/AD/camcan995/
projectdir=$basedir/derivatives/scd/gm_roi/
t1dir=$projectdir/freesurfer/
dwidir=$projectdir/dwi_processed/

# export FREESURFER_HOME=/gpfs/share/apps/freesurfer/7.4.1/raw/freesurfer/
# module load freesurfer/7.4.1
export SUBJECTS_DIR=$t1dir
# module load miniconda3/gpu/4.9.2
# source activate ~/.conda/envs/fsl_eddy/
# export FSLDIR=$CONDA_PREFIX
# source $FSLDIR/etc/fslconf/fsl.sh

subj_list=$(cut -f1 -d$'\t' $projectdir/dwi_over_55.tsv)
subj_list=($subj_list)

cd $t1dir
mkdir -p $t1dir/group_qc/recon
echo '<HTML><TITLE>recon</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/recon/index.html

cd $t1dir/group_qc/1_recon
for subj in "${subj_list[@]}"; do
	if [ ! -f group_qc/recon/$subj.png ]; then
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'x' -slice 102 128 128 -nocursor -screenshot group_qc/recon/grota.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'x' -slice 128 128 128 -nocursor -screenshot group_qc/recon/grotb.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'x' -slice 154 128 128 -nocursor -screenshot group_qc/recon/grotc.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'y' -slice 128 128 102 -nocursor -screenshot group_qc/recon/grotd.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'y' -slice 128 128 128 -nocursor -screenshot group_qc/recon/grote.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'y' -slice 128 128 154 -nocursor -screenshot group_qc/recon/grotf.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'z' -slice 128 102 128 -nocursor -screenshot group_qc/recon/grotg.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'z' -slice 128 128 128 -nocursor -screenshot group_qc/recon/groth.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'z' -slice 128 154 128 -nocursor -screenshot group_qc/recon/groti.png
		pngappend group_qc/recon/grota.png + group_qc/recon/grotb.png + group_qc/recon/grotc.png + group_qc/recon/grotd.png + group_qc/recon/grote.png + group_qc/recon/grotf.png + group_qc/recon/grotg.png + group_qc/recon/groth.png + group_qc/recon/groti.png group_qc/recon/$subj.png
	fi
	echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> group_qc/recon/index.html
done
echo '</BODY></HTML>' >> group_qc/recon/index.html

cd $dwidir
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

mkdir -p $t1dir/group_qc/aparc+aseg2diff
echo '<HTML><TITLE>aparc2diff</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/aparc+aseg2diff/index.html
cd $t1dir/group_qc/aparc+aseg2diff
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
echo '</BODY></HTML>' >> $t1dir/group_qc/aparc+aseg2diff/index.html

cd $t1dir
rm -rf group_qc/1_mti_reg
mkdir -p group_qc/1_mti_reg
for subj in "${subj_list[@]}"; do
# for subj in "${(f)subj_list}"; do
	echo $subj
	cp $subj/mti/intermediate/1_mti_reg.nii.gz group_qc/1_mti_reg/${subj}_1_mti_reg.nii.gz
	cp $subj/mti/mti_bl_raw.nii.gz group_qc/1_mti_reg/${subj}_2_mti_bl_raw.nii.gz
done
slicesdir -o group_qc/1_mti_reg/*.nii.gz
rm -rf group_qc/1_mti_reg
mv slicesdir group_qc/1_mti_reg

slicesdir */mti/intermediate/2_mti_bl_degibbs.nii.gz
rm -rf group_qc/2_mti_bl_degibbs
mv slicesdir group_qc/2_mti_bl_degibbs

slicesdir */mti/intermediate/2_mti_degibbs.nii.gz
rm -rf group_qc/2_mti_degibbs
mv slicesdir group_qc/2_mti_degibbs

slicesdir */mti/mtr.nii.gz
rm -rf group_qc/3_mtr
mv slicesdir group_qc/3_mtr

mkdir -p group_qc/aparc+aseg2mtr
echo '<HTML><TITLE>aparc+aseg2mtr</TITLE><BODY BGCOLOR="#aaaaff">' > group_qc/aparc+aseg2mtr/index.html
for subj in "${subj_list[@]}"; do
# for subj in "${(f)subj_list}"; do
	if [ ! -f group_qc/aparc+aseg2mtr/$subj.png ]; then
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 42 64 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grota.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 52 64 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grotb.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 62 64 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grotc.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 51 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grotd.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 64 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grote.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 77 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grotf.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 51 -nocursor -screenshot group_qc/aparc+aseg2mtr/grotg.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/groth.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 77 -nocursor -screenshot group_qc/aparc+aseg2mtr/groti.png
		pngappend group_qc/aparc+aseg2mtr/grota.png + group_qc/aparc+aseg2mtr/grotb.png + group_qc/aparc+aseg2mtr/grotc.png + group_qc/aparc+aseg2mtr/grotd.png + group_qc/aparc+aseg2mtr/grote.png + group_qc/aparc+aseg2mtr/grotf.png + group_qc/aparc+aseg2mtr/grotg.png + group_qc/aparc+aseg2mtr/groth.png + group_qc/aparc+aseg2mtr/groti.png group_qc/aparc+aseg2mtr/$subj.png
	fi
	echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> group_qc/aparc+aseg2mtr/index.html
done

echo '</BODY></HTML>' >> group_qc/aparc+aseg2mtr/index.html