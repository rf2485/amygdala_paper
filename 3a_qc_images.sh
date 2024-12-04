# basedir=/gpfs/data/lazarlab/CamCan995/
basedir=/Volumes/Research/lazarm03lab/labspace/AD/camcan995/
projectdir=$basedir/derivatives/scd/main/
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

#generate QC webpage for freesurfer recon (WM and pial edges)
if [ ! -f $t1dir/group_qc/recon/index.html ]; then
	mkdir -p $t1dir/group_qc/recon
	echo '<HTML><TITLE>recon</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/recon/index.html
	for subj in "${subj_list[@]}"; do
		if [ ! -f $t1dir/group_qc/recon/$subj.png ]; then
			freeview -v $t1dir/$subj/mri/T1.mgz -f $t1dir/$subj/surf/lh.white:edgecolor=blue $t1dir/$subj/surf/rh.white:edgecolor=blue $t1dir/$subj/surf/lh.pial:edgecolor=red $t1dir/$subj/surf/rh.pial:edgecolor=red -viewport 'x' -slice 102 128 128 -nocursor -screenshot $t1dir/group_qc/recon/grota.png
			freeview -v $t1dir/$subj/mri/T1.mgz -f $t1dir/$subj/surf/lh.white:edgecolor=blue $t1dir/$subj/surf/rh.white:edgecolor=blue $t1dir/$subj/surf/lh.pial:edgecolor=red $t1dir/$subj/surf/rh.pial:edgecolor=red -viewport 'x' -slice 128 128 128 -nocursor -screenshot $t1dir/group_qc/recon/grotb.png
			freeview -v $t1dir/$subj/mri/T1.mgz -f $t1dir/$subj/surf/lh.white:edgecolor=blue $t1dir/$subj/surf/rh.white:edgecolor=blue $t1dir/$subj/surf/lh.pial:edgecolor=red $t1dir/$subj/surf/rh.pial:edgecolor=red -viewport 'x' -slice 154 128 128 -nocursor -screenshot $t1dir/group_qc/recon/grotc.png
			freeview -v $t1dir/$subj/mri/T1.mgz -f $t1dir/$subj/surf/lh.white:edgecolor=blue $t1dir/$subj/surf/rh.white:edgecolor=blue $t1dir/$subj/surf/lh.pial:edgecolor=red $t1dir/$subj/surf/rh.pial:edgecolor=red -viewport 'y' -slice 128 128 102 -nocursor -screenshot $t1dir/group_qc/recon/grotd.png
			freeview -v $t1dir/$subj/mri/T1.mgz -f $t1dir/$subj/surf/lh.white:edgecolor=blue $t1dir/$subj/surf/rh.white:edgecolor=blue $t1dir/$subj/surf/lh.pial:edgecolor=red $t1dir/$subj/surf/rh.pial:edgecolor=red -viewport 'y' -slice 128 128 128 -nocursor -screenshot $t1dir/group_qc/recon/grote.png
			freeview -v $t1dir/$subj/mri/T1.mgz -f $t1dir/$subj/surf/lh.white:edgecolor=blue $t1dir/$subj/surf/rh.white:edgecolor=blue $t1dir/$subj/surf/lh.pial:edgecolor=red $t1dir/$subj/surf/rh.pial:edgecolor=red -viewport 'y' -slice 128 128 154 -nocursor -screenshot $t1dir/group_qc/recon/grotf.png
			freeview -v $t1dir/$subj/mri/T1.mgz -f $t1dir/$subj/surf/lh.white:edgecolor=blue $t1dir/$subj/surf/rh.white:edgecolor=blue $t1dir/$subj/surf/lh.pial:edgecolor=red $t1dir/$subj/surf/rh.pial:edgecolor=red -viewport 'z' -slice 128 102 128 -nocursor -screenshot $t1dir/group_qc/recon/grotg.png
			freeview -v $t1dir/$subj/mri/T1.mgz -f $t1dir/$subj/surf/lh.white:edgecolor=blue $t1dir/$subj/surf/rh.white:edgecolor=blue $t1dir/$subj/surf/lh.pial:edgecolor=red $t1dir/$subj/surf/rh.pial:edgecolor=red -viewport 'z' -slice 128 128 128 -nocursor -screenshot $t1dir/group_qc/recon/groth.png
			freeview -v $t1dir/$subj/mri/T1.mgz -f $t1dir/$subj/surf/lh.white:edgecolor=blue $t1dir/$subj/surf/rh.white:edgecolor=blue $t1dir/$subj/surf/lh.pial:edgecolor=red $t1dir/$subj/surf/rh.pial:edgecolor=red -viewport 'z' -slice 128 154 128 -nocursor -screenshot $t1dir/group_qc/recon/groti.png
			pngappend $t1dir/group_qc/recon/grota.png + $t1dir/group_qc/recon/grotb.png + $t1dir/group_qc/recon/grotc.png + $t1dir/group_qc/recon/grotd.png + $t1dir/group_qc/recon/grote.png + $t1dir/group_qc/recon/grotf.png + $t1dir/group_qc/recon/grotg.png + $t1dir/group_qc/recon/groth.png + $t1dir/group_qc/recon/groti.png $t1dir/group_qc/recon/$subj.png
			echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> $t1dir/group_qc/recon/index.html
		fi
	done
	echo '</BODY></HTML>' >> $t1dir/group_qc/recon/index.html
fi

#generate QC webpage for JHU in T1 subject space
if [ ! -f $t1dir/group_qc/jhu2fs/index.html ]; then
	mkdir -p $t1dir/group_qc/jhu2fs
	echo '<HTML><TITLE>jhu2fs</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/jhu2fs/index.html
	for subj in "${subj_list[@]}"; do
		if [ ! -f $t1dir/group_qc/jhu2fs/$subj.png ]; then
			freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/jhu2fs.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 102 128 128 -nocursor -screenshot $t1dir/group_qc/jhu2fs/grota.png
			freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/jhu2fs.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 128 128 128 -nocursor -screenshot $t1dir/group_qc/jhu2fs/grotb.png
			freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/jhu2fs.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 154 128 128 -nocursor -screenshot $t1dir/group_qc/jhu2fs/grotc.png
			freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/jhu2fs.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 128 128 102 -nocursor -screenshot $t1dir/group_qc/jhu2fs/grotd.png
			freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/jhu2fs.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 128 128 128 -nocursor -screenshot $t1dir/group_qc/jhu2fs/grote.png
			freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/jhu2fs.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 128 128 154 -nocursor -screenshot $t1dir/group_qc/jhu2fs/grotf.png
			freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/jhu2fs.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 128 102 128 -nocursor -screenshot $t1dir/group_qc/jhu2fs/grotg.png
			freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/jhu2fs.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 128 128 128 -nocursor -screenshot $t1dir/group_qc/jhu2fs/groth.png
			freeview -v $t1dir/$subj/mri/T1.mgz $t1dir/$subj/mri/jhu2fs.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 128 154 128 -nocursor -screenshot $t1dir/group_qc/jhu2fs/groti.png
			pngappend $t1dir/group_qc/jhu2fs/grota.png + $t1dir/group_qc/jhu2fs/grotb.png + $t1dir/group_qc/jhu2fs/grotc.png + $t1dir/group_qc/jhu2fs/grotd.png + $t1dir/group_qc/jhu2fs/grote.png + $t1dir/group_qc/jhu2fs/grotf.png + $t1dir/group_qc/jhu2fs/grotg.png + $t1dir/group_qc/jhu2fs/groth.png + $t1dir/group_qc/jhu2fs/groti.png $t1dir/group_qc/jhu2fs/$subj.png
			echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> $t1dir/group_qc/jhu2fs/index.html
		fi
	done
	echo '</BODY></HTML>' >> $t1dir/group_qc/jhu2fs/index.html
fi

#generate QC webpages for diffusion processing
mkdir -p $dwidir/group_qc/intermediate_nifti
mkdir -p $dwidir/group_qc/metrics

nii_list=( dwi_raw noisemap residual intermediate_nifti/1_dwi_denoised intermediate_nifti/2_dwi_degibbs intermediate_nifti/2_dwi_undistorted intermediate_nifti/3_dwi_smoothed intermediate_nifti/4_dwi_rician B0 B1000 B2000 metrics/dti_md metrics/dti_rd metrics/dti_ad metrics/dti_fa metrics/dki_mk metrics/dki_rk metrics/dki_ak metrics/dki_kfa metrics/dki_mkt metrics/smi_matlab_Da metrics/smi_matlab_DePar metrics/smi_matlab_DePerp metrics/smi_matlab_f metrics/smi_matlab_p2 )
for i in "${nii_list[@]}"; do
	if [ ! -f $dwidir/group_qc/$i/index.html ]; then
		slicesdir $dwidir/*/${i}.nii
		rm -rf $dwidir/group_qc/$i
		mv slicesdir $dwidir/group_qc/$i
	fi
done

#generate QC webpages for diffusion masks
mask_list=( brain_mask csf_mask FWF_mask ODI_mask )
for mask in "${mask_list[@]}"; do
	if [ ! -f $dwidir/group_qc/$mask/index.html ]; then
		mkdir -p $dwidir/group_qc/$mask
		for subj in "${subj_list[@]}"; do
			echo $subj
			cp $dwidir/$subj/B0.nii $dwidir/group_qc/$mask/${subj}_1_B0.nii
			cp $dwidir/$subj/$mask.nii $dwidir/group_qc/$mask/${subj}_2_${mask}.nii
		done
		slicesdir -o $dwidir/group_qc/$mask/*.nii
		rm -rf $dwidir/group_qc/$mask/
		mv slicesdir $dwidir/group_qc/$mask
	fi
done

#generate QC webpages for NODDI metrics (nii.gz)
noddi_list=( AMICO/NODDI/fit_FWF AMICO/NODDI/fit_NDI AMICO/NODDI/fit_ODI )
mkdir -p $dwidir/group_qc/AMICO/NODDI
for i in "${noddi_list[@]}"; do
	if [ ! -f $dwidir/group_qc/$i/index.html ]; then
		slicesdir $dwidir/*/${i}.nii.gz
		rm -rf $dwidir/group_qc/$i
		mv slicesdir $dwidir/group_qc/$i
	fi
done

#generate QC webpage for aparc+aseg2diff
if [ ! -f $t1dir/group_qc/aparc+aseg2diff/index.html ]; then
	mkdir -p $t1dir/group_qc/aparc+aseg2diff
	echo '<HTML><TITLE>aparc+aseg2diff</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/aparc+aseg2diff/index.html
	for subj in "${subj_list[@]}"; do
		if [ ! -f $t1dir/group_qc/aparc+aseg2diff/$subj.png ]; then
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 38 48 33 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2diff/grota.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2diff/grotb.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 58 48 33 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2diff/grotc.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 38 33 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2diff/grotd.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2diff/grote.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 58 33 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2diff/grotf.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 26 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2diff/grotg.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2diff/groth.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 40 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2diff/groti.png
			pngappend $t1dir/group_qc/aparc+aseg2diff/grota.png + $t1dir/group_qc/aparc+aseg2diff/grotb.png + $t1dir/group_qc/aparc+aseg2diff/grotc.png + $t1dir/group_qc/aparc+aseg2diff/grotd.png + $t1dir/group_qc/aparc+aseg2diff/grote.png + $t1dir/group_qc/aparc+aseg2diff/grotf.png + $t1dir/group_qc/aparc+aseg2diff/grotg.png + $t1dir/group_qc/aparc+aseg2diff/groth.png + $t1dir/group_qc/aparc+aseg2diff/groti.png $t1dir/group_qc/aparc+aseg2diff/$subj.png
			echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> $t1dir/group_qc/aparc+aseg2diff/index.html
		fi
	done
	echo '</BODY></HTML>' >> $t1dir/group_qc/aparc+aseg2diff/index.html
fi

#generate QC webpage for jhu2diff
if [ ! -f $t1dir/group_qc/jhu2diff/index.html ]; then
	mkdir -p $t1dir/group_qc/jhu2diff
	echo '<HTML><TITLE>jhu2diff</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/jhu2diff/index.html
	for subj in "${subj_list[@]}"; do
		if [ ! -f $t1dir/group_qc/jhu2diff/$subj.png ]; then
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/jhu2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 38 48 33 -nocursor -screenshot $t1dir/group_qc/jhu2diff/grota.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/jhu2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/jhu2diff/grotb.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/jhu2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 58 48 33 -nocursor -screenshot $t1dir/group_qc/jhu2diff/grotc.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/jhu2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 38 33 -nocursor -screenshot $t1dir/group_qc/jhu2diff/grotd.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/jhu2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/jhu2diff/grote.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/jhu2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 58 33 -nocursor -screenshot $t1dir/group_qc/jhu2diff/grotf.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/jhu2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 26 -nocursor -screenshot $t1dir/group_qc/jhu2diff/grotg.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/jhu2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/jhu2diff/groth.png
			freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/jhu2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 40 -nocursor -screenshot $t1dir/group_qc/jhu2diff/groti.png
			pngappend $t1dir/group_qc/jhu2diff/grota.png + $t1dir/group_qc/jhu2diff/grotb.png + $t1dir/group_qc/jhu2diff/grotc.png + $t1dir/group_qc/jhu2diff/grotd.png + $t1dir/group_qc/jhu2diff/grote.png + $t1dir/group_qc/jhu2diff/grotf.png + $t1dir/group_qc/jhu2diff/grotg.png + $t1dir/group_qc/jhu2diff/groth.png + $t1dir/group_qc/jhu2diff/groti.png $t1dir/group_qc/jhu2diff/$subj.png
			echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> $t1dir/group_qc/jhu2diff/index.html
		fi
	done
	echo '</BODY></HTML>' >> $t1dir/group_qc/jhu2diff/index.html
fi

#generate QC webpages for MTI processing
if [ ! -f $t1dir/group_qc/1_mti_reg/index.html ]; then
	mkdir -p $t1dir/group_qc/1_mti_reg
	for subj in "${subj_list[@]}"; do
		echo $subj
		cp $t1dir/$subj/mti/intermediate/1_mti_reg.nii.gz $t1dir/group_qc/1_mti_reg/${subj}_1_mti_reg.nii.gz
		cp $t1dir/$subj/mti/mti_bl_raw.nii.gz $t1dir/group_qc/1_mti_reg/${subj}_2_mti_bl_raw.nii.gz
	done
	slicesdir -o $t1dir/group_qc/1_mti_reg/*.nii.gz
	rm -rf $t1dir/group_qc/1_mti_reg
	mv slicesdir $t1dir/group_qc/1_mti_reg
fi

if [ ! -f $t1dir/group_qc/2_mti_bl_degibbs/index.html ]; then
	slicesdir $t1dir/*/mti/intermediate/2_mti_bl_degibbs.nii.gz
	rm -rf $t1dir/group_qc/2_mti_bl_degibbs
	mv slicesdir $t1dir/group_qc/2_mti_bl_degibbs
fi

if [ ! -f $t1dir/group_qc/2_mti_degibbs/index.html ];then
	slicesdir $t1dir/*/mti/intermediate/2_mti_degibbs.nii.gz
	rm -rf $t1dir/group_qc/2_mti_degibbs
	mv slicesdir $t1dir/group_qc/2_mti_degibbs
fi

if [ ! -f $t1dir/group_qc/3_mtr/index.html ]; then
	slicesdir $t1dir/*/mti/mtr.nii.gz
	rm -rf $t1dir/group_qc/3_mtr
	mv slicesdir $t1dir/group_qc/3_mtr
fi

#generate QC webpage for aparc+aseg2mtr
if [ ! -f $t1dir/group_qc/aparc+aseg2mtr/index.html ]; then
	mkdir -p $t1dir/group_qc/aparc+aseg2mtr
	echo '<HTML><TITLE>aparc+aseg2mtr</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/aparc+aseg2mtr/index.html
	for subj in "${subj_list[@]}"; do
		if [ ! -f $t1dir/group_qc/aparc+aseg2mtr/$subj.png ]; then
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 42 64 64 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2mtr/grota.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 52 64 64 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2mtr/grotb.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 62 64 64 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2mtr/grotc.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 51 64 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2mtr/grotd.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 64 64 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2mtr/grote.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 77 64 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2mtr/grotf.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 51 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2mtr/grotg.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 64 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2mtr/groth.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 77 -nocursor -screenshot $t1dir/group_qc/aparc+aseg2mtr/groti.png
			pngappend $t1dir/group_qc/aparc+aseg2mtr/grota.png + $t1dir/group_qc/aparc+aseg2mtr/grotb.png + $t1dir/group_qc/aparc+aseg2mtr/grotc.png + $t1dir/group_qc/aparc+aseg2mtr/grotd.png + $t1dir/group_qc/aparc+aseg2mtr/grote.png + $t1dir/group_qc/aparc+aseg2mtr/grotf.png + $t1dir/group_qc/aparc+aseg2mtr/grotg.png + $t1dir/group_qc/aparc+aseg2mtr/groth.png + $t1dir/group_qc/aparc+aseg2mtr/groti.png $t1dir/group_qc/aparc+aseg2mtr/$subj.png
			echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> $t1dir/group_qc/aparc+aseg2mtr/index.html
		fi
	done
	echo '</BODY></HTML>' >> group_qc/aparc+aseg2mtr/index.html
fi

#generate QC webpage for jhu2mtr
if [ ! -f $t1dir/group_qc/jhu2mtr/index.html ]; then
	mkdir -p $t1dir/group_qc/jhu2mtr
	echo '<HTML><TITLE>jhu2mtr</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/jhu2mtr/index.html
	for subj in "${subj_list[@]}"; do
		if [ ! -f $t1dir/group_qc/jhu2mtr/$subj.png ]; then
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/jhu2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 42 64 64 -nocursor -screenshot $t1dir/group_qc/jhu2mtr/grota.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/jhu2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 52 64 64 -nocursor -screenshot $t1dir/group_qc/jhu2mtr/grotb.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/jhu2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 62 64 64 -nocursor -screenshot $t1dir/group_qc/jhu2mtr/grotc.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/jhu2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 51 64 -nocursor -screenshot $t1dir/group_qc/jhu2mtr/grotd.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/jhu2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 64 64 -nocursor -screenshot $t1dir/group_qc/jhu2mtr/grote.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/jhu2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 77 64 -nocursor -screenshot $t1dir/group_qc/jhu2mtr/grotf.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/jhu2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 51 -nocursor -screenshot $t1dir/group_qc/jhu2mtr/grotg.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/jhu2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 64 -nocursor -screenshot $t1dir/group_qc/jhu2mtr/groth.png
			freeview -v $t1dir/$subj/mti/mtr.nii.gz $t1dir/$subj/mti/jhu2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 77 -nocursor -screenshot $t1dir/group_qc/jhu2mtr/groti.png
			pngappend $t1dir/group_qc/jhu2mtr/grota.png + $t1dir/group_qc/jhu2mtr/grotb.png + $t1dir/group_qc/jhu2mtr/grotc.png + $t1dir/group_qc/jhu2mtr/grotd.png + $t1dir/group_qc/jhu2mtr/grote.png + $t1dir/group_qc/jhu2mtr/grotf.png + $t1dir/group_qc/jhu2mtr/grotg.png + $t1dir/group_qc/jhu2mtr/groth.png + $t1dir/group_qc/jhu2mtr/groti.png $t1dir/group_qc/jhu2mtr/$subj.png
			echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> $t1dir/group_qc/jhu2mtr/index.html
		fi
	done
	echo '</BODY></HTML>' >> $t1dir/group_qc/jhu2mtr/index.html
fi

#generate QC webpage for mtr2diff
subj_list=$(cut -f1 -d$'\t' $projectdir/mti_over_55.tsv)
subj_list=($subj_list)

if [ ! -f $t1dir/group_qc/mtr2diff/index.html ]; then
	mkdir -p $t1dir/group_qc/mtr2diff
	for subj in "${subj_list[@]}"; do
		echo $subj
		cp $t1dir/$subj/diffusion/B0.nii $t1dir/group_qc/mtr2diff/${subj}_1_B0.nii
		cp $t1dir/$subj/diffusion/mtr2diff.nii $t1dir/group_qc/mtr2diff/${subj}_2_mtr2diff.nii
	done
	cd $t1dir
	slicesdir -o group_qc/mtr2diff/*.nii
	rm -rf $t1dir/group_qc/mtr2diff
	mv slicesdir $t1dir/group_qc/mtr2diff
fi

#generate QC webpage for g-ratio
if [ ! -f $t1dir/group_qc/g_ratio/index.html ]; then
	slicesdir $t1dir/*/diffusion/g_ratio.nii
	rm -rf $t1dir/group_qc/g_ratio
	mv slicesdir $t1dir/group_qc/g_ratio
fi

