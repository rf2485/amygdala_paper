#!/bin/bash
#SBATCH --mail-user=rf2485@nyulangone.org
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --array=1-326
#SBATCH -o ./slurm_output/3_image_pipeline/slurm-%A_%a.out

basedir=/gpfs/data/lazarlab/CamCan995/
projectdir=$basedir/derivatives/scd/main/

subj_list=$(cut -f1 -d$'\t' $projectdir/dwi_mti_over_55.tsv)
subj_list=($subj_list)
subj_num=$(($SLURM_ARRAY_TASK_ID))
j=${subj_list[$subj_num]}

rawdwi=$basedir/raw/$j/dwi/${j}_dwi
rawanat=$basedir/raw/$j/anat/${j}
rawmti=$basedir/raw/$j/mti/${j}
designerdir=$projectdir/dwi_processed/$j/
intdir=$designerdir/intermediate_nifti/
mkdir -p $intdir
synb0in=$designerdir/synb0/INPUTS/
mkdir -p $synb0in
synb0out=$designerdir/synb0/OUTPUTS/
mkdir -p $synb0out
freesurferdir=$projectdir/freesurfer
mkdir -p $freesurferdir

module load freesurfer/7.4.1
export SUBJECTS_DIR=$freesurferdir
module load singularity/3.9.8
module load matlab/R2024a
module load miniconda3/gpu/4.9.2
source activate ~/.conda/envs/fsl_eddy/
export FSLDIR=$CONDA_PREFIX
source $FSLDIR/etc/fslconf/fsl.sh


#recon-all-clinical to generate segmentations
if [ ! -f $freesurferdir/$j/mri/wmparc.mgz ] && [ -f ${rawanat}_T1w.nii.gz ]; then
	rm $freesurferdir/$j/scripts/IsRunning.lh+rh
	recon-all-clinical.sh ${rawanat}_T1w.nii.gz $j 4 $freesurferdir
fi

if [ ! -f $freesurferdir/$j/label/BA_exvivo.thresh.ctab ] && [ -f ${rawanat}_T1w.nii.gz ]; then
	mkdir -p $freesurferdir/$j/mri/orig
	mri_convert ${rawanat}_T1w.nii.gz $freesurferdir/$j/mri/orig/001.mgz
	##create files necessary for bbregister
	recon-all -motioncor -talairach -nuintensitycor -normalization -gcareg -pctsurfcon -segstats -wmparc -balabels -subjid ${j} -sd ${freesurferdir} -threads 4 -no-isrunning
fi

if [ ! -f $freesurferdir/$j/stats/rh.AD_sig_thickness.stats ] && [ -f ${rawanat}_T1w.nii.gz ]; then
	##AD signature regions 
	cd $freesurferdir/$j
	side=( lh rh )
	for s in "${side[@]}"; do
		mri_annotation2label --hemi ${s} --subject ${j} --annotation aparc.a2009s --outdir a2009s
		mkdir ${s}_AD_sig
		AD_sig_list=( G_front_sup G_oc-temp_med-Parahip G_pariet_inf-Angular G_pariet_inf-Supramar G_parietal_sup G_precuneus G_temporal_inf Pole_temporal S_front_inf )
		for roi in "${AD_sig_list[@]}"; do
			cp a2009s/$s.$roi.label ${s}_AD_sig
		done
		mri_mergelabels -d ${s}_AD_sig -o label/$s.AD_sig.label
		mri_segstats --slabel ${j} ${s} label/$s.AD_sig.label --i surf/$s.thickness --excludeid 0 --sum stats/$s.AD_sig_thickness.stats
	done
	rm -rf a2009s
	rm -rf *h_AD_sig
fi

#diffusion preproc step 1
if [ ! -f $designerdir/residual.nii.gz ] && [ -f $rawdwi.nii.gz ]; then
	##concatenate and convert raw dwi
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert -json_import $rawdwi.json -fslgrad $rawdwi.bvec $rawdwi.bval $rawdwi.nii.gz $designerdir/working.mif -force
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert $designerdir/working.mif -json_export $designerdir/dwi_raw.json -export_grad_fsl $designerdir/dwi_raw.bvec $designerdir/dwi_raw.bval $designerdir/dwi_raw.nii -force

	##denoise, degibbs, and normalize with designer
	rm -rf $designerdir/tmp
	singularity exec --nv --bind $basedir $basedir/designer2_v2.0.10.sif designer -denoise -shrinkage frob -adaptive_patch -degibbs -nthreads 4 -nocleanup -scratch $designerdir/tmp $designerdir/dwi_raw.nii $intdir/2_dwi_degibbs.nii

	##move denoised files and calculate residual
	mv $designerdir/tmp/sigma.nii $designerdir/noisemap.nii
	mv $designerdir/tmp/tmp_dwidn.nii $intdir/1_dwi_denoised.nii
	cp $designerdir/dwi_raw.json $intdir/1_dwi_denoised.json
	cp $designerdir/dwi_raw.bvec $intdir/1_dwi_denoised.bvec
	cp $designerdir/dwi_raw.bval $intdir/1_dwi_denoised.bval
	singularity exec --bind $basedir $basedir/neurodock_v1.0.0.sif mrcalc $designerdir/dwi_raw.nii $designerdir/intermediate_nifti/1_dwi_denoised.nii -sub $designerdir/residual.nii.gz
fi

#diffusion preproc step 2
if [ ! -f $synb0out/b0_all.mif ] && [ -f $rawdwi.nii.gz ]; then
	##move and convert degibbsed files
	cp $designerdir/dwi_raw.json $intdir/2_dwi_degibbs.json
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert -json_import $intdir/2_dwi_degibbs.json -fslgrad $intdir/2_dwi_degibbs.bvec $intdir/2_dwi_degibbs.bval $intdir/2_dwi_degibbs.nii $designerdir/working.mif -force

	##prepare for synb0
	###create b0 and acqparams
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif dwiextract -bzero $designerdir/working.mif - | singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert - -coord 3 0 -axes 0,1,2 $intdir/b0_1_working.mif -force
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert -export_pe_table $synb0in/acqparams.txt $intdir/b0_1_working.mif $synb0in/b0.nii.gz -force
	rm $intdir/b0_1_working.mif
	echo " 0 1 0 0.000" >> $synb0in/acqparams.txt
	###copy T1
	cp ${rawanat}_T1w.nii.gz $synb0in/T1.nii.gz

	##synb0
	singularity run -e -B $synb0in:/INPUTS -B $synb0out:/OUTPUTS -B $basedir/license.txt:/extra/freesurfer/license.txt $basedir/synb0-disco_v3.1.sif
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert -import_pe_table $synb0in/acqparams.txt $synb0out/b0_all.nii.gz $synb0out/b0_all.mif -force
fi

#diffusion preproc step 3
if [ ! -f $intdir/2_dwi_undistorted.nii ] && [ -f $rawdwi.nii.gz ]; then
	##gather files for eddy_correct
	mkdir -p $designerdir/tmp/
	cp $synb0out/b0_all.mif $designerdir/tmp/se_epi.mif
	cd $designerdir/tmp
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert -json_import $intdir/2_dwi_degibbs.json -fslgrad $intdir/2_dwi_degibbs.bvec $intdir/2_dwi_degibbs.bval $intdir/2_dwi_degibbs.nii $designerdir/working.mif -force
	cp $designerdir/working.mif dwi.mif
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert dwi.mif dwi.nii -strides -1,+2,+3,+4 -json_export dwi.json -export_grad_fsl dwi.bvec dwi.bval

	##eddy_correct, using first volume as reference, with no blurring
	eddy_correct dwi dwi_post_eddy 0 -trilinear -noresampleblur
	###rotate bvecs using eddy_correct transformation log
	fdt_rotate_bvecs dwi.bvec dwi_post_eddy.bvec dwi_post_eddy.ecclog
	###convert eddy_correct output
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert dwi_post_eddy.nii.gz dwi_post_eddy.mif -json_import dwi.json -fslgrad dwi_post_eddy.bvec dwi.bval

	##quantify motion
	logfile=dwi_post_eddy.ecclog
	basenm=`basename $logfile .ecclog`
	nums=`grep -n 'Final' $logfile | sed 's/:.*//'`

	touch grot_ts.txt
	touch grot.mat

	firsttime=yes;
	m=1;
	for n in $nums ; do 
	    echo "Timepoint $m"
	    n1=`echo $n + 1 | bc` ; 
	    n2=`echo $n + 5 | bc` ;
	    sed -n  "$n1,${n2}p" $logfile > grot.mat ; 
	    if [ $firsttime = yes ] ; then firsttime=no; cp grot.mat grot.refmat ; cp grot.mat grot.oldmat ; fi
	    absval=`rmsdiff grot.mat grot.refmat $basenm`;
	    relval=`rmsdiff grot.mat grot.oldmat $basenm`;
	    cp grot.mat grot.oldmat
	    echo $absval $relval >> ec_disp.txt ;
	    avscale --allparams grot.mat $basenm | grep 'Rotation Angles' | sed 's/.* = //' >> ec_rot.txt ;
	    avscale --allparams grot.mat $basenm | grep 'Translations' | sed 's/.* = //' >> ec_trans.txt ;
	    m=`echo $m + 1 | bc`;
	done

	echo "absolute" > grot_labels.txt
	echo "relative" >> grot_labels.txt

	fsl_tsplot -i ec_disp.txt -t 'Eddy Current estimated mean displacement (mm)' -l grot_labels.txt -o ec_disp.png

	echo "x" > grot_labels.txt
	echo "y" >> grot_labels.txt
	echo "z" >> grot_labels.txt

	fsl_tsplot -i ec_rot.txt -t 'Eddy Current estimated rotations (radians)' -l grot_labels.txt -o ec_rot.png
	fsl_tsplot -i ec_trans.txt -t 'Eddy Current estimated translations (mm)' -l grot_labels.txt -o ec_trans.png

	mkdir -p $designerdir/metrics_qc/eddy
	cp ec* $designerdir/metrics_qc/eddy/
	cp dwi_post_eddy.ecclog $designerdir/metrics_qc/eddy/

	##topup
	###replace first b0 of se_epi with first b0 of dwi
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert dwi_post_eddy.mif dwi_first_bzero.mif -coord 3 0 -axes 0,1,2
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert se_epi.mif - -coord 3 1 | singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrcat dwi_first_bzero.mif - se_epi_firstdwibzero.mif -axis 3
	###create files for topup input
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert se_epi_firstdwibzero.mif topup_in.nii -strides -1,+2,+3,+4 -export_pe_table topup_datain.txt
	###run topup
	topup --imain=topup_in.nii --datain=topup_datain.txt --out=field --fout=field_map.nii.gz --config=${FSLDIR}/etc/flirtsch/b02b0.cnf --verbose
	
	##applytopup
	###create files for applytopup
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrinfo dwi_post_eddy.mif -export_pe_eddy applytopup_config.txt applytopup_indices.txt
	##run applytopup
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif applytopup --imain=dwi_post_eddy.nii.gz --datain=applytopup_config.txt --inindex=1 --topup=field --out=dwi_applytopup --method=jac
	###move and convert applytopup output
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert dwi_applytopup.nii.gz dwi_applytopup.mif -json_import dwi.json -fslgrad dwi_post_eddy.bvec dwi.bval
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrconvert -export_grad_fsl $intdir/2_dwi_undistorted.bvec $intdir/2_dwi_undistorted.bval -json_export $intdir/2_dwi_undistorted.json dwi_applytopup.mif $intdir/2_dwi_undistorted.nii -force
	cd $designerdir
	rm -rf $designerdir/tmp/
fi

#diffusion preproc step 4
if [ ! -f $designerdir/metrics/dki_dsistudio.fib ] && [ -f $rawdwi.nii.gz ]; then
	singularity exec --bind $basedir/ $basedir/neurodock_v1.0.0.sif mrconvert -json_import $intdir/2_dwi_undistorted.json -fslgrad $intdir/2_dwi_undistorted.bvec $intdir/2_dwi_undistorted.bval $intdir/2_dwi_undistorted.nii $designerdir/working.mif -force
	##brain mask
	mkdir -p $designerdir/tmp/
	cd $designerdir/tmp
	singularity exec --bind $basedir/ $basedir/neurodock_v1.0.0.sif dwiextract $designerdir/working.mif - -bzero | singularity exec --bind $basedir/ $basedir/neurodock_v1.0.0.sif mrmath - mean - -axis 3 | singularity exec --bind $basedir/ $basedir/neurodock_v1.0.0.sif mrconvert - meanb0.nii
	mri_synthstrip -i meanb0.nii -m new_brain_mask.nii
	##csf mask, csf-excluded smoothing, rician bias correction, and model fitting with pydesigner
	mkdir -p $designerdir/metrics_qc/fitting
	singularity exec --bind $basedir/ $basedir/neurodock_v1.0.0.sif pydesigner --resume -n -u -cd 2 -z -r --user_mask $designerdir/tmp/new_brain_mask.nii --verbose --noqc --nthreads 4 -o $designerdir/ $rawdwi.nii.gz
	cd $designerdir
	rm -rf $designerdir/tmp/
fi
if [ ! -f $designerdir/metrics/dti_V1.nii.gz ] && [ -f $rawdwi.nii.gz ]; then
	##calculate colormap from diffusion tensor
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif tensor2metric $designerdir/metrics/DT.nii -vector $designerdir/metrics/dti_V1.nii.gz -modulate none -force
fi
if [ ! -f $designerdir/AMICO/NODDI/fit_FWF.nii.gz ] && [ -f $rawdwi.nii.gz ]; then
	##fit NODDI model using the AMICO program
	conda deactivate
	source activate ~/.conda/envs/amico/ #switch to amico conda env
	python3 $projectdir/amico_noddi.py $j
	conda deactivate
	source activate ~/.conda/envs/fsl_eddy/ #switch back to FSL eddy conda env
	export FSLDIR=$CONDA_PREFIX
	source $FSLDIR/etc/fslconf/fsl.sh
fi
if [ ! -f $designerdir/FWF_mask.nii.gz ] && [ -f $rawdwi.nii.gz ]; then
	##mask voxels with orientation dispersion > 0.8 or free water > 0.7 (Correia et al.)
	cd $designerdir
	mri_binarize --i AMICO/NODDI/fit_ODI.nii.gz --max 0.8 --o ODI_mask.nii.gz
	mri_binarize --i AMICO/NODDI/fit_FWF.nii.gz --max 0.7 --o FWF_mask.nii.gz
fi

#generate diffusion Freesurfer stats
if [ ! -f $freesurferdir/$j/diffusion/b02fs.lta ] && [ -f $rawdwi.nii.gz ]; then
	cd $freesurferdir/$j/
	mkdir -p diffusion
	## copy B0 to freesurferdir
	mri_convert $designerdir/B0.nii diffusion/B0.nii.gz
	## register B0 to t1
	bbregister --s ${j} --mov diffusion/B0.nii.gz --reg diffusion/b02fs.lta --dti --init-fsl --9
fi
if [ ! -f $freesurferdir/$j/diffusion/mni2diff_field.mgz ] && [ -f $rawdwi.nii.gz ]; then
	## register MNI to B0
	cd $freesurferdir/$j/diffusion
	mri_easyreg --ref $FREESURFER_HOME/subjects/cvs_avg35_inMNI152/mri/norm.mgz --flo B0.nii.gz --ref_seg $FREESURFER_HOME/subjects/cvs_avg35_inMNI152/mri/aparc+aseg.mgz --flo_seg aparc+aseg2diff.mgz --bak_field mni2diff_field.mgz --threads 4
fi
if [ ! -f $freesurferdir/$j/diffusion/aparc+aseg2diff_NODDI_mask.mgz ] && [ -f $rawdwi.nii.gz ]; then
	## use inverse of b02fs registration to register aparc+aseg to diffusion, then apply NODDI masks	
	cd $freesurferdir/$j/
	mri_vol2vol --mov diffusion/B0.nii.gz --targ mri/aparc+aseg.mgz --inv --interp nearest --o diffusion/aparc+aseg2diff.mgz --reg diffusion/b02fs.lta --no-save-reg
	mri_mask diffusion/aparc+aseg2diff.mgz $designerdir/FWF_mask.nii.gz diffusion/aparc+aseg2diff_FWF_mask.mgz
	mri_mask diffusion/aparc+aseg2diff_FWF_mask.mgz $designerdir/ODI_mask.nii.gz diffusion/aparc+aseg2diff_NODDI_mask.mgz
fi
if [ ! -f $freesurferdir/$j/diffusion/lh_diff_AD_sig_NODDI_mask.mgz ] && [ -f $rawdwi.nii.gz ]; then
	## use aparc.a2009s+aseg to generate Dickerson AD signature regions mask, then apply NODDI masks
	cd $freesurferdir/$j/	
	mri_vol2vol --mov diffusion/B0.nii.gz --targ mri/aparc.a2009s+aseg.mgz --inv --interp nearest --o diffusion/aparc.a2009s+aseg2diff.mgz --reg diffusion/b02fs.lta --no-save-reg
	side=( rh lh )
	for s in "${side[@]}"; do
		mri_binarize --i diffusion/aparc.a2009s+aseg2diff.mgz --match-ctab $projectdir/${s}_AD_signature_LUT.txt --o diffusion/${s}_diff_AD_sig.mgz
		mri_mask diffusion/${s}_diff_AD_sig.mgz $designerdir/FWF_mask.nii.gz diffusion/${s}_diff_AD_sig_FWF_mask.mgz
		mri_mask diffusion/${s}_diff_AD_sig_FWF_mask.mgz $designerdir/ODI_mask.nii.gz diffusion/${s}_diff_AD_sig_NODDI_mask.mgz
	done  	
fi
if [ ! -f $freesurferdir/$j/diffusion/fit_ODI.nii.gz ] && [ -f $rawdwi.nii.gz ]; then
	## copy diffusion metrics to freesurferdir
	cd $freesurferdir/$j/
	nii_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd )
	for meas in "${nii_list[@]}"; do
		mri_convert $designerdir/metrics/${meas}.nii diffusion/${meas}.nii.gz
	done
	noddi_list=( fit_FWF fit_NDI fit_ODI )
	for meas in "${noddi_list[@]}"; do
		cp $designerdir/AMICO/NODDI/${meas}.nii.gz diffusion/
	done
fi
	
if [ ! -f $freesurferdir/$j/stats/lh_AD_sig2smi_matlab_p2.stats ] && [ -f $rawdwi.nii.gz ]; then	
	## apply registered aparc+aseg, AD signature mask, and JHU to diffusion metrics to generate stats
	cd $freesurferdir/$j/
	diff_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd fit_FWF fit_NDI fit_ODI )
	for meas in "${pyd_list[@]}"; do
		### generate stats
		mri_segstats --seg diffusion/aparc+aseg2diff_NODDI_mask.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i diffusion/${meas}.nii.gz --sum stats/aparc+aseg2${meas}.stats --excludeid 0 30 31 62 63 72 77 78 79 80 81 82 85 2 7 41 46 251 252 253 254 255
		side=( rh lh )
		for s in "${side[@]}"; do
			mri_segstats --seg diffusion/${s}_diff_AD_sig_NODDI_mask.mgz --mask diffusion/${s}_diff_AD_sig_NODDI_mask.mgz --excludeid 0 --ctab $projectdir/${s}_AD_signature_LUT.txt --i diffusion/${meas}.nii.gz --sum stats/${s}_AD_sig2${meas}.stats
		done
	done	
fi

# mtr analysis
if [ ! -f $freesurferdir/$j/mti/mti_raw.nii.gz ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	##copy raw mti images
	mkdir -p $freesurferdir/$j/mti/intermediate/
	cd $freesurferdir/$j/mti
	cp ${rawmti}_mti_bl.nii.gz mti_bl_raw.nii.gz
	cp  ${rawmti}_mti.nii.gz mti_raw.nii.gz
fi
if [ ! -f $freesurferdir/$j/mti/intermediate/1_mti_reg.nii.gz ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	##register mti_bl and mti
	cd $freesurferdir/$j/mti
	flirt -in mti_raw -ref mti_bl_raw -dof 6 -out intermediate/1_mti_reg
fi
if [ ! -f $freesurferdir/$j/mti/intermediate/2_mti_degibbs.nii.gz ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	#degibbs mti_bl and registered mti
	cd $freesurferdir/$j/mti
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrdegibbs mti_bl_raw.nii.gz intermediate/2_mti_bl_degibbs.nii.gz
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrdegibbs intermediate/1_mti_reg.nii.gz intermediate/2_mti_degibbs.nii.gz
fi
if [ ! -f $freesurferdir/$j/mti/mtr.nii.gz ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	cd $freesurferdir/$j/mti
	##calculate mtr
	fslmaths intermediate/2_mti_bl_degibbs.nii.gz -sub intermediate/2_mti_degibbs.nii.gz -div intermediate/2_mti_bl_degibbs.nii.gz -thr 0 -uthr 1 mtr.nii.gz
fi

#generate mtr Freesurfer stats
if [ ! -f $freesurferdir/$j/mti/fs2mtr_field.mgz ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	##register mtr to t1
	cd $freesurferdir/$j
	mri_easyreg --ref mri/T1.mgz --flo mti/mtr.nii.gz --ref_seg mri/synthseg.mgz --flo_seg mti/mtr_synthseg.mgz --bak_field mti/fs2mtr_field.mgz --threads 4 --affine_only
fi
if [ ! -f $freesurferdir/$j/mti/mni2mtr_field.mgz ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	## register MNI to MTR
	cd $freesurferdir/$j/mti
	mri_easyreg --ref $FREESURFER_HOME/subjects/cvs_avg35_inMNI152/mri/norm.mgz --flo mtr.nii.gz --ref_seg $FREESURFER_HOME/subjects/cvs_avg35_inMNI152/mri/aparc+aseg.mgz --flo_seg aparc+aseg2mtr.mgz --bak_field mni2mtr_field.mgz --threads 4
fi
if [ ! -f $freesurferdir/$j/mti/aparc+aseg2mtr.mgz ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	##register aparc+aseg to mtr, using the inverse of the mtr to t1 registration
	cd $freesurferdir/$j
	mri_easywarp --i mri/aparc+aseg.mgz --o mti/aparc+aseg2mtr.mgz --field mti/fs2mtr_field.mgz --threads 4 --nearest
fi
if [ ! -f $freesurferdir/$j/mti/lh_mtr_AD_sig.mgz ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then	
	## use aparc.a2009s+aseg to generate Dickerson AD signature regions mask
	cd $freesurferdir/$j
	mri_easywarp --i mri/aparc.a2009s+aseg.mgz --o mti/aparc.a2009s+aseg2mtr.mgz --field mti/fs2mtr_field.mgz --threads 4 --nearest
	side=( rh lh )
	for s in "${side[@]}"; do
		mri_binarize --i mti/aparc.a2009s+aseg2mtr.mgz --match-ctab $projectdir/${s}_AD_signature_LUT.txt --o mti/${s}_mtr_AD_sig.mgz
	done
fi
if [ ! -f $freesurferdir/$j/stats/lh_AD_sig2mtr.stats ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	## apply registered aparc+aseg, AD signature mask to mtr to generate stats
	cd $freesurferdir/$j
	mri_segstats --seg mti/aparc+aseg2mtr.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i mti/mtr.nii.gz --sum stats/aparc+aseg2mtr.stats --excludeid 0 30 31 62 63 72 77 78 79 80 81 82 85
	side=( rh lh )
	for s in "${side[@]}"; do
		mri_segstats --seg mti/${s}_mtr_AD_sig.mgz --mask mti/${s}_mtr_AD_sig.mgz --excludeid 0 --ctab $projectdir/${s}_AD_signature_LUT.txt --i mti/mtr.nii.gz --sum stats/${s}_AD_sig2mtr.stats
	done	
fi
if [ ! -f $freesurferdir/$j/stats/mtr_gm.stats ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	##whole brain WM and GM in MTR (to choose which TR has better contrast)
	cd $freesurferdir/$j
	mri_binarize --i mti/aparc+aseg2mtr.mgz --all-wm --o mti/mtr_wm.mgz #all wm
	mri_binarize --i mti/aparc+aseg2mtr.mgz --gm --o mti/mtr_gm.mgz #all gm
	mri_segstats --seg mti/mtr_wm.mgz --mask mti/mtr_wm.mgz --excludeid 0 --i mti/mtr.nii.gz --sum stats/mtr_wm.stats
	mri_segstats --seg mti/mtr_gm.mgz --mask mti/mtr_gm.mgz --excludeid 0 --i mti/mtr.nii.gz --sum stats/mtr_gm.stats
fi