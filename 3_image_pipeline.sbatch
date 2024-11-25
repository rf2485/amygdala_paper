#!/bin/bash
#SBATCH --mail-user=rf2485@nyulangone.org
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --array=1-325
#SBATCH -o ./slurm_output/3_image_pipeline/slurm-%A_%a.out

basedir=/gpfs/data/lazarlab/CamCan995/
projectdir=$basedir/derivatives/scd/main/

subj_list=$(cut -f1 -d$'\t' $projectdir/dwi_over_55.tsv)
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
module load miniconda3/gpu/4.9.2
source activate ~/.conda/envs/fsl_eddy/
export FSLDIR=$CONDA_PREFIX
source $FSLDIR/etc/fslconf/fsl.sh
module load matlab/R2024a

#recon-all-clinical to generate segmentations
if [ ! -f $freesurferdir/$j/mri/wmparc.mgz ]; then
	rm $freesurferdir/$j/scripts/IsRunning.lh+rh
	recon-all-clinical.sh ${rawanat}_T1w.nii.gz $j 4 $freesurferdir
fi

if [ ! -f $freesurferdir/$j/label/BA_exvivo.thresh.ctab ]; then
	mkdir -p $freesurferdir/$j/mri/orig
	mri_convert ${rawanat}_T1w.nii.gz $freesurferdir/$j/mri/orig/001.mgz
	##create files necessary for bbregister
	recon-all -motioncor -talairach -nuintensitycor -normalization -gcareg -pctsurfcon -segstats -wmparc -balabels -subjid ${j} -sd ${freesurferdir} -threads 4 -no-isrunning
fi
##AD signature regions 
if [ ! -f $freesurferdir/$j/stats/rh.AD_sig_thickness.stats ]; then
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
##JHU volumes
if [ ! -f $freesurferdir/$j/stats/jhu_volume.stats ]; then
	cd $freesurferdir/$j
	mri_easyreg --ref $FREESURFER_HOME/subjects/cvs_avg35_inMNI152/mri/norm.mgz --flo mri/norm.mgz --ref_seg $FREESURFER_HOME/subjects/cvs_avg35_inMNI152/mri/aparc+aseg.mgz --flo_seg mri/aparc+aseg.mgz --bak_field mri/mni2fs_field.mgz --threads 4
	mri_easywarp --i $FSLDIR/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz --o mri/jhu2fs.mgz --field mri/mni2fs_field.mgz --nearest --threads 4
	cd $freesurferdir/$j/
	mri_segstats --seg mri/jhu2fs.mgz --pv mri/norm.mgz --mask mri/jhu2fs.mgz --excludeid 0 --in mri/norm.mgz --sum stats/jhu_volume.stats
fi

#diffusion preproc step 1
if [ ! -f $designerdir/residual.nii ]; then
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
	singularity exec --bind $basedir $basedir/neurodock_v1.0.0.sif mrcalc $designerdir/dwi_raw.nii $designerdir/intermediate_nifti/1_dwi_denoised.nii -sub $designerdir/residual.nii
fi

#diffusion preproc step 2
if [ ! -f $synb0out/b0_all.mif ]; then
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
if [ ! -f $intdir/2_dwi_undistorted.nii ]; then
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
if [ ! -f $designerdir/metrics/dki_dsistudio.fib ]; then
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
if [ ! -f $designerdir/metrics/dti_V1.nii ]; then
	##calculate colormap from diffusion tensor
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif tensor2metric $designerdir/metrics/DT.nii -vector $designerdir/metrics/dti_V1.nii -modulate none -force
fi
if [ ! -f $designerdir/AMICO/NODDI/fit_FWF.nii.gz ]; then
	##fit NODDI model using the AMICO program
	conda deactivate
	source activate ~/.conda/envs/amico/
	python3 $projectdir/amico_noddi.py $j
fi
if [ ! -f $designerdir/metrics/smi_matlab_p2.nii ]; then
	##fit SMI model using matlab
	matlab -sd $designerdir -batch "addpath('${basedir}/SMI'); nii = niftiread('dwi_preprocessed.nii'); dwi = abs(double(nii)); bval = load('dwi_preprocessed.bval'); dirs = load('dwi_preprocessed.bvec'); nii_mask = niftiread('brain_mask.nii'); sigma_nii = niftiread('noisemap.nii'); options.b = bval; options.beta = []; options.dirs = dirs; options.TE   = []; options.MergeDistance = 0.1; options.mask = logical(nii_mask); options.sigma = abs(sigma_nii); options.compartments = {'IAS','EAS'}; options.NoiseBias = 'Rician'; options.MLTraining.bounds = [0.05, 1, 1, 0.1, 0, 50, 50, 0.05; 0.95, 3, 3, 1.2, 0.5, 150, 120, 0.99]; [out] = SMI.fit(dwi,options); niftiwrite(out.kernel(:,:,:,1), 'metrics/smi_matlab_f.nii'); niftiwrite(out.kernel(:,:,:,2), 'metrics/smi_matlab_Da.nii'); niftiwrite(out.kernel(:,:,:,3), 'metrics/smi_matlab_DePar.nii'); niftiwrite(out.kernel(:,:,:,4), 'metrics/smi_matlab_DePerp.nii'); niftiwrite(out.kernel(:,:,:,6), 'metrics/smi_matlab_p2.nii');"
	###copy header to SMI output images
	smi_list=( smi_matlab_f smi_matlab_Da smi_matlab_DePar smi_matlab_DePerp smi_matlab_p2 )
	for meas in "${smi_list[@]}"; do
		fslcpgeom $designerdir/dwi_preprocessed.nii $designerdir/metrics/${meas}.nii -d
	done	
fi

#mask voxels with orientation dispersion > 0.8 or free water > 0.7 (Correia et al.)
if [ ! -f $designerdir/FWF_mask.nii ]; then
	cd $designerdir
	mri_binarize --i AMICO/NODDI/fit_ODI.nii.gz --max 0.8 --o ODI_mask.nii
	mri_binarize --i AMICO/NODDI/fit_FWF.nii.gz --max 0.7 --o FWF_mask.nii
fi

#generate diffusion Freesurfer stats
if [ ! -f $freesurferdir/$j/diffusion/b02fs.lta ]; then
	cd $freesurferdir/$j/
	mkdir -p diffusion
	## copy B0 to freesurferdir
	cp $designerdir/B0.nii diffusion/
	## register B0 to t1
	bbregister --s ${j} --mov diffusion/B0.nii --reg diffusion/b02fs.lta --dti --init-fsl --9
fi
	## use inverse of b02fs registration to register aparc+aseg to diffusion, then apply NODDI masks
if [ ! -f $freesurferdir/$j/diffusion/aparc+aseg2diff_NODDI_mask.mgz ]; then	
	cd $freesurferdir/$j/
	mri_vol2vol --mov diffusion/B0.nii --targ mri/aparc+aseg.mgz --inv --interp nearest --o diffusion/aparc+aseg2diff.mgz --reg diffusion/b02fs.lta --no-save-reg
	mri_mask diffusion/aparc+aseg2diff.mgz $designerdir/FWF_mask.nii diffusion/aparc+aseg2diff_FWF_mask.mgz
	mri_mask diffusion/aparc+aseg2diff_FWF_mask.mgz $designerdir/ODI_mask.nii diffusion/aparc+aseg2diff_NODDI_mask.mgz
fi
	## use aparc.a2009s+aseg to generate Dickerson AD signature regions mask, then apply NODDI masks
if [ ! -f $freesurferdir/$j/diffusion/lh_diff_AD_sig_NODDI_mask.mgz ]; then
	cd $freesurferdir/$j/	
	mri_vol2vol --mov diffusion/B0.nii --targ mri/aparc.a2009s+aseg.mgz --inv --interp nearest --o diffusion/aparc.a2009s+aseg2diff.mgz --reg diffusion/b02fs.lta --no-save-reg
	side=( rh lh )
	for s in "${side[@]}"; do
		mri_binarize --i diffusion/aparc.a2009s+aseg2diff.mgz --match-ctab $projectdir/${s}_AD_signature_LUT.txt --o diffusion/${s}_diff_AD_sig.mgz
		mri_mask diffusion/${s}_diff_AD_sig.mgz $designerdir/FWF_mask.nii diffusion/${s}_diff_AD_sig_FWF_mask.mgz
		mri_mask diffusion/${s}_diff_AD_sig_FWF_mask.mgz $designerdir/ODI_mask.nii diffusion/${s}_diff_AD_sig_NODDI_mask.mgz
	done  	
fi	

	## apply registered aparc+aseg and AD signature mask to diffusion metrics to generate stats
if [ ! -f $freesurferdir/$j/stats/lh_AD_sig2fit_ODI.stats ]; then	
	cd $freesurferdir/$j/
	pyd_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd )
	for meas in "${pyd_list[@]}"; do
		### copy to freesurferdir
		cp $designerdir/metrics/${meas}.nii diffusion/${meas}.nii
		### generate stats
		mri_segstats --seg diffusion/aparc+aseg2diff_NODDI_mask.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i diffusion/${meas}.nii --sum stats/aparc+aseg2${meas}.stats --excludeid 0 30 31 62 63 72 77 78 79 80 81 82 85
		side=( rh lh )
		for s in "${side[@]}"; do
			mri_segstats --seg diffusion/${s}_diff_AD_sig_NODDI_mask.mgz --mask diffusion/${s}_diff_AD_sig_NODDI_mask.mgz --excludeid 0 --ctab $projectdir/${s}_AD_signature_LUT.txt --i diffusion/${meas}.nii --sum stats/${s}_AD_sig2${meas}.stats
		done
	done
	### repeat for NODDI (GM only)	
	noddi_list=( fit_FWF fit_NDI fit_ODI )
	for meas in "${noddi_list[@]}"; do
		cp $designerdir/AMICO/NODDI/${meas}.nii.gz diffusion/
		mri_segstats --seg diffusion/aparc+aseg2diff_NODDI_mask.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i diffusion/${meas}.nii.gz --sum stats/aparc+aseg2${meas}.stats --excludeid 0 30 31 62 63 72 77 78 79 80 81 82 85 2 7 41 46 251 252 253 254 255
		side=( rh lh )
		for s in "${side[@]}"; do
			mri_segstats --seg diffusion/${s}_diff_AD_sig_NODDI_mask.mgz --mask diffusion/${s}_diff_AD_sig_NODDI_mask.mgz --excludeid 0 --ctab $projectdir/${s}_AD_signature_LUT.txt --i diffusion/${meas}.nii.gz --sum stats/${s}_AD_sig2${meas}.stats
		done
	done
fi
### repeat for SMI (WM only)
if [ ! -f $freesurferdir/$j/stats/aparc+aseg2smi_matlab_p2.stats ]; then
	cd $freesurferdir/$j/
	smi_list=( smi_matlab_f smi_matlab_Da smi_matlab_DePar smi_matlab_DePerp smi_matlab_p2 )
	for meas in "${smi_list[@]}"; do
		cp $designerdir/metrics/${meas}.nii diffusion/${meas}.nii
		mri_segstats --seg diffusion/aparc+aseg2diff_NODDI_mask.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i diffusion/${meas}.nii --sum stats/aparc+aseg2${meas}.stats --id 2 7 41 46 251 252 253 254 255
	done
fi
#register JHU MNI atlas to diffusion images
if [ ! -f $freesurferdir/$j/stats/jhu2smi_matlab_p2.stats ]; then
	cd $freesurferdir/$j/diffusion
	mri_easyreg --ref $FREESURFER_HOME/subjects/cvs_avg35_inMNI152/mri/norm.mgz --flo B0.nii --ref_seg $FREESURFER_HOME/subjects/cvs_avg35_inMNI152/mri/aparc+aseg.mgz --flo_seg aparc+aseg2diff.mgz --bak_field mni2diff_field.mgz --threads 4
	mri_easywarp --i $FSLDIR/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz --o jhu2diff.mgz --field mni2diff_field.mgz --nearest --threads 4
	mri_mask jhu2diff.mgz $designerdir/FWF_mask.nii jhu2diff_FWF_mask.mgz
	mri_mask jhu2diff_FWF_mask.mgz $designerdir/ODI_mask.nii jhu2diff_NODDI_mask.mgz
	cd $freesurferdir/$j/
	diff_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd smi_matlab_f smi_matlab_Da smi_matlab_DePar smi_matlab_DePerp smi_matlab_p2 )
	for meas in "${diff_list[@]}"; do
		mri_segstats --seg diffusion/jhu2diff_NODDI_mask.mgz --mask diffusion/jhu2diff_NODDI_mask.mgz --i diffusion/$meas.nii --sum stats/jhu2${meas}.stats --excludeid 0
	done
fi

#mtr analysis
if [ ! -f $freesurferdir/$j/mti/mtr.nii.gz ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	#copy raw mti images
	mkdir -p $freesurferdir/$j/mti/intermediate/
	cd $freesurferdir/$j/mti
	cp ${rawmti}_mti_bl.nii.gz mti_bl_raw.nii.gz
	cp  ${rawmti}_mti.nii.gz mti_raw.nii.gz
	#register mti_bl and mti
	flirt -in mti_raw -ref mti_bl_raw -dof 6 -out intermediate/1_mti_reg
	#degibbs mti_bl and registered mti
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrdegibbs mti_bl_raw.nii.gz intermediate/2_mti_bl_degibbs.nii.gz
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrdegibbs intermediate/1_mti_reg.nii.gz intermediate/2_mti_degibbs.nii.gz
	#calculate mtr
	fslmaths intermediate/2_mti_bl_degibbs.nii.gz -sub intermediate/2_mti_degibbs.nii.gz -div intermediate/2_mti_bl_degibbs.nii.gz -thr 0 -uthr 1 mtr.nii.gz
fi

if [ ! -f $freesurferdir/$j/stats/aparc+aseg2mtr.stats ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	#register mtr to t1
	cd $freesurferdir/$j
	mri_easyreg --ref mri/T1.mgz --flo mti/mtr.nii.gz --ref_seg mri/synthseg.mgz --flo_seg mti/mtr_synthseg.mgz --bak_field mti/fs2mtr_field.mgz --threads 4 --affine_only
	#register aparc+aseg to mtr, using the inverse of the mtr to t1 registration
	mri_easywarp --i mri/aparc+aseg.mgz --o mti/aparc+aseg2mtr.mgz --field mti/fs2mtr_field.mgz --threads 4 --nearest
	mri_easywarp --i mri/aparc.a2009s+aseg.mgz --o mti/aparc.a2009s+aseg2mtr.mgz --field mti/fs2mtr_field.mgz --threads 4 --nearest
	mri_segstats --seg mti/aparc+aseg2mtr.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i mti/mtr.nii.gz --sum stats/aparc+aseg2mtr.stats --excludeid 0 30 31 62 63 72 77 78 79 80 81 82 85
fi
##MTR AD sig
if [ ! -f $freesurferdir/$j/stats/lh_AD_sig2mtr.stats ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	cd $freesurferdir/$j/
	side=( rh lh )
	for s in "${side[@]}"; do
		mri_binarize --i mti/aparc.a2009s+aseg2mtr.mgz --match-ctab $projectdir/${s}_AD_signature_LUT.txt --o mti/${s}_mtr_AD_sig.mgz
		mri_segstats --seg mti/${s}_mtr_AD_sig.mgz --mask mti/${s}_mtr_AD_sig.mgz --excludeid 0 --ctab $projectdir/${s}_AD_signature_LUT.txt --i mti/mtr.nii.gz --sum stats/${s}_AD_sig2mtr.stats
	done
fi
##MTR JHU
if [ ! -f $freesurferdir/$j/stats/jhu2mtr.stats ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	cd $freesurferdir/$j/mti
	mri_easyreg --ref $FREESURFER_HOME/subjects/cvs_avg35_inMNI152/mri/norm.mgz --flo mtr.nii.gz --ref_seg $FREESURFER_HOME/subjects/cvs_avg35_inMNI152/mri/aparc+aseg.mgz --flo_seg aparc+aseg2mtr.mgz --bak_field mni2mtr_field.mgz --threads 4
	mri_easywarp --i $FSLDIR/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz --o jhu2mtr.mgz --field mni2mtr_field.mgz --nearest --threads 4
	cd $freesurferdir/$j/
	mri_segstats --seg mti/jhu2mtr.mgz --mask mti/jhu2mtr.mgz --i mti/mtr.nii.gz --sum stats/jhu2mtr.stats --exludeid 0
fi

#G ratio (WM only)
if [ ! -f $freesurferdir/$j/stats/aparc+aseg2g_ratio.stats ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	#register MTR to diffusion
	cd $freesurferdir/$j
	mri_easyreg --ref mti/mtr.nii.gz --flo diffusion/B0.nii --ref_seg mti/mtr_synthseg.mgz --flo_seg diffusion/aparc+aseg2diff.mgz --bak_field diffusion/mti2diff_field.mgz --threads 4 --affine_only
	mri_easywarp --i mti/mtr.nii.gz --o diffusion/mtr2diff.nii --field diffusion/mti2diff_field.mgz --threads 4
	#calculate g-ratio in diffusion space (white matter only)
	singularity exec --nv --bind $basedir $basedir/neurodock_v1.0.0.sif mrcalc diffusion/mtr2diff.nii -neg diffusion/smi_matlab_f.nii -div 1 -add -sqrt diffusion/g_ratio.nii
	mri_segstats --seg diffusion/aparc+aseg2diff_NODDI_mask.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i diffusion/g_ratio.nii --sum stats/aparc+aseg2g_ratio.stats --id 2 7 41 46 251 252 253 254 255
fi
##JHU g-ratio
if [ ! -f $freesurferdir/$j/stats/jhu2g_ratio.stats ] && [ -f ${rawmti}_mti_bl.nii.gz ]; then
	cd $freesurferdir/$j
	mri_segstats --seg diffusion/jhu2diff_NODDI_mask.mgz --mask diffusion/jhu2diff_NODDI_mask.mgz --i diffusion/g_ratio.nii --sum stats/jhu2g_ratio.stats --excludeid 0
fi
