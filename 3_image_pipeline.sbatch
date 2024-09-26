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
projectdir=$basedir/derivatives/gm_roi/

subj_list=$(cut -f1 -d$'\t' $projectdir/dwi_over_55.tsv)
subj_list=($subj_list)
subj_num=$(($SLURM_ARRAY_TASK_ID))
j=${subj_list[$subj_num]}

rawdwi=$basedir/raw/$j/dwi/${j}_dwi
rawanat=$basedir/raw/$j/anat/${j}
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

#recon-all-clinical to generate segmentations
recon-all-clinical.sh ${rawanat}_T1w.nii.gz $j 4 $freesurferdir
##create files necessary for including T2 in the recon
mkdir -p $freesurferdir/$j/mri/orig
mri_convert ${rawanat}_T1w.nii.gz $freesurferdir/$j/mri/orig/001.mgz
mv $freesurferdir/$j/surf/lh.pial $freesurferdir/$j/surf/lh.pial.T1
mv $freesurferdir/$j/surf/rh.pial $freesurferdir/$j/surf/rh.pial.T1
##include T2 in the recon
recon-all -T2 ${rawanat}_T2w.nii.gz -motioncor -talairach -nuintensitycor -normalization -gcareg -maskbfs -T2pial -cortribbon -parcstats -cortparc2 -parcstats2 -cortparc3 -parcstats3 -pctsurfcon -hyporelabel -aparc2aseg -apas2aseg -segstats -wmparc -balabels -subjid ${j} -sd ${freesurferdir} -threads 4

#concatenate and convert raw dwi
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert -json_import $rawdwi.json -fslgrad $rawdwi.bvec $rawdwi.bval $rawdwi.nii.gz $designerdir/working.mif -force
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert $designerdir/working.mif -json_export $designerdir/dwi_raw.json -export_grad_fsl $designerdir/dwi_raw.bvec $designerdir/dwi_raw.bval $designerdir/dwi_raw.nii -force

#denoise, degibbs, and normalize with designer
rm -rf $designerdir/tmp
singularity exec --nv --bind $basedir $basedir/designer2_latest.sif designer -denoise -shrinkage frob -adaptive_patch -degibbs -nthreads 4 -nocleanup -scratch $designerdir/tmp $designerdir/dwi_raw.nii $intdir/2_dwi_degibbs.nii

#move denoised files and calculate residual
mv $designerdir/tmp/sigma.nii $designerdir/noisemap.nii
mv $designerdir/tmp/tmp_dwidn.nii $intdir/1_dwi_denoised.nii
cp $designerdir/dwi_raw.json $intdir/1_dwi_denoised.json
cp $designerdir/dwi_raw.bvec $intdir/1_dwi_denoised.bvec
cp $designerdir/dwi_raw.bval $intdir/1_dwi_denoised.bval
singularity exec --bind $basedir $basedir/neurodock_latest.sif mrcalc $designerdir/dwi_raw.nii $designerdir/intermediate_nifti/1_dwi_denoised.nii -sub $designerdir/residual.nii

#move and convert degibbsed files
cp $designerdir/dwi_raw.json $intdir/2_dwi_degibbs.json
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert -json_import $intdir/2_dwi_degibbs.json -fslgrad $intdir/2_dwi_degibbs.bvec $intdir/2_dwi_degibbs.bval $intdir/2_dwi_degibbs.nii $designerdir/working.mif -force

#prepare for synb0
##create b0 and acqparams
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif dwiextract -bzero $designerdir/working.mif - | singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert - -coord 3 0 -axes 0,1,2 $intdir/b0_1_working.mif -force
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert -export_pe_table $synb0in/acqparams.txt $intdir/b0_1_working.mif $synb0in/b0.nii.gz -force
rm $intdir/b0_1_working.mif
echo " 0 1 0 0.000" >> $synb0in/acqparams.txt
##copy T1
cp ${rawanat}_T1w.nii.gz $synb0in/T1.nii.gz


#synb0
singularity run -e -B $synb0in:/INPUTS -B $synb0out:/OUTPUTS -B $basedir/license.txt:/extra/freesurfer/license.txt $basedir/synb0-disco_v3.1.sif
##copy and convert synb0 output
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert -import_pe_table $synb0in/acqparams.txt $synb0out/b0_all.nii.gz $synb0out/b0_all.mif -force

#copy and convert files for eddy_correct
cp $synb0out/b0_all.mif $designerdir/tmp/se_epi.mif
cd $designerdir/tmp
cp $designerdir/working.mif dwi.mif
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert dwi.mif dwi.nii -strides -1,+2,+3,+4 -json_export dwi.json -export_grad_fsl dwi.bvec dwi.bval

#eddy_correct, using first volume as reference, with no blurring
eddy_correct dwi dwi_post_eddy 0 -trilinear -noresampleblur
##rotate bvecs using eddy_correct transformation log
fdt_rotate_bvecs dwi.bvec dwi_post_eddy.bvec dwi_post_eddy.ecclog
##convert eddy_correct output
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert dwi_post_eddy.nii.gz dwi_post_eddy.mif -json_import dwi.json -fslgrad dwi_post_eddy.bvec dwi.bval

#prepare for topup
##replace first b0 of se_epi with first b0 of dwi
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert dwi_post_eddy.mif dwi_first_bzero.mif -coord 3 0 -axes 0,1,2
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert se_epi.mif - -coord 3 1 | singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrcat dwi_first_bzero.mif - se_epi_firstdwibzero.mif -axis 3
##create files for topup input
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert se_epi_firstdwibzero.mif topup_in.nii -strides -1,+2,+3,+4 -export_pe_table topup_datain.txt

#topup
topup --imain=topup_in.nii --datain=topup_datain.txt --out=field --fout=field_map.nii.gz --config=${FSLDIR}/etc/flirtsch/b02b0.cnf --verbose
##create files for applytopup
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrinfo dwi_post_eddy.mif -export_pe_eddy applytopup_config.txt applytopup_indices.txt
##applytopup
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif applytopup --imain=dwi_post_eddy.nii.gz --datain=applytopup_config.txt --inindex=1 --topup=field --out=dwi_applytopup --method=jac
##convert applytopup output
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert dwi_applytopup.nii.gz dwi_applytopup.mif -json_import dwi.json -fslgrad dwi_post_eddy.bvec dwi.bval
##move and convert applytopup output
intermediate=$designerdir/intermediate_nifti/2_dwi_undistorted
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif mrconvert -export_grad_fsl $intermediate.bvec $intermediate.bval -json_export $intermediate.json dwi_applytopup.mif $intermediate.nii -force
mv dwi_applytopup.mif $designerdir/working.mif

#brain mask
singularity exec --bind $basedir/ $basedir/neurodock_latest.sif dwiextract $designerdir/working.mif - -bzero | singularity exec --bind $basedir/ $basedir/neurodock_latest.sif mrmath - mean - -axis 3 | singularity exec --bind $basedir/ $basedir/neurodock_latest.sif mrconvert - meanb0.nii
mri_synthstrip -i meanb0.nii -m new_brain_mask.nii

#csf mask, csf-excluded smoothing, rician bias correction, and model fitting with pydesigner
cd $designerdir
mkdir -p $designerdir/metrics_qc/fitting
singularity exec --bind $basedir/ $basedir/neurodock_latest.sif pydesigner --resume -n -u -cd 2 -z -r --fit_constraints 0,0,0 --akc --user_mask $designerdir/tmp/new_brain_mask.nii --verbose --noqc --nthreads 4 -o $designerdir/ $rawdwi.nii.gz
rm -rf $designerdir/tmp/

##calculate colormap from diffusion tensor
singularity exec --nv --bind $basedir $basedir/neurodock_latest.sif tensor2metric $designerdir/metrics/DT.nii -vector $designerdir/metrics/dti_V1.nii -modulate none -force

#fit NODDI model using the AMICO program
conda deactivate
source activate ~/.conda/envs/amico/
python3 $projectdir/amico_noddi.py $j

#register aparc to diffusion
cd $freesurferdir/$j/
mkdir -p diffusion
## copy B0 to freesurferdir
cp $designerdir/B0.nii diffusion/
## register B0 to freesurfer output
bbregister --s ${j} --mov diffusion/B0.nii --reg diffusion/b02fs.lta --dti --init-fsl
## use inverse of b02fs registration to register aparc+aseg to diffusion
mri_vol2vol --mov diffusion/B0.nii --targ mri/aparc+aseg.mgz --inv --interp nearest --o diffusion/aparc+aseg2diff.mgz --reg diffusion/b02fs.lta --no-save-reg
## apply registered aparc+aseg to diffusion metrics to generate stats
pyd_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd )
for meas in "${pyd_list[@]}"; do
	## copy metric to freesurferdir
	cp $designerdir/metrics/${meas}.nii diffusion/
	## generate stats
	mri_segstats --seg diffusion/aparc+aseg2diff.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i diffusion/${meas}.nii --sum stats/aparc+aseg2${meas}.stats
done
## repeat for NODDI
noddi_list=( fit_FWF fit_NDI fit_ODI )
for meas in "${noddi_list[@]}"; do
	cp $designerdir/AMICO/NODDI/${meas}.nii.gz diffusion/
	mri_segstats --seg diffusion/aparc+aseg2diff.mgz --ctab $FREESURFER_HOME/FreeSurferColorLUT.txt --i diffusion/${meas}.nii.gz --sum stats/aparc+aseg2${meas}.stats
done
