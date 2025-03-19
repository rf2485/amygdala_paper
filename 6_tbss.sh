basedir=/gpfs/data/lazarlab/CamCan995/
projectdir=$basedir/derivatives/scd/camcan_scd/

module load miniconda3/gpu/4.9.2
source activate ~/.conda/envs/fsl_eddy/
export FSLDIR=$CONDA_PREFIX
source $FSLDIR/etc/fslconf/fsl.sh
module load r/4.3.2

problem_subjs=( scd_sub-CC510255 scd_sub-CC620821 ctl_sub-CC510438 ctl_sub-CC621011 ctl_sub-CC710551 ctl_sub-CC721292 )
# 510255 has pathology in the L temporal lobe causing errors in registration to template
# sub-CC510438 has pathology in L frontal lobe
# sub-CC710551 has motion artifacts in DWI
# 620821, 621011, and 721292 have big ventricles, causing errors in registration to template

#create stats contrasts and matrices
mkdir -p $projectdir/tbss/stats
Rscript fsl_glm_matrices.R

#copy FA files to TBSS folder
for j in $(cut -f1 $projectdir/dwi_over_55_ctl.tsv); do
	echo "copying FA for ${j}"
	fslmaths $projectdir/freesurfer/$j/diffusion/dti_fa.nii.gz -nan $projectdir/tbss/ctl_${j}.nii.gz
done
for j in $(cut -f1 $projectdir/dwi_over_55_scd.tsv); do
	echo "copying FA for ${j}"
	fslmaths $projectdir/freesurfer/$j/diffusion/dti_fa.nii.gz -nan $projectdir/tbss/scd_${j}.nii.gz
done

#copy non-FA diffusion files to TBSS folder
meas_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd smi_matlab_f smi_matlab_Da smi_matlab_DePar smi_matlab_DePerp smi_matlab_p2 wm_fit_FWF wm_fit_NDI wm_fit_ODI )
for meas in "${meas_list[@]}"; do
	mkdir -p $projectdir/tbss/$meas
	for j in $(cut -f1 $projectdir/dwi_over_55_ctl.tsv); do
		echo "copying ${meas} files for ${j}"
		fslmaths $projectdir/freesurfer/$j/diffusion/${meas}_masked_wm.nii.gz -nan $projectdir/tbss/$meas/ctl_${j}.nii.gz
	done
	for j in $(cut -f1 $projectdir/dwi_over_55_scd.tsv); do
		echo "copying ${meas} files for ${j}"
		fslmaths $projectdir/freesurfer/$j/diffusion/${meas}_masked_wm.nii.gz -nan $projectdir/tbss/$meas/scd_${j}.nii.gz
	done
done

#remove problem subjects
for subj in "${problem_subjs[@]}"; do
	echo "removing ${subj}"
	rm $projectdir/tbss/${subj}.nii.gz
	rm $projectdir/tbss/*/${subj}.nii.gz
done

#complete TBSS pipeline
cd $projectdir/tbss
tbss_1_preproc *.nii.gz
tbss_2_reg -T
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
fsl_sub -T 239 -R 32 -j $most_recent_job -l tbss_logs tbss_3_postreg -S
# fsleyes $tbssdir/stats/all_FA -dr 0 1 $tbssdir/stats/mean_FA_skeleton -dr 0.2 1 -cm green
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
fsl_sub -T 239 -R 32 -j $most_recent_job -l tbss_logs tbss_4_prestats 0.3
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_non_FA_array
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)

#stats
cd $projectdir/tbss/stats
design_ttest2 design 194 125
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_randomise_array
int_test_list=(memory_int anxiety_int depression_int age_int bmi_int memory anxiety depression age bmi)
for test in "${int_test_list[@]}"; do
	Text2Vest ${test}_mat.txt ${test}.mat
	Text2Vest ${test}_con.txt ${test}.con
	mkdir ${test}
done
fsl_sub  -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_randomise_int_array
fsl_sub  -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_randomise_across_groups_array

#TBSS for each MTR TR
# mkdir -p $projectdir/tbss_30/stats
# cd $projectdir/tbss_30/stats
# design_ttest2 design 101 75
#
# mkdir -p $projectdir/tbss_50/stats
# cd $projectdir/tbss_30/stats
# design_ttest2 design 81 43
#
# #copy FA files to TBSS folder
# tr_list=( tr30 tr50 )
# for tr in "${tr_list[@]}"; do
# 	mkdir -p $projectdir/tbss_${tr}/stats
# 	for j in $(cut -f1 $projectdir/mti_over_55_${tr}_ctl.tsv); do
# 		echo "copying FA for ${j}"
# 		fslmaths $projectdir/freesurfer/$j/diffusion/dti_fa_masked_wm.nii.gz -nan $projectdir/tbss_${tr}/ctl_${j}.nii.gz
# 	done
# 	for j in $(cut -f1 $projectdir/mti_over_55_${tr}_scd.tsv); do
# 		echo "copying FA for ${j}"
# 		fslmaths $projectdir/freesurfer/$j/diffusion/dti_fa_masked_wm.nii.gz -nan $projectdir/tbss_${tr}/scd_${j}.nii.gz
# 	done
#
# 	#copy non-FA diffusion files to TBSS folder
# 	non_fa_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_md dti_rd smi_matlab_f smi_matlab_Da smi_matlab_DePar smi_matlab_DePerp smi_matlab_p2 wm_fit_FWF wm_fit_NDI wm_fit_ODI mtr2diff g_ratio )
# 	for meas in "${non_fa_list[@]}"; do
# 		mkdir -p $projectdir/tbss_${tr}/$meas
# 		for j in $(cut -f1 $projectdir/mti_over_55_${tr}_ctl.tsv); do
# 			echo "copying ${meas} files for ${j}"
# 			fslmaths $projectdir/freesurfer/$j/diffusion/${meas}_masked_wm.nii.gz -nan $projectdir/tbss_${tr}/$meas/ctl_${j}.nii.gz
# 		done
# 		for j in $(cut -f1 $projectdir/mti_over_55_${tr}_scd.tsv); do
# 			echo "copying ${meas} files for ${j}"
# 			fslmaths $projectdir/freesurfer/$j/diffusion/${meas}_masked_wm.nii.gz -nan $projectdir/tbss_${tr}/$meas/scd_${j}.nii.gz
# 		done
# 	done
#
# 	for subj in "${problem_subjs[@]}"; do
# 		echo "removing ${subj}"
# 		rm $projectdir/tbss_${tr}/${subj}.nii.gz
# 		rm $projectdir/tbss_${tr}/*/${subj}.nii.gz
# 	done
#
# 	cd $projectdir/tbss_${tr}
# 	tbss_1_preproc *.nii.gz
# 	tbss_2_reg -T
# 	most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
# 	fsl_sub -T 239 -R 32 -j $most_recent_job -l tbss_logs tbss_3_postreg -S
# 	# fsleyes $tbssdir/stats/all_FA -dr 0 1 $tbssdir/stats/mean_FA_skeleton -dr 0.2 1 -cm green
# 	most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
# 	fsl_sub -T 239 -R 32 -j $most_recent_job -l tbss_logs tbss_4_prestats 0.3
# 	most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
# 	fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/mtr_non_FA_array
# 	fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_non_FA_array
# 	most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
# 	cd $projectdir/tbss_${tr}/stats
# 	fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/mtr_randomise_array
# 	fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_randomise_array
# done
