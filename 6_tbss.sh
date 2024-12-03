basedir=/gpfs/data/lazarlab/CamCan995/
projectdir=$basedir/derivatives/scd/main/

module load miniconda3/gpu/4.9.2
source activate ~/.conda/envs/fsl_eddy/
export FSLDIR=$CONDA_PREFIX
source $FSLDIR/etc/fslconf/fsl.sh

problem_subjs=( scd_sub-CC510255 scd_sub-CC620821 ctl_sub-CC621011 ctl_sub-CC721292 )
# 510255 has pathology in the L temporal lobe causing errors in registration to template
# 620821, 621011, and 721292 have big ventricles, causing errors in registration to template

#TBSS for non-FA diffusion
mkdir -p $projectdir/tbss/stats
mkdir $projectdir/tbss/MD
mkdir $projectdir/tbss/RD
mkdir $projectdir/tbss/AxD
mkdir $projectdir/tbss/KFA
mkdir $projectdir/tbss/MK
mkdir $projectdir/tbss/RK
mkdir $projectdir/tbss/AK
mkdir $projectdir/tbss/Da_smi
mkdir $projectdir/tbss/DePar_smi
mkdir $projectdir/tbss/DePerp_smi
mkdir $projectdir/tbss/f_smi
mkdir $projectdir/tbss/p2_smi
 
for j in $(cut -f1 $projectdir/dwi_over_55_ctl.tsv); do
	echo "copying files for ${j}"
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_fa.nii -nan $projectdir/tbss/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_md.nii -nan $projectdir/tbss/MD/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_rd.nii -nan $projectdir/tbss/RD/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_ad.nii -nan $projectdir/tbss/AxD/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_kfa.nii -nan $projectdir/tbss/KFA/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_mk.nii -nan $projectdir/tbss/MK/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_rk.nii -nan $projectdir/tbss/RK/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_ak.nii -nan $projectdir/tbss/AK/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_Da.nii -nan $projectdir/tbss/Da_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics//smi_matlab_DePar.nii -nan $projectdir/tbss/DePar_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_DePerp.nii -nan $projectdir/tbss/DePerp_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_f.nii -nan $projectdir/tbss/f_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_p2.nii -nan $projectdir/tbss/p2_smi/ctl_${j}.nii.gz
done

for j in $(cut -f1 $projectdir/dwi_over_55_scd.tsv); do
	echo "copying files for ${j}"
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_fa.nii -nan $projectdir/tbss/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_md.nii -nan $projectdir/tbss/MD/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_rd.nii -nan $projectdir/tbss/RD/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_ad.nii -nan $projectdir/tbss/AxD/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_kfa.nii -nan $projectdir/tbss/KFA/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_mk.nii -nan $projectdir/tbss/MK/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_rk.nii -nan $projectdir/tbss/RK/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_ak.nii -nan $projectdir/tbss/AK/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_Da.nii -nan $projectdir/tbss/Da_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics//smi_matlab_DePar.nii -nan $projectdir/tbss/DePar_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_DePerp.nii -nan $projectdir/tbss/DePerp_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_f.nii -nan $projectdir/tbss/f_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_p2.nii -nan $projectdir/tbss/p2_smi/scd_${j}.nii.gz
done

for subj in "${problem_subjs[@]}"; do
	echo "removing ${subj}"
	rm $projectdir/tbss/${subj}.nii.gz
	rm $projectdir/tbss/*/${subj}.nii.gz
done

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
cd $projectdir/tbss/stats
design_ttest2 design 196 125
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_randomise_array

#TBSS for MTR TR=30ms
mkdir -p $projectdir/tbss_tr30/stats
mkdir $projectdir/tbss_tr30/MD
mkdir $projectdir/tbss_tr30/RD
mkdir $projectdir/tbss_tr30/AxD
mkdir $projectdir/tbss_tr30/KFA
mkdir $projectdir/tbss_tr30/MK
mkdir $projectdir/tbss_tr30/RK
mkdir $projectdir/tbss_tr30/AK
mkdir $projectdir/tbss_tr30/Da_smi
mkdir $projectdir/tbss_tr30/DePar_smi
mkdir $projectdir/tbss_tr30/DePerp_smi
mkdir $projectdir/tbss_tr30/f_smi
mkdir $projectdir/tbss_tr30/p2_smi
mkdir $projectdir/tbss_tr30/mtr
mkdir $projectdir/tbss_tr30/g_ratio
 
for j in $(cut -f1 $projectdir/mti_over_55_tr30_ctl.tsv); do
	echo "copying files for ${j}"
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_fa.nii -nan $projectdir/tbss_tr30/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_md.nii -nan $projectdir/tbss_tr30/MD/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_rd.nii -nan $projectdir/tbss_tr30/RD/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_ad.nii -nan $projectdir/tbss_tr30/AxD/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_kfa.nii -nan $projectdir/tbss_tr30/KFA/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_mk.nii -nan $projectdir/tbss_tr30/MK/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_rk.nii -nan $projectdir/tbss_tr30/RK/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_ak.nii -nan $projectdir/tbss_tr30/AK/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_Da.nii -nan $projectdir/tbss_tr30/Da_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics//smi_matlab_DePar.nii -nan $projectdir/tbss_tr30/DePar_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_DePerp.nii -nan $projectdir/tbss_tr30/DePerp_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_f.nii -nan $projectdir/tbss_tr30/f_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_p2.nii -nan $projectdir/tbss_tr30/p2_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/freesurfer/$j/diffusion/mtr2diff.nii -nan $projectdir/tbss_tr30/mtr/ctl_${j}.nii.gz
	fslmaths $projectdir/freesurfer/$j/diffusion/g_ratio.nii -nan $projectdir/tbss_tr30/g_ratio/ctl_${j}.nii.gz
done

for j in $(cut -f1 $projectdir/mti_over_55_tr30_scd.tsv); do
	echo "copying files for ${j}"
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_fa.nii -nan $projectdir/tbss_tr30/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_md.nii -nan $projectdir/tbss_tr30/MD/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_rd.nii -nan $projectdir/tbss_tr30/RD/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_ad.nii -nan $projectdir/tbss_tr30/AxD/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_kfa.nii -nan $projectdir/tbss_tr30/KFA/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_mk.nii -nan $projectdir/tbss_tr30/MK/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_rk.nii -nan $projectdir/tbss_tr30/RK/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_ak.nii -nan $projectdir/tbss_tr30/AK/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_Da.nii -nan $projectdir/tbss_tr30/Da_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics//smi_matlab_DePar.nii -nan $projectdir/tbss_tr30/DePar_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_DePerp.nii -nan $projectdir/tbss_tr30/DePerp_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_f.nii -nan $projectdir/tbss_tr30/f_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_p2.nii -nan $projectdir/tbss_tr30/p2_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_fa.nii -nan $projectdir/tbss_tr30/scd_${j}.nii.gz
	fslmaths $projectdir/freesurfer/$j/diffusion/mtr2diff.nii -nan $projectdir/tbss_tr30/mtr/scd_${j}.nii.gz
	fslmaths $projectdir/freesurfer/$j/diffusion/g_ratio.nii -nan $projectdir/tbss_tr30/g_ratio/scd_${j}.nii.gz
done

for subj in "${problem_subjs[@]}"; do
	echo "removing ${subj}"
	rm $projectdir/tbss_tr30/${subj}.nii.gz
	rm $projectdir/tbss_tr30/*/${subj}.nii.gz
done

cd $projectdir/tbss_tr30
tbss_1_preproc *.nii.gz
tbss_2_reg -T
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
fsl_sub -T 239 -R 32 -j $most_recent_job -l tbss_logs tbss_3_postreg -S
# fsleyes $tbssdir/stats/all_FA -dr 0 1 $tbssdir/stats/mean_FA_skeleton -dr 0.2 1 -cm green
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
fsl_sub -T 239 -R 32 -j $most_recent_job -l tbss_logs tbss_4_prestats 0.3
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/mtr_non_FA_array
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_non_FA_array
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
cd $projectdir/tbss_tr30/stats
design_ttest2 design 102 75
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/mtr_randomise_array
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_randomise_array


#TBSS for MTR TR=50ms
mkdir -p $projectdir/tbss_tr50/stats
mkdir $projectdir/tbss_tr50/MD
mkdir $projectdir/tbss_tr50/RD
mkdir $projectdir/tbss_tr50/AxD
mkdir $projectdir/tbss_tr50/KFA
mkdir $projectdir/tbss_tr50/MK
mkdir $projectdir/tbss_tr50/RK
mkdir $projectdir/tbss_tr50/AK
mkdir $projectdir/tbss_tr50/Da_smi
mkdir $projectdir/tbss_tr50/DePar_smi
mkdir $projectdir/tbss_tr50/DePerp_smi
mkdir $projectdir/tbss_tr50/f_smi
mkdir $projectdir/tbss_tr50/p2_smi
mkdir $projectdir/tbss_tr50/mtr
mkdir $projectdir/tbss_tr50/g_ratio
 
for j in $(cut -f1 $projectdir/mti_over_55_tr50_ctl.tsv); do
	echo "copying files for ${j}"
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_fa.nii -nan $projectdir/tbss_tr50/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_md.nii -nan $projectdir/tbss_tr50/MD/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_rd.nii -nan $projectdir/tbss_tr50/RD/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_ad.nii -nan $projectdir/tbss_tr50/AxD/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_kfa.nii -nan $projectdir/tbss_tr50/KFA/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_mk.nii -nan $projectdir/tbss_tr50/MK/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_rk.nii -nan $projectdir/tbss_tr50/RK/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_ak.nii -nan $projectdir/tbss_tr50/AK/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_Da.nii -nan $projectdir/tbss_tr50/Da_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics//smi_matlab_DePar.nii -nan $projectdir/tbss_tr50/DePar_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_DePerp.nii -nan $projectdir/tbss_tr50/DePerp_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_f.nii -nan $projectdir/tbss_tr50/f_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_p2.nii -nan $projectdir/tbss_tr50/p2_smi/ctl_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_fa.nii -nan $projectdir/tbss_tr50/ctl_${j}.nii.gz
	fslmaths $projectdir/freesurfer/$j/diffusion/mtr2diff.nii -nan $projectdir/tbss_tr50/mtr/ctl_${j}.nii.gz
	fslmaths $projectdir/freesurfer/$j/diffusion/g_ratio.nii -nan $projectdir/tbss_tr50/g_ratio/ctl_${j}.nii.gz
done

for j in $(cut -f1 $projectdir/mti_over_55_tr50_scd.tsv); do
	echo "copying files for ${j}"
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_fa.nii -nan $projectdir/tbss_tr50/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_md.nii -nan $projectdir/tbss_tr50/MD/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_rd.nii -nan $projectdir/tbss_tr50/RD/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_ad.nii -nan $projectdir/tbss_tr50/AxD/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_kfa.nii -nan $projectdir/tbss_tr50/KFA/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_mk.nii -nan $projectdir/tbss_tr50/MK/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_rk.nii -nan $projectdir/tbss_tr50/RK/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dki_ak.nii -nan $projectdir/tbss_tr50/AK/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_Da.nii -nan $projectdir/tbss_tr50/Da_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics//smi_matlab_DePar.nii -nan $projectdir/tbss_tr50/DePar_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_DePerp.nii -nan $projectdir/tbss_tr50/DePerp_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_f.nii -nan $projectdir/tbss_tr50/f_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/smi_matlab_p2.nii -nan $projectdir/tbss_tr50/p2_smi/scd_${j}.nii.gz
	fslmaths $projectdir/dwi_processed/$j/metrics/dti_fa.nii -nan $projectdir/tbss_tr50/scd_${j}.nii.gz
	fslmaths $projectdir/freesurfer/$j/diffusion/mtr2diff.nii -nan $projectdir/tbss_tr50/mtr/scd_${j}.nii.gz
	fslmaths $projectdir/freesurfer/$j/diffusion/g_ratio.nii -nan $projectdir/tbss_tr50/g_ratio/scd_${j}.nii.gz
done

for subj in "${problem_subjs[@]}"; do
	echo "removing ${subj}"
	rm $projectdir/tbss_tr50/${subj}.nii.gz
	rm $projectdir/tbss_tr50/*/${subj}.nii.gz
done

cd $projectdir/tbss_tr50
tbss_1_preproc *.nii.gz
tbss_2_reg -T
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
fsl_sub -T 239 -R 32 -j $most_recent_job -l tbss_logs tbss_3_postreg -S
# fsleyes $tbssdir/stats/all_FA -dr 0 1 $tbssdir/stats/mean_FA_skeleton -dr 0.2 1 -cm green
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
fsl_sub -T 239 -R 32 -j $most_recent_job -l tbss_logs tbss_4_prestats 0.3
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/mtr_non_FA_array
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_non_FA_array
most_recent_job=$(squeue -u rf2485 --nohead --format %F | head -n 1)
cd $projectdir/tbss_tr50/stats
design_ttest2 design 82 42
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/mtr_randomise_array
fsl_sub -T 239 -R 64 -j $most_recent_job -l tbss_logs -t $projectdir/dwi_randomise_array
