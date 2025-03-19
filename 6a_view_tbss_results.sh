metric_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd smi_matlab_f smi_matlab_Da smi_matlab_DePar smi_matlab_DePerp smi_matlab_p2 wm_fit_FWF wm_fit_NDI wm_fit_ODI)
# basedir=/gpfs/data/lazarlab/CamCan995/
basedir=/Volumes/Research/lazarm03lab/labspace/AD/camcan995/
projectdir=$basedir/derivatives/scd/main/

# module load miniconda3/gpu/4.9.2
# source activate ~/.conda/envs/fsl_eddy/
# export FSLDIR=$CONDA_PREFIX
# source $FSLDIR/etc/fslconf/fsl.sh

cd $projectdir/tbss/stats
for metric in "${metric_list[@]}"; do
	fsleyes -std1mm all_${metric} mean_FA_skeleton_mask -cm green ${metric}_clusterm_corrp_tstat1 -cm blue-lightblue -dr 0.949 1 ${metric}_clusterm_corrp_tstat2 -cm red -dr 0.949 1 #red means higher in SCD (tstat2)
done

#interaction tests only done on metrics with group differences
int_test_list=(memory_int anxiety_int depression_int age_int bmi_int)
group_diffs_list=( dti_fa dti_md dti_rd smi_matlab_DePerp smi_matlab_p2 )
for test in "${int_test_list[@]}"; do
	for metric in "${group_diffs_list[@]}"; do
		fsleyes -std1mm all_${metric} mean_FA_skeleton_mask -cm green ${test}/${test}_${metric}_clusterm_corrp_tstat1 -cm blue-lightblue -dr 0.949 1 ${test}/${test}_${metric}_clusterm_corrp_tstat2 -cm red -dr 0.949 1 #red means SCD has a steeper slope than controls (tstat2)
	done
done

#across group tests only done on metrics without interaction effects
