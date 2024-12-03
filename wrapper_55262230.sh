#!/usr/bin/bash

#SBATCH --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,MKL_DOMAIN_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,GOTO_NUM_THREADS=1,FSLSUB_PARALLEL=1,FSLSUB_JOB_ID_VAR=SLURM_JOB_ID,FSLSUB_ARRAYTASKID_VAR=SLURM_ARRAY_TASK_ID,FSLSUB_ARRAYSTARTID_VAR=SLURM_ARRAY_TASK_MIN,FSLSUB_ARRAYENDID_VAR=SLURM_ARRAY_TASK_MAX,FSLSUB_ARRAYSTEPSIZE_VAR=SLURM_ARRAY_TASK_STEP,FSLSUB_ARRAYCOUNT_VAR=SLURM_ARRAY_TASK_COUNT,FSLSUB_NSLOTS=SLURM_NPROCS
#SBATCH -o tbss_logs/dwi_randomise_array.o%A.%a
#SBATCH -e tbss_logs/dwi_randomise_array.e%A.%a
#SBATCH --dependency=afterany:55262168
#SBATCH --mem=64G
#SBATCH -t 239
#SBATCH --job-name=dwi_randomise_array
#SBATCH --chdir=/gpfs/data/lazarlab/CamCan995/derivatives/scd/main
#SBATCH -p cpu_dev
#SBATCH --parsable
#SBATCH --requeue
#SBATCH --array=1-13
#SBATCH --ntasks=1
# Built by fsl_sub v.2.8.4 and fsl_sub_plugin_slurm v.1.6.5
# Command line: /gpfs/home/rf2485/.conda/envs/fsl_eddy/bin/fsl_sub -T 239 -R 64 -j 55262168 -l tbss_logs -t /gpfs/data/lazarlab/CamCan995//derivatives/scd/main//dwi_randomise_array
# Submission time (H:M:S DD/MM/YYYY): 17:52:57 26/11/2024


the_command=$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" /gpfs/data/lazarlab/CamCan995//derivatives/scd/main//dwi_randomise_array)

exec /usr/bin/bash -c "$the_command"

