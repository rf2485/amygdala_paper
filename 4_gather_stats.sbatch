#!/bin/bash
#SBATCH --mail-user=rf2485@nyulangone.org
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH -o ./slurm_output/4_gather_stats/slurm-%A_%a.out

basedir=/gpfs/data/lazarlab/CamCan995
# basedir=/Volumes/Research/lazarm03lab/labspace/AD/camcan995/
projectdir=$basedir/derivatives/scd/main/
t1dir=$projectdir/freesurfer/

module load freesurfer/7.4.1
export SUBJECTS_DIR=$t1dir

cut -f1 $projectdir/dwi_mti_over_55.tsv > $t1dir/subjectsfile.txt
cd $t1dir
sed -i '' '1d' subjectsfile.txt

#generate stats tables with Freesurfer
side=( rh lh )
for s in "${side[@]}"; do
	aparcstats2table --subjectsfile=subjectsfile.txt --hemi ${s} --tablefile=${s}_aparctable_thickness.tsv --measure=thickness --skip
done

asegstats2table --subjectsfile=subjectsfile.txt --tablefile=asegtable.tsv --common-segs --skip --no-segno 0 30 31 62 63 72 77 78 79 80 81 82 85

meas_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd mtr )
for meas in "${meas_list[@]}"; do
  asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=aparc+aseg2${meas}.stats --tablefile=aparc+aseg2${meas}.tsv --common-segs --skip --no-segno 0 30 31 62 63 72 77 78 79 80 81 82 85
done

gm_meas_list=( fit_FWF fit_NDI fit_ODI )
for meas in "${gm_meas_list[@]}"; do
	asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=aparc+aseg2${meas}.stats --tablefile=aparc+aseg2${meas}.tsv --common-segs --skip --no-segno 0 30 31 62 63 72 77 78 79 80 81 82 85 2 7 41 46 251 252 253 254 255
done
	
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_wm.stats --skip --tablefile=mtr_wm.tsv --common-segs
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_gm.stats --skip --tablefile=mtr_gm.tsv --common-segs

