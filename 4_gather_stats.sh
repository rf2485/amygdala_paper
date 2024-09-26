basedir=/gpfs/data/lazarlab/CamCan995
projectdir=$basedir/derivatives/gm_roi/
t1dir=$projectdir/freesurfer/

module load freesurfer/7.4.1
export SUBJECTS_DIR=$t1dir

cut -f1 $projectdir/dwi_over_55.tsv > $t1dir/subjectsfile.txt
cd $t1dir
sed -i '' '1d' subjectsfile.txt

#generate stats tables with Freesurfer
aparcstats2table --subjectsfile=subjectsfile.txt --hemi lh --tablefile=lh_aparctable.tsv --measure=thickness --common-parcs --skip
aparcstats2table --subjectsfile=subjectsfile.txt --hemi rh --tablefile=rh_aparctable.tsv --measure=thickness --common-parcs --skip
asegstats2table --subjectsfile=subjectsfile.txt --tablefile=asegtable.tsv --common-segs --skip
asegstats2table --subjectsfile=subjectsfile.txt --stats=wmparc.stats --tablefile=wmparctable.tsv --common-segs --skip

meas_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd fit_FWF fit_NDI fit_ODI )
for meas in "${meas_list[@]}"; do
  asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=aparc+aseg2${meas}.stats --tablefile=aparc+aseg2${meas}.tsv --common-segs
done