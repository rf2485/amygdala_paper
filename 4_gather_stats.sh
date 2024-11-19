# basedir=/gpfs/data/lazarlab/CamCan995
basedir=/Volumes/Research/lazarm03lab/labspace/AD/camcan995/
projectdir=$basedir/derivatives/scd/gm_roi/
t1dir=$projectdir/freesurfer/

# module load freesurfer/7.4.1
export SUBJECTS_DIR=$t1dir

cut -f1 $projectdir/dwi_over_55.tsv > $t1dir/subjectsfile.txt
cd $t1dir
sed -i '' '1d' subjectsfile.txt

#generate stats tables with Freesurfer
side=( rh lh )
for s in "${side[@]}"; do
	aparcstats2table --subjectsfile=subjectsfile.txt --hemi ${s} --tablefile=${s}_aparctable_thickness.tsv --measure=thickness --common-parcs --skip
	aparcstats2table --subjectsfile=subjectsfile.txt --hemi ${s} --tablefile=${s}_aparctable_volume.tsv --measure=volume --common-parcs --skip
	asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=$s.AD_sig_thickness.stats --tablefile=${s}_AD_sig_thickness.tsv
	asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=$s.AD_sig_volume.stats --tablefile=${s}_AD_sig_volume.tsv
done
asegstats2table --subjectsfile=subjectsfile.txt --tablefile=asegtable.tsv --common-segs --skip
asegstats2table --subjectsfile=subjectsfile.txt --stats=wmparc.stats --tablefile=wmparctable.tsv --common-segs --skip
asegstats2table --subjectsfile=subjectsfile.txt --stats=aparc+aseg2dki_ak.stats --tablefile=aparc+aseg_volume.tsv --common-segs
meas_list=( dki_ak dki_kfa dki_mk dki_mkt dki_rk dti_ad dti_fa dti_md dti_rd fit_FWF fit_NDI fit_ODI )
for meas in "${meas_list[@]}"; do
  asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=aparc+aseg2${meas}.stats --tablefile=aparc+aseg2${meas}.tsv --common-segs
  side=( rh lh )
  for s in "${side[@]}"; do
	  asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=${s}_AD_sig2${meas}.stats --tablefile=${s}_AD_sig2${meas}.tsv --common-segs
  done
done