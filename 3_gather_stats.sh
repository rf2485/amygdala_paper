basedir=/gpfs/data/lazarlab/CamCan995
projectdir=$basedir/derivatives/mti_whole_wm_gm
t1dir=$projectdir/freesurfer/

module load freesurfer/7.4.1
export SUBJECTS_DIR=$t1dir

cut -f1 $projectdir/mti_over_55.tsv > $t1dir/subjectsfile.txt
cd $t1dir
sed -i '1d' subjectsfile.txt

#generate stats tables with Freesurfer
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_aparc+aseg.stats --tablefile=mtr_aparc+aseg.tsv --common-segs
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_wm.stats --tablefile=mtr_wm.tsv --common-segs
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_gm.stats --tablefile=mtr_gm.tsv --common-segs
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_ctx_wm.stats --tablefile=mtr_ctx_wm.tsv --common-segs
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_ctx_gm.stats --tablefile=mtr_ctx_gm.tsv --common-segs
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_lh_ctx_gm.stats --tablefile=mtr_lh_ctx_gm.tsv --common-segs
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_rh_ctx_gm.stats --tablefile=mtr_rh_ctx_gm.tsv --common-segs
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_subcort_gm.stats --tablefile=mtr_subcort_gm.tsv --common-segs
asegstats2table --subjectsfile=subjectsfile.txt --meas mean --stats=mtr_AD_sig.stats --tablefile=mtr_AD_sig.tsv --common-segs
