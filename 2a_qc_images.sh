# basedir=/gpfs/data/lazarlab/CamCan995/
# projectdir=$basedir/derivatives/mti_whole_wm_gm/
#
# module load freesurfer/7.4.1
# module load fsl/6.0.7

projectdir=/Volumes/Research/lazarm03lab/labspace/AD/camcan995/derivatives/scd/mti_whole_wm_gm

subj_list=$(cut -f1 -d$'\t' $projectdir/mti_over_55.tsv)
subj_list=($subj_list)
freesurferdir=$projectdir/freesurfer
export SUBJECTS_DIR=$freesurferdir
cd $freesurferdir

mkdir -p $freesurferdir/group_qc/recon
echo '<HTML><TITLE>recon</TITLE><BODY BGCOLOR="#aaaaff">' > group_qc/recon/index.html
# for subj in "${(f)subj_list}"; do
for subj in "${subj_list[@]}"; do
	if [ ! -f group_qc/recon/$subj.png ]; then
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'x' -slice 102 128 128 -nocursor -screenshot group_qc/recon/grota.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'x' -slice 128 128 128 -nocursor -screenshot group_qc/recon/grotb.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'x' -slice 154 128 128 -nocursor -screenshot group_qc/recon/grotc.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'y' -slice 128 128 102 -nocursor -screenshot group_qc/recon/grotd.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'y' -slice 128 128 128 -nocursor -screenshot group_qc/recon/grote.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'y' -slice 128 128 154 -nocursor -screenshot group_qc/recon/grotf.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'z' -slice 128 102 128 -nocursor -screenshot group_qc/recon/grotg.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'z' -slice 128 128 128 -nocursor -screenshot group_qc/recon/groth.png
		freeview -v $subj/mri/T1.mgz -f $subj/surf/lh.white:edgecolor=blue $subj/surf/rh.white:edgecolor=blue $subj/surf/lh.pial:edgecolor=red $subj/surf/rh.pial:edgecolor=red -viewport 'z' -slice 128 154 128 -nocursor -screenshot group_qc/recon/groti.png
		pngappend group_qc/recon/grota.png + group_qc/recon/grotb.png + group_qc/recon/grotc.png + group_qc/recon/grotd.png + group_qc/recon/grote.png + group_qc/recon/grotf.png + group_qc/recon/grotg.png + group_qc/recon/groth.png + group_qc/recon/groti.png group_qc/recon/$subj.png
	fi
	echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> group_qc/recon/index.html
done

echo '</BODY></HTML>' >> group_qc/recon/index.html

rm -rf group_qc/1_mti_reg
mkdir -p group_qc/1_mti_reg

for subj in "${subj_list[@]}"; do
# for subj in "${(f)subj_list}"; do
	echo $subj
	cp $subj/mti/intermediate/1_mti_reg.nii.gz group_qc/1_mti_reg/${subj}_1_mti_reg.nii.gz
	cp $subj/mti/mti_bl_raw.nii.gz group_qc/1_mti_reg/${subj}_2_mti_bl_raw.nii.gz
done
slicesdir -o group_qc/1_mti_reg/*.nii.gz
rm -rf group_qc/1_mti_reg
mv slicesdir group_qc/1_mti_reg

slicesdir */mti/intermediate/2_mti_bl_degibbs.nii.gz
rm -rf group_qc/2_mti_bl_degibbs
mv slicesdir group_qc/2_mti_bl_degibbs

slicesdir */mti/intermediate/2_mti_degibbs.nii.gz
rm -rf group_qc/2_mti_degibbs
mv slicesdir group_qc/2_mti_degibbs

slicesdir */mti/mtr.nii.gz
rm -rf group_qc/mtr
mv slicesdir group_qc/mtr

mkdir -p group_qc/aparc+aseg2mtr
echo '<HTML><TITLE>aparc+aseg2mtr</TITLE><BODY BGCOLOR="#aaaaff">' > group_qc/aparc+aseg2mtr/index.html
for subj in "${subj_list[@]}"; do
# for subj in "${(f)subj_list}"; do
	if [ ! -f group_qc/aparc+aseg2mtr/$subj.png ]; then
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 42 64 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grota.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 52 64 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grotb.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 62 64 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grotc.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 51 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grotd.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 64 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grote.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 52 77 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/grotf.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 51 -nocursor -screenshot group_qc/aparc+aseg2mtr/grotg.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 64 -nocursor -screenshot group_qc/aparc+aseg2mtr/groth.png
		freeview -v $subj/mti/mtr.nii.gz $subj/mti/aparc+aseg2mtr.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 52 64 77 -nocursor -screenshot group_qc/aparc+aseg2mtr/groti.png
		pngappend group_qc/aparc+aseg2mtr/grota.png + group_qc/aparc+aseg2mtr/grotb.png + group_qc/aparc+aseg2mtr/grotc.png + group_qc/aparc+aseg2mtr/grotd.png + group_qc/aparc+aseg2mtr/grote.png + group_qc/aparc+aseg2mtr/grotf.png + group_qc/aparc+aseg2mtr/grotg.png + group_qc/aparc+aseg2mtr/groth.png + group_qc/aparc+aseg2mtr/groti.png group_qc/aparc+aseg2mtr/$subj.png
	fi
	echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> group_qc/aparc+aseg2mtr/index.html
done

echo '</BODY></HTML>' >> group_qc/aparc+aseg2mtr/index.html