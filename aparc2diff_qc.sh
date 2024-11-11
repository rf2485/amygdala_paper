basedir=/Volumes/Research/lazarm03lab/labspace/AD/camcan995/
projectdir=$basedir/derivatives/scd/gm_roi/
t1dir=$projectdir/freesurfer/
dwidir=$projectdir/dwi_processed/

subj_list=$(cut -f1 $projectdir/dwi_over_55.tsv)
subj_list=($subj_list)

rm -rf $t1dir/group_qc/2_aparc2diff
mkdir -p $t1dir/group_qc/2_aparc2diff
echo '<HTML><TITLE>aparc2diff</TITLE><BODY BGCOLOR="#aaaaff">' > $t1dir/group_qc/2_aparc2diff/index.html

cd $t1dir/group_qc/2_aparc2diff
for subj in "${subj_list[@]}"; do
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 38 48 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grota.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grotb.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'x' -slice 58 48 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grotc.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 38 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grotd.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grote.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'y' -slice 48 58 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grotf.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 26 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/grotg.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 33 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/groth.png
	freeview -v $t1dir/$subj/diffusion/B0.nii $t1dir/$subj/diffusion/aparc+aseg2diff_NODDI_mask.mgz:colormap=lut:opacity=0.2 -viewport 'z' -slice 48 48 40 -nocursor -screenshot $t1dir/group_qc/2_aparc2diff/groti.png
	pngappend grota.png + grotb.png + grotc.png + grotd.png + grote.png + grotf.png + grotg.png + groth.png + groti.png $subj.png
	echo '<a href="'${subj}'.png"><img src="'${subj}'.png" WIDTH='1000' >' ${subj}'</a><br>' >> $t1dir/group_qc/2_aparc2diff/index.html
done

echo '</BODY></HTML>' >> $t1dir/group_qc/2_aparc2diff/index.html