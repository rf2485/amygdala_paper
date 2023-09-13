#!/bin/bash

dirS=/Volumes/Research/lazarm03lab/labspace/AD/camcan995/raw

mkdir -p dwi_preprocessing/
cut -f1 dwi_over_55.tsv > dwi_preprocessing/subjectsfile.txt
sed -i '' '1d' dwi_preprocessing/subjectsfile.txt

for j in $(cut -f1 dwi_preprocessing/subjectsfile.txt); do
  if [ ! -f dwi_preprocessing/${j}/b0_brain_mask.nii.gz ]; then
  	mkdir -p dwi_preprocessing/${j}
    cp $dirS/${j}/dwi/${j}_dwi.* dwi_preprocessing/${j}/
    cd dwi_preprocessing/${j}/
  	dwidenoise ${j}_dwi.nii.gz dwi_denoised.nii
	
  	mrcalc ${j}_dwi.nii.gz dwi_denoised.nii -subtract dwi_residuals.nii
  	mrdegibbs dwi_denoised.nii dwi_degibbs.nii

  	mcflirt -in dwi_degibbs.nii -out dwi_corr.nii.gz
    fslroi dwi_corr.nii.gz b0 0 1 #pull b0 image from corrected dwi
    fslmaths b0.nii.gz -Tmean b0_mean.nii.gz #mean across time
    bet b0_mean b0_brain -f 0.2 -g 0 -n -m #generate brain mask
    cd ../..
  fi
done

