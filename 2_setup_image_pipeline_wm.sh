basedir=/gpfs/data/lazarlab/CamCan995/
projectdir=$basedir/derivatives/designer_tbss/

module load singularity/3.9.8
module load miniconda3/gpu/4.9.2

cd $basedir
singularity pull docker://dmri/neurodock:1.0.0
singularity pull docker://nyudiffusionmri/designer2:v2.0.10
singularity pull docker://leonyichencai/synb0-disco:v3.1

conda create -c https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public/ -c conda-forge -n fsl_eddy fsl-topup==2203.5 fsl-avwutils==2209.2 fsl-fdt==2202.10 fsl-tbss==2111.2 fsl-bet2==2111.8 fsl-eddy-cuda-10.2==2401.2 fsl-miscvis==2201.1 fsl-data_standard==2208.0 fsl-sub==2.8.4 fsl-sub-plugin-slurm==1.6.5 fsl-randomise==2203.3
cp $projectdir/fsl_sub.yml ~/.fsl_sub.yml

conda create -n amico python=3.11
source activate ~/.conda/envs/amico/
pip install dmri-amico==2.0.3
conda deactivate

module load git
git clone https://github.com/NYU-DiffusionMRI/SMI.git