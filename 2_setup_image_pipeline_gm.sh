basedir=/gpfs/data/lazarlab/CamCan995/
projectdir=$basedir/derivatives/scd/main/

module load singularity/3.9.8
module load miniconda3/gpu/4.9.2

singularity pull docker://dmri/neurodock:v1.0.0
singularity pull docker://nyudiffusionmri/designer2:v2.0.10
singularity pull docker://leonyichencai/synb0-disco:v3.1
mv *.sif $basedir
conda create -c https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public/ -c conda-forge -n fsl_eddy fsl-topup==2203.5 fsl-avwutils==2209.2 fsl-fdt==2202.10 fsl-tbss==2111.2 fsl-bet2==2111.8 fsl-eddy-cuda-10.2==2401.2

conda create -n amico python=3.11
source activate ~/.conda/envs/amico/
pip install dmri-amico==2.0.3
conda deactivate