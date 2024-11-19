import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import amico
import sys
study="/gpfs/data/lazarlab/CamCan995/derivatives/scd/gm_roi/dwi_processed/"
subject=sys.argv[1]
amico.core.setup()
ae = amico.Evaluation(study, subject)
bvals = os.path.join(study, subject, "dwi_preprocessed.bval")
bvecs = os.path.join(study, subject, "dwi_preprocessed.bvec")
amico.util.fsl2scheme(bvals, bvecs)

ae.load_data(dwi_filename = "dwi_preprocessed.nii", scheme_filename = "dwi_preprocessed.scheme", mask_filename = "brain_mask.nii", b0_thr = 0)
ae.set_model("NODDI")
ae.model.set(dPar=1.1E-3)
ae.generate_kernels()
ae.load_kernels()
ae.fit()
ae.save_results()