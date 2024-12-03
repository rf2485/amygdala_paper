source("1_data_preparation.R")
library(arsenal)

failed_qc <- c('sub-CC510255', #abnormality in left temporal pole
               'sub-CC510438', #abnormality in left frontal lobe
               # 'sub-CC610308', #parietal lobe cutoff
               # 'sub-CC610469', #parietal lobe cutoff
               # 'sub-CC620466', #parietal lobe cutoff
               'sub-CC620821', #segmentation errors from large ventricles
               'sub-CC621011', #segmentation errors from large ventricles
               'sub-CC621080', #segmentation errors
               # 'sub-CC710214', #parietal lobe cutoff
               'sub-CC710551', #motion artifacts in DWI
               'sub-CC711027', #severe motion artifacts in T1
               # 'sub-CC712027', #parietal lobe cutoff
               'sub-CC721434' #segmentation errors from large ventricles
)

#define remove_outliers function
remove_outliers <- function(x, na.rm = TRUE) 
{
  ## Find 25% and 75% Quantiles using inbuild function
  quant <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  
  ## Find Interquantile range and multiply it by 1.5 
  ## to derive factor for range calculation
  H <- 1.5 * IQR(x, na.rm = na.rm)
  
  y <- x
  
  ## fill the outlier elements with NA
  y[x < (quant[1] - H)] <- NA
  y[x > (quant[2] + H)] <- NA
  
  y
}

#extract SCD status from dwi_over_55 table
scd_status <- dwi_over_55 %>% select(participant_id, SCD) %>%
  filter(!participant_id %in% failed_qc)
scd_status$SCD <- factor(scd_status$SCD,
                         levels = c(1,0),
                         labels = c('SCD', 'Control'))

#import aseg stats table (subcortical volumes)
aseg = read_tsv("freesurfer/asegtable.tsv") %>%
  rename(participant_id=`Measure:volume`) #%>%
names(aseg) <- make.names(names(aseg))
# aseg_table <- tableby(formulize('SCD', names(aseg)[3:54]), 
#                          data = aseg, numeric.test="wt", total = F) 
# summary(aseg_table, text = TRUE)

#import JHU stats table (WM volumes)
jhu_volume <- read_tsv("freesurfer/jhu_volume.tsv") %>%
  rename(participant_id=`Measure:volume`,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
         ) %>%
  select(!starts_with("Seg00"))
# jhu_volume[jhu_volume == 0] <- NA
volumes <- left_join(scd_status, jhu_volume) %>% #join JHU volumes with SCD status
  left_join(., aseg) %>%
  mutate(across(c(3:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) %>% #normalize by intracranial volume
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
volumes[volumes == 0] <- NA
volumes_table <- tableby(formulize('SCD', names(volumes)[3:90]), 
                      data = volumes, numeric.test="wt", total = F) 
summary(volumes_table, text = TRUE)

#import left aparc stats tables (cortical thickness and volume)
lh_AD_sig_thickness <- read_tsv("freesurfer/lh_AD_sig_thickness.tsv") %>%
  rename(participant_id=`Measure:mean`, lh_AD_sig_thickness=Seg0001)
lh_aparc_thickness = read_tsv("freesurfer/lh_aparctable_thickness.tsv") %>%
  rename(participant_id=lh.aparc.thickness) %>%
  left_join(., lh_AD_sig_thickness) %>%
  select(!(lh_MeanThickness_thickness:eTIV))
#import right aparc stats tables (cortical thickness and volume)
rh_AD_sig_thickness <- read_tsv("freesurfer/rh_AD_sig_thickness.tsv") %>%
  rename(participant_id=`Measure:mean`, rh_AD_sig_thickness=Seg0001)
rh_aparc_thickness = read_tsv("freesurfer/rh_aparctable_thickness.tsv") %>%
  rename(participant_id=rh.aparc.thickness) %>%
  left_join(., rh_AD_sig_thickness) %>%
  select(!(rh_MeanThickness_thickness:eTIV))
#create aparc table
thickness <- left_join(scd_status, lh_aparc_thickness) %>% #join left cortical thickness with SCD stats
  left_join(., rh_aparc_thickness) %>% #add right cortical thickness
  mutate(across(where(is.double), remove_outliers))
thickness_table <- tableby(formulize('SCD', names(thickness)[3:72]), 
                           data = thickness, numeric.test="wt", total = F)
summary(thickness_table, text = T)

####diffusion and MTR group means
#read in diffusion and MTR tables
aparc2meas_files <- list.files(path = "freesurfer", 
                               pattern = "aparc.*aseg.*\\.*tsv", 
                               full.names = T)
for (i in 1:length(aparc2meas_files)) {
  assign(gsub(".tsv", "", 
              gsub("freesurfer.aparc.", "", make.names(aparc2meas_files[i]))), 
         read.delim(aparc2meas_files[i]))
}
AD_sig2meas_files <- list.files(path = "freesurfer", 
                                pattern = "?h_AD_sig2.*\\.*tsv", 
                                full.names = T)
for (i in 1:length(AD_sig2meas_files)) {
  assign(gsub(".tsv", "", 
              gsub("freesurfer.", "", make.names(AD_sig2meas_files[i]))), 
         read.delim(AD_sig2meas_files[i]))
}
jhu2meas_files <- list.files(path = "freesurfer",
                             pattern = "jhu2.*\\.*tsv",
                             full.names = T)
for (i in 1:length(jhu2meas_files)) {
  assign(gsub(".tsv", "",
              gsub("freesurfer.", "", make.names(jhu2meas_files[i]))),
         read.delim(jhu2meas_files[i]))
}

fit_FWF <- aseg2fit_FWF %>%
  left_join(., lh_AD_sig2fit_FWF) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2fit_FWF) %>%
  rename(rh_AD_signature = AD_signature) %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
fit_FWF[fit_FWF == 0] <- NA
FWF_table <- tableby(formulize('SCD', names(fit_FWF)[3:91]),
                     data = fit_FWF, numeric.test="wt", total = FALSE)
summary(FWF_table, text = T)

fit_NDI <- aseg2fit_NDI %>%
  left_join(., lh_AD_sig2fit_NDI) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2fit_NDI) %>%
  rename(rh_AD_signature = AD_signature) %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
fit_NDI[fit_NDI == 0] <- NA
NDI_table <- tableby(formulize('SCD', names(fit_NDI)[3:91]),
                     data = fit_NDI, numeric.test="wt", total = FALSE)
summary(NDI_table, text = T)

fit_ODI <- aseg2fit_ODI %>%
  left_join(., lh_AD_sig2fit_ODI) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2fit_ODI) %>%
  rename(rh_AD_signature = AD_signature) %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
fit_ODI[fit_ODI == 0] <- NA
ODI_table <- tableby(formulize('SCD', names(fit_ODI)[3:91]),
                     data = fit_ODI, numeric.test="wt", total = FALSE)
summary(ODI_table, text = T)

dti_fa <- aseg2dti_fa %>%
  left_join(., lh_AD_sig2dti_fa) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dti_fa) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_fa) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
dti_fa[dti_fa == 0] <- NA
fa_table <- tableby(formulize('SCD', names(dti_fa)[3:147]),
                     data = dti_fa, numeric.test="wt", total = FALSE)
summary(fa_table, text = T)

dti_md <- aseg2dti_md %>%
  left_join(., lh_AD_sig2dti_md) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dti_md) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_md) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
dti_md[dti_md == 0] <- NA
md_table <- tableby(formulize('SCD', names(dti_md)[3:147]),
                    data = dti_md, numeric.test="wt", total = FALSE)
summary(md_table, text = T)

dti_rd <- aseg2dti_rd %>%
  left_join(., lh_AD_sig2dti_rd) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dti_rd) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_rd) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
dti_rd[dti_rd == 0] <- NA
rd_table <- tableby(formulize('SCD', names(dti_rd)[3:147]),
                    data = dti_rd, numeric.test="wt", total = FALSE)
summary(rd_table, text = T)

dti_ad <- aseg2dti_ad %>%
  left_join(., lh_AD_sig2dti_ad) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dti_ad) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_ad) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
dti_ad[dti_ad == 0] <- NA
ad_table <- tableby(formulize('SCD', names(dti_ad)[3:147]),
                    data = dti_ad, numeric.test="wt", total = FALSE)
summary(ad_table, text = T)

dki_kfa <- aseg2dki_kfa %>%
  left_join(., lh_AD_sig2dki_kfa) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dki_kfa) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_kfa) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
dki_kfa[dki_kfa == 0] <- NA
kfa_table <- tableby(formulize('SCD', names(dki_kfa)[3:147]),
                    data = dki_kfa, numeric.test="wt", total = FALSE)
summary(kfa_table, text = T)

dki_mk <- aseg2dki_mk %>%
  left_join(., lh_AD_sig2dki_mk) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dki_mk) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_mk) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
dki_mk[dki_mk == 0] <- NA
mk_table <- tableby(formulize('SCD', names(dki_mk)[3:147]),
                     data = dki_mk, numeric.test="wt", total = FALSE)
summary(mk_table, text = T)

dki_rk <- aseg2dki_rk %>%
  left_join(., lh_AD_sig2dki_rk) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dki_rk) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_rk) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
dki_rk[dki_rk == 0] <- NA
rk_table <- tableby(formulize('SCD', names(dki_rk)[3:147]),
                     data = dki_rk, numeric.test="wt", total = FALSE)
summary(rk_table, text = T)

dki_ak <- aseg2dki_ak %>%
  left_join(., lh_AD_sig2dki_ak) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dki_ak) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_ak) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
dki_ak[dki_ak == 0] <- NA
ak_table <- tableby(formulize('SCD', names(dki_ak)[3:147]),
                     data = dki_ak, numeric.test="wt", total = FALSE)
summary(ak_table, text = T)

smi_matlab_Da <- aseg2smi_matlab_Da %>%
  left_join(., jhu2smi_matlab_Da) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
smi_matlab_Da[smi_matlab_Da == 0] <- NA
Da_table <- tableby(formulize('SCD', names(smi_matlab_Da)[3:56]),
                     data = smi_matlab_Da, numeric.test="wt", total = FALSE)
summary(Da_table, text = T)

smi_matlab_DePar <- aseg2smi_matlab_DePar %>%
  left_join(., jhu2smi_matlab_DePar) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
smi_matlab_DePar[smi_matlab_DePar == 0] <- NA
DePar_table <- tableby(formulize('SCD', names(smi_matlab_DePar)[3:56]),
                    data = smi_matlab_DePar, numeric.test="wt", total = FALSE)
summary(DePar_table, text = T)

smi_matlab_DePerp <- aseg2smi_matlab_DePerp %>%
  left_join(., jhu2smi_matlab_DePerp) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
smi_matlab_DePerp[smi_matlab_DePerp == 0] <- NA
DePerp_table <- tableby(formulize('SCD', names(smi_matlab_DePerp)[3:56]),
                       data = smi_matlab_DePerp, numeric.test="wt", total = FALSE)
summary(DePerp_table, text = T)

smi_matlab_f <- aseg2smi_matlab_f %>%
  left_join(., jhu2smi_matlab_f) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
smi_matlab_f[smi_matlab_f == 0] <- NA
f_table <- tableby(formulize('SCD', names(smi_matlab_f)[3:56]),
                       data = smi_matlab_f, numeric.test="wt", total = FALSE)
summary(f_table, text = T)

smi_matlab_p2 <- aseg2smi_matlab_p2 %>%
  left_join(., jhu2smi_matlab_p2) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
smi_matlab_p2[smi_matlab_p2 == 0] <- NA
p2_table <- tableby(formulize('SCD', names(smi_matlab_p2)[3:56]),
                       data = smi_matlab_p2, numeric.test="wt", total = FALSE)
summary(p2_table, text = T)

#extract SCD status from mti_over_55 table for each TR
scd_status_tr30 <- mti_over_55_tr30 %>% select(participant_id, SCD) %>%
  filter(!participant_id %in% failed_qc)
scd_status_tr30$SCD <- factor(scd_status_tr30$SCD,
                              levels = c(1,0),
                              labels = c('SCD', 'Control'))

scd_status_tr50 <- mti_over_55_tr50 %>% select(participant_id, SCD) %>%
  filter(!participant_id %in% failed_qc)
scd_status_tr50$SCD <- factor(scd_status_tr50$SCD,
                              levels = c(1,0),
                              labels = c('SCD', 'Control'))

#mtr stats for each TR
mtr_tr30 <- aseg2mtr %>%
  left_join(., lh_AD_sig2mtr) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2mtr) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2mtr) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status_tr30, .) %>%
  select(!CSF & !ends_with("Ventricle") & !ends_with("Vent") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
mtr_tr30[mtr_tr30 == 0] <- NA
mtr_tr30_table <- tableby(formulize('SCD', names(mtr_tr30)[3:145]),
                          data = mtr_tr30, numeric.test="wt", total = FALSE)
summary(mtr_tr30_table, text = T)

g_ratio_tr30 <- aseg2g_ratio %>%
  left_join(., jhu2g_ratio) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status_tr30, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
g_ratio_tr30[g_ratio_tr30 == 0] <- NA
g_ratio_tr30$SCD <- factor(g_ratio_tr30$SCD,
                         levels = c('SCD', 'Control'),
                         labels = c('SCD', 'Control'))
g_ratio_tr30_table <- tableby(formulize('SCD', names(g_ratio_tr30)[3:56]),
                          data = g_ratio_tr30, numeric.test="wt", total = FALSE
                          )
summary(g_ratio_tr30_table, text = T)

mtr_tr50 <- aseg2mtr %>%
  left_join(., lh_AD_sig2mtr) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2mtr) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2mtr) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status_tr50, .) %>%
  select(!CSF & !ends_with("Ventricle") & !ends_with("Vent") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
mtr_tr50[mtr_tr50 == 0] <- NA
mtr_tr50_table <- tableby(formulize('SCD', names(mtr_tr50)[3:145]),
                          data = mtr_tr50, numeric.test="wt", total = FALSE)
summary(mtr_tr50_table, text = T)

g_ratio_tr50 <- aseg2g_ratio %>%
  left_join(., jhu2g_ratio) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         fornix_column_body=Seg0006,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049,
         tapetum_L=Seg0050
  ) %>%
  left_join(scd_status_tr50, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
g_ratio_tr50[g_ratio_tr50 == 0] <- NA
g_ratio_tr50$SCD <- factor(g_ratio_tr50$SCD,
                           levels = c('SCD', 'Control'),
                           labels = c('SCD', 'Control'))
g_ratio_tr50_table <- tableby(formulize('SCD', names(g_ratio_tr50)[3:56]),
                              data = g_ratio_tr50, numeric.test="kwt", total = FALSE
)
summary(g_ratio_tr50_table, text = T)

# results_tables <- c(aseg_table, thickness_table, FWF_table, NDI_table, ODI_table,
#                     fa_table, md_table, rd_table, ad_table, kfa_table, mk_table, 
#                     rk_table, ak_table, mtr_tr30_table, mtr_tr50_table)
# write2word(results_tables, "results.docx")
