source("1_data_preparation.R")
library(gtsummary)

failed_qc <- c('sub-CC510255', #SCD abnormality in left temporal pole
               'sub-CC510438', #CTL abnormality in left frontal lobe
               # 'sub-CC610308', #parietal lobe cutoff
               # 'sub-CC610469', #parietal lobe cutoff
               # 'sub-CC620466', #parietal lobe cutoff
               'sub-CC620821', #SCD segmentation errors from large ventricles
               'sub-CC621011', #CTL segmentation errors from large ventricles
               'sub-CC621080', #SCD segmentation errors
               # 'sub-CC710214', #parietal lobe cutoff
               'sub-CC710551', #CTL motion artifacts in DWI
               'sub-CC711027', #SCD severe motion artifacts in T1
               # 'sub-CC712027', #parietal lobe cutoff
               'sub-CC721434' #CTL segmentation errors from large ventricles
)

subcort_gm <- c("SCD", "Left.Thalamus", "Right.Thalamus", "Left.Caudate", "Right.Caudate",
                "Left.Putamen", "Right.Putamen", "Left.Pallidum", "Right.Pallidum",
                "Left.Hippocampus", "Right.Hippocampus", "Left.Amygdala", "Right.Amygdala",
                "Left.Accumbens.area", "Right.Accumbens.area", "Left.VentralDC", "Right.VentralDC")
# ctx_gm_exclude <- c("participant_id", "lh_AD_signature", "rh_AD_signature",
#                     "Left.Thalamus", "Right.Thalamus", "Left.Caudate", "Right.Caudate",
#                 "Left.Putamen", "Right.Putamen", "Left.Pallidum", "Right.Pallidum",
#                 "Left.Hippocampus", "Right.Hippocampus", "Left.Amygdala", "Right.Amygdala",
#                 "Left.Accumbens.area", "Right.Accumbens.area", "Left.VentralDC", "Right.VentralDC",
#                 "CC_Posterior", "CC_Mid_Posterior", "CC_Central", "CC_Mid_Anterior", "CC_Anterior")
jhu_lookup <- c(participant_id="Measure:volume",
                participant_id="Measure.mean",
                middle_cerebellar_peduncle="Seg0001",
                genu_corpus_callosum="Seg0003",
                body_corpus_callosum="Seg0004",
                splenium_corpus_callosum="Seg0005",
                fornix_column_body="Seg0006",
                inferior_cerebellar_peduncle_R="Seg0011",
                inferior_cerebellar_peduncle_L="Seg0012",
                superior_cerebellar_peduncle_R="Seg0013",
                superior_cerebellar_peduncle_L="Seg0014",
                cerebral_peduncle_R="Seg0015",
                cerebral_peduncle_L="Seg0016",
                anterior_limb_internal_capsule_R="Seg0017",
                anterior_limb_internal_capsule_L="Seg0018",
                posterior_limb_internal_capsule_R="Seg0019",
                posterior_limb_internal_capsule_L="Seg0020",
                retrolenticular_part_internal_capsule_R="Seg0021",
                retrolenticular_part_internal_capsule_L="Seg0022",
                anterior_corona_radiata_R="Seg0023",
                anterior_corona_radiata_L="Seg0024",
                superior_corona_radiata_R="Seg0025",
                superior_corona_radiata_L="Seg0026",
                posterior_corona_radiata_R="Seg0027",
                posterior_corona_radiata_L="Seg0028",
                posterior_thalamic_radiation_R="Seg0029",
                posterior_thalamic_radiation_L="Seg0030",
                sagittal_stratum_R="Seg0031",
                sagittal_stratum_L="Seg0032",
                external_capsule_R="Seg0033",
                external_capsule_L="Seg0034",
                upper_cingulum_R="Seg0035",
                upper_cingulum_L="Seg0036",
                lower_cingulum_R="Seg0037",
                lower_cingulum_L="Seg0038",
                fornix_cres_R="Seg0039",
                fornix_cres_L="Seg0040",
                superior_longitudinal_fasciculus_R="Seg0041",
                superior_longitudinal_fasciculus_L="Seg0042",
                superior_fronto_occipital_fasciculus_R="Seg0043",
                superior_fronto_occipital_fasciculus_L="Seg0044",
                inferior_fronto_occipital_fasciculus_R="Seg0045",
                inferior_fronto_occipital_fasciculus_L="Seg0046",
                uncinate_fasciculus_R="Seg0047",
                uncinate_fasciculus_L="Seg0048",
                tapetum_R="Seg0049",
                tapetum_L="Seg0050")

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
  rename(participant_id=`Measure:volume`)
names(aseg) <- make.names(names(aseg))

#import JHU stats table (WM volumes)
jhu_volume <- read_tsv("freesurfer/jhu_volume.tsv") %>%
  rename(any_of(jhu_lookup)) %>%
  select(!starts_with("Seg00"))
volumes <- left_join(scd_status, jhu_volume) %>% #join JHU volumes with SCD status
  left_join(., aseg) %>%
  mutate(across(c(3:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) %>% #normalize by intracranial volume
  mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  rename(Left.Cerebral.White.Matter	= lhCerebralWhiteMatterVol, 
         Right.Cerebral.White.Matter	= rhCerebralWhiteMatterVol) %>%
  select(!c((BrainSegVol:CortexVol), (CerebralWhiteMatterVol:EstimatedTotalIntraCranialVol)))
subcort_gm_volumes_table <- volumes %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_volumes_table <- volumes %>% 
  select(c(SCD:middle_cerebellar_peduncle, 
           inferior_cerebellar_peduncle_R:tapetum_L, 
           CC_Posterior:CC_Anterior)) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

#import left aparc stats tables (cortical thickness and volume)
lh_AD_sig_thickness <- read_tsv("freesurfer/lh_AD_sig_thickness.tsv") %>%
  rename(participant_id=`Measure:mean`, lh_AD_signature=Seg0001)
lh_aparc_thickness = read_tsv("freesurfer/lh_aparctable_thickness.tsv") %>%
  rename(participant_id=lh.aparc.thickness) %>%
  select(!(lh_MeanThickness_thickness:eTIV))
#import right aparc stats tables (cortical thickness and volume)
rh_AD_sig_thickness <- read_tsv("freesurfer/rh_AD_sig_thickness.tsv") %>%
  rename(participant_id=`Measure:mean`, rh_AD_signature=Seg0001)
rh_aparc_thickness = read_tsv("freesurfer/rh_aparctable_thickness.tsv") %>%
  rename(participant_id=rh.aparc.thickness) %>%
  select(!(rh_MeanThickness_thickness:eTIV))
#create aparc table
thickness <- left_join(scd_status, lh_aparc_thickness) %>% #join left cortical thickness with SCD stats
  left_join(., rh_aparc_thickness) %>% #add right cortical thickness
  mutate(across(where(is.double), remove_outliers)) %>%
  rename_with(., ~ gsub("_", ".", .x, fixed = T), ends_with("thickness")) %>%
  rename_with(., ~ paste0("ctx.", .x, recycle0 = T), ends_with("thickness")) %>%
  rename_with(., ~ gsub(".thickness", "", .x, fixed = T)) %>%
  left_join(., lh_AD_sig_thickness) %>%
  left_join(., rh_AD_sig_thickness)
ctx_gm_thickness_table <- thickness %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

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
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_FWF_table <- fit_FWF %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_FWF_table <- fit_FWF %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

fit_NDI <- aseg2fit_NDI %>%
  left_join(., lh_AD_sig2fit_NDI) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2fit_NDI) %>%
  rename(rh_AD_signature = AD_signature) %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_NDI_table <- fit_NDI %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_NDI_table <- fit_NDI %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

fit_ODI <- aseg2fit_ODI %>%
  left_join(., lh_AD_sig2fit_ODI) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2fit_ODI) %>%
  rename(rh_AD_signature = AD_signature) %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_ODI_table <- fit_ODI %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_ODI_table <- fit_ODI %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

dti_fa <- aseg2dti_fa %>%
  left_join(., lh_AD_sig2dti_fa) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dti_fa) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_fa) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_fa_table <- dti_fa %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_fa_table <- dti_fa %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_fa_table <- dti_fa %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
           )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

dti_md <- aseg2dti_md %>%
  left_join(., lh_AD_sig2dti_md) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dti_md) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_md) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_md_table <- dti_md %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_md_table <- dti_md %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_md_table <- dti_md %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

dti_rd <- aseg2dti_rd %>%
  left_join(., lh_AD_sig2dti_rd) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dti_rd) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_rd) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_rd_table <- dti_rd %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_rd_table <- dti_rd %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_rd_table <- dti_rd %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

dti_ad <- aseg2dti_ad %>%
  left_join(., lh_AD_sig2dti_ad) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dti_ad) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_ad) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_ad_table <- dti_ad %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_ad_table <- dti_ad %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_ad_table <- dti_ad %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing_text = "Excluded Outliers") %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

dki_kfa <- aseg2dki_kfa %>%
  left_join(., lh_AD_sig2dki_kfa) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dki_kfa) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_kfa) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_kfa_table <- dki_kfa %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_kfa_table <- dki_kfa %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_kfa_table <- dki_kfa %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

dki_mk <- aseg2dki_mk %>%
  left_join(., lh_AD_sig2dki_mk) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dki_mk) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_mk) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_mk_table <- dki_mk %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_mk_table <- dki_mk %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_mk_table <- dki_mk %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

dki_rk <- aseg2dki_rk %>%
  left_join(., lh_AD_sig2dki_rk) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dki_rk) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_rk) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_rk_table <- dki_rk %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_rk_table <- dki_rk %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_rk_table <- dki_rk %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

dki_ak <- aseg2dki_ak %>%
  left_join(., lh_AD_sig2dki_ak) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2dki_ak) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_ak) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_ak_table <- dki_ak %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_ak_table <- dki_ak %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_ak_table <- dki_ak %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

smi_matlab_Da <- aseg2smi_matlab_Da %>%
  left_join(., jhu2smi_matlab_Da) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_Da_table <- smi_matlab_Da %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

smi_matlab_DePar <- aseg2smi_matlab_DePar %>%
  left_join(., jhu2smi_matlab_DePar) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_DePar_table <- smi_matlab_DePar %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

smi_matlab_DePerp <- aseg2smi_matlab_DePerp %>%
  left_join(., jhu2smi_matlab_DePerp) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_DePerp_table <- smi_matlab_DePerp %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

smi_matlab_f <- aseg2smi_matlab_f %>%
  left_join(., jhu2smi_matlab_f) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_f_table <- smi_matlab_f %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

smi_matlab_p2 <- aseg2smi_matlab_p2 %>%
  left_join(., jhu2smi_matlab_p2) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_p2_table <- smi_matlab_p2 %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

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
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status_tr30, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_mtr_tr30_table <- mtr_tr30 %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_mtr_tr30_table <- mtr_tr30 %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_mtr_tr30_table <- mtr_tr30 %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

g_ratio_tr30 <- aseg2g_ratio %>%
  left_join(., jhu2g_ratio) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status_tr30, .) %>%
  select(!starts_with("Seg00")) #%>%
  # mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_g_ratio_tr30_table <- g_ratio_tr30 %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

mtr_tr50 <- aseg2mtr %>%
  left_join(., lh_AD_sig2mtr) %>%
  rename(lh_AD_signature = AD_signature) %>%
  left_join(., rh_AD_sig2mtr) %>%
  rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2mtr) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status_tr50, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_mtr_tr50_table <- mtr_tr50 %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
ctx_gm_mtr_tr50_table <- mtr_tr50 %>% 
  select(SCD, ctx.lh.bankssts:ctx.rh.insula) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",               
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
wm_mtr_tr50_table <- mtr_tr50 %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
              ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

g_ratio_tr50 <- aseg2g_ratio %>%
  left_join(., jhu2g_ratio) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status_tr50, .) %>%
  select(!starts_with("Seg00")) %>%
  filter(participant_id!="sub-CC520083") %>% #no diffusion available
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_g_ratio_tr50_table <- g_ratio_tr50 %>%
  select(c(SCD, CC_Posterior:CC_Anterior, middle_cerebellar_peduncle,
           fornix_column_body:tapetum_L
  )) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no"
              # missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

mtr_wm <- read_tsv("freesurfer/mtr_wm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_wm=Seg0001)
mtr_gm <- read_tsv("freesurfer/mtr_gm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_gm=Seg0001)
mtr_wm_gm_ratio <- mti_over_55 %>%
  select(participant_id, coil, mt_tr) %>%
  left_join(., mtr_wm) %>%
  left_join(., mtr_gm)
mtr_wm_gm_ratio$mtr_wm_gm_ratio <- mtr_wm_gm_ratio$mtr_wm / mtr_wm_gm_ratio$mtr_gm
mtr_wm_gm_ratio$coil <- as.factor(mtr_wm_gm_ratio$coil)
# mtr_wm_gm_ratio$mt_tr <- as.factor(mtr_wm_gm_ratio$mt_tr)
mtr_wm_gm_ratio$mt_tr <- factor(mtr_wm_gm_ratio$mt_tr,
                                levels = c(30, 50),
                                labels = c("TR=30ms", "TR=50ms"))
mtr_wm_gm_ratio %>% select(mt_tr, mtr_wm_gm_ratio) %>% tbl_summary(by = mt_tr, missing = "no") %>% add_p()

tbl_merge(list(subcort_gm_volumes_table,
               subcort_gm_fa_table, subcort_gm_md_table, subcort_gm_ad_table, subcort_gm_rd_table,
               subcort_gm_kfa_table, subcort_gm_mk_table, subcort_gm_rk_table,
               subcort_gm_FWF_table, subcort_gm_ODI_table,
               subcort_gm_mtr_tr50_table), 
          tab_spanner = c("**Volume**",
                          "**FA**", "**MD**", "**AD**", "**RD**",
                          "**KFA**", "**MK**", "**RK**",
                          "**FWF**", "**ODI**",
                          "**MTR TR=50ms**"
          )) %>% 
  as_gt() %>% 
  gt::gtsave(filename = "subcort_gm_table.html")

tbl_merge(list(ctx_gm_thickness_table, 
               ctx_gm_fa_table, ctx_gm_md_table, ctx_gm_ad_table, ctx_gm_rd_table,
               ctx_gm_kfa_table, ctx_gm_mk_table,
               ctx_gm_FWF_table, ctx_gm_NDI_table,
               ctx_gm_mtr_tr50_table, ctx_gm_mtr_tr30_table), 
          tab_spanner = c("**Thickness**", 
                          "**FA**", "**MD**", "**AD**", "**RD**",
                          "**KFA**", "**MK**",
                          "**FWF**", "**NDI**",
                          "**MTR TR=50ms**", "**MTR TR=30ms**"
          )) %>% 
  as_gt() %>% 
  gt::gtsave(filename = "ctx_gm_table.html")

tbl_merge(list(wm_volumes_table,
               wm_fa_table, wm_md_table, wm_ad_table, wm_rd_table,
               wm_kfa_table, wm_mk_table, wm_ak_table, wm_rk_table,
               wm_Da_table, wm_DePar_table, wm_DePerp_table, wm_f_table, wm_p2_table,
               wm_mtr_tr30_table, wm_g_ratio_tr30_table, wm_mtr_tr50_table, wm_g_ratio_tr50_table),
          tab_spanner = c("**Volume**",
                          "**FA**", "**MD**", "**AD**", "**RD**",
                          "**KFA**", "**MK**", "**AK**", "**RK**",
                          "**Da**", "**DePar**", "**DePerp**", "**f**", "**p2**",
                          "**MTR TR=30ms**", "**g-ratio TR=30ms**", "**MTR TR=50ms**", "**g-ratio TR=50ms**"
                          )) %>% 
  as_gt() %>% 
  gt::gtsave(filename = "wm_table.html")

wm_dti_table <- tbl_merge(list(
  wm_fa_table, wm_md_table
),
tab_spanner = c("**FA**", "**MD**")
)

wm_dki_table <- tbl_merge(list(
  wm_kfa_table, wm_mk_table
),
tab_spanner = c("**KFA**", "**MK**")
)

wm_smi_table_ax <- tbl_merge(list(
  wm_Da_table, wm_f_table
),
tab_spanner = c("**Axonal Diffusion**", "**Axonal Fraction**")
)

wm_smi_table_De <- tbl_merge(list(
  wm_DePar_table, wm_DePerp_table
),
tab_spanner = c("**Parallel Extracellular Diffusion**", "**Perpendicular Extracellular Diffusion**")
)

wm_mti_table <- tbl_merge(list(
  wm_mtr_tr30_table, wm_mtr_tr50_table
),
tab_spanner = c("**MTR TR=30ms**", "**MTR TR=50ms**")
)

wm_g_ratio_table <- tbl_merge(list(
  wm_g_ratio_tr30_table, wm_g_ratio_tr50_table
),
tab_spanner = c("**g-ratio TR=30ms**", "**g-ratio TR=50ms**")
)

gm_vol_thick_table <- tbl_merge(list(
  subcort_gm_volumes_table, ctx_gm_thickness_table
),
tab_spanner = c("**Subcortical Volume**", "**Cortical Thickness**")
)

gm_fa_table <- tbl_stack(list(subcort_gm_fa_table, ctx_gm_fa_table))
gm_md_table <- tbl_stack(list(subcort_gm_md_table, ctx_gm_md_table))
gm_dti_table <- tbl_merge(list(
  gm_fa_table, gm_md_table
),
tab_spanner = c("**FA**", "**MD**")
)

subcort_gm_dki_table <- tbl_merge(list(
  subcort_gm_kfa_table, subcort_gm_mk_table
),
tab_spanner = c("**KFA**", "**MK**")
)

ctx_gm_dki_table <- tbl_merge(list(
  ctx_gm_kfa_table, ctx_gm_mk_table
),
tab_spanner = c("**KFA**", "**MK**")
)

gm_FWF_table <- tbl_stack(list(subcort_gm_FWF_table, ctx_gm_FWF_table))
gm_NDI_table <- tbl_stack(list(subcort_gm_NDI_table, ctx_gm_NDI_table))
gm_ODI_table <- tbl_stack(list(subcort_gm_ODI_table, ctx_gm_ODI_table))
gm_NDI_ODI_table <- tbl_merge(list(
  gm_ODI_table, gm_NDI_table
),
tab_spanner = c("**ODI**", "**NDI**")
)

gm_mtr_tr30_table <- tbl_stack(list(subcort_gm_mtr_tr30_table, ctx_gm_mtr_tr30_table))
gm_mtr_tr50_table <- tbl_stack(list(subcort_gm_mtr_tr50_table, ctx_gm_mtr_tr50_table))
gm_mtr_table <- tbl_merge(list(
  gm_mtr_tr50_table, gm_mtr_tr30_table
),
tab_spanner = c("**MTR TR=50ms**", "**MTR TR=30ms**")
)

library(ggpmisc)
amygdala_df <- dki_kfa %>% select(participant_id, SCD, Left.Amygdala, Right.Amygdala) %>%
  rename(kfa_lh_amygdala = Left.Amygdala, kfa_rh_amygdala = Right.Amygdala)
amygdala_df <- fit_ODI %>% select(participant_id, SCD, Left.Amygdala, Right.Amygdala) %>%
  rename(odi_lh_amygdala = Left.Amygdala, odi_rh_amygdala = Right.Amygdala) %>%
  left_join(amygdala_df, .)
amygdala_df <- dwi_over_55 %>% select(participant_id, homeint_storyrecall_d, 
                                      bp_dia_mean_cardio, bp_sys_mean_cardio,
                                      additional_hads_anxiety, additional_hads_depression) %>%
  left_join(amygdala_df, .)


summary(lm(kfa_lh_amygdala ~ homeint_storyrecall_d))
amygdala_lh_kfa_odi <- lm(kfa_lh_amygdala ~ odi_lh_amygdala, data = amygdala_df)
summary(amygdala_lh_kfa_odi)
ggplot(amygdala_df, aes(kfa_lh_amygdala, odi_lh_amygdala)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x) +
  geom_richtext(aes(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                    label = paste0(
                      "p < 0.001",
                      ", adj-R<sup>2</sup> = ", signif(summary(amygdala_lh_kfa_odi)$adj.r.squared, 2),
                      ", \u03B2 = ", signif(summary(amygdala_lh_kfa_odi)$coefficients[2,1], 2))), 
                fill = NA, label.color = NA) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
