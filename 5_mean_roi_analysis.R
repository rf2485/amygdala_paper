source("1_data_preparation.R")
library(gtsummary)
library(ggeffects)
library(ggtext)

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

subcort_labels <- list(
  Left.Thalamus ~ "Left Thalamus",
  Right.Thalamus ~ "Right Thalamus",
  Left.Caudate ~ "Left Caudate",
  Right.Caudate ~ "Right Caudate",
  Left.Putamen ~ "Left Putamen",
  Right.Putamen ~ "Right Putamen",
  Left.Pallidum ~ "Left Pallidum",
  Right.Pallidum ~ "Right Pallidum",
  Left.Hippocampus ~ "Left Hippocampus",
  Right.Hippocampus ~ "Right Hippocampus",
  Left.Amygdala ~ "Left Amygdala",
  Right.Amygdala ~ "Right Amygdala",
  Left.Accumbens.area ~ "Left Accumbens Area",
  Right.Accumbens.area ~ "Right Accumbens Area",
  Left.VentralDC ~ "Left Ventral Diencephalon",
  Right.VentralDC ~ "Right Ventral Diencephalon"
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

my_ES_test <- function(data, variable, by, ...) {
  rstatix::cohens_d(data, as.formula(glue::glue("{variable} ~ {by}")))$effsize
}
my_cramer_v <- function(data, variable, by, ...) {
  table(data[[variable]], data[[by]]) %>%
    rstatix::cramer_v()
}

#extract SCD status and demographics for participants
scd_status <- dwi_over_55 %>% select(participant_id, SCD, mt_tr, Income, 
                                     Ethnicity, Sex, age, age_education_completed,
                                     homeint_storyrecall_d, bp_dia_mean_cardio, 
                                     bp_sys_mean_cardio, pulse_mean_cardio,
                                     height_cardio, weight_cardio,
                                     additional_hads_anxiety, additional_hads_depression) %>%
  filter(!participant_id %in% failed_qc) %>%
  mutate(bmi_cardio = weight_cardio / (height_cardio/100)^2) 
scd_status <- set_label(scd_status, bmi_cardio = "BMI (kg/m<sup>2</sup>)")
scd_status$additional_hads_depression[scd_status$additional_hads_depression > 21] <- NA

#demographics table
scd_status %>% select(SCD, age, age_education_completed, Sex, Income, Ethnicity) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})") %>%
  # add_difference(test = list(all_continuous() ~ 'cohens_d')) %>%
  # modify_column_hide(conf.low) %>%
  add_stat(
    fns = list(all_continuous() ~ my_ES_test,
               all_categorical() ~ my_cramer_v)) %>%  add_p() %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**",
                add_stat_1 ~ "**Effect size**",
                label ~ "") %>%
  bold_labels() %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Table 1:</div> <div style='text-align: left'> 
  Demographics of the sample. Group differences in age and the age the participant completed their education were assessed 
                 with Wilcoxon rank sum tests and Cohen's D. Group differences in sex and income 
                 were assessed with Pearson’s Chi-squared test and Cramer's V. Group differences 
                 in ethnicity were assessed with Fisher’s exact test and Cramer's V, because some 
                 ethnicity categories had counts less than 5.  
                 SCD: Subjective Cognitive Decline. SD: Standard Deviation.</div>") %>%
  modify_footnote(add_stat_1 ~ "Cohen's D; Cramer's V")

#cog and psych table
scd_status %>% select(SCD, homeint_storyrecall_d, additional_hads_anxiety, additional_hads_depression) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})") %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "",
                estimate ~ "**Effect Size**") %>%
  bold_labels() %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Table 2:</div> <div style='text-align: left'> 
                 Group differences in memory, anxiety, and depression scores. HADS: Hospital Anxiety and Depression
                 Scale. SCD: Subjective Cognitive Decline. SD: Standard Deviation.</div>")

#cardiovascular health table
scd_status %>% select(SCD, bp_sys_mean_cardio, bp_dia_mean_cardio, pulse_mean_cardio, bmi_cardio) %>%
  filter(!is.na(bmi_cardio)) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})") %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "",
                estimate ~ "**Effect Size**") %>%
  bold_labels() %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Table 3:</div> <div style='text-align: left'> 
                 Group differences in cardiovascular health measures. 
                 SCD: Subjective Cognitive Decline. BMI: Body Mass Index. SD: Standard Deviation.</div>") %>%
  as_gt() %>%
  gt::fmt_markdown(columns = vars(label))

#import aseg stats table (subcortical volumes)
aseg <- read_tsv("freesurfer/asegtable.tsv") %>%
  rename(participant_id=`Measure:volume`)
names(aseg) <- make.names(names(aseg))
aparc2meas_files <- list.files(path = "freesurfer", 
                               pattern = "aparc.*aseg.*\\.*tsv", 
                               full.names = T)

for (i in 1:length(aparc2meas_files)) {
  assign(gsub(".tsv", "", 
              gsub("freesurfer.aparc.", "", make.names(aparc2meas_files[i]))), 
         read.delim(aparc2meas_files[i]))
}

# #import JHU stats table (WM volumes)
# jhu_volume <- read_tsv("freesurfer/jhu_volume.tsv") %>%
#   rename(any_of(jhu_lookup)) %>%
#   select(!starts_with("Seg00"))
volumes <- left_join(scd_status, aseg) %>% #join volumes with SCD status
  # left_join(., jhu_volume) %>%
  mutate(across(c(Left.Lateral.Ventricle:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) %>% #normalize by intracranial volume
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename(Left.Cerebral.White.Matter	= lhCerebralWhiteMatterVol, 
         Right.Cerebral.White.Matter	= rhCerebralWhiteMatterVol) %>%
  select(!c((BrainSegVol:CortexVol), (CerebralWhiteMatterVol:EstimatedTotalIntraCranialVol)))
subcort_gm_volumes_table <- volumes %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Supplemental Table 1:</div> <div style='text-align: left'> 
                 Group differences in subcortical regional volumes, normalized by estimated total intracranial volume, and mean NDI. 
                 NDI: Neurite Density Index. SCD: Subjective Cognitive Decline. SD: Standard Deviation.</div>")

fit_NDI <- aseg2fit_NDI %>%
  # left_join(., lh_AD_sig2fit_NDI) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2fit_NDI) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  rename(participant_id = Measure.mean) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_NDI_table <- fit_NDI %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

tbl_merge(list(subcort_gm_volumes_table, subcort_gm_NDI_table), 
          tab_spanner = c("**Volume**", "**NDI**"))

fit_FWF <- aseg2fit_FWF %>%
  # left_join(., lh_AD_sig2fit_FWF) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2fit_FWF) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  rename(participant_id = Measure.mean) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_FWF_table <- fit_FWF %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Table 4:</div> <div style='text-align: left'> 
                 Statistically significant group differences in subcortical regional mean FWF and ODI. 
                 FWF: Free Water Fraction. ODI: Orientation Dispersion Index. 
                 SCD: Subjective Cognitive Decline. SD: Standard Deviation.</div>")
subcort_gm_FWF_table_all <- fit_FWF %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% bold_p() %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Supplemental Table 2:</div> <div style='text-align: left'> 
                 All group differences in subcortical regional mean FWF and ODI. 
                 FWF: Free Water Fraction. ODI: Orientation Dispersion Index. 
                 SCD: Subjective Cognitive Decline. SD: Standard Deviation.</div>")

fit_ODI <- aseg2fit_ODI %>%
  # left_join(., lh_AD_sig2fit_ODI) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2fit_ODI) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  rename(participant_id = Measure.mean) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_ODI_table <- fit_ODI %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
subcort_gm_ODI_table_all <- fit_ODI %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% bold_p () %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

tbl_merge(list(subcort_gm_FWF_table, subcort_gm_ODI_table), 
          tab_spanner = c("**FWF**", "**ODI**"))
tbl_merge(list(subcort_gm_FWF_table_all, subcort_gm_ODI_table_all), 
          tab_spanner = c("**FWF**", "**ODI**"))

dti_fa <- aseg2dti_fa %>%
  # left_join(., lh_AD_sig2dti_fa) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2dti_fa) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_fa) %>%
  rename(any_of(jhu_lookup)) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_fa_table_all <- dti_fa %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% bold_p() %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Supplemental Table 3:</div> <div style='text-align: left'> 
                 All group differences in subcortical regional mean FA and MD. 
                 FA: Fractional Anisotropy. MD: Mean Diffusivity. 
                 SCD: Subjective Cognitive Decline. SD: Standard Deviation.</div>")

dti_md <- aseg2dti_md %>%
  # left_join(., lh_AD_sig2dti_md) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2dti_md) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_md) %>%
  rename(any_of(jhu_lookup)) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_md_table_all <- dti_md %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% bold_p() %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

tbl_merge(list(subcort_gm_fa_table_all, subcort_gm_md_table_all), 
          tab_spanner = c("**FA**", "**MD**"))

dki_kfa <- aseg2dki_kfa %>%
  # left_join(., lh_AD_sig2dki_kfa) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2dki_kfa) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_kfa) %>%
  rename(any_of(jhu_lookup)) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_kfa_table <- dki_kfa %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Table 5:</div> <div style='text-align: left'> 
                 Statistically significant group differences in subcortical regional mean KFA and MD. 
                 KFA: Kurtosis Fractional Anisotropy. MK: Mean Kurtosis. 
                 SCD: Subjective Cognitive Decline. SD: Standard Deviation.</div>")
subcort_gm_kfa_table_all <- dki_kfa %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% bold_p() %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Supplemental Table 4:</div> <div style='text-align: left'> 
                 All group differences in subcortical regional mean KFA and MD. 
                 KFA: Kurtosis Fractional Anisotropy. MK: Mean Kurtosis. 
                 SCD: Subjective Cognitive Decline. SD: Standard Deviation.</div>")

dki_mk <- aseg2dki_mk %>%
  # left_join(., lh_AD_sig2dki_mk) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2dki_mk) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_mk) %>%
  rename(any_of(jhu_lookup)) %>%
  #remove biologically implausible values (Veraat, Hecke, and Sijbers 2011)
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  mutate(across(where(is.double), ~ replace(., . > 4, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
subcort_gm_mk_table <- dki_mk %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
subcort_gm_mk_table_all <- dki_mk %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% bold_p() %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

tbl_merge(list(subcort_gm_kfa_table, subcort_gm_mk_table), 
          tab_spanner = c("**KFA**", "**MK**"))
tbl_merge(list(subcort_gm_kfa_table_all, subcort_gm_mk_table_all), 
          tab_spanner = c("**KFA**", "**MK**"))

mtr_wm <- read_tsv("freesurfer/mtr_wm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_wm=Seg0001)
mtr_gm <- read_tsv("freesurfer/mtr_gm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_gm=Seg0001)
mtr_wm_gm_ratio <- scd_status %>%
  left_join(., mtr_wm) %>%
  left_join(., mtr_gm)
mtr_wm_gm_ratio$mtr_wm_gm_ratio <- mtr_wm_gm_ratio$mtr_wm / mtr_wm_gm_ratio$mtr_gm
# mtr_wm_gm_ratio$coil <- as.factor(mtr_wm_gm_ratio$coil)
# mtr_wm_gm_ratio$mt_tr <- as.factor(mtr_wm_gm_ratio$mt_tr)
# mtr_wm_gm_ratio$mt_tr <- factor(mtr_wm_gm_ratio$mt_tr,
#                                 levels = c(30, 50),
#                                 labels = c("TR=30ms", "TR=50ms"))
mtr_wm_gm_ratio %>% select(mt_tr, mtr_wm_gm_ratio) %>% 
  tbl_summary(by = mt_tr, statistic = all_continuous() ~ "{mean} ({sd})",
              label = mtr_wm_gm_ratio ~ "MTR white matter-gray matter ratio") %>% 
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Table 6:</div> <div style='text-align: left'> 
                 Differences in MTR gray matter-white matter ratio between images with TR=30ms and TR=50ms, 
                 collapsed across SCD and controls.
                 MTR: Magnetization Transfer Ratio. TR: repetition time, an image acquisition parameter. 
                 SCD: Subjective Cognitive Decline. SD: Standard Deviation.</div>")


mtr_tr30 <- aseg2mtr %>%
  # left_join(., lh_AD_sig2mtr) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2mtr) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2mtr) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  # mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
  filter(mt_tr=="TR=30ms")
subcort_gm_mtr_tr30_table_all <- mtr_tr30 %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p() %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Supplemental Table 5:</div> <div style='text-align: left'> 
                 All group differences in subcortical regional mean MTR for images with TR=30ms and TR=50ms.
                 MTR: Magnetization Transfer Ratio. TR: repetition time, an image acquisition parameter. 
                 SCD: Subjective Cognitive Decline. SD: Standard Deviation.</div>")

mtr_tr50 <- aseg2mtr %>%
  # left_join(., lh_AD_sig2mtr) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2mtr) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2mtr) %>%
  rename(any_of(jhu_lookup)) %>%
  left_join(scd_status, .) %>%
  # mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
  filter(mt_tr=="TR=50ms")
subcort_gm_mtr_tr50_table_all <- mtr_tr50 %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              # missing_text = "Excluded Outliers",
              label = subcort_labels
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p() %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")  

tbl_merge(list(subcort_gm_mtr_tr30_table_all, subcort_gm_mtr_tr50_table_all), 
          tab_spanner = c("**MTR TR=30ms**", "**MTR TR=50ms**"))

left_putamen <- volumes %>% dplyr::select(participant_id:bmi_cardio, Left.Putamen) %>%
  rename(volume = "Left.Putamen")
left_putamen <- fit_NDI %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(fit_NDI = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- fit_FWF %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(fit_FWF = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- fit_ODI %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(fit_ODI = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- dti_fa %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(dti_fa = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- dti_md %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(dti_md = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- dki_kfa %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(dki_kfa = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- dki_mk %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(dki_mk = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- mtr_tr30 %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(mtr_tr30 = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- mtr_tr50 %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(mtr_tr50 = "Left.Putamen") %>% full_join(left_putamen, .)

left_putamen <- set_label(left_putamen,
                          fit_NDI = "NDI",
                          fit_FWF = "FWF",
                          fit_ODI = "ODI",
                          dti_fa = "FA",
                          dti_md = "MD",
                          dki_kfa = "KFA",
                          dki_mk = "MK",
                          mtr_tr30 = "MTR TR=30ms",
                          mtr_tr50 = "MTR TR=50ms"
)


left_putamen_FWF_by_volume_all <- 
  lm(fit_FWF ~ volume * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_volume_all) #age significant
left_putamen_FWF_by_NDI_all <- 
  lm(fit_FWF ~ fit_NDI * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_NDI_all) #age significant
left_putamen_FWF_by_ODI_all <- 
  lm(fit_FWF ~ fit_ODI * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_ODI_all) #age significant
left_putamen_FWF_by_FA_all <- 
  lm(fit_FWF ~ dti_fa * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_FA_all) #age significant
left_putamen_FWF_by_MD_all <- 
  lm(fit_FWF ~ dti_md * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_MD_all) #none significant
left_putamen_FWF_by_KFA_all <- 
  lm(fit_FWF ~ dki_kfa * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_KFA_all) #age and sex significant
left_putamen_FWF_by_MK_all <- 
  lm(fit_FWF ~ dki_mk * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_MK_all) #age significant
left_putamen_FWF_by_mtr_tr30_all <- 
  lm(fit_FWF ~ mtr_tr30 * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_mtr_tr30_all) #age significant
left_putamen_FWF_by_mtr_tr50_all <- 
  lm(fit_FWF ~ mtr_tr50 * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_mtr_tr50_all) #age significant

left_amygdala <- volumes %>% dplyr::select(participant_id:bmi_cardio, Left.Amygdala) %>%
  rename(volume = "Left.Amygdala")
left_amygdala <- fit_NDI %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(fit_NDI = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- fit_FWF %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(fit_FWF = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- fit_ODI %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(fit_ODI = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dti_fa %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(dti_fa = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dti_md %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(dti_md = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_kfa %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(dki_kfa = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_mk %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(dki_mk = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- mtr_tr30 %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(mtr_tr30 = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- mtr_tr50 %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(mtr_tr50 = "Left.Amygdala") %>% full_join(left_amygdala, .)

left_amygdala <- set_label(left_amygdala,
                           fit_NDI = "NDI",
                           fit_FWF = "FWF",
                           fit_ODI = "ODI",
                           dti_fa = "FA",
                           dti_md = "MD",
                           dki_kfa = "KFA",
                           dki_mk = "MK",
                           mtr_tr30 = "MTR TR=30ms",
                           mtr_tr50 = "MTR TR=50ms"
)

left_amygdala_KFA_by_volume_all <- 
  lm(dki_kfa ~ volume * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_amygdala)
anova(left_amygdala_KFA_by_volume_all) #Income significant
left_amygdala_KFA_by_NDI_all <- 
  lm(dki_kfa ~ fit_NDI * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_amygdala)
anova(left_amygdala_KFA_by_NDI_all) #Income significant
left_amygdala_KFA_by_FWF_all <- 
  lm(dki_kfa ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_amygdala)
anova(left_amygdala_KFA_by_FWF_all) #income and education significant
left_amygdala_KFA_by_ODI_all <- 
  lm(dki_kfa ~ fit_ODI * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_amygdala)
anova(left_amygdala_KFA_by_ODI_all) #income significant
left_amygdala_KFA_by_MK_all <- 
  lm(dki_kfa ~ dki_mk * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_amygdala)
anova(left_amygdala_KFA_by_MK_all) #income significant
left_amygdala_KFA_by_FA_all <- 
  lm(dki_kfa ~ dti_fa * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_amygdala)
anova(left_amygdala_KFA_by_FA_all) #income significant
left_amygdala_KFA_by_MD_all <- 
  lm(dki_kfa ~ dti_md * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_amygdala)
anova(left_amygdala_KFA_by_MD_all) #income significant
left_amygdala_KFA_by_mtr_tr30_all <- 
  lm(dki_kfa ~ mtr_tr30 * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_amygdala)
anova(left_amygdala_KFA_by_mtr_tr30_all) #none significant
left_amygdala_KFA_by_mtr_tr50_all <- 
  lm(dki_kfa ~ mtr_tr50 * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_amygdala)
anova(left_amygdala_KFA_by_mtr_tr50_all) #none significant

right_amygdala <- volumes %>% dplyr::select(participant_id:bmi_cardio, Right.Amygdala) %>%
  rename(volume = "Right.Amygdala")
right_amygdala <- fit_NDI %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(fit_NDI = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- fit_FWF %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(fit_FWF = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- fit_ODI %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(fit_ODI = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dti_fa %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(dti_fa = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dti_md %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(dti_md = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_kfa %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(dki_kfa = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_mk %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(dki_mk = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- mtr_tr30 %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(mtr_tr30 = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- mtr_tr50 %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(mtr_tr50 = "Right.Amygdala") %>% full_join(right_amygdala, .)

right_amygdala <- set_label(right_amygdala,
                            fit_NDI = "NDI",
                            fit_FWF = "FWF",
                            fit_ODI = "ODI",
                            dti_fa = "FA",
                            dti_md = "MD",
                            dki_kfa = "KFA",
                            dki_mk = "MK",
                            mtr_tr30 = "MTR TR=30ms",
                            mtr_tr50 = "MTR TR=50ms"
)
right_amygdala_KFA_by_volume_all <- 
  lm(dki_kfa ~ volume * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     right_amygdala)
anova(right_amygdala_KFA_by_volume_all) #sex significant
right_amygdala_KFA_by_NDI_all <- 
  lm(dki_kfa ~ fit_NDI * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     right_amygdala)
anova(right_amygdala_KFA_by_NDI_all) #age and sex significant
right_amygdala_KFA_by_FWF_all <- 
  lm(dki_kfa ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     right_amygdala)
anova(right_amygdala_KFA_by_FWF_all) #sex significant
right_amygdala_KFA_by_ODI_all <- 
  lm(dki_kfa ~ fit_ODI * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     right_amygdala)
anova(right_amygdala_KFA_by_ODI_all) #sex significant
right_amygdala_KFA_by_MK_all <- 
  lm(dki_kfa ~ dki_mk * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     right_amygdala)
anova(right_amygdala_KFA_by_MK_all) #age and sex significant
right_amygdala_KFA_by_FA_all <- 
  lm(dki_kfa ~ dti_fa * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     right_amygdala)
anova(right_amygdala_KFA_by_FA_all) #none significant
right_amygdala_KFA_by_MD_all <- 
  lm(dki_kfa ~ dti_md * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     right_amygdala)
anova(right_amygdala_KFA_by_MD_all) #sex significant
right_amygdala_KFA_by_mtr_tr30_all <- 
  lm(dki_kfa ~ mtr_tr30 * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     right_amygdala)
anova(right_amygdala_KFA_by_mtr_tr30_all) #none significant
right_amygdala_KFA_by_mtr_tr50_all <- 
  lm(dki_kfa ~ mtr_tr50 * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     right_amygdala)
anova(right_amygdala_KFA_by_mtr_tr50_all) #none significant

left_accumbens_area <- volumes %>% dplyr::select(participant_id:bmi_cardio, Left.Accumbens.area) %>%
  rename(volume = "Left.Accumbens.area")
left_accumbens_area <- fit_NDI %>% dplyr::select(participant_id, Left.Accumbens.area) %>%
  rename(fit_NDI = "Left.Accumbens.area") %>% full_join(left_accumbens_area, .)
left_accumbens_area <- fit_FWF %>% dplyr::select(participant_id, Left.Accumbens.area) %>%
  rename(fit_FWF = "Left.Accumbens.area") %>% full_join(left_accumbens_area, .)
left_accumbens_area <- fit_ODI %>% dplyr::select(participant_id, Left.Accumbens.area) %>%
  rename(fit_ODI = "Left.Accumbens.area") %>% full_join(left_accumbens_area, .)
left_accumbens_area <- dti_fa %>% dplyr::select(participant_id, Left.Accumbens.area) %>%
  rename(dti_fa = "Left.Accumbens.area") %>% full_join(left_accumbens_area, .)
left_accumbens_area <- dti_md %>% dplyr::select(participant_id, Left.Accumbens.area) %>%
  rename(dti_md = "Left.Accumbens.area") %>% full_join(left_accumbens_area, .)
left_accumbens_area <- dki_kfa %>% dplyr::select(participant_id, Left.Accumbens.area) %>%
  rename(dki_kfa = "Left.Accumbens.area") %>% full_join(left_accumbens_area, .)
left_accumbens_area <- dki_mk %>% dplyr::select(participant_id, Left.Accumbens.area) %>%
  rename(dki_mk = "Left.Accumbens.area") %>% full_join(left_accumbens_area, .)
left_accumbens_area <- mtr_tr30 %>% dplyr::select(participant_id, Left.Accumbens.area) %>%
  rename(mtr_tr30 = "Left.Accumbens.area") %>% full_join(left_accumbens_area, .)
left_accumbens_area <- mtr_tr50 %>% dplyr::select(participant_id, Left.Accumbens.area) %>%
  rename(mtr_tr50 = "Left.Accumbens.area") %>% full_join(left_accumbens_area, .)

left_accumbens_area <- set_label(left_accumbens_area,
                           fit_NDI = "NDI",
                           fit_FWF = "FWF",
                           fit_ODI = "ODI",
                           dti_fa = "FA",
                           dti_md = "MD",
                           dki_kfa = "KFA",
                           dki_mk = "MK",
                           mtr_tr30 = "MTR TR=30ms",
                           mtr_tr50 = "MTR TR=50ms"
)

## old code
wm_volumes_table <- volumes %>% 
  select(c(SCD, CC_Posterior:Right.Cerebral.White.Matter, 
           Left.Cerebellum.White.Matter, Right.Cerebellum.White.Matter)) %>%
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
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  # mutate(across(where(is.double), remove_outliers)) %>%
  rename_with(., ~ gsub("_", ".", .x, fixed = T), ends_with("thickness")) %>%
  rename_with(., ~ paste0("ctx.", .x, recycle0 = T), ends_with("thickness")) %>%
  rename_with(., ~ gsub(".thickness", "", .x, fixed = T)) #%>%
# left_join(., lh_AD_sig_thickness) %>%
# left_join(., rh_AD_sig_thickness)
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


AD_sig2meas_files <- list.files(path = "freesurfer", 
                                pattern = "?h_AD_sig2.*\\.*tsv", 
                                full.names = T)

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
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  # left_join(., lh_AD_sig2dti_rd) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2dti_rd) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_rd) %>%
  rename(any_of(jhu_lookup)) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
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
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  # left_join(., lh_AD_sig2dti_ad) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2dti_ad) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dti_ad) %>%
  rename(any_of(jhu_lookup)) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
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
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              missing_text = "Excluded Outliers") %>%
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
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  # left_join(., lh_AD_sig2dki_rk) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2dki_rk) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_rk) %>%
  rename(any_of(jhu_lookup)) %>%
  #remove biologically implausible values (Veraat, Hecke, and Sijbers 2011)
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  mutate(across(where(is.double), ~ replace(., . > 4, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
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
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  # left_join(., lh_AD_sig2dki_ak) %>%
  # rename(lh_AD_signature = AD_signature) %>%
  # left_join(., rh_AD_sig2dki_ak) %>%
  # rename(rh_AD_signature = AD_signature) %>%
  left_join(., jhu2dki_ak) %>%
  rename(any_of(jhu_lookup)) %>%
  #remove biologically implausible values (Veraat, Hecke, and Sijbers 2011)
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  mutate(across(where(is.double), ~ replace(., . > 4, NA))) %>%
  left_join(scd_status, .) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
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
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_Da_table <- smi_matlab_Da %>%
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_DePar_table <- smi_matlab_DePar %>%
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_DePerp_table <- smi_matlab_DePerp %>%
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  select(!starts_with("Seg00")) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_f_table <- smi_matlab_f %>%
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  select(!starts_with("Seg00")) #%>%
# mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
wm_p2_table <- smi_matlab_p2 %>%
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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

#mtr stats for each TR
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
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  # mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
  filter(mt_tr=="TR=30ms")
wm_g_ratio_tr30_table <- g_ratio_tr30 %>%
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  # mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
  filter(mt_tr=="TR=50ms")
wm_g_ratio_tr50_table <- g_ratio_tr50 %>%
  select(c(SCD, middle_cerebellar_peduncle:tapetum_L)) %>%
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

mtr <- aseg2mtr %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)


right_amygdala <- fit_NDI %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(fit_NDI = "Right.Amygdala")
right_amygdala <- fit_ODI %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(fit_ODI = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_kfa %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dki_kfa = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_mk %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dki_mk = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_ak %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dki_ak = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_rk %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dki_rk = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- mtr %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(mtr = "Right.Amygdala") %>% full_join(right_amygdala, .) %>%
  filter(!participant_id %in% failed_qc) %>%
  filter(mt_tr=="TR=30ms")

right_amygdala <- set_label(right_amygdala,
                            fit_NDI = "NODDI Neurite Density (ND)",
                            fit_ODI = "NODDI Orientation Dispersion (OD)",
                            dki_kfa = "DKI Kurtosis Fractional Anisotropy (KFA)",
                            dki_mk = "DKI Mean Kurtosis (MK)",
                            dki_ak = "DKI Axial Kurtosis (AK)",
                            dki_rk = "DKI Radial Kurtosis (RK)",
                            mtr = "Magnetization Transfer Ratio (MTR)"
)

#group means of mtr tr30 subset
left_amygdala %>%
  select(SCD, Sex, Income, Ethnicity, age, age_education_completed) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})") %>% add_p()

left_amygdala_table <- left_amygdala %>% 
  select(SCD, fit_NDI, fit_ODI, dki_kfa, dki_mk, dki_ak, dki_rk, mtr) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Metric**",
                estimate ~ "**Effect Size**")

right_amygdala_table <- right_amygdala %>% 
  select(SCD, fit_NDI, fit_ODI, dki_kfa, dki_mk, dki_ak, dki_rk, mtr) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>%  bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Metric**",
                estimate ~ "**Effect Size**")

tbl_merge(list(left_amygdala_table, right_amygdala_table),
          tab_spanner = c("**Left Amygdala**", "**Right Amygdala**"))

#so education can be compared to other covariates:
right_amygdala_na_omit <- mutate(right_amygdala,
                                 across(where(is.factor),~fct_explicit_na(.,'Unknown'))) %>%
  na.omit(.)

#scatterplots: neurodegen and demyelination/protein aggregation
#first assess if each covariate improves the model. Discard if ANOVA is not significant.
KFA_by_NDI <- lm(dki_kfa ~ fit_NDI, right_amygdala_na_omit)
KFA_by_NDI_scd <- lm(dki_kfa ~ fit_NDI + SCD, right_amygdala_na_omit)
anova(KFA_by_NDI, KFA_by_NDI_scd)
KFA_by_NDI_age <- lm(dki_kfa ~ fit_NDI + age, right_amygdala_na_omit)
anova(KFA_by_NDI, KFA_by_NDI_age) #improvement
KFA_by_NDI_age_sex <- lm(dki_kfa ~ fit_NDI + age + Sex, right_amygdala_na_omit)
anova(KFA_by_NDI_age, KFA_by_NDI_age_sex) #improvement
KFA_by_NDI_age_sex_education <- lm(dki_kfa ~ fit_NDI + age + Sex + age_education_completed, right_amygdala_na_omit)
anova(KFA_by_NDI_age_sex, KFA_by_NDI_age_sex_education)
KFA_by_NDI_age_sex_income <- lm(dki_kfa ~ fit_NDI + age + Sex + Income, right_amygdala_na_omit)
anova(KFA_by_NDI_age_sex, KFA_by_NDI_age_sex_income)
KFA_by_NDI_age_sex_ethnicity <- lm(dki_kfa ~ fit_NDI + age + Sex + Ethnicity, right_amygdala_na_omit)
anova(KFA_by_NDI_age_sex, KFA_by_NDI_age_sex_ethnicity)

KFA_by_NDI_scd_age_sex_education_income_ethnicity <- 
  lm(dki_kfa ~ fit_NDI * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     right_amygdala_na_omit)
anova(KFA_by_NDI_scd_age_sex_education_income_ethnicity)

#image width 570 height 481
#only education has NAs and we are not using it, so switch back to full dataset
KFA_by_NDI_age_sex <- lm(dki_kfa ~ fit_NDI + age + Sex, right_amygdala)
adj_r_squared <- summary(KFA_by_NDI_age_sex)$adj.r.squared
model_p_value <- pf(summary(KFA_by_NDI_age_sex)$fstatistic[1], 
                    summary(KFA_by_NDI_age_sex)$fstatistic[2], 
                    summary(KFA_by_NDI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(KFA_by_NDI_age_sex, c( "fit_NDI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "darkorange") + 
  labs(
    title = "A. Mean Right Amygdala KFA and ND, Corrected by Age and Sex",
    subtitle = paste0("** p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5, face = "bold"))

MK_by_NDI_age_sex <- lm(dki_mk ~ fit_NDI + age + Sex, right_amygdala)
adj_r_squared <- summary(MK_by_NDI_age_sex)$adj.r.squared
model_p_value <- pf(summary(MK_by_NDI_age_sex)$fstatistic[1], 
                    summary(MK_by_NDI_age_sex)$fstatistic[2], 
                    summary(MK_by_NDI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(MK_by_NDI_age_sex, c( "fit_NDI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "darkorange") + 
  labs(
    title = "B. Mean Right Amygdala MK and ND, Corrected by Age and Sex",
    subtitle = paste0("*** p < 0.001",
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5, face = "bold"))

RK_by_NDI_age_sex <- lm(dki_rk ~ fit_NDI + age + Sex, right_amygdala)
adj_r_squared <- summary(RK_by_NDI_age_sex)$adj.r.squared
model_p_value <- pf(summary(RK_by_NDI_age_sex)$fstatistic[1], 
                    summary(RK_by_NDI_age_sex)$fstatistic[2], 
                    summary(RK_by_NDI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(RK_by_NDI_age_sex, c( "fit_NDI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "darkorange") + 
  labs(
    title = "C. Mean Right Amygdala RK and ND, Corrected by Age and Sex",
    subtitle = paste0("*** p < 0.001",
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5, face = "bold"))

KFA_by_ODI_age_sex <- lm(dki_kfa ~ fit_ODI + age + Sex, right_amygdala)
adj_r_squared <- summary(KFA_by_ODI_age_sex)$adj.r.squared
model_p_value <- pf(summary(KFA_by_ODI_age_sex)$fstatistic[1], 
                    summary(KFA_by_ODI_age_sex)$fstatistic[2], 
                    summary(KFA_by_ODI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(KFA_by_ODI_age_sex, c( "fit_ODI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "purple") + 
  labs(
    title = "D. Mean Right Amygdala KFA and OD, Corrected by Age and Sex",
    subtitle = paste0("p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

MK_by_ODI_age_sex <- lm(dki_mk ~ fit_ODI + age + Sex, right_amygdala)
adj_r_squared <- summary(MK_by_ODI_age_sex)$adj.r.squared
model_p_value <- pf(summary(MK_by_ODI_age_sex)$fstatistic[1], 
                    summary(MK_by_ODI_age_sex)$fstatistic[2], 
                    summary(MK_by_ODI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(MK_by_ODI_age_sex, c( "fit_ODI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "purple") + 
  labs(
    title = "E. Mean Right Amygdala MK and OD, Corrected by Age and Sex",
    subtitle = paste0("\\* p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5, face = "bold"))

RK_by_ODI_age_sex <- lm(dki_rk ~ fit_ODI + age + Sex, right_amygdala)
adj_r_squared <- summary(RK_by_ODI_age_sex)$adj.r.squared
model_p_value <- pf(summary(RK_by_ODI_age_sex)$fstatistic[1], 
                    summary(RK_by_ODI_age_sex)$fstatistic[2], 
                    summary(RK_by_ODI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(RK_by_ODI_age_sex, c( "fit_ODI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "purple") + 
  labs(
    title = "F. Mean Right Amygdala RK and OD, Corrected by Age and Sex",
    subtitle = paste0("p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

KFA_by_MTR_age_sex <- lm(dki_kfa ~ mtr + age + Sex, right_amygdala)
adj_r_squared <- summary(KFA_by_MTR_age_sex)$adj.r.squared
model_p_value <- pf(summary(KFA_by_MTR_age_sex)$fstatistic[1], 
                    summary(KFA_by_MTR_age_sex)$fstatistic[2], 
                    summary(KFA_by_MTR_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(KFA_by_MTR_age_sex, c( "mtr [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "turquoise3") + 
  labs(
    title = "G. Mean Right Amygdala KFA and MTR, Corrected by Age and Sex",
    subtitle = paste0("p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

MK_by_MTR_age_sex <- lm(dki_mk ~ mtr + age + Sex, right_amygdala)
adj_r_squared <- summary(MK_by_MTR_age_sex)$adj.r.squared
model_p_value <- pf(summary(MK_by_MTR_age_sex)$fstatistic[1], 
                    summary(MK_by_MTR_age_sex)$fstatistic[2], 
                    summary(MK_by_MTR_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(MK_by_MTR_age_sex, c( "mtr [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "turquoise3") + 
  labs(
    title = "H. Mean Right Amygdala MK and MTR, Corrected by Age and Sex",
    subtitle = paste0("p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))


RK_by_MTR_age_sex <- lm(dki_rk ~ mtr + age + Sex, right_amygdala)
adj_r_squared <- summary(RK_by_MTR_age_sex)$adj.r.squared
model_p_value <- pf(summary(RK_by_MTR_age_sex)$fstatistic[1], 
                    summary(RK_by_MTR_age_sex)$fstatistic[2], 
                    summary(RK_by_MTR_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(RK_by_MTR_age_sex, c( "mtr [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "turquoise3") + 
  labs(
    title = "I. Mean Right Amygdala RK and MTR, Corrected by Age and Sex",
    subtitle = paste0("p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))




library(ggpmisc)
amygdala_df <- fit_NDI %>% select(participant_id, SCD, Left.Amygdala, Right.Amygdala) %>%
  rename(kfa_lh_amygdala = Left.Amygdala, kfa_rh_amygdala = Right.Amygdala)
amygdala_df <- dwi_over_55 %>% select(participant_id, homeint_storyrecall_d, 
                                      bp_dia_mean_cardio, bp_sys_mean_cardio,
                                      additional_hads_anxiety, additional_hads_depression) %>%
  left_join(amygdala_df, .) %>% mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
filter()


summary(lm(kfa_rh_amygdala ~ homeint_storyrecall_d, amygdala_df))
summary(lm(kfa_rh_amygdala ~ additional_hads_anxiety, amygdala_df))
summary(lm(kfa_rh_amygdala ~ additional_hads_depression, amygdala_df))
summary(lm(kfa_rh_amygdala ~ bp_sys_mean_cardio, amygdala_df))


amygdala_lh_kfa_odi <- lm(kfa_lh_amygdala ~ additional_hads_depression, data = amygdala_df)
summary(amygdala_lh_kfa_odi)
ggplot(amygdala_df, aes(kfa_lh_amygdala, additional_hads_depression)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x) +
  geom_richtext(aes(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                    label = paste0(
                      "p = ", signif(summary(amygdala_lh_kfa_odi)$coeff)
                      ", adj-R<sup>2</sup> = ", signif(summary(amygdala_lh_kfa_odi)$adj.r.squared, 2),
                      ", \u03B2 = ", signif(summary(amygdala_lh_kfa_odi)$coefficients[2,1], 2))), 
                fill = NA, label.color = NA) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
