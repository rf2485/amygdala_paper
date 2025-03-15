source("1_data_preparation_gm.R")
library(gtsummary)
library(ggeffects)
library(ggtext)

######### subcort GM #################
failed_qc <- c('sub-CC510255', #SCD abnormality in left temporal pole
               'sub-CC510438', #CTL abnormality in left frontal lobe
               'sub-CC620821', #SCD segmentation errors from large ventricles
               'sub-CC621011', #CTL segmentation errors from large ventricles
               'sub-CC621080', #SCD segmentation errors
               'sub-CC710551', #CTL motion artifacts in DWI
               'sub-CC711027', #SCD severe motion artifacts in T1
               'sub-CC721434'  #CTL segmentation errors from large ventricles
)

subcort_gm <- c("SCD", "Left.Thalamus", "Right.Thalamus", "Left.Caudate", "Right.Caudate",
                "Left.Putamen", "Right.Putamen", "Left.Pallidum", "Right.Pallidum",
                "Left.Hippocampus", "Right.Hippocampus", "Left.Amygdala", "Right.Amygdala",
                "Left.Accumbens.area", "Right.Accumbens.area", "Left.VentralDC", "Right.VentralDC")

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

my_ES_test <- function(data, variable, by, ...) {
  rstatix::cohens_d(data, as.formula(glue::glue("{variable} ~ {by}")))$effsize
}
my_cramer_v <- function(data, variable, by, ...) {
  table(data[[variable]], data[[by]]) %>%
    rstatix::cramer_v()
}

### 3.1	Group differences between SCD and controls ###

#extract SCD status and demographics for participants
scd_status <- dwi_over_55 %>% select(participant_id, SCD, mt_tr, Income, 
                                     Ethnicity, Sex, age, age_education_completed,
                                     homeint_storyrecall_d, bp_dia_mean_cardio, 
                                     bp_sys_mean_cardio, pulse_mean_cardio,
                                     height_cardio, weight_cardio,
                                     additional_hads_anxiety, additional_hads_depression) %>%
  filter(!participant_id %in% failed_qc) %>%
  mutate(bmi_cardio = weight_cardio / (height_cardio/100)^2) %>%
  select(!height_cardio, !weight_cardio)
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
aseg <- read_tsv("freesurfer/asegtable_gm.tsv") %>%
  rename(participant_id=`Measure:volume`)
names(aseg) <- make.names(names(aseg))
aparc2meas_files <- list.files(path = "freesurfer", 
                               pattern = "aparc.*aseg.*_gm\\.*tsv", 
                               full.names = T)

for (i in 1:length(aparc2meas_files)) {
  assign(gsub(".tsv", "", 
              gsub("freesurfer.aparc.", "", make.names(aparc2meas_files[i]))), 
         read.delim(aparc2meas_files[i]))
}

volumes <- left_join(scd_status, aseg) %>% #join volumes with SCD status
  mutate(across(c(Left.Lateral.Ventricle:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) %>% #normalize by intracranial volume
  rename(Left.Cerebral.White.Matter	= lhCerebralWhiteMatterVol, 
         Right.Cerebral.White.Matter	= rhCerebralWhiteMatterVol) %>%
  select(!c((BrainSegVol:CortexVol), (CerebralWhiteMatterVol:EstimatedTotalIntraCranialVol)))

fit_NDI <- aseg2gm_fit_NDI_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) 

fit_FWF <- aseg2gm_fit_FWF_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) 

fit_ODI <- aseg2gm_fit_ODI_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) 

dti_fa <- aseg2dti_fa_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) 

dti_md <- aseg2dti_md_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) 

dti_ad <- aseg2dti_ad_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) 

dti_rd <- aseg2dti_rd_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) 

dki_kfa <- aseg2dki_kfa_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

dki_mk <- aseg2dki_mk_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) 

dki_ak <- aseg2dki_ak_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) 

dki_rk <- aseg2dki_rk_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) 

mtr_wm <- read_tsv("freesurfer/mtr_wm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_wm=Seg0001)
mtr_gm <- read_tsv("freesurfer/mtr_gm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_gm=Seg0001)
mtr_wm_gm_ratio <- scd_status %>%
  left_join(., mtr_wm) %>%
  left_join(., mtr_gm)
mtr_wm_gm_ratio$mtr_wm_gm_ratio <- mtr_wm_gm_ratio$mtr_wm / mtr_wm_gm_ratio$mtr_gm

mtr_tr30 <- aseg2mtr %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  filter(mt_tr=="TR=30ms")

mtr_tr50 <- aseg2mtr %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  filter(mt_tr=="TR=50ms")


#### focus on MTL ####
left_amygdala %>% 
  select(-participant_id) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(dti_fa, dti_md, dki_kfa, dki_mk)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p() %>% bold_p(q=T) %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Left Amygdala </div>")

right_amygdala %>% 
  select(-participant_id) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(dti_fa, dti_md, dki_kfa, dki_mk)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p() %>% bold_p(q=T) %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Right Amygdala </div>")

mtl_imaging <- aseg %>% 
  mutate(across(c(Left.Lateral.Ventricle:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) %>% #normalize by intracranial volume
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename_with(~paste0(., "_volume"), !participant_id)
mtl_imaging <- aseg2dki_ak %>%
  rename(participant_id = Measure.mean) %>%
  #remove biologically implausible values (Veraat, Hecke, and Sijbers 2011)
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  mutate(across(where(is.double), ~ replace(., . > 4, NA))) %>%
  rename_with(~paste0(., "_AK"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2dki_kfa %>%
  rename(participant_id = Measure.mean) %>%
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename_with(~paste0(., "_KFA"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2dki_mk %>%
  rename(participant_id = Measure.mean) %>%
  #remove biologically implausible values (Veraat, Hecke, and Sijbers 2011)
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  mutate(across(where(is.double), ~ replace(., . > 4, NA))) %>%
  rename_with(~paste0(., "_MK"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2dki_mkt %>%
  rename(participant_id = Measure.mean) %>%
  #remove biologically implausible values (Veraat, Hecke, and Sijbers 2011)
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  mutate(across(where(is.double), ~ replace(., . > 4, NA))) %>%
  rename_with(~paste0(., "_MKT"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2dki_rk %>%
  rename(participant_id = Measure.mean) %>%
  #remove biologically implausible values (Veraat, Hecke, and Sijbers 2011)
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  mutate(across(where(is.double), ~ replace(., . > 4, NA))) %>%
  rename_with(~paste0(., "_RK"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2dti_ad %>%
  rename(participant_id = Measure.mean) %>%
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename_with(~paste0(., "_AxD"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2dti_fa %>%
  rename(participant_id = Measure.mean) %>%
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename_with(~paste0(., "_FA"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2dti_md %>%
  rename(participant_id = Measure.mean) %>%
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename_with(~paste0(., "_MD"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2dti_rd %>%
  rename(participant_id = Measure.mean) %>%
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename_with(~paste0(., "_RD"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2fit_FWF %>%
  rename(participant_id = Measure.mean) %>%
  #remove biologically implausible values (Correia et al.)
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  mutate(across(where(is.double), ~ replace(., . > 0.7, NA))) %>%
  rename_with(~paste0(., "_FWF"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2fit_NDI %>%
  rename(participant_id = Measure.mean) %>%
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename_with(~paste0(., "_NDI"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2fit_ODI %>%
  rename(participant_id = Measure.mean) %>%
  #remove biologically implausible values (Correia et al.)
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  mutate(across(where(is.double), ~ replace(., . > 0.8, NA))) %>%
  rename_with(~paste0(., "_ODI"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2mtr %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  filter(mt_tr=="TR=30ms") %>%
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename_with(~paste0(., "_MTR_TR30"), !participant_id) %>%
  left_join(mtl_imaging, .)
mtl_imaging <- aseg2mtr %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  filter(mt_tr=="TR=50ms") %>%
  select(participant_id, Left.Amygdala, Left.Hippocampus, Right.Amygdala, Right.Hippocampus) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename_with(~paste0(., "_MTR_TR50"), !participant_id) %>%
  left_join(mtl_imaging, .)


#across groups
mtl_imaging_across_groups <- scd_status %>%
  select(participant_id, age, age_education_completed, 
         bp_dia_mean_cardio, bp_sys_mean_cardio, pulse_mean_cardio, bmi_cardio,
         homeint_storyrecall_d, additional_hads_anxiety, additional_hads_depression) %>%
  left_join(., mtl_imaging)
mtl_imaging_matrix <- mtl_imaging_across_groups %>% select(Left.Amygdala_volume:ncol(.))
scd_status_matrix <- mtl_imaging_across_groups %>% select(age:additional_hads_depression)
corr_matrix <- psych::corr.test(scd_status_matrix, mtl_imaging_matrix, use = "pairwise.complete.obs", 
                                adjust = "fdr")
corrplot::corrplot(corr_matrix$r, p.mat = corr_matrix$p.adj, method = 'color',
                   sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', pch.cex = 0.7)

#scd only
mtl_imaging_scd <- scd_status %>%
  dplyr::filter(SCD == "SCD") %>%
  select(participant_id, age, age_education_completed, 
         bp_dia_mean_cardio, bp_sys_mean_cardio, pulse_mean_cardio, bmi_cardio,
         homeint_storyrecall_d, additional_hads_anxiety, additional_hads_depression) %>%
  left_join(., mtl_imaging)
mtl_imaging_matrix <- mtl_imaging_scd %>% select(Left.Amygdala_volume:ncol(.))
scd_status_matrix <- mtl_imaging_scd %>% select(age:additional_hads_depression)
corr_matrix <- psych::corr.test(scd_status_matrix, mtl_imaging_matrix, use = "pairwise.complete.obs", 
                                adjust = "fdr")
corrplot::corrplot(corr_matrix$r, p.mat = corr_matrix$p.adj, method = 'color',
                   sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', pch.cex = 0.7)

#control only
mtl_imaging_control <- scd_status %>%
  dplyr::filter(SCD == "Control") %>%
  select(participant_id, age, age_education_completed, 
         bp_dia_mean_cardio, bp_sys_mean_cardio, pulse_mean_cardio, bmi_cardio,
         homeint_storyrecall_d, additional_hads_anxiety, additional_hads_depression) %>%
  left_join(., mtl_imaging)
mtl_imaging_matrix <- mtl_imaging_control %>% select(Left.Amygdala_volume:ncol(.))
scd_status_matrix <- mtl_imaging_control %>% select(age:additional_hads_depression)
corr_matrix <- psych::corr.test(scd_status_matrix, mtl_imaging_matrix, use = "pairwise.complete.obs", 
                                adjust = "fdr")
corrplot::corrplot(corr_matrix$r, p.mat = corr_matrix$p.adj, method = 'color',
                   sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', pch.cex = 0.7)

mtl_imaging <- left_join(scd_status, mtl_imaging)
mtl_imaging %>%
  select(-participant_id) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(Left.Amygdala_MK, Left.Amygdala_RK, Left.Amygdala_KFA, Left.Amygdala_AK)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p() %>% bold_p(q=T) %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

### 3.2	Correlations Between Imaging Metrics ####
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


##models 
left_amygdala_KFA_by_volume_age_sex_income_education <- 
  lm(volume ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_volume_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_volume_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_volume_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_volume_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Amygdala KFA and Volume, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** volume p < 0.001",
      "volume p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_KFA_by_FWF_age_sex_income_education <- 
  lm(fit_FWF ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_FWF_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_FWF_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_FWF_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_FWF_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Amygdala KFA and FWF, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** FWF p < 0.001",
      # "FWF p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_KFA_by_NDI_age_sex_income_education <- 
  lm(fit_NDI ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_NDI_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_NDI_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_NDI_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_NDI_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Left Amygdala KFA and NDI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** model p < 0.001",
      "NDI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_KFA_by_ODI_age_sex_income_education <- 
  lm(fit_ODI ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_ODI_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_ODI_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_ODI_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_ODI_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Left Amygdala KFA and ODI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** ODI p < 0.001",
      "\\* ODI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_KFA_by_FA_age_sex_income_education <- 
  lm(dti_fa ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_FA_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_FA_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_FA_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_FA_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Left Amygdala KFA and FA, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** FA p < 0.001",
      # "FA p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> \\* group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_KFA_by_MD_age_sex_income_education <- 
  lm(dti_md ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_MD_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_MD_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_MD_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_MD_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Amygdala KFA and MD, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** MD p < 0.001",
      # "MD p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_KFA_by_MK_age_sex_income_education <- 
  lm(dki_mk ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_MK_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_MK_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_MK_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_MK_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "E. Mean Left Amygdala KFA and MK, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** MK p < 0.001",
      # "\\* MK p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_KFA_by_mtr_tr30_age_sex_income_education <- 
  lm(mtr_tr30 ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_mtr_tr30_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_mtr_tr30_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_mtr_tr30_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_mtr_tr30_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Left Amygdala KFA and MTR TR=30ms, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** MTR TR=30ms p < 0.001",
      "MTR TR=30ms p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_KFA_by_mtr_tr50_age_sex_income_education <- 
  lm(mtr_tr50 ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_mtr_tr50_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_mtr_tr50_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_mtr_tr50_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_mtr_tr50_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Amygdala KFA and MTR TR=50ms, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** MTR TR=50ms p < 0.001",
      "MTR TR=50ms p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_volume_age_sex_income_education <- 
  lm(volume ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_volume_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_volume_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_volume_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_volume_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Right Amygdala KFA and Volume, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** volume p < 0.001",
      "volume p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_FWF_age_sex_income_education <- 
  lm(fit_FWF ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_FWF_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_FWF_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_FWF_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_FWF_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Right Amygdala KFA and FWF, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** FWF p < 0.001",
      # "FWF p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_NDI_age_sex_income_education <- 
  lm(fit_NDI ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_NDI_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_NDI_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_NDI_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_NDI_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Right Amygdala KFA and NDI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** model p < 0.001",
      "** NDI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_ODI_age_sex_income_education <- 
  lm(fit_ODI ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_ODI_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_ODI_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_ODI_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_ODI_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Right Amygdala KFA and ODI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** ODI p < 0.001",
      "ODI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_FA_age_sex_income_education <- 
  lm(dti_fa ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_FA_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_FA_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_FA_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_FA_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Right Amygdala KFA and FA, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** FA p < 0.001",
      # "FA p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_MD_age_sex_income_education <- 
  lm(dti_md ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_MD_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_MD_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_MD_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_MD_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Right Amygdala KFA and MD, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** MD p < 0.001",
      # "MD p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_MK_age_sex_income_education <- 
  lm(dki_mk ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_MK_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_MK_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_MK_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_MK_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "E. Mean Right Amygdala KFA and MK, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** MK p < 0.001",
      # "\\* MK p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_mtr_tr30_age_sex_income_education <- 
  lm(mtr_tr30 ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_mtr_tr30_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_mtr_tr30_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_mtr_tr30_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_mtr_tr30_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Right Amygdala KFA and MTR TR=30ms, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** MTR TR=30ms p < 0.001",
      "MTR TR=30ms p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_mtr_tr50_age_sex_income_education <- 
  lm(mtr_tr50 ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_mtr_tr50_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_mtr_tr50_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_mtr_tr50_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_mtr_tr50_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Right Amygdala KFA and MTR TR=50ms, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** MTR TR=50ms p < 0.001",
      "MTR TR=50ms p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

### 3.3	Correlations of Significant Imaging Metrics with Memory, Anxiety, and Depression ###
left_amygdala_KFA_by_story_age_sex_income_education <- 
  lm(homeint_storyrecall_d ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_story_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_story_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_story_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_story_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Amygdala KFA and Memory, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** story p < 0.001",
      "story p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_KFA_by_anxiety_age_sex_income_education <- 
  lm(additional_hads_anxiety ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_anxiety_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_anxiety_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_anxiety_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_anxiety_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "E. Mean Left Amygdala KFA and Anxiety, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** anxiety p < 0.001",
      "anxiety p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_KFA_by_depression_age_sex_income_education <- 
  lm(additional_hads_depression ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_KFA_by_depression_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_KFA_by_depression_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_KFA_by_depression_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_KFA_by_depression_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Amygdala KFA and Depression, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** depression p < 0.001",
      "depression p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_story_age_sex_income_education <- 
  lm(homeint_storyrecall_d ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_story_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_story_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_story_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_story_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Right Amygdala KFA and Memory, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** story p < 0.001",
      "story p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_anxiety_age_sex_income_education <- 
  lm(additional_hads_anxiety ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_anxiety_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_anxiety_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_anxiety_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_anxiety_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "F. Mean Right Amygdala KFA and Anxiety, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** anxiety p < 0.001",
      "anxiety p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_KFA_by_depression_age_sex_income_education <- 
  lm(additional_hads_depression ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_KFA_by_depression_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_KFA_by_depression_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_KFA_by_depression_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_KFA_by_depression_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Right Amygdala KFA and Depression, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** depression p < 0.001",
      "depression p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

### 3.4	Correlations of Significant Imaging Metrics with Measures of Cardiovascular Health ###
left_amygdala_kfa_by_bp_dia_age_sex_income_education <- 
  lm(bp_dia_mean_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_kfa_by_bp_dia_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_kfa_by_bp_dia_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_kfa_by_bp_dia_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_kfa_by_bp_dia_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Amygdala KFA and Diastolic BP, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** diastolic p < 0.001",
      "diastolic p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_kfa_by_bp_sys_age_sex_income_education <- 
  lm(bp_sys_mean_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_kfa_by_bp_sys_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_kfa_by_bp_sys_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_kfa_by_bp_sys_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_kfa_by_bp_sys_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Left Amygdala KFA and systolic BP, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** systolic p < 0.001",
      "systolic p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_kfa_by_pulse_age_sex_income_education <- 
  lm(pulse_mean_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_kfa_by_pulse_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_kfa_by_pulse_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_kfa_by_pulse_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_kfa_by_pulse_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Left Amygdala KFA and Pulse, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** pulse p < 0.001",
      "pulse p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_kfa_by_bmi_age_sex_income_education <- 
  lm(bmi_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_amygdala)
corr_p_value <- summary(left_amygdala_kfa_by_bmi_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_amygdala_kfa_by_bmi_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_amygdala_kfa_by_bmi_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_amygdala_kfa_by_bmi_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    y = BMI~(kg/m^2),
    title = "A. Mean Left Amygdala KFA and BMI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** BMI p < 0.001",
      "\\* BMI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) 
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_kfa_by_bp_dia_age_sex_income_education <- 
  lm(bp_dia_mean_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_kfa_by_bp_dia_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_kfa_by_bp_dia_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_kfa_by_bp_dia_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_kfa_by_bp_dia_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Right Amygdala KFA and Diastolic BP, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** diastolic p < 0.001",
      "diastolic p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_kfa_by_bp_sys_age_sex_income_education <- 
  lm(bp_sys_mean_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_kfa_by_bp_sys_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_kfa_by_bp_sys_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_kfa_by_bp_sys_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_kfa_by_bp_sys_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Right Amygdala KFA and systolic BP, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** systolic p < 0.001",
      "systolic p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_kfa_by_pulse_age_sex_income_education <- 
  lm(pulse_mean_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_kfa_by_pulse_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_kfa_by_pulse_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_kfa_by_pulse_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_kfa_by_pulse_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Right Amygdala KFA and Pulse, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** pulse p < 0.001",
      "pulse p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

right_amygdala_kfa_by_bmi_age_sex_income_education <- 
  lm(bmi_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     right_amygdala)
corr_p_value <- summary(right_amygdala_kfa_by_bmi_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(right_amygdala_kfa_by_bmi_age_sex_income_education)$adj.r.squared
interaction_p <- summary(right_amygdala_kfa_by_bmi_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(right_amygdala_kfa_by_bmi_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    y = BMI~(kg/m^2),
    title = "F. Mean Right Amygdala KFA and BMI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** BMI p < 0.001",
      "BMI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
