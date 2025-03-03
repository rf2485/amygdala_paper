source("1_data_preparation.R")
library(gtsummary)
library(ggeffects)
library(ggtext)

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

volumes <- left_join(scd_status, aseg) %>% #join volumes with SCD status
  mutate(across(c(Left.Lateral.Ventricle:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) %>% #normalize by intracranial volume
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  rename(Left.Cerebral.White.Matter	= lhCerebralWhiteMatterVol, 
         Right.Cerebral.White.Matter	= rhCerebralWhiteMatterVol) %>%
  select(!c((BrainSegVol:CortexVol), (CerebralWhiteMatterVol:EstimatedTotalIntraCranialVol)))
subcort_gm_volumes_table <- volumes %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
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
  rename(participant_id = Measure.mean) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) 
subcort_gm_NDI_table <- fit_NDI %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              label = subcort_labels) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

tbl_merge(list(subcort_gm_volumes_table, subcort_gm_NDI_table), 
          tab_spanner = c("**Volume**", "**NDI**"))

fit_FWF <- aseg2fit_FWF %>%
  rename(participant_id = Measure.mean) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) 
subcort_gm_FWF_table <- fit_FWF %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
               label = subcort_labels) %>%
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
              label = subcort_labels) %>%
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
  rename(participant_id = Measure.mean) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) 
subcort_gm_ODI_table <- fit_ODI %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
               label = subcort_labels) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
subcort_gm_ODI_table_all <- fit_ODI %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              label = subcort_labels) %>%
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
  rename(participant_id = Measure.mean) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) 
subcort_gm_fa_table_all <- dti_fa %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
               label = subcort_labels) %>%
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
  rename(participant_id = Measure.mean) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .) 
subcort_gm_md_table_all <- dti_md %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              label = subcort_labels) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% bold_p() %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

tbl_merge(list(subcort_gm_fa_table_all, subcort_gm_md_table_all), 
          tab_spanner = c("**FA**", "**MD**"))

dki_kfa <- aseg2dki_kfa %>%
  rename(participant_id = Measure.mean) %>%
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  left_join(scd_status, .)
subcort_gm_kfa_table <- dki_kfa %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              label = subcort_labels) %>%
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
              label = subcort_labels) %>%
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
  rename(participant_id = Measure.mean) %>%
  #remove biologically implausible values (Veraat, Hecke, and Sijbers 2011)
  mutate(across(where(is.double), ~ replace(., . < 0, NA))) %>%
  mutate(across(where(is.double), ~ replace(., . > 4, NA))) %>%
  left_join(scd_status, .) 
subcort_gm_mk_table <- dki_mk %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              label = subcort_labels) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")
subcort_gm_mk_table_all <- dki_mk %>% 
  select(any_of(subcort_gm)) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              label = subcort_labels) %>%
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
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
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
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
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

### 3.2	Correlations Between Imaging Metrics ####

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
  lm(volume ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_volume_all) #age significant
left_putamen_FWF_by_NDI_all <- 
  lm(fit_NDI ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_NDI_all) #age and sex significant
left_putamen_FWF_by_ODI_all <- 
  lm(fit_ODI ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_ODI_all) #age significant
left_putamen_FWF_by_FA_all <- 
  lm(dti_fa ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_FA_all) #age and sex significant
left_putamen_FWF_by_MD_all <- 
  lm(dti_md ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_MD_all) #age significant
left_putamen_FWF_by_KFA_all <- 
  lm(dki_kfa ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_KFA_all) #sex and ethnicity significant
left_putamen_FWF_by_MK_all <- 
  lm(dki_mk ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_MK_all) #age significant
left_putamen_FWF_by_mtr_tr30_all <- 
  lm(mtr_tr30 ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_mtr_tr30_all) #none significant
left_putamen_FWF_by_mtr_tr50_all <- 
  lm(mtr_tr50 ~ fit_FWF * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_putamen)
anova(left_putamen_FWF_by_mtr_tr50_all) #age and sex significant

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
  lm(volume ~ dki_kfa * SCD + age + Sex + age_education_completed + Income + Ethnicity,
     left_amygdala)
anova(left_amygdala_KFA_by_volume_all) #age, education, Income significant

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

##models 
left_putamen_FWF_by_volume_age_sex_income_education <- 
  lm(volume ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_volume_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_volume_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_volume_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_volume_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "Mean Left Putamen FWF and Volume, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** volume p < 0.001",
      "volume p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_NDI_age_sex_income_education <- 
  lm(fit_NDI ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_NDI_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_NDI_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_NDI_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_NDI_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Putamen FWF and NDI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** model p < 0.001",
      "\\* NDI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> * group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_ODI_age_sex_income_education <- 
  lm(fit_ODI ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_ODI_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_ODI_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_ODI_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_ODI_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Left Putamen FWF and ODI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** ODI p < 0.001",
      "ODI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> * group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_FA_age_sex_income_education <- 
  lm(dti_fa ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_FA_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_FA_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_FA_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_FA_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Left Putamen FWF and FA, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** FA p < 0.001",
      # "* FA p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_MD_age_sex_income_education <- 
  lm(dti_md ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_MD_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_MD_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_MD_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_MD_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Putamen FWF and MD, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** MD p < 0.001",
      # "* MD p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> \\* group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_KFA_age_sex_income_education <- 
  lm(dki_kfa ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_KFA_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_KFA_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_KFA_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_KFA_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "E. Mean Left Putamen FWF and KFA, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** KFA p < 0.001",
      # "* KFA p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_MK_age_sex_income_education <- 
  lm(dki_mk ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_MK_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_MK_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_MK_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_MK_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "F. Mean Left Putamen FWF and MK, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** MK p < 0.001",
      # "* MK p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> \\* group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_mtr_tr30_age_sex_income_education <- 
  lm(mtr_tr30 ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_mtr_tr30_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_mtr_tr30_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_mtr_tr30_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_mtr_tr30_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "G. Mean Left Putamen FWF and MTR TR=30ms, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** MTR TR=30ms p < 0.001",
      "\\* MTR TR=30ms p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_mtr_tr50_age_sex_income_education <- 
  lm(mtr_tr50 ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_mtr_tr50_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_mtr_tr50_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_mtr_tr50_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_mtr_tr50_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "H. Mean Left Putamen FWF and MTR TR=50ms, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** MTR TR=50ms p < 0.001",
      "MTR TR=50ms p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

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

left_accumbens_area_KFA_by_volume_age_sex_income_education <- 
  lm(volume ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_volume_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_volume_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_volume_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_volume_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Accumbens Area KFA and Volume, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** volume p < 0.001",
      "volume p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> * group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_KFA_by_FWF_age_sex_income_education <- 
  lm(fit_FWF ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_FWF_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_FWF_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_FWF_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_FWF_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Accumbens Area KFA and FWF, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** FWF p < 0.001",
      # "FWF p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_KFA_by_NDI_age_sex_income_education <- 
  lm(fit_NDI ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_NDI_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_NDI_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_NDI_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_NDI_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Left Accumbens Area KFA and NDI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** model p < 0.001",
      "** NDI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_KFA_by_ODI_age_sex_income_education <- 
  lm(fit_ODI ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_ODI_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_ODI_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_ODI_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_ODI_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Left Accumbens Area KFA and ODI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** ODI p < 0.001",
      "ODI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_KFA_by_FA_age_sex_income_education <- 
  lm(dti_fa ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_FA_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_FA_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_FA_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_FA_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Left Accumbens Area KFA and FA, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** FA p < 0.001",
      # "FA p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_KFA_by_MD_age_sex_income_education <- 
  lm(dti_md ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_MD_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_MD_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_MD_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_MD_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Accumbens Area KFA and MD, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** MD p < 0.001",
      # "MD p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_KFA_by_MK_age_sex_income_education <- 
  lm(dki_mk ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_MK_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_MK_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_MK_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_MK_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "E. Mean Left Accumbens Area KFA and MK, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** MK p < 0.001",
      # "\\* MK p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_KFA_by_mtr_tr30_age_sex_income_education <- 
  lm(mtr_tr30 ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_mtr_tr30_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_mtr_tr30_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_mtr_tr30_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_mtr_tr30_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Left Accumbens Area KFA and MTR TR=30ms, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** MTR TR=30ms p < 0.001",
      "MTR TR=30ms p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_KFA_by_mtr_tr50_age_sex_income_education <- 
  lm(mtr_tr50 ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_mtr_tr50_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_mtr_tr50_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_mtr_tr50_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_mtr_tr50_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = ". Mean Left Accumbens Area KFA and MTR TR=50ms, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** MTR TR=50ms p < 0.001",
      "MTR TR=50ms p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_MK_by_volume_age_sex_income_education <- 
  lm(volume ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_MK_by_volume_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_MK_by_volume_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_MK_by_volume_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_MK_by_volume_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Accumbens Area MK and Volume, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** volume p < 0.001",
      "volume p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_MK_by_FWF_age_sex_income_education <- 
  lm(fit_FWF ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_MK_by_FWF_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_MK_by_FWF_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_MK_by_FWF_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_MK_by_FWF_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Accumbens Area MK and FWF, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** FWF p < 0.001",
      # "FWF p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> \\* group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_MK_by_NDI_age_sex_income_education <- 
  lm(fit_NDI ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_MK_by_NDI_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_MK_by_NDI_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_MK_by_NDI_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_MK_by_NDI_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Left Accumbens Area MK and NDI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** model p < 0.001",
      # "** NDI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_MK_by_ODI_age_sex_income_education <- 
  lm(fit_ODI ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_MK_by_ODI_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_MK_by_ODI_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_MK_by_ODI_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_MK_by_ODI_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Left Accumbens Area MK and ODI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** ODI p < 0.001",
      " \\* ODI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_MK_by_FA_age_sex_income_education <- 
  lm(dti_fa ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_MK_by_FA_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_MK_by_FA_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_MK_by_FA_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_MK_by_FA_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Left Accumbens Area MK and FA, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** FA p < 0.001",
      "FA p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_MK_by_MD_age_sex_income_education <- 
  lm(dti_md ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_MK_by_MD_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_MK_by_MD_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_MK_by_MD_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_MK_by_MD_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Accumbens Area MK and MD, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** MD p < 0.001",
      # "MD p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_MK_by_KFA_age_sex_income_education <- 
  lm(dki_kfa ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_MK_by_KFA_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_MK_by_KFA_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_MK_by_KFA_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_MK_by_KFA_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "E. Mean Left Accumbens Area MK and KFA, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      "*** MK p < 0.001",
      # "\\* MK p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_MK_by_mtr_tr30_age_sex_income_education <- 
  lm(mtr_tr30 ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_MK_by_mtr_tr30_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_MK_by_mtr_tr30_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_MK_by_mtr_tr30_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_MK_by_mtr_tr30_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Left Accumbens Area MK and MTR TR=30ms, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** MTR TR=30ms p < 0.001",
      "MTR TR=30ms p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_MK_by_mtr_tr50_age_sex_income_education <- 
  lm(mtr_tr50 ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_MK_by_mtr_tr50_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_MK_by_mtr_tr50_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_MK_by_mtr_tr50_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_MK_by_mtr_tr50_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Accumbens Area MK and MTR TR=50ms, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** MTR TR=50ms p < 0.001",
      "MTR TR=50ms p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())


### 3.3	Correlations of Significant Imaging Metrics with Memory, Anxiety, and Depression ###
left_putamen_FWF_by_story_age_sex_income_education <- 
  lm(homeint_storyrecall_d ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_story_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_story_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_story_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_story_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Putamen FWF and Memory, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** story p < 0.001",
      "story p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> \\* group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_anxiety_age_sex_income_education <- 
  lm(additional_hads_anxiety ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_anxiety_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_anxiety_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_anxiety_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_anxiety_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Putamen FWF and Anxiety, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** anxiety p < 0.001",
      "anxiety p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_depression_age_sex_income_education <- 
  lm(additional_hads_depression ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_depression_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_depression_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_depression_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_depression_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "B. Mean Left Putamen FWF and Depression, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** depression p < 0.001",
      "\\* depression p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

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

left_accumbens_area_KFA_by_story_age_sex_income_education <- 
  lm(homeint_storyrecall_d ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_story_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_story_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_story_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_story_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Left Accumbens Area KFA and Memory, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** story p < 0.001",
      "story p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_KFA_by_anxiety_age_sex_income_education <- 
  lm(additional_hads_anxiety ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_anxiety_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_anxiety_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_anxiety_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_anxiety_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Accumbens Area KFA and Anxiety, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** anxiety p < 0.001",
      "anxiety p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_KFA_by_depression_age_sex_income_education <- 
  lm(additional_hads_depression ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_KFA_by_depression_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_KFA_by_depression_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_KFA_by_depression_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_KFA_by_depression_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "G. Mean Left Accumbens Area KFA and Depression, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** depression p < 0.001",
      "\\* depression p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> \\* group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_mk_by_story_age_sex_income_education <- 
  lm(homeint_storyrecall_d ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_mk_by_story_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_mk_by_story_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_mk_by_story_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_mk_by_story_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Accumbens Area MK and Memory, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** story p < 0.001",
      "story p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_mk_by_anxiety_age_sex_income_education <- 
  lm(additional_hads_anxiety ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_mk_by_anxiety_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_mk_by_anxiety_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_mk_by_anxiety_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_mk_by_anxiety_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "H. Mean Left Accumbens Area MK and Anxiety, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** anxiety p < 0.001",
      "anxiety p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_mk_by_depression_age_sex_income_education <- 
  lm(additional_hads_depression ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_mk_by_depression_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_mk_by_depression_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_mk_by_depression_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_mk_by_depression_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Left Accumbens Area MK and Depression, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** depression p < 0.001",
      "depression p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

### 3.4	Correlations of Significant Imaging Metrics with Measures of Cardiovascular Health ###
left_putamen_FWF_by_bp_dia_age_sex_income_education <- 
  lm(bp_dia_mean_cardio ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_bp_dia_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_bp_dia_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_bp_dia_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_bp_dia_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Putamen FWF and Diastolic BP, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** diastolic p < 0.001",
      "\\* diastolic p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_bp_sys_age_sex_income_education <- 
  lm(bp_sys_mean_cardio ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_bp_sys_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_bp_sys_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_bp_sys_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_bp_sys_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Putamen FWF and systolic BP, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** systolic p < 0.001",
      "systolic p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_pulse_age_sex_income_education <- 
  lm(pulse_mean_cardio ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_pulse_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_pulse_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_pulse_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_pulse_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "A. Mean Left Putamen FWF and Pulse, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** pulse p < 0.001",
      "pulse p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_putamen_FWF_by_bmi_age_sex_income_education <- 
  lm(bmi_cardio ~ fit_FWF * SCD + age + Sex + age_education_completed + Income,
     left_putamen)
corr_p_value <- summary(left_putamen_FWF_by_bmi_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_putamen_FWF_by_bmi_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_putamen_FWF_by_bmi_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_putamen_FWF_by_bmi_age_sex_income_education, c( "fit_FWF[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    y = BMI~(kg/m^2),
    title = "E. Mean Left Putamen FWF and BMI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** BMI p < 0.001",
      "BMI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

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

left_accumbens_area_kfa_by_bp_dia_age_sex_income_education <- 
  lm(bp_dia_mean_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_kfa_by_bp_dia_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_kfa_by_bp_dia_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_kfa_by_bp_dia_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_kfa_by_bp_dia_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "C. Mean Left Accumbens Area KFA and Diastolic BP, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** diastolic p < 0.001",
      "diastolic p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_kfa_by_bp_sys_age_sex_income_education <- 
  lm(bp_sys_mean_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_kfa_by_bp_sys_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_kfa_by_bp_sys_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_kfa_by_bp_sys_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_kfa_by_bp_sys_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Accumbens Area KFA and systolic BP, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** systolic p < 0.001",
      "systolic p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_kfa_by_pulse_age_sex_income_education <- 
  lm(pulse_mean_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_kfa_by_pulse_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_kfa_by_pulse_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_kfa_by_pulse_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_kfa_by_pulse_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Accumbens Area KFA and Pulse, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** pulse p < 0.001",
      "pulse p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_kfa_by_bmi_age_sex_income_education <- 
  lm(bmi_cardio ~ dki_kfa * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_kfa_by_bmi_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_kfa_by_bmi_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_kfa_by_bmi_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_kfa_by_bmi_age_sex_income_education, c( "dki_kfa[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    y = BMI~(kg/m^2),
    title = "G. Mean Left Accumbens Area KFA and BMI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** BMI p < 0.001",
      "BMI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_mk_by_bp_dia_age_sex_income_education <- 
  lm(bp_dia_mean_cardio ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_mk_by_bp_dia_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_mk_by_bp_dia_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_mk_by_bp_dia_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_mk_by_bp_dia_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "D. Mean Left Accumbens Area MK and Diastolic BP, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** diastolic p < 0.001",
      "diastolic p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_mk_by_bp_sys_age_sex_income_education <- 
  lm(bp_sys_mean_cardio ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_mk_by_bp_sys_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_mk_by_bp_sys_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_mk_by_bp_sys_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_mk_by_bp_sys_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "E. Mean Left Accumbens Area MK and systolic BP, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** systolic p < 0.001",
      "systolic p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_mk_by_pulse_age_sex_income_education <- 
  lm(pulse_mean_cardio ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_mk_by_pulse_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_mk_by_pulse_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_mk_by_pulse_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_mk_by_pulse_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    title = "E. Mean Left Accumbens Area MK and Pulse, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** pulse p < 0.001",
      "pulse p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_accumbens_area_mk_by_bmi_age_sex_income_education <- 
  lm(bmi_cardio ~ dki_mk * SCD + age + Sex + age_education_completed + Income,
     left_accumbens_area)
corr_p_value <- summary(left_accumbens_area_mk_by_bmi_age_sex_income_education)$coefficients[2,4]
adj_r_squared <- summary(left_accumbens_area_mk_by_bmi_age_sex_income_education)$adj.r.squared
interaction_p <- summary(left_accumbens_area_mk_by_bmi_age_sex_income_education)$coefficients[12,4]
pr <- predict_response(left_accumbens_area_mk_by_bmi_age_sex_income_education, c( "dki_mk[all]", "SCD"))
plot(pr, show_data = T, dot_alpha = 1) + 
  labs(
    y = BMI~(kg/m^2),
    title = "H. Mean Left Accumbens Area MK and BMI, \n Corrected by Age, Sex, Education, and Income",
    subtitle = paste0(
      # "*** BMI p < 0.001",
      "BMI p = ", signif(corr_p_value, 2),
      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2),
      "<br> group interaction p = ", signif(interaction_p, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
