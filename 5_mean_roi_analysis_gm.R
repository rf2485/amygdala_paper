source("1_data_preparation_gm.R")
library(tidyverse)
library(labelVector)
library(gtsummary)
library(ggeffects)
library(ggtext)
library(patchwork)
library(corrplot)


failed_qc <- c('sub-CC510255', #SCD abnormality in left temporal pole
               'sub-CC510438', #CTL abnormality in left frontal lobe
               'sub-CC620821', #SCD segmentation errors from large ventricles
               'sub-CC621011', #CTL segmentation errors from large ventricles
               'sub-CC621080', #SCD segmentation errors
               'sub-CC710551', #CTL motion artifacts in DWI
               'sub-CC711027', #SCD severe motion artifacts in T1
               'sub-CC721434'  #CTL segmentation errors from large ventricles
)

my_ES_test <- function(data, variable, by, ...) {
  rstatix::cohens_d(data, as.formula(glue::glue("{variable} ~ {by}")))$effsize
}
my_cramer_v <- function(data, variable, by, ...) {
  table(data[[variable]], data[[by]]) %>%
    rstatix::cramer_v()
}
theme_gtsummary_mean_sd()

### Memory Self-report within SCD ###
scd_10mq <- dwi_over_55 %>% filter(SCD == "SCD") %>%
  filter(!participant_id %in% failed_qc) %>%
  select(participant_id, homeint_v231:homeint_v240) %>%
  naniar::replace_with_na_if(.predicate = is.numeric, condition = ~.x > 4) %>%
  mutate(ten_mq_total = rowSums(across(homeint_v231:homeint_v240), na.rm = T))

hist(scd_10mq$ten_mq_total)
median(scd_10mq$ten_mq_total)

### 3.1	Group differences between SCD and controls ###

#extract SCD status and demographics for participants
scd_status <- dwi_over_55 %>% select(participant_id, SCD, mt_tr, Income, Ethnicity, 
                                     Sex, age, age_education_completed,homeint_storyrecall_d,
                                     additional_hads_anxiety, additional_hads_depression,
                                     objprpos_emotional_mem, objprneu_emotional_mem, 
                                     objprneg_emotional_mem) %>%
  filter(!participant_id %in% failed_qc) %>%
  rename("Group" = "SCD") %>%
  mutate(objmem_emotion_enhance = 
           ((objprpos_emotional_mem + objprneg_emotional_mem)/2) - objprneu_emotional_mem,
         objmem_pos_bias = objprpos_emotional_mem - objprneg_emotional_mem) %>%
  select(!objprpos_emotional_mem, !objprneu_emotional_mem, !objprneg_emotional_mem)
scd_status$additional_hads_depression[scd_status$additional_hads_depression > 21] <- NA
scd_status <- scd_10mq %>% select(participant_id, ten_mq_total) %>%
  left_join(scd_status, .)

scd_status_10mq_thresh <- scd_status %>%
  filter(Group == "Control" | ten_mq_total > 10)

#demographics table
scd_status %>% select(Group, age, age_education_completed, Sex, Income, Ethnicity) %>%
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({min}-{max})") %>%
  # add_difference(test = list(all_continuous() ~ 'cohens_d')) %>%
  # modify_column_hide(conf.low) %>%
  add_stat(
    fns = list(all_continuous() ~ my_ES_test,
               all_categorical() ~ my_cramer_v)) %>%  add_p() %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**",
                add_stat_1 ~ "**Effect size**",
                label ~ "") %>%
  bold_labels() %>%
  modify_footnote(add_stat_1 ~ "Cohen's D; Cramer's V") %>%
  as_gt() %>%
  gt::gtsave(filename = "demographics_table.docx")

scd_status_10mq_thresh %>% select(Group, age, age_education_completed, Sex, Income, Ethnicity) %>%
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({min}-{max})") %>%
  # add_difference(test = list(all_continuous() ~ 'cohens_d')) %>%
  # modify_column_hide(conf.low) %>%
  add_stat(
    fns = list(all_continuous() ~ my_ES_test,
               all_categorical() ~ my_cramer_v)) %>%  add_p() %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**",
                add_stat_1 ~ "**Effect size**",
                label ~ "") %>%
  bold_labels() %>%
  modify_footnote(add_stat_1 ~ "Cohen's D; Cramer's V")
  
#cog and psych table
scd_status %>% select(Group, homeint_storyrecall_d, 
                      additional_hads_anxiety, additional_hads_depression,
                      objmem_emotion_enhance, objmem_pos_bias) %>%
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})") %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "",
                estimate ~ "**Effect Size**") %>%
  bold_labels() #%>%
  # as_gt() %>%
  # gt::gtsave(filename = "cog_and_psych_table.docx")

scd_status_10mq_thresh %>% select(Group, homeint_storyrecall_d, additional_hads_anxiety, additional_hads_depression) %>%
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})") %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "",
                estimate ~ "**Effect Size**") %>%
  bold_labels()

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

#### focus on amygdala ####
left_amygdala <- volumes %>% dplyr::select(participant_id:ten_mq_total, Left.Amygdala) %>%
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
left_amygdala <- dti_ad %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(dti_ad = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dti_rd %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(dti_rd = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_kfa %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(dki_kfa = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_mk %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(dki_mk = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_ak %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(dki_ak = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_rk %>% dplyr::select(participant_id, Left.Amygdala) %>%
  rename(dki_rk = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala$region <- "Left Amygdala"

left_amygdala <- set_label(left_amygdala,
                           fit_NDI = "NDI",
                           fit_FWF = "FWF",
                           fit_ODI = "ODI",
                           dti_fa = "FA",
                           dti_md = "MD",
                           dti_ad = "AxD",
                           dti_rd = "RD",
                           dki_kfa = "KFA",
                           dki_mk = "MK",
                           dki_ak = "AK",
                           dki_rk = "RK"
)

left_amygdala_10mq_thresh <- left_amygdala %>% 
  filter(Group == "Control" | ten_mq_total > 10)

right_amygdala <- volumes %>% dplyr::select(participant_id:ten_mq_total, Right.Amygdala) %>%
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
right_amygdala <- dti_ad %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(dti_ad = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dti_rd %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(dti_rd = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_kfa %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(dki_kfa = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_mk %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(dki_mk = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_ak %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(dki_ak = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_rk %>% dplyr::select(participant_id, Right.Amygdala) %>%
  rename(dki_rk = "Right.Amygdala") %>% full_join(right_amygdala, .)

right_amygdala$region <- "Right Amygdala"
right_amygdala <- set_label(right_amygdala,
                            fit_NDI = "NDI",
                            fit_FWF = "FWF",
                            fit_ODI = "ODI",
                            dti_fa = "FA",
                            dti_md = "MD",
                            dti_ad = "AxD",
                            dti_rd = "RD",
                            dki_kfa = "KFA",
                            dki_mk = "MK",
                            dki_ak = "AK",
                            dki_rk = "RK"
)

right_amygdala_10mq_thresh <- right_amygdala %>% 
  filter(Group == "Control" | ten_mq_total > 10)

#### between groups DTI and DKI differences ####
left_amygdala %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(dti_fa, dti_md, dki_kfa, dki_mk)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p(t=0.025) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**") #%>%
  # as_gt() %>%
  # gt::gtsave(filename = "left_amygdala_dti_dki_table.docx")

left_amygdala_10mq_thresh %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(dti_fa, dti_md, dki_kfa, dki_mk)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p(t=0.025) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**")

right_amygdala %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(dti_fa, dti_md, dki_kfa, dki_mk)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p(t=0.025) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**") #%>%
  # as_gt() %>%
  # gt::gtsave(filename = "right_amygdala_dti_dki_table.docx")

right_amygdala_10mq_thresh %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(dti_fa, dti_md, dki_kfa, dki_mk)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p(t=0.025) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**")

#### Within SCD correlations between DTI, DKI, volume, and NODDI ####
trace(corrplot, edit = T) #change line 446 to text(pos.pNew[, 1][sig.locs], (pos.pNew[, 2][sig.locs])+0.25, 

left_dti_dki_matrix <- left_amygdala %>% 
  filter(Group == "SCD") %>%
  select(dti_fa, dti_md, dki_mk, dki_kfa) %>%
  rename_with( ~ paste0(.x, "_left"))
right_dti_dki_matrix <- right_amygdala %>% 
  filter(Group == "SCD") %>%
  select(dti_fa, dti_md, dki_mk, dki_kfa) %>%
  rename_with( ~ paste0(.x, "_right"))
dti_dki_matrix <- cbind(left_dti_dki_matrix, right_dti_dki_matrix)

left_vol_noddi_matrix <- left_amygdala %>%
  filter(Group == "SCD") %>%
  select(volume, fit_NDI, fit_FWF, fit_ODI) %>%
  rename_with( ~ paste0(.x, "_left"))
right_vol_noddi_matrix <- right_amygdala %>%
  filter(Group == "SCD") %>%
  select(volume, fit_NDI, fit_FWF, fit_ODI) %>%
  rename_with( ~ paste0(.x, "_right"))
vol_noddi_matrix <- cbind(left_vol_noddi_matrix, right_vol_noddi_matrix)

corr_matrix_pearson <- psych::corr.test(dti_dki_matrix, vol_noddi_matrix)
right_corr_matrix_pearson <- psych::corr.test(right_dti_dki_matrix, right_vol_noddi_matrix)
right_corr_matrix_spearman <- psych::corr.test(right_dti_dki_matrix, right_vol_noddi_matrix,
                                               method = "spearman")
left_corr_matrix_pearson <- psych::corr.test(left_dti_dki_matrix, left_vol_noddi_matrix)
left_corr_matrix_spearman <- psych::corr.test(left_dti_dki_matrix, left_vol_noddi_matrix,
                                               method = "spearman")

# colnames(corr_matrix_pearson$r) <- c("Left Volume", "Left NDI", "Left FWF", "Left ODI",
                                   # "Right Volume", "Right NDI", "Right FWF", "Right ODI")

corrplot(left_corr_matrix_pearson$r, p.mat = left_corr_matrix_pearson$p, method = 'color',
         addCoef.col = "black",
                   sig.level = c(0.001, 0.01, 0.025), insig = 'label_sig', pch.cex = 0.9)

corrplot(left_corr_matrix_spearman$r, p.mat = left_corr_matrix_spearman$p, method = 'color',
         addCoef.col = "black",
         sig.level = c(0.001, 0.01, 0.025), insig = 'label_sig', pch.cex = 0.9)

# left_dti_dki_matrix <- left_amygdala_10mq_thresh %>% 
#   filter(Group == "SCD") %>%
#   select(dti_fa, dki_kfa) %>%
#   rename_with( ~ paste0(.x, "_left"))
# right_dti_dki_matrix <- right_amygdala_10mq_thresh %>% 
#   filter(Group == "SCD") %>%
#   select(dki_kfa) %>%
#   rename_with( ~ paste0(.x, "_right"))
# dti_dki_matrix <- cbind(left_dti_dki_matrix, right_dti_dki_matrix)
# 
# left_vol_noddi_matrix <- left_amygdala_10mq_thresh %>%
#   filter(Group == "SCD") %>%
#   select(volume, fit_NDI, fit_FWF, fit_ODI) %>%
#   rename_with( ~ paste0(.x, "_left"))
# right_vol_noddi_matrix <- right_amygdala_10mq_thresh %>%
#   filter(Group == "SCD") %>%
#   select(volume, fit_NDI, fit_FWF, fit_ODI) %>%
#   rename_with( ~ paste0(.x, "_right"))
# vol_noddi_matrix <- cbind(left_vol_noddi_matrix, right_vol_noddi_matrix)
# 
# corr_matrix_pearson <- psych::corr.test(dti_dki_matrix, vol_noddi_matrix)
# colnames(corr_matrix_pearson$r) <- c("Left Volume", "Left NDI", "Left FWF", "Left ODI",
#                                      "Right Volume", "Right NDI", "Right FWF", "Right ODI")
# 
# corrplot(corr_matrix_pearson$r, p.mat = corr_matrix_pearson$p, method = 'color',
#          addCoef.col = "black",
#          sig.level = c(0.001, 0.01, 0.025), insig = 'label_sig', pch.cex = 0.9)

### Between Groups Volume and NODDI Differences ###
left_amygdala %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(volume, fit_NDI, fit_FWF, fit_ODI)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p(t=0.025) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**") #%>%
  # as_gt() %>%
  # gt::gtsave(filename = "left_amygdala_vol_noddi_table.docx")

left_amygdala_10mq_thresh %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(volume, fit_NDI, fit_FWF, fit_ODI)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p(t=0.025) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**")

right_amygdala %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(volume, fit_NDI, fit_FWF, fit_ODI)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p(t=0.025) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**") #%>%
  # as_gt() %>%
  # gt::gtsave(filename = "right_amygdala_vol_noddi_table.docx")

right_amygdala_10mq_thresh %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(volume, fit_NDI, fit_FWF, fit_ODI)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p(t=0.025) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**")

##### Linear Models ######
glm_left_amygdala_dti_fa_age_int <- lm(dti_fa ~ age * Group, left_amygdala)
summary(glm_left_amygdala_dti_fa_age_int)
glm_left_amygdala_dti_fa_age_int <- lm(dti_fa ~ age * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_dti_fa_age_int)

tbl_regression(glm_left_amygdala_dti_fa_age_int, 
               estimate_fun = label_style_sigfig(digits = 3),
               show_single_row = c(Group, `age:Group`)) %>% 
  bold_p(t=0.025) %>%
  add_glance_table(include = c(r.squared, p.value), 
                   label = list(p.value = "overall p-value")) %>%
  modify_spanning_header(c(estimate, conf.low, conf.high, p.value) ~ "**Left Amygdala FA**") #%>%
  # as_gt() %>%
  # gt::gtsave(filename = "left_amygdala_dti_fa_age_int.docx")

glm_left_amygdala_fit_FWF_age_int <- lm(fit_FWF ~ age * Group, left_amygdala)
summary(glm_left_amygdala_fit_FWF_age_int)
glm_left_amygdala_fit_FWF_age_int <- lm(fit_FWF ~ age * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_fit_FWF_age_int)

tbl_regression(glm_left_amygdala_fit_FWF_age_int, 
               estimate_fun = label_style_sigfig(digits = 3),
               show_single_row = c(Group, `age:Group`)) %>% 
  bold_p(t=0.025) %>%
  add_glance_table(include = c(r.squared, p.value), 
                   label = list(p.value = "overall p-value")) %>%
  modify_spanning_header(c(estimate, conf.low, conf.high, p.value) ~ "**Left Amygdala FWF**") #%>%
  # as_gt() %>%
  # gt::gtsave(filename = "left_amygdala_fit_FWF_age_int.docx")

glm_left_amygdala_dki_kfa_age_int <- lm(dki_kfa ~ age * Group, left_amygdala)
summary(glm_left_amygdala_dki_kfa_age_int)
glm_left_amygdala_dki_kfa_age_int <- lm(dki_kfa ~ age * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_dki_kfa_age_int)

tbl_regression(glm_left_amygdala_dki_kfa_age_int, 
                 estimate_fun = label_style_sigfig(digits = 3),
                 show_single_row = c(Group, `age:Group`)) %>% 
  bold_p(t=0.025) %>%
  add_glance_table(include = c(r.squared, p.value), 
                   label = list(p.value = "overall p-value")) %>%
  modify_spanning_header(c(estimate, conf.low, conf.high, p.value) ~ "**Left Amygdala KFA**") #%>%
  # as_gt() %>%
  # gt::gtsave(filename = "left_amygdala_dki_kfa_age_int.docx") 

glm_right_amygdala_dki_kfa_age_int <- lm(dki_kfa ~ age * Group, right_amygdala)
summary(glm_right_amygdala_dki_kfa_age_int)
glm_right_amygdala_dki_kfa_age_int <- lm(dki_kfa ~ age * Group, right_amygdala_10mq_thresh)
summary(glm_right_amygdala_dki_kfa_age_int)

tbl_regression(glm_right_amygdala_dki_kfa_age_int, 
               estimate_fun = label_style_sigfig(digits = 3),
               show_single_row = c(Group, `age:Group`)) %>% 
  bold_p(t=0.025) %>%
  add_glance_table(include = c(r.squared, p.value), 
                         label = list(p.value = "overall p-value")) %>%
  modify_spanning_header(c(estimate, conf.low, conf.high, p.value) ~ "**Right Amygdala KFA**") #%>%
  # as_gt() %>%
  # gt::gtsave(filename = "right_amygdala_dki_kfa_age_int.docx")
 
pr <- predict_response(glm_left_amygdala_dti_fa_age_int, c( "age[all]"))
plot_left_amygdala_dti_fa_age <- plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.1) + 
   theme_bw(base_size = 10) +
  labs(y = "Left Amygdala FA",
       title = paste0(
         symnum(summary(glm_left_amygdala_dti_fa_age_int)$coefficients[2,4], 
                corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
         " slope p = ",
         format.pval(summary(glm_left_amygdala_dti_fa_age_int)$coefficients[2,4], digits = 2, eps = 0.001)
         )) +
   theme(
     plot.title = element_markdown(),
     legend.title=element_blank())

pr <- predict_response(glm_left_amygdala_fit_FWF_age_int, c( "age[all]"))
plot_left_amygdala_fit_FWF_age <- plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.1) + 
  theme_bw(base_size = 10) +
  labs(y = "Left Amygdala FWF",
       title = paste0(
         symnum(summary(glm_left_amygdala_fit_FWF_age_int)$coefficients[2,4], 
                corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
         " slope p = ",
         format.pval(summary(glm_left_amygdala_fit_FWF_age_int)$coefficients[2,4], digits = 1, eps = 0.001)
       )) +
  theme(
    plot.title = element_markdown(),
    legend.title=element_blank())

diffusion_age_plots <- plot_left_amygdala_dti_fa_age + 
  plot_left_amygdala_fit_FWF_age + 
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave(filename = "diffusion_age_plots.tif", diffusion_age_plots,
       width = 6.5, height = 3, dpi = 600, units = "in", device='tiff')

plot_right_amygdala_dki_kfa_age_int <- interactions::interact_plot(glm_right_amygdala_dki_kfa_age_int, pred = age, modx = Group, point.alpha = 1,
                                                          plot.points = T, interval = T, point.size = 0.1, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "Right Amygdala KFA", x = "Age (Years)",
    subtitle = paste0(
                      symnum(summary(glm_right_amygdala_dki_kfa_age_int)$coefficients[4,4], corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " interaction p: ", format.pval(summary(glm_right_amygdala_dki_kfa_age_int)$coefficients[4,4], digits = 2, eps = 0.001))) +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_markdown())

ggsave(filename = "kfa_age_int_plot.tif", plot_right_amygdala_dki_kfa_age_int,
       width = 4, height = 3, dpi = 600, units = "in", device = "tiff")

right_amygdala_ctl <- right_amygdala %>% filter(Group == "Control")
right_amygdala_scd <- right_amygdala %>% filter(Group == "SCD")

glm_right_amygdala_dki_kfa_age_scd <- lm(dki_kfa ~ age, right_amygdala_scd)
summary(glm_right_amygdala_dki_kfa_age_scd)
tbl_regression(glm_right_amygdala_dki_kfa_age_scd, 
               estimate_fun = label_style_sigfig(digits = 3)) %>% 
  bold_p(t=0.025) %>%
  add_glance_table(include = c(r.squared), 
                   label = list(p.value = "overall p-value")) %>%
  modify_spanning_header(c(estimate, conf.low, conf.high, p.value) ~ "**SCD**")

glm_right_amygdala_dki_kfa_age_ctl <- lm(dki_kfa ~ age, right_amygdala_ctl)
summary(glm_right_amygdala_dki_kfa_age_ctl)
tbl_regression(glm_right_amygdala_dki_kfa_age_ctl, 
               estimate_fun = label_style_sigfig(digits = 3)) %>% 
  bold_p(t=0.025) %>%
  add_glance_table(include = c(r.squared), 
                   label = list(p.value = "overall p-value")) %>%
  modify_spanning_header(c(estimate, conf.low, conf.high, p.value) ~ "**Control**")

glm_left_amygdala_dti_fa_additional_hads_anxiety_int <- lm(dti_fa ~ additional_hads_anxiety * Group, left_amygdala)
summary(glm_left_amygdala_dti_fa_additional_hads_anxiety_int)
glm_left_amygdala_dti_fa_additional_hads_anxiety_int <- lm(dti_fa ~ additional_hads_anxiety * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_dti_fa_additional_hads_anxiety_int)
glm_left_amygdala_fit_FWF_additional_hads_anxiety_int <- lm(fit_FWF ~ additional_hads_anxiety * Group, left_amygdala)
summary(glm_left_amygdala_fit_FWF_additional_hads_anxiety_int)
glm_left_amygdala_fit_FWF_additional_hads_anxiety_int <- lm(fit_FWF ~ additional_hads_anxiety * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_fit_FWF_additional_hads_anxiety_int)
glm_left_amygdala_dki_kfa_additional_hads_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * Group, left_amygdala)
summary(glm_left_amygdala_dki_kfa_additional_hads_anxiety_int)
glm_left_amygdala_dki_kfa_additional_hads_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_dki_kfa_additional_hads_anxiety_int)
glm_right_amygdala_dki_kfa_additional_hads_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * Group, right_amygdala)
summary(glm_right_amygdala_dki_kfa_additional_hads_anxiety_int)
glm_right_amygdala_dki_kfa_additional_hads_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * Group, right_amygdala_10mq_thresh)
summary(glm_right_amygdala_dki_kfa_additional_hads_anxiety_int)

plot_left_amygdala_dti_fa_additional_hads_anxiety_int <- interactions::interact_plot(glm_left_amygdala_dti_fa_additional_hads_anxiety_int, pred = additional_hads_anxiety, modx = Group, point.alpha = 1,
                                                                   plot.points = T, interval = T, point.size = 0.1, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "Left Amygdala FA", x = "Anxiety (HADS)",
    subtitle = paste0(
      symnum(summary(glm_left_amygdala_dti_fa_additional_hads_anxiety_int)$coefficients[4,4], corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
      " interaction p: ", format.pval(summary(glm_left_amygdala_dti_fa_additional_hads_anxiety_int)$coefficients[4,4], digits = 2, eps = 0.001))) +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_markdown())

ggsave(filename = "fa_anxiety_int_plot.tif", plot_left_amygdala_dti_fa_additional_hads_anxiety_int,
       width = 4, height = 3, dpi = 600, units = "in", device = "tiff")

tbl_regression(glm_left_amygdala_dti_fa_additional_hads_anxiety_int, 
               estimate_fun = label_style_sigfig(digits = 3),
               show_single_row = c(Group, `additional_hads_anxiety:Group`)) %>% 
  bold_p(t=0.05) %>%
  add_glance_table(include = c(r.squared, p.value), 
                   label = list(p.value = "overall p-value")) %>%
  modify_spanning_header(c(estimate, conf.low, conf.high, p.value) ~ "**Left Amygdala FA**") #%>%
# as_gt() %>%
# gt::gtsave(filename = "right_amygdala_dki_kfa_age_int.docx")

glm_left_amygdala_dti_fa_additional_hads_depression_int <- lm(dti_fa ~ additional_hads_depression * Group, left_amygdala)
summary(glm_left_amygdala_dti_fa_additional_hads_depression_int)
glm_left_amygdala_dti_fa_additional_hads_depression_int <- lm(dti_fa ~ additional_hads_depression * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_dti_fa_additional_hads_depression_int)
glm_left_amygdala_fit_FWF_additional_hads_depression_int <- lm(fit_FWF ~ additional_hads_depression * Group, left_amygdala)
summary(glm_left_amygdala_fit_FWF_additional_hads_depression_int)
glm_left_amygdala_fit_FWF_additional_hads_depression_int <- lm(fit_FWF ~ additional_hads_depression * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_fit_FWF_additional_hads_depression_int)
glm_left_amygdala_dki_kfa_additional_hads_depression_int <- lm(dki_kfa ~ additional_hads_depression * Group, left_amygdala)
summary(glm_left_amygdala_dki_kfa_additional_hads_depression_int)
glm_left_amygdala_dki_kfa_additional_hads_depression_int <- lm(dki_kfa ~ additional_hads_depression * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_dki_kfa_additional_hads_depression_int)
glm_right_amygdala_dki_kfa_additional_hads_depression_int <- lm(dki_kfa ~ additional_hads_depression * Group, right_amygdala)
summary(glm_right_amygdala_dki_kfa_additional_hads_depression_int)
glm_right_amygdala_dki_kfa_additional_hads_depression_int <- lm(dki_kfa ~ additional_hads_depression * Group, right_amygdala_10mq_thresh)
summary(glm_right_amygdala_dki_kfa_additional_hads_depression_int)

glm_left_amygdala_dti_fa_homeint_storyrecall_d_int <- lm(dti_fa ~ homeint_storyrecall_d * Group, left_amygdala)
summary(glm_left_amygdala_dti_fa_homeint_storyrecall_d_int)
glm_left_amygdala_dti_fa_homeint_storyrecall_d_int <- lm(dti_fa ~ homeint_storyrecall_d * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_dti_fa_homeint_storyrecall_d_int)
glm_left_amygdala_fit_FWF_homeint_storyrecall_d_int <- lm(fit_FWF ~ homeint_storyrecall_d * Group, left_amygdala)
summary(glm_left_amygdala_fit_FWF_homeint_storyrecall_d_int)
glm_left_amygdala_fit_FWF_homeint_storyrecall_d_int <- lm(fit_FWF ~ homeint_storyrecall_d * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_fit_FWF_homeint_storyrecall_d_int)
glm_left_amygdala_dki_kfa_homeint_storyrecall_d_int <- lm(dki_kfa ~ homeint_storyrecall_d * Group, left_amygdala)
summary(glm_left_amygdala_dki_kfa_homeint_storyrecall_d_int)
glm_left_amygdala_dki_kfa_homeint_storyrecall_d_int <- lm(dki_kfa ~ homeint_storyrecall_d * Group, left_amygdala_10mq_thresh)
summary(glm_left_amygdala_dki_kfa_homeint_storyrecall_d_int)
glm_right_amygdala_dki_kfa_homeint_storyrecall_d_int <- lm(dki_kfa ~ homeint_storyrecall_d * Group, right_amygdala)
summary(glm_right_amygdala_dki_kfa_homeint_storyrecall_d_int)
glm_right_amygdala_dki_kfa_homeint_storyrecall_d_int <- lm(dki_kfa ~ homeint_storyrecall_d * Group, right_amygdala_10mq_thresh)
summary(glm_right_amygdala_dki_kfa_homeint_storyrecall_d_int)

summary(lm(dti_fa ~ ten_mq_total, left_amygdala))
summary(lm(dti_md ~ ten_mq_total, left_amygdala))
summary(lm(fit_FWF ~ ten_mq_total, left_amygdala))
summary(lm(dki_kfa ~ ten_mq_total, left_amygdala))
summary(lm(dki_mk ~ ten_mq_total, left_amygdala))
summary(lm(dki_kfa ~ ten_mq_total, right_amygdala))
summary(lm(dki_mk ~ ten_mq_total, right_amygdala))

glm_left_amygdala_dti_fa_objmem_emotion_enhance_int <- lm(dti_fa ~ objmem_emotion_enhance * Group, left_amygdala)
summary(glm_left_amygdala_dti_fa_objmem_emotion_enhance_int)
glm_left_amygdala_fit_FWF_objmem_emotion_enhance_int <- lm(fit_FWF ~ objmem_emotion_enhance * Group, left_amygdala)
summary(glm_left_amygdala_fit_FWF_objmem_emotion_enhance_int)
glm_left_amygdala_dki_kfa_objmem_emotion_enhance_int <- lm(dki_kfa ~ objmem_emotion_enhance * Group, left_amygdala)
summary(glm_left_amygdala_dki_kfa_objmem_emotion_enhance_int)
glm_right_amygdala_dki_kfa_objmem_emotion_enhance_int <- lm(dki_kfa ~ objmem_emotion_enhance * Group, right_amygdala)
summary(glm_right_amygdala_dki_kfa_objmem_emotion_enhance_int)

glm_left_amygdala_dti_fa_objmem_pos_bias_int <- lm(dti_fa ~ objmem_pos_bias * Group, left_amygdala)
summary(glm_left_amygdala_dti_fa_objmem_pos_bias_int)
glm_left_amygdala_fit_FWF_objmem_pos_bias_int <- lm(fit_FWF ~ objmem_pos_bias * Group, left_amygdala)
summary(glm_left_amygdala_fit_FWF_objmem_pos_bias_int)
glm_left_amygdala_dki_kfa_objmem_pos_bias_int <- lm(dki_kfa ~ objmem_pos_bias * Group, left_amygdala)
summary(glm_left_amygdala_dki_kfa_objmem_pos_bias_int)
glm_right_amygdala_dki_kfa_objmem_pos_bias_int <- lm(dki_kfa ~ objmem_pos_bias * Group, right_amygdala)
summary(glm_right_amygdala_dki_kfa_objmem_pos_bias_int)

plot_right_amygdala_dki_kfa_emo_mem_bias_int <- interactions::interact_plot(glm_right_amygdala_dki_kfa_objmem_pos_bias_int, pred = objmem_pos_bias, modx = Group, point.alpha = 1,
                            plot.points = T, interval = T, point.size = 0.1, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "Right Amygdala KFA", x = "Positivity Memory Bias",
    subtitle = paste0(
      symnum(summary(glm_right_amygdala_dki_kfa_objmem_pos_bias_int)$coefficients[4,4], corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
      " interaction p: ", format.pval(summary(glm_right_amygdala_dki_kfa_objmem_pos_bias_int)$coefficients[4,4], digits = 2, eps = 0.001))) +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_markdown())

ggsave(filename = "kfa_emo_mem_bias_int_plot.tif", plot_right_amygdala_dki_kfa_emo_mem_bias_int,
       width = 5, height = 3, dpi = 600, units = "in", device = "tiff")

right_amygdala_ctl <- right_amygdala %>% filter(Group == "Control")
right_amygdala_scd <- right_amygdala %>% filter(Group == "SCD")

glm_right_amygdala_dki_kfa_objmem_pos_bias_scd <- lm(dki_kfa ~ objmem_pos_bias, right_amygdala_scd)
summary(glm_right_amygdala_dki_kfa_objmem_pos_bias_scd)
tbl_regression(glm_right_amygdala_dki_kfa_objmem_pos_bias_scd, 
               estimate_fun = label_style_sigfig(digits = 3)) %>% 
  bold_p(t=0.025) %>%
  add_glance_table(include = c(r.squared), 
                   label = list(p.value = "overall p-value")) %>%
  modify_spanning_header(c(estimate, conf.low, conf.high, p.value) ~ "**SCD**")

glm_right_amygdala_dki_kfa_objmem_pos_bias_ctl <- lm(dki_kfa ~ objmem_pos_bias, right_amygdala_ctl)
summary(glm_right_amygdala_dki_kfa_objmem_pos_bias_ctl)
tbl_regression(glm_right_amygdala_dki_kfa_objmem_pos_bias_ctl, 
               estimate_fun = label_style_sigfig(digits = 3)) %>% 
  bold_p(t=0.025) %>%
  add_glance_table(include = c(r.squared), 
                   label = list(p.value = "overall p-value")) %>%
  modify_spanning_header(c(estimate, conf.low, conf.high, p.value) ~ "**Control**")


#shapiro-wilk test for normality. if violated, can argue data is sufficiently large
#so violation is trivial or transform dependent variable (ie box-cox transform or log transform)
#homogeneity assumption plot residual and y-hat. can use square root or log transform if violated
#levene's test, aka equal variance
#proc transreg procedure: explore different transformation methods

library(car)
leveneTest(dti_fa ~ Group, left_amygdala)
leveneTest(fit_FWF ~ Group, left_amygdala)
leveneTest(dki_kfa ~ Group, left_amygdala)
leveneTest(dki_kfa ~ Group, right_amygdala)
