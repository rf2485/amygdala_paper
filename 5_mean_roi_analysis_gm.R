source("1_data_preparation_gm.R")
library(tidyverse)
library(labelVector)
library(gtsummary)
library(ggeffects)
library(ggtext)
library(patchwork)

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

mtr_tr30 <- aseg2mtr_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  filter(mt_tr=="TR=30ms")

mtr_tr50 <- aseg2mtr_gm %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  filter(mt_tr=="TR=50ms")


#### focus on amygdala ####
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
# left_amygdala <- mtr_tr30 %>% dplyr::select(participant_id, Left.Amygdala) %>%
#   rename(mtr_tr30 = "Left.Amygdala") %>% full_join(left_amygdala, .)
# left_amygdala <- mtr_tr50 %>% dplyr::select(participant_id, Left.Amygdala) %>%
#   rename(mtr_tr50 = "Left.Amygdala") %>% full_join(left_amygdala, .)
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
                           # mtr_tr30 = "MTR TR=30ms",
                           # mtr_tr50 = "MTR TR=50ms"
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
# right_amygdala <- mtr_tr30 %>% dplyr::select(participant_id, Right.Amygdala) %>%
#   rename(mtr_tr30 = "Right.Amygdala") %>% full_join(right_amygdala, .)
# right_amygdala <- mtr_tr50 %>% dplyr::select(participant_id, Right.Amygdala) %>%
#   rename(mtr_tr50 = "Right.Amygdala") %>% full_join(right_amygdala, .)

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
                           # mtr_tr30 = "MTR TR=30ms",
                           # mtr_tr50 = "MTR TR=50ms"
)


left_amygdala %>% 
  select(-participant_id) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(dti_fa, dti_md, dki_kfa, dki_mk, 
                           volume, fit_FWF, fit_NDI, fit_ODI)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p() %>% bold_p(q=T) %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Left Amygdala </div>")

right_amygdala %>% 
  select(-participant_id) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(dti_fa, dti_md, dki_kfa, dki_mk, 
                          volume, fit_FWF, fit_NDI, fit_ODI)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p() %>% bold_p(q=T) %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**") %>%
  modify_caption("<div style='text-align: left; font-weight: bold;'> Right Amygdala </div>")
#### ROC and AUC ####
library(pROC)
pROC::roc.test(response = right_amygdala$SCD, predictor1 = right_amygdala$volume, 
               predictor2 = right_amygdala$dki_kfa, levels = c("Control", "SCD"))

rocobj1 <- plot.roc(right_amygdala$SCD, right_amygdala$volume, print.auc = T,
                    
                    main="AUC Volume vs KFA, Right Amygdala", percent=TRUE, col="#1c61b6")

rocobj2 <- plot.roc(right_amygdala$SCD, right_amygdala$dki_kfa, print.auc = T, print.auc.adj = c(0, 3), add = T,
                     percent=TRUE, col="#008600")

testobj <- roc.test(rocobj1, rocobj2)

text(20, 45, labels=paste("p-value =", format.pval(testobj$p.value, digits = 2, eps = 0.001)), adj=c(0, .5))

legend("bottomright", legend=c("Volume", "KFA"), col=c("#1c61b6", "#008600"), lwd=2)

##### Scatterplots of Linear Models ######
left_amygdala_ctl <- left_amygdala %>% filter(SCD == "Control")
left_amygdala_scd <- left_amygdala %>% filter(SCD == "SCD")
right_amygdala_ctl <- right_amygdala %>% filter(SCD == "Control")
right_amygdala_scd <- right_amygdala %>% filter(SCD == "SCD")

### DTI FA correlations with bioloigcally interpretable imaging metrics
glm_left_amygdala_dti_fa_volume_scd <- lm(dti_fa ~ volume, left_amygdala_scd)
summary(glm_left_amygdala_dti_fa_volume_scd)
corr_p_value <- summary(glm_left_amygdala_dti_fa_volume_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_dti_fa_volume_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_dti_fa_volume_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_dti_fa_volume_scd, c( "volume[all]"))
plot_glm_left_amygdala_dti_fa_volume_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    x = expression(Volume~(mm^3)),
    title = "Left Amygdala Volume",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_left_amygdala_dti_fa_fit_NDI_scd <- lm(dti_fa ~ fit_NDI, left_amygdala_scd)
summary(glm_left_amygdala_dti_fa_fit_NDI_scd)
corr_p_value <- summary(glm_left_amygdala_dti_fa_fit_NDI_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_dti_fa_fit_NDI_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_dti_fa_fit_NDI_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_dti_fa_fit_NDI_scd, c( "fit_NDI[all]"))
plot_glm_left_amygdala_dti_fa_fit_NDI_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Left Amygdala Mean NDI",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_left_amygdala_dti_fa_fit_ODI_scd <- lm(dti_fa ~ fit_ODI, left_amygdala_scd)
summary(glm_left_amygdala_dti_fa_fit_ODI_scd)
corr_p_value <- summary(glm_left_amygdala_dti_fa_fit_ODI_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_dti_fa_fit_ODI_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_dti_fa_fit_ODI_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_dti_fa_fit_ODI_scd, c( "fit_ODI[all]"))
plot_glm_left_amygdala_dti_fa_fit_ODI_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Left Amygdala Mean ODI",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_left_amygdala_dti_fa_fit_FWF_scd <- lm(dti_fa ~ fit_FWF, left_amygdala_scd)
summary(glm_left_amygdala_dti_fa_fit_FWF_scd)
corr_p_value <- summary(glm_left_amygdala_dti_fa_fit_FWF_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_dti_fa_fit_FWF_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_dti_fa_fit_FWF_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_dti_fa_fit_FWF_scd, c( "fit_FWF[all]"))
plot_glm_left_amygdala_dti_fa_fit_FWF_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Left Amygdala Mean FWF",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
left_amygdala_dti_fa_bio_plots <- plot_glm_left_amygdala_dti_fa_volume_scd + 
  plot_glm_left_amygdala_dti_fa_fit_NDI_scd + 
  plot_glm_left_amygdala_dti_fa_fit_ODI_scd + 
  plot_glm_left_amygdala_dti_fa_fit_FWF_scd +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(title = "Correlations between Left Amygdala FA and Biologically Interpretable 
  \nImaging Metrics within the SCD group. * p < 0.05, ** p < 0.01, *** p < 0.001",
                  tag_levels = "A")
ggsave(filename = "left_amygdala_dti_fa_bio_plots.tif", left_amygdala_dti_fa_bio_plots,
       width = 6.5, height = 5.5, dpi = 600, units = "in", device='tiff')

glm_left_amygdala_dki_kfa_volume_scd <- lm(dki_kfa ~ volume, left_amygdala_scd)
summary(glm_left_amygdala_dki_kfa_volume_scd)
corr_p_value <- summary(glm_left_amygdala_dki_kfa_volume_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_dki_kfa_volume_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_dki_kfa_volume_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_dki_kfa_volume_scd, c( "volume[all]"))
plot_glm_left_amygdala_dki_kfa_volume_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    x = expression(Volume~(mm^3)),
    title = "Left Amygdala Volume",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_left_amygdala_dki_kfa_fit_NDI_scd <- lm(dki_kfa ~ fit_NDI, left_amygdala_scd)
summary(glm_left_amygdala_dki_kfa_fit_NDI_scd)
corr_p_value <- summary(glm_left_amygdala_dki_kfa_fit_NDI_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_dki_kfa_fit_NDI_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_dki_kfa_fit_NDI_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_dki_kfa_fit_NDI_scd, c( "fit_NDI[all]"))
plot_glm_left_amygdala_dki_kfa_fit_NDI_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Left Amygdala Mean NDI",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_left_amygdala_dki_kfa_fit_ODI_scd <- lm(dki_kfa ~ fit_ODI, left_amygdala_scd)
summary(glm_left_amygdala_dki_kfa_fit_ODI_scd)
corr_p_value <- summary(glm_left_amygdala_dki_kfa_fit_ODI_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_dki_kfa_fit_ODI_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_dki_kfa_fit_ODI_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_dki_kfa_fit_ODI_scd, c( "fit_ODI[all]"))
plot_glm_left_amygdala_dki_kfa_fit_ODI_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Left Amygdala Mean ODI",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_left_amygdala_dki_kfa_fit_FWF_scd <- lm(dki_kfa ~ fit_FWF, left_amygdala_scd)
summary(glm_left_amygdala_dki_kfa_fit_FWF_scd)
corr_p_value <- summary(glm_left_amygdala_dki_kfa_fit_FWF_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_dki_kfa_fit_FWF_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_dki_kfa_fit_FWF_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_dki_kfa_fit_FWF_scd, c( "fit_FWF[all]"))
plot_glm_left_amygdala_dki_kfa_fit_FWF_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Left Amygdala Mean FWF",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
left_amygdala_dki_kfa_bio_plots <- plot_glm_left_amygdala_dki_kfa_volume_scd + 
  plot_glm_left_amygdala_dki_kfa_fit_NDI_scd + 
  plot_glm_left_amygdala_dki_kfa_fit_ODI_scd + 
  plot_glm_left_amygdala_dki_kfa_fit_FWF_scd +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(title = "Correlations between Left Amygdala KFA and Biologically Interpretable 
  \nImaging Metrics within the SCD group. * p < 0.05, ** p < 0.01, *** p < 0.001",
                  tag_levels = "A")
ggsave(filename = "left_amygdala_dki_kfa_bio_plots.tif", left_amygdala_dki_kfa_bio_plots,
       width = 6.5, height = 5.5, dpi = 600, units = "in", device='tiff')

glm_right_amygdala_dki_kfa_volume_scd <- lm(dki_kfa ~ volume, right_amygdala_scd)
summary(glm_right_amygdala_dki_kfa_volume_scd)
corr_p_value <- summary(glm_right_amygdala_dki_kfa_volume_scd)$coefficients[2,4]
corr_beta <- summary(glm_right_amygdala_dki_kfa_volume_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_right_amygdala_dki_kfa_volume_scd)$adj.r.squared
pr <- predict_response(glm_right_amygdala_dki_kfa_volume_scd, c( "volume[all]"))
plot_glm_right_amygdala_dki_kfa_volume_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    x = expression(Volume~(mm^3)),
    title = "Right Amygdala Volume",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_right_amygdala_dki_kfa_fit_NDI_scd <- lm(dki_kfa ~ fit_NDI, right_amygdala_scd)
summary(glm_right_amygdala_dki_kfa_fit_NDI_scd)
corr_p_value <- summary(glm_right_amygdala_dki_kfa_fit_NDI_scd)$coefficients[2,4]
corr_beta <- summary(glm_right_amygdala_dki_kfa_fit_NDI_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_right_amygdala_dki_kfa_fit_NDI_scd)$adj.r.squared
pr <- predict_response(glm_right_amygdala_dki_kfa_fit_NDI_scd, c( "fit_NDI[all]"))
plot_glm_right_amygdala_dki_kfa_fit_NDI_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Right Amygdala Mean NDI",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_right_amygdala_dki_kfa_fit_ODI_scd <- lm(dki_kfa ~ fit_ODI, right_amygdala_scd)
summary(glm_right_amygdala_dki_kfa_fit_ODI_scd)
corr_p_value <- summary(glm_right_amygdala_dki_kfa_fit_ODI_scd)$coefficients[2,4]
corr_beta <- summary(glm_right_amygdala_dki_kfa_fit_ODI_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_right_amygdala_dki_kfa_fit_ODI_scd)$adj.r.squared
pr <- predict_response(glm_right_amygdala_dki_kfa_fit_ODI_scd, c( "fit_ODI[all]"))
plot_glm_right_amygdala_dki_kfa_fit_ODI_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Right Amygdala Mean ODI",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_right_amygdala_dki_kfa_fit_FWF_scd <- lm(dki_kfa ~ fit_FWF, right_amygdala_scd)
summary(glm_right_amygdala_dki_kfa_fit_FWF_scd)
corr_p_value <- summary(glm_right_amygdala_dki_kfa_fit_FWF_scd)$coefficients[2,4]
corr_beta <- summary(glm_right_amygdala_dki_kfa_fit_FWF_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_right_amygdala_dki_kfa_fit_FWF_scd)$adj.r.squared
pr <- predict_response(glm_right_amygdala_dki_kfa_fit_FWF_scd, c( "fit_FWF[all]"))
plot_glm_right_amygdala_dki_kfa_fit_FWF_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Right Amygdala Mean FWF",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
right_amygdala_dki_kfa_bio_plots <- plot_glm_right_amygdala_dki_kfa_volume_scd + 
  plot_glm_right_amygdala_dki_kfa_fit_NDI_scd + 
  plot_glm_right_amygdala_dki_kfa_fit_ODI_scd + 
  plot_glm_right_amygdala_dki_kfa_fit_FWF_scd +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(title = "Correlations between Right Amygdala KFA and Biologically Interpretable 
  \nImaging Metrics within the SCD group. * p < 0.05, ** p < 0.01, *** p < 0.001",
                  tag_levels = "A")
ggsave(filename = "right_amygdala_dki_kfa_bio_plots.tif", right_amygdala_dki_kfa_bio_plots,
       width = 6.5, height = 5.5, dpi = 600, units = "in", device='tiff')

glm_left_amygdala_dti_fa_age_int <- lm(dti_fa ~ age * SCD, left_amygdala)
summary(glm_left_amygdala_dti_fa_age_int)
glm_left_amygdala_dki_kfa_age_int <- lm(dki_kfa ~ age * SCD, left_amygdala)
summary(glm_left_amygdala_dki_kfa_age_int)
glm_left_amygdala_fit_FWF_age_int <- lm(fit_FWF ~ age * SCD, left_amygdala)
summary(glm_left_amygdala_fit_FWF_age_int)
glm_right_amygdala_dki_kfa_age_int <- lm(dki_kfa ~ age * SCD, right_amygdala)
summary(glm_right_amygdala_dki_kfa_age_int)

glm_left_amygdala_dti_fa_additional_hads_anxiety_int <- lm(dti_fa ~ additional_hads_anxiety * SCD, left_amygdala)
summary(glm_left_amygdala_dti_fa_additional_hads_anxiety_int)
glm_left_amygdala_dki_kfa_additional_hads_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * SCD, left_amygdala)
summary(glm_left_amygdala_dki_kfa_additional_hads_anxiety_int)
glm_left_amygdala_fit_FWF_additional_hads_anxiety_int <- lm(fit_FWF ~ additional_hads_anxiety * SCD, left_amygdala)
summary(glm_left_amygdala_fit_FWF_additional_hads_anxiety_int)
glm_right_amygdala_dki_kfa_additional_hads_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * SCD, right_amygdala)
summary(glm_right_amygdala_dki_kfa_additional_hads_anxiety_int)

###age interaction
glm_left_amygdala_dti_fa_age_int <- lm(dti_fa ~ age * SCD, left_amygdala)
summary(glm_left_amygdala_dti_fa_age_int)
across_group_p_value <- summary(glm_left_amygdala_dti_fa_age_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_dti_fa_age_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_dti_fa_age_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_dti_fa_age_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_dti_fa_age_int)$coefficients[4,1]                
glm_left_amygdala_dti_fa_age_scd <- lm(dti_fa ~ age, left_amygdala_scd)
summary(glm_left_amygdala_dti_fa_age_scd)
scd_p_value <- summary(glm_left_amygdala_dti_fa_age_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_dti_fa_age_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_dti_fa_age_scd)$adj.r.squared
glm_left_amygdala_dti_fa_age_ctl <- lm(dti_fa ~ age, left_amygdala_ctl)
summary(glm_left_amygdala_dti_fa_age_ctl)
ctl_p_value <- summary(glm_left_amygdala_dti_fa_age_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_dti_fa_age_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_dti_fa_age_ctl)$adj.r.squared
plot_left_amygdala_dti_fa <- interactions::interact_plot(glm_left_amygdala_dti_fa_age_int, pred = age, modx = SCD, 
                            plot.points = T, interval = T, point.size = 0.5,
                            legend.main = 'Cohort') +
  theme(legend.position = 'none') +
  labs(
    y = "FA", x = "Age (Years)", title = "Left Amygdala Mean FA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
                      )
    ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                    label = paste0(
                      symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(scd_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(scd_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(scd_adj_r_squared, 2)
                      ),
                    color = "SCD"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
                ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 2.5, hjust = 1.01,
                    label = paste0(
                      symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(ctl_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(ctl_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(ctl_adj_r_squared, 2)
                    ),
                    color = "Control"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
                ) +
  ylim(NA,0.27) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

glm_left_amygdala_fit_FWF_age_int <- lm(fit_FWF ~ age * SCD, left_amygdala)
summary(glm_left_amygdala_fit_FWF_age_int)
across_group_p_value <- summary(glm_left_amygdala_fit_FWF_age_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_fit_FWF_age_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_age_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_fit_FWF_age_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_fit_FWF_age_int)$coefficients[4,1]                
glm_left_amygdala_fit_FWF_age_scd <- lm(fit_FWF ~ age, left_amygdala_scd)
summary(glm_left_amygdala_fit_FWF_age_scd)
scd_p_value <- summary(glm_left_amygdala_fit_FWF_age_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_fit_FWF_age_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_age_scd)$adj.r.squared
glm_left_amygdala_fit_FWF_age_ctl <- lm(fit_FWF ~ age, left_amygdala_ctl)
summary(glm_left_amygdala_fit_FWF_age_ctl)
ctl_p_value <- summary(glm_left_amygdala_fit_FWF_age_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_fit_FWF_age_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_age_ctl)$adj.r.squared
plot_left_amygdala_fit_FWF <- interactions::interact_plot(glm_left_amygdala_fit_FWF_age_int, pred = age, modx = SCD, 
                                                         plot.points = T, interval = T, point.size = 0.5,
                                                         legend.main = 'Cohort') +
  theme(legend.position = 'none') +
  labs(
    y = "FWF", x = "Age (Years)", title = "Left Amygdala Mean FWF",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                    label = paste0(
                      symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(scd_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(scd_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(scd_adj_r_squared, 2)
                    ),
                    color = "SCD"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 2.5, hjust = 1.01,
                    label = paste0(
                      symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(ctl_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(ctl_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(ctl_adj_r_squared, 2)
                    ),
                    color = "Control"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  ylim(NA,0.23) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

glm_left_amygdala_dki_kfa_age_int <- lm(dki_kfa ~ age * SCD, left_amygdala)
summary(glm_left_amygdala_dki_kfa_age_int)
across_group_p_value <- summary(glm_left_amygdala_dki_kfa_age_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_dki_kfa_age_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_age_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_dki_kfa_age_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_dki_kfa_age_int)$coefficients[4,1]                
glm_left_amygdala_dki_kfa_age_scd <- lm(dki_kfa ~ age, left_amygdala_scd)
summary(glm_left_amygdala_dki_kfa_age_scd)
scd_p_value <- summary(glm_left_amygdala_dki_kfa_age_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_dki_kfa_age_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_age_scd)$adj.r.squared
glm_left_amygdala_dki_kfa_age_ctl <- lm(dki_kfa ~ age, left_amygdala_ctl)
summary(glm_left_amygdala_dki_kfa_age_ctl)
ctl_p_value <- summary(glm_left_amygdala_dki_kfa_age_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_dki_kfa_age_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_age_ctl)$adj.r.squared
plot_left_amygdala_dki_kfa <- interactions::interact_plot(glm_left_amygdala_dki_kfa_age_int, pred = age, modx = SCD, 
                                                          plot.points = T, interval = T, point.size = 0.5,
                                                          legend.main = 'Cohort') +
  theme(legend.position = 'none') +
  labs(
    y = "KFA", x = "Age (Years)", title = "Left Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                    label = paste0(
                      symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(scd_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(scd_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(scd_adj_r_squared, 2)
                    ),
                    color = "SCD"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 2.5, hjust = 1.01,
                    label = paste0(
                      symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(ctl_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(ctl_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(ctl_adj_r_squared, 2)
                    ),
                    color = "Control"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  ylim(0.2, 0.6) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

glm_right_amygdala_dki_kfa_age_int <- lm(dki_kfa ~ age * SCD, right_amygdala)
summary(glm_right_amygdala_dki_kfa_age_int)
across_group_p_value <- summary(glm_right_amygdala_dki_kfa_age_int)$coefficients[2,4]
across_group_beta <- summary(glm_right_amygdala_dki_kfa_age_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_age_int)$adj.r.squared
int_p_value <- summary(glm_right_amygdala_dki_kfa_age_int)$coefficients[4,4]
int_beta <-  summary(glm_right_amygdala_dki_kfa_age_int)$coefficients[4,1]                
glm_right_amygdala_dki_kfa_age_scd <- lm(dki_kfa ~ age, right_amygdala_scd)
summary(glm_right_amygdala_dki_kfa_age_scd)
scd_p_value <- summary(glm_right_amygdala_dki_kfa_age_scd)$coefficients[2,4]
scd_beta <- summary(glm_right_amygdala_dki_kfa_age_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_age_scd)$adj.r.squared
glm_right_amygdala_dki_kfa_age_ctl <- lm(dki_kfa ~ age, right_amygdala_ctl)
summary(glm_right_amygdala_dki_kfa_age_ctl)
ctl_p_value <- summary(glm_right_amygdala_dki_kfa_age_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_right_amygdala_dki_kfa_age_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_age_ctl)$adj.r.squared
plot_right_amygdala_dki_kfa <- interactions::interact_plot(glm_right_amygdala_dki_kfa_age_int, pred = age, modx = SCD, 
                                                          plot.points = T, interval = T, point.size = 0.5,
                                                          legend.main = 'Cohort') +
  theme(legend.position = 'none') +
  labs(
    y = "KFA", x = "Age (Years)", title = "Right Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                    label = paste0(
                      symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(scd_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(scd_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(scd_adj_r_squared, 2)
                    ),
                    color = "SCD"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 2.5, hjust = 1.01,
                    label = paste0(
                      symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(ctl_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(ctl_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(ctl_adj_r_squared, 2)
                    ),
                    color = "Control"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  ylim(0.2, 0.6) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

diffusion_age_plots <- plot_left_amygdala_dti_fa + 
  plot_left_amygdala_fit_FWF + 
  plot_left_amygdala_dki_kfa + 
  plot_right_amygdala_dki_kfa +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(title = "Correlations Between Age and Mean Amygdala Diffusion Metrics. \n* p < 0.05, ** p < 0.01, *** p < 0.001",
                  tag_levels = "A")

ggsave(filename = "diffusion_age_plots.tif", diffusion_age_plots,
       width = 7.5, height = 7.5, dpi = 600, units = "in", device='tiff')

###anxiety interaction
glm_left_amygdala_dti_fa_anxiety_int <- lm(dti_fa ~ additional_hads_anxiety * SCD, left_amygdala)
summary(glm_left_amygdala_dti_fa_anxiety_int)
across_group_p_value <- summary(glm_left_amygdala_dti_fa_anxiety_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_dti_fa_anxiety_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_dti_fa_anxiety_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_dti_fa_anxiety_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_dti_fa_anxiety_int)$coefficients[4,1]                
glm_left_amygdala_dti_fa_anxiety_scd <- lm(dti_fa ~ additional_hads_anxiety, left_amygdala_scd)
summary(glm_left_amygdala_dti_fa_anxiety_scd)
scd_p_value <- summary(glm_left_amygdala_dti_fa_anxiety_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_dti_fa_anxiety_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_dti_fa_anxiety_scd)$adj.r.squared
glm_left_amygdala_dti_fa_anxiety_ctl <- lm(dti_fa ~ additional_hads_anxiety, left_amygdala_ctl)
summary(glm_left_amygdala_dti_fa_anxiety_ctl)
ctl_p_value <- summary(glm_left_amygdala_dti_fa_anxiety_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_dti_fa_anxiety_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_dti_fa_anxiety_ctl)$adj.r.squared
plot_left_amygdala_dti_fa <- interactions::interact_plot(glm_left_amygdala_dti_fa_anxiety_int, pred = additional_hads_anxiety, modx = SCD, 
                                                         plot.points = T, interval = T, point.size = 0.5,
                                                         legend.main = 'Cohort') +
  theme(legend.position = 'none') +
  labs(
    y = "FA", x = "Anxiety", title = "Left Amygdala Mean FA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                    label = paste0(
                      symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(scd_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(scd_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(scd_adj_r_squared, 2)
                    ),
                    color = "SCD"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 2.5, hjust = 1.01,
                    label = paste0(
                      symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(ctl_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(ctl_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(ctl_adj_r_squared, 2)
                    ),
                    color = "Control"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  ylim(NA,0.27) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))
plot_left_amygdala_dti_fa
glm_left_amygdala_fit_FWF_anxiety_int <- lm(fit_FWF ~ additional_hads_anxiety * SCD, left_amygdala)
summary(glm_left_amygdala_fit_FWF_anxiety_int)
across_group_p_value <- summary(glm_left_amygdala_fit_FWF_anxiety_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_fit_FWF_anxiety_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_anxiety_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_fit_FWF_anxiety_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_fit_FWF_anxiety_int)$coefficients[4,1]                
glm_left_amygdala_fit_FWF_anxiety_scd <- lm(fit_FWF ~ additional_hads_anxiety, left_amygdala_scd)
summary(glm_left_amygdala_fit_FWF_anxiety_scd)
scd_p_value <- summary(glm_left_amygdala_fit_FWF_anxiety_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_fit_FWF_anxiety_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_anxiety_scd)$adj.r.squared
glm_left_amygdala_fit_FWF_anxiety_ctl <- lm(fit_FWF ~ additional_hads_anxiety, left_amygdala_ctl)
summary(glm_left_amygdala_fit_FWF_anxiety_ctl)
ctl_p_value <- summary(glm_left_amygdala_fit_FWF_anxiety_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_fit_FWF_anxiety_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_anxiety_ctl)$adj.r.squared
plot_left_amygdala_fit_FWF <- interactions::interact_plot(glm_left_amygdala_fit_FWF_anxiety_int, pred = additional_hads_anxiety, modx = SCD, 
                                                         plot.points = T, interval = T, point.size = 0.5,
                                                         legend.main = 'Cohort') +
  theme(legend.position = 'none') +
  labs(
    y = "FWF", x = "Anxiety", title = "Left Amygdala Mean FWF",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                    label = paste0(
                      symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(scd_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(scd_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(scd_adj_r_squared, 2)
                    ),
                    color = "SCD"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 2.5, hjust = 1.01,
                    label = paste0(
                      symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(ctl_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(ctl_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(ctl_adj_r_squared, 2)
                    ),
                    color = "Control"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  ylim(NA,0.27) +
  # scale_y_continuous(expand = expansion(mult = c(0.02,0.1))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

glm_left_amygdala_dki_kfa_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * SCD, left_amygdala)
summary(glm_left_amygdala_dki_kfa_anxiety_int)
across_group_p_value <- summary(glm_left_amygdala_dki_kfa_anxiety_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_dki_kfa_anxiety_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_anxiety_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_dki_kfa_anxiety_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_dki_kfa_anxiety_int)$coefficients[4,1]                
glm_left_amygdala_dki_kfa_anxiety_scd <- lm(dki_kfa ~ additional_hads_anxiety, left_amygdala_scd)
summary(glm_left_amygdala_dki_kfa_anxiety_scd)
scd_p_value <- summary(glm_left_amygdala_dki_kfa_anxiety_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_dki_kfa_anxiety_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_anxiety_scd)$adj.r.squared
glm_left_amygdala_dki_kfa_anxiety_ctl <- lm(dki_kfa ~ additional_hads_anxiety, left_amygdala_ctl)
summary(glm_left_amygdala_dki_kfa_anxiety_ctl)
ctl_p_value <- summary(glm_left_amygdala_dki_kfa_anxiety_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_dki_kfa_anxiety_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_anxiety_ctl)$adj.r.squared
plot_left_amygdala_dki_kfa <- interactions::interact_plot(glm_left_amygdala_dki_kfa_anxiety_int, pred = additional_hads_anxiety, modx = SCD, 
                                                          plot.points = T, interval = T, point.size = 0.5,
                                                          legend.main = 'Cohort') +
  theme(legend.position = 'none') +
  labs(
    y = "KFA", x = "Anxiety", title = "Left Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                    label = paste0(
                      symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(scd_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(scd_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(scd_adj_r_squared, 2)
                    ),
                    color = "SCD"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 2.5, hjust = 1.01,
                    label = paste0(
                      symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(ctl_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(ctl_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(ctl_adj_r_squared, 2)
                    ),
                    color = "Control"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  ylim(0.2, 0.6) +
  # scale_y_continuous(expand = expansion(mult = c(0.02,0.1))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))
glm_right_amygdala_dki_kfa_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * SCD, right_amygdala)
summary(glm_right_amygdala_dki_kfa_anxiety_int)
across_group_p_value <- summary(glm_right_amygdala_dki_kfa_anxiety_int)$coefficients[2,4]
across_group_beta <- summary(glm_right_amygdala_dki_kfa_anxiety_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_anxiety_int)$adj.r.squared
int_p_value <- summary(glm_right_amygdala_dki_kfa_anxiety_int)$coefficients[4,4]
int_beta <-  summary(glm_right_amygdala_dki_kfa_anxiety_int)$coefficients[4,1]                
glm_right_amygdala_dki_kfa_anxiety_scd <- lm(dki_kfa ~ additional_hads_anxiety, right_amygdala_scd)
summary(glm_right_amygdala_dki_kfa_anxiety_scd)
scd_p_value <- summary(glm_right_amygdala_dki_kfa_anxiety_scd)$coefficients[2,4]
scd_beta <- summary(glm_right_amygdala_dki_kfa_anxiety_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_anxiety_scd)$adj.r.squared
glm_right_amygdala_dki_kfa_anxiety_ctl <- lm(dki_kfa ~ additional_hads_anxiety, right_amygdala_ctl)
summary(glm_right_amygdala_dki_kfa_anxiety_ctl)
ctl_p_value <- summary(glm_right_amygdala_dki_kfa_anxiety_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_right_amygdala_dki_kfa_anxiety_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_anxiety_ctl)$adj.r.squared
plot_right_amygdala_dki_kfa <- interactions::interact_plot(glm_right_amygdala_dki_kfa_anxiety_int, pred = additional_hads_anxiety, modx = SCD, 
                                                          plot.points = T, interval = T, point.size = 0.5,
                                                          legend.main = 'Cohort') +
  theme(legend.position = 'none') +
  labs(
    y = "KFA", x = "Anxiety", title = "Right Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                    label = paste0(
                      symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(scd_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(scd_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(scd_adj_r_squared, 2)
                    ),
                    color = "SCD"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 2.5, hjust = 1.01,
                    label = paste0(
                      symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(ctl_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(ctl_beta, 2),
                      ", adj-R<sup>2</sup>: ", signif(ctl_adj_r_squared, 2)
                    ),
                    color = "Control"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  ylim(0.2, 0.6) +
  # scale_y_continuous(expand = expansion(mult = c(0.02,0.1))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

diffusion_anxiety_plots <- plot_left_amygdala_dti_fa + 
  plot_left_amygdala_fit_FWF + 
  plot_left_amygdala_dki_kfa + 
  plot_right_amygdala_dki_kfa +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(title = "Correlations Between Anxiety and Mean Amygdala Diffusion Metrics. \n* p < 0.05, ** p < 0.01, *** p < 0.001",
                  tag_levels = "A")

ggsave(filename = "diffusion_anxiety_plots.tif", diffusion_anxiety_plots,
       width = 7.5, height = 7.5, dpi = 600, units = "in", device='tiff')

##### Pearson and Spearman Tables for significant diffusion measures ######
#across groups
scd_status_matrix <- scd_status %>% select(age, homeint_storyrecall_d, additional_hads_anxiety, additional_hads_depression)
left_imaging_matrix <- left_amygdala %>% 
  select(dti_fa, fit_FWF, dki_kfa) %>%
  rename_with( ~ paste0(.x, "_left"))
right_imaging_matrix <- right_amygdala %>%
  select(dki_kfa) %>%
  rename_with( ~ paste0(.x, "_right"))
imaging_matrix <- cbind(left_imaging_matrix, right_imaging_matrix)
corr_matrix_pearson <- psych::corr.test(scd_status_matrix, imaging_matrix, use = "pairwise.complete.obs", 
                                adjust = "fdr")
corr_pearson <- as.data.frame(corr_matrix_pearson$r)
corr_pearson <- tibble::rownames_to_column(corr_pearson, "behav_demo_metric")
corr_pearson_long <- corr_pearson %>%
  pivot_longer( cols = !behav_demo_metric,
    names_to = "diffusion_metric", values_to = "Pearson")
corr_pearson_p <- as.data.frame(corr_matrix_pearson$p) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_pearson_p <- tibble::rownames_to_column(corr_pearson_p, "behav_demo_metric")
corr_pearson_p_long <- corr_pearson_p %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Pearson p value")
corr_pearson_long <- left_join(corr_pearson_long, corr_pearson_p_long)
corr_pearson_p_adj <- as.data.frame(corr_matrix_pearson$p.adj) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_pearson_p_adj <- tibble::rownames_to_column(corr_pearson_p_adj, "behav_demo_metric")
corr_pearson_p_adj_long <- corr_pearson_p_adj %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Pearson q value")
corr_pearson_long <- left_join(corr_pearson_long, corr_pearson_p_adj_long)

corr_matrix_spearman <- psych::corr.test(scd_status_matrix, imaging_matrix, use = "pairwise.complete.obs",
                                              method = "spearman", adjust = "fdr")
corr_spearman <- as.data.frame(corr_matrix_spearman$r)
corr_spearman <- tibble::rownames_to_column(corr_spearman, "behav_demo_metric")
corr_spearman_long <- corr_spearman %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Spearman")
corr_spearman_p <- as.data.frame(corr_matrix_spearman$p) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_spearman_p <- tibble::rownames_to_column(corr_spearman_p, "behav_demo_metric")
corr_spearman_p_long <- corr_spearman_p %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Spearman p value")
corr_spearman_long <- left_join(corr_spearman_long, corr_spearman_p_long)
corr_spearman_p_adj <- as.data.frame(corr_matrix_spearman$p.adj) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_spearman_p_adj <- tibble::rownames_to_column(corr_spearman_p_adj, "behav_demo_metric")
corr_spearman_p_adj_long <- corr_spearman_p_adj %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Spearman q value")
corr_spearman_long <- left_join(corr_spearman_long, corr_spearman_p_adj_long)

corr_table <- left_join(corr_pearson_long, corr_spearman_long)
tinytable::tt(corr_table, digits = 2, caption = "Across Groups")

corrplot::corrplot(corr_matrix_spearman$r, p.mat = corr_matrix_spearman$p, method = 'color',
                   sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', pch.cex = 0.9)

#scd_only
scd_status_matrix <- scd_status %>% 
  filter(SCD == "SCD") %>%
  select(age, homeint_storyrecall_d, additional_hads_anxiety, additional_hads_depression)
left_imaging_matrix <- left_amygdala %>% 
  filter(SCD == "SCD") %>%
  select(dti_fa, fit_FWF, dki_kfa) %>%
  rename_with( ~ paste0(.x, "_left"))
right_imaging_matrix <- right_amygdala %>% 
  filter(SCD == "SCD") %>%
  select(dki_kfa) %>%
  rename_with( ~ paste0(.x, "_right"))
imaging_matrix <- cbind(left_imaging_matrix, right_imaging_matrix)
corr_matrix_pearson <- psych::corr.test(scd_status_matrix, imaging_matrix, use = "pairwise.complete.obs", 
                                        adjust = "fdr")
corr_pearson <- as.data.frame(corr_matrix_pearson$r)
corr_pearson <- tibble::rownames_to_column(corr_pearson, "behav_demo_metric")
corr_pearson_long <- corr_pearson %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Pearson")
corr_pearson_p <- as.data.frame(corr_matrix_pearson$p) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_pearson_p <- tibble::rownames_to_column(corr_pearson_p, "behav_demo_metric")
corr_pearson_p_long <- corr_pearson_p %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Pearson p value")
corr_pearson_long <- left_join(corr_pearson_long, corr_pearson_p_long)
corr_pearson_p_adj <- as.data.frame(corr_matrix_pearson$p.adj) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_pearson_p_adj <- tibble::rownames_to_column(corr_pearson_p_adj, "behav_demo_metric")
corr_pearson_p_adj_long <- corr_pearson_p_adj %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Pearson q value")
corr_pearson_long <- left_join(corr_pearson_long, corr_pearson_p_adj_long)

corr_matrix_spearman <- psych::corr.test(scd_status_matrix, imaging_matrix, use = "pairwise.complete.obs",
                                         method = "spearman", adjust = "fdr")
corr_spearman <- as.data.frame(corr_matrix_spearman$r)
corr_spearman <- tibble::rownames_to_column(corr_spearman, "behav_demo_metric")
corr_spearman_long <- corr_spearman %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Spearman")
corr_spearman_p <- as.data.frame(corr_matrix_spearman$p) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_spearman_p <- tibble::rownames_to_column(corr_spearman_p, "behav_demo_metric")
corr_spearman_p_long <- corr_spearman_p %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Spearman p value")
corr_spearman_long <- left_join(corr_spearman_long, corr_spearman_p_long)
corr_spearman_p_adj <- as.data.frame(corr_matrix_spearman$p.adj) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_spearman_p_adj <- tibble::rownames_to_column(corr_spearman_p_adj, "behav_demo_metric")
corr_spearman_p_adj_long <- corr_spearman_p_adj %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Spearman q value")
corr_spearman_long <- left_join(corr_spearman_long, corr_spearman_p_adj_long)

corr_table <- left_join(corr_pearson_long, corr_spearman_long)
tinytable::tt(corr_table, digits = 2, caption = "SCD Only")

corrplot::corrplot(corr_matrix_spearman$r, p.mat = corr_matrix_spearman$p, method = 'color',
                   sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', pch.cex = 0.9)


#ctl_only
scd_status_matrix <- scd_status %>% 
  filter(SCD == "Control") %>%
  select(age, homeint_storyrecall_d, additional_hads_anxiety, additional_hads_depression)
left_imaging_matrix <- left_amygdala %>% 
  filter(SCD == "Control") %>%
  select(dti_fa, fit_FWF, dki_kfa) %>%
  rename_with( ~ paste0(.x, "_left"))
right_imaging_matrix <- right_amygdala %>% 
  filter(SCD == "Control") %>%
  select(dki_kfa) %>%
  rename_with( ~ paste0(.x, "_right"))
imaging_matrix <- cbind(left_imaging_matrix, right_imaging_matrix)
corr_matrix_pearson <- psych::corr.test(scd_status_matrix, imaging_matrix, use = "pairwise.complete.obs", 
                                        adjust = "fdr")
corr_pearson <- as.data.frame(corr_matrix_pearson$r)
corr_pearson <- tibble::rownames_to_column(corr_pearson, "behav_demo_metric")
corr_pearson_long <- corr_pearson %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Pearson")
corr_pearson_p <- as.data.frame(corr_matrix_pearson$p) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_pearson_p <- tibble::rownames_to_column(corr_pearson_p, "behav_demo_metric")
corr_pearson_p_long <- corr_pearson_p %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Pearson p value")
corr_pearson_long <- left_join(corr_pearson_long, corr_pearson_p_long)
corr_pearson_p_adj <- as.data.frame(corr_matrix_pearson$p.adj) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_pearson_p_adj <- tibble::rownames_to_column(corr_pearson_p_adj, "behav_demo_metric")
corr_pearson_p_adj_long <- corr_pearson_p_adj %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Pearson q value")
corr_pearson_long <- left_join(corr_pearson_long, corr_pearson_p_adj_long)

corr_matrix_spearman <- psych::corr.test(Control_status_matrix, imaging_matrix, use = "pairwise.complete.obs",
                                         method = "spearman", adjust = "fdr")
corr_spearman <- as.data.frame(corr_matrix_spearman$r)
corr_spearman <- tibble::rownames_to_column(corr_spearman, "behav_demo_metric")
corr_spearman_long <- corr_spearman %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Spearman")
corr_spearman_p <- as.data.frame(corr_matrix_spearman$p) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_spearman_p <- tibble::rownames_to_column(corr_spearman_p, "behav_demo_metric")
corr_spearman_p_long <- corr_spearman_p %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Spearman p value")
corr_spearman_long <- left_join(corr_spearman_long, corr_spearman_p_long)
corr_spearman_p_adj <- as.data.frame(corr_matrix_spearman$p.adj) %>%
  mutate(across(everything(), ~ format.pval(.x, digits = 2, eps = 0.001)))
corr_spearman_p_adj <- tibble::rownames_to_column(corr_spearman_p_adj, "behav_demo_metric")
corr_spearman_p_adj_long <- corr_spearman_p_adj %>%
  pivot_longer( cols = !behav_demo_metric,
                names_to = "diffusion_metric", values_to = "Spearman q value")
corr_spearman_long <- left_join(corr_spearman_long, corr_spearman_p_adj_long)

corr_table <- left_join(corr_pearson_long, corr_spearman_long)
tinytable::tt(corr_table, digits = 2, caption = "Control Only")

corrplot::corrplot(corr_matrix_spearman$r, p.mat = corr_matrix_spearman$p, method = 'color',
                   sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', pch.cex = 0.9)

