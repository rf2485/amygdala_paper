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
### 3.1	Group differences between SCD and controls ###

#extract SCD status and demographics for participants
scd_status <- dwi_over_55 %>% select(participant_id, SCD, mt_tr, Income, 
                                     Ethnicity, Sex, age, age_education_completed,
                                     homeint_storyrecall_d, bp_dia_mean_cardio, 
                                     bp_sys_mean_cardio, pulse_mean_cardio,
                                     height_cardio, weight_cardio,
                                     additional_hads_anxiety, additional_hads_depression) %>%
  filter(!participant_id %in% failed_qc) %>%
  rename("Group" = "SCD") %>%
  mutate(bmi_cardio = weight_cardio / (height_cardio/100)^2) %>%
  select(!height_cardio, !weight_cardio)
scd_status$additional_hads_depression[scd_status$additional_hads_depression > 21] <- NA

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

#cog and psych table
scd_status %>% select(Group, homeint_storyrecall_d, additional_hads_anxiety, additional_hads_depression) %>%
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})") %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "",
                estimate ~ "**Effect Size**") %>%
  bold_labels() %>%
  as_gt() %>%
  gt::gtsave(filename = "cog_and_psych_table.docx")

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
                estimate ~ "**Effect Size**") %>%
  as_gt() %>%
  gt::gtsave(filename = "left_amygdala_dti_dki_table.docx")

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
                estimate ~ "**Effect Size**") %>%
  as_gt() %>%
  gt::gtsave(filename = "right_amygdala_dti_dki_table.docx")

#### Within SCD correlations between DTI, DKI, volume, and NODDI ####
left_dti_dki_matrix <- left_amygdala %>% 
  filter(Group == "SCD") %>%
  select(dti_fa, dki_kfa) %>%
  rename_with( ~ paste0(.x, "_left"))
right_dti_dki_matrix <- right_amygdala %>% 
  filter(Group == "SCD") %>%
  select(dki_kfa) %>%
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

trace(corrplot, edit = T) #change line 447 to text(pos.pNew[, 1][sig.locs], (pos.pNew[, 2][sig.locs])+0.25, 
corr_matrix_pearson <- psych::corr.test(dti_dki_matrix, vol_noddi_matrix)
colnames(corr_matrix_pearson$r) <- c("Left Volume", "Left NDI", "Left FWF", "Left ODI",
                                   "Right Volume", "Right NDI", "Right FWF", "Right ODI")

corrplot(corr_matrix_pearson$r, p.mat = corr_matrix_pearson$p, method = 'color',
         addCoef.col = "black",
                   sig.level = c(0.001, 0.01, 0.025), insig = 'label_sig', pch.cex = 0.9)


##### Scatterplots of Linear Models ######
left_amygdala_ctl <- left_amygdala %>% filter(Group == "Control")
left_amygdala_scd <- left_amygdala %>% filter(Group == "SCD")
right_amygdala_ctl <- right_amygdala %>% filter(Group == "Control")
right_amygdala_scd <- right_amygdala %>% filter(Group == "SCD")

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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
  plot_annotation(tag_levels = "a")
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
  plot_annotation(tag_levels = "a")
ggsave(filename = "left_amygdala_dki_kfa_bio_plots.tif", left_amygdala_dki_kfa_bio_plots,
       width = 6.5, height = 5.5, dpi = 600, units = "in", device='tiff')

glm_left_amygdala_fit_FWF_volume_scd <- lm(fit_FWF ~ volume, left_amygdala_scd)
summary(glm_left_amygdala_fit_FWF_volume_scd)
corr_p_value <- summary(glm_left_amygdala_fit_FWF_volume_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_fit_FWF_volume_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_fit_FWF_volume_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_fit_FWF_volume_scd, c( "volume[all]"))
plot_glm_left_amygdala_fit_FWF_volume_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    x = expression(Volume~(mm^3)),
    title = "Left Amygdala Volume",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_left_amygdala_fit_FWF_fit_NDI_scd <- lm(fit_FWF ~ fit_NDI, left_amygdala_scd)
summary(glm_left_amygdala_fit_FWF_fit_NDI_scd)
corr_p_value <- summary(glm_left_amygdala_fit_FWF_fit_NDI_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_fit_FWF_fit_NDI_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_fit_FWF_fit_NDI_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_fit_FWF_fit_NDI_scd, c( "fit_NDI[all]"))
plot_glm_left_amygdala_fit_FWF_fit_NDI_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Left Amygdala Mean NDI",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())
glm_left_amygdala_fit_FWF_fit_ODI_scd <- lm(fit_FWF ~ fit_ODI, left_amygdala_scd)
summary(glm_left_amygdala_fit_FWF_fit_ODI_scd)
corr_p_value <- summary(glm_left_amygdala_fit_FWF_fit_ODI_scd)$coefficients[2,4]
corr_beta <- summary(glm_left_amygdala_fit_FWF_fit_ODI_scd)$coefficients[2,1]
adj_r_squared <- summary(glm_left_amygdala_fit_FWF_fit_ODI_scd)$adj.r.squared
pr <- predict_response(glm_left_amygdala_fit_FWF_fit_ODI_scd, c( "fit_ODI[all]"))
plot_glm_left_amygdala_fit_FWF_fit_ODI_scd <- 
  plot(pr, show_data = T, dot_alpha = 1, dot_size = 0.5) + 
  labs(
    title = "Left Amygdala Mean ODI",
    subtitle = paste0(
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
      " p: ", format.pval(corr_p_value, digits = 2, eps = 0.001),
      ", \u03B2: ", signif(corr_beta, 2),
      ", adj-R<sup>2</sup>: ", signif(adj_r_squared, 2))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        legend.title=element_blank())

left_amygdala_fit_FWF_bio_plots <- plot_glm_left_amygdala_fit_FWF_volume_scd + 
  plot_glm_left_amygdala_fit_FWF_fit_NDI_scd + 
  plot_glm_left_amygdala_fit_FWF_fit_ODI_scd + 
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a")
ggsave(filename = "left_amygdala_fit_FWF_bio_plots.tif", left_amygdala_fit_FWF_bio_plots,
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
      symnum(corr_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
  plot_annotation(tag_levels = "a")
ggsave(filename = "right_amygdala_dki_kfa_bio_plots.tif", right_amygdala_dki_kfa_bio_plots,
       width = 6.5, height = 5.5, dpi = 600, units = "in", device='tiff')

###age interaction
glm_left_amygdala_dti_fa_age_int <- lm(dti_fa ~ age * Group, left_amygdala)
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
plot_left_amygdala_dti_fa <- interactions::interact_plot(glm_left_amygdala_dti_fa_age_int, pred = age, modx = Group, point.alpha = 1,
                                                         plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "FA", x = "Age (Years)", title = "Left Amygdala Mean FA", fill = "Group",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                       " p: ", format.pval(ctl_p_value, digits = 2, eps = 0.001),
                       ", \u03B2: ", signif(ctl_beta, 2),
                       ", adj-R<sup>2</sup>: ", signif(ctl_adj_r_squared, 2)
                     ),
                     color = "Control"), 
                show.legend = F, size = 3, fill = NA, label.color = NA, 
                label.padding = grid::unit(rep(0,4), "pt")
  ) +
  ylim(NA,0.27) +
  theme_bw(base_size = 10, ) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

glm_left_amygdala_fit_FWF_age_int <- lm(fit_FWF ~ age * Group, left_amygdala)
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
plot_left_amygdala_fit_FWF <- interactions::interact_plot(glm_left_amygdala_fit_FWF_age_int, pred = age, modx = Group, point.alpha = 1,
                                                          plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "FWF", x = "Age (Years)", title = "Left Amygdala Mean FWF",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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

glm_left_amygdala_dki_kfa_age_int <- lm(dki_kfa ~ age * Group, left_amygdala)
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
plot_left_amygdala_dki_kfa <- interactions::interact_plot(glm_left_amygdala_dki_kfa_age_int, pred = age, modx = Group, point.alpha = 1,
                                                          plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "KFA", x = "Age (Years)", title = "Left Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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

glm_right_amygdala_dki_kfa_age_int <- lm(dki_kfa ~ age * Group, right_amygdala)
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
plot_right_amygdala_dki_kfa <- interactions::interact_plot(glm_right_amygdala_dki_kfa_age_int, pred = age, modx = Group, point.alpha = 1,
                                                           plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "KFA", x = "Age (Years)", title = "Right Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
  plot_annotation(tag_levels = "a")

ggsave(filename = "diffusion_age_plots.tif", diffusion_age_plots,
       width = 7.5, height = 7.5, dpi = 600, units = "in", device='tiff')

###anxiety interaction
clm_right_amygdala_dki_kfa_depression_int <- ordinal::clm(ordered(additional_hads_depression) ~ dki_kfa * Group, data = right_amygdala)
summary(clm_right_amygdala_dki_kfa_depression_int)
across_group_p_value <- summary(clm_right_amygdala_dki_kfa_depression_int)$coefficients["dki_kfa", "Pr(>|z|)"]
across_group_beta <- summary(clm_right_amygdala_dki_kfa_depression_int)$coefficients["dki_kfa", "Estimate"]
across_group_odds_ratio <- exp(across_group_beta)
int_p_value <- summary(clm_right_amygdala_dki_kfa_depression_int)$coefficients["dki_kfa:GroupControl", "Pr(>|z|)"]
int_beta <- summary(clm_right_amygdala_dki_kfa_depression_int)$coefficients["dki_kfa:GroupControl", "Estimate"]
int_odds_ratio <- exp(int_beta)
clm_right_amygdala_dki_kfa_depression_scd <- ordinal::clm(ordered(additional_hads_depression) ~ dki_kfa, data = right_amygdala_scd)
summary(clm_right_amygdala_dki_kfa_depression_scd)
scd_p_value <- summary(clm_right_amygdala_dki_kfa_depression_scd)$coefficients["dki_kfa", "Pr(>|z|)"]
scd_beta <- summary(clm_right_amygdala_dki_kfa_depression_scd)$coefficients["dki_kfa", "Estimate"]
scd_odds_ratio <- exp(scd_beta)
clm_right_amygdala_dki_kfa_depression_ctl <- ordinal::clm(ordered(additional_hads_depression) ~ dki_kfa, data = right_amygdala_ctl)
summary(clm_right_amygdala_dki_kfa_depression_ctl)
ctl_p_value <- summary(clm_right_amygdala_dki_kfa_depression_ctl)$coefficients["dki_kfa", "Pr(>|z|)"]
ctl_beta <- summary(clm_right_amygdala_dki_kfa_depression_ctl)$coefficients["dki_kfa", "Estimate"]
ctl_odds_ratio <- exp(ctl_beta)
marginaleffects::plot_predictions(clm_right_amygdala_dki_kfa_depression_int, type = "prob", condition = "dki_kfa")
predict_response(clm_right_amygdala_dki_kfa_depression_int, c("dki_kfa", "Group")) %>% plot()

ggplot(right_amygdala, aes(dti_fa, additional_hads_anxiety, color=Group)) +
  geom_point()

glm_left_amygdala_dti_fa_anxiety_int <- lm(dti_fa ~ additional_hads_anxiety * Group, left_amygdala)
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
plot_left_amygdala_dti_fa <- interactions::interact_plot(glm_left_amygdala_dti_fa_anxiety_int, pred = additional_hads_anxiety, modx = Group, point.alpha = 1,
                                                         plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "FA", x = "Anxiety (HADS)", title = "Left Amygdala Mean FA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
glm_left_amygdala_fit_FWF_anxiety_int <- lm(fit_FWF ~ additional_hads_anxiety * Group, left_amygdala)
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
plot_left_amygdala_fit_FWF <- interactions::interact_plot(glm_left_amygdala_fit_FWF_anxiety_int, pred = additional_hads_anxiety, modx = Group, point.alpha = 1,
                                                          plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "FWF", x = "Anxiety (HADS)", title = "Left Amygdala Mean FWF",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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

glm_left_amygdala_dki_kfa_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * Group, left_amygdala)
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
plot_left_amygdala_dki_kfa <- interactions::interact_plot(glm_left_amygdala_dki_kfa_anxiety_int, pred = additional_hads_anxiety, modx = Group, point.alpha = 1,
                                                          plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "KFA", x = "Anxiety (HADS)", title = "Left Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
glm_right_amygdala_dki_kfa_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * Group, right_amygdala)
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
plot_right_amygdala_dki_kfa <- interactions::interact_plot(glm_right_amygdala_dki_kfa_anxiety_int, pred = additional_hads_anxiety, modx = Group, point.alpha = 1,
                                                           plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "KFA", x = "Anxiety (HADS)", title = "Right Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
  plot_annotation(tag_levels = "a")

ggsave(filename = "diffusion_anxiety_plots.tif", diffusion_anxiety_plots,
       width = 7.5, height = 7.5, dpi = 600, units = "in", device='tiff')

###depression interaction
glm_left_amygdala_dti_fa_depression_int <- lm(dti_fa ~ additional_hads_depression * Group, left_amygdala)
summary(glm_left_amygdala_dti_fa_depression_int)
across_group_p_value <- summary(glm_left_amygdala_dti_fa_depression_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_dti_fa_depression_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_dti_fa_depression_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_dti_fa_depression_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_dti_fa_depression_int)$coefficients[4,1]                
glm_left_amygdala_dti_fa_depression_scd <- lm(dti_fa ~ additional_hads_depression, left_amygdala_scd)
summary(glm_left_amygdala_dti_fa_depression_scd)
scd_p_value <- summary(glm_left_amygdala_dti_fa_depression_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_dti_fa_depression_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_dti_fa_depression_scd)$adj.r.squared
glm_left_amygdala_dti_fa_depression_ctl <- lm(dti_fa ~ additional_hads_depression, left_amygdala_ctl)
summary(glm_left_amygdala_dti_fa_depression_ctl)
ctl_p_value <- summary(glm_left_amygdala_dti_fa_depression_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_dti_fa_depression_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_dti_fa_depression_ctl)$adj.r.squared
plot_left_amygdala_dti_fa <- interactions::interact_plot(glm_left_amygdala_dti_fa_depression_int, pred = additional_hads_depression, modx = Group, point.alpha = 1,
                                                         plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "FA", x = "Depression (HADS)", title = "Left Amygdala Mean FA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
glm_left_amygdala_fit_FWF_depression_int <- lm(fit_FWF ~ additional_hads_depression * Group, left_amygdala)
summary(glm_left_amygdala_fit_FWF_depression_int)
across_group_p_value <- summary(glm_left_amygdala_fit_FWF_depression_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_fit_FWF_depression_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_depression_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_fit_FWF_depression_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_fit_FWF_depression_int)$coefficients[4,1]                
glm_left_amygdala_fit_FWF_depression_scd <- lm(fit_FWF ~ additional_hads_depression, left_amygdala_scd)
summary(glm_left_amygdala_fit_FWF_depression_scd)
scd_p_value <- summary(glm_left_amygdala_fit_FWF_depression_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_fit_FWF_depression_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_depression_scd)$adj.r.squared
glm_left_amygdala_fit_FWF_depression_ctl <- lm(fit_FWF ~ additional_hads_depression, left_amygdala_ctl)
summary(glm_left_amygdala_fit_FWF_depression_ctl)
ctl_p_value <- summary(glm_left_amygdala_fit_FWF_depression_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_fit_FWF_depression_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_depression_ctl)$adj.r.squared
plot_left_amygdala_fit_FWF <- interactions::interact_plot(glm_left_amygdala_fit_FWF_depression_int, pred = additional_hads_depression, modx = Group, point.alpha = 1,
                                                          plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "FWF", x = "Depression (HADS)", title = "Left Amygdala Mean FWF",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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

glm_left_amygdala_dki_kfa_depression_int <- lm(dki_kfa ~ additional_hads_depression * Group, left_amygdala)
summary(glm_left_amygdala_dki_kfa_depression_int)
across_group_p_value <- summary(glm_left_amygdala_dki_kfa_depression_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_dki_kfa_depression_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_depression_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_dki_kfa_depression_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_dki_kfa_depression_int)$coefficients[4,1]                
glm_left_amygdala_dki_kfa_depression_scd <- lm(dki_kfa ~ additional_hads_depression, left_amygdala_scd)
summary(glm_left_amygdala_dki_kfa_depression_scd)
scd_p_value <- summary(glm_left_amygdala_dki_kfa_depression_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_dki_kfa_depression_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_depression_scd)$adj.r.squared
glm_left_amygdala_dki_kfa_depression_ctl <- lm(dki_kfa ~ additional_hads_depression, left_amygdala_ctl)
summary(glm_left_amygdala_dki_kfa_depression_ctl)
ctl_p_value <- summary(glm_left_amygdala_dki_kfa_depression_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_dki_kfa_depression_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_depression_ctl)$adj.r.squared
plot_left_amygdala_dki_kfa <- interactions::interact_plot(glm_left_amygdala_dki_kfa_depression_int, pred = additional_hads_depression, modx = Group, point.alpha = 1,
                                                          plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "KFA", x = "Depression (HADS)", title = "Left Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
glm_right_amygdala_dki_kfa_depression_int <- lm(dki_kfa ~ additional_hads_depression * Group, right_amygdala)
summary(glm_right_amygdala_dki_kfa_depression_int)
across_group_p_value <- summary(glm_right_amygdala_dki_kfa_depression_int)$coefficients[2,4]
across_group_beta <- summary(glm_right_amygdala_dki_kfa_depression_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_depression_int)$adj.r.squared
int_p_value <- summary(glm_right_amygdala_dki_kfa_depression_int)$coefficients[4,4]
int_beta <-  summary(glm_right_amygdala_dki_kfa_depression_int)$coefficients[4,1]                
glm_right_amygdala_dki_kfa_depression_scd <- lm(dki_kfa ~ additional_hads_depression, right_amygdala_scd)
summary(glm_right_amygdala_dki_kfa_depression_scd)
scd_p_value <- summary(glm_right_amygdala_dki_kfa_depression_scd)$coefficients[2,4]
scd_beta <- summary(glm_right_amygdala_dki_kfa_depression_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_depression_scd)$adj.r.squared
glm_right_amygdala_dki_kfa_depression_ctl <- lm(dki_kfa ~ additional_hads_depression, right_amygdala_ctl)
summary(glm_right_amygdala_dki_kfa_depression_ctl)
ctl_p_value <- summary(glm_right_amygdala_dki_kfa_depression_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_right_amygdala_dki_kfa_depression_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_depression_ctl)$adj.r.squared
plot_right_amygdala_dki_kfa <- interactions::interact_plot(glm_right_amygdala_dki_kfa_depression_int, pred = additional_hads_depression, modx = Group, point.alpha = 1,
                                                           plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "KFA", x = "Depression (HADS)", title = "Right Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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

diffusion_depression_plots <- plot_left_amygdala_dti_fa + 
  plot_left_amygdala_fit_FWF + 
  plot_left_amygdala_dki_kfa + 
  plot_right_amygdala_dki_kfa +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave(filename = "diffusion_depression_plots.tif", diffusion_depression_plots,
       width = 7.5, height = 7.5, dpi = 600, units = "in", device='tiff')

#memory
glm_left_amygdala_dti_fa_memory_int <- lm(dti_fa ~ homeint_storyrecall_d * Group, left_amygdala)
summary(glm_left_amygdala_dti_fa_memory_int)
across_group_p_value <- summary(glm_left_amygdala_dti_fa_memory_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_dti_fa_memory_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_dti_fa_memory_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_dti_fa_memory_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_dti_fa_memory_int)$coefficients[4,1]                
glm_left_amygdala_dti_fa_memory_scd <- lm(dti_fa ~ homeint_storyrecall_d, left_amygdala_scd)
summary(glm_left_amygdala_dti_fa_memory_scd)
scd_p_value <- summary(glm_left_amygdala_dti_fa_memory_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_dti_fa_memory_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_dti_fa_memory_scd)$adj.r.squared
glm_left_amygdala_dti_fa_memory_ctl <- lm(dti_fa ~ homeint_storyrecall_d, left_amygdala_ctl)
summary(glm_left_amygdala_dti_fa_memory_ctl)
ctl_p_value <- summary(glm_left_amygdala_dti_fa_memory_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_dti_fa_memory_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_dti_fa_memory_ctl)$adj.r.squared
plot_left_amygdala_dti_fa <- interactions::interact_plot(glm_left_amygdala_dti_fa_memory_int, pred = homeint_storyrecall_d, modx = Group, point.alpha = 1,
                                                         plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "FA", x = "Memory (# Details)", title = "Left Amygdala Mean FA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
glm_left_amygdala_fit_FWF_memory_int <- lm(fit_FWF ~ homeint_storyrecall_d * Group, left_amygdala)
summary(glm_left_amygdala_fit_FWF_memory_int)
across_group_p_value <- summary(glm_left_amygdala_fit_FWF_memory_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_fit_FWF_memory_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_memory_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_fit_FWF_memory_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_fit_FWF_memory_int)$coefficients[4,1]                
glm_left_amygdala_fit_FWF_memory_scd <- lm(fit_FWF ~ homeint_storyrecall_d, left_amygdala_scd)
summary(glm_left_amygdala_fit_FWF_memory_scd)
scd_p_value <- summary(glm_left_amygdala_fit_FWF_memory_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_fit_FWF_memory_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_memory_scd)$adj.r.squared
glm_left_amygdala_fit_FWF_memory_ctl <- lm(fit_FWF ~ homeint_storyrecall_d, left_amygdala_ctl)
summary(glm_left_amygdala_fit_FWF_memory_ctl)
ctl_p_value <- summary(glm_left_amygdala_fit_FWF_memory_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_fit_FWF_memory_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_fit_FWF_memory_ctl)$adj.r.squared
plot_left_amygdala_fit_FWF <- interactions::interact_plot(glm_left_amygdala_fit_FWF_memory_int, pred = homeint_storyrecall_d, modx = Group, point.alpha = 1,
                                                          plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "FWF", x = "Memory (# Details)", title = "Left Amygdala Mean FWF",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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

glm_left_amygdala_dki_kfa_memory_int <- lm(dki_kfa ~ homeint_storyrecall_d * Group, left_amygdala)
summary(glm_left_amygdala_dki_kfa_memory_int)
across_group_p_value <- summary(glm_left_amygdala_dki_kfa_memory_int)$coefficients[2,4]
across_group_beta <- summary(glm_left_amygdala_dki_kfa_memory_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_memory_int)$adj.r.squared
int_p_value <- summary(glm_left_amygdala_dki_kfa_memory_int)$coefficients[4,4]
int_beta <-  summary(glm_left_amygdala_dki_kfa_memory_int)$coefficients[4,1]                
glm_left_amygdala_dki_kfa_memory_scd <- lm(dki_kfa ~ homeint_storyrecall_d, left_amygdala_scd)
summary(glm_left_amygdala_dki_kfa_memory_scd)
scd_p_value <- summary(glm_left_amygdala_dki_kfa_memory_scd)$coefficients[2,4]
scd_beta <- summary(glm_left_amygdala_dki_kfa_memory_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_memory_scd)$adj.r.squared
glm_left_amygdala_dki_kfa_memory_ctl <- lm(dki_kfa ~ homeint_storyrecall_d, left_amygdala_ctl)
summary(glm_left_amygdala_dki_kfa_memory_ctl)
ctl_p_value <- summary(glm_left_amygdala_dki_kfa_memory_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_left_amygdala_dki_kfa_memory_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_left_amygdala_dki_kfa_memory_ctl)$adj.r.squared
plot_left_amygdala_dki_kfa <- interactions::interact_plot(glm_left_amygdala_dki_kfa_memory_int, pred = homeint_storyrecall_d, modx = Group, point.alpha = 1,
                                                          plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "KFA", x = "Memory (# Details)", title = "Left Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
glm_right_amygdala_dki_kfa_memory_int <- lm(dki_kfa ~ homeint_storyrecall_d * Group, right_amygdala)
summary(glm_right_amygdala_dki_kfa_memory_int)
across_group_p_value <- summary(glm_right_amygdala_dki_kfa_memory_int)$coefficients[2,4]
across_group_beta <- summary(glm_right_amygdala_dki_kfa_memory_int)$coefficients[2,1]
across_group_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_memory_int)$adj.r.squared
int_p_value <- summary(glm_right_amygdala_dki_kfa_memory_int)$coefficients[4,4]
int_beta <-  summary(glm_right_amygdala_dki_kfa_memory_int)$coefficients[4,1]                
glm_right_amygdala_dki_kfa_memory_scd <- lm(dki_kfa ~ homeint_storyrecall_d, right_amygdala_scd)
summary(glm_right_amygdala_dki_kfa_memory_scd)
scd_p_value <- summary(glm_right_amygdala_dki_kfa_memory_scd)$coefficients[2,4]
scd_beta <- summary(glm_right_amygdala_dki_kfa_memory_scd)$coefficients[2,1]
scd_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_memory_scd)$adj.r.squared
glm_right_amygdala_dki_kfa_memory_ctl <- lm(dki_kfa ~ homeint_storyrecall_d, right_amygdala_ctl)
summary(glm_right_amygdala_dki_kfa_memory_ctl)
ctl_p_value <- summary(glm_right_amygdala_dki_kfa_memory_ctl)$coefficients[2,4]
ctl_beta <- summary(glm_right_amygdala_dki_kfa_memory_ctl)$coefficients[2,1]
ctl_adj_r_squared <- summary(glm_right_amygdala_dki_kfa_memory_ctl)$adj.r.squared
plot_right_amygdala_dki_kfa <- interactions::interact_plot(glm_right_amygdala_dki_kfa_memory_int, pred = homeint_storyrecall_d, modx = Group, point.alpha = 1,
                                                           plot.points = T, interval = T, point.size = 0.5, point.shape = T, colors = c("black", "gray50")) +
  labs(
    y = "KFA", x = "Memory (# Details)", title = "Right Amygdala Mean KFA",
    subtitle = paste0("across group ", 
                      symnum(across_group_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(across_group_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(across_group_beta, 2),
                      "<br> adj-R<sup>2</sup>: ", signif(across_group_adj_r_squared, 2),
                      "<br> interaction ",
                      symnum(int_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
                      " p: ", format.pval(int_p_value, digits = 2, eps = 0.001),
                      ", \u03B2: ", signif(int_beta, 2)
    )
  ) +
  geom_richtext(aes_(x = Inf, y = Inf, vjust = 1.1, hjust = 1.01,
                     label = paste0(
                       symnum(scd_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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
                       symnum(ctl_p_value, corr = F, cutpoints = c(0, .001, .01, .025, 1), symbols = c("***", "**", "\\*", "")),
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

diffusion_memory_plots <- plot_left_amygdala_dti_fa + 
  plot_left_amygdala_fit_FWF + 
  plot_left_amygdala_dki_kfa + 
  plot_right_amygdala_dki_kfa +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave(filename = "diffusion_memory_plots.tif", diffusion_memory_plots,
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

#shapiro-wilk test for normality. if violated, can argue data is sufficiently large
#so violation is trivial or transform dependent variable (ie box-cox transform or log transform)
#homogeneity assumption plot residual and y-hat. can use square root or log transform if violated
#levene's test, aka equal variance
#proc transreg procedure: explore different transformation methods

shapiro.test(log(left_amygdala$dti_fa))
shapiro.test(log(left_amygdala$fit_FWF))
shapiro.test(log(left_amygdala$dki_kfa))
shapiro.test(1/right_amygdala$dki_kfa)

new_model <- lm(sqrt(additional_hads_anxiety) ~ dti_fa * Group, data = left_amygdala)
qqnorm(glm_right_amygdala_dki_kfa_anxiety_int$residuals)
qqline(glm_right_amygdala_dki_kfa_anxiety_int$residuals)
qqnorm(new_model$residuals)
qqline(new_model$residuals)

library(car)
leveneTest(dti_fa ~ Group, left_amygdala)
leveneTest(fit_FWF ~ Group, left_amygdala)
leveneTest(dki_kfa ~ Group, left_amygdala)
leveneTest(dki_kfa ~ Group, right_amygdala)
leveneTest(sqrt(additional_hads_anxiety) ~ Group, left_amygdala)
