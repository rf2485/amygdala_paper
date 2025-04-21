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

trace(corrplot, edit = T) #change line 446 to text(pos.pNew[, 1][sig.locs], (pos.pNew[, 2][sig.locs])+0.25, 
corr_matrix_pearson <- psych::corr.test(dti_dki_matrix, vol_noddi_matrix)
colnames(corr_matrix_pearson$r) <- c("Left Volume", "Left NDI", "Left FWF", "Left ODI",
                                   "Right Volume", "Right NDI", "Right FWF", "Right ODI")

corrplot(corr_matrix_pearson$r, p.mat = corr_matrix_pearson$p, method = 'color',
         addCoef.col = "black",
                   sig.level = c(0.001, 0.01, 0.025), insig = 'label_sig', pch.cex = 0.9)

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
                estimate ~ "**Effect Size**") %>%
  as_gt() %>%
  gt::gtsave(filename = "left_amygdala_vol_noddi_table.docx")

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
                estimate ~ "**Effect Size**") %>%
  as_gt() %>%
  gt::gtsave(filename = "right_amygdala_vol_noddi_table.docx")

##### Linear Models ######
glm_left_amygdala_dti_fa_age_int <- lm(dti_fa ~ age * Group, left_amygdala)
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
       width = 5, height = 3, dpi = 600, units = "in", device = "tiff")

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
