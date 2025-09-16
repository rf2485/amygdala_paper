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
scd_status <- dwi_over_55 %>% select(participant_id, SCD, mt_tr, Income, Ethnicity, 
                                     Sex, age, age_education_completed,homeint_storyrecall_d,
                                     additional_hads_anxiety, additional_hads_depression) %>%
  filter(!participant_id %in% failed_qc) %>%
  rename("Group" = "SCD")
scd_status$additional_hads_depression[scd_status$additional_hads_depression > 21] <- NA

#demographics table
demo_table <- scd_status %>% select(Group, age, age_education_completed, Sex, Income, Ethnicity) %>%
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
  modify_footnote(add_stat_1 ~ "Cohen's D and Cramer's V; Positive effect sizes 
                  indicate larger values in SCD.",
                  stat_1 ~ "Mean (Min-Max); n (%); Mean (Standard Deviation)",
                  stat_2 ~ "Mean (Min-Max); n (%); Mean (Standard Deviation)")
  
#cog and psych table
cog_psych_table <- scd_status %>% select(Group, homeint_storyrecall_d, 
                      additional_hads_anxiety, additional_hads_depression) %>%
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})") %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "",
                estimate ~ "**Effect Size**") %>%
  bold_labels()

tbl_stack(tbls = list(demo_table, cog_psych_table)) %>%
  as_gt() %>%
  gt::gtsave(filename = "demo_cog_psych_table.docx")

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
left_amygdala <- volumes %>% dplyr::select(participant_id:additional_hads_depression, Left.Amygdala) %>%
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


right_amygdala <- volumes %>% dplyr::select(participant_id:additional_hads_depression, Right.Amygdala) %>%
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
              include = c(volume, dti_fa, dti_md, dki_kfa, dki_mk, 
                          fit_FWF, fit_NDI, fit_ODI)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p(t=0.025) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**") %>%
  as_gt() %>%
  gt::gtsave(filename = "left_amygdala_table.docx")

right_amygdala %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(volume, dti_fa, dti_md, dki_kfa, dki_mk, 
                          fit_FWF, fit_NDI, fit_ODI)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% bold_p(t=0.025) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**") %>%
  as_gt() %>%
  gt::gtsave(filename = "right_amygdala_table.docx")

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

left_corr_matrix_pearson <- psych::corr.test(left_dti_dki_matrix, left_vol_noddi_matrix)
colnames(left_corr_matrix_pearson$r) <- c("Volume", "NDI", "FWF", "ODI")
rownames(left_corr_matrix_pearson$r) <- c("FA", "MD", "MK", "KFA")
right_corr_matrix_pearson <- psych::corr.test(right_dti_dki_matrix, right_vol_noddi_matrix)
colnames(right_corr_matrix_pearson$r) <- c("Volume", "NDI", "FWF", "ODI")
rownames(right_corr_matrix_pearson$r) <- c("FA", "MD", "MK", "KFA")

corrplot(left_corr_matrix_pearson$r, p.mat = left_corr_matrix_pearson$p, method = 'color',
         addCoef.col = "black", tl.col = "black", tl.srt = 0, tl.offset = 1,
         title = "(a) Left Amygdala", mar = c(0,0,1,0),
         sig.level = c(0.001, 0.01, 0.03), insig = 'label_sig', pch.cex = 0.9)

corrplot(right_corr_matrix_pearson$r, p.mat = right_corr_matrix_pearson$p, method = 'color',
         addCoef.col = "black", tl.col = "black", tl.srt = 0, tl.offset = 1, 
         title = "(b) Right Amygdala", mar = c(0,0,1,0),
         sig.level = c(0.001, 0.01, 0.03), insig = 'label_sig', pch.cex = 0.9)

##### Linear Models ######
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
  plot_annotation(tag_levels = "a")

ggsave(filename = "diffusion_anxiety_plots_sig.tif", diffusion_anxiety_plots_sig,
       width = 7.5, height = 3.75, dpi = 600, units = "in", device='tiff')

diffusion_anxiety_plots_nonsig <- plot_left_amygdala_dki_kfa + 
  plot_right_amygdala_dki_kfa + 
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave(filename = "diffusion_anxiety_plots_nonsig.tif", diffusion_anxiety_plots_nonsig,
       width = 7.5, height = 3.75, dpi = 600, units = "in", device='tiff')

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

diffusion_memory_plots <- plot_left_amygdala_dti_fa + 
  plot_left_amygdala_fit_FWF + 
  plot_left_amygdala_dki_kfa + 
  plot_right_amygdala_dki_kfa +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave(filename = "diffusion_memory_plots.tif", diffusion_memory_plots,
       width = 7.5, height = 7.5, dpi = 600, units = "in", device='tiff')