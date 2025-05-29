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
left_putamen <- dti_ad %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(dti_ad = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- dti_rd %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(dti_rd = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- dki_kfa %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(dki_kfa = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- dki_mk %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(dki_mk = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- dki_ak %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(dki_ak = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen <- dki_rk %>% dplyr::select(participant_id, Left.Putamen) %>%
  rename(dki_rk = "Left.Putamen") %>% full_join(left_putamen, .)
# left_putamen <- mtr_tr30 %>% dplyr::select(participant_id, Left.Putamen) %>%
#   rename(mtr_tr30 = "Left.Putamen") %>% full_join(left_putamen, .)
# left_putamen <- mtr_tr50 %>% dplyr::select(participant_id, Left.Putamen) %>%
#   rename(mtr_tr50 = "Left.Putamen") %>% full_join(left_putamen, .)
left_putamen$region <- "Left putamen"

left_putamen <- set_label(left_putamen,
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

right_putamen <- volumes %>% dplyr::select(participant_id:bmi_cardio, Right.Putamen) %>%
  rename(volume = "Right.Putamen")
right_putamen <- fit_NDI %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(fit_NDI = "Right.Putamen") %>% full_join(right_putamen, .)
right_putamen <- fit_FWF %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(fit_FWF = "Right.Putamen") %>% full_join(right_putamen, .)
right_putamen <- fit_ODI %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(fit_ODI = "Right.Putamen") %>% full_join(right_putamen, .)
right_putamen <- dti_fa %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(dti_fa = "Right.Putamen") %>% full_join(right_putamen, .)
right_putamen <- dti_md %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(dti_md = "Right.Putamen") %>% full_join(right_putamen, .)
right_putamen <- dti_ad %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(dti_ad = "Right.Putamen") %>% full_join(right_putamen, .)
right_putamen <- dti_rd %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(dti_rd = "Right.Putamen") %>% full_join(right_putamen, .)
right_putamen <- dki_kfa %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(dki_kfa = "Right.Putamen") %>% full_join(right_putamen, .)
right_putamen <- dki_mk %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(dki_mk = "Right.Putamen") %>% full_join(right_putamen, .)
right_putamen <- dki_ak %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(dki_ak = "Right.Putamen") %>% full_join(right_putamen, .)
right_putamen <- dki_rk %>% dplyr::select(participant_id, Right.Putamen) %>%
  rename(dki_rk = "Right.Putamen") %>% full_join(right_putamen, .)
# right_putamen <- mtr_tr30 %>% dplyr::select(participant_id, Right.Putamen) %>%
#   rename(mtr_tr30 = "Right.Putamen") %>% full_join(right_putamen, .)
# right_putamen <- mtr_tr50 %>% dplyr::select(participant_id, Right.Putamen) %>%
#   rename(mtr_tr50 = "Right.Putamen") %>% full_join(right_putamen, .)

right_putamen$region <- "Right putamen"
right_putamen <- set_label(right_putamen,
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

left_putamen_ctl <- left_putamen %>% filter(Group == "Control")
left_putamen_scd <- left_putamen %>% filter(Group == "SCD")
right_putamen_ctl <- right_putamen %>% filter(Group == "Control")
right_putamen_scd <- right_putamen %>% filter(Group == "SCD")

left_putamen %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(volume, dti_fa, dti_md, dki_kfa, dki_mk, 
                          fit_FWF, fit_NDI, fit_ODI)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p() %>% bold_p(q=T) %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**")

right_putamen %>% 
  select(-participant_id) %>% 
  tbl_summary(by = Group, statistic = all_continuous() ~ "{mean} ({sd})",
              include = c(dti_fa, dti_md, dki_kfa, dki_mk, 
                          volume, fit_FWF, fit_NDI, fit_ODI)
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p() %>% bold_p(q=T) %>% #filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Imaging Metric**",
                estimate ~ "**Effect Size**") 

###age
glm_left_putamen_volume_age_int <- lm(volume ~ age * Group, left_putamen)
summary(glm_left_putamen_volume_age_int) #associated with age

glm_left_putamen_dti_fa_age_int <- lm(dti_fa ~ age * Group, left_putamen)
summary(glm_left_putamen_dti_fa_age_int)

glm_left_putamen_dti_md_age_int <- lm(dti_md ~ age * Group, left_putamen)
summary(glm_left_putamen_dti_md_age_int) #associated with age

glm_left_putamen_dki_kfa_age_int <- lm(dki_kfa ~ age * Group, left_putamen)
summary(glm_left_putamen_dki_kfa_age_int) #associated with age

glm_left_putamen_dki_mk_age_int <- lm(dki_mk ~ age * Group, left_putamen)
summary(glm_left_putamen_dki_mk_age_int)

glm_left_putamen_fit_FWF_age_int <- lm(fit_FWF ~ age * Group, left_putamen)
summary(glm_left_putamen_fit_FWF_age_int) #associated with age

glm_left_putamen_fit_NDI_age_int <- lm(fit_NDI ~ age * Group, left_putamen)
summary(glm_left_putamen_fit_NDI_age_int)

glm_left_putamen_fit_ODI_age_int <- lm(fit_ODI ~ age * Group, left_putamen)
summary(glm_left_putamen_fit_ODI_age_int) #associated with age

###homeint_storyrecall_d
glm_left_putamen_volume_homeint_storyrecall_d_int <- lm(volume ~ homeint_storyrecall_d * Group, left_putamen)
summary(glm_left_putamen_volume_homeint_storyrecall_d_int) #associated with memory

glm_left_putamen_dti_fa_homeint_storyrecall_d_int <- lm(dti_fa ~ homeint_storyrecall_d * Group, left_putamen)
summary(glm_left_putamen_dti_fa_homeint_storyrecall_d_int)

glm_left_putamen_dti_md_homeint_storyrecall_d_int <- lm(dti_md ~ homeint_storyrecall_d * Group, left_putamen)
summary(glm_left_putamen_dti_md_homeint_storyrecall_d_int) #associated with memory, group diff

glm_left_putamen_dki_kfa_homeint_storyrecall_d_int <- lm(dki_kfa ~ homeint_storyrecall_d * Group, left_putamen)
summary(glm_left_putamen_dki_kfa_homeint_storyrecall_d_int) #group diff

glm_left_putamen_dki_mk_homeint_storyrecall_d_int <- lm(dki_mk ~ homeint_storyrecall_d * Group, left_putamen)
summary(glm_left_putamen_dki_mk_homeint_storyrecall_d_int)

glm_left_putamen_fit_FWF_homeint_storyrecall_d_int <- lm(fit_FWF ~ homeint_storyrecall_d * Group, left_putamen)
summary(glm_left_putamen_fit_FWF_homeint_storyrecall_d_int) #group diff

glm_left_putamen_fit_NDI_homeint_storyrecall_d_int <- lm(fit_NDI ~ homeint_storyrecall_d * Group, left_putamen)
summary(glm_left_putamen_fit_NDI_homeint_storyrecall_d_int)

glm_left_putamen_fit_ODI_homeint_storyrecall_d_int <- lm(fit_ODI ~ homeint_storyrecall_d * Group, left_putamen)
summary(glm_left_putamen_fit_ODI_homeint_storyrecall_d_int) #associated with memory

###additional_hads_anxiety
glm_left_putamen_volume_additional_hads_anxiety_int <- lm(volume ~ additional_hads_anxiety * Group, left_putamen)
summary(glm_left_putamen_volume_additional_hads_anxiety_int) #associated with anxiety

glm_left_putamen_dti_fa_additional_hads_anxiety_int <- lm(dti_fa ~ additional_hads_anxiety * Group, left_putamen)
summary(glm_left_putamen_dti_fa_additional_hads_anxiety_int)

glm_left_putamen_dti_md_additional_hads_anxiety_int <- lm(dti_md ~ additional_hads_anxiety * Group, left_putamen)
summary(glm_left_putamen_dti_md_additional_hads_anxiety_int) #group diff

glm_left_putamen_dki_kfa_additional_hads_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * Group, left_putamen)
summary(glm_left_putamen_dki_kfa_additional_hads_anxiety_int)

glm_left_putamen_dki_mk_additional_hads_anxiety_int <- lm(dki_mk ~ additional_hads_anxiety * Group, left_putamen)
summary(glm_left_putamen_dki_mk_additional_hads_anxiety_int)

glm_left_putamen_fit_FWF_additional_hads_anxiety_int <- lm(fit_FWF ~ additional_hads_anxiety * Group, left_putamen)
summary(glm_left_putamen_fit_FWF_additional_hads_anxiety_int) #group diff

glm_left_putamen_fit_NDI_additional_hads_anxiety_int <- lm(fit_NDI ~ additional_hads_anxiety * Group, left_putamen)
summary(glm_left_putamen_fit_NDI_additional_hads_anxiety_int)

glm_left_putamen_fit_ODI_additional_hads_anxiety_int <- lm(fit_ODI ~ additional_hads_anxiety * Group, left_putamen)
summary(glm_left_putamen_fit_ODI_additional_hads_anxiety_int)

###additional_hads_depression
glm_left_putamen_volume_additional_hads_depression_int <- lm(volume ~ additional_hads_depression * Group, left_putamen)
summary(glm_left_putamen_volume_additional_hads_depression_int)

glm_left_putamen_dti_fa_additional_hads_depression_int <- lm(dti_fa ~ additional_hads_depression * Group, left_putamen)
summary(glm_left_putamen_dti_fa_additional_hads_depression_int)

glm_left_putamen_dti_md_additional_hads_depression_int <- lm(dti_md ~ additional_hads_depression * Group, left_putamen)
summary(glm_left_putamen_dti_md_additional_hads_depression_int) #associated with depression

glm_left_putamen_dki_kfa_additional_hads_depression_int <- lm(dki_kfa ~ additional_hads_depression * Group, left_putamen)
summary(glm_left_putamen_dki_kfa_additional_hads_depression_int)

glm_left_putamen_dki_mk_additional_hads_depression_int <- lm(dki_mk ~ additional_hads_depression * Group, left_putamen)
summary(glm_left_putamen_dki_mk_additional_hads_depression_int)

glm_left_putamen_fit_FWF_additional_hads_depression_int <- lm(fit_FWF ~ additional_hads_depression * Group, left_putamen)
summary(glm_left_putamen_fit_FWF_additional_hads_depression_int) #associated with depression

glm_left_putamen_fit_NDI_additional_hads_depression_int <- lm(fit_NDI ~ additional_hads_depression * Group, left_putamen)
summary(glm_left_putamen_fit_NDI_additional_hads_depression_int)

glm_left_putamen_fit_ODI_additional_hads_depression_int <- lm(fit_ODI ~ additional_hads_depression * Group, left_putamen)
summary(glm_left_putamen_fit_ODI_additional_hads_depression_int)

###bmi_cardio
glm_left_putamen_volume_bmi_cardio_int <- lm(volume ~ bmi_cardio * Group, left_putamen)
summary(glm_left_putamen_volume_bmi_cardio_int)

glm_left_putamen_dti_fa_bmi_cardio_int <- lm(dti_fa ~ bmi_cardio * Group, left_putamen)
summary(glm_left_putamen_dti_fa_bmi_cardio_int)

glm_left_putamen_dti_md_bmi_cardio_int <- lm(dti_md ~ bmi_cardio * Group, left_putamen)
summary(glm_left_putamen_dti_md_bmi_cardio_int) 

glm_left_putamen_dki_kfa_bmi_cardio_int <- lm(dki_kfa ~ bmi_cardio * Group, left_putamen)
summary(glm_left_putamen_dki_kfa_bmi_cardio_int)

glm_left_putamen_dki_mk_bmi_cardio_int <- lm(dki_mk ~ bmi_cardio * Group, left_putamen)
summary(glm_left_putamen_dki_mk_bmi_cardio_int) #associated with BMI

glm_left_putamen_fit_FWF_bmi_cardio_int <- lm(fit_FWF ~ bmi_cardio * Group, left_putamen)
summary(glm_left_putamen_fit_FWF_bmi_cardio_int) 

glm_left_putamen_fit_NDI_bmi_cardio_int <- lm(fit_NDI ~ bmi_cardio * Group, left_putamen)
summary(glm_left_putamen_fit_NDI_bmi_cardio_int)

glm_left_putamen_fit_ODI_bmi_cardio_int <- lm(fit_ODI ~ bmi_cardio * Group, left_putamen)
summary(glm_left_putamen_fit_ODI_bmi_cardio_int)

###age
glm_right_putamen_volume_age_int <- lm(volume ~ age * Group, right_putamen)
summary(glm_right_putamen_volume_age_int) #associated with age

glm_right_putamen_dti_fa_age_int <- lm(dti_fa ~ age * Group, right_putamen)
summary(glm_right_putamen_dti_fa_age_int)

glm_right_putamen_dti_md_age_int <- lm(dti_md ~ age * Group, right_putamen)
summary(glm_right_putamen_dti_md_age_int) #associated with age

glm_right_putamen_dki_kfa_age_int <- lm(dki_kfa ~ age * Group, right_putamen)
summary(glm_right_putamen_dki_kfa_age_int) #associated with age

glm_right_putamen_dki_mk_age_int <- lm(dki_mk ~ age * Group, right_putamen)
summary(glm_right_putamen_dki_mk_age_int)

glm_right_putamen_fit_FWF_age_int <- lm(fit_FWF ~ age * Group, right_putamen)
summary(glm_right_putamen_fit_FWF_age_int) #associated with age

glm_right_putamen_fit_NDI_age_int <- lm(fit_NDI ~ age * Group, right_putamen)
summary(glm_right_putamen_fit_NDI_age_int)

glm_right_putamen_fit_ODI_age_int <- lm(fit_ODI ~ age * Group, right_putamen)
summary(glm_right_putamen_fit_ODI_age_int) #associated with age

###homeint_storyrecall_d
glm_right_putamen_volume_homeint_storyrecall_d_int <- lm(volume ~ homeint_storyrecall_d * Group, right_putamen)
summary(glm_right_putamen_volume_homeint_storyrecall_d_int) #associated with memory, interaction effect
glm_right_putamen_volume_homeint_storyrecall_d_scd <- lm(volume ~ homeint_storyrecall_d, right_putamen_scd)
summary(glm_right_putamen_volume_homeint_storyrecall_d_scd) #sig
glm_right_putamen_volume_homeint_storyrecall_d_ctl <- lm(volume ~ homeint_storyrecall_d, right_putamen_ctl)
summary(glm_right_putamen_volume_homeint_storyrecall_d_ctl)


glm_right_putamen_dti_fa_homeint_storyrecall_d_int <- lm(dti_fa ~ homeint_storyrecall_d * Group, right_putamen)
summary(glm_right_putamen_dti_fa_homeint_storyrecall_d_int)

glm_right_putamen_dti_md_homeint_storyrecall_d_int <- lm(dti_md ~ homeint_storyrecall_d * Group, right_putamen)
summary(glm_right_putamen_dti_md_homeint_storyrecall_d_int) #associated with memory, group diff

glm_right_putamen_dki_kfa_homeint_storyrecall_d_int <- lm(dki_kfa ~ homeint_storyrecall_d * Group, right_putamen)
summary(glm_right_putamen_dki_kfa_homeint_storyrecall_d_int) #group diff

glm_right_putamen_dki_mk_homeint_storyrecall_d_int <- lm(dki_mk ~ homeint_storyrecall_d * Group, right_putamen)
summary(glm_right_putamen_dki_mk_homeint_storyrecall_d_int)

glm_right_putamen_fit_FWF_homeint_storyrecall_d_int <- lm(fit_FWF ~ homeint_storyrecall_d * Group, right_putamen)
summary(glm_right_putamen_fit_FWF_homeint_storyrecall_d_int) #group diff

glm_right_putamen_fit_NDI_homeint_storyrecall_d_int <- lm(fit_NDI ~ homeint_storyrecall_d * Group, right_putamen)
summary(glm_right_putamen_fit_NDI_homeint_storyrecall_d_int)

glm_right_putamen_fit_ODI_homeint_storyrecall_d_int <- lm(fit_ODI ~ homeint_storyrecall_d * Group, right_putamen)
summary(glm_right_putamen_fit_ODI_homeint_storyrecall_d_int) #associated with memory

###additional_hads_anxiety
glm_right_putamen_volume_additional_hads_anxiety_int <- lm(volume ~ additional_hads_anxiety * Group, right_putamen)
summary(glm_right_putamen_volume_additional_hads_anxiety_int)

glm_right_putamen_dti_fa_additional_hads_anxiety_int <- lm(dti_fa ~ additional_hads_anxiety * Group, right_putamen)
summary(glm_right_putamen_dti_fa_additional_hads_anxiety_int)

glm_right_putamen_dti_md_additional_hads_anxiety_int <- lm(dti_md ~ additional_hads_anxiety * Group, right_putamen)
summary(glm_right_putamen_dti_md_additional_hads_anxiety_int) #group diff

glm_right_putamen_dki_kfa_additional_hads_anxiety_int <- lm(dki_kfa ~ additional_hads_anxiety * Group, right_putamen)
summary(glm_right_putamen_dki_kfa_additional_hads_anxiety_int)

glm_right_putamen_dki_mk_additional_hads_anxiety_int <- lm(dki_mk ~ additional_hads_anxiety * Group, right_putamen)
summary(glm_right_putamen_dki_mk_additional_hads_anxiety_int)

glm_right_putamen_fit_FWF_additional_hads_anxiety_int <- lm(fit_FWF ~ additional_hads_anxiety * Group, right_putamen)
summary(glm_right_putamen_fit_FWF_additional_hads_anxiety_int) #group diff

glm_right_putamen_fit_NDI_additional_hads_anxiety_int <- lm(fit_NDI ~ additional_hads_anxiety * Group, right_putamen)
summary(glm_right_putamen_fit_NDI_additional_hads_anxiety_int)

glm_right_putamen_fit_ODI_additional_hads_anxiety_int <- lm(fit_ODI ~ additional_hads_anxiety * Group, right_putamen)
summary(glm_right_putamen_fit_ODI_additional_hads_anxiety_int)

###additional_hads_depression
glm_right_putamen_volume_additional_hads_depression_int <- lm(volume ~ additional_hads_depression * Group, right_putamen)
summary(glm_right_putamen_volume_additional_hads_depression_int)

glm_right_putamen_dti_fa_additional_hads_depression_int <- lm(dti_fa ~ additional_hads_depression * Group, right_putamen)
summary(glm_right_putamen_dti_fa_additional_hads_depression_int)

glm_right_putamen_dti_md_additional_hads_depression_int <- lm(dti_md ~ additional_hads_depression * Group, right_putamen)
summary(glm_right_putamen_dti_md_additional_hads_depression_int) #associated with depression

glm_right_putamen_dki_kfa_additional_hads_depression_int <- lm(dki_kfa ~ additional_hads_depression * Group, right_putamen)
summary(glm_right_putamen_dki_kfa_additional_hads_depression_int) #associated with depression

glm_right_putamen_dki_mk_additional_hads_depression_int <- lm(dki_mk ~ additional_hads_depression * Group, right_putamen)
summary(glm_right_putamen_dki_mk_additional_hads_depression_int)

glm_right_putamen_fit_FWF_additional_hads_depression_int <- lm(fit_FWF ~ additional_hads_depression * Group, right_putamen)
summary(glm_right_putamen_fit_FWF_additional_hads_depression_int) #associated with depression

glm_right_putamen_fit_NDI_additional_hads_depression_int <- lm(fit_NDI ~ additional_hads_depression * Group, right_putamen)
summary(glm_right_putamen_fit_NDI_additional_hads_depression_int)

glm_right_putamen_fit_ODI_additional_hads_depression_int <- lm(fit_ODI ~ additional_hads_depression * Group, right_putamen)
summary(glm_right_putamen_fit_ODI_additional_hads_depression_int)

###bmi_cardio
glm_right_putamen_volume_bmi_cardio_int <- lm(volume ~ bmi_cardio * Group, right_putamen)
summary(glm_right_putamen_volume_bmi_cardio_int)

glm_right_putamen_dti_fa_bmi_cardio_int <- lm(dti_fa ~ bmi_cardio * Group, right_putamen)
summary(glm_right_putamen_dti_fa_bmi_cardio_int)

glm_right_putamen_dti_md_bmi_cardio_int <- lm(dti_md ~ bmi_cardio * Group, right_putamen)
summary(glm_right_putamen_dti_md_bmi_cardio_int) #associated with depression

glm_right_putamen_dki_kfa_bmi_cardio_int <- lm(dki_kfa ~ bmi_cardio * Group, right_putamen)
summary(glm_right_putamen_dki_kfa_bmi_cardio_int) #associated with depression

glm_right_putamen_dki_mk_bmi_cardio_int <- lm(dki_mk ~ bmi_cardio * Group, right_putamen)
summary(glm_right_putamen_dki_mk_bmi_cardio_int)

glm_right_putamen_fit_FWF_bmi_cardio_int <- lm(fit_FWF ~ bmi_cardio * Group, right_putamen)
summary(glm_right_putamen_fit_FWF_bmi_cardio_int) #associated with depression

glm_right_putamen_fit_NDI_bmi_cardio_int <- lm(fit_NDI ~ bmi_cardio * Group, right_putamen)
summary(glm_right_putamen_fit_NDI_bmi_cardio_int)

glm_right_putamen_fit_ODI_bmi_cardio_int <- lm(fit_ODI ~ bmi_cardio * Group, right_putamen)
summary(glm_right_putamen_fit_ODI_bmi_cardio_int) #associated with BMI
