library(tidyverse)
library(labelVector)
#replace with location for your CamCan data
basedir = "/Volumes/labspace/AD/camcan995/"
data_dir = file.path(basedir, "source_materials")
raw_dir = file.path(basedir, "raw")

### pulling all spreadsheets into BIDS compliant participants.tsv ###

# standard_data = read_csv(file.path(data_dir, "standard_data.csv"))
# approved_data = read_tsv(file.path(data_dir, "behavioral/approved_data.tsv"))
# 
# #summary txt files had header and footer removed and were converted to tsv in a text editor
# vstm = read_tsv(file.path(data_dir, "behavioral/VSTMcolour_summary.tsv")) %>%
#   rename_at(vars(-CCID), ~ paste0(., "_vstm"))
# rtsimple = read_tsv(file.path(data_dir, "behavioral/RTsimple_summary.tsv")) %>%
#   rename_at(vars(-CCID), ~ paste0(., "_rtsimple"))
# rtchoice = read_tsv(file.path(data_dir, "behavioral/RTchoice_summary.tsv")) %>%
#   rename_at(vars(-CCID), ~ paste0(., "_rtchoice"))
# emotional_memory = read_tsv(file.path(data_dir, "behavioral/EmotionalMemory_summary.tsv")) %>%
#   rename_at(vars(-CCID), ~ paste0(., "_emotional_mem"))
# cattell = read_tsv(file.path(data_dir, "behavioral/Cattell_summary.tsv")) %>%
#   rename_at(vars(-CCID), ~ paste0(., "_cattell"))
# cardio_measures = read_tsv(file.path(data_dir, "behavioral/cardio_measures.tsv")) %>%
#   rename(CCID = Observations) %>%   rename_at(vars(-CCID), ~ paste0(., "_cardio"))
# 
# participants = full_join(standard_data, approved_data, by='CCID') %>%
#   left_join(., vstm, by='CCID') %>%
#   left_join(., rtsimple, by='CCID') %>% left_join(., rtchoice, by='CCID') %>%
#   left_join(., emotional_memory, by='CCID') %>% full_join(., cattell, by='CCID') %>%
#   left_join(., cardio_measures, by='CCID') %>%
#   mutate(CCID = str_replace(CCID, "CC", "sub-CC")) %>% rename(participant_id=CCID)
# names(participants) <- tolower(names(participants))
# write_tsv(participants, file.path(raw_dir, "participants.tsv"))

#replace with location of your participants.tsv
participants = read.delim(file.path(raw_dir, "participants.tsv"), tryLogical = FALSE)
participants$SCD <- participants$homeint_v230
participants$SCD[participants$SCD == 1] <- FALSE
participants$SCD[participants$SCD == 2] <- TRUE
participants$SCD[is.na(participants$SCD)] <- FALSE
participants$SCD <- factor(participants$SCD,
                           levels = c(1,0),
                           labels = c('SCD', 'Control'))
participants[participants$participant_id=="sub-CC610050", "mt_tr"] <- 30 #from json
participants[participants$participant_id=="sub-CC620821", "mt_tr"] <- 50 #from json
participants$mt_tr <- factor(participants$mt_tr,
                             levels = c(30, 50),
                             labels = c("TR=30ms", "TR=50ms"))
participants$Income <- factor(participants$homeint_v15,
                              levels = c("D", "B", "C", "A", "F", "E"),
                              labels = c("Less than  £18000", 
                                         "£18000 to 30999", 
                                         "£31000 to 51999",
                                         "£52000 to 100000",
                                         "Greater than £100000",
                                         "Prefer not to answer"))
participants$Ethnicity <- factor(participants$homeint_v24,
                                 levels = c(1,2,3,4,6),
                                 labels = c("White", "Mixed", "Asian", "Black", "Other"))
participants$Sex <- str_to_title(participants$sex)
participants$age_education_completed <- participants$homeint_v74
participants <- set_label(participants,
                          age = "Age (Years)",
                          age_education_completed = "Age Education Completed (Years)",
                          homeint_storyrecall_d = "Delayed Story Recall (Number of details)",
                          bp_dia_mean_cardio = "Diastolic Blood Pressure (mmHg)",
                          bp_sys_mean_cardio = "Systolic Blood Pressure (mmHg)",
                          pulse_mean_cardio = "Pulse (beats/min)",
                          additional_hads_anxiety = "HADS Anxiety Score",
                          additional_hads_depression = "HADS Depression Score")

#all DWI participants
#replace with location of your dwi participants.tsv
dwi_participants = read_tsv(file.path(data_dir, "imaging/dwi/participants.tsv")) %>%
  select(participant_id) %>%
  left_join(., participants, by='participant_id')
#DWI participants over age 55
dwi_over_55 = dwi_participants %>% filter(age > 55)
### split dwi_over_55 into SCD and controls ###
scd_dwi = dwi_over_55 %>% filter(SCD == "SCD")
write_tsv(scd_dwi, "dwi_over_55_scd.tsv")
ctl_dwi = dwi_over_55 %>% filter(SCD == "Control")
write_tsv(ctl_dwi, "dwi_over_55_ctl.tsv")

#all MTI participants
#replace with location of your mti participants.tsv
mti_participants = read_tsv(file.path(data_dir, "imaging/mti/participants.tsv")) %>%
  select(participant_id) %>%
  mutate(participant_id = str_replace(participant_id, "CC", "sub-CC")) %>%
  left_join(., participants, by='participant_id')
#MTI participants over age 55
mti_over_55 = mti_participants %>% filter(age > 55) %>% 
  filter(participant_id != "sub-CC410129") #error in scanning protocol, mti TR=34ms bl TR=30ms
#separate by TR
mti_over_55_tr50 <- mti_over_55 %>% filter(mt_tr == "TR=50ms")
scd_tr50 <- mti_over_55_tr50 %>% filter(SCD == "SCD")
write_tsv(scd_tr50, "mti_over_55_tr50_scd.tsv")
ctl_tr50 <- mti_over_55_tr50 %>% filter(SCD == "Control")
write_tsv(ctl_tr50, "mti_over_55_tr30_ctl.tsv")
mti_over_55_tr30 <- mti_over_55 %>% filter(mt_tr == "TR=30ms")
scd_tr30 <- mti_over_55_tr30 %>% filter(SCD == "SCD")
write_tsv(scd_tr30, "mti_over_55_tr30_scd.tsv")
ctl_tr30 <- mti_over_55_tr30 %>% filter(SCD == "Control")
write_tsv(ctl_tr30, "mti_over_55_tr30_ctl.tsv")

#dwi and mti participants
dwi_mti_over_55 <- full_join(dwi_over_55, mti_over_55)
write_tsv(dwi_mti_over_55, "dwi_mti_over_55.tsv")
