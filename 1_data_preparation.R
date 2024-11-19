library(tidyverse)
#replace with location for your CamCan data
basedir = "/Volumes/Research/lazarm03lab/labspace/AD/camcan995/"
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

#all DWI participants
#replace with location of your dwi participants.tsv
dwi_participants = read_tsv(file.path(data_dir, "imaging/dwi/participants.tsv")) %>%
  select(participant_id) %>%
  left_join(., participants, by='participant_id')
#DWI participants over age 55
dwi_over_55 = dwi_participants %>% filter(age > 55)

#all MTI participants
#replace with location of your mti participants.tsv
mti_participants = read_tsv(file.path(data_dir, "imaging/mti/participants.tsv")) %>%
  select(participant_id) %>%
  mutate(participant_id = str_replace(participant_id, "CC", "sub-CC")) %>%
  left_join(., participants, by='participant_id')
#MTI participants over age 55
mti_over_55 = mti_participants %>% filter(age > 55) %>% 
  filter(participant_id != "sub-CC410129") #error in scanning protocol, mti TR=34ms bl TR=30ms
mti_over_55[mti_over_55$participant_id=="sub-CC610050", "mt_tr"] <- 30 #from json
mti_over_55[mti_over_55$participant_id=="sub-CC620821", "mt_tr"] <- 50 #from json
write_tsv(mti_over_55, "mti_over_55.tsv")
#separate by TR
mti_over_55_tr50 <- mti_over_55 %>% filter(mt_tr == 50)
write_tsv(mti_over_55_tr50, "mti_over_55_tr50.tsv")
mti_over_55_tr30 <- mti_over_55 %>% filter(mt_tr == 30)
write_tsv(mti_over_55_tr30, "mti_over_55_tr30.tsv")

write_tsv(dwi_over_55, "dwi_over_55.tsv") #write to file