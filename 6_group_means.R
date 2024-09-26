source("1_data_preparation.R")
library(arsenal)

failed_qc <- c('sub-CC510255', #abnormality in left temporal pole
               'sub-CC620821', #segmentation errors from large ventricles
               'sub-CC621011', #segmentation errors from large ventricles
               'sub-CC711027', #severe motion artifacts in T1
               'sub-CC721434', #segmentation errors from large ventricles
               'sub-CC710551' #motion artifacts in DWI
)

#define remove_outliers function
remove_outliers <- function(x, na.rm = TRUE) 
{
  ## Find 25% and 75% Quantiles using inbuild function
  quant <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  
  ## Find Interquantile range and multiply it by 1.5 
  ## to derive factor for range calculation
  H <- 1.5 * IQR(x, na.rm = na.rm)
  
  y <- x
  
  ## fill the outlier elements with NA
  y[x < (quant[1] - H)] <- NA
  y[x > (quant[2] + H)] <- NA
  
  y
}

#extract SCD status from dwi_over_55 table
scd_status <- dwi_over_55 %>% select(participant_id, SCD) %>%
  filter(!participant_id %in% failed_qc)
scd_status$SCD <- factor(scd_status$SCD,
                         levels = c(1,0),
                         labels = c('SCD', 'Control'))

#import aseg stats table (gray matter volumes)
aseg = read_tsv("freesurfer/asegtable.tsv") %>%
  rename(participant_id=`Measure:volume`) %>%
  mutate(across(c(2:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) #normalize by intracranial volume
#import wmparc stats table (white matter volumes)
wmparc = read_tsv("freesurfer/wmparctable.tsv") %>%
  rename(participant_id=`Measure:volume`) %>%
  mutate(across(c(2:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) #normalize by intracranial volume
#create volumes table
volumes <- left_join(scd_status, aseg) %>% #join gray matter volumes with SCD status
  left_join(., wmparc) %>% #add white matter volumes
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
names(volumes) <- make.names(names(volumes))
volumes_table <- tableby(formulize('SCD', names(volumes)[3:134]), 
                         data = volumes, numeric.test="wt", total = F) 
summary(volumes_table, text = TRUE)

#import left aparc stats table (cortical thickness)
lh_aparc = read_tsv("freesurfer/lh_aparctable.tsv") %>%
  rename(participant_id=lh.aparc.thickness)
#import right aparc stats table (cortical thickness)
rh_aparc = read_tsv("freesurfer/rh_aparctable.tsv") %>%
  rename(participant_id=rh.aparc.thickness)
#create thickness table
thickness <- left_join(scd_status, lh_aparc) %>% #join left cortical thickness with SCD stats
  left_join(., rh_aparc) %>% #add right cortical thickness
  #remove "_thickness" from column names
  rename_with(~ str_remove(., "_thickness")) %>%
  mutate(across(where(is.double), remove_outliers))
thickness_table <- tableby(formulize('SCD', names(thickness)[3:74]), 
                           data = thickness, numeric.test="wt", total = F)
summary(thickness_table, text = T)

####diffusion group means
#read in diffusion tables
aparc2diff_files <- list.files(path = "freesurfer", 
                               pattern = "aparc.*aseg.*\\.*tsv", 
                               full.names = T)
for (i in 1:length(aparc2diff_files)) {
  assign(gsub(".tsv", "", 
              gsub("freesurfer.aparc.aseg2", "", make.names(aparc2diff_files[i]))), 
         read.delim(aparc2diff_files[i]))
}

fit_FWF <- fit_FWF %>% 
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
FWF_table <- tableby(formulize('SCD', names(fit_FWF)[4:106]),
                     data = fit_FWF, numeric.test="wt", total = FALSE)
summary(FWF_table, text = T)

fit_NDI <- fit_NDI %>% 
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
NDI_table <- tableby(formulize('SCD', names(fit_NDI)[4:106]),
                     data = fit_NDI, numeric.test="wt", total = FALSE)
summary(NDI_table, text = T)

fit_ODI <- fit_ODI %>% 
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
ODI_table <- tableby(formulize('SCD', names(fit_ODI)[4:106]),
                     data = fit_ODI, numeric.test="wt", total = FALSE)
summary(ODI_table, text = T)

dti_fa <- dti_fa %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers))
fa_table <- tableby(formulize('SCD', names(dti_fa)[4:106]),
                    data = dti_fa, numeric.test="wt", total=F)
summary(fa_table, text=T)

dti_md <- dti_md %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers))
md_table <- tableby(formulize('SCD', names(dti_md)[4:106]),
                    data = dti_md, numeric.test="wt", total=F)
summary(md_table, text=T)

dti_rd <- dti_rd %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers))
rd_table <- tableby(formulize('SCD', names(dti_rd)[4:106]),
                    data = dti_rd, numeric.test="wt", total=F)
summary(rd_table, text=T)

dti_ad <- dti_ad %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers))
ad_table <- tableby(formulize('SCD', names(dti_ad)[4:106]),
                    data = dti_ad, numeric.test="wt", total=F)
summary(ad_table, text=T)

dki_kfa <- dki_kfa %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers))
kfa_table <- tableby(formulize('SCD', names(dki_kfa)[4:106]),
                    data = dki_kfa, numeric.test="wt", total=F)
summary(kfa_table, text=T)

dki_mk <- dki_mk %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers))
mk_table <- tableby(formulize('SCD', names(dki_mk)[4:106]),
                     data = dki_mk, numeric.test="wt", total=F)
summary(mk_table, text=T)

dki_rk <- dki_rk %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers))
rk_table <- tableby(formulize('SCD', names(dki_rk)[4:106]),
                    data = dki_rk, numeric.test="wt", total=F)
summary(rk_table, text=T)

dki_ak <- dki_ak %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .) %>%
  mutate(across(where(is.double), remove_outliers))
ak_table <- tableby(formulize('SCD', names(dki_ak)[4:106]),
                    data = dki_ak, numeric.test="wt", total=F)
summary(ak_table, text=T)
